#include "../DISCON_MPC_CPP.cpp"

#include <iostream>

namespace
{
    void check(bool condition, const std::string& message)
    {
        if (!condition) throw std::runtime_error(message);
    }

    std::string readText(const std::filesystem::path& path)
    {
        std::ifstream in(path);
        check(static_cast<bool>(in), "Cannot read " + path.string());
        return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
    }

    std::string changed(std::string text, const std::string& from, const std::string& to)
    {
        const auto pos = text.find(from);
        check(pos != std::string::npos, "Missing fixture text: " + from);
        text.replace(pos, from.size(), to);
        return text;
    }

    void writeText(const std::filesystem::path& path, const std::string& text)
    {
        std::ofstream out(path);
        check(static_cast<bool>(out), "Cannot write test fixture");
        out << text;
    }

    struct Session
    {
        std::array<float, 3000> swap{};
        std::array<char, 1024> message{};
        int fail = 0;
        std::string input;
        std::string output;

        Session(const std::filesystem::path& config, const std::filesystem::path& out)
            : input(config.string()), output(out.string())
        {
            swap[48] = static_cast<float>(message.size());
            swap[49] = static_cast<float>(input.size());
            swap[50] = static_cast<float>(output.size());
            swap[3] = 0.1f;
            swap[19] = 120.0f;
            swap[20] = 1.2f;
            swap[26] = 12.0f;
            swap[1018] = 0.05f;
            swap[1020] = 0.1f;
        }

        void call(int status, float time)
        {
            swap[0] = static_cast<float>(status);
            swap[1] = time;
            DISCON(swap.data(), &fail, input.c_str(), output.c_str(), message.data());
        }
    };

    void checkLimits(float oldTorque, float oldPitch, float dt, const Session& session)
    {
        const float torque = session.swap[46], pitch = session.swap[44];
        check(std::isfinite(torque) && std::isfinite(pitch), "Nonfinite applied output");
        check(torque >= VS_MIN_TQ && torque <= VS_MAX_TQ, "Torque magnitude limit");
        check(pitch >= PC_MIN_PIT && pitch <= PC_MAX_PIT, "Pitch magnitude limit");
        check(std::fabs(torque - oldTorque) <= VS_MAX_TQ_RATE * dt + 0.05f, "Torque rate limit");
        check(std::fabs(pitch - oldPitch) <= PC_MAX_RAT * dt + 1.0e-6f, "Pitch rate limit");
    }

    void testRiccati()
    {
        gQOmega = gQX = gQV = gQTg = gQBeta = 1.0f;
        gRT = gRB = 1.0f;
        MatNN A{}, P{};
        MatNU B{};
        MatUN K{};
        for (int i = 0; i < N_STATE; ++i) A[i * N_STATE + i] = 0.7f;
        A[1] = 0.3f;
        A[N_STATE] = -0.15f;
        A[N_STATE + 2] = 0.2f;
        B[0] = 1.0f;
        B[N_CTRL + 1] = 0.5f;
        check(computeDiscreteLqrTerminal(A, B, P, K), "Stable LQR fixture failed to converge");
        MatNN closed = A;
        for (int i = 0; i < N_STATE; ++i)
            for (int j = 0; j < N_STATE; ++j)
                for (int u = 0; u < N_CTRL; ++u)
                    closed[i * N_STATE + j] += B[i * N_CTRL + u] * K[u * N_STATE + j];
        float residual = 0.0f;
        for (int i = 0; i < N_STATE; ++i)
            for (int j = 0; j < N_STATE; ++j)
            {
                float rhs = (i == j) ? 1.0f : 0.0f;
                for (int u = 0; u < N_CTRL; ++u) rhs += K[u * N_STATE + i] * K[u * N_STATE + j];
                for (int a = 0; a < N_STATE; ++a)
                    for (int b = 0; b < N_STATE; ++b)
                        rhs += closed[a * N_STATE + i] * P[a * N_STATE + b] * closed[b * N_STATE + j];
                residual = std::max(residual, std::fabs(P[i * N_STATE + j] - rhs));
            }
        check(residual < 0.002f, "LQR Bellman residual exposes an incorrect matrix product");
    }
}

int main(int argc, char** argv)
{
    try
    {
        check(argc == 2, "Pass the repository path");
        const auto root = std::filesystem::absolute(argv[1]);
        const auto temp = root / ".verification" / "testcases";
        std::filesystem::create_directories(temp);
        const std::string base = changed(readText(root / "DISCON_MPC_CPP.IN"),
            "TABLE_DIR = .", "TABLE_DIR = " + root.generic_string());
        const auto config = temp / "controller.IN";
        writeText(config, base);
        std::string error;
        check(loadGainTables(config.string(), error), error);
        check(gNPred == 330 && gNCtrlH == 70 && gTable.GTomega.size() == 29, "Default configuration");
        check(gTable.GTomega.front().size() == 23, "Default speed grid");
        int invalidCount = 0;
        auto invalid = [&](const std::string& text) {
            writeText(config, text);
            check(!loadGainTables(config.string(), error), "Invalid configuration accepted");
            check(!gTable.loaded && !error.empty(), "Invalid configuration retained loaded state");
            ++invalidCount;
        };
        invalid(base + "DEBUG = 1\n");
        invalid(base + "UNKNOWN = 1\n");
        invalid(changed(base, "MPC_DT = 0.010", "MPC_DT = nan"));
        invalid(changed(base, "MPC_DT = 0.010", "MPC_DT = 0.01garbage"));
        invalid(changed(base, "CTRL_DT = 0.010", "CTRL_DT = 0"));
        invalid(changed(base, "N_CTRL = 70", "N_CTRL = 331"));
        invalid(changed(base, "N_PRED = 330", "N_PRED = 1.5"));
        invalid(changed(base, "Q = 30.0", "Q = -30.0"));
        invalid(changed(base, "R = 5.0, 60.0", "R = 0, 60"));
        invalid(changed(base, "QP_START = cold", "QP_START = hot"));
        invalid(changed(base, "1.0,1.5,2.0", "1.0,1.0,2.0"));
        invalid(changed(base, "TABLE_DIR = " + root.generic_string(), "TABLE_DIR = missing-directory"));
        invalid("23\n29\n1,2\n");
        std::cout << "PASS: " << invalidCount << " invalid configurations rejected\n";

        Session bad(config, temp / "bad.out");
        bad.swap[46] = 12345.0f;
        bad.swap[44] = 0.123f;
        bad.call(0, 0.0f);
        check(bad.fail < 0 && !gState.initialized, "Invalid config must stop initialization");
        check(bad.swap[46] == 12345.0f && bad.swap[44] == 0.123f, "Failure must not write commands");
        bad.call(1, 1.0f);
        check(bad.fail < 0 && !gState.initialized, "Failed init must not become active on the next callback");
        std::cout << "PASS: initialization failure leaves output records untouched\n";

        testRiccati();
        std::cout << "PASS: coupled LQR Bellman residual\n";
        std::array<std::array<float, 2>, 2> finalOutputs{};
        for (int warm = 0; warm < 2; ++warm)
        {
            std::string fixture = changed(base, "N_PRED = 330", "N_PRED = 30");
            fixture = changed(fixture, "N_CTRL = 70", "N_CTRL = 5");
            fixture = changed(fixture, "DEBUG = 0", "DEBUG = 1");
            if (warm) fixture = changed(fixture, "QP_START = cold", "QP_START = warm");
            writeText(config, fixture);
            Session session(config, temp / (warm ? "warm.out" : "cold.out"));
            session.call(0, 0.0f);
            check(gState.initialized && session.fail >= 0, session.message.data());
            session.call(1, 29.99f);
            check(!gState.fractureActive, "Premature MPC switch");
            float oldTime = 29.99f;
            for (int step = 0; step < 41; ++step)
            {
                const float torque = session.swap[46], pitch = session.swap[44];
                const float time = 30.0f + 0.001f * step;
                session.call(1, time);
                check(session.fail == 0, session.message.data());
                check(gState.fractureActive && !gState.lastStepUsedFallback, "MPC solve failed");
                checkLimits(torque, pitch, time - oldTime, session);
                oldTime = time;
            }
            check(std::fabs(gState.lastMpcSolveTime - 30.04f) < 1.0e-5f, "10 ms solve gate slipped");
            check(gLastQpSolveMode == warm, "Cold/warm solver mode not used");
            const float torque = session.swap[46], pitch = session.swap[44];
            session.call(1, oldTime);
            check(session.swap[46] == torque && session.swap[44] == pitch, "Same-time callback advanced output");
            finalOutputs[warm] = {torque, pitch};
            gNWSR = 0;
            session.call(1, 30.05f);
            check(session.fail == 1 && gState.lastStepUsedFallback, "Forced QP failure did not invoke fallback");
            checkLimits(torque, pitch, 30.05f - oldTime, session);
            const float finalTorque = session.swap[46], finalPitch = session.swap[44];
            const float finalSolveTime = gState.lastMpcSolveTime;
            session.call(-1, 30.1f);
            check(!gState.initialized && session.fail == 0, "Final callback failed");
            check(gState.lastMpcSolveTime == finalSolveTime, "Final callback solved another QP");
            check(session.swap[46] == finalTorque && session.swap[44] == finalPitch, "Final callback changed outputs");
            check(std::filesystem::exists(temp / (warm ? "warm.mpc_command_trace.csv" : "cold.mpc_command_trace.csv")), "Debug trace missing");
            std::cout << "PASS: " << (warm ? "warm" : "cold") << " mode, solve gate, actuator limits, fallback and final callback\n";
        }
        check(std::fabs(finalOutputs[0][0] - finalOutputs[1][0]) < 1.0f &&
            std::fabs(finalOutputs[0][1] - finalOutputs[1][1]) < 1.0e-4f, "Cold and warm solutions disagree");

        writeText(config, base);
        Session full(config, temp / "default.out");
        full.call(0, 0.0f);
        check(full.fail >= 0, full.message.data());
        full.call(1, 29.99f);
        full.call(1, 30.0f);
        check(full.fail == 0 && !gState.lastStepUsedFallback, "Supplied 330/70 horizon failed");
        full.call(-1, 30.0f);
        std::cout << "PASS: supplied 330/70 horizon and cold/warm agreement\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}

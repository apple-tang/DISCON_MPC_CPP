#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include <qpOASES.hpp>

#ifdef _WIN32
#define DLL_EXPORT extern "C" __declspec(dllexport)
#else
#define DLL_EXPORT extern "C"
#endif

namespace
{
    USING_NAMESPACE_QPOASES

    struct ControllerState
    {
        bool initialized = false;
        bool fractureActive = false;
        float lastTime = 0.0f;
        float lastGenTorque = 0.0f;
        float pitchCmd = 0.0f;
        float lastPitchRate = 0.0f;
        float genSpeedF = 0.0f;
        float intSpdErr = 0.0f;
        float lastTimeVS = 0.0f;
        float lastTimePC = 0.0f;
        float VS_Slope15 = 0.0f;
        float VS_Slope25 = 0.0f;
        float VS_SySp = 0.0f;
        float VS_TrGnSp = 0.0f;
    };

    struct GainTableData
    {
        bool loaded = false;
        std::vector<float> speedPtsRpm;
        std::vector<float> windPtsMs;
        std::vector<std::vector<float>> GTomega;
        std::vector<std::vector<float>> GTbeta;
        std::vector<std::vector<float>> GTu;
        std::vector<std::vector<float>> GFomega;
        std::vector<std::vector<float>> GFbeta;
        std::vector<std::vector<float>> GFu;
    };

    ControllerState gState;
    GainTableData gTable;
    float gFractureTime = 30.0f;

    constexpr float D2R = 0.017453292f;
    constexpr float R2D = 57.295780f;
    constexpr float RPS2RPM = 9.5492966f;
    constexpr float PC_MAX_PIT = 1.570796f;
    constexpr float PC_MIN_PIT = 0.0f;
    constexpr float PC_MAX_RAT = 0.1396263f;   // rad/s
    constexpr float PC_DT      = 0.000125f;
    constexpr float PC_KP      = 0.01882681f;
    constexpr float PC_KI      = 0.008068634f;
    constexpr float PC_KK      = 0.1099965f;
    constexpr float PC_REFSPD  = 122.9096f;
    constexpr float CORNER_FREQ = 1.570796f;
    constexpr float ONE_PLUS_EPS = 1.0f + 1.1920929e-07f;
    constexpr float VS_MAX_TQ  = 47402.91f;    // N-m
    constexpr float VS_MIN_TQ  = 0.0f;
    constexpr float VS_MAX_TQ_RATE = 15000.0f; // N-m/s
    constexpr float VS_DT      = 0.000125f;
    constexpr float VS_CtInSp  = 70.16224f;
    constexpr float VS_Rgn2Sp  = 91.21091f;
    constexpr float VS_Rgn2K   = 2.332287f;
    constexpr float VS_SlPc    = 10.0f;
    constexpr float VS_Rgn3MP  = 0.01745329f;
    constexpr float VS_RtGnSp  = 121.6805f;
    constexpr float VS_RtPwr   = 5296610.0f;
    constexpr float OMEGA_REF = 0.0f;          // shutdown target rotor speed

    // Continuous-time physical parameters for the simplified shutdown MPC model.
    // Updated from the user's latest identified values.
    constexpr float J_RF = 3.19609962e7f;      // kg m^2  (identified rotor-equivalent inertia at 70%)
    constexpr float M_T  = 3.822772e5f;        // kg
    constexpr float C_T  = 6.858330e3f;        // N s/m
    constexpr float K_T  = 1.184008e6f;        // N/m
    constexpr float OMEGA_T = 1.759900f;       // rad/s
    constexpr float ZETA_T  = 5.097086e-3f;    // -
    constexpr float G_TOMEGA_DEFAULT = 0.0f;   // fallback until tables are loaded
    constexpr float G_FOMEGA_DEFAULT = 0.0f;
    constexpr float G_TBETA_DEFAULT = -2.0e6f;
    constexpr float G_FBETA_DEFAULT = 4.133385e4f;
    constexpr float G_TU_DEFAULT = 3.5e5f;
    constexpr float G_FU_DEFAULT = 2.0e4f;

    // Tunable MPC weights and state-constraint limits.
    float gQOmega = 22.9f;
    float gQX     = 40.0f;
    float gQV     = 80.0f;
    float gRT     = 1.0e-9f;
    float gRB     = 120.0f;
    int   gNPred  = 5;
    int   gNCtrlH = 5;
    float gOmegaErrMax  = 2.0f;   // rad/s
    float gTowerDispMax = 0.5f;   // m
    float gTowerVelMax  = 0.5f;   // m/s

    constexpr int   N_STATE = 5;
    constexpr int   N_CTRL  = 2;
    constexpr float BIG_NEG = -1.0e20f;

    inline std::string cArrayToString(const char* data, int n)
    {
        if (data == nullptr || n <= 0) return {};
        std::string s(data, data + n);
        auto pos = s.find('\0');
        if (pos != std::string::npos) s.resize(pos);
        return s;
    }

    inline void writeMessage(char* dst, int n, const std::string& msg)
    {
        if (dst == nullptr || n <= 0) return;
        std::fill(dst, dst + n, '\0');
        const int m = static_cast<int>(std::min<std::size_t>(msg.size(), static_cast<std::size_t>(n - 1)));
        std::memcpy(dst, msg.data(), static_cast<std::size_t>(m));
    }

    inline std::string trim(const std::string& s)
    {
        const auto b = s.find_first_not_of(" \t\r\n");
        if (b == std::string::npos) return {};
        const auto e = s.find_last_not_of(" \t\r\n");
        return s.substr(b, e - b + 1);
    }

    inline std::string stripComment(const std::string& s)
    {
        const auto p = s.find('!');
        return (p == std::string::npos) ? s : s.substr(0, p);
    }

    inline std::vector<float> splitCsvLine(const std::string& line)
    {
        std::vector<float> vals;
        std::stringstream ss(line);
        std::string item;
        while (std::getline(ss, item, ','))
        {
            item = trim(item);
            if (!item.empty())
                vals.push_back(std::stof(item));
        }
        return vals;
    }

    inline bool loadCsvTable(const std::string& path, int nWind, int nSpeed, std::vector<std::vector<float>>& table)
    {
        std::ifstream in(path);
        if (!in) return false;
        table.assign(nWind, std::vector<float>(nSpeed, 0.0f));
        std::string line;
        int row = 0;
        while (std::getline(in, line) && row < nWind)
        {
            line = trim(stripComment(line));
            if (line.empty()) continue;
            const auto vals = splitCsvLine(line);
            if (static_cast<int>(vals.size()) < nSpeed) return false;
            for (int j = 0; j < nSpeed; ++j) table[row][j] = vals[j];
            ++row;
        }
        return row == nWind;
    }

    inline bool loadGainTables(const std::string& inFilePath, std::string& err)
    {
        std::ifstream in(inFilePath);
        if (!in)
        {
            err = "Could not open MPC gain-table config file.";
            return false;
        }

        std::vector<std::string> lines;
        std::string line;
        while (std::getline(in, line))
        {
            line = trim(stripComment(line));
            if (!line.empty()) lines.push_back(line);
        }

        if (lines.size() < 10)
        {
            err = "MPC gain-table config file has too few non-empty lines.";
            return false;
        }

        const int nSpeed = std::stoi(lines[0]);
        const int nWind = std::stoi(lines[1]);
        gTable.speedPtsRpm = splitCsvLine(lines[2]);
        gTable.windPtsMs = splitCsvLine(lines[3]);

        if (static_cast<int>(gTable.speedPtsRpm.size()) != nSpeed || static_cast<int>(gTable.windPtsMs.size()) != nWind)
        {
            err = "Speed/wind point counts do not match config.";
            return false;
        }

        std::string dir = inFilePath;
        const auto pos = dir.find_last_of("\\/");
        dir = (pos == std::string::npos) ? "." : dir.substr(0, pos);
        auto makePath = [&](const std::string& file) { return dir + "\\" + trim(file); };

        if (!loadCsvTable(makePath(lines[4]), nWind, nSpeed, gTable.GTomega)) { err = "Failed loading GTomega table."; return false; }
        if (!loadCsvTable(makePath(lines[5]), nWind, nSpeed, gTable.GTbeta))  { err = "Failed loading GTbeta table.";  return false; }
        if (!loadCsvTable(makePath(lines[6]), nWind, nSpeed, gTable.GTu))     { err = "Failed loading GTu table.";     return false; }
        if (!loadCsvTable(makePath(lines[7]), nWind, nSpeed, gTable.GFomega)) { err = "Failed loading GFomega table."; return false; }
        if (!loadCsvTable(makePath(lines[8]), nWind, nSpeed, gTable.GFbeta))  { err = "Failed loading GFbeta table.";  return false; }
        if (!loadCsvTable(makePath(lines[9]), nWind, nSpeed, gTable.GFu))     { err = "Failed loading GFu table.";     return false; }

        if (lines.size() >= 11) gFractureTime = std::stof(lines[10]);

        // Backward-compatible parsing:
        // old format:  19 data lines, horizons fixed at 5/5 in code
        // new format:  21 data lines, lines[11:12] are N_PRED/N_CTRL_H
        std::size_t idx = 11;
        gNPred = 5;
        gNCtrlH = 5;
        if (lines.size() >= 21)
        {
            gNPred = std::stoi(lines[idx++]);
            gNCtrlH = std::stoi(lines[idx++]);
        }

        if (gNPred < 1)
        {
            err = "N_PRED must be at least 1.";
            return false;
        }
        if (gNCtrlH < 1)
        {
            err = "N_CTRL_H must be at least 1.";
            return false;
        }
        if (gNCtrlH > gNPred)
        {
            err = "N_CTRL_H must be less than or equal to N_PRED.";
            return false;
        }

        if (lines.size() > idx) gQOmega       = std::stof(lines[idx++]);
        if (lines.size() > idx) gQX           = std::stof(lines[idx++]);
        if (lines.size() > idx) gQV           = std::stof(lines[idx++]);
        if (lines.size() > idx) gRT           = std::stof(lines[idx++]);
        if (lines.size() > idx) gRB           = std::stof(lines[idx++]);
        if (lines.size() > idx) gOmegaErrMax  = std::stof(lines[idx++]);
        if (lines.size() > idx) gTowerDispMax = std::stof(lines[idx++]);
        if (lines.size() > idx) gTowerVelMax  = std::stof(lines[idx++]);

        gTable.loaded = true;
        return true;
    }

    inline float clampToRange(float x, float lo, float hi)
    {
        return std::max(lo, std::min(x, hi));
    }

    inline float interp2d(
        const std::vector<float>& windPts,
        const std::vector<float>& speedPts,
        const std::vector<std::vector<float>>& table,
        float windNow,
        float speedNow)
    {
        if (windPts.empty() || speedPts.empty() || table.empty()) return 0.0f;

        windNow = clampToRange(windNow, windPts.front(), windPts.back());
        speedNow = clampToRange(speedNow, speedPts.front(), speedPts.back());

        int iw = 0;
        while (iw + 1 < static_cast<int>(windPts.size()) && windPts[iw + 1] < windNow) ++iw;
        int is = 0;
        while (is + 1 < static_cast<int>(speedPts.size()) && speedPts[is + 1] < speedNow) ++is;

        const int iw2 = std::min(iw + 1, static_cast<int>(windPts.size()) - 1);
        const int is2 = std::min(is + 1, static_cast<int>(speedPts.size()) - 1);

        const float w1 = windPts[iw];
        const float w2 = windPts[iw2];
        const float s1 = speedPts[is];
        const float s2 = speedPts[is2];

        const float a = (iw2 == iw || std::fabs(w2 - w1) < 1e-8f) ? 0.0f : (windNow - w1) / (w2 - w1);
        const float b = (is2 == is || std::fabs(s2 - s1) < 1e-8f) ? 0.0f : (speedNow - s1) / (s2 - s1);

        const float g11 = table[iw ][is ];
        const float g21 = table[iw2][is ];
        const float g12 = table[iw ][is2];
        const float g22 = table[iw2][is2];

        return (1.0f - a) * (1.0f - b) * g11
             + a * (1.0f - b) * g21
             + (1.0f - a) * b * g12
             + a * b * g22;
    }

    inline void initializeBaselineStates(float time, float genSpeed, float bladePitch1)
    {
        gState.genSpeedF = genSpeed;
        gState.pitchCmd = bladePitch1;
        const float GK = 1.0f / (1.0f + gState.pitchCmd / PC_KK);
        gState.intSpdErr = gState.pitchCmd / (GK * PC_KI);
        gState.lastTime = time;
        gState.lastTimePC = time - PC_DT;
        gState.lastTimeVS = time - VS_DT;
        gState.VS_SySp = VS_RtGnSp / (1.0f + 0.01f * VS_SlPc);
        gState.VS_Slope15 = (VS_Rgn2K * VS_Rgn2Sp * VS_Rgn2Sp) / (VS_Rgn2Sp - VS_CtInSp);
        gState.VS_Slope25 = (VS_RtPwr / VS_RtGnSp) / (VS_RtGnSp - gState.VS_SySp);
        if (VS_Rgn2K == 0.0f)
            gState.VS_TrGnSp = gState.VS_SySp;
        else
            gState.VS_TrGnSp = (gState.VS_Slope25 - std::sqrt(gState.VS_Slope25 * (gState.VS_Slope25 - 4.0f * VS_Rgn2K * gState.VS_SySp))) / (2.0f * VS_Rgn2K);
    }

    inline void runBaselineDISCON(
        int iStatus,
        float time,
        float dt,
        float bladePitch1,
        float genSpeed,
        float& demandedGenTorque,
        float& demandedPitchCmd)
    {
        const float alpha = std::exp((gState.lastTime - time) * CORNER_FREQ);
        gState.genSpeedF = (1.0f - alpha) * genSpeed + alpha * gState.genSpeedF;

        // Baseline variable-speed torque control
        if ((time * ONE_PLUS_EPS - gState.lastTimeVS) >= VS_DT)
        {
            float elapTime = std::max(time - gState.lastTimeVS, 1.0e-6f);
            float genTrq = 0.0f;
            if ((gState.genSpeedF >= VS_RtGnSp) || (gState.pitchCmd >= VS_Rgn3MP))
                genTrq = VS_RtPwr / gState.genSpeedF;
            else if (gState.genSpeedF <= VS_CtInSp)
                genTrq = 0.0f;
            else if (gState.genSpeedF < VS_Rgn2Sp)
                genTrq = gState.VS_Slope15 * (gState.genSpeedF - VS_CtInSp);
            else if (gState.genSpeedF < gState.VS_TrGnSp)
                genTrq = VS_Rgn2K * gState.genSpeedF * gState.genSpeedF;
            else
                genTrq = gState.VS_Slope25 * (gState.genSpeedF - gState.VS_SySp);

            genTrq = std::min(genTrq, VS_MAX_TQ);
            if (iStatus == 0) gState.lastGenTorque = genTrq;
            float trqRate = (genTrq - gState.lastGenTorque) / elapTime;
            trqRate = std::clamp(trqRate, -VS_MAX_TQ_RATE, VS_MAX_TQ_RATE);
            gState.lastGenTorque = gState.lastGenTorque + trqRate * elapTime;
            gState.lastTimeVS = time;
        }

        // Baseline collective pitch PI control
        if ((time * ONE_PLUS_EPS - gState.lastTimePC) >= PC_DT)
        {
            float elapTime = std::max(time - gState.lastTimePC, 1.0e-6f);
            float GK = 1.0f / (1.0f + gState.pitchCmd / PC_KK);
            float spdErr = gState.genSpeedF - PC_REFSPD;
            gState.intSpdErr += spdErr * elapTime;
            gState.intSpdErr = std::clamp(gState.intSpdErr, PC_MIN_PIT / (GK * PC_KI), PC_MAX_PIT / (GK * PC_KI));
            float pitComP = GK * PC_KP * spdErr;
            float pitComI = GK * PC_KI * gState.intSpdErr;
            float pitComT = std::clamp(pitComP + pitComI, PC_MIN_PIT, PC_MAX_PIT);
            float pitRate = (pitComT - bladePitch1) / elapTime;
            pitRate = std::clamp(pitRate, -PC_MAX_RAT, PC_MAX_RAT);
            gState.pitchCmd = std::clamp(bladePitch1 + pitRate * elapTime, PC_MIN_PIT, PC_MAX_PIT);
            gState.lastPitchRate = pitRate;
            gState.lastTimePC = time;
        }

        demandedGenTorque = std::clamp(gState.lastGenTorque, VS_MIN_TQ, VS_MAX_TQ);
        demandedPitchCmd = std::clamp(gState.pitchCmd, PC_MIN_PIT, PC_MAX_PIT);
        (void)dt;
    }

    inline bool solveMultiStepMPC(
        float dt,
        float rotSpeed,
        float horWindV,
        float towerDispFA,
        float towerVelFA,
        float prevGenTorque,
        float prevPitchCmd,
        float& demandedGenTorque,
        float& demandedPitchCmd)
    {
        const int nPred = gNPred;
        const int nCtrlH = gNCtrlH;
        const float rotSpeedRPM = rotSpeed * 9.5492966f;
        const float windRef = gTable.loaded ? gTable.windPtsMs[std::min_element(gTable.windPtsMs.begin(), gTable.windPtsMs.end(),
            [&](float a, float b){ return std::fabs(a - horWindV) < std::fabs(b - horWindV); }) - gTable.windPtsMs.begin()] : 11.4f;
        const float dWind = horWindV - windRef;

        const float gTOmegaNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GTomega, horWindV, rotSpeedRPM) : G_TOMEGA_DEFAULT;
        const float gTBetaNow  = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GTbeta,  horWindV, rotSpeedRPM) : G_TBETA_DEFAULT;
        const float gTUNow     = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GTu,     horWindV, rotSpeedRPM) : G_TU_DEFAULT;
        const float gFOmegaNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GFomega, horWindV, rotSpeedRPM) : G_FOMEGA_DEFAULT;
        const float gFBetaNow  = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GFbeta,  horWindV, rotSpeedRPM) : G_FBETA_DEFAULT;
        const float gFUNow     = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GFu,     horWindV, rotSpeedRPM) : G_FU_DEFAULT;

        // Augmented state:
        // xbar = [ dOmega, x_t, v_t, u_prev_Tg, u_prev_beta ]^T
        const float x0 = rotSpeed - OMEGA_REF;
        const float x1 = towerDispFA;
        const float x2 = towerVelFA;
        const float x3 = prevGenTorque;
        const float x4 = prevPitchCmd;

        std::array<float, N_STATE> xbar0 = { x0, x1, x2, x3, x4 };

        // One-step Euler-discretized augmented model
        std::array<float, N_STATE * N_STATE> Abar{};
        std::array<float, N_STATE * N_CTRL> Bbar{};
        std::array<float, N_STATE> Ebar{};
        auto Aat = [&](int r, int c) -> float& { return Abar[r * N_STATE + c]; };
        auto Bat = [&](int r, int c) -> float& { return Bbar[r * N_CTRL + c]; };

        const float a11 = 1.0f + dt * gTOmegaNow / J_RF;
        const float a22 = 1.0f;
        const float a23 = dt;
        const float a31 = dt * gFOmegaNow / M_T;
        const float a32 = -dt * K_T / M_T;
        const float a33 = 1.0f - dt * C_T / M_T;
        const float b11 = -dt / J_RF;
        const float b12 =  dt * gTBetaNow / J_RF;
        const float b32 =  dt * gFBetaNow / M_T;
        const float e11 =  dt * gTUNow / J_RF;
        const float e31 =  dt * gFUNow / M_T;

        Aat(0,0) = a11;  Aat(0,3) = b11;  Aat(0,4) = b12;
        Aat(1,1) = a22;  Aat(1,2) = a23;
        Aat(2,0) = a31;  Aat(2,1) = a32;  Aat(2,2) = a33;  Aat(2,4) = b32;
        Aat(3,3) = 1.0f;
        Aat(4,4) = 1.0f;

        Bat(0,0) = b11;  Bat(0,1) = b12;
        Bat(1,0) = 0.0f; Bat(1,1) = 0.0f;
        Bat(2,0) = 0.0f; Bat(2,1) = b32;
        Bat(3,0) = 1.0f; Bat(3,1) = 0.0f;
        Bat(4,0) = 0.0f; Bat(4,1) = 1.0f;

        Ebar[0] = e11;
        Ebar[1] = 0.0f;
        Ebar[2] = e31;
        Ebar[3] = 0.0f;
        Ebar[4] = 0.0f;

        // Build prediction matrices X = F*x0 + G*U
        const int NX = N_STATE * nPred;
        const int NU = N_CTRL * nCtrlH;
        std::vector<float> F(NX * N_STATE, 0.0f);
        std::vector<float> G(NX * NU, 0.0f);

        auto matMulSq = [&](const std::array<float, N_STATE * N_STATE>& M,
                            const std::array<float, N_STATE * N_STATE>& N) {
            std::array<float, N_STATE * N_STATE> R{};
            for (int i = 0; i < N_STATE; ++i)
                for (int j = 0; j < N_STATE; ++j)
                    for (int k = 0; k < N_STATE; ++k)
                        R[i * N_STATE + j] += M[i * N_STATE + k] * N[k * N_STATE + j];
            return R;
        };

        auto matMulAB = [&](const std::array<float, N_STATE * N_STATE>& M,
                            const std::array<float, N_STATE * N_CTRL>& N) {
            std::array<float, N_STATE * N_CTRL> R{};
            for (int i = 0; i < N_STATE; ++i)
                for (int j = 0; j < N_CTRL; ++j)
                    for (int k = 0; k < N_STATE; ++k)
                        R[i * N_CTRL + j] += M[i * N_STATE + k] * N[k * N_CTRL + j];
            return R;
        };

        std::array<float, N_STATE * N_STATE> Apow{};
        for (int i = 0; i < N_STATE; ++i) Apow[i * N_STATE + i] = 1.0f;

        std::vector<std::array<float, N_STATE * N_STATE>> powers(static_cast<std::size_t>(nPred));
        for (int p = 0; p < nPred; ++p)
        {
            Apow = matMulSq(Abar, Apow);
            powers[p] = Apow;
            for (int i = 0; i < N_STATE; ++i)
                for (int j = 0; j < N_STATE; ++j)
                    F[(p * N_STATE + i) * N_STATE + j] = Apow[i * N_STATE + j];
        }

        for (int p = 0; p < nPred; ++p)
        {
            for (int c = 0; c <= p && c < nCtrlH; ++c)
            {
                std::array<float, N_STATE * N_CTRL> block{};
                if (p == c)
                {
                    block = Bbar;
                }
                else
                {
                    block = matMulAB(powers[p - c - 1], Bbar);
                }

                for (int i = 0; i < N_STATE; ++i)
                    for (int j = 0; j < N_CTRL; ++j)
                        G[(p * N_STATE + i) * NU + (c * N_CTRL + j)] = block[i * N_CTRL + j];
            }
        }

        // c = predicted state stack under zero control increments but
        //     with constant preview disturbance dWind over the horizon.
        std::vector<float> c(NX, 0.0f);
        std::array<float, N_STATE> xpred = xbar0;
        for (int p = 0; p < nPred; ++p)
        {
            std::array<float, N_STATE> xnext{};
            for (int i = 0; i < N_STATE; ++i)
            {
                for (int j = 0; j < N_STATE; ++j)
                    xnext[i] += Abar[i * N_STATE + j] * xpred[j];
                xnext[i] += Ebar[i] * dWind;
                c[p * N_STATE + i] = xnext[i];
            }
            xpred = xnext;
        }

        // Block-diagonal Q and R
        std::vector<float> Qblk(NX * NX, 0.0f);
        std::vector<float> Rblk(NU * NU, 0.0f);

        for (int p = 0; p < nPred; ++p)
        {
            const int base = p * N_STATE;
            Qblk[(base + 0) * NX + (base + 0)] = gQOmega;
            Qblk[(base + 1) * NX + (base + 1)] = gQX;
			Qblk[(base + 2) * NX + (base + 2)] = gQV;
        }
        for (int p = 0; p < nCtrlH; ++p)
        {
            const int base = p * N_CTRL;
            Rblk[(base + 0) * NU + (base + 0)] = gRT;
            Rblk[(base + 1) * NU + (base + 1)] = gRB;
        }

        // H = 2*(G'QG + R), g = 2*G'Qc
        std::vector<real_t> H(NU * NU, 0.0);
        std::vector<real_t> g(NU, 0.0);

        for (int i = 0; i < NU; ++i)
        {
            for (int j = 0; j < NU; ++j)
            {
                float val = 0.0f;
                for (int k = 0; k < NX; ++k)
                {
                    const float qkk = Qblk[k * NX + k];
                    if (qkk != 0.0f)
                        val += G[k * NU + i] * qkk * G[k * NU + j];
                }
                val += Rblk[i * NU + j];
                H[i * NU + j] = static_cast<real_t>(2.0f * val);
            }
            float gv = 0.0f;
            for (int k = 0; k < NX; ++k)
            {
                const float qkk = Qblk[k * NX + k];
                if (qkk != 0.0f)
                    gv += G[k * NU + i] * qkk * c[k];
            }
            g[i] = static_cast<real_t>(2.0f * gv);
        }

        // Bounds on input increments
        const float dTgRate = VS_MAX_TQ_RATE * dt;
        const float dBetaRate = PC_MAX_RAT * dt;

        std::vector<real_t> lb(NU, 0.0), ub(NU, 0.0);
        for (int p = 0; p < nCtrlH; ++p)
        {
            lb[p * N_CTRL + 0] = static_cast<real_t>(-dTgRate);
            ub[p * N_CTRL + 0] = static_cast<real_t>( dTgRate);
            lb[p * N_CTRL + 1] = static_cast<real_t>(-dBetaRate);
            ub[p * N_CTRL + 1] = static_cast<real_t>( dBetaRate);
        }

        // Absolute-input constraints through cumulative-sum matrix Tu
        const int NC_INPUT = 2 * NU;
        const int NC_STATE = 2 * 3 * nPred; // upper/lower bounds for [dOmega, x_t, v_t]
        const int NC = NC_INPUT + NC_STATE;

        std::vector<real_t> Acon(NC * NU, 0.0), lbA(NC, BIG_NEG), ubA(NC, 0.0);
        for (int rowBlk = 0; rowBlk < nCtrlH; ++rowBlk)
        {
            for (int colBlk = 0; colBlk <= rowBlk; ++colBlk)
            {
                Acon[(rowBlk * N_CTRL + 0) * NU + (colBlk * N_CTRL + 0)] = 1.0;
                Acon[(rowBlk * N_CTRL + 1) * NU + (colBlk * N_CTRL + 1)] = 1.0;
                Acon[((rowBlk + nCtrlH) * N_CTRL + 0) * NU + (colBlk * N_CTRL + 0)] = -1.0;
                Acon[((rowBlk + nCtrlH) * N_CTRL + 1) * NU + (colBlk * N_CTRL + 1)] = -1.0;
            }

            ubA[rowBlk * N_CTRL + 0] = static_cast<real_t>(VS_MAX_TQ - prevGenTorque);
            ubA[rowBlk * N_CTRL + 1] = static_cast<real_t>(PC_MAX_PIT - prevPitchCmd);
            ubA[(rowBlk + nCtrlH) * N_CTRL + 0] = static_cast<real_t>(prevGenTorque - VS_MIN_TQ);
            ubA[(rowBlk + nCtrlH) * N_CTRL + 1] = static_cast<real_t>(prevPitchCmd - PC_MIN_PIT);
        }

        // State constraints on predicted [dOmega, x_t, v_t]
        const int stateRowBase = NC_INPUT;
        for (int p = 0; p < nPred; ++p)
        {
            const int xBase = p * N_STATE;
            const int upRow = stateRowBase + p * 3;
            const int lowRow = stateRowBase + 3 * nPred + p * 3;

            // Upper-bound rows: G_state * U <= xmax - c_state
            for (int j = 0; j < NU; ++j)
            {
                Acon[(upRow + 0) * NU + j] = static_cast<real_t>(G[(xBase + 0) * NU + j]); // dOmega
                Acon[(upRow + 1) * NU + j] = static_cast<real_t>(G[(xBase + 1) * NU + j]); // x_t
                Acon[(upRow + 2) * NU + j] = static_cast<real_t>(G[(xBase + 2) * NU + j]); // v_t

                Acon[(lowRow + 0) * NU + j] = static_cast<real_t>(-G[(xBase + 0) * NU + j]);
                Acon[(lowRow + 1) * NU + j] = static_cast<real_t>(-G[(xBase + 1) * NU + j]);
                Acon[(lowRow + 2) * NU + j] = static_cast<real_t>(-G[(xBase + 2) * NU + j]);
            }

            const float cOmega = c[xBase + 0];
            const float cDisp  = c[xBase + 1];
            const float cVel   = c[xBase + 2];

            ubA[upRow + 0] = static_cast<real_t>( gOmegaErrMax - cOmega );
            ubA[upRow + 1] = static_cast<real_t>( gTowerDispMax - cDisp );
            ubA[upRow + 2] = static_cast<real_t>( gTowerVelMax  - cVel );

            ubA[lowRow + 0] = static_cast<real_t>( gOmegaErrMax + cOmega );
            ubA[lowRow + 1] = static_cast<real_t>( gTowerDispMax + cDisp );
            ubA[lowRow + 2] = static_cast<real_t>( gTowerVelMax  + cVel );
        }

        // create and solve QP
        QProblem qp(NU, NC);
        Options options;
        options.printLevel = PL_NONE;
        qp.setOptions(options);

        int_t nWSR = 200;
        returnValue rv = qp.init(
            H.data(),
            g.data(),
            Acon.data(),
            lb.data(),
            ub.data(),
            lbA.data(),
            ubA.data(),
            nWSR
        );
        if (rv != SUCCESSFUL_RETURN)
            return false;

        std::vector<real_t> xOpt(NU, 0.0);
        qp.getPrimalSolution(xOpt.data());

        const float dTg = static_cast<float>(xOpt[0]);
        const float dBeta = static_cast<float>(xOpt[1]);

        demandedGenTorque = std::clamp(prevGenTorque + dTg, VS_MIN_TQ, VS_MAX_TQ);
        demandedPitchCmd  = std::clamp(prevPitchCmd + dBeta, PC_MIN_PIT, PC_MAX_PIT);
        return true;
    }
}

// Minimal C++ DISCON shell.
// This version is intentionally simple:
// - reads the same avrSWAP layout as the Fortran DISCON
// - keeps DLL linkage compatible with ServoDyn/Bladed interface
// - provides a clean location to insert qpOASES-based MPC later
DLL_EXPORT void DISCON(float* avrSWAP, int* aviFAIL, const char* accINFILE, const char* avcOUTNAME, char* avcMSG)
{
    if (avrSWAP == nullptr || aviFAIL == nullptr || avcMSG == nullptr)
        return;

    const int msgLen = std::max(1, static_cast<int>(std::lround(avrSWAP[48])));   // avrSWAP(49) in Fortran
    const int inFileLen = std::max(0, static_cast<int>(std::lround(avrSWAP[49]))); // avrSWAP(50)
    const int outNameLen = std::max(0, static_cast<int>(std::lround(avrSWAP[50])));// avrSWAP(51)

    const int iStatus = static_cast<int>(std::lround(avrSWAP[0]));  // avrSWAP(1)
    const float time = avrSWAP[1];                                  // avrSWAP(2)
    const float bladePitch1 = avrSWAP[3];                           // avrSWAP(4)
    const float genSpeed = avrSWAP[19];                             // avrSWAP(20)
    const float rotSpeed = avrSWAP[20];                             // avrSWAP(21)
    const float horWindV = avrSWAP[26];                             // avrSWAP(27)
    const float towerVelFA = avrSWAP[1018];                         // avrSWAP(1019)
    const float towerVelSS = avrSWAP[1019];                         // avrSWAP(1020)
    const float towerDispFA = avrSWAP[1020];                        // avrSWAP(1021)

    (void)genSpeed;
    (void)rotSpeed;
    (void)horWindV;
    (void)towerVelFA;
    (void)towerVelSS;
    (void)towerDispFA;
    (void)cArrayToString(accINFILE, inFileLen);
    (void)cArrayToString(avcOUTNAME, outNameLen);

    if (iStatus == 0 || !gState.initialized)
    {
        gState.initialized = true;
        gState.fractureActive = false;
        gState.lastTime = time;
        gState.lastGenTorque = 0.0f;
        gState.pitchCmd = bladePitch1;
        gState.lastPitchRate = 0.0f;

        std::string loadErr;
        const std::string inFile = cArrayToString(accINFILE, inFileLen);
        const bool tablesLoaded = loadGainTables(inFile, loadErr);

        if (tablesLoaded)
            initializeBaselineStates(time, genSpeed, bladePitch1);

        *aviFAIL = tablesLoaded ? 1 : -1;
        if (tablesLoaded)
        {
            std::ostringstream oss;
            oss << "Running C++ DISCON shell: baseline control before fracture, MPC after fracture"
                << " (N_PRED=" << gNPred << ", N_CTRL_H=" << gNCtrlH << ").";
            writeMessage(avcMSG, msgLen, oss.str());
        }
        else
            writeMessage(avcMSG, msgLen, loadErr);
    }
    else
    {
        *aviFAIL = 0;
        writeMessage(avcMSG, msgLen, "");
    }

    const float dt = std::max(time - gState.lastTime, 1.0e-4f);

    if (!gState.fractureActive && time >= gFractureTime)
        gState.fractureActive = true;

    float demandedGenTorque = 0.0f;
    float demandedPitch = bladePitch1;

    if (!gState.fractureActive)
    {
        runBaselineDISCON(
            iStatus,
            time,
            dt,
            bladePitch1,
            genSpeed,
            demandedGenTorque,
            demandedPitch
        );
    }
    else
    {
        const bool qpSolved = solveMultiStepMPC(
            dt,
            rotSpeed,
            horWindV,
            towerDispFA,
            towerVelFA,
            gState.lastGenTorque,
            gState.pitchCmd,
            demandedGenTorque,
            demandedPitch
        );

        if (!qpSolved)
        {
            *aviFAIL = -1;
            writeMessage(avcMSG, msgLen, "qpOASES multi-step MPC QP failed inside DISCON.");
            return;
        }
    }

    const float demandedPitchRate = std::clamp((demandedPitch - gState.pitchCmd) / dt, -PC_MAX_RAT, PC_MAX_RAT);

    avrSWAP[34] = 1.0f;            // avrSWAP(35)  generator contactor: main variable-speed generator
    avrSWAP[46] = demandedGenTorque; // avrSWAP(47) demanded generator torque
    avrSWAP[55] = 0.0f;            // avrSWAP(56) torque override = yes

    avrSWAP[41] = demandedPitch;   // avrSWAP(42) blade 1 pitch command
    avrSWAP[42] = demandedPitch;   // avrSWAP(43) blade 2 pitch command
    avrSWAP[43] = demandedPitch;   // avrSWAP(44) blade 3 pitch command
    avrSWAP[44] = demandedPitch;   // avrSWAP(45) collective pitch command
    avrSWAP[54] = 0.0f;            // avrSWAP(55) pitch override = yes
    avrSWAP[45] = demandedPitchRate; // avrSWAP(46) demanded collective pitch rate

    gState.lastTime = time;
    gState.lastGenTorque = demandedGenTorque;
    gState.pitchCmd = demandedPitch;
    gState.lastPitchRate = demandedPitchRate;

    std::ostringstream oss;
    if (!gState.fractureActive)
    {
        oss << "baseline control active: Tg=" << demandedGenTorque
            << " Nm, beta=" << demandedPitch * R2D
            << " deg, dBeta=" << demandedPitchRate * R2D << " deg/s";
    }
    else
    {
        oss << "MPC shutdown active: Tg=" << demandedGenTorque
            << " Nm, beta=" << demandedPitch * R2D
            << " deg, dBeta=" << demandedPitchRate * R2D
            << " deg/s, wind=" << horWindV << " m/s";
    }
    writeMessage(avcMSG, msgLen, oss.str());
}

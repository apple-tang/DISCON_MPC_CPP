#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <ctime>

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
        float targetGenTorque = 0.0f;
        float targetPitchCmd = 0.0f;
        float lastPitchRate = 0.0f;
        float genSpeedF = 0.0f;
        float intSpdErr = 0.0f;
        float lastTimeVS = 0.0f;
        float lastTimePC = 0.0f;
        float lastLidarLogTime = -1.0f;
        int fallbackCount = 0;
        bool lastStepUsedFallback = false;
        float VS_Slope15 = 0.0f;
        float VS_Slope25 = 0.0f;
        float VS_SySp = 0.0f;
        float VS_TrGnSp = 0.0f;
        bool operatingWindRefInitialized = false;
        float operatingWindRef = 0.0f;
        int turbulentWindBand = -1;
        float lastTuningTraceTime = -1.0f;
        bool lastTuningTraceFallbackState = false;
        float lastCommandTraceTime = -1.0f;
        bool mpcCommandInitialized = false;
        float lastMpcSolveTime = 0.0f;
        bool hasOneStepPrediction = false;
        float oneStepPredictionIssuedTime = 0.0f;
        float oneStepPredictionTime = 0.0f;
        float oneStepPredictionDt = 0.0f;
        float oneStepPredOmegaRef = 0.0f;
        float oneStepPredDOmega = 0.0f;
        float oneStepPredTowerDispFA = 0.0f;
        float oneStepPredTowerVelFA = 0.0f;
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

    struct ReferenceCurveTableData
    {
        bool loaded = false;
        std::vector<float> windPtsMs;
        std::vector<float> tgHiKnm;
        std::vector<float> omegaHiRpm;
        std::vector<float> omegaLoRpm;
        std::vector<float> shapeP;
    };

    struct ShutdownReferenceTrajectoryData
    {
        bool loaded = false;
        std::vector<float> timeS;
        std::vector<float> omegaRefRpm;
        std::vector<float> betaRefDeg;
        std::vector<float> tgRefKnm;
    };

    struct LidarPreviewData
    {
        bool available = false;
        int sensorType = 0;
        int numBeams = 0;
        int numPulseGates = 0;
        float referenceWind = 0.0f;
        std::vector<float> measuredSpeeds;
        std::vector<float> posX;
        std::vector<float> posY;
        std::vector<float> posZ;
    };

    struct LidarHistoryPoint
    {
        float tMeas = 0.0f;
        float convectionWind = 0.0f;
        float measuredSpeed = 0.0f;
        float posX = 0.0f;
    };

    ControllerState gState;
    GainTableData gTable;
    ReferenceCurveTableData gRefCurve;
    ShutdownReferenceTrajectoryData gShutdownRef;
    LidarPreviewData gLidar;
    std::deque<LidarHistoryPoint> gLidarHistory;
    float gFractureTime = 30.0f;
    std::string gDebugLogPath;
    std::string gPredictionLogPath;
    std::string gInputTracePath;
    std::string gTuningTracePath;
    std::string gCommandTracePath;
#ifdef MPC_ENABLE_TIMING
    std::string gTimingLogPath;
#endif
    bool gEnableTraceFiles = true;
    bool gEnablePredictionTrace = false;
    bool gEnableInputTrace = false;
    bool gEnableTuningTrace = false;
    int gCurrentTerminalIndex = -1;

    constexpr float D2R = 0.017453292f;
    constexpr float R2D = 57.295780f;
    constexpr float RPS2RPM = 9.5492966f;
    constexpr float PC_MAX_PIT = 1.570796f;
    constexpr float PC_MIN_PIT = 0.0f;
    constexpr float PC_MAX_RAT = 0.1396263f;   // rad/s
    constexpr float PC_DT = 0.000125f;
    constexpr float PC_KP = 0.01882681f;
    constexpr float PC_KI = 0.008068634f;
    constexpr float PC_KK = 0.1099965f;
    constexpr float PC_REFSPD = 122.9096f;
    constexpr float CORNER_FREQ = 1.570796f;
    constexpr float ONE_PLUS_EPS = 1.0f + 1.1920929e-07f;
    constexpr float VS_MAX_TQ = 47402.91f;    // N-m
    constexpr float VS_MIN_TQ = 0.0f;
    constexpr float VS_MAX_TQ_RATE = 15000.0f; // N-m/s
    constexpr float VS_DT = 0.000125f;
    constexpr float VS_CtInSp = 70.16224f;
    constexpr float VS_Rgn2Sp = 91.21091f;
    constexpr float VS_Rgn2K = 2.332287f;
    constexpr float VS_SlPc = 10.0f;
    constexpr float VS_Rgn3MP = 0.01745329f;
    constexpr float VS_RtGnSp = 121.6805f;
    constexpr float VS_RtPwr = 5296610.0f;
    constexpr float OMEGA_REF = 0.0f;          // shutdown target rotor speed
    constexpr float OMEGA_MIN_RPM = 0.0f;      // predicted rotor-speed lower bound [rpm]
    constexpr float OMEGA_MIN = OMEGA_MIN_RPM / RPS2RPM; // predicted rotor-speed lower bound [rad/s]
    constexpr float TG_REF = VS_MIN_TQ;        // shutdown target generator torque [N-m]
    constexpr float BETA_REF = PC_MAX_PIT;     // shutdown target collective pitch [rad]
    constexpr float TERMINAL_UNLOAD_TRIGGER_RPM_DEFAULT = 0.2f;

    // Continuous-time physical parameters for the simplified shutdown MPC model.
    // Updated from the user's latest identified values.
    constexpr float J_RF_DEFAULT = 2.56488355e7f; // kg m^2, fallback until OpenFAST publishes fracture inertia
    constexpr float J_generator = 534.116f;    // Generator inertia about HSS (kg m^2)

    constexpr float N_gear = 97.0f;            // gearbox ratio
    constexpr float J_EQ_DEFAULT = J_RF_DEFAULT + J_generator * N_gear * N_gear; // 等效转动惯量 fallback
    constexpr float M_T = 3.822772e5f;        // kg
    constexpr float C_T = 6.858330e3f;        // N s/m
    constexpr float K_T = 1.184008e6f;        // N/m
    constexpr float OMEGA_T = 1.759900f;       // rad/s
    constexpr float ZETA_T = 5.097086e-3f;    // -
    constexpr float G_TOMEGA_DEFAULT = 0.0f;   // fallback until tables are loaded
    constexpr float G_FOMEGA_DEFAULT = 0.0f;
    constexpr float G_TBETA_DEFAULT = -2.0e6f;
    constexpr float G_FBETA_DEFAULT = 4.133385e4f;
    constexpr float G_TU_DEFAULT = 3.5e5f;
    constexpr float G_FU_DEFAULT = 2.0e4f;
    constexpr int   FRACTURE_J_RF_IDX = 1022;           // C index for avrSWAP(1023)
    constexpr int   FRACTURE_MASS_IMBALANCE_IDX = 1023; // C index for avrSWAP(1024)
    constexpr int   FRACTURE_LOCATION_IDX = 1024;       // C index for avrSWAP(1025)
    constexpr int   FRACTURE_SIGMA_IDX = 1025;          // C index for avrSWAP(1026)
    constexpr int   FRACTURE_STATUS_IDX = 1026;         // C index for avrSWAP(1027): 0 inactive, 1 ramping, 2 final
    constexpr int   LIDAR_MSR_START = 2000;  // C index for avrSWAP(2001)
    constexpr int   LIDAR_MAX_CHAN  = 500;
    constexpr float LIDAR_LOG_DT = 0.1f;     // seconds between lidar debug log entries
    constexpr float LIDAR_MIN_CONV_WIND = 1.0f;   // minimum convection speed used in Taylor mapping [m/s]
    constexpr float LIDAR_EVOLUTION_LENGTH = 300.0f; // decay length for preview confidence [m]
    constexpr float LIDAR_HISTORY_MAX_AGE = 30.0f;   // seconds of lidar history kept for preview mapping
    constexpr float OPERATING_WIND_REF_TAU = 1.0f / 1.570796f; // matched to baseline DISCON CornerFreq [s]
    constexpr float TUNING_TRACE_DT = 0.5f;          // seconds between ordinary tuning-trace rows
    constexpr float COMMAND_TRACE_DT = 0.01f;        // dataset-aligned command-trace sampling interval [s]
    constexpr float TURB_WIND_BAND_HYST = 0.5f;      // [m/s] hysteresis on filtered mean wind band switching
    constexpr float TURB_WIND_BAND_0_MAX = 9.0f;     // 8 m/s anchor
    constexpr float TURB_WIND_BAND_1_MAX = 13.0f;    // 10/12 m/s anchor block
    constexpr float TURB_WIND_BAND_2_MAX = 17.0f;    // 14/16 m/s anchor block
    constexpr float TURB_WIND_BAND_3_MAX = 21.0f;    // 18/20 m/s anchor block
    constexpr float TURB_OU_BAND0 = 0.05f;           // 8 m/s
    constexpr float TURB_OU_BAND1 = 0.20f;           // 10/12 m/s
    constexpr float TURB_OU_BAND2 = 1.50f;           // 14/16 m/s
    constexpr float TURB_OU_BAND3 = 0.20f;           // 18/20 m/s
    constexpr float TURB_OU_BAND4 = 2.00f;           // 22/24 m/s

    float gJRFRuntime = J_RF_DEFAULT;
    float gJEQRuntime = J_EQ_DEFAULT;
    float gFractureMassImbalance = 0.0f;
    float gFractureLocationRR = 0.0f;
    float gFractureSigma = 1.0f;
    int gFractureStatusFromFAST = 0;
    bool gFractureInertiaLocked = false;

    // Tunable MPC weights and state-constraint limits.
    float gQOmega = 22.9f;
    float gQX = 40.0f;
    float gQV = 80.0f;
    float gQTg = 1.0e-6f;
    float gQBeta = 10.0f;
    float gRT = 1.0e-9f;
    float gRB = 120.0f;
    int   gNPred = 5;
    int   gNCtrlH = 5;
    int   gNWSR = 200;
    float gPredictionDt = 0.000125f;
    float gTerminalCostScale = 1.0f;
    float gTerminalUnloadTriggerRpm = TERMINAL_UNLOAD_TRIGGER_RPM_DEFAULT;
    int   gActuatorOrder = 0;
    float gPitchActuatorTau = 0.25f;
    float gTorqueActuatorTau = 0.25f;
    float gControlDt = 0.010f;
    bool  gEnableQpHotstart = false;

    constexpr int   N_STATE = 5;
    constexpr int   N_CTRL = 2;
    constexpr float BIG_NEG = -1.0e20f;
    constexpr float BIG_POS = 1.0e20f;

    using MatNN = std::array<float, N_STATE * N_STATE>;
    using MatNU = std::array<float, N_STATE * N_CTRL>;
    using MatUN = std::array<float, N_CTRL * N_STATE>;
    using MatUU = std::array<float, N_CTRL * N_CTRL>;

    inline float getGeneratorTorqueReference(float time, float rotSpeed, float windNow);
    inline float getRotorSpeedReference(float time);
    inline float getPitchReference(float time);
    inline float actuatorRetention(float dtController, float tau);

    struct TerminalScheduleData
    {
        bool ready = false;
        int fallbackPoints = 0;
        std::vector<float> windPtsMs;
        std::vector<float> speedPtsRpm;
        std::vector<MatNN> PTable;
        std::vector<MatUN> KTable;
    };

    TerminalScheduleData gTerminal;

    struct SolverWorkspace
    {
        int nx = 0;
        int nu = 0;
        int nc = 0;
        std::vector<float> G;
        std::vector<float> c;
        std::vector<real_t> H;
        std::vector<real_t> g;
        std::vector<real_t> lb;
        std::vector<real_t> ub;
        std::vector<real_t> Acon;
        std::vector<real_t> lbA;
        std::vector<real_t> ubA;
        std::vector<real_t> xOpt;
        std::vector<float> WG;
        std::vector<MatNN> powers;
    };

    SolverWorkspace gSolverWs;
    std::unique_ptr<SQProblem> gHotstartQp;
    bool gHotstartQpReady = false;
    int gLastQpSolveMode = 0; // 0=cold, 1=hot, 2=hot-failed-then-cold
    int gLastQpNwsrUsed = 0;

#ifdef MPC_ENABLE_TIMING
    using TimingClock = std::chrono::steady_clock;
    constexpr int TIMING_VALUE_COUNT = 12;
    std::array<double, TIMING_VALUE_COUNT> gLastTimingUs{};
    bool gLastTimingValid = false;
    bool gTimingLogReady = false;
    struct TimingLogRow
    {
        double simulationTimeS = 0.0;
        double demandedGenTorqueNm = 0.0;
        double demandedPitchRad = 0.0;
        int qpSolveMode = 0;
        int nwsrUsed = 0;
        std::array<double, TIMING_VALUE_COUNT> timingUs{};
    };
    std::vector<TimingLogRow> gTimingRows;
    constexpr std::size_t TIMING_FLUSH_ROWS = 256;
#endif

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

    inline std::string replaceExtension(const std::string& path, const std::string& newExt)
    {
        const auto slashPos = path.find_last_of("\\/");
        const auto dotPos = path.find_last_of('.');
        if (dotPos == std::string::npos || (slashPos != std::string::npos && dotPos < slashPos))
            return path + newExt;
        return path.substr(0, dotPos) + newExt;
    }

#ifdef MPC_ENABLE_TIMING
    inline void flushTimingRows()
    {
        if (!gTimingLogReady || gTimingLogPath.empty() || gTimingRows.empty()) return;
        std::ofstream out(gTimingLogPath, std::ios::app);
        if (!out)
        {
            gTimingLogReady = false;
            gTimingRows.clear();
            return;
        }
        out << std::fixed << std::setprecision(6);
        for (const auto& row : gTimingRows)
        {
            out << row.simulationTimeS << ','
                << row.demandedGenTorqueNm << ','
                << row.demandedPitchRad << ','
                << row.qpSolveMode << ','
                << row.nwsrUsed;
            for (double value : row.timingUs) out << ',' << value;
            out << '\n';
        }
        gTimingRows.clear();
    }

    inline bool resetTimingLog(const std::string& path)
    {
        gTimingRows.clear();
        gTimingLogReady = false;
        if (path.empty()) return false;

        std::error_code removeError;
        std::filesystem::remove(std::filesystem::path(path), removeError);
        if (removeError) return false;

        std::ofstream out(path, std::ios::out | std::ios::trunc);
        if (!out) return false;
        out << "simulation_time_s,demanded_gen_torque_Nm,demanded_pitch_rad,qp_solve_mode,nwsr_used,"
               "start_timestamp_us,end_timestamp_us,preprocess_us,model_workspace_us,"
               "prediction_us,cost_us,constraints_us,qpoases_us,extract_u_us,"
               "post_diagnostics_us,total_to_u_us,total_function_us\n";
        out.flush();
        gTimingLogReady = out.good();
        return gTimingLogReady;
    }

    inline void queueTimingRow(float simulationTime, float demandedGenTorque, float demandedPitch)
    {
        if (!gTimingLogReady || !gLastTimingValid) return;
        gTimingRows.push_back({
            static_cast<double>(simulationTime),
            static_cast<double>(demandedGenTorque),
            static_cast<double>(demandedPitch),
            gLastQpSolveMode,
            gLastQpNwsrUsed,
            gLastTimingUs
        });
        if (gTimingRows.size() >= TIMING_FLUSH_ROWS) flushTimingRows();
    }
#endif

    inline bool invert2x2(const MatUU& M, MatUU& Minv)
    {
        const float a = M[0];
        const float b = M[1];
        const float c = M[2];
        const float d = M[3];
        const float det = a * d - b * c;
        if (std::fabs(det) < 1.0e-12f) return false;

        const float invDet = 1.0f / det;
        Minv = { d * invDet, -b * invDet,
                -c * invDet,  a * invDet };
        return true;
    }

    inline MatNN makeDiagonalQ()
    {
        MatNN Q{};
        Q[0 * N_STATE + 0] = gQOmega;
        Q[1 * N_STATE + 1] = gQX;
        Q[2 * N_STATE + 2] = gQV;
        Q[3 * N_STATE + 3] = gQTg;
        Q[4 * N_STATE + 4] = gQBeta;
        return Q;
    }

    inline MatUU makeDiagonalR()
    {
        MatUU R{};
        R[0 * N_CTRL + 0] = gRT;
        R[1 * N_CTRL + 1] = gRB;
        return R;
    }

    inline bool solveRiccati2x2(const MatUU& S, const MatUN& BtPA, MatUN& K)
    {
        MatUU Sinv{};
        if (!invert2x2(S, Sinv)) return false;
        for (int r = 0; r < N_CTRL; ++r)
        {
            for (int c = 0; c < N_STATE; ++c)
            {
                const float val =
                    Sinv[r * N_CTRL + 0] * BtPA[0 * N_STATE + c] +
                    Sinv[r * N_CTRL + 1] * BtPA[1 * N_STATE + c];
                K[r * N_STATE + c] = -val;
            }
        }
        return true;
    }

    inline bool computeDiscreteLqrTerminal(const MatNN& A, const MatNU& B, MatNN& P, MatUN& K)
    {
        const MatNN Q = makeDiagonalQ();
        const MatUU R = makeDiagonalR();
        P = Q;

        for (int iter = 0; iter < 500; ++iter)
        {
            MatUN BtPA{};
            MatUU BtPB{};
            MatNN AtPA{};
            MatNN AtPBK{};
            MatUN Knew{};
            MatNN Pnext{};

            for (int i = 0; i < N_CTRL; ++i)
            {
                for (int j = 0; j < N_STATE; ++j)
                {
                    float acc = 0.0f;
                    for (int k = 0; k < N_STATE; ++k)
                        acc += B[k * N_CTRL + i] * P[k * N_STATE + j];
                    BtPA[i * N_STATE + j] = acc;
                }
            }

            for (int i = 0; i < N_CTRL; ++i)
            {
                for (int j = 0; j < N_CTRL; ++j)
                {
                    float acc = 0.0f;
                    for (int k = 0; k < N_STATE; ++k)
                        acc += BtPA[i * N_STATE + k] * B[k * N_CTRL + j];
                    BtPB[i * N_CTRL + j] = acc + R[i * N_CTRL + j];
                }
            }

            for (int i = 0; i < N_CTRL; ++i)
            {
                for (int j = 0; j < N_STATE; ++j)
                {
                    float acc = 0.0f;
                    for (int k = 0; k < N_STATE; ++k)
                        acc += BtPA[i * N_STATE + k] * A[k * N_STATE + j];
                    BtPA[i * N_STATE + j] = acc;
                }
            }

            if (!solveRiccati2x2(BtPB, BtPA, Knew)) return false;

            for (int i = 0; i < N_STATE; ++i)
            {
                for (int j = 0; j < N_STATE; ++j)
                {
                    float acc = 0.0f;
                    for (int k = 0; k < N_STATE; ++k)
                        acc += A[k * N_STATE + i] * P[k * N_STATE + j];
                    AtPA[i * N_STATE + j] = acc;
                }
            }

            for (int i = 0; i < N_STATE; ++i)
            {
                for (int j = 0; j < N_STATE; ++j)
                {
                    float accAPA = 0.0f;
                    for (int k = 0; k < N_STATE; ++k)
                        accAPA += AtPA[i * N_STATE + k] * A[k * N_STATE + j];

                    Pnext[i * N_STATE + j] = Q[i * N_STATE + j] + accAPA;
                }
            }

            MatNU PB{};
            MatNU APB{};
            for (int i = 0; i < N_STATE; ++i)
            {
                for (int j = 0; j < N_CTRL; ++j)
                {
                    float acc = 0.0f;
                    for (int k = 0; k < N_STATE; ++k)
                        acc += P[i * N_STATE + k] * B[k * N_CTRL + j];
                    PB[i * N_CTRL + j] = acc;
                }
            }
            for (int i = 0; i < N_STATE; ++i)
            {
                for (int j = 0; j < N_CTRL; ++j)
                {
                    float acc = 0.0f;
                    for (int k = 0; k < N_STATE; ++k)
                        acc += A[k * N_STATE + i] * PB[k * N_CTRL + j];
                    APB[i * N_CTRL + j] = acc;
                }
            }
            for (int i = 0; i < N_STATE; ++i)
            {
                for (int j = 0; j < N_STATE; ++j)
                {
                    float correction = 0.0f;
                    for (int m = 0; m < N_CTRL; ++m)
                        correction += APB[i * N_CTRL + m] * (-Knew[m * N_STATE + j]);
                    Pnext[i * N_STATE + j] -= correction;
                }
            }

            float maxDiff = 0.0f;
            for (int i = 0; i < N_STATE * N_STATE; ++i)
                maxDiff = std::max(maxDiff, std::fabs(Pnext[i] - P[i]));

            P = Pnext;
            K = Knew;
            if (maxDiff < 1.0e-4f) return true;
        }
        return false;
    }

    inline bool buildRuntimeTerminalSchedule(std::string& err)
    {
        gTerminal = {};
        if (!gTable.loaded)
        {
            err = "Gain tables must be loaded before building terminal schedule.";
            return false;
        }

        gTerminal.windPtsMs = gTable.windPtsMs;
        gTerminal.speedPtsRpm = gTable.speedPtsRpm;
        const std::size_t nWind = gTerminal.windPtsMs.size();
        const std::size_t nSpeed = gTerminal.speedPtsRpm.size();
        gTerminal.PTable.resize(nWind * nSpeed);
        gTerminal.KTable.resize(nWind * nSpeed);

        for (std::size_t iw = 0; iw < nWind; ++iw)
        {
            for (std::size_t is = 0; is < nSpeed; ++is)
            {
                const float windNow = gTerminal.windPtsMs[iw];
                const float speedNow = gTerminal.speedPtsRpm[is];

                const float gTOmegaNow = gTable.GTomega[iw][is];
                const float gTBetaNow = gTable.GTbeta[iw][is];
                const float gFOmegaNow = gTable.GFomega[iw][is];
                const float gFBetaNow = gTable.GFbeta[iw][is];

                MatNN A{};
                MatNU B{};
                const float dt = gPredictionDt;
                const float a11 = 1.0f + dt * gTOmegaNow / gJEQRuntime;
                const float a23 = dt;
                const float a31 = dt * gFOmegaNow / M_T;
                const float a32 = -dt * K_T / M_T;
                const float a33 = 1.0f - dt * C_T / M_T;
                const float b11 = -dt * N_gear / gJEQRuntime;
                const float b12 = dt * gTBetaNow / gJEQRuntime;
                const float b32 = dt * gFBetaNow / M_T;
                const float tgRetain = actuatorRetention(dt, gTorqueActuatorTau);
                const float betaRetain = actuatorRetention(dt, gPitchActuatorTau);
                const float tgCmdGain = 1.0f - tgRetain;
                const float betaCmdGain = 1.0f - betaRetain;

                A[0 * N_STATE + 0] = a11;  A[0 * N_STATE + 3] = b11 * tgRetain;  A[0 * N_STATE + 4] = b12 * betaRetain;
                A[1 * N_STATE + 1] = 1.0f; A[1 * N_STATE + 2] = a23;
                A[2 * N_STATE + 0] = a31;  A[2 * N_STATE + 1] = a32;  A[2 * N_STATE + 2] = a33; A[2 * N_STATE + 4] = b32 * betaRetain;
                A[3 * N_STATE + 3] = tgRetain;
                A[4 * N_STATE + 4] = betaRetain;

                B[0 * N_CTRL + 0] = b11 * tgCmdGain;  B[0 * N_CTRL + 1] = b12 * betaCmdGain;
                B[1 * N_CTRL + 0] = 0.0f; B[1 * N_CTRL + 1] = 0.0f;
                B[2 * N_CTRL + 0] = 0.0f; B[2 * N_CTRL + 1] = b32 * betaCmdGain;
                B[3 * N_CTRL + 0] = tgCmdGain; B[3 * N_CTRL + 1] = 0.0f;
                B[4 * N_CTRL + 0] = 0.0f; B[4 * N_CTRL + 1] = betaCmdGain;

                MatNN P{};
                MatUN K{};
                if (!computeDiscreteLqrTerminal(A, B, P, K))
                {
                    P = makeDiagonalQ();
                    K = {};
                    gTerminal.fallbackPoints += 1;
                }

                // The Riccati iteration is evaluated in single precision.
                // Symmetrizing avoids numerical asymmetry in the QP Hessian.
                for (int r = 0; r < N_STATE; ++r)
                {
                    for (int c = r + 1; c < N_STATE; ++c)
                    {
                        const float sym = 0.5f * (P[r * N_STATE + c] + P[c * N_STATE + r]);
                        P[r * N_STATE + c] = sym;
                        P[c * N_STATE + r] = sym;
                    }
                }

                const std::size_t idx = iw * nSpeed + is;
                gTerminal.PTable[idx] = P;
                gTerminal.KTable[idx] = K;
            }
        }

        gTerminal.ready = true;
        return true;
    }

    inline int lookupTerminalScheduleIndex(float windNow, float speedNowRpm)
    {
        auto nearestIndex = [](const std::vector<float>& arr, float x) -> int
        {
            int idx = 0;
            float best = std::fabs(arr[0] - x);
            for (int i = 1; i < static_cast<int>(arr.size()); ++i)
            {
                const float err = std::fabs(arr[i] - x);
                if (err < best)
                {
                    best = err;
                    idx = i;
                }
            }
            return idx;
        };

        const int iw = nearestIndex(gTerminal.windPtsMs, windNow);
        const int is = nearestIndex(gTerminal.speedPtsRpm, speedNowRpm);
        return iw * static_cast<int>(gTerminal.speedPtsRpm.size()) + is;
    }

    inline MatNN getScheduledTerminalP(float windNow, float speedNowRpm)
    {
        MatNN P{};
        if (!gTerminal.ready || gTerminal.windPtsMs.empty() || gTerminal.speedPtsRpm.empty())
            return makeDiagonalQ();

        windNow = std::clamp(windNow, gTerminal.windPtsMs.front(), gTerminal.windPtsMs.back());
        speedNowRpm = std::clamp(speedNowRpm, gTerminal.speedPtsRpm.front(), gTerminal.speedPtsRpm.back());

        int iw = 0;
        while (iw + 1 < static_cast<int>(gTerminal.windPtsMs.size()) && gTerminal.windPtsMs[iw + 1] < windNow) ++iw;
        int is = 0;
        while (is + 1 < static_cast<int>(gTerminal.speedPtsRpm.size()) && gTerminal.speedPtsRpm[is + 1] < speedNowRpm) ++is;

        const int iw2 = std::min(iw + 1, static_cast<int>(gTerminal.windPtsMs.size()) - 1);
        const int is2 = std::min(is + 1, static_cast<int>(gTerminal.speedPtsRpm.size()) - 1);
        const float w1 = gTerminal.windPtsMs[iw];
        const float w2 = gTerminal.windPtsMs[iw2];
        const float s1 = gTerminal.speedPtsRpm[is];
        const float s2 = gTerminal.speedPtsRpm[is2];
        const float a = (iw2 == iw || std::fabs(w2 - w1) < 1.0e-8f) ? 0.0f : (windNow - w1) / (w2 - w1);
        const float b = (is2 == is || std::fabs(s2 - s1) < 1.0e-8f) ? 0.0f : (speedNowRpm - s1) / (s2 - s1);

        const std::size_t nSpeed = gTerminal.speedPtsRpm.size();
        const MatNN& p11 = gTerminal.PTable[static_cast<std::size_t>(iw) * nSpeed + static_cast<std::size_t>(is)];
        const MatNN& p21 = gTerminal.PTable[static_cast<std::size_t>(iw2) * nSpeed + static_cast<std::size_t>(is)];
        const MatNN& p12 = gTerminal.PTable[static_cast<std::size_t>(iw) * nSpeed + static_cast<std::size_t>(is2)];
        const MatNN& p22 = gTerminal.PTable[static_cast<std::size_t>(iw2) * nSpeed + static_cast<std::size_t>(is2)];
        for (int i = 0; i < N_STATE * N_STATE; ++i)
        {
            P[i] = (1.0f - a) * (1.0f - b) * p11[i]
                + a * (1.0f - b) * p21[i]
                + (1.0f - a) * b * p12[i]
                + a * b * p22[i];
        }
        return P;
    }

    inline MatUN getScheduledTerminalK(float windNow, float speedNowRpm)
    {
        const int idx = lookupTerminalScheduleIndex(windNow, speedNowRpm);
        return gTerminal.KTable[static_cast<std::size_t>(idx)];
    }

    inline std::array<float, N_STATE> buildAugmentedState(
        float time,
        float rotSpeed,
        float windNow,
        float towerDispFA,
        float towerVelFA,
        float prevGenTorque,
        float prevPitchCmd)
    {
        const float omegaRefNow = getRotorSpeedReference(time);
        const float tgRefNow = getGeneratorTorqueReference(time, rotSpeed, windNow);
        const float betaRefNow = getPitchReference(time);
        return {
            rotSpeed - omegaRefNow,
            towerDispFA,
            towerVelFA,
            prevGenTorque - tgRefNow,
            prevPitchCmd - betaRefNow
        };
    }

    inline bool updateRuntimeFractureMetricsFromFAST(const float* avrSWAP)
    {
        const float jRfFromFAST = avrSWAP[FRACTURE_J_RF_IDX];
        const int statusFromFAST = static_cast<int>(std::lround(avrSWAP[FRACTURE_STATUS_IDX]));

        gFractureStatusFromFAST = statusFromFAST;
        gFractureMassImbalance = avrSWAP[FRACTURE_MASS_IMBALANCE_IDX];
        gFractureLocationRR = avrSWAP[FRACTURE_LOCATION_IDX];
        gFractureSigma = avrSWAP[FRACTURE_SIGMA_IDX];

        const bool hasValidInertia = std::isfinite(jRfFromFAST) && jRfFromFAST > 1.0e5f;
        if (statusFromFAST <= 0 || !hasValidInertia)
            return false;

        if (!gFractureInertiaLocked)
        {
            const float inertiaTol = std::max(1.0f, std::fabs(gJRFRuntime) * 1.0e-5f);
            if (std::fabs(jRfFromFAST - gJRFRuntime) > inertiaTol)
            {
                gJRFRuntime = jRfFromFAST;
                gJEQRuntime = gJRFRuntime + J_generator * N_gear * N_gear;
            }

            if (statusFromFAST >= 2)
            {
                gFractureInertiaLocked = true;
                return true;
            }
        }

        return false;
    }

    inline void applyScheduledTerminalFallback(
        float time,
        float windRef,
        float rotSpeed,
        float towerDispFA,
        float towerVelFA,
        float prevAppliedGenTorque,
        float prevAppliedPitch,
        float prevCommandGenTorque,
        float prevCommandPitch,
        float& demandedGenTorque,
        float& demandedPitchCmd)
    {
        const float rotSpeedRPM = rotSpeed * RPS2RPM;
        const MatUN K = getScheduledTerminalK(windRef, rotSpeedRPM);
        const auto z = buildAugmentedState(time, rotSpeed, windRef, towerDispFA, towerVelFA, prevAppliedGenTorque, prevAppliedPitch);

        float cmdTgErr = 0.0f;
        float cmdBetaErr = 0.0f;
        for (int j = 0; j < N_STATE; ++j)
        {
            cmdTgErr += K[0 * N_STATE + j] * z[j];
            cmdBetaErr += K[1 * N_STATE + j] * z[j];
        }

        const float tgRefNow = getGeneratorTorqueReference(time, rotSpeed, windRef);
        const float betaRefNow = getPitchReference(time);
        const float dTgMax = VS_MAX_TQ_RATE * std::max(gPredictionDt, 1.0e-4f);
        const float dBetaMax = PC_MAX_RAT * std::max(gPredictionDt, 1.0e-4f);
        demandedGenTorque = std::clamp(tgRefNow + cmdTgErr, VS_MIN_TQ, VS_MAX_TQ);
        demandedPitchCmd = std::clamp(betaRefNow + cmdBetaErr, PC_MIN_PIT, PC_MAX_PIT);
        demandedGenTorque = std::clamp(demandedGenTorque, prevCommandGenTorque - dTgMax, prevCommandGenTorque + dTgMax);
        demandedPitchCmd = std::clamp(demandedPitchCmd, prevCommandPitch - dBetaMax, prevCommandPitch + dBetaMax);
        demandedGenTorque = std::clamp(demandedGenTorque, VS_MIN_TQ, VS_MAX_TQ);
        demandedPitchCmd = std::clamp(demandedPitchCmd, PC_MIN_PIT, PC_MAX_PIT);
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

    inline bool loadReferenceCurveTable(const std::string& path)
    {
        std::ifstream in(path);
        if (!in) return false;

        ReferenceCurveTableData table{};
        std::string line;
        bool firstRow = true;
        while (std::getline(in, line))
        {
            line = trim(stripComment(line));
            if (line.empty()) continue;

            if (firstRow)
            {
                firstRow = false;
                continue;
            }

            std::stringstream ss(line);
            std::string item;
            std::vector<std::string> fields;
            while (std::getline(ss, item, ','))
                fields.push_back(trim(item));

            if (fields.size() < 5) return false;
            table.windPtsMs.push_back(std::stof(fields[0]));
            table.tgHiKnm.push_back(std::stof(fields[1]));
            table.omegaHiRpm.push_back(std::stof(fields[2]));
            table.omegaLoRpm.push_back(std::stof(fields[3]));
            table.shapeP.push_back(std::stof(fields[4]));
        }

        if (table.windPtsMs.empty()) return false;

        std::vector<std::size_t> order(table.windPtsMs.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return table.windPtsMs[a] < table.windPtsMs[b];
            });

        ReferenceCurveTableData sorted{};
        sorted.loaded = true;
        for (std::size_t idx : order)
        {
            sorted.windPtsMs.push_back(table.windPtsMs[idx]);
            sorted.tgHiKnm.push_back(table.tgHiKnm[idx]);
            sorted.omegaHiRpm.push_back(table.omegaHiRpm[idx]);
            sorted.omegaLoRpm.push_back(table.omegaLoRpm[idx]);
            sorted.shapeP.push_back(table.shapeP[idx]);
        }

        gRefCurve = std::move(sorted);
        return true;
    }

    inline bool loadShutdownReferenceTrajectory(const std::string& path)
    {
        std::ifstream in(path);
        if (!in) return false;

        ShutdownReferenceTrajectoryData table{};
        std::string line;
        bool firstRow = true;
        while (std::getline(in, line))
        {
            line = trim(stripComment(line));
            if (line.empty()) continue;

            if (firstRow)
            {
                firstRow = false;
                continue;
            }

            std::stringstream ss(line);
            std::string item;
            std::vector<std::string> fields;
            while (std::getline(ss, item, ','))
                fields.push_back(trim(item));

            if (fields.size() < 4) return false;
            table.timeS.push_back(std::stof(fields[0]));
            table.omegaRefRpm.push_back(std::stof(fields[1]));
            table.betaRefDeg.push_back(std::stof(fields[2]));
            table.tgRefKnm.push_back(std::stof(fields[3]));
        }

        if (table.timeS.empty()) return false;

        std::vector<std::size_t> order(table.timeS.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
            return table.timeS[a] < table.timeS[b];
            });

        ShutdownReferenceTrajectoryData sorted{};
        sorted.loaded = true;
        for (std::size_t idx : order)
        {
            sorted.timeS.push_back(table.timeS[idx]);
            sorted.omegaRefRpm.push_back(table.omegaRefRpm[idx]);
            sorted.betaRefDeg.push_back(table.betaRefDeg[idx]);
            sorted.tgRefKnm.push_back(table.tgRefKnm[idx]);
        }

        gShutdownRef = std::move(sorted);
        return true;
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
        auto makePath = [&](const std::string& file) { return dir + "/" + trim(file); };

        if (!loadCsvTable(makePath(lines[4]), nWind, nSpeed, gTable.GTomega)) { err = "Failed loading GTomega table."; return false; }
        if (!loadCsvTable(makePath(lines[5]), nWind, nSpeed, gTable.GTbeta)) { err = "Failed loading GTbeta table.";  return false; }
        if (!loadCsvTable(makePath(lines[6]), nWind, nSpeed, gTable.GTu)) { err = "Failed loading GTu table.";     return false; }
        if (!loadCsvTable(makePath(lines[7]), nWind, nSpeed, gTable.GFomega)) { err = "Failed loading GFomega table."; return false; }
        if (!loadCsvTable(makePath(lines[8]), nWind, nSpeed, gTable.GFbeta)) { err = "Failed loading GFbeta table.";  return false; }
        if (!loadCsvTable(makePath(lines[9]), nWind, nSpeed, gTable.GFu)) { err = "Failed loading GFu table.";     return false; }
        gRefCurve = {};
        if (!loadReferenceCurveTable(makePath("TgRefCurveTable_runtime.csv")))
            loadReferenceCurveTable(makePath("TgRefCurveTable.csv"));

        // The paper-study path keeps fixed post-fracture references.
        // Time-varying shutdown trajectories remain disabled unless they are
        // explicitly re-enabled in code.
        gShutdownRef = {};

        if (lines.size() >= 11) gFractureTime = std::stof(lines[10]);

        // Backward-compatible parsing:
        // old format:  19 data lines, horizons fixed at 5/5 in code
        // mid format:  21 data lines, lines[11:12] are N_PRED/N_CTRL_H
        // new format:  23 data lines, adds Q_TG and Q_BETA before R_T/R_B
        // legacy newest format: 27 data lines, lines 21:23 are unused hard limits,
        // line 24 is MPC_DT, line 25 toggles trace-file generation,
        // line 26 sets the qpOASES nWSR iteration budget, and line 27 scales the terminal cost
        // current format: 25+ data lines, line 12 is TERMINAL_UNLOAD_OMEGA_RPM and the unused hard limits are removed
        const bool hasTerminalUnloadLine =
            (lines.size() >= 25 && std::fabs(std::stof(lines[11]) - std::round(std::stof(lines[11]))) > 1.0e-6f);
        std::size_t idx = 11;
        gNPred = 5;
        gNCtrlH = 5;
        gNWSR = 200;
        gPredictionDt = 0.000125f;
        gTerminalCostScale = 1.0f;
        gActuatorOrder = 0;
        gPitchActuatorTau = 0.25f;
        gTorqueActuatorTau = 0.25f;
        gControlDt = 0.010f;
        gEnableQpHotstart = false;
        gQTg = 1.0e-6f;
        gQBeta = 10.0f;
        gEnableTraceFiles = true;
        gTerminalUnloadTriggerRpm = TERMINAL_UNLOAD_TRIGGER_RPM_DEFAULT;
        if (hasTerminalUnloadLine && lines.size() > idx)
            gTerminalUnloadTriggerRpm = std::stof(lines[idx++]);

        if (lines.size() >= idx + 2)
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

        if (lines.size() >= idx + 7)
        {
            if (lines.size() > idx) gQOmega = std::stof(lines[idx++]);
            if (lines.size() > idx) gQX = std::stof(lines[idx++]);
            if (lines.size() > idx) gQV = std::stof(lines[idx++]);
            if (lines.size() > idx) gQTg = std::stof(lines[idx++]);
            if (lines.size() > idx) gQBeta = std::stof(lines[idx++]);
            if (lines.size() > idx) gRT = std::stof(lines[idx++]);
            if (lines.size() > idx) gRB = std::stof(lines[idx++]);
        }
        else
        {
            if (lines.size() > idx) gQOmega = std::stof(lines[idx++]);
            if (lines.size() > idx) gQX = std::stof(lines[idx++]);
            if (lines.size() > idx) gQV = std::stof(lines[idx++]);
            if (lines.size() > idx) gRT = std::stof(lines[idx++]);
            if (lines.size() > idx) gRB = std::stof(lines[idx++]);
        }

        if (lines.size() > idx)
        {
            gPredictionDt = std::stof(lines[idx++]);
        }

        if (gPredictionDt <= 0.0f)
        {
            err = "MPC_DT must be positive.";
            return false;
        }

        if (lines.size() > idx)
        {
            gEnableTraceFiles = (std::stoi(lines[idx++]) != 0);
        }

        if (lines.size() > idx)
        {
            gNWSR = std::stoi(lines[idx++]);
        }

        if (gNWSR < 1)
        {
            err = "nWSR must be at least 1.";
            return false;
        }

        if (lines.size() > idx)
        {
            gTerminalCostScale = std::stof(lines[idx++]);
        }
        if (gTerminalCostScale < 0.0f)
        {
            err = "TERMINAL_COST_SCALE must be nonnegative.";
            return false;
        }

        if (lines.size() > idx)
        {
            gActuatorOrder = std::stoi(lines[idx++]);
        }
        if (lines.size() > idx)
        {
            gPitchActuatorTau = std::stof(lines[idx++]);
        }
        if (lines.size() > idx)
        {
            gTorqueActuatorTau = std::stof(lines[idx++]);
        }
        if (gActuatorOrder != 0 && gActuatorOrder != 1)
        {
            err = "ACTUATOR_ORDER must be 0 or 1.";
            return false;
        }
        if (gPitchActuatorTau < 0.0f || gTorqueActuatorTau < 0.0f)
        {
            err = "Actuator time constants must be nonnegative.";
            return false;
        }
        if (lines.size() > idx)
        {
            gControlDt = std::stof(lines[idx++]);
        }
        if (gControlDt <= 0.0f)
        {
            err = "CTRL_DT must be positive.";
            return false;
        }
        if (lines.size() > idx)
        {
            gEnableQpHotstart = (std::stoi(lines[idx++]) != 0);
        }

        gTable.loaded = true;
        return true;
    }

    inline float clampToRange(float x, float lo, float hi)
    {
        return std::max(lo, std::min(x, hi));
    }

    inline float updateOperatingWindReference(float currentWind, float dtController)
    {
        if (!gState.operatingWindRefInitialized)
        {
            gState.operatingWindRef = currentWind;
            gState.operatingWindRefInitialized = true;
        }
        else
        {
            const float tau = std::max(OPERATING_WIND_REF_TAU, dtController);
            const float alpha = std::exp(-dtController / tau);
            gState.operatingWindRef = alpha * gState.operatingWindRef + (1.0f - alpha) * currentWind;
        }

        if (gTable.loaded && !gTable.windPtsMs.empty())
            gState.operatingWindRef = clampToRange(gState.operatingWindRef, gTable.windPtsMs.front(), gTable.windPtsMs.back());

        return gState.operatingWindRef;
    }

    inline float applyFirstOrderActuator(float previousApplied, float target, float dtController, float tau)
    {
        if (gActuatorOrder == 0 || tau <= 0.0f)
            return target;

        const float dt = std::max(dtController, 0.0f);
        const float alpha = 1.0f - std::exp(-dt / std::max(tau, 1.0e-6f));
        return previousApplied + alpha * (target - previousApplied);
    }

    inline float applyHardRateLimit(
        float previousApplied,
        float requestedApplied,
        float maxAbsRate,
        float dtController)
    {
        if (dtController <= 0.0f)
            return previousApplied;
        const float maxStep = std::max(maxAbsRate, 0.0f) * dtController;
        return previousApplied + std::clamp(
            requestedApplied - previousApplied,
            -maxStep,
            maxStep
        );
    }

    inline float actuatorRetention(float dtController, float tau)
    {
        if (gActuatorOrder == 0 || tau <= 0.0f)
            return 0.0f;
        const float dt = std::max(dtController, 0.0f);
        return std::exp(-dt / std::max(tau, 1.0e-6f));
    }

    inline int classifyTurbulentWindBand(float windRef)
    {
        if (windRef < TURB_WIND_BAND_0_MAX) return 0;
        if (windRef < TURB_WIND_BAND_1_MAX) return 1;
        if (windRef < TURB_WIND_BAND_2_MAX) return 2;
        if (windRef < TURB_WIND_BAND_3_MAX) return 3;
        return 4;
    }

    inline int updateTurbulentWindBand(float windRef)
    {
        if (gState.turbulentWindBand < 0)
        {
            gState.turbulentWindBand = classifyTurbulentWindBand(windRef);
            return gState.turbulentWindBand;
        }

        switch (gState.turbulentWindBand)
        {
        case 0:
            if (windRef >= TURB_WIND_BAND_0_MAX + TURB_WIND_BAND_HYST) gState.turbulentWindBand = 1;
            break;
        case 1:
            if (windRef < TURB_WIND_BAND_0_MAX - TURB_WIND_BAND_HYST) gState.turbulentWindBand = 0;
            else if (windRef >= TURB_WIND_BAND_1_MAX + TURB_WIND_BAND_HYST) gState.turbulentWindBand = 2;
            break;
        case 2:
            if (windRef < TURB_WIND_BAND_1_MAX - TURB_WIND_BAND_HYST) gState.turbulentWindBand = 1;
            else if (windRef >= TURB_WIND_BAND_2_MAX + TURB_WIND_BAND_HYST) gState.turbulentWindBand = 3;
            break;
        case 3:
            if (windRef < TURB_WIND_BAND_2_MAX - TURB_WIND_BAND_HYST) gState.turbulentWindBand = 2;
            else if (windRef >= TURB_WIND_BAND_3_MAX + TURB_WIND_BAND_HYST) gState.turbulentWindBand = 4;
            break;
        default:
            if (windRef < TURB_WIND_BAND_3_MAX - TURB_WIND_BAND_HYST) gState.turbulentWindBand = 3;
            break;
        }

        return gState.turbulentWindBand;
    }

    inline float getScheduledTerminalUnloadTriggerRpm(float windRef)
    {
        switch (updateTurbulentWindBand(windRef))
        {
        case 0: return TURB_OU_BAND0;
        case 1: return TURB_OU_BAND1;
        case 2: return TURB_OU_BAND2;
        case 3: return TURB_OU_BAND3;
        default: return TURB_OU_BAND4;
        }
    }

    inline void readLidarPreviewData(const float* avrSWAP, LidarPreviewData& lidar)
    {
        lidar = {};
        if (avrSWAP == nullptr) return;

        const int sensorType = static_cast<int>(std::lround(avrSWAP[LIDAR_MSR_START + 0]));
        const int numBeams = static_cast<int>(std::lround(avrSWAP[LIDAR_MSR_START + 1]));
        const int numPulseGates = static_cast<int>(std::lround(avrSWAP[LIDAR_MSR_START + 2]));
        const int nPts = numBeams * numPulseGates;
        if (sensorType == 0 || nPts <= 0) return;

        const int dataStart = LIDAR_MSR_START + 4;
        const int maxPts = (LIDAR_MAX_CHAN - 4) / 4;
        const int usedPts = std::min(nPts, maxPts);

        lidar.available = true;
        lidar.sensorType = sensorType;
        lidar.numBeams = numBeams;
        lidar.numPulseGates = numPulseGates;
        lidar.referenceWind = avrSWAP[LIDAR_MSR_START + 3];
        lidar.measuredSpeeds.resize(static_cast<std::size_t>(usedPts));
        lidar.posX.resize(static_cast<std::size_t>(usedPts));
        lidar.posY.resize(static_cast<std::size_t>(usedPts));
        lidar.posZ.resize(static_cast<std::size_t>(usedPts));

        for (int i = 0; i < usedPts; ++i)
        {
            lidar.measuredSpeeds[static_cast<std::size_t>(i)] = avrSWAP[dataStart + i];
            lidar.posX[static_cast<std::size_t>(i)] = avrSWAP[dataStart + usedPts + i];
            lidar.posY[static_cast<std::size_t>(i)] = avrSWAP[dataStart + 2 * usedPts + i];
            lidar.posZ[static_cast<std::size_t>(i)] = avrSWAP[dataStart + 3 * usedPts + i];
        }
    }

    inline void updateLidarHistory(
        float time,
        float currentWind,
        const LidarPreviewData& lidar)
    {
        while (!gLidarHistory.empty() && (time - gLidarHistory.front().tMeas) > LIDAR_HISTORY_MAX_AGE)
            gLidarHistory.pop_front();

        if (!lidar.available || lidar.measuredSpeeds.empty()) return;

        const float convectionWind = std::max(std::fabs(currentWind), LIDAR_MIN_CONV_WIND);
        for (std::size_t i = 0; i < lidar.measuredSpeeds.size(); ++i)
        {
            LidarHistoryPoint point;
            point.tMeas = time;
            point.convectionWind = convectionWind;
            point.measuredSpeed = lidar.measuredSpeeds[i];
            point.posX = (i < lidar.posX.size()) ? lidar.posX[i] : 0.0f;
            gLidarHistory.push_back(point);
        }
    }

    inline std::vector<float> buildWindPreviewDeltas(
        float currentTime,
        float currentWind,
        int nPred,
        float dtPred,
        bool& usedLidarHistory)
    {
        std::vector<float> dWindPreview(static_cast<std::size_t>(std::max(nPred, 0)), 0.0f);
        usedLidarHistory = false;
        if (nPred <= 0 || dtPred <= 0.0f || gLidarHistory.empty()) return dWindPreview;

        std::vector<float> weightSum(static_cast<std::size_t>(nPred), 0.0f);

        for (const auto& pt : gLidarHistory)
        {
            const float xPos = pt.posX;
            const float distanceAhead = std::max(0.0f, -xPos);
            const float convectionWind = std::max(std::fabs(pt.convectionWind), LIDAR_MIN_CONV_WIND);
            const float tArrive = pt.tMeas + distanceAhead / convectionWind;
            const float futureTime = tArrive - currentTime;
            if (futureTime < 0.0f) continue;

            const int stepOffset = static_cast<int>(std::lround(futureTime / dtPred));
            if (stepOffset < 0) continue;
            const int p = std::min(stepOffset, std::max(nPred - 1, 0));
            const float decay = std::exp(-distanceAhead / LIDAR_EVOLUTION_LENGTH);
            const float delta = decay * (pt.measuredSpeed - currentWind);
            dWindPreview[static_cast<std::size_t>(p)] += delta;
            weightSum[static_cast<std::size_t>(p)] += decay;
        }

        float lastValid = 0.0f;
        bool anyAssigned = false;
        for (int p = 0; p < nPred; ++p)
        {
            if (weightSum[static_cast<std::size_t>(p)] > 1.0e-8f)
            {
                dWindPreview[static_cast<std::size_t>(p)] /= weightSum[static_cast<std::size_t>(p)];
                lastValid = dWindPreview[static_cast<std::size_t>(p)];
                anyAssigned = true;
            }
            else
            {
                dWindPreview[static_cast<std::size_t>(p)] = lastValid;
            }
        }
        usedLidarHistory = anyAssigned;
        return dWindPreview;
    }

    inline float getLidarMeanWind(const LidarPreviewData& lidar, float fallbackWind)
    {
        if (!lidar.available || lidar.measuredSpeeds.empty()) return fallbackWind;

        float sumWind = 0.0f;
        for (float v : lidar.measuredSpeeds) sumWind += v;
        return sumWind / static_cast<float>(lidar.measuredSpeeds.size());
    }

    inline void appendDebugLog(const std::string& path, const std::string& line)
    {
        if (!gEnableTraceFiles) return;
        if (path.empty()) return;
        std::ofstream out(path, std::ios::app);
        if (!out) return;
        out << line << '\n';
    }

    inline void resetDebugLog(const std::string& path)
    {
        if (!gEnableTraceFiles) return;
        if (path.empty()) return;
        std::ofstream out(path, std::ios::trunc);
        if (!out) return;
    }

    inline std::string makeTimestampString()
    {
        std::time_t now = std::time(nullptr);
        std::tm tmNow{};
#ifdef _WIN32
        localtime_s(&tmNow, &now);
#else
        tmNow = *std::localtime(&now);
#endif
        std::ostringstream oss;
        oss << (tmNow.tm_year + 1900) << "-"
            << std::setw(2) << std::setfill('0') << (tmNow.tm_mon + 1) << "-"
            << std::setw(2) << std::setfill('0') << tmNow.tm_mday << " "
            << std::setw(2) << std::setfill('0') << tmNow.tm_hour << ":"
            << std::setw(2) << std::setfill('0') << tmNow.tm_min << ":"
            << std::setw(2) << std::setfill('0') << tmNow.tm_sec;
        return oss.str();
    }

    inline void writeDebugLogHeader(const std::string& path)
    {
        if (path.empty()) return;
        resetDebugLog(path);
        appendDebugLog(path, "# DISCON_MPC_CPP lidar debug log");
        appendDebugLog(path, "# generated_at=" + makeTimestampString());

        std::ostringstream cfg;
        cfg << "# N_PRED=" << gNPred
            << ", N_CTRL_H=" << gNCtrlH
            << ", MPC_DT=" << gPredictionDt
            << ", nWSR=" << gNWSR
            << ", terminal_cost_scale=" << gTerminalCostScale
            << ", Q=[" << gQOmega << "," << gQX << "," << gQV << "," << gQTg << "," << gQBeta << "]"
            << ", R=[" << gRT << "," << gRB << "]"
            << ", terminalUnloadOmegaRpm=" << gTerminalUnloadTriggerRpm;
        appendDebugLog(path, cfg.str());

        std::ostringstream sched;
        sched << "# terminal_schedule_points=" << gTerminal.PTable.size()
              << ", wind_grid=" << gTerminal.windPtsMs.size()
              << ", speed_grid=" << gTerminal.speedPtsRpm.size()
              << ", runtime_generated=" << (gTerminal.ready ? 1 : 0)
              << ", fallback_points=" << gTerminal.fallbackPoints;
        appendDebugLog(path, sched.str());
    }

    inline void writePredictionLogHeader(const std::string& path)
    {
        if (!gEnableTraceFiles) return;
        if (!gEnablePredictionTrace) return;
        if (path.empty()) return;
        resetDebugLog(path);
        std::ofstream out(path, std::ios::app);
        if (!out) return;
        out << "time,horizon_step,pred_time,dOmega_pred,x_t_pred,v_t_pred,"
               "dOmega_meas,x_t_meas,v_t_meas,wind_meas\n";
    }

    inline void appendPredictionLogRow(
        const std::string& path,
        float time,
        int horizonStep,
        float predTime,
        float dOmegaPred,
        float xPred,
        float vPred,
        float dOmegaMeas,
        float xMeas,
        float vMeas,
        float windMeas)
    {
        if (!gEnableTraceFiles) return;
        if (!gEnablePredictionTrace) return;
        if (path.empty()) return;
        std::ofstream out(path, std::ios::app);
        if (!out) return;
        out << std::fixed << std::setprecision(6)
            << time << ','
            << horizonStep << ','
            << predTime << ','
            << dOmegaPred << ','
            << xPred << ','
            << vPred << ','
            << dOmegaMeas << ','
            << xMeas << ','
            << vMeas << ','
            << windMeas << '\n';
    }

    inline void writeCommandTraceHeader(const std::string& path)
    {
        if (!gEnableTraceFiles) return;
        if (path.empty()) return;
        resetDebugLog(path);
        std::ofstream out(path, std::ios::app);
        if (!out) return;
        out << "time_s,iStatus,fracture_active,qp_solved,fallback_used,fallback_count,terminal_idx,"
               "rotSpeed_rpm,genSpeed_rpm,towerDispFA_m,towerVelFA_mps,horWindV_mps,operatingWindRef_mps,"
               "tg_ref_Nm,beta_ref_rad,beta_ref_deg,prevGenTorque_Nm,prevPitchCmd_rad,prevPitchCmd_deg,"
               "demandedGenTorque_Nm,demandedPitch_rad,demandedPitch_deg,demandedPitchRate_radps,demandedPitchRate_degps,"
               "appliedGenTorque_Nm,appliedPitch_rad,appliedPitch_deg,appliedPitchRate_radps,appliedPitchRate_degps,"
               "targetMinusAppliedGenTorque_Nm,targetMinusAppliedPitch_rad,targetMinusAppliedPitch_deg,"
               "deltaTg_cmd_Nm,deltaBeta_cmd_rad,has_pred_1step,pred_issued_time_s,pred_time_s,actual_time_s,"
               "pred_dOmega_radps,actual_dOmega_radps,err_dOmega_radps,"
               "pred_rotSpeed_rpm,actual_rotSpeed_rpm,err_rotSpeed_rpm,"
               "pred_towerDispFA_m,actual_towerDispFA_m,err_towerDispFA_m,"
               "pred_towerVelFA_mps,actual_towerVelFA_mps,err_towerVelFA_mps,pred_dt_s,actual_dt_s\n";
    }

    inline void appendCommandTraceRow(
        const std::string& path,
        float time,
        int iStatus,
        bool fractureActive,
        bool qpSolved,
        bool fallbackUsed,
        int fallbackCount,
        int terminalIdx,
        float rotSpeedRpm,
        float genSpeedRpm,
        float towerDispFA,
        float towerVelFA,
        float horWindV,
        float operatingWindRef,
        float tgRefNow,
        float betaRefNow,
        float prevGenTorque,
        float prevPitchCmd,
        float demandedGenTorque,
        float demandedPitch,
        float demandedPitchRate,
        float appliedGenTorque,
        float appliedPitch,
        float appliedPitchRate,
        bool hasPred1Step,
        float predIssuedTime,
        float predTime,
        float predDOmega,
        float actualDOmega,
        float predRotSpeedRpm,
        float actualRotSpeedRpm,
        float predTowerDispFA,
        float actualTowerDispFA,
        float predTowerVelFA,
        float actualTowerVelFA,
        float predDt,
        float actualDt)
    {
        if (!gEnableTraceFiles) return;
        if (path.empty()) return;
        const bool periodicDue =
            (gState.lastCommandTraceTime < 0.0f) ||
            ((time - gState.lastCommandTraceTime) >= (COMMAND_TRACE_DT - 1.0e-6f));
        if (!periodicDue && !fallbackUsed && !hasPred1Step)
            return;

        const float nan = std::numeric_limits<float>::quiet_NaN();
        if (!hasPred1Step)
        {
            predIssuedTime = nan;
            predTime = nan;
            predDOmega = nan;
            actualDOmega = nan;
            predRotSpeedRpm = nan;
            actualRotSpeedRpm = nan;
            predTowerDispFA = nan;
            actualTowerDispFA = nan;
            predTowerVelFA = nan;
            actualTowerVelFA = nan;
            predDt = nan;
            actualDt = nan;
        }

        std::ofstream out(path, std::ios::app);
        if (!out) return;
        gState.lastCommandTraceTime = time;
        out << std::fixed << std::setprecision(6)
            << time << ','
            << iStatus << ','
            << (fractureActive ? 1 : 0) << ','
            << (qpSolved ? 1 : 0) << ','
            << (fallbackUsed ? 1 : 0) << ','
            << fallbackCount << ','
            << terminalIdx << ','
            << rotSpeedRpm << ','
            << genSpeedRpm << ','
            << towerDispFA << ','
            << towerVelFA << ','
            << horWindV << ','
            << operatingWindRef << ','
            << tgRefNow << ','
            << betaRefNow << ','
            << betaRefNow * R2D << ','
            << prevGenTorque << ','
            << prevPitchCmd << ','
            << prevPitchCmd * R2D << ','
            << demandedGenTorque << ','
            << demandedPitch << ','
            << demandedPitch * R2D << ','
            << demandedPitchRate << ','
            << demandedPitchRate * R2D << ','
            << appliedGenTorque << ','
            << appliedPitch << ','
            << appliedPitch * R2D << ','
            << appliedPitchRate << ','
            << appliedPitchRate * R2D << ','
            << demandedGenTorque - appliedGenTorque << ','
            << demandedPitch - appliedPitch << ','
            << (demandedPitch - appliedPitch) * R2D << ','
            << demandedGenTorque - prevGenTorque << ','
            << demandedPitch - prevPitchCmd << ','
            << (hasPred1Step ? 1 : 0) << ','
            << predIssuedTime << ','
            << predTime << ','
            << time << ','
            << predDOmega << ','
            << actualDOmega << ','
            << predDOmega - actualDOmega << ','
            << predRotSpeedRpm << ','
            << actualRotSpeedRpm << ','
            << predRotSpeedRpm - actualRotSpeedRpm << ','
            << predTowerDispFA << ','
            << actualTowerDispFA << ','
            << predTowerDispFA - actualTowerDispFA << ','
            << predTowerVelFA << ','
            << actualTowerVelFA << ','
            << predTowerVelFA - actualTowerVelFA << ','
            << predDt << ','
            << actualDt << '\n';
    }

    inline void writeInputTraceHeader(const std::string& path)
    {
        if (!gEnableTraceFiles) return;
        if (!gEnableInputTrace) return;
        if (path.empty()) return;
        resetDebugLog(path);
        std::ofstream out(path, std::ios::app);
        if (!out) return;
        out << "time,iStatus,fracture_active,rotSpeed,horWindV,"
               "avr_1019_towerVelFA,avr_1020_towerVelSS,avr_1021_towerDispFA,"
               "avr_1023_J_RF,avr_1024_mass_imbalance,avr_1025_frac_r_R,avr_1026_sigma,avr_1027_status,"
               "runtime_J_RF,runtime_J_EQ,"
               "out_rotSpeed_rads,out_towerVelFA,out_towerDispFA\n";
    }

    inline void appendInputTraceRow(
        const std::string& path,
        float time,
        int iStatus,
        bool fractureActive,
        float rotSpeed,
        float horWindV,
        float avrTowerVelFA,
        float avrTowerVelSS,
        float avrTowerDispFA,
        float avrJRF,
        float avrMassImbalance,
        float avrFractureLocation,
        float avrFractureSigma,
        int avrFractureStatus,
        float towerVelFA,
        float towerDispFA)
    {
        if (!gEnableTraceFiles) return;
        if (!gEnableInputTrace) return;
        if (path.empty()) return;
        std::ofstream out(path, std::ios::app);
        if (!out) return;
        out << std::fixed << std::setprecision(6)
            << time << ','
            << iStatus << ','
            << (fractureActive ? 1 : 0) << ','
            << rotSpeed << ','
            << horWindV << ','
            << avrTowerVelFA << ','
            << avrTowerVelSS << ','
            << avrTowerDispFA << ','
            << avrJRF << ','
            << avrMassImbalance << ','
            << avrFractureLocation << ','
            << avrFractureSigma << ','
            << avrFractureStatus << ','
            << gJRFRuntime << ','
            << gJEQRuntime << ','
            << rotSpeed << ','
            << towerVelFA << ','
            << towerDispFA << '\n';
    }

    inline void writeTuningTraceHeader(const std::string& path)
    {
        if (!gEnableTraceFiles) return;
        if (!gEnableTuningTrace) return;
        if (path.empty()) return;
        resetDebugLog(path);
        std::ofstream out(path, std::ios::app);
        if (!out) return;
        out << "time,fracture_active,qp_solved,fallback_used,fallback_count,terminal_idx,"
               "rotSpeed_rpm,rotSpeed_rads,towerDispFA,towerVelFA,horWindV,dWind,"
               "tg_ref,beta_ref,prevGenTorque,prevPitchCmd,demandedGenTorque,demandedPitch,"
               "deltaTg_cmd,deltaBeta_cmd,"
               "x0_dOmega,x1_xt,x2_vt,x3_tgerr,x4_betaerr,"
               "J_total,J_omega,J_x,J_v,J_tg,J_beta,J_du_tg,J_du_beta,J_terminal,"
               "Qomega_eff,Qx_eff,Qv_eff,QTg_eff,QBeta_eff,RT_eff,RB_eff,"
               "pred_max_abs_dOmega,pred_max_abs_xt,pred_max_abs_vt,"
               "trace_note\n";
    }

    inline void appendTuningTraceRow(
        const std::string& path,
        float time,
        bool fractureActive,
        bool qpSolved,
        bool fallbackUsed,
        int fallbackCount,
        int terminalIdx,
        float rotSpeedRpm,
        float rotSpeedRads,
        float towerDispFA,
        float towerVelFA,
        float horWindV,
        float dWind,
        float tgRefNow,
        float betaRefNow,
        float prevGenTorque,
        float prevPitchCmd,
        float demandedGenTorque,
        float demandedPitch,
        float deltaTgCmd,
        float deltaBetaCmd,
        const std::array<float, N_STATE>& xbar0,
        float Jtotal,
        float Jomega,
        float Jx,
        float Jv,
        float Jtg,
        float Jbeta,
        float JduTg,
        float JduBeta,
        float Jterminal,
        float qOmegaEff,
        float qXEff,
        float qVEff,
        float qTgEff,
        float qBetaEff,
        float rTEff,
        float rBEff,
        float predMaxAbsDOmega,
        float predMaxAbsXt,
        float predMaxAbsVt,
        const std::string& note)
    {
        if (!gEnableTraceFiles) return;
        if (!gEnableTuningTrace) return;
        if (path.empty()) return;
        const bool fallbackTransition = (fallbackUsed != gState.lastTuningTraceFallbackState);
        const bool noteIsSpecial = (note != "qp");
        const bool periodicDue =
            (gState.lastTuningTraceTime < 0.0f) ||
            ((time - gState.lastTuningTraceTime) >= TUNING_TRACE_DT);
        if (!(fallbackUsed || fallbackTransition || noteIsSpecial || periodicDue))
            return;

        std::ofstream out(path, std::ios::app);
        if (!out) return;
        gState.lastTuningTraceTime = time;
        gState.lastTuningTraceFallbackState = fallbackUsed;
        out << std::fixed << std::setprecision(6)
            << time << ','
            << (fractureActive ? 1 : 0) << ','
            << (qpSolved ? 1 : 0) << ','
            << (fallbackUsed ? 1 : 0) << ','
            << fallbackCount << ','
            << terminalIdx << ','
            << rotSpeedRpm << ','
            << rotSpeedRads << ','
            << towerDispFA << ','
            << towerVelFA << ','
            << horWindV << ','
            << dWind << ','
            << tgRefNow << ','
            << betaRefNow << ','
            << prevGenTorque << ','
            << prevPitchCmd << ','
            << demandedGenTorque << ','
            << demandedPitch << ','
            << deltaTgCmd << ','
            << deltaBetaCmd << ','
            << xbar0[0] << ','
            << xbar0[1] << ','
            << xbar0[2] << ','
            << xbar0[3] << ','
            << xbar0[4] << ','
            << Jtotal << ','
            << Jomega << ','
            << Jx << ','
            << Jv << ','
            << Jtg << ','
            << Jbeta << ','
            << JduTg << ','
            << JduBeta << ','
            << Jterminal << ','
            << qOmegaEff << ','
            << qXEff << ','
            << qVEff << ','
            << qTgEff << ','
            << qBetaEff << ','
            << rTEff << ','
            << rBEff << ','
            << predMaxAbsDOmega << ','
            << predMaxAbsXt << ','
            << predMaxAbsVt << ','
            << note << '\n';
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

        const float g11 = table[iw][is];
        const float g21 = table[iw2][is];
        const float g12 = table[iw][is2];
        const float g22 = table[iw2][is2];

        return (1.0f - a) * (1.0f - b) * g11
            + a * (1.0f - b) * g21
            + (1.0f - a) * b * g12
            + a * b * g22;
    }

    inline const char* qpReturnValueToString(returnValue rv)
    {
        switch (rv)
        {
        case SUCCESSFUL_RETURN:                  return "SUCCESSFUL_RETURN";
        case RET_INVALID_ARGUMENTS:             return "RET_INVALID_ARGUMENTS";
        case RET_INIT_FAILED:                   return "RET_INIT_FAILED";
        case RET_INIT_FAILED_TQ:                return "RET_INIT_FAILED_TQ";
        case RET_INIT_FAILED_CHOLESKY:          return "RET_INIT_FAILED_CHOLESKY";
        case RET_INIT_FAILED_HOTSTART:          return "RET_INIT_FAILED_HOTSTART";
        case RET_INIT_FAILED_INFEASIBILITY:     return "RET_INIT_FAILED_INFEASIBILITY";
        case RET_INIT_FAILED_UNBOUNDEDNESS:     return "RET_INIT_FAILED_UNBOUNDEDNESS";
        case RET_QP_UNBOUNDED:                  return "RET_QP_UNBOUNDED";
        case RET_QP_INFEASIBLE:                 return "RET_QP_INFEASIBLE";
        case RET_QP_NOT_SOLVED:                 return "RET_QP_NOT_SOLVED";
        case RET_UNABLE_TO_SOLVE_QP:            return "RET_UNABLE_TO_SOLVE_QP";
        case RET_HOTSTART_FAILED:               return "RET_HOTSTART_FAILED";
        case RET_STEPDIRECTION_DETERMINATION_FAILED:
            return "RET_STEPDIRECTION_DETERMINATION_FAILED";
        case RET_STEPLENGTH_DETERMINATION_FAILED:
            return "RET_STEPLENGTH_DETERMINATION_FAILED";
        case RET_HOMOTOPY_STEP_FAILED:          return "RET_HOMOTOPY_STEP_FAILED";
        case RET_HOTSTART_STOPPED_INFEASIBILITY:
            return "RET_HOTSTART_STOPPED_INFEASIBILITY";
        case RET_HOTSTART_STOPPED_UNBOUNDEDNESS:
            return "RET_HOTSTART_STOPPED_UNBOUNDEDNESS";
        case RET_MAX_NWSR_REACHED:              return "RET_MAX_NWSR_REACHED";
        case RET_ADDCONSTRAINT_FAILED_INFEASIBILITY:
            return "RET_ADDCONSTRAINT_FAILED_INFEASIBILITY";
        case RET_ADDBOUND_FAILED_INFEASIBILITY: return "RET_ADDBOUND_FAILED_INFEASIBILITY";
        case RET_ENSURELI_FAILED:               return "RET_ENSURELI_FAILED";
        case RET_HESSIAN_NOT_SPD:               return "RET_HESSIAN_NOT_SPD";
        case RET_HESSIAN_INDEFINITE:            return "RET_HESSIAN_INDEFINITE";
        case RET_MATRIX_FACTORISATION_FAILED:   return "RET_MATRIX_FACTORISATION_FAILED";
        case RET_USING_REGULARISATION:          return "RET_USING_REGULARISATION";
        default:                                return "RET_UNKNOWN_OR_UNMAPPED";
        }
    }

    inline const char* qpSimpleStatusToString(returnValue rv)
    {
        switch (getSimpleStatus(rv, BT_FALSE))
        {
        case 0:  return "solved";
        case 1:  return "iteration_limit";
        case -1: return "internal_error";
        case -2: return "infeasible";
        case -3: return "unbounded";
        default: return "unknown";
        }
    }

    inline float interp1dClamped(const std::vector<float>& xPts, const std::vector<float>& yPts, float x)
    {
        if (xPts.empty() || yPts.empty()) return 0.0f;
        if (xPts.size() == 1 || yPts.size() == 1) return yPts.front();

        x = clampToRange(x, xPts.front(), xPts.back());
        int idx = 0;
        while (idx + 1 < static_cast<int>(xPts.size()) && xPts[idx + 1] < x) ++idx;
        const int idx2 = std::min(idx + 1, static_cast<int>(xPts.size()) - 1);
        const float x1 = xPts[idx];
        const float x2 = xPts[idx2];
        const float y1 = yPts[idx];
        const float y2 = yPts[idx2];
        if (idx2 == idx || std::fabs(x2 - x1) < 1.0e-8f) return y1;
        const float a = (x - x1) / (x2 - x1);
        return (1.0f - a) * y1 + a * y2;
    }

    inline float evaluateReferenceCurveKnm(float rotSpeedRpm, float windNow)
    {
        if (!gRefCurve.loaded) return TG_REF * 1.0e-3f;

        const float tgHi = interp1dClamped(gRefCurve.windPtsMs, gRefCurve.tgHiKnm, windNow);
        const float omegaHi = interp1dClamped(gRefCurve.windPtsMs, gRefCurve.omegaHiRpm, windNow);
        const float omegaLo = interp1dClamped(gRefCurve.windPtsMs, gRefCurve.omegaLoRpm, windNow);
        const float shapeP = std::max(0.1f, interp1dClamped(gRefCurve.windPtsMs, gRefCurve.shapeP, windNow));

        if (rotSpeedRpm >= omegaHi) return tgHi;
        if (rotSpeedRpm <= omegaLo) return 0.0f;

        const float xi = (rotSpeedRpm - omegaLo) / std::max(omegaHi - omegaLo, 1.0e-6f);
        return tgHi * std::pow(clampToRange(xi, 0.0f, 1.0f), shapeP);
    }

    inline float getShutdownRefTime(float time)
    {
        return std::max(0.0f, time - gFractureTime);
    }

    inline float getRotorSpeedReference(float time)
    {
        if (!gShutdownRef.loaded) return OMEGA_REF;
        const float tau = getShutdownRefTime(time);
        return interp1dClamped(gShutdownRef.timeS, gShutdownRef.omegaRefRpm, tau) / RPS2RPM;
    }

    inline float getPitchReference(float time)
    {
        if (!gShutdownRef.loaded) return BETA_REF;
        const float tau = getShutdownRefTime(time);
        return interp1dClamped(gShutdownRef.timeS, gShutdownRef.betaRefDeg, tau) * D2R;
    }

    inline float getGeneratorTorqueReference(float time, float rotSpeed, float windNow)
    {
        if (gShutdownRef.loaded)
        {
            const float tau = getShutdownRefTime(time);
            return interp1dClamped(gShutdownRef.timeS, gShutdownRef.tgRefKnm, tau) * 1000.0f;
        }
        // Current paper-aligned path uses the fixed generator-torque reference.
        (void)time;
        (void)rotSpeed;
        (void)windNow;
        return TG_REF;
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
        float time,
        float dtController,
        float rotSpeed,
        float horWindV,
        float towerDispFA,
        float towerVelFA,
        const LidarPreviewData& lidar,
        float prevAppliedGenTorque,
        float prevAppliedPitch,
        float prevCommandGenTorque,
        float prevCommandPitch,
        std::string& solveErr,
        float& demandedGenTorque,
        float& demandedPitchCmd)
    {
#ifdef MPC_ENABLE_TIMING
        gLastTimingUs.fill(0.0);
        gLastTimingValid = false;
        const auto timingStart = TimingClock::now();
        auto timingLast = timingStart;
        auto markTiming = [&](int index) {
            const auto now = TimingClock::now();
            gLastTimingUs[static_cast<std::size_t>(index)] =
                std::chrono::duration<double, std::micro>(now - timingLast).count();
            timingLast = now;
            return now;
        };
        gLastTimingUs[0] =
            std::chrono::duration<double, std::micro>(timingStart.time_since_epoch()).count();
#endif
        const int nPred = gNPred;
        const int nCtrlH = gNCtrlH;
        const float dtPred = gPredictionDt;
        const float rotSpeedRPM = rotSpeed * 9.5492966f;
        const float windRef = gState.operatingWindRef;
        gCurrentTerminalIndex = lookupTerminalScheduleIndex(windRef, rotSpeedRPM);
        const float omegaRefNow = getRotorSpeedReference(time);
        const float tgRefNow = getGeneratorTorqueReference(time, rotSpeed, windRef);
        const float betaRefNow = getPitchReference(time);
        const float dWind = horWindV - windRef;
        bool usedLidarHistory = false;
        const std::vector<float> dWindPreview = buildWindPreviewDeltas(time, horWindV, nPred, dtPred, usedLidarHistory);

        const float gTOmegaNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GTomega, windRef, rotSpeedRPM) : G_TOMEGA_DEFAULT;
        const float gTBetaNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GTbeta, windRef, rotSpeedRPM) : G_TBETA_DEFAULT;
        const float gTUNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GTu, windRef, rotSpeedRPM) : G_TU_DEFAULT;
        const float gFOmegaNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GFomega, windRef, rotSpeedRPM) : G_FOMEGA_DEFAULT;
        const float gFBetaNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GFbeta, windRef, rotSpeedRPM) : G_FBETA_DEFAULT;
        const float gFUNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GFu, windRef, rotSpeedRPM) : G_FU_DEFAULT;
#ifdef MPC_ENABLE_TIMING
        markTiming(2); // references, preview construction, and gain-table interpolation
#endif

        // Augmented state with first-order actuator memory:
        // xbar = [ dOmega, x_t, v_t, Tg_applied - Tg_ref, beta_applied - beta_ref ]^T.
        // The QP decision variables are absolute command errors
        // [Tg_cmd - Tg_ref, beta_cmd - beta_ref], not command increments.
        const float x0 = rotSpeed - omegaRefNow;
        const float x1 = towerDispFA;
        const float x2 = towerVelFA;
        const float x3 = prevAppliedGenTorque - tgRefNow;
        const float x4 = prevAppliedPitch - betaRefNow;
        const float prevCommandTgErr = prevCommandGenTorque - tgRefNow;
        const float prevCommandBetaErr = prevCommandPitch - betaRefNow;

        std::array<float, N_STATE> xbar0 = { x0, x1, x2, x3, x4 };

        // One-step Euler-discretized augmented model
        MatNN Abar{};
        MatNU Bbar{};
        std::array<float, N_STATE> Ebar{};
        auto Aat = [&](int r, int c) -> float& { return Abar[r * N_STATE + c]; };
        auto Bat = [&](int r, int c) -> float& { return Bbar[r * N_CTRL + c]; };

        const float a11 = 1.0f + dtPred * gTOmegaNow / gJEQRuntime;
        const float a22 = 1.0f;
        const float a23 = dtPred;
        const float a31 = dtPred * gFOmegaNow / M_T;
        const float a32 = -dtPred * K_T / M_T;
        const float a33 = 1.0f - dtPred * C_T / M_T;
        const float b11 = -dtPred * N_gear / gJEQRuntime;
        const float b12 = dtPred * gTBetaNow / gJEQRuntime;
        const float b32 = dtPred * gFBetaNow / M_T;
        const float e11 = dtPred * gTUNow / gJEQRuntime;
        const float e31 = dtPred * gFUNow / M_T;
        const float tgRetain = actuatorRetention(dtPred, gTorqueActuatorTau);
        const float betaRetain = actuatorRetention(dtPred, gPitchActuatorTau);
        const float tgCmdGain = 1.0f - tgRetain;
        const float betaCmdGain = 1.0f - betaRetain;

        Aat(0, 0) = a11;  Aat(0, 3) = b11 * tgRetain;  Aat(0, 4) = b12 * betaRetain;
        Aat(1, 1) = a22;  Aat(1, 2) = a23;
        Aat(2, 0) = a31;  Aat(2, 1) = a32;  Aat(2, 2) = a33;  Aat(2, 4) = b32 * betaRetain;
        Aat(3, 3) = tgRetain;
        Aat(4, 4) = betaRetain;

        Bat(0, 0) = b11 * tgCmdGain;  Bat(0, 1) = b12 * betaCmdGain;
        Bat(1, 0) = 0.0f; Bat(1, 1) = 0.0f;
        Bat(2, 0) = 0.0f; Bat(2, 1) = b32 * betaCmdGain;
        Bat(3, 0) = tgCmdGain; Bat(3, 1) = 0.0f;
        Bat(4, 0) = 0.0f; Bat(4, 1) = betaCmdGain;

        Ebar[0] = e11;
        Ebar[1] = 0.0f;
        Ebar[2] = e31;
        Ebar[3] = 0.0f;
        Ebar[4] = 0.0f;

        MatNN terminalP = getScheduledTerminalP(horWindV, rotSpeedRPM);
        for (float& value : terminalP) value *= gTerminalCostScale;
        const MatUN terminalK = getScheduledTerminalK(horWindV, rotSpeedRPM);
        (void)terminalK;

        // Build prediction matrices X = F*x0 + G*U
        const int NX = N_STATE * nPred;
        const int NU = N_CTRL * nCtrlH;
        const int NC_INPUT = NU;
        const int NC_STATE = nPred;
        const int NC = NC_INPUT + NC_STATE;
        if (gSolverWs.nx != NX || gSolverWs.nu != NU || gSolverWs.nc != NC)
        {
            gSolverWs.nx = NX;
            gSolverWs.nu = NU;
            gSolverWs.nc = NC;
            gSolverWs.G.assign(static_cast<std::size_t>(NX * NU), 0.0f);
            gSolverWs.c.assign(static_cast<std::size_t>(NX), 0.0f);
            gSolverWs.H.assign(static_cast<std::size_t>(NU * NU), 0.0);
            gSolverWs.g.assign(static_cast<std::size_t>(NU), 0.0);
            gSolverWs.lb.assign(static_cast<std::size_t>(NU), 0.0);
            gSolverWs.ub.assign(static_cast<std::size_t>(NU), 0.0);
            gSolverWs.Acon.assign(static_cast<std::size_t>(NC * NU), 0.0);
            gSolverWs.lbA.assign(static_cast<std::size_t>(NC), BIG_NEG);
            gSolverWs.ubA.assign(static_cast<std::size_t>(NC), BIG_POS);
            gSolverWs.xOpt.assign(static_cast<std::size_t>(NU), 0.0);
            gSolverWs.WG.assign(static_cast<std::size_t>(N_STATE * NU), 0.0f);
            gSolverWs.powers.resize(static_cast<std::size_t>(nPred));
        }
        else
        {
            std::fill(gSolverWs.G.begin(), gSolverWs.G.end(), 0.0f);
            std::fill(gSolverWs.c.begin(), gSolverWs.c.end(), 0.0f);
            std::fill(gSolverWs.H.begin(), gSolverWs.H.end(), 0.0);
            std::fill(gSolverWs.g.begin(), gSolverWs.g.end(), 0.0);
            std::fill(gSolverWs.lb.begin(), gSolverWs.lb.end(), 0.0);
            std::fill(gSolverWs.ub.begin(), gSolverWs.ub.end(), 0.0);
            std::fill(gSolverWs.Acon.begin(), gSolverWs.Acon.end(), 0.0);
            std::fill(gSolverWs.lbA.begin(), gSolverWs.lbA.end(), BIG_NEG);
            std::fill(gSolverWs.ubA.begin(), gSolverWs.ubA.end(), BIG_POS);
        }
#ifdef MPC_ENABLE_TIMING
        markTiming(3); // local model, terminal lookup, workspace allocation/reset
#endif

        auto& G = gSolverWs.G;
        auto& c = gSolverWs.c;
        auto& H = gSolverWs.H;
        auto& g = gSolverWs.g;
        auto& lb = gSolverWs.lb;
        auto& ub = gSolverWs.ub;
        auto& Acon = gSolverWs.Acon;
        auto& lbA = gSolverWs.lbA;
        auto& ubA = gSolverWs.ubA;
        auto& xOpt = gSolverWs.xOpt;
        auto& WG = gSolverWs.WG;

        auto matMulSq = [&](const std::array<float, N_STATE* N_STATE>& M,
            const std::array<float, N_STATE* N_STATE>& N) {
                std::array<float, N_STATE* N_STATE> R{};
                for (int i = 0; i < N_STATE; ++i)
                    for (int j = 0; j < N_STATE; ++j)
                        for (int k = 0; k < N_STATE; ++k)
                            R[i * N_STATE + j] += M[i * N_STATE + k] * N[k * N_STATE + j];
                return R;
            };

        auto matMulAB = [&](const std::array<float, N_STATE* N_STATE>& M,
            const std::array<float, N_STATE* N_CTRL>& N) {
                std::array<float, N_STATE* N_CTRL> R{};
                for (int i = 0; i < N_STATE; ++i)
                    for (int j = 0; j < N_CTRL; ++j)
                        for (int k = 0; k < N_STATE; ++k)
                            R[i * N_CTRL + j] += M[i * N_STATE + k] * N[k * N_CTRL + j];
                return R;
            };

        std::vector<float> sensitivity(static_cast<std::size_t>(N_STATE * NU), 0.0f);
        std::vector<float> nextSensitivity(static_cast<std::size_t>(N_STATE * NU), 0.0f);
        for (int p = 0; p < nPred; ++p)
        {
            std::fill(nextSensitivity.begin(), nextSensitivity.end(), 0.0f);
            for (int i = 0; i < N_STATE; ++i)
            {
                for (int j = 0; j < NU; ++j)
                {
                    float acc = 0.0f;
                    for (int k = 0; k < N_STATE; ++k)
                        acc += Abar[i * N_STATE + k] * sensitivity[k * NU + j];
                    nextSensitivity[i * NU + j] = acc;
                }
            }

            const int cmdBlk = std::min(p, nCtrlH - 1);
            for (int i = 0; i < N_STATE; ++i)
            {
                for (int j = 0; j < N_CTRL; ++j)
                    nextSensitivity[i * NU + (cmdBlk * N_CTRL + j)] += Bbar[i * N_CTRL + j];
            }

            for (int i = 0; i < N_STATE; ++i)
                for (int j = 0; j < NU; ++j)
                    G[(p * N_STATE + i) * NU + j] = nextSensitivity[i * NU + j];

            sensitivity.swap(nextSensitivity);
        }

        // c = predicted state stack under zero absolute command errors.
        // If lidar preview is available, each prediction step uses the
        // corresponding preview disturbance dWindPreview[p]. Otherwise, fall
        // back to the legacy constant-over-horizon disturbance dWind.
        std::array<float, N_STATE> xpred = xbar0;
        for (int p = 0; p < nPred; ++p)
        {
            const float dWindStep = gLidar.available ? dWindPreview[static_cast<std::size_t>(p)] : dWind;
            std::array<float, N_STATE> xnext{};
            for (int i = 0; i < N_STATE; ++i)
            {
                for (int j = 0; j < N_STATE; ++j)
                    xnext[i] += Abar[i * N_STATE + j] * xpred[j];
                xnext[i] += Ebar[i] * dWindStep;
            }

            if (gShutdownRef.loaded)
            {
                const float tPrev = time + static_cast<float>(p) * dtPred;
                const float tNext = time + static_cast<float>(p + 1) * dtPred;
                xnext[0] -= (getRotorSpeedReference(tNext) - getRotorSpeedReference(tPrev));
                xnext[3] -= (getGeneratorTorqueReference(tNext, rotSpeed, horWindV) - getGeneratorTorqueReference(tPrev, rotSpeed, horWindV));
                xnext[4] -= (getPitchReference(tNext) - getPitchReference(tPrev));
            }

            for (int i = 0; i < N_STATE; ++i)
                c[p * N_STATE + i] = xnext[i];
            xpred = xnext;
        }
#ifdef MPC_ENABLE_TIMING
        markTiming(4); // condensed prediction G and zero-input trajectory c
#endif

        MatNN stageQ = makeDiagonalQ();

        // H = 2*(G'Qbar*G + Rbar), g = 2*G'Qbar*c. The terminal
        // block is the full scheduled Riccati matrix, not just its diagonal.
        for (int p = 0; p < nPred; ++p)
        {
            const int base = p * N_STATE;
            const MatNN& blockQ = (p == nPred - 1) ? terminalP : stageQ;
            std::array<float, N_STATE> Qc{};
            std::fill(WG.begin(), WG.end(), 0.0f);

            for (int r = 0; r < N_STATE; ++r)
            {
                for (int s = 0; s < N_STATE; ++s)
                {
                    const float qrs = blockQ[r * N_STATE + s];
                    Qc[r] += qrs * c[base + s];
                    for (int j = 0; j < NU; ++j)
                        WG[r * NU + j] += qrs * G[(base + s) * NU + j];
                }
            }

            for (int i = 0; i < NU; ++i)
            {
                float linearTerm = 0.0f;
                for (int r = 0; r < N_STATE; ++r)
                {
                    const float gri = G[(base + r) * NU + i];
                    linearTerm += gri * Qc[r];
                    for (int j = 0; j < NU; ++j)
                        H[i * NU + j] += static_cast<real_t>(2.0f * gri * WG[r * NU + j]);
                }
                g[i] += static_cast<real_t>(2.0f * linearTerm);
            }
        }
        for (int p = 0; p < nCtrlH; ++p)
        {
            for (int ctrl = 0; ctrl < N_CTRL; ++ctrl)
            {
                const int idx = p * N_CTRL + ctrl;
                const float weight = (ctrl == 0) ? gRT : gRB;
                H[idx * NU + idx] += static_cast<real_t>(2.0f * weight);
                if (p == 0)
                {
                    const float prevCmdErr = (ctrl == 0) ? prevCommandTgErr : prevCommandBetaErr;
                    g[idx] += static_cast<real_t>(-2.0f * weight * prevCmdErr);
                }
                else
                {
                    const int prevIdx = (p - 1) * N_CTRL + ctrl;
                    H[prevIdx * NU + prevIdx] += static_cast<real_t>(2.0f * weight);
                    H[idx * NU + prevIdx] -= static_cast<real_t>(2.0f * weight);
                    H[prevIdx * NU + idx] -= static_cast<real_t>(2.0f * weight);
                }
            }
        }
#ifdef MPC_ENABLE_TIMING
        markTiming(5); // Hessian H and gradient g
#endif

        // Bounds on absolute command errors
        const float dTgRatePred = VS_MAX_TQ_RATE * dtPred;
        const float dBetaRatePred = PC_MAX_RAT * dtPred;

        for (int p = 0; p < nCtrlH; ++p)
        {
            lb[p * N_CTRL + 0] = static_cast<real_t>(VS_MIN_TQ - tgRefNow);
            ub[p * N_CTRL + 0] = static_cast<real_t>(VS_MAX_TQ - tgRefNow);
            lb[p * N_CTRL + 1] = static_cast<real_t>(PC_MIN_PIT - betaRefNow);
            ub[p * N_CTRL + 1] = static_cast<real_t>(PC_MAX_PIT - betaRefNow);
        }

        // Rate constraints on consecutive absolute commands.
        for (int rowBlk = 0; rowBlk < nCtrlH; ++rowBlk)
        {
            const int rowTg = rowBlk * N_CTRL + 0;
            const int rowBeta = rowBlk * N_CTRL + 1;
            const int idxTg = rowBlk * N_CTRL + 0;
            const int idxBeta = rowBlk * N_CTRL + 1;
            Acon[rowTg * NU + idxTg] = 1.0;
            Acon[rowBeta * NU + idxBeta] = 1.0;
            if (rowBlk == 0)
            {
                lbA[rowTg] = static_cast<real_t>(prevCommandTgErr - dTgRatePred);
                ubA[rowTg] = static_cast<real_t>(prevCommandTgErr + dTgRatePred);
                lbA[rowBeta] = static_cast<real_t>(prevCommandBetaErr - dBetaRatePred);
                ubA[rowBeta] = static_cast<real_t>(prevCommandBetaErr + dBetaRatePred);
            }
            else
            {
                Acon[rowTg * NU + ((rowBlk - 1) * N_CTRL + 0)] = -1.0;
                Acon[rowBeta * NU + ((rowBlk - 1) * N_CTRL + 1)] = -1.0;
                lbA[rowTg] = static_cast<real_t>(-dTgRatePred);
                ubA[rowTg] = static_cast<real_t>(dTgRatePred);
                lbA[rowBeta] = static_cast<real_t>(-dBetaRatePred);
                ubA[rowBeta] = static_cast<real_t>(dBetaRatePred);
            }
        }

        // Predicted rotor-speed lower bound. Once the rotor is essentially
        // stopped, this lower bound can conflict with the shutdown target.
        const float omegaConstraintUnloadRpm = std::max(
            gTerminalUnloadTriggerRpm,
            getScheduledTerminalUnloadTriggerRpm(windRef)
        );
        if (rotSpeedRPM > omegaConstraintUnloadRpm)
        {
            for (int p = 0; p < nPred; ++p)
            {
                const int row = NC_INPUT + p;
                const int omegaStateIndex = p * N_STATE + 0;
                for (int j = 0; j < NU; ++j)
                    Acon[row * NU + j] = G[omegaStateIndex * NU + j];

                const float omegaRefPred = gShutdownRef.loaded
                    ? getRotorSpeedReference(time + static_cast<float>(p + 1) * dtPred)
                    : OMEGA_REF;
                lbA[row] = static_cast<real_t>(OMEGA_MIN - omegaRefPred - c[omegaStateIndex]);
                ubA[row] = static_cast<real_t>(BIG_POS);
            }
        }
#ifdef MPC_ENABLE_TIMING
        markTiming(6); // command bounds, slew constraints, and state constraints
#endif

        // Create and solve the QP. SQProblem is required for hotstarting when
        // H and A change between control updates.
        Options options;
        options.printLevel = PL_NONE;
        int_t nWSR = static_cast<int_t>(gNWSR);
        returnValue rv = RET_QP_NOT_SOLVED;
        std::unique_ptr<QProblem> coldQp;
        QProblem* activeQp = nullptr;
        gLastQpSolveMode = 0;

        if (gEnableQpHotstart)
        {
            const bool dimensionsChanged =
                !gHotstartQp || gSolverWs.nu != NU || gSolverWs.nc != NC;
            if (dimensionsChanged)
            {
                gHotstartQp = std::make_unique<SQProblem>(NU, NC);
                gHotstartQpReady = false;
            }
            gHotstartQp->setOptions(options);
            activeQp = gHotstartQp.get();

            if (gHotstartQpReady)
            {
                gLastQpSolveMode = 1;
                rv = gHotstartQp->hotstart(
                    H.data(), g.data(), Acon.data(),
                    lb.data(), ub.data(), lbA.data(), ubA.data(), nWSR
                );
                if (rv != SUCCESSFUL_RETURN)
                {
                    gLastQpSolveMode = 2;
                    gHotstartQp = std::make_unique<SQProblem>(NU, NC);
                    gHotstartQp->setOptions(options);
                    activeQp = gHotstartQp.get();
                    nWSR = static_cast<int_t>(gNWSR);
                    rv = gHotstartQp->init(
                        H.data(), g.data(), Acon.data(),
                        lb.data(), ub.data(), lbA.data(), ubA.data(), nWSR
                    );
                }
            }
            else
            {
                rv = gHotstartQp->init(
                    H.data(), g.data(), Acon.data(),
                    lb.data(), ub.data(), lbA.data(), ubA.data(), nWSR
                );
            }
            gHotstartQpReady = (rv == SUCCESSFUL_RETURN);
        }
        else
        {
            gHotstartQp.reset();
            gHotstartQpReady = false;
            coldQp = std::make_unique<QProblem>(NU, NC);
            coldQp->setOptions(options);
            activeQp = coldQp.get();
            rv = activeQp->init(
                H.data(), g.data(), Acon.data(),
                lb.data(), ub.data(), lbA.data(), ubA.data(), nWSR
            );
        }
        gLastQpNwsrUsed = static_cast<int>(nWSR);
#ifdef MPC_ENABLE_TIMING
        markTiming(7); // qpOASES object setup and active-set solve
#endif
        if (rv != SUCCESSFUL_RETURN)
        {
            std::ostringstream oss;
            oss << "qpOASES init failed: "
                << qpReturnValueToString(rv)
                << " (code=" << static_cast<int>(rv)
                << ", status=" << qpSimpleStatusToString(rv)
                << ", nWSR_used=" << nWSR
                << ", NU=" << NU
                << ", NC=" << NC
                << ", TgRef=" << tgRefNow
                << ")";
            solveErr = oss.str();
            return false;
        }

        rv = activeQp->getPrimalSolution(xOpt.data());
        if (rv != SUCCESSFUL_RETURN)
        {
            std::ostringstream oss;
            oss << "qpOASES getPrimalSolution failed: "
                << qpReturnValueToString(rv)
                << " (code=" << static_cast<int>(rv)
                << ", status=" << qpSimpleStatusToString(rv)
                << ", TgRef=" << tgRefNow
                << ")";
            solveErr = oss.str();
            return false;
        }

        const float cmdTgErr = static_cast<float>(xOpt[0]);
        const float cmdBetaErr = static_cast<float>(xOpt[1]);

        const float dTgRateApply = VS_MAX_TQ_RATE * dtController;
        const float dBetaRateApply = PC_MAX_RAT * dtController;
        const float rawDemandedGenTorque = std::clamp(tgRefNow + cmdTgErr, VS_MIN_TQ, VS_MAX_TQ);
        const float rawDemandedPitchCmd = std::clamp(betaRefNow + cmdBetaErr, PC_MIN_PIT, PC_MAX_PIT);
        demandedGenTorque = std::clamp(rawDemandedGenTorque, prevCommandGenTorque - dTgRateApply, prevCommandGenTorque + dTgRateApply);
        demandedPitchCmd = std::clamp(rawDemandedPitchCmd, prevCommandPitch - dBetaRateApply, prevCommandPitch + dBetaRateApply);
        demandedGenTorque = std::clamp(demandedGenTorque, VS_MIN_TQ, VS_MAX_TQ);
        demandedPitchCmd = std::clamp(demandedPitchCmd, PC_MIN_PIT, PC_MAX_PIT);
        const float dTgApplied = demandedGenTorque - prevCommandGenTorque;
        const float dBetaApplied = demandedPitchCmd - prevCommandPitch;
#ifdef MPC_ENABLE_TIMING
        const auto uReadyTime = markTiming(8); // primal extraction and first command limiting
        gLastTimingUs[10] =
            std::chrono::duration<double, std::micro>(uReadyTime - timingStart).count();
#endif

        float Jomega = 0.0f;
        float Jx = 0.0f;
        float Jv = 0.0f;
        float Jtg = 0.0f;
        float Jbeta = 0.0f;
        float JduTg = 0.0f;
        float JduBeta = 0.0f;
        float Jterminal = 0.0f;
        float predMaxAbsDOmega = 0.0f;
        float predMaxAbsXt = 0.0f;
        float predMaxAbsVt = 0.0f;

        for (int p = 0; p < nPred; ++p)
        {
            const int xBase = p * N_STATE;
            float dOmegaPred = c[xBase + 0];
            float xPred = c[xBase + 1];
            float vPred = c[xBase + 2];
            for (int j = 0; j < NU; ++j)
            {
                dOmegaPred += G[(xBase + 0) * NU + j] * static_cast<float>(xOpt[j]);
                xPred += G[(xBase + 1) * NU + j] * static_cast<float>(xOpt[j]);
                vPred += G[(xBase + 2) * NU + j] * static_cast<float>(xOpt[j]);
            }

            predMaxAbsDOmega = std::max(predMaxAbsDOmega, std::fabs(dOmegaPred));
            predMaxAbsXt = std::max(predMaxAbsXt, std::fabs(xPred));
            predMaxAbsVt = std::max(predMaxAbsVt, std::fabs(vPred));

            if (p == 0 && !gState.hasOneStepPrediction)
            {
                gState.hasOneStepPrediction = true;
                gState.oneStepPredictionIssuedTime = time;
                gState.oneStepPredictionTime = time + dtPred;
                gState.oneStepPredictionDt = dtPred;
                gState.oneStepPredOmegaRef = getRotorSpeedReference(gState.oneStepPredictionTime);
                gState.oneStepPredDOmega = dOmegaPred;
                gState.oneStepPredTowerDispFA = xPred;
                gState.oneStepPredTowerVelFA = vPred;
            }

            float tgErrPred = 0.0f;
            float betaErrPred = 0.0f;
            for (int j = 0; j < NU; ++j)
            {
                tgErrPred += G[(xBase + 3) * NU + j] * static_cast<float>(xOpt[j]);
                betaErrPred += G[(xBase + 4) * NU + j] * static_cast<float>(xOpt[j]);
            }
            tgErrPred += c[xBase + 3];
            betaErrPred += c[xBase + 4];

            if (p == nPred - 1)
            {
                const std::array<float, N_STATE> xTerminal = {
                    dOmegaPred, xPred, vPred, tgErrPred, betaErrPred
                };
                for (int r = 0; r < N_STATE; ++r)
                    for (int s = 0; s < N_STATE; ++s)
                        Jterminal += xTerminal[r] * terminalP[r * N_STATE + s] * xTerminal[s];
            }
            else
            {
                Jomega += gQOmega * dOmegaPred * dOmegaPred;
                Jx += gQX * xPred * xPred;
                Jv += gQV * vPred * vPred;
                Jtg += gQTg * tgErrPred * tgErrPred;
                Jbeta += gQBeta * betaErrPred * betaErrPred;
            }

            appendPredictionLogRow(
                gPredictionLogPath,
                time,
                p + 1,
                time + static_cast<float>(p + 1) * dtPred,
                dOmegaPred,
                xPred,
                vPred,
                x0,
                x1,
                x2,
                horWindV
            );
        }

        for (int p = 0; p < nCtrlH; ++p)
        {
            const float tgCmdErr = static_cast<float>(xOpt[p * N_CTRL + 0]);
            const float betaCmdErr = static_cast<float>(xOpt[p * N_CTRL + 1]);
            const float prevTgCmdErr = (p == 0)
                ? prevCommandTgErr
                : static_cast<float>(xOpt[(p - 1) * N_CTRL + 0]);
            const float prevBetaCmdErr = (p == 0)
                ? prevCommandBetaErr
                : static_cast<float>(xOpt[(p - 1) * N_CTRL + 1]);
            const float duTg = tgCmdErr - prevTgCmdErr;
            const float duBeta = betaCmdErr - prevBetaCmdErr;
            JduTg += gRT * duTg * duTg;
            JduBeta += gRB * duBeta * duBeta;
        }

        const float Jtotal = Jomega + Jx + Jv + Jtg + Jbeta + JduTg + JduBeta + Jterminal;
        appendTuningTraceRow(
            gTuningTracePath,
            time,
            gState.fractureActive,
            true,
            false,
            gState.fallbackCount,
            gCurrentTerminalIndex,
            rotSpeedRPM,
            rotSpeed,
            towerDispFA,
            towerVelFA,
            horWindV,
            dWind,
            tgRefNow,
            betaRefNow,
            prevAppliedGenTorque,
            prevAppliedPitch,
            demandedGenTorque,
            demandedPitchCmd,
            dTgApplied,
            dBetaApplied,
            xbar0,
            Jtotal,
            Jomega,
            Jx,
            Jv,
            Jtg,
            Jbeta,
            JduTg,
            JduBeta,
            Jterminal,
            gQOmega,
            gQX,
            gQV,
            gQTg,
            gQBeta,
            gRT,
            gRB,
            predMaxAbsDOmega,
            predMaxAbsXt,
            predMaxAbsVt,
            "qp"
        );

#ifdef MPC_ENABLE_TIMING
        const auto timingEnd = markTiming(9); // optional prediction/tuning diagnostics
        gLastTimingUs[1] =
            std::chrono::duration<double, std::micro>(timingEnd.time_since_epoch()).count();
        gLastTimingUs[11] =
            std::chrono::duration<double, std::micro>(timingEnd - timingStart).count();
        gLastTimingValid = true;
#endif

        return true;
    }
}

// Returns timestamps and the most recent successful MPC timing breakdown in microseconds.
// Layout: start, end, preprocess, model/workspace, prediction, cost, constraints,
// qpOASES, extract-U, post-diagnostics, total-to-U, total-function.
DLL_EXPORT int MPC_GET_LAST_TIMING(double* values, int capacity)
{
#ifdef MPC_ENABLE_TIMING
    if (values == nullptr || capacity < TIMING_VALUE_COUNT || !gLastTimingValid)
        return 0;
    std::copy(gLastTimingUs.begin(), gLastTimingUs.end(), values);
    return TIMING_VALUE_COUNT;
#else
    (void)values;
    (void)capacity;
    return 0;
#endif
}

DLL_EXPORT int MPC_GET_LAST_SOLVE_INFO(int* values, int capacity)
{
    if (values == nullptr || capacity < 2) return 0;
    values[0] = gLastQpSolveMode;
    values[1] = gLastQpNwsrUsed;
    return 2;
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

#ifdef MPC_ENABLE_TIMING
    gLastTimingValid = false;
#endif

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
    readLidarPreviewData(avrSWAP, gLidar);

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
        gState.targetGenTorque = 0.0f;
        gState.targetPitchCmd = bladePitch1;
        gState.lastPitchRate = 0.0f;
        gState.lastLidarLogTime = -1.0f;
        gState.fallbackCount = 0;
        gState.lastStepUsedFallback = false;
        gState.lastTuningTraceTime = -1.0f;
        gState.lastTuningTraceFallbackState = false;
        gState.lastCommandTraceTime = -1.0f;
        gState.mpcCommandInitialized = false;
        gState.lastMpcSolveTime = time;
        gState.hasOneStepPrediction = false;
        gState.oneStepPredictionIssuedTime = 0.0f;
        gState.oneStepPredictionTime = 0.0f;
        gState.oneStepPredictionDt = 0.0f;
        gState.oneStepPredOmegaRef = 0.0f;
        gState.oneStepPredDOmega = 0.0f;
        gState.oneStepPredTowerDispFA = 0.0f;
        gState.oneStepPredTowerVelFA = 0.0f;
        gState.operatingWindRefInitialized = false;
        gState.operatingWindRef = 0.0f;
        gState.turbulentWindBand = -1;
        gJRFRuntime = J_RF_DEFAULT;
        gJEQRuntime = J_EQ_DEFAULT;
        gFractureMassImbalance = 0.0f;
        gFractureLocationRR = 0.0f;
        gFractureSigma = 1.0f;
        gFractureStatusFromFAST = 0;
        gFractureInertiaLocked = false;
        gLidarHistory.clear();
        gHotstartQp.reset();
        gHotstartQpReady = false;

        std::string loadErr;
        const std::string inFile = cArrayToString(accINFILE, inFileLen);
        const std::string outRoot = cArrayToString(avcOUTNAME, outNameLen);
        gDebugLogPath = replaceExtension(outRoot.empty() ? "DISCON_MPC_CPP" : outRoot, ".lidar_debug.log");
        gPredictionLogPath = replaceExtension(outRoot.empty() ? "DISCON_MPC_CPP" : outRoot, ".mpc_prediction.csv");
        gInputTracePath = replaceExtension(outRoot.empty() ? "DISCON_MPC_CPP" : outRoot, ".mpc_input_trace.csv");
        gTuningTracePath = replaceExtension(outRoot.empty() ? "DISCON_MPC_CPP" : outRoot, ".mpc_tuning_trace.csv");
        gCommandTracePath = replaceExtension(outRoot.empty() ? "DISCON_MPC_CPP" : outRoot, ".mpc_command_trace.csv");
#ifdef MPC_ENABLE_TIMING
        gTimingLogPath = replaceExtension(outRoot.empty() ? "DISCON_MPC_CPP" : outRoot, ".mpc_u_timing.csv");
        const bool timingLogReady = resetTimingLog(gTimingLogPath);
#endif
        const bool tablesLoaded = loadGainTables(inFile, loadErr);
        bool terminalBuilt = false;
        if (tablesLoaded)
            terminalBuilt = buildRuntimeTerminalSchedule(loadErr);
        writeDebugLogHeader(gDebugLogPath);
        writePredictionLogHeader(gPredictionLogPath);
        writeInputTraceHeader(gInputTracePath);
        writeTuningTraceHeader(gTuningTracePath);
        writeCommandTraceHeader(gCommandTracePath);

        if (tablesLoaded && terminalBuilt)
            initializeBaselineStates(time, genSpeed, bladePitch1);

        *aviFAIL = (tablesLoaded && terminalBuilt) ? 1 : -1;
        if (tablesLoaded && terminalBuilt)
        {
            std::ostringstream oss;
            oss << "Running C++ DISCON shell: baseline control before fracture, MPC after fracture"
                << " (N_PRED=" << gNPred << ", N_CTRL_H=" << gNCtrlH
                << ", MPC_DT=" << gPredictionDt
                << ", CTRL_DT=" << gControlDt
                << ", QP_HOTSTART=" << (gEnableQpHotstart ? 1 : 0)
                << ", nWSR=" << gNWSR
                << ", terminalScale=" << gTerminalCostScale
                << ", actuatorOrder=" << gActuatorOrder
                << ", pitchTau=" << gPitchActuatorTau
                << ", torqueTau=" << gTorqueActuatorTau
#ifdef MPC_ENABLE_TIMING
                << ", timingCsv=" << (timingLogReady ? "overwritten" : "disabled-file-busy")
#endif
                << ", terminalFallbackPts=" << gTerminal.fallbackPoints << ").";
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

    const bool fractureInertiaFinalized = updateRuntimeFractureMetricsFromFAST(avrSWAP);
    if (fractureInertiaFinalized && gTable.loaded)
    {
        std::string terminalErr;
        if (!buildRuntimeTerminalSchedule(terminalErr))
        {
            appendDebugLog(gDebugLogPath, "terminal schedule rebuild failed after fracture inertia lock: " + terminalErr);
        }
        else
        {
            appendDebugLog(gDebugLogPath, "locked fracture inertia: J_RF=" + std::to_string(gJRFRuntime) +
                                       ", J_EQ=" + std::to_string(gJEQRuntime) +
                                       ", mass_imbalance=" + std::to_string(gFractureMassImbalance) +
                                       ", fracture_r_R=" + std::to_string(gFractureLocationRR));
        }
    }

    const float dt = std::max(time - gState.lastTime, 1.0e-4f);
    updateOperatingWindReference(horWindV, dt);
    updateLidarHistory(time, horWindV, gLidar);

    if (!gState.fractureActive && time >= gFractureTime)
    {
        gState.fractureActive = true;
    }

    const float prevGenTorqueForTrace = gState.lastGenTorque;
    const float prevPitchCmdForTrace = gState.pitchCmd;
    const float prevAppliedGenTorqueForActuator = gState.lastGenTorque;
    const float prevAppliedPitchForActuator = gState.pitchCmd;
    const float prevTargetGenTorqueForControl =
        (gActuatorOrder == 1) ? gState.targetGenTorque : gState.lastGenTorque;
    const float prevTargetPitchForControl =
        (gActuatorOrder == 1) ? gState.targetPitchCmd : gState.pitchCmd;
    const bool hasPredForTrace =
        gState.hasOneStepPrediction &&
        (time + 1.0e-6f >= gState.oneStepPredictionTime);
    const float predIssuedTimeForTrace = gState.oneStepPredictionIssuedTime;
    const float predTimeForTrace = gState.oneStepPredictionTime;
    const float predDtForTrace = gState.oneStepPredictionDt;
    const float predDOmegaForTrace = gState.oneStepPredDOmega;
    const float actualDOmegaForTrace = rotSpeed - getRotorSpeedReference(time);
    const float predRotSpeedRpmForTrace =
        (gState.oneStepPredDOmega + gState.oneStepPredOmegaRef) * RPS2RPM;
    const float actualRotSpeedRpmForTrace = rotSpeed * RPS2RPM;
    const float predTowerDispFAForTrace = gState.oneStepPredTowerDispFA;
    const float actualTowerDispFAForTrace = towerDispFA;
    const float predTowerVelFAForTrace = gState.oneStepPredTowerVelFA;
    const float actualTowerVelFAForTrace = towerVelFA;
    const float actualDtForTrace = time - predIssuedTimeForTrace;
    if (hasPredForTrace)
        gState.hasOneStepPrediction = false;

    appendInputTraceRow(
        gInputTracePath,
        time,
        iStatus,
        gState.fractureActive,
        rotSpeed,
        horWindV,
        avrSWAP[1018],
        avrSWAP[1019],
        avrSWAP[1020],
        avrSWAP[FRACTURE_J_RF_IDX],
        avrSWAP[FRACTURE_MASS_IMBALANCE_IDX],
        avrSWAP[FRACTURE_LOCATION_IDX],
        avrSWAP[FRACTURE_SIGMA_IDX],
        gFractureStatusFromFAST,
        towerVelFA,
        towerDispFA
    );

    float demandedGenTorque = 0.0f;
    float demandedPitch = bladePitch1;
    bool qpSolvedForTrace = false;
    bool fallbackUsedForTrace = false;

    if (!gState.fractureActive)
    {
        const float savedAppliedGenTorque = gState.lastGenTorque;
        const float savedAppliedPitchCmd = gState.pitchCmd;
        if (gActuatorOrder == 1)
        {
            gState.lastGenTorque = prevTargetGenTorqueForControl;
            gState.pitchCmd = prevTargetPitchForControl;
        }

        runBaselineDISCON(
            iStatus,
            time,
            dt,
            bladePitch1,
            genSpeed,
            demandedGenTorque,
            demandedPitch
        );

        if (gActuatorOrder == 1)
        {
            gState.targetGenTorque = demandedGenTorque;
            gState.targetPitchCmd = demandedPitch;
            gState.lastGenTorque = savedAppliedGenTorque;
            gState.pitchCmd = savedAppliedPitchCmd;
        }
    }
    else
    {
        // avrSWAP time is single precision; at tens of seconds its quantization
        // is already a few microseconds. A 0.1% period tolerance prevents a
        // nominal 10 ms update from slipping to the next 1 ms simulation step.
        const float controlTimeTolerance = std::max(1.0e-6f, 1.0e-3f * gControlDt);
        const bool mpcUpdateDue =
            !gState.mpcCommandInitialized ||
            (time - gState.lastMpcSolveTime) >= (gControlDt - controlTimeTolerance);

        if (!mpcUpdateDue)
        {
            demandedGenTorque = gState.targetGenTorque;
            demandedPitch = gState.targetPitchCmd;
            gState.lastStepUsedFallback = false;
        }
        else
        {
            const float elapsedMpcTime = gState.mpcCommandInitialized
                ? std::max(time - gState.lastMpcSolveTime, 1.0e-4f)
                : gControlDt;
            std::string qpErr;
            const bool qpSolved = solveMultiStepMPC(
                time,
                elapsedMpcTime,
                rotSpeed,
                horWindV,
                towerDispFA,
                towerVelFA,
                gLidar,
                prevAppliedGenTorqueForActuator,
                prevAppliedPitchForActuator,
                prevTargetGenTorqueForControl,
                prevTargetPitchForControl,
                qpErr,
                demandedGenTorque,
                demandedPitch
            );
            qpSolvedForTrace = qpSolved;
#ifdef MPC_ENABLE_TIMING
            if (qpSolved) queueTimingRow(time, demandedGenTorque, demandedPitch);
#endif

            if (!qpSolved)
            {
                applyScheduledTerminalFallback(
                    time,
                    gState.operatingWindRef,
                    rotSpeed,
                    towerDispFA,
                    towerVelFA,
                    prevAppliedGenTorqueForActuator,
                    prevAppliedPitchForActuator,
                    prevTargetGenTorqueForControl,
                    prevTargetPitchForControl,
                    demandedGenTorque,
                    demandedPitch
                );
                gState.lastStepUsedFallback = true;
                fallbackUsedForTrace = true;
                gState.fallbackCount += 1;
                *aviFAIL = 1;
                writeMessage(avcMSG, msgLen, qpErr.empty() ? "qpOASES failed; fallback K used." : (qpErr + " | fallback K used"));

                const float tgRefNow = getGeneratorTorqueReference(time, rotSpeed, gState.operatingWindRef);
                const float betaRefNow = getPitchReference(time);
                const std::array<float, N_STATE> z = buildAugmentedState(
                    time, rotSpeed, gState.operatingWindRef, towerDispFA, towerVelFA, prevAppliedGenTorqueForActuator, prevAppliedPitchForActuator
                );
                appendTuningTraceRow(
                    gTuningTracePath,
                    time,
                    gState.fractureActive,
                    false,
                    true,
                    gState.fallbackCount,
                    gCurrentTerminalIndex,
                    rotSpeed * RPS2RPM,
                    rotSpeed,
                    towerDispFA,
                    towerVelFA,
                    horWindV,
                    0.0f,
                    tgRefNow,
                    betaRefNow,
                    prevAppliedGenTorqueForActuator,
                    prevAppliedPitchForActuator,
                    demandedGenTorque,
                    demandedPitch,
                    demandedGenTorque - prevTargetGenTorqueForControl,
                    demandedPitch - prevTargetPitchForControl,
                    z,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    0.0f,
                    gQOmega,
                    gQX,
                    gQV,
                    gQTg,
                    gQBeta,
                    gRT,
                    gRB,
                    0.0f,
                    0.0f,
                    0.0f,
                    "fallback"
                );
            }
            else
            {
                gState.lastStepUsedFallback = false;
            }

            gState.mpcCommandInitialized = true;
            gState.lastMpcSolveTime = time;
        }
    }

    const float demandedPitchRate = std::clamp((demandedPitch - prevTargetPitchForControl) / dt, -PC_MAX_RAT, PC_MAX_RAT);
    const float rawActuatorDt = time - gState.lastTime;
    const float actuatorDt = (iStatus == 0 || rawActuatorDt <= 0.0f) ? 0.0f : rawActuatorDt;
    const bool actuatorLagActive = (gActuatorOrder == 1 && iStatus != 0);
    const float requestedAppliedGenTorque = actuatorLagActive
        ? applyFirstOrderActuator(prevAppliedGenTorqueForActuator, demandedGenTorque, actuatorDt, gTorqueActuatorTau)
        : demandedGenTorque;
    const float requestedAppliedPitch = actuatorLagActive
        ? applyFirstOrderActuator(prevAppliedPitchForActuator, demandedPitch, actuatorDt, gPitchActuatorTau)
        : demandedPitch;
    const float appliedGenTorque = std::clamp(
        (iStatus == 0)
            ? requestedAppliedGenTorque
            : applyHardRateLimit(
                prevAppliedGenTorqueForActuator,
                requestedAppliedGenTorque,
                VS_MAX_TQ_RATE,
                actuatorDt
            ),
        VS_MIN_TQ,
        VS_MAX_TQ
    );
    const float appliedPitch = std::clamp(
        (iStatus == 0)
            ? requestedAppliedPitch
            : applyHardRateLimit(
                prevAppliedPitchForActuator,
                requestedAppliedPitch,
                PC_MAX_RAT,
                actuatorDt
            ),
        PC_MIN_PIT,
        PC_MAX_PIT
    );
    const float appliedPitchRate = (actuatorDt > 0.0f)
        ? (appliedPitch - prevAppliedPitchForActuator) / actuatorDt
        : 0.0f;
    const float tgRefForTrace = getGeneratorTorqueReference(time, rotSpeed, gState.operatingWindRef);
    const float betaRefForTrace = getPitchReference(time);
    appendCommandTraceRow(
        gCommandTracePath,
        time,
        iStatus,
        gState.fractureActive,
        qpSolvedForTrace,
        fallbackUsedForTrace,
        gState.fallbackCount,
        gCurrentTerminalIndex,
        rotSpeed * RPS2RPM,
        genSpeed * RPS2RPM,
        towerDispFA,
        towerVelFA,
        horWindV,
        gState.operatingWindRef,
        tgRefForTrace,
        betaRefForTrace,
        prevTargetGenTorqueForControl,
        prevTargetPitchForControl,
        demandedGenTorque,
        demandedPitch,
        demandedPitchRate,
        appliedGenTorque,
        appliedPitch,
        appliedPitchRate,
        hasPredForTrace,
        predIssuedTimeForTrace,
        predTimeForTrace,
        predDOmegaForTrace,
        actualDOmegaForTrace,
        predRotSpeedRpmForTrace,
        actualRotSpeedRpmForTrace,
        predTowerDispFAForTrace,
        actualTowerDispFAForTrace,
        predTowerVelFAForTrace,
        actualTowerVelFAForTrace,
        predDtForTrace,
        actualDtForTrace
    );

    avrSWAP[34] = 1.0f;            // avrSWAP(35)  generator contactor: main variable-speed generator
    avrSWAP[46] = appliedGenTorque; // avrSWAP(47) demanded generator torque
    avrSWAP[55] = 0.0f;            // avrSWAP(56) torque override = yes

    avrSWAP[41] = appliedPitch;    // avrSWAP(42) blade 1 pitch command
    avrSWAP[42] = appliedPitch;    // avrSWAP(43) blade 2 pitch command
    avrSWAP[43] = appliedPitch;    // avrSWAP(44) blade 3 pitch command
    avrSWAP[44] = appliedPitch;    // avrSWAP(45) collective pitch command
    avrSWAP[54] = 0.0f;            // avrSWAP(55) pitch override = yes
    avrSWAP[45] = appliedPitchRate; // avrSWAP(46) demanded collective pitch rate

    gState.lastTime = time;
    gState.lastGenTorque = appliedGenTorque;
    gState.pitchCmd = appliedPitch;
    gState.targetGenTorque = demandedGenTorque;
    gState.targetPitchCmd = demandedPitch;
    gState.lastPitchRate = appliedPitchRate;

#ifdef MPC_ENABLE_TIMING
    if (iStatus < 0) flushTimingRows();
#endif

    if (gState.fractureActive)
    {
        const float lidarMeanWind = getLidarMeanWind(gLidar, horWindV);
        const float lidarPreviewDelta = lidarMeanWind - horWindV;
        const bool usingLidarPreview = !gLidarHistory.empty();

        if (gState.lastLidarLogTime < 0.0f || (time - gState.lastLidarLogTime) >= LIDAR_LOG_DT)
        {
            std::ostringstream dbg;
            dbg << "time=" << time
                << ", wind=" << horWindV
                << ", windRef=" << gState.operatingWindRef
                << ", rotSpeedRpm=" << rotSpeed * RPS2RPM
                << ", lidarAvail=" << (gLidar.available ? 1 : 0)
                << ", lidarPts=" << gLidar.measuredSpeeds.size()
                << ", lidarMean=" << lidarMeanWind
                << ", lidarDelta=" << lidarPreviewDelta
                << ", previewSource=" << (usingLidarPreview ? "lidar-history" : "fallback-constant")
                << ", terminalIdx=" << gCurrentTerminalIndex
                << ", fallbackUsed=" << (gState.lastStepUsedFallback ? 1 : 0)
                << ", fallbackCount=" << gState.fallbackCount
                << ", fractureStatusFAST=" << gFractureStatusFromFAST
                << ", J_RF=" << gJRFRuntime
                << ", J_EQ=" << gJEQRuntime
                << ", massImbalance=" << gFractureMassImbalance
                << ", fracture_r_R=" << gFractureLocationRR
                << ", fractureSigma=" << gFractureSigma
                << ", Tg=" << demandedGenTorque
                << ", betaDeg=" << demandedPitch * R2D;
            appendDebugLog(gDebugLogPath, dbg.str());
            gState.lastLidarLogTime = time;
        }
    }
}

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <deque>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <ctime>

#include <qpOASES.hpp>
#include "terminal_mpc_schedule.hpp"

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
        float lastLidarLogTime = -1.0f;
        int fallbackCount = 0;
        bool lastStepUsedFallback = false;
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
    LidarPreviewData gLidar;
    std::deque<LidarHistoryPoint> gLidarHistory;
    float gFractureTime = 30.0f;
    std::string gDebugLogPath;
    std::string gPredictionLogPath;
    std::string gInputTracePath;
    bool gEnableTraceFiles = true;
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
    constexpr float OMEGA_MIN = -0.2f;         // hard lower bound to avoid reverse rotation
    constexpr float TG_REF = VS_MIN_TQ;        // shutdown target generator torque [N-m]
    constexpr float BETA_REF = PC_MAX_PIT;     // shutdown target collective pitch [rad]

    // Continuous-time physical parameters for the simplified shutdown MPC model.
    // Updated from the user's latest identified values.
    constexpr float J_RF = 3.19609962e7f;      // kg m^2  (identified rotor-equivalent inertia at 70%)
    constexpr float J_generator = 534.116f;    // Generator inertia about HSS (kg m^2)

    constexpr float N_gear = 97.0f;            // gearbox ratio
    constexpr float J_EQ = J_RF + J_generator * N_gear * N_gear; // 等效转动惯量
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
    constexpr int   LIDAR_MSR_START = 2000;  // C index for avrSWAP(2001)
    constexpr int   LIDAR_MAX_CHAN  = 500;
    constexpr float LIDAR_LOG_DT = 0.1f;     // seconds between lidar debug log entries
    constexpr float LIDAR_MIN_CONV_WIND = 1.0f;   // minimum convection speed used in Taylor mapping [m/s]
    constexpr float LIDAR_EVOLUTION_LENGTH = 300.0f; // decay length for preview confidence [m]
    constexpr float LIDAR_HISTORY_MAX_AGE = 30.0f;   // seconds of lidar history kept for preview mapping

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
    float gOmegaErrMax = 2.0f;   // rad/s
    float gTowerDispMax = 0.5f;   // m
    float gTowerVelMax = 0.5f;   // m/s

    constexpr int   N_STATE = 5;
    constexpr int   N_CTRL = 2;
    constexpr float BIG_NEG = -1.0e20f;

    using MatNN = std::array<float, N_STATE * N_STATE>;
    using MatNU = std::array<float, N_STATE * N_CTRL>;
    using MatUN = std::array<float, N_CTRL * N_STATE>;
    using MatUU = std::array<float, N_CTRL * N_CTRL>;

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

    inline int lookupTerminalScheduleIndex(float windNow, float speedNowRpm)
    {
        auto nearestIndex = [](const auto& arr, float x) -> int
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

        const int iw = nearestIndex(kTerminalWindPoints, windNow);
        const int is = nearestIndex(kTerminalSpeedPoints, speedNowRpm);
        const int tableIdx = iw * kTerminalNumSpeed + is;
        return kTerminalPointToIndex[tableIdx];
    }

    inline MatNN getScheduledTerminalP(float windNow, float speedNowRpm)
    {
        MatNN P{};
        const int idx = lookupTerminalScheduleIndex(windNow, speedNowRpm);
        const int offset = idx * (N_STATE * N_STATE);
        for (int i = 0; i < N_STATE * N_STATE; ++i)
            P[i] = kTerminalPTable[static_cast<std::size_t>(offset + i)];
        return P;
    }

    inline MatUN getScheduledTerminalK(float windNow, float speedNowRpm)
    {
        MatUN K{};
        const int idx = lookupTerminalScheduleIndex(windNow, speedNowRpm);
        const int offset = idx * (N_CTRL * N_STATE);
        for (int i = 0; i < N_CTRL * N_STATE; ++i)
            K[i] = kTerminalKTable[static_cast<std::size_t>(offset + i)];
        return K;
    }

    inline std::array<float, N_STATE> buildAugmentedState(
        float time,
        float rotSpeed,
        float towerDispFA,
        float towerVelFA,
        float prevGenTorque,
        float prevPitchCmd)
    {
        (void)time;
        (void)rotSpeed;
        const float tgRefNow = TG_REF;
        return {
            rotSpeed - OMEGA_REF,
            towerDispFA,
            towerVelFA,
            prevGenTorque - tgRefNow,
            prevPitchCmd - BETA_REF
        };
    }

    inline void applyScheduledTerminalFallback(
        float time,
        float horWindV,
        float rotSpeed,
        float towerDispFA,
        float towerVelFA,
        float prevGenTorque,
        float prevPitchCmd,
        float& demandedGenTorque,
        float& demandedPitchCmd)
    {
        const float rotSpeedRPM = rotSpeed * RPS2RPM;
        const MatUN K = getScheduledTerminalK(horWindV, rotSpeedRPM);
        const auto z = buildAugmentedState(time, rotSpeed, towerDispFA, towerVelFA, prevGenTorque, prevPitchCmd);

        float deltaTg = 0.0f;
        float deltaBeta = 0.0f;
        for (int j = 0; j < N_STATE; ++j)
        {
            deltaTg += K[0 * N_STATE + j] * z[j];
            deltaBeta += K[1 * N_STATE + j] * z[j];
        }

        demandedGenTorque = std::clamp(prevGenTorque + deltaTg, VS_MIN_TQ, VS_MAX_TQ);
        demandedPitchCmd = std::clamp(prevPitchCmd + deltaBeta, PC_MIN_PIT, PC_MAX_PIT);
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
        if (!loadCsvTable(makePath(lines[5]), nWind, nSpeed, gTable.GTbeta)) { err = "Failed loading GTbeta table.";  return false; }
        if (!loadCsvTable(makePath(lines[6]), nWind, nSpeed, gTable.GTu)) { err = "Failed loading GTu table.";     return false; }
        if (!loadCsvTable(makePath(lines[7]), nWind, nSpeed, gTable.GFomega)) { err = "Failed loading GFomega table."; return false; }
        if (!loadCsvTable(makePath(lines[8]), nWind, nSpeed, gTable.GFbeta)) { err = "Failed loading GFbeta table.";  return false; }
        if (!loadCsvTable(makePath(lines[9]), nWind, nSpeed, gTable.GFu)) { err = "Failed loading GFu table.";     return false; }

        if (lines.size() >= 11) gFractureTime = std::stof(lines[10]);

        // Backward-compatible parsing:
        // old format:  19 data lines, horizons fixed at 5/5 in code
        // mid format:  21 data lines, lines[11:12] are N_PRED/N_CTRL_H
        // new format:  23 data lines, adds Q_TG and Q_BETA before R_T/R_B
        // newest format: 24 data lines, line 24 toggles trace-file generation
        std::size_t idx = 11;
        gNPred = 5;
        gNCtrlH = 5;
        gQTg = 1.0e-6f;
        gQBeta = 10.0f;
        gEnableTraceFiles = true;
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

        if (lines.size() >= 23)
        {
            if (lines.size() > idx) gQOmega = std::stof(lines[idx++]);
            if (lines.size() > idx) gQX = std::stof(lines[idx++]);
            if (lines.size() > idx) gQV = std::stof(lines[idx++]);
            if (lines.size() > idx) gQTg = std::stof(lines[idx++]);
            if (lines.size() > idx) gQBeta = std::stof(lines[idx++]);
            if (lines.size() > idx) gRT = std::stof(lines[idx++]);
            if (lines.size() > idx) gRB = std::stof(lines[idx++]);
            if (lines.size() > idx) gOmegaErrMax = std::stof(lines[idx++]);
            if (lines.size() > idx) gTowerDispMax = std::stof(lines[idx++]);
            if (lines.size() > idx) gTowerVelMax = std::stof(lines[idx++]);
        }
        else
        {
            if (lines.size() > idx) gQOmega = std::stof(lines[idx++]);
            if (lines.size() > idx) gQX = std::stof(lines[idx++]);
            if (lines.size() > idx) gQV = std::stof(lines[idx++]);
            if (lines.size() > idx) gRT = std::stof(lines[idx++]);
            if (lines.size() > idx) gRB = std::stof(lines[idx++]);
            if (lines.size() > idx) gOmegaErrMax = std::stof(lines[idx++]);
            if (lines.size() > idx) gTowerDispMax = std::stof(lines[idx++]);
            if (lines.size() > idx) gTowerVelMax = std::stof(lines[idx++]);
        }

        if (lines.size() > idx)
        {
            gEnableTraceFiles = (std::stoi(lines[idx++]) != 0);
        }

        gTable.loaded = true;
        return true;
    }

    inline float clampToRange(float x, float lo, float hi)
    {
        return std::max(lo, std::min(x, hi));
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
            << ", Q=[" << gQOmega << "," << gQX << "," << gQV << "," << gQTg << "," << gQBeta << "]"
            << ", R=[" << gRT << "," << gRB << "]"
            << ", limits=[omegaErrMax=" << gOmegaErrMax
            << ", towerDispMax=" << gTowerDispMax
            << ", towerVelMax=" << gTowerVelMax << "]";
        appendDebugLog(path, cfg.str());

        std::ostringstream sched;
        sched << "# terminal_schedule_points=" << kTerminalNumPoints
              << ", wind_grid=" << kTerminalNumWind
              << ", speed_grid=" << kTerminalNumSpeed;
        appendDebugLog(path, sched.str());
    }

    inline void writePredictionLogHeader(const std::string& path)
    {
        if (!gEnableTraceFiles) return;
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

    inline void writeInputTraceHeader(const std::string& path)
    {
        if (!gEnableTraceFiles) return;
        if (path.empty()) return;
        resetDebugLog(path);
        std::ofstream out(path, std::ios::app);
        if (!out) return;
        out << "time,iStatus,fracture_active,rotSpeed,horWindV,"
               "avr_1019_towerVelFA,avr_1020_towerVelSS,avr_1021_towerDispFA,"
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
        float towerVelFA,
        float towerDispFA)
    {
        if (!gEnableTraceFiles) return;
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
            << rotSpeed << ','
            << towerVelFA << ','
            << towerDispFA << '\n';
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

    inline float getGeneratorTorqueReference(float time, float rotSpeed)
    {
        const float tau = time - gFractureTime;
        if (tau <= 0.0f)
            return TG_REF;
        const float rotSpeedRPM = rotSpeed * RPS2RPM;
        if (rotSpeedRPM > 3.0f)
            return VS_MAX_TQ;
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
        float dt,
        float rotSpeed,
        float horWindV,
        float towerDispFA,
        float towerVelFA,
        const LidarPreviewData& lidar,
        float prevGenTorque,
        float prevPitchCmd,
        std::string& solveErr,
        float& demandedGenTorque,
        float& demandedPitchCmd)
    {
        const int nPred = gNPred;
        const int nCtrlH = gNCtrlH;
        const float rotSpeedRPM = rotSpeed * 9.5492966f;
        gCurrentTerminalIndex = lookupTerminalScheduleIndex(horWindV, rotSpeedRPM);
        const float tgRefNow = getGeneratorTorqueReference(time, rotSpeed);
        const float windRef = gTable.loaded ? gTable.windPtsMs[std::min_element(gTable.windPtsMs.begin(), gTable.windPtsMs.end(),
            [&](float a, float b) { return std::fabs(a - horWindV) < std::fabs(b - horWindV); }) - gTable.windPtsMs.begin()] : 11.4f;
        const float dWind = horWindV - windRef;
        bool usedLidarHistory = false;
        const std::vector<float> dWindPreview = buildWindPreviewDeltas(time, horWindV, nPred, dt, usedLidarHistory);

        const float gTOmegaNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GTomega, horWindV, rotSpeedRPM) : G_TOMEGA_DEFAULT;
        const float gTBetaNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GTbeta, horWindV, rotSpeedRPM) : G_TBETA_DEFAULT;
        const float gTUNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GTu, horWindV, rotSpeedRPM) : G_TU_DEFAULT;
        const float gFOmegaNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GFomega, horWindV, rotSpeedRPM) : G_FOMEGA_DEFAULT;
        const float gFBetaNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GFbeta, horWindV, rotSpeedRPM) : G_FBETA_DEFAULT;
        const float gFUNow = gTable.loaded ? interp2d(gTable.windPtsMs, gTable.speedPtsRpm, gTable.GFu, horWindV, rotSpeedRPM) : G_FU_DEFAULT;

        // Augmented state:
        // xbar = [ dOmega, x_t, v_t, Tg - Tg_ref(time, omega), beta - beta_ref ]^T
        const float x0 = rotSpeed - OMEGA_REF;
        const float x1 = towerDispFA;
        const float x2 = towerVelFA;
        const float x3 = prevGenTorque - tgRefNow;
        const float x4 = prevPitchCmd - BETA_REF;

        std::array<float, N_STATE> xbar0 = { x0, x1, x2, x3, x4 };

        // One-step Euler-discretized augmented model
        MatNN Abar{};
        MatNU Bbar{};
        std::array<float, N_STATE> Ebar{};
        auto Aat = [&](int r, int c) -> float& { return Abar[r * N_STATE + c]; };
        auto Bat = [&](int r, int c) -> float& { return Bbar[r * N_CTRL + c]; };

        const float a11 = 1.0f + dt * gTOmegaNow / J_EQ;
        const float a22 = 1.0f;
        const float a23 = dt;
        const float a31 = dt * gFOmegaNow / M_T;
        const float a32 = -dt * K_T / M_T;
        const float a33 = 1.0f - dt * C_T / M_T;
        const float b11 = -dt * N_gear / J_EQ;
        const float b12 = dt * gTBetaNow / J_EQ;
        const float b32 = dt * gFBetaNow / M_T;
        const float e11 = dt * gTUNow / J_EQ;
        const float e31 = dt * gFUNow / M_T;

        Aat(0, 0) = a11;  Aat(0, 3) = b11;  Aat(0, 4) = b12;
        Aat(1, 1) = a22;  Aat(1, 2) = a23;
        Aat(2, 0) = a31;  Aat(2, 1) = a32;  Aat(2, 2) = a33;  Aat(2, 4) = b32;
        Aat(3, 3) = 1.0f;
        Aat(4, 4) = 1.0f;

        Bat(0, 0) = b11;  Bat(0, 1) = b12;
        Bat(1, 0) = 0.0f; Bat(1, 1) = 0.0f;
        Bat(2, 0) = 0.0f; Bat(2, 1) = b32;
        Bat(3, 0) = 1.0f; Bat(3, 1) = 0.0f;
        Bat(4, 0) = 0.0f; Bat(4, 1) = 1.0f;

        Ebar[0] = e11;
        Ebar[1] = 0.0f;
        Ebar[2] = e31;
        Ebar[3] = 0.0f;
        Ebar[4] = 0.0f;

        const MatNN terminalP = getScheduledTerminalP(horWindV, rotSpeedRPM);
        const MatUN terminalK = getScheduledTerminalK(horWindV, rotSpeedRPM);
        (void)terminalK;

        // Build prediction matrices X = F*x0 + G*U
        const int NX = N_STATE * nPred;
        const int NU = N_CTRL * nCtrlH;
        std::vector<float> F(NX * N_STATE, 0.0f);
        std::vector<float> G(NX * NU, 0.0f);

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

        std::array<float, N_STATE* N_STATE> Apow{};
        for (int i = 0; i < N_STATE; ++i) Apow[i * N_STATE + i] = 1.0f;

        std::vector<std::array<float, N_STATE* N_STATE>> powers(static_cast<std::size_t>(nPred));
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
                std::array<float, N_STATE* N_CTRL> block{};
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

        // c = predicted state stack under zero control increments.
        // If lidar preview is available, each prediction step uses the
        // corresponding preview disturbance dWindPreview[p]. Otherwise, fall
        // back to the legacy constant-over-horizon disturbance dWind.
        std::vector<float> c(NX, 0.0f);
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
                c[p * N_STATE + i] = xnext[i];
            }
            xpred = xnext;
        }

        // Block-diagonal Q and R
        std::vector<float> Qblk(NX * NX, 0.0f);
        std::vector<float> Rblk(NU * NU, 0.0f);

        const MatNN stageQ = makeDiagonalQ();
        for (int p = 0; p < nPred; ++p)
        {
            const int base = p * N_STATE;
            const MatNN& blockQ = (p == nPred - 1) ? terminalP : stageQ;
            for (int i = 0; i < N_STATE; ++i)
                for (int j = 0; j < N_STATE; ++j)
                    Qblk[(base + i) * NX + (base + j)] = blockQ[i * N_STATE + j];
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
            ub[p * N_CTRL + 0] = static_cast<real_t>(dTgRate);
            lb[p * N_CTRL + 1] = static_cast<real_t>(-dBetaRate);
            ub[p * N_CTRL + 1] = static_cast<real_t>(dBetaRate);
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
        // dOmega uses an asymmetric bound:
        //   OMEGA_MIN <= rotSpeed = dOmega + OMEGA_REF <= gOmegaErrMax + OMEGA_REF
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
            const float cDisp = c[xBase + 1];
            const float cVel = c[xBase + 2];

            ubA[upRow + 0] = static_cast<real_t>(gOmegaErrMax - cOmega);
            ubA[upRow + 1] = static_cast<real_t>(gTowerDispMax - cDisp);
            ubA[upRow + 2] = static_cast<real_t>(gTowerVelMax - cVel);

            //ubA[lowRow + 0] = static_cast<real_t>( cOmega - OMEGA_MIN );
            ubA[lowRow + 0] = static_cast<real_t>(gOmegaErrMax + cOmega);
            ubA[lowRow + 1] = static_cast<real_t>(gTowerDispMax + cDisp);
            ubA[lowRow + 2] = static_cast<real_t>(gTowerVelMax + cVel);
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

        std::vector<real_t> xOpt(NU, 0.0);
        rv = qp.getPrimalSolution(xOpt.data());
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

        const float dTg = static_cast<float>(xOpt[0]);
        const float dBeta = static_cast<float>(xOpt[1]);

        demandedGenTorque = std::clamp(prevGenTorque + dTg, VS_MIN_TQ, VS_MAX_TQ);
        demandedPitchCmd = std::clamp(prevPitchCmd + dBeta, PC_MIN_PIT, PC_MAX_PIT);

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

            appendPredictionLogRow(
                gPredictionLogPath,
                time,
                p + 1,
                time + static_cast<float>(p + 1) * dt,
                dOmegaPred,
                xPred,
                vPred,
                x0,
                x1,
                x2,
                horWindV
            );
        }

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
        gState.lastPitchRate = 0.0f;
        gState.lastLidarLogTime = -1.0f;
        gState.fallbackCount = 0;
        gState.lastStepUsedFallback = false;
        gLidarHistory.clear();

        std::string loadErr;
        const std::string inFile = cArrayToString(accINFILE, inFileLen);
        const std::string outRoot = cArrayToString(avcOUTNAME, outNameLen);
        gDebugLogPath = replaceExtension(outRoot.empty() ? "DISCON_MPC_CPP" : outRoot, ".lidar_debug.log");
        gPredictionLogPath = replaceExtension(outRoot.empty() ? "DISCON_MPC_CPP" : outRoot, ".mpc_prediction.csv");
        gInputTracePath = replaceExtension(outRoot.empty() ? "DISCON_MPC_CPP" : outRoot, ".mpc_input_trace.csv");
        const bool tablesLoaded = loadGainTables(inFile, loadErr);
        writeDebugLogHeader(gDebugLogPath);
        writePredictionLogHeader(gPredictionLogPath);
        writeInputTraceHeader(gInputTracePath);

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
    updateLidarHistory(time, horWindV, gLidar);

    if (!gState.fractureActive && time >= gFractureTime)
    {
        gState.fractureActive = true;
    }

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
        towerVelFA,
        towerDispFA
    );

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
        std::string qpErr;
        const bool qpSolved = solveMultiStepMPC(
            time,
            dt,
            rotSpeed,
            horWindV,
            towerDispFA,
            towerVelFA,
            gLidar,
            gState.lastGenTorque,
            gState.pitchCmd,
            qpErr,
            demandedGenTorque,
            demandedPitch
        );

        if (!qpSolved)
        {
            applyScheduledTerminalFallback(
                time,
                horWindV,
                rotSpeed,
                towerDispFA,
                towerVelFA,
                gState.lastGenTorque,
                gState.pitchCmd,
                demandedGenTorque,
                demandedPitch
            );
            gState.lastStepUsedFallback = true;
            gState.fallbackCount += 1;
            *aviFAIL = 1;
            writeMessage(avcMSG, msgLen, qpErr.empty() ? "qpOASES failed; fallback K used." : (qpErr + " | fallback K used"));
        }
        else
        {
            gState.lastStepUsedFallback = false;
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
                << ", rotSpeedRpm=" << rotSpeed * RPS2RPM
                << ", lidarAvail=" << (gLidar.available ? 1 : 0)
                << ", lidarPts=" << gLidar.measuredSpeeds.size()
                << ", lidarMean=" << lidarMeanWind
                << ", lidarDelta=" << lidarPreviewDelta
                << ", previewSource=" << (usingLidarPreview ? "lidar-history" : "fallback-constant")
                << ", terminalIdx=" << gCurrentTerminalIndex
                << ", fallbackUsed=" << (gState.lastStepUsedFallback ? 1 : 0)
                << ", fallbackCount=" << gState.fallbackCount
                << ", Tg=" << demandedGenTorque
                << ", betaDeg=" << demandedPitch * R2D;
            appendDebugLog(gDebugLogPath, dbg.str());
            gState.lastLidarLogTime = time;
        }
    }
}

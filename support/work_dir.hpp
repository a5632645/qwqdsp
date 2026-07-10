#pragma once

#include <filesystem>

namespace qwqdsp_support {

// ------------------------------------------------------------
// folder
// ------------------------------------------------------------

inline std::filesystem::path GetWorkDir() {
#ifdef QWQDSP_WORK_DIR
    return std::filesystem::absolute(QWQDSP_WORK_DIR);
#else
    #error "cmake config wrongly"
#endif
}

inline std::filesystem::path GetInputDir() {
    return GetWorkDir() / "input";
}

inline std::filesystem::path GetOutputDir() {
    return GetWorkDir() / "output";
}

// ------------------------------------------------------------
// file
// ------------------------------------------------------------

inline std::string WormholeWav() {
    return (GetInputDir() / "wormhole.wav").string();
}

inline std::string SweepWav() {
    return (GetInputDir() / "sweep.wav").string();
}

inline std::string InputFile(std::string str) {
    return (GetInputDir() / str).string();
}

inline std::string OutputFile(std::string str) {
    return (GetOutputDir() / str).string();
}

} // namespace qwqdsp_support

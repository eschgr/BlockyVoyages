#pragma once

#include <vector>
#include <string>
#include <fstream>

#include "../types.h"

namespace BlockyVoyages {
namespace Debugability {

#define LOG(LVL, STR, ...)     BlockyVoyages::Debugability::logger->writeln(__FILE__, __LINE__, __FUNCTION__, LVL, STR, __VA_ARGS__)

enum DebugLevel {
    DBG_LVL_MESSAGE = 0,     // Normal Passive logging, will always log the message
    DBG_LVL_NORMAL,          // Normal Passive logging, may be cut out by minimum level
    DBG_LVL_LOW,
	DBG_LVL_WARNING,         // Logs the message with a warning tag
	DBG_LVL_MEDIUM,
	DBG_LVL_HIGH,
    DBG_LVL_CRITICAL,        // Logs the message with a critical tag
	DBG_LVL_FATAL,           // Logs the message as a fatal error. Probably will crash just after this
	DBG_LVL_ALWAYS_TERMINAL, // fatal error that should prevent the program from ever proceeding
    DBG_LVL_COUNT            // the number of log levels supported. Always the last one
};

class Logger {
public:
    Logger();
    Logger(const std::string& logfile);
    ~Logger();

    // Logging functionality
    void write(DebugLevel dbgLevel, char* format, ...);
    void writeln(DebugLevel dbgLevel, char* format, ...);

    void write(char* file, int32 line, char* func, DebugLevel dbgLevel, char* format, ...);
    void writeln(char* file, int32 line, char* func, DebugLevel dbgLevel, char* format, ...);

    void setMinimumRecordLevel(DebugLevel newMin) { m_minLogLevel = newMin; }
    void setCriticalDebugLevel(DebugLevel newMax) { m_critDbgLvl = newMax; }

    bool criticalMessageSeen(void) { return m_criticalStop; }
private:
    std::string m_log_file_name;
    std::fstream m_logFile;
    bool m_criticalStop;

    DebugLevel m_minLogLevel;
    DebugLevel m_critDbgLvl;

    void Initialize();
    inline void WriteToLog(DebugLevel level, const std::string& text) {
        WriteToLog(level, text, nullptr, -1, nullptr);
    }
    void WriteToLog(DebugLevel level,
                    const std::string& text,
                    const char* file,
                    int32 line,
                    const char* func);
};

extern Logger* logger;

}
}
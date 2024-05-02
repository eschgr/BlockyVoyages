#include "stdafx.h"

#include "Logger.h"

#include <Windows.h>
#include <sstream>
#include <cctype>

namespace BlockyVoyages {
namespace Debugability {
Logger* logger = nullptr;

Logger::Logger()
    : m_minLogLevel(DBG_LVL_MESSAGE),
      m_critDbgLvl(DBG_LVL_CRITICAL),
      m_criticalStop(false),
      m_log_file_name("Logs/main_log.txt") {
    Initialize();
}

Logger::Logger(const std::string& log_file)
    : m_minLogLevel(DBG_LVL_MESSAGE),
      m_critDbgLvl(DBG_LVL_CRITICAL),
      m_criticalStop(false),
      m_log_file_name(log_file) {
    Initialize();
}

Logger::~Logger() {
    WriteToLog(DBG_LVL_NORMAL, "LOGGER: Log file closed.");

    m_logFile.close();
    m_logFile.clear();
}

void Logger::write(DebugLevel dbgLevel, char* format, ...) {
    va_list args;
    char buffer[255];

    // create a string from the arguments sent
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    WriteToLog(dbgLevel, buffer, nullptr, -1, nullptr);
}

void Logger::writeln(DebugLevel dbgLevel, char* format, ...) {
    // create a string from the arguments sent
    va_list args;
    char buffer[255];
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    WriteToLog(dbgLevel, buffer, nullptr, -1, nullptr);
}

void Logger::write(char* file, int32 line, char* func, DebugLevel dbgLevel, char* format, ...) {
    // create a string from the arguments sent
    va_list args;
    char buffer[255];
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    WriteToLog(dbgLevel, buffer, file, line, func);
}

void Logger::writeln(char* file, int32 line, char* func,DebugLevel dbgLevel, char* format, ...) {
    // create a string from the arguments sent
    va_list args;
    char buffer[1024];
    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    // send out the message
    WriteToLog(dbgLevel, buffer, file, line, func);
}

void Logger::WriteToLog(DebugLevel level, const std::string& text, const char* file, int32 line, const char* func) {
    if (level != DBG_LVL_MESSAGE && level < m_minLogLevel) {
        return;
    }
    std::stringstream sStream;
    switch(level) {
        case DBG_LVL_WARNING:
        case DBG_LVL_MEDIUM:
            m_logFile << "Warning: ";
            break;
        case DBG_LVL_HIGH:
        case DBG_LVL_CRITICAL:
            m_logFile << "Critical Error: ";
            break;
        case DBG_LVL_FATAL:
        case DBG_LVL_ALWAYS_TERMINAL:
            m_logFile << "Fatal Error: ";
            break;
    }
    if (level > DBG_LVL_MEDIUM) {
        if (nullptr != file) {
            m_logFile << file;
            if (line >= 0) {
                m_logFile << "(" << line << ")";
            }
            m_logFile << ": ";
        }
        if (nullptr != func) {
            m_logFile << func << ": ";
        }
    }

    m_logFile << text << std::endl;
    m_logFile.flush();
    m_criticalStop = m_criticalStop || level >= m_critDbgLvl;
}

void Logger::Initialize() {
    m_logFile.open(m_log_file_name, std::ios_base::out | std::ios_base::trunc);
    WriteToLog(DBG_LVL_NORMAL, "LOGGER: Log file opened.");
}
}
}
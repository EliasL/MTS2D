#pragma once
#ifndef LOGGING_H
#define LOGGING_H

#include <iostream>

#ifdef ENABLE_LOGGING
    #ifdef WRITE_LOG_TO_FILE
        #include <fstream>
        extern std::ofstream logFile;
        #define LOG(x) logToFile(x)
    #else
        #define LOG(x) log(x)
    #endif
#else
    #define LOG(x)
#endif

// Forward declaration for logToFile
template <typename T, typename... Args>
void logToFile(const T &firstArg, const Args&... args);

template <typename T>
void logValue(const T &value) {
    std::cout << value;
}

template <typename T, typename... Args>
void logValue(const T &value, const Args&... args) {
    std::cout << value << " ";
    logValue(args...);
}

template <typename T, typename... Args>
void log(const T &firstArg, const Args&... args) {
    logValue(firstArg, args...);
    std::cout << std::endl;
}

#ifdef WRITE_LOG_TO_FILE
extern std::ofstream logFile;

template <typename T>
void logValueToFile(const T &value) {
    logFile << value;
}

template <typename T, typename... Args>
void logValueToFile(const T &value, const Args&... args) {
    logFile << value << " ";
    logValueToFile(args...);
}

template <typename T, typename... Args>
void logToFile(const T &firstArg, const Args&... args) {
    logValueToFile(firstArg, args...);
    logFile << std::endl;
}
#endif

#endif

#ifndef LOGGING_H
#define LOGGING_H
#pragma once

#include <string>

#include "settings.h"
#include "Data/dataExport.h"

#include <chrono>
#include <iostream>

class Timer
{
public:
    Timer();
    void Start();
    void Stop();
    std::string CurrentTime();
    // CurrentTime in milli seconds
    long long CTms() const;
    void Reset();
    static std::string FormatDuration(long long milliseconds);

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_point_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_time_point_;
    bool running_;
};

#endif
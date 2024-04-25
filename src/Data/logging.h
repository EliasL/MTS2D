#ifndef LOGGING_H
#define LOGGING_H
#pragma once

#include <string>

#include "settings.h"
#include "Data/dataExport.h"

#include <chrono>
#include <iostream>
#include <cereal/types/chrono.hpp> // Include Cereal support for std::chrono types

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

    std::chrono::time_point<std::chrono::steady_clock> startTime;
    std::chrono::time_point<std::chrono::steady_clock> endTime;
    bool running;

    template <class Archive>
    void serialize(Archive &ar)
    {
        ar(running,
           startTime,
           endTime);
    }
};

#endif

#ifndef LOGGING_H
#define LOGGING_H
#pragma once

#include <chrono>
#include <iostream>
#include <sstream>
#include <iomanip>                 // For std::setprecision
#include <cereal/types/chrono.hpp> // Include Cereal support for std::chrono types
#include <cereal/access.hpp>

class Timer
{
public:
    Timer();
    void Start();
    void Stop();

    std::chrono::milliseconds RunTime() const;
    std::string RunTimeString() const;
    void Reset();

    static std::string FormatDuration(std::chrono::milliseconds duration);

private:
    std::chrono::time_point<std::chrono::steady_clock> lastCheckpoint;
    std::chrono::milliseconds runtime;
    bool running;

    friend class cereal::access;
    template <class Archive>
    void serialize(Archive &ar)
    {
        ar(running,
           lastCheckpoint,
           runtime);
    }
};

#endif
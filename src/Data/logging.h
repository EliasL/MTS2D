#ifndef LOGGING_H
#define LOGGING_H
#pragma once

#include <string>

#include "settings.h"
#include "Data/dataExport.h"

#include <chrono>
#include <iostream>
#include <boost/serialization/access.hpp>

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

    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive &ar, const unsigned int version);
};

#endif

template <class Archive>
void Timer::serialize(Archive &ar, const unsigned int version)
{
    using namespace std::chrono;
    // Serialize running_ state
    ar & running_;

    // Serialize the time points as long long after converting them from time_point to duration
    if (Archive::is_saving::value)
    {
        auto start_dur = duration_cast<microseconds>(start_time_point_.time_since_epoch()).count();
        auto end_dur = duration_cast<microseconds>(end_time_point_.time_since_epoch()).count();
        ar & start_dur;
        ar & end_dur;
    }
    else
    {
        long long start_dur, end_dur;
        ar & start_dur;
        ar & end_dur;
        start_time_point_ = time_point<high_resolution_clock>(microseconds(start_dur));
        end_time_point_ = time_point<high_resolution_clock>(microseconds(end_dur));
    }
}
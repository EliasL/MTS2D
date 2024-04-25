#include "logging.h"
#include <sstream>
#include <iomanip> // For std::setprecision

Timer::Timer() : running(false) {}

void Timer::Start()
{
    startTime = std::chrono::steady_clock::now();
    running = true;
}

void Timer::Stop()
{
    endTime = std::chrono::steady_clock::now();
    running = false;
}

// CurrentTime in milli seconds
long long Timer::CTms() const
{
    auto currentTimePoint = running ? std::chrono::steady_clock::now() : endTime;
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(currentTimePoint - startTime).count();
    return duration;
}

std::string Timer::CurrentTime()
{
    return FormatDuration(CTms());
}

void Timer::Reset()
{
    startTime = std::chrono::steady_clock::now();
    running = true;
}

std::string Timer::FormatDuration(long long milliseconds)
{
    std::ostringstream stream;
    long long total_seconds = milliseconds / 1000;
    milliseconds %= 1000;
    long long hours = total_seconds / 3600;
    total_seconds %= 3600;
    long long minutes = total_seconds / 60;
    double seconds = total_seconds % 60 + milliseconds / 1000.0;

    bool displayHigherUnits = false; // Used to track whether higher units were displayed

    // Days
    if (hours >= 24)
    {
        long long days = hours / 24;
        hours %= 24;
        stream << days << "d ";
        displayHigherUnits = true;
    }

    // Hours
    if (hours > 0 || displayHigherUnits)
    {
        stream << hours << "h ";
        displayHigherUnits = true;
    }

    // Minutes
    if (minutes > 0 || displayHigherUnits)
    {
        stream << minutes << "m ";
    }

    // Seconds with three decimal places
    // Save the current format state of the stream
    std::streamsize prec = stream.precision();
    std::ios_base::fmtflags f(stream.flags());

    stream << std::fixed << std::setprecision(3) << seconds << "s";
    // Save the current format state of the stream

    // Restore the saved precision state
    stream.precision(prec);
    stream.flags(f);
    return stream.str();
}

#include "logging.h"

Timer::Timer() : runtime(0), running(false) {}

void Timer::Start()
{
    if (!running)
    {
        lastCheckpoint = std::chrono::steady_clock::now();
        running = true;
    }
}

void Timer::Stop()
{
    if (running)
    {
        auto now = std::chrono::steady_clock::now();
        runtime += std::chrono::duration_cast<std::chrono::milliseconds>(now - lastCheckpoint);
        running = false;
    }
}

std::string Timer::RunTimeString() const
{
    return FormatDuration(RunTime());
}

std::chrono::milliseconds Timer::RunTime() const
{
    if (running)
    {
        auto now = std::chrono::steady_clock::now();
        auto currentDuration = std::chrono::duration_cast<std::chrono::milliseconds>(now - lastCheckpoint);
        return (runtime + currentDuration);
    }
    return runtime;
}

void Timer::Reset()
{
    runtime = std::chrono::milliseconds(0);
    running = false;
}

#include <chrono>
#include <sstream>
#include <iomanip>

std::string Timer::FormatDuration(std::chrono::milliseconds duration)
{
    std::ostringstream stream;
    auto hours = std::chrono::duration_cast<std::chrono::hours>(duration);
    duration -= hours;
    auto minutes = std::chrono::duration_cast<std::chrono::minutes>(duration);
    duration -= minutes;
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration);
    auto milliseconds = duration.count() % 1000; // Remainder in milliseconds

    bool displayHigherUnits = false;

    if (hours.count() >= 24)
    {
        auto days = hours.count() / 24;
        hours = std::chrono::hours(hours.count() % 24); // This ensures compatibility without C++20
        stream << days << "d ";
        displayHigherUnits = true;
    }

    if (hours.count() > 0 || displayHigherUnits)
    {
        stream << hours.count() << "h ";
        displayHigherUnits = true;
    }

    if (minutes.count() > 0 || displayHigherUnits)
    {
        stream << minutes.count() << "m ";
    }

    // Save the current format state of the stream
    std::streamsize prec = stream.precision();
    std::ios_base::fmtflags f(stream.flags());

    // Seconds with three decimal places for milliseconds
    stream << std::fixed << std::setprecision(3) << seconds.count() + milliseconds / 1000.0 << "s";

    // Restore the saved precision state
    stream.precision(prec);
    stream.flags(f);
    return stream.str();
}

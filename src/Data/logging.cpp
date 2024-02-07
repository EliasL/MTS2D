#include "logging.h"
#include <sstream>
#include <iomanip> // For std::setprecision

Timer::Timer() : running_(false) {}

void Timer::Start()
{
    start_time_point_ = std::chrono::high_resolution_clock::now();
    running_=true;
}

void Timer::Stop() {
    end_time_point_ = std::chrono::high_resolution_clock::now();
    running_ = false;
}


// CurrentTime in milli seconds
long long Timer::CTms() {
    auto currentTimePoint = running_ ? std::chrono::high_resolution_clock::now() : end_time_point_;
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(currentTimePoint - start_time_point_).count();
    return duration;
}   

std::string Timer::CurrentTime() {
    return FormatDuration(CTms());
}


void Timer::Reset() {
    start_time_point_ = std::chrono::high_resolution_clock::now();
    running_ = true;
}


std::string Timer::FormatDuration(long long milliseconds) {
    std::ostringstream stream;
    long long total_seconds = milliseconds / 1000;
    milliseconds %= 1000;
    long long hours = total_seconds / 3600;
    total_seconds %= 3600;
    long long minutes = total_seconds / 60;
    double seconds = total_seconds % 60 + milliseconds / 1000.0;

    bool displayHigherUnits = false; // Used to track whether higher units were displayed

    // Days
    if (hours >= 24) {
        long long days = hours / 24;
        hours %= 24;
        stream << days << "d ";
        displayHigherUnits = true;
    }

    // Hours
    if (hours > 0 || displayHigherUnits) {
        stream << hours << "h ";
        displayHigherUnits = true;
    }

    // Minutes
    if (minutes > 0 || displayHigherUnits) {
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


void setLogFile(std::string simulationName, std::string dataPath)
{

    std::string path = getOutputPath(simulationName, dataPath);
    std::string filePath = path + simulationName + ".log";
    // We set the logging settings
    auto file_logger = spdlog::basic_logger_mt(LOGNAME, filePath);
    spdlog::set_default_logger(file_logger);
    spdlog::info("Starting simulation.");
}



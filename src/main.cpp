#include <iostream>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include <indicators/cursor_control.hpp>
#include "settings.h"
#include "Simulation/simulation.h"
#include "Simulation/energyPlotting.h"

void init()
{

    // We fix the random seed to get reproducable results
    srand(0);

    // We set the logging settings
    auto file_logger = spdlog::basic_logger_mt(LOGNAME, "test.log");
    spdlog::set_default_logger(file_logger);
    spdlog::info("Starting simulation.");

    // Hide cursor
    indicators::show_console_cursor(false);
}

void exit()
{

    // Show cursor
    indicators::show_console_cursor(true);

    // Close and flush logger
    spdlog::drop(LOGNAME);
}

int main()
{
    init();

    Simulation s = Simulation();
    s.run_simulation();
    // drawPicture();

    exit();
}
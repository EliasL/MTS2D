#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#include "Data/param_parser.h"
#include <filesystem>

// Log settings
// #define ELPP_DISABLE_LOGS
// #define ELPP_DISABLE_DEBUG_LOGS

int main(int argc, char *argv[])
{
    std::filesystem::current_path(MTS2D_SOURCE_DIR);
    doctest::Context context;
    context.applyCommandLine(argc, argv);
    setQuiet(true);
    int res = context.run(); // Run tests
    return res;
}

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

// Log settings
// #define ELPP_DISABLE_LOGS
// #define ELPP_DISABLE_DEBUG_LOGS

int main(int argc, char *argv[])
{
    doctest::Context context;
    context.applyCommandLine(argc, argv);
    int res = context.run(); // Run tests
    return res;
}
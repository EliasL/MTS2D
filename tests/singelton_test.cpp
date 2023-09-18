#include "run/doctest.h" // Include Doctest
#include "../src/singelton.h"

// Include the Singleton class declaration here

TEST_CASE("Singleton Instance") {
    // Get two instances of the Singleton class
    S& instance1 = S::getInstance();
    S& instance2 = S::getInstance();

    // Check that both instances are the same
    CHECK(&instance1 == &instance2);
}

TEST_CASE("Setting Grid Size") {
    S& instance = S::getInstance();
    // Set grid size for the first time
    instance.setGridSize(3, 3);
    // Try to set grid size again, it should throw an exception
    CHECK_THROWS_WITH(
        instance.setGridSize(4, 4),
        "The size of the grid has already been set."
    );
}

// Add more test cases as needed

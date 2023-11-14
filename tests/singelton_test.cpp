#include "run/doctest.h" // Include Doctest
#include "../src/Utility/singelton.h"

// Include the Singleton class declaration here

TEST_CASE("Singleton Instance") {
    // Get two instances of the Singleton class
    Singelton &instance1 = Singelton::getInstance();
    Singelton &instance2 = Singelton::getInstance();

    // Check that both instances are the same
    CHECK(&instance1 == &instance2);
}

TEST_CASE("Setting Surface Size") {
    Singelton &instance = Singelton::getInstance();
    // Set surface size for the first time
    instance.setSurfaceSize(3, 3);
    // Try to set surface size again, it should throw an exception
    CHECK_THROWS_WITH(
        instance.setSurfaceSize(4, 4),
        "The size of the surface has already been set."
    );
}

// Add more test cases as needed

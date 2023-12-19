#ifndef SETTINGS_H
#define SETTINGS_H
#pragma once
/*
        SETTINGS
*/

// When you run the simulation, the output files will be stored in this folder
#define OUTPUTFOLDERPATH "output/"

// It can be usefull to have a subfolder for different dates/modes/simulations
// This is the default subfolder, but you can assign a different folder in the
// function argument instead of changing this values.
#define DEFAULTSUBFOLDER "testing/"

// Number of threads. Must be between 1 and nr of cpus on machine.
#define NUMEROFTHREADS 8

#endif
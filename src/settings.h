#ifndef SETTINGS_H
#define SETTINGS_H
#pragma once


/*
        SIMULATION SETTINGS
*/
// Number of threads. Must be between 1 and nr of cpus on machine.
#define NUMEROFTHREADS 8




/*
        FOLDER SETTINGS

        Example layout
        output/
        ├── shearSimulation/
        │   ├── data/
        │   ├── frames/
        │   ├── animation.mp4
        │   └── energy_plot.pdf
*/
// When you run the simulation, the output files will be stored in this folder
#define OUTPUTFOLDERPATH "output/"

// It can be usefull to have a subfolder for different dates/modes/simulations
// This is the default subfolder, but you can assign a different folder in the
// function argument instead of changing this values.
#define SUBFOLDERPATH "testing/"

// Each subfolder will have two raw data folders, one for .vtu, and one for .png
#define DATAFOLDERPATH "data/"
#define FRAMEFOLDERPATH "frames/" 

#endif
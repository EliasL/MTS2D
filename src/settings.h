#ifndef SETTINGS_H
#define SETTINGS_H
#pragma once
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
// A collection file to easily find all relevant vtu files
#define COLLECTIONNAME "collection"
// A small cvs file where each line holds procceced data about a frame
#define MACRODATAFILE "macroData"

// LOG SETTINGS
#define LOGNAME "infoLog"

#endif
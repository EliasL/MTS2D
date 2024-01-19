#ifndef SETTINGS_H
#define SETTINGS_H
#pragma once

/*
        FOLDER SETTINGS

        Example layout
        output/
        ├── shearSimulation/             # simulationName
        │   ├── data/
        │   ├── frames/
        │   ├── animation.mp4
        │   └── energy_plot.pdf
*/
// Remember to include a trailing '/' in the path variables
// When you run the simulation, the output files will be stored in this folder
#define OUTPUTFOLDERPATH "/media/elias/Data/Output/" // Local
//#define OUTPUTFOLDERPATH "/media/elias/Data/Output/" // Cluster

// Each subfolder will have two raw data folders, one for .vtu, and one for .png
#define DATAFOLDERPATH "data/"
#define FRAMEFOLDERPATH "frames/"
// A collection file to easily find all relevant vtu files
#define COLLECTIONNAME "collection"
// A small cvs file where each line holds procceced data about a frame
#define MACRODATAFILE "macroData"

// LOG SETTINGS
// The log name is insignificant and does not determine any output file name
#define LOGNAME "infoLog"

#endif
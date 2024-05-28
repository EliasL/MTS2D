#ifndef SETTINGS_H
#define SETTINGS_H
#pragma once

/*
        FOLDER SETTINGS

        Example layout
        output/                         # dataPath
        ├── shearSimulation/            # simulationName / folderName
        │   ├── data/
        │   ├── frames/
        │   ├── snapshots/
        │   ├── animation.mp4
        │   └── energy_plot.pdf
*/

// Each subfolder will have two raw data folders, one for .vtu, and one for .png
#define OUTPUTFOLDERPATH "MTS2D_output/"
#define DATAFOLDERPATH "data/"
#define FRAMEFOLDERPATH "frames/"
#define DUMPFOLDERPATH "dumps/"
#define BACKUPFOLDERPATH "backups/"
// A collection file to easily find all relevant vtu files
#define COLLECTIONNAME "collection"
// A small cvs file where each line holds procceced data about a frame
#define MACRODATANAME "macroData"
// The name of the config file that is saved at the start of the simulation
#define CONFIGNAME "config.conf"
#endif
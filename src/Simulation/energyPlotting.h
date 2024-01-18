
#ifndef ENERGYPLOTTING_H
#define ENERGYPLOTTING_H
#pragma once

#include "settings.h"
#include "Mesh/tElement.h"
#include "Matrix/matrix.h"
#include "Data/dataExport.h"

// exports image data in the form of a csv file.
void drawPicture(std::string simulationName, int resolution = 500)
{
    /*

    The idea here will be to fill a grid of x, y coordinates with the energy
    associated with that point. We first need to find what sort of element is
    associated with each x,y point. We assume that x,y has been calculated using

    (Homogeneous nucleation of dislocations as a pattern formation phenomenon, page 5)
    a = c12/c22, b = 1/c22,
    x = (a^2 + b^2 - 1)/(a^2 + (b + 1)^2),
    y = (2a)/(a^2 + (b + 1)^2)

    which solves to (almost always, see wolfram link)

    https://www.wolframalpha.com/input?i=solve+x+%3D+%28a%5E2+%2B+b%5E2+-+1%29%2F%28a%5E2+%2B+%28b+%2B+1%29%5E2%29%2C+y+%3D+%282a%29%2F%28a%5E2+%2B+%28b+%2B+1%29%5E2%29+for+a%2C+b
    a = (2 y)/(x^2 - 2 x + y^2 + 1),
    b = -(x^2 + y^2 - 1)/(x^2 - 2 x + y^2 + 1)

    finally,
    c12 = a/b, c22 = 1/b
    and using det(C)=1, we find
    c11 = 1/c12 + c22

    now we use the calculate energy function to find the energy,
    and then export the data to a csv file using the export function.
    */

    // Construct the full file path
    std::string filePath = getDataPath(simulationName) + "energy_grid.csv";
    std::ofstream outputFile(filePath);

    // Define the range for x and y based on the unit circle
    double radius = 1.0;
    double minX = -radius, maxX = radius;
    double minY = -radius, maxY = radius;
    double step = 2 * radius / resolution;

    for (double x = minX; x <= maxX; x += step)
    {
        for (double y = minY; y <= maxY; y += step)
        {

            // Skip points outside the unit circle
            if (x * x + y * y > 0.999999999999)
            {
                outputFile << x << "," << y << ","
                           << "nan" //<< std::endl;
                           << "," << 0 << "," << 0 << "," << 0 << std::endl;
                continue;
            }

            // Calculate a and b from x and y
            double a = (2 * y) / (x * x - 2 * x + y * y + 1);
            double b = -(x * x + y * y - 1) / (x * x - 2 * x + y * y + 1);

            // Matrix2x2<double> F = {{1+x,0+y},{0,1}};
            // Matrix2x2<double> C = F.transpose()*F;
            // C = C * (1/sqrt(C.det()));

            // // Calculate c12, c22, c11
            // double c12 = C[0][1];
            // double c22 = C[1][2];
            // double c11 = C[0][0];//(1+c12*c12)/c22;

            double c12 = a / b;
            double c22 = 1 / b;
            double c11 = (1 + c12 * c12) / c22;

            // Calculate the energy at this point
            double energy = TElement::calculateEnergy(c11, c22, c12);

            // Export the data to CSV
            outputFile << x << "," << y << "," << energy //<< std::endl;
                       << "," << c11 << "," << c22 << "," << c12 << std::endl;
        }
    }

    outputFile.close();
}
#endif
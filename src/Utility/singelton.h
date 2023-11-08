#pragma once
#ifndef SingeltonINGELTON_H
#define SingeltonINGELTON_H

#include "../Grid/grid2D.h"


// https://stackoverflow.com/questions/1008019/how-do-you-implement-the-singleton-design-pattern
class Singelton
{
public:
    static Singelton& getInstance()
    {
        static Singelton instance; // Guaranteed to be destroyed.
                                // Instantiated on first use.
        return instance;
    }
private:
    Singelton(){}

    bool _grid_has_been_given_size = false;

public:
    Singelton(Singelton const&) = delete;
    void operator=(Singelton const&)  = delete;

    Grid g;
    double zeroing_energy; // Used to calculate the energy. 
    void setGridSize(int n, int m, double a=1){
        if(_grid_has_been_given_size){
            throw std::logic_error("The size of the grid has already been set.");
        } else {
            g = Grid(n,m,a);
            _grid_has_been_given_size = true;
        }
    }
};

#endif
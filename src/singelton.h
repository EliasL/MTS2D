#pragma once
#ifndef SINGELTON_H
#define SINGELTON_H

#include "grid2D.h"


// https://stackoverflow.com/questions/1008019/how-do-you-implement-the-singleton-design-pattern
class S
{
public:
    static S& getInstance()
    {
        static S    instance; // Guaranteed to be destroyed.
                                // Instantiated on first use.
        return instance;
    }
private:
    S() {}                   // Constructor? (the {} brackets) are needed here.

    bool _grid_has_been_given_size = false;

public:
    S(S const&)               = delete;
    void operator=(S const&)  = delete;

    Grid g;
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
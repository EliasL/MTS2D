#pragma once
#ifndef SINGELTON_H
#define SINGELTON_H

#include "../Surface/surface.h"


// https://stackoverflow.com/questions/1008019/how-do-you-implement-the-singleton-design-pattern
class Singelton
{
public:
    static Singelton &getInstance()
    {
        static Singelton instance; // Guaranteed to be destroyed.
                                // Instantiated on first use.
        return instance;
    }
private:
    Singelton(){}

    bool _surface_has_been_given_size = false;

public:
    Singelton(Singelton const&) = delete;
    void operator=(Singelton const&)  = delete;

    Surface g;
    double zeroing_energy; // Used to calculate the energy. 
    void setSurfaceSize(int n, int m, double a=1){
        if(_surface_has_been_given_size){
            throw std::logic_error("The size of the surface has already been set.");
        } else {
            g = Surface(n,m,a);
            _surface_has_been_given_size = true;
        }
    }
};

#endif
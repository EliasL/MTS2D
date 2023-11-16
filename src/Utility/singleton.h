#pragma once
#ifndef SINGELTON_H
#define SINGELTON_H

#include "Mesh/mesh.h"


// https://stackoverflow.com/questions/1008019/how-do-you-implement-the-singleton-design-pattern
class Singleton
{
public:
    static Singleton &getInstance()
    {
        static Singleton instance; // Guaranteed to be destroyed.
                                // Instantiated on first use.
        return instance;
    }
private:
    Singleton(){}

    bool _surface_has_been_given_size = false;

public:
    Singleton(Singleton const&) = delete;
    void operator=(Singleton const&)  = delete;

    Mesh mesh;
    double zeroing_energy; // Used to calculate the energy. 
    void setSurfaceSize(int n, int m, double a=1){
        if(_surface_has_been_given_size){
            throw std::logic_error("The size of the surface has already been set.");
        } else {
            mesh = Mesh(n,m,a);
            _surface_has_been_given_size = true;
        }
    }
};

#endif
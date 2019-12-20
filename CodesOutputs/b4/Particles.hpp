//
//  Particles.hpp
//  
//
//  Created by Miranda Holmes-Cerfon on 2019-07-18.
//

#ifndef Particles_hpp
#define Particles_hpp

#include <stdio.h>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>



using namespace std;
using namespace Eigen;


class Particles {
  
public:
    // Constructor
    Particles(int n, int dim, VectorXd params0);
    
    // Must set all the things below
    void setup(void);   // initialize the variables below
    
    // -------   Functions that change with system  ------ //
    // Particle parameters (initialized in setup
    int n;     // number of particles
    int dim;   // dimension
    VectorXd x0;  // initial condition
    VectorXd params;  // various parameters we might wish to set in script
    VectorXd rad;  // radii of particles
    double rho;   // (inverse) width parameter for interactions
    VectorXd KapParams; // sticky parameters
    VectorXd EParams;   // energy depth parameters
    double Ewca;   // repulsive parameters
    int mwca;     
    double kspring;   // spring parameter (spring constant)
    
    // Energy & gradient
    double U(VectorXd& x);  // energy
    void dU(VectorXd& x, VectorXd& df);  // gradient of energy
    
    // Data to save along the way
    int nstats;  // number of different quantities to compute
    VectorXd computestats(VectorXd& x);
    
    
    // -------   Fixed functions; don't change with system  ------ //
    double dist(VectorXd& x, int i, int j);  // distance between particles i,j
    // Specific energy functions
    double Morse(double r, double r0, double E, double rho);
    double dMorse(double r, double r0, double E, double rho);
    double Spring(double r, double r0, double k);  // spring potential for permanent interactions
    double dSpring(double r, double r0, double k);
    double HalfSpring(double r, double r0, double k);  // for repulsive interactions
    double dHalfSpring(double r, double r0, double k);
    double WCA(double r, double r0, double E, int m);  // WCA potential for repulsive interactions
    double dWCA(double r, double r0, double E, int m);
    
    // Convert & compute parameters
    double E_morse(double kap, double rho, double e0=7);  // computes E parameter, from sticky etc
    double Sticky_morse(double E, double rho);  // computes sticky parameter

    // Other useful math functions
    double costh(VectorXd& x, int i, int j, int k);  // cos angle between i,j,k
    
    // mathematical constants
    constexpr static const double pi = 3.141592653589793238462643383279502884;
    constexpr static const double sqrtpi = 1.77245385090551602729816748334115;
    constexpr static const double sqrt2pi = 2.50662827463100050241576528481105;

};

#endif /* Particles_hpp */

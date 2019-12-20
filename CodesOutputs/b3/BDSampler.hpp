//
//  BDSampler.hpp
//  
//
//  Created by Miranda Holmes-Cerfon on 2019-07-17.
//

#ifndef BDSampler_hpp
#define BDSampler_hpp

#include <stdio.h>
#include <iostream>
#include <iomanip>           // for setw
#include <fstream>           // for writing to files
#include <cmath>
#include <ctime>             // for timing operations
#include <Eigen/Dense>
#include <random>
#include "Particles.hpp"


using namespace std;
using namespace Eigen;


class BDSampler {
private:
    // Random number generators and variables
    mt19937                       mrng;        // defines the generator
    int                           mseed;       // seed for rng
    normal_distribution<double>   mZdist;      // standard random numbers
    void         initializeRandom(void);       // initialize random number generators
    
    // Info about particles
    Particles* mp;         // holds particle system information
    MatrixXd mPts;         // list of points sampled
    ArrayXXd mData;        // holds derived statistics
    int mnsave;            // how many data points are saved
    double mtotalTime;     // how many points generated
    double mdisptime;      // how often to display time
    
    // Clock
    double   mtime;            // time it took to sample (in seconds)
    
    
public:
    // Constructor
    BDSampler(Particles *parts,
              int seed=0);
    
    
    // sampling function
    int sampleEM(double dt, double totaltime, int dsave=1);  // sampling algorithm

    // Return internal variables
    const MatrixXd* pts()   const { return &mPts; }
    const ArrayXXd* data()   const { return &mData; }
    ArrayXd stats(void){ return mData.colwise().mean();}  // return vector of statistics
    static ArrayXd std(ArrayXXd& a) { return (a.square().colwise().sum()/a.rows() - a.colwise().sum().square()/(a.rows()*a.rows())).sqrt(); }
    ArrayXd statsstd(void) { return std(mData); }
    int write_pts(string ptsfile, int dwrite=1);  // write points and/or flags to file
    int write_data(string datafile, int dwrite=1);  // write data to file
    int prec = 8;    // precision for writing output to file
    
    
    // Generate random numbers and set internal variables
    double normRand(void) { return mZdist(mrng);}  // output a standard normal random number
    void setSeed  (int    seed)  { mrng.seed(seed); mseed = seed; }
    int  seed(void) const  { return mseed;  }
    double  time(void) const  { return mtime;  }
    void printtime(void);
    void setdisptime(double);
};

#endif /* BDSampler_hpp */

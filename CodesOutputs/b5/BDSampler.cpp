//
//  BDSampler.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 2019-07-17.
//

#include "BDSampler.hpp"


// Constructor
BDSampler::BDSampler(Particles *parts,
                     int seed)
: mp(parts), mseed(seed) {
    initializeRandom();     // initialize random number generators
    setdisptime(1e2);
}


// Initialize random number generator
void  BDSampler::BDSampler::initializeRandom(void) {
    if(mseed==0) {
        mseed = random_device()(); // choose seed randomly based on device configurations
    }
    mrng.seed(mseed);    // initialize generator for random number stream
}


// Set how often to display time
void BDSampler::setdisptime(double dtime) {
    mdisptime = dtime;
}


// Sampling algorithm
int BDSampler::sampleEM(double dt, double totalTime, int dsave) {
    
    mtotalTime = totalTime;
    
    int n = mp->n;     // number of particles
    int dim = mp->dim;     // dimension they're in
    VectorXd x = mp->x0;   // configuration; set to initial condition
    VectorXd dx(n*dim);   // holds change in configuration
    double F;    // holds force
    VectorXd gradF(n*dim);  // gradient of force
    VectorXd noise(n*dim);  // noise vector
    clock_t start_time, end_time;     // start,end times
    
    // Initialize matrices to hold data
    mnsave = floor(totalTime / dt / dsave) + 1;  // prev: floor(npts / dsave) + 1;
    mPts   = MatrixXd::Zero(mnsave,n*dim);
    mData  = MatrixXd::Zero(mnsave,mp->nstats);
    mPts.row(0) = x; // save data initially
    mData.row(0) = mp->computestats(x);
    
    
    // *************  Loop through MCMC steps  ***************
    start_time = clock();            // start timing now
    double currentTime = 0;
    double dispTime = 0;
    int ipt = 1;
    int ctr = 0;
    while(currentTime < totalTime){
        
        currentTime = currentTime + dt;
        dispTime = dispTime + dt;
        ctr++;   // index for saving points

        // Get random noise
        for(int j=0; j<n*dim; j++) {
            noise(j) = normRand();
        }
        
        // Update configuration
        mp->dU(x,gradF);
        dx = -gradF*dt + sqrt(2*dt)*noise;
        x = x+dx;
        
        // save x, statistics
        if(ctr == dsave && ipt < mnsave) {
            mPts.row(ipt) = x;
            mData.row(ipt) = mp->computestats(x);
            ipt++;
            ctr = 0;
        }
        
        // display output, every 100 simulation seconds
        if(dispTime >= mdisptime) {
            cout << "  currentTime = " << currentTime << ", simulation time = " << (clock() - start_time ) / (double) CLOCKS_PER_SEC / 60. << " minutes." << endl;
            dispTime = 0;
        }
        
    }
    // save time
    mtime = (clock() - start_time ) / (double) CLOCKS_PER_SEC;
    
    return 0;
}



// Write points to file, at intervals of dwrite points
int BDSampler::write_pts(string ptsfile,
                           int dwrite) {
    
    ofstream myfile1 (ptsfile); // SHOULD check for problems
    myfile1.precision(prec);   // Set precision
    
    // Write points to file
    int i=0;
    while(i<mnsave) {
        myfile1 << mPts.row(i) << endl;
        i += dwrite;
    }
    myfile1.close();
    
    // return
    return 0;
}

// Write data to file, at intervals of dwrite points
int BDSampler::write_data(string datafile,
                         int dwrite) {
    
    ofstream myfile1 (datafile); // SHOULD check for problems
    myfile1.precision(prec);   // Set precision
    
    // Write data to file
    int i=0;
    while(i<mnsave) {
        myfile1 << mData.row(i) << endl;
        i += dwrite;
    }
    myfile1.close();
    
    // return
    return 0;
}


void BDSampler::printtime(void) {
    cout << "-----  Timing  -----" << endl;
    cout << "Total time    : " << mtime << " seconds" << endl;
    cout << "You can simulate 1000 time units in  " << 1000*mtime/mtotalTime  << " seconds"<< endl;
}

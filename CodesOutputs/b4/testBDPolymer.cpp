//
//  testWormChain.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 2019-07-18.
//

#include <stdio.h>
#include <Eigen/Dense>
#include <chrono>   // for high-precision timing
#include <cmath>

#include "Particles.hpp"
#include "BDSampler.hpp"


typedef std::chrono::high_resolution_clock Clock;

using namespace Eigen;
using namespace std;



// Main loop
int main(int argc,char *argv[])
{
    
    /*  ------  Parameters  ------ */
    int n = 6;
    int dim = 3;
    double kappa = 4.;  // sticky parameter for interactions
    VectorXd params(1);
    params << kappa;
    
    // filename
    string datafile = "polymerBD_b4_data.txt";
    string ptsfile = "polymerBD_pts.txt";
    
    // Sampling parameters
    double dt = 1e-6;         // time step
    double totalTime = 1e4;   // how long to simulate for
    int dsave = floor(0.05/dt); //1000*10 * 5;  // how often to save points
    int seed = 0;     // = 0 for random seed (based on clock time)
    

    /*  ------  Begin code  ------ */
    
    // Set up objects
    Particles p(n,dim,params);
    BDSampler mysampler(&p,seed);
    mysampler.setdisptime(500.);
    cout << "seed = " << mysampler.seed() << endl;
    
    
    // Sample!
    mysampler.sampleEM(dt, totalTime, dsave);
    
    
    
    // check for numerical instability, by computing energy
    //const MatrixXd* data = mysampler.data();
    //cout << "Mean energy: " << (*data).col(0).sum()/(double)npts << endl;
    
    
    /*  ------  Print Output  ------ */
    // Print out statistics
    ArrayXd stats = mysampler.stats();
    //cout << "Mean energy: " << stats(0) << endl;
    cout << "Statistics: \n           " << stats.transpose() << endl;
    cout << "Std dev: \n              " << mysampler.statsstd().transpose() << endl;
    
    
    mysampler.printtime();
    
    
    // --------   Write Points & Data to File  --------
    
    // Save points and flags to file
    cout << "writing data to file " << endl;
    auto t1 = Clock::now();
    //mysampler.write_pts(ptsfile);
    mysampler.write_data(datafile);
    auto t2 = Clock::now();
    cout << "Writing to file took "
    << chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9
    << " seconds" << std::endl;
    
    

    
    return 0;
}






// DEBUG
/* // each potential energy function & its gradient
 double r = 1.2;
 double r0 = 1;
 double v = 1e-8*mysampler.normRand();
 double E = 400;
 double rho = 30;
 int m = 12;
 double du1 = (p.Spring(r+v,r0,E) - p.Spring(r,r0,E))/v;
 double du2 = p.dSpring(r,r0,E);
 cout << "du1 = " << du1 << ", du2 = " << du2 << ", diff = " << du1-du2 << endl;
 */

/*
 VectorXd v(p.n*p.dim), df(p.n*p.dim), w(p.n*p.dim);
 double sig = 1e-7;
 double sig2 = 1;
 for(int i=0;i<p.n*p.dim;i++) {
 v(i) = sig*mysampler.normRand();
 w(i) = sig2*mysampler.normRand();
 }
 VectorXd x0 = p.x0 + w;
 VectorXd xv = x0+v;
 double du1 = (p.U(xv) - p.U(x0));
 p.dU(x0,df);
 double du2=df.transpose()*v;
 cout << "du1 = " << du1 << ", du2 = " << du2 << ", diff = " << du1-du2 << endl;
 */





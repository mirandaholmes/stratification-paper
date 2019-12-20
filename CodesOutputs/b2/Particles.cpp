//
//  Particles.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 2019-07-18.
//

#include "Particles.hpp"


// Constructor
Particles::Particles(int n0, int dim0, VectorXd params0)
: n(n0), dim(dim0), params(params0)
{
    setup();    // initialize variables
}


// ===========  Standard, keep in every code  =========== //

// distance between particles i,j
double Particles::dist(VectorXd& x, int i, int j){
    double d = 0;
    for(int k=0;k<dim; k++) {
        d += (x(dim*i+k) - x(dim*j+k)) * (x(dim*i+k) - x(dim*j+k)) ;
    }
    return sqrt(d);
}

// cosine of angle between particles i,j,k
double Particles::costh(VectorXd& x, int i, int j, int k) {
    VectorXd p1(dim);
    VectorXd p2(dim);
    for (int c=0; c<dim; c++) {
        p1(c) = x(dim*i+c) - x(dim*j+c);
        p2(c) = x(dim*k+c) - x(dim*j+c);
    }
    return p1.dot(p2)/p1.norm()/p2.norm();
}

// Morse potential
double Particles::Morse(double r, double r0, double E, double rho) {
    return E*(1.-exp(-rho*(r-r0)))*(1.-exp(-rho*(r-r0))) - E;
}

// gradient of Morse potential
double Particles::dMorse(double r, double r0, double E, double rho) {
    return 2*rho*E*(exp(-rho*(r-r0)) - exp(-2*rho*(r-r0)) );
}
// 2nd deriv is 2*rho^2*E*(-exp(-rho*(r-r0)) + 2*exp(-2*rho*(r-r0))).

// spring potential
double Particles::Spring(double r, double r0, double k){
    return 0.5*k*(r-r0)*(r-r0);
}
double Particles::dSpring(double r, double r0, double k){
    return k*(r-r0);
}

// Spring potential; repulsive part only
double Particles::HalfSpring(double r, double r0, double k) {
    if(r > r0) return 0.;
    else return Spring(r, r0, k);
}

double Particles::dHalfSpring(double r, double r0, double k) {
    if(r > r0) return 0.;
    else return dSpring(r, r0, k);
}

// WCA potential for repulsive interactions
double Particles::WCA(double r, double r0, double E, int m) {
    if(r > r0) return 0.;
    return E*( pow(r0/r,2*m) - 2*pow(r0/r,m) ) + E;
}

double Particles::dWCA(double r, double r0, double E, int m) {
    if(r > r0) return 0.;
    return -2*m*E/r*( pow(r0/r,2*m) - pow(r0/r,m) );
}



// Compute kappa from E, rho
double Particles::Sticky_morse(double E, double rho) {
    return (sqrtpi / rho) * exp(E) / sqrt(E);
}


// Compute E from kappa, rho
double Particles::E_morse(double k0, double rho, double e0) {
    
    int maxiter = 15;
    double tol = 1e-5;
    
    double e = e0;  // initial guess
    double de, fp, k; // change in e, d(sticky)/de, kappa
    
    // run newton to get a better guess
    for(int i=0; i<maxiter; i++) {
        k = Sticky_morse(e,rho);
        if(fabs(k - k0) < tol) {
            return e;
        }
        fp = k * (1-0.5/e);
        de = -(k-k0) / fp;
        e = e+de;
        
        /* //debug
        cout << "i = " << i << endl;
        cout << "k = " << k << endl;
        cout << "e = " << e << endl;
        cout << "de = " << de << endl;
         */ 
    }
    return NAN;
}

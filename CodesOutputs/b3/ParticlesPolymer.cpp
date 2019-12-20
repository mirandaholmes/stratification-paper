//
//  Particles_WormChain.cpp
//  
//
//  Created by Miranda Holmes-Cerfon on 2019-07-18.
//

#include "Particles.hpp"



// Initialize all the variables that depend on particular system.
void Particles::setup(void) {
    
    nstats = 2;

    rad = 0.5*VectorXd::Ones(n*dim); // radii (all ones)

    // potential energy parameters
    double kap_reg = params(0);       // sticky parameter for 1-n interaction
    rho = 60;       // width of interaction
    kspring = 4*rho*rho * 1.5;  // spring constant for backbone interaction
    double E_reg = E_morse(kap_reg,rho, 7.);
    KapParams.resize(1);
    EParams.resize(1);
    KapParams << kap_reg;
    EParams << E_reg;
    

    // Initial condition: Straight or Bent line
    x0 = VectorXd::Zero(n*dim);
    for(int i=0; i<n; i++) {     // straight line
        x0(i*dim) = i-(n-2.)/2.;
    }
    x0((n-1)*dim) = n-2.-(n-2.)/2.; // bent at end
    x0((n-1)*dim+1) = 1;
    
    
    //debug
    cout << "kappa (approx) = " << kap_reg << ", E = "<< E_reg << endl;
    
}

VectorXd Particles::computestats(VectorXd& x) {
    VectorXd data(nstats);
    // Energy
    //data(0) = U(x);
    
    // number of equations
    double rcut = 2.5;  // cutoff for interactions (diamter+rcut/rho)
    int neqns = n-1;  // start with backbone
    for(int ii=0;ii<n-2;ii++) {
        for(int jj=ii+2;jj<n;jj++) {
            if(dist(x,ii,jj) < rad(ii)+rad(jj)+rcut/rho) neqns++;
        }
    }
    data(0) = neqns;
    return data;
}


double Particles::U(VectorXd& x) {
    
    double e = 0.;
    double r,r0;
    
    for(int i=0; i< n-1; i++) {
        for (int j=i+1; j< n; j++){
            r = dist(x,i,j);
            r0 = rad(i)+rad(j);
            if(j == i+1) {
                e += Spring(r,r0,kspring);
            } else  {
                e += Morse(r,r0,EParams(0), rho);
            }
        }
    }
    return e;
}

void Particles::dU(VectorXd& x, VectorXd& dF) {
    
    dF = VectorXd::Zero(n*dim);
    double up;
    double r,r0;
    
    for(int i=0; i< n-1; i++) {
        for (int j=i+1; j< n; j++){
            r = dist(x,i,j);
            r0 = rad(i)+rad(j);
            if(j == i+1) {
                up = dSpring(r,r0,kspring);
            } else {
                up = dMorse(r,r0,EParams(0), rho);
            }
            for(int k=0;k<dim;k++) {
                dF(dim*i + k) += up/r * (x(dim*i+k)-x(dim*j+k));
                dF(dim*j + k) += -up/r * (x(dim*i+k)-x(dim*j+k));
            }
        }
    }
}


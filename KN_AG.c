//
//  intro.c
//  intro
//
//  Created by Adithan Kathirgamaraju on 12/3/18.
//  Copyright © 2018 Adithan Kathirgamaraju. All rights reserved.
//

#include "intro.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <sys/time.h>

#define imax  60
#define jmax  60
#define smax  50  // Number of timesteps between tmin and tmax over which lightcurve is plotted
double hypg2f1 (double a, double b, double c, double x)
{
    double y;
    if (x>=-1 && x<=1) {
        y = gsl_sf_hyperg_2F1 (a,b,c,x);
        
    }
    else if (x<-1){
        double small=(c-b)/100000;
        y=pow((1-x),-a)*gsl_sf_hyperg_2F1(a, c-b+small, c, x/(x-1));
    }
    /*else if (x){
     y=pow(1-x, -a)*((gsl_sf_gamma(c)*gsl_sf_gamma(b-a))/(gsl_sf_gamma(b)*gsl_sf_gamma(c-a)))*gsl_sf_hyperg_2F1(a, c-b, a-b+1, 1/(1-x))+pow(1-x, -b)*((gsl_sf_gamma(c)*gsl_sf_gamma(a-b))/(gsl_sf_gamma(a)*gsl_sf_gamma(c-b)))*gsl_sf_hyperg_2F1(b, c-a, b-a+1, 1/(1-x));}*/
    else{
        y=pow(x, -a)*((gsl_sf_gamma(c)*gsl_sf_gamma(c-b-a))/(gsl_sf_gamma(c-b)*gsl_sf_gamma(c-a)))*gsl_sf_hyperg_2F1(a, a-c+1, a+b-c+1, 1-1/(x))+pow(1-x, c-a-b)*pow(x,a-c)*((gsl_sf_gamma(c)*gsl_sf_gamma(a+b-c))/(gsl_sf_gamma(a)*gsl_sf_gamma(b)))*gsl_sf_hyperg_2F1(c-a, 1-a, c-b-a+1, 1-1/(x));}
    
    return y;
}

double const thetaobs=0.5236;
double const tmin=8.; // The minimum time of the light curve 10^tmin
double const tmax=10.; // The max time
double const rmin= 16.;
double const rmax= 20.; // max and min radius to interpolate within 10^(rmin) to 10^(rmax)
double const clight=3.e10;
double const sigmat=6.65e-25;
double const mp=1.67e-24;
double const me=9.11e-28;
double const e=4.8e-10;
double const beta0=0.35;  //min velocity of KN blastwave
double const theta0=15.*M_PI/180;  // polar extent
double const thetamax=60.*M_PI/180;// of blast wave
double const thetacore=thetamax-theta0;
double const thetadiv=thetacore/imax;
double const phidiv=2.*M_PI/jmax;
double const E=1e51; //true energy
double const n=0.001;
double const alpha=5;
double const eps_e=0.1;
double const eps_b=0.001;
double const p=2.05;
double const nu=3.0e9;   //observing frequency (Hz)
double const Gamma_DN=1+(0.0005455*(p-1))/(eps_e*(p-2));
double const distance=1.23e26;

double theta[imax];
double phi[jmax];
double Tobs[smax+1];
double r_0[smax+1];
double T_array[smax+1][imax][jmax];
int l[smax+1][imax][jmax];
double r_obs[smax+1][imax][jmax];
double Fnu_obs[smax+1][imax][jmax];
double E_obs[smax+1][imax][jmax];


void initialize(){
    
    int i,j,s;
    
    for(i = 0; i <= imax-1; i++){
        theta[i]=theta0+(i+0.5)*thetadiv;
    }
    for(j = 0; j <= jmax-1; j++){
        phi[j]=(j+0.5)*phidiv;
    }
    for(s = 0; s <= smax; s++){
        Tobs[s]=0.3*pow(10,tmin+s*(tmax-tmin)/smax);
        r_0[s]=pow(10,rmin+s*(rmax-rmin)/smax);
        //Tobs[s]=pow(10,tmin)+(s/smax)*(pow(10, tmax)-pow(10, tmin));
        //r_0[s]=pow(10,rmin)+(s/smax)*(pow(10, rmax)-pow(10, rmin));
    }
}


double E_iso (void){
    double result;
    result=2*E/(cos(theta0)-cos(thetamax));
    return result;
}

double g0 (void){
    double result;
    result=1/(sqrt(1-pow(beta0,2.)));
    return result;
}
double R_dec(){
    double rdec;
    rdec=pow(3*E_iso()/(4*M_PI*mp*clight*clight*g0()*(g0()-1)*n), 0.333333);
    return rdec;
}
double betagamma(double r){
    double bg;
    if (r<R_dec()) {
        bg=beta0*g0()/(pow(r/R_dec(), 3./(2+alpha)));
    }
    else{
        bg=beta0*g0()/(pow(r/R_dec(), 1.5));
    }
        return bg;
    
}
double beta(double r){
    double result;
    result=sqrt(1-1/(pow(betagamma(r),2.)+1));
    return result;
}

double Gamma(double r){
    double result;
    result=sqrt(pow(betagamma(r),2.)+1);
    return result;
}

double Einjected(double r){
    double e;
    if (r<=R_dec()) {
        e=E*(pow((beta0*g0())/betagamma(r), alpha));
    }
    else{
        e=E;
    }
    return e;
}

//lab frame time
double t(double r){
    double result;
    if (r<=R_dec()) {
        result=(r*(2+alpha)*(pow(r/R_dec(), 3./(2+alpha))))/((5.+alpha)*beta0*g0()*clight);
    }
    else{
        result=((R_dec()*(2+alpha))/(5+alpha)+0.4*r*(pow(r/R_dec(), 1.5)-1))/(clight*beta0*g0());
    }
    return result;
}
 
 //Defining hypergeometric2F1 for domain outside the unit circle as well as inside, caution that this function will fail (as defined below) when Gamma function is infinity (at 0,-1,-2,....)
 double hypg2f1 (double a, double b, double c, double x)
 {
 double y;
 if (x>=-1 && x<=1) {
 y = gsl_sf_hyperg_2F1 (a,b,c,x);
 
 }
 else{
 y=pow(1-x, -a)*((gsl_sf_gamma(c)*gsl_sf_gamma(b-a))/(gsl_sf_gamma(b)*gsl_sf_gamma(c-a)))*gsl_sf_hyperg_2F1(a, c-b, a-b+1, 1/(1-x))+pow(1-x, -b)*((gsl_sf_gamma(c)*gsl_sf_gamma(a-b))/(gsl_sf_gamma(a)*gsl_sf_gamma(c-b)))*gsl_sf_hyperg_2F1(b, c-a, b-a+1, 1/(1-x));}
 
 return y;
 }

double doppler(double r, int i, int j){
    double result;
    result= 1/(Gamma(r)*(1-beta(r)*(cos(thetaobs)*cos(theta[i])+sin(thetaobs)*sin(theta[i])*cos(phi[j]))));
    return result;
}

double tobs(double r, int i, int j){
    double result;
    if (r<=R_dec()) {
        result=(-(r*(cos(theta[i])*cos(thetaobs) + cos(phi[j])*sin(theta[i])*sin(thetaobs))) + (r*sqrt((pow(g0(),2)*pow(R_dec(),6./(2 + alpha))*pow(beta0,2.))/
                                                                                                             (pow(r,6./(2 + alpha)) + pow(g0(),2.)*pow(R_dec(),6./(2 + alpha))*pow(beta0,2.)))*
                                                                                                      ((2 + alpha)*(pow(r,6/(2 + alpha)) + pow(g0(),2.)*pow(R_dec(),6/(2 + alpha))*pow(beta0,2.)) +
                                                                                                       3*pow(g0(),2.)*hypg2f1(0.5,(2 + alpha)/6.,(8 + alpha)/6.,-(pow(r,6/(2 + alpha))/(pow(g0(),2.)*pow(R_dec(),6/(2 + alpha))*pow(beta0,2.))))*
                                                                                                       pow(R_dec(),6/(2 + alpha))*pow(beta0,2.)*sqrt((pow(r,6/(2 + alpha)) + pow(g0(),2.)*pow(R_dec(),6/(2 + alpha))*pow(beta0,2.))/
                                                                                                            (pow(g0(),2.)*pow(R_dec(),6/(2 + alpha))*pow(beta0,2.)))))/((5 + alpha)*pow(g0(),2.)*pow(R_dec(),6/(2 + alpha))*pow(beta0,2.)))/clight;
    }
    else{
        result=(r*hypg2f1 (-0.5,0.3333333333333333,1.3333333333333333,-(pow(r,3.)/(pow(g0(),2.)*pow(R_dec(),3.)*pow(beta0,2.)))) -
                hypg2f1 (-0.5,0.3333333333333333,1.3333333333333333,-(1/(pow(g0(),2.)*pow(beta0,2.))))*R_dec() +
                (-r + R_dec())*(cos(theta[i])*cos(thetaobs) + cos(phi[j])*sin(theta[i])*sin(thetaobs)))/clight +
        (-(R_dec()*(cos(theta[i])*cos(thetaobs) + cos(phi[j])*sin(theta[i])*sin(thetaobs))) + (pow(R_dec(),(-4. + alpha)/(2. + alpha))*
                                                                                                              sqrt((pow(g0(),2.)*pow(beta0,2.))/(1 + pow(g0(),2.)*pow(beta0,2.)))*
                                                                                                              ((2 + alpha)*pow(R_dec(),6/(2 + alpha))*(1 + pow(g0(),2.)*pow(beta0,2.)) +
                                                                                                               3*pow(g0(),2.)*hypg2f1 (0.5,(2 + alpha)/6.,(8 + alpha)/6.,-(1/(pow(g0(),2.)*pow(beta0,2.))))*pow(R_dec(),6/(2 + alpha))*pow(beta0,2.)*
                                                                                                               sqrt((1 + pow(g0(),2.)*pow(beta0,2.))/(pow(g0(),2.)*pow(beta0,2.)))))/((5 + alpha)*pow(g0(),2.)*pow(beta0,2.)))/clight;
    }
    return result;
}
/*double F_numax(double r, int i, int j){
    double result;
    result=1/(pow(distance, 2.))*4.775*pow(10, 2.)*pow(doppler(r, i, j), 3.)*((p - 1)/(3*p - 1))*pow((Gamma(r)-1), 0.5)*(pow(Gamma(r), 0.5)*pow((eps_b), 0.5)*pow(n, 1.5))*sin(theta[i]);
    return result;
    
}*/
// nu_m in observer frame
double nu_m(double r, int i, int j){
    double result;
    result=5.5*pow(10, 12.)*doppler(r, i, j)*pow((p - 2)/(p - 1), 2.)*pow((Gamma(r)-1), 2.5)*(pow(Gamma(r), 0.5)*pow((eps_b*n), 0.5)*pow(eps_e, 2.));
    return result;
    
}
//nu_c in observer frame
double nu_c(double r, int i, int j){
    double result;
    result=(1.225*pow(10, 25)*doppler(r, i, j)*pow((Gamma(r)-1), -1.5))/(pow(Gamma(r), 1.5)*pow((eps_b*n), 1.5)*pow(t(r), 2.));
    return result;
    
}

//nu_c with modifications to include synchrotron and self synchrotron compton energy loss (beniamini et al 2015)

double nu_c2(double r, int i, int j){
    double result;
result=nu_c(r, i, j)/pow(1 +
0.5*(pow(1 + (4*eps_e*
               pow((4*pow(10, 12.)*pow(sigmat, 2.)*
                     pow(32*M_PI*mp, 1.5)*pow(clight, 2.)*
                     pow((p - 2)/(p - 1), 2.)*pow(eps_e, 2.)*
                     pow(eps_b, 2.)*pow(n, 2.)*
                     pow(Gamma(r) - 1, 4.)*pow(Gamma(r), 2.))/
                    ((0.286*27.*M_PI*me*e)/pow(t (r),
                                                2.)), (p - 2)/2.))/((3 - p)*eps_b), 0.5) - 1),
                         2);
    return result;
}

// in future more accurate to replace pow(2*pow(10, 6.) with pow(2.34*pow(10, 6.)
double F_nu(double r, int i, int j){
    
    double fnu;
    //if(nu_c < nu_m)
    if(5.5*pow(10, 12.)*doppler(r, i, j)*pow((p - 2)/(p - 1), 2.)*pow((Gamma(r)-1), 2.5)*(pow(Gamma(r), 0.5)*pow((eps_b*n), 0.5)*pow(eps_e, 2.)) > (1.225*pow(10, 25)*doppler(r, i, j)*pow((Gamma(r)-1), -1.5))/(pow(Gamma(r), 1.5)*pow((eps_b*n), 1.5)*pow(t(r), 2.))){
        printf("nu_m>nu_c ordering is not taken into account\n");
        return -1;
        
    }
    //if(nu_m(r, i, j) > nu && Gamma(r) < Gamma_DN )
   else if(5.5*pow(10, 12.)*doppler(r, i, j)*pow((p - 2)/(p - 1), 2.)*pow((Gamma(r)-1), 2.5)*(pow(Gamma(r), 0.5)*pow((eps_b*n), 0.5)*pow(eps_e, 2.)) > nu && Gamma(r) < Gamma_DN ){
        printf("DN phase when nu < nu_m is not accounted for\n");
        return -1;
        
    }
    // if (nu<nu_m)
  else if(nu < 5.5*pow(10, 12.)*doppler(r, i, j)*pow((p - 2)/(p - 1), 2.)*pow((Gamma(r)-1), 2.5)*(pow(Gamma(r), 0.5)*pow((eps_b*n), 0.5)*pow(eps_e, 2.))){
       
       fnu=1/(pow(distance, 2.))*4.775*pow(10, 2.)*pow(doppler(r, i, j), 3.)*((p - 1)/(3*p - 1))*pow((Gamma(r)-1), 0.5)*(pow(Gamma(r), 0.5)*pow((eps_b), 0.5)*pow(n, 1.5))*sin(theta[i])*pow(5.5*pow(10, 12.)*doppler(r, i, j)*pow((p - 2)/(p - 1), 2.)*pow((Gamma(r)-1), 2.5)*(pow(Gamma(r), 0.5)*pow((eps_b*n), 0.5)*pow(eps_e, 2.)),-0.3333)*pow(nu, 0.3333);
       
   }
     // if (nu<nu_c & Gamma>Gamma_DN), nu_c is modified here to take into account synchrotron and self synchrotron compton energy loss (beniamini et al 2015)
   else if (Gamma(r) > Gamma_DN && nu <
        nu_c2(r, i, j))
    {
        fnu =
        1/(pow(distance, 2.)*(3*p - 1))*
        pow(doppler(r, i, j), (5 + p)/2.)*(4.775*pow(10, 2.))*(p - 1)*
        pow(eps_b, (p + 1)/4.)*pow(n, (p + 5)/4.)*
        pow(r, 3.)*pow(Gamma(r), (p + 1)/4.)*
        pow(Gamma(r) - 1, (5*p - 3)/4.)*
        sin(theta[i])*pow(2*pow(10, 6.), p - 1)*
        pow((p - 2)/(p - 1), p - 1)*
        pow(eps_e, p - 1)*pow(nu, (1 - p)/2.);}
    
    
    // if (nu>nu_c & Gamma>Gamma_DN)
    else if (Gamma(r) > Gamma_DN && nu>=nu_c2(r, i, j))
    {
        fnu = (pow(doppler(r, i, j), (6 + p)/
                    2.)*(1.67*pow(10, 15.))*(p - 1)*pow(eps_b, (p - 2)/4.)*
               pow(n, (p + 2)/4.)*pow(r, 3.)*
               pow(Gamma(r), (p - 2)/4.)*
               pow(Gamma(r) - 1, (5*p - 6)/4.)*
               sin(theta[i])*pow(2*pow(10, 6.), p - 1)*
               pow((p - 2)/(p - 1), p - 1)*
               pow(eps_e, p - 1))/(pow(distance, 2.)*(3*p - 1)*(t (r)*pow(nu, p/2.))*
                                    (1 +
                                     0.5*(pow(1 + (4*eps_e*
                                                    pow((4*pow(10, 12.)*pow(sigmat, 2.)*
                                                          pow(32*M_PI*mp, 1.5)*pow(clight, 2.)*
                                                          pow((p - 2)/(p - 1), 2.)*pow(eps_e, 2.)*
                                                          pow(eps_b, 2.)*pow(n, 2.)*
                                                          pow(Gamma(r) - 1, 4.)*pow(Gamma(r), 2.))/
                                                         ((0.286*27.*M_PI*me*e)/
                                                          pow(t (r), 2.)), (p - 2)/2.))/((3 - p)*eps_b),
                                               0.5) - 1)));}
    
    else if (Gamma(r) <= Gamma_DN && nu < nu_c2(r, i, j))
    {
        fnu =
        1/(pow(distance, 2)*(3*p - 1)*(Gamma_DN - 1))*
        pow(doppler(r, i, j), (p + 5)/
             2.)*(4.775*pow(10, 2.)*pow(1.19*pow(10, 6.), (p - 1)/2.))*
        (p - 1)*pow(eps_b, (p + 1)/4.)*pow(n, (p + 5)/4.)*
        pow(r, 3.)*pow(Gamma(r), (p + 1)/4.)*
        pow(Gamma(r) - 1, (p + 5)/4.)*
        sin(theta[i])*pow(nu, (1 - p)/2.);}
    
    else
    {
        fnu = (pow(doppler(r, i, j), (p + 6)/
                    2.)*(1.671*pow(10, 15.)*
                         pow(1.19*pow(10, 6.), (p - 1)/2.))*(p - 1)*
               pow(eps_b, (p - 2)/4.)*
               pow(n, (p + 2)/4.)*pow(r, 3.)*
               pow(Gamma(r), (p - 2)/4.)*
               pow(Gamma(r) - 1, (p + 2)/4.)*
               sin(theta[i]))/
        (pow(nu, p/2.)*(pow(distance, 2)*(3*p - 1)*(Gamma_DN - 1)*t (r)*
                         (1 + 
                          0.5*(pow(1 + (4*eps_e*
                                         pow((4*pow(10, 12.)*pow(sigmat, 2.)*
                                               pow(32*M_PI*mp, 1.5)*pow(clight, 2.)*
                                               pow((p - 2)/(p - 1), 2.)*pow(eps_e, 2.)*
                                               pow(eps_b, 2.)*pow(n, 2.)*
                                               pow(Gamma(r) - 1, 4.)*pow(Gamma(r), 2.))/
                                              ((0.286*27.*M_PI*me*e)/
                                               pow(t (r), 2.)), (p - 2)/2.))/((3 - p)*eps_b),
                                    0.5) - 1))));}
    
    
    return fnu;
}

void initialize();

void Tarray(){
    int i,j,s;
    
    for(s = 0; s <= smax; s++){
        for(i = 0; i <= imax-1; i++){
            for(j = 0; j <= jmax-1; j++){
                T_array[s][i][j]= tobs(r_0[s], i, j);
                
            }
        }
    }
    
}
void Tarray();

int main (){
   

    //double r=6e17;
    //int itest=1;
    //int jtest=1;
    //double thetaobs=0;
    //int stest=smax-1;
    struct timeval start_time, stop_time, elapsed_time;
    gettimeofday(&start_time,NULL);
    
    initialize();
    Tarray();
    int i,j,s,k;
    /*for(s = 0; s <= smax; s++){
     for(i = 0; i <= imax-1; i++){
     for(j = 0; j <= jmax-1; j++){
     T_array[s][i][j]= tobs(r_0[s], i, j);
     
     }
     }
     }*/

    int Nmax=ceil(sqrt(smax));
    int kmax=ceil(smax/Nmax);
    //int l[smax+1][imax][jmax];
    //double Fnu_obs[smax+1][imax][jmax];
    double F_nutotal[smax+1];
    double E_tot[smax+1];
    //double nu_mobs[smax+1];
    //double nu_cobs[smax+1];
    double gamma_obs[smax+1];
    double gamma_lab[smax+1];
    double t_lab[smax+1];
    for(i=0; i<=imax-1; i++){
        for(j=0; j<=jmax-1; j++){
         for (s=0; s<=smax; s++) {
             if (T_array[0][i][j]>Tobs[s] || T_array[smax][i][j]<Tobs[s]) {
                 printf("Tobs[%d,%d,%d] is not within bounds of r_0[0] and r_0[%d], decrease (increase) rmin (rmax)\n",s,i,j,smax);
                 return -1;
             }

             else{
                for (k=0; k<=kmax-1; k++) {
                    if (T_array[k*Nmax][i][j]<=Tobs[s] && Tobs[s]<=T_array[smax ^ ((Nmax*(k+1) ^ smax) & -(Nmax*(n+1) < smax))][i][j]) {
                        int dmax,d;
                        if (k==kmax-1) {
                            dmax=smax-Nmax*(kmax-1);
                        }
                        else{dmax=Nmax;}
                        for (d=0; d<=dmax-1; d++) {
                            if (T_array[k*Nmax+d][i][j]<=Tobs[s] && Tobs[s]<=T_array[Nmax*k+d+1][i][j]) {
                                l[s][i][j]=Nmax*k+d;
                                //printf("l[%d,%d,%d]=%d",s,i,j,l[s][i][j]);
                                goto next_s;
                            }
                            //break;
                        } //break;
                    }
                }
             }
             next_s: continue;}
            
           
                
            
            //else{continue;}
            
        }
    }
    
    for(s = 0; s <= smax; s++){
        for(i = 0; i <= imax-1; i++){
            for(j = 0; j <= jmax-1; j++){
                r_obs[s][i][j]=r_0[l[s][i][j]]+(Tobs[s]-T_array[l[s][i][j]][i][j])*((r_0[l[s][i][j]+1]-r_0[l[s][i][j]])/(T_array[l[s][i][j]+1][i][j]-T_array[l[s][i][j]][i][j]));
                //printf("r_obs[%d,%d,%d]=%.3e\n",s,i,j,r_obs[s][i][j]);
                
            }
        }
    }
  /*  for(s = smax-10; s <= smax; s++){
       for(i = 1; i <= 1; i++){
            for(j = 1; j <= 1; j++){
        printf("nucobs[%d,%d,%d]=%.5e\n",s,i,j,nu_c(r_obs[s][i][j],i,j));
            }
    
        }
   }*/
    
    for(s = 0; s <= smax; s++){
        for(i = 0; i <= imax-1; i++){
            for(j = 0; j <= jmax-1; j++){
                Fnu_obs[s][i][j]=F_nu(r_obs[s][i][j], i, j);
                E_obs[s][i][j]=Einjected(r_obs[s][i][j]);

            
            }
        }
    }



    s=0;
    while (s<=smax){
        double sum=0;
        double sum2=0;
        for (i=0; i<=imax-1; i++) {
            for (j=0; j<=jmax/2-1; j++) {
                sum+=Fnu_obs[s][i][j];
                sum2+=E_obs[s][i][j];
                
            }
        }
        F_nutotal[s]=2*sum*thetadiv*phidiv*1000; //multiply by 2 for one half in phi,multiply by 1000 to convert to micro Jy
        E_tot[s]=2*sum2/(imax*jmax);
        s++;
    }
    //get values of nu_m and nu_c for one portion of blastwave as a representative value
   for(s = 0; s <= smax; s++){
                //nu_mobs[s]=nu_m(r_obs[s][59][0], 59, 0);
                //nu_cobs[s]=nu_c2(r_obs[s][59][0], 59, 0);
                //gamma_obs[s]=Gamma(r_obs[s][30][0]);
                gamma_lab[s]=Gamma(r_0[s]);
                t_lab[s]=t(r_0[s]);
                

            }
    
    char buf[0x100];
    snprintf(buf, sizeof(buf), "/Users/Adi/Documents/intro/lowp/KN_AG_%.0e_%g_%g_%0.2eHz_test.txt",n,alpha,beta0,nu);
    FILE *filePtr= fopen(buf,"w");
    //filePtr = fopen(buf,"w");
    fprintf(filePtr, "alpha=%.3g,n=%.3g,beta0=%.3g, thetaobs=%.3g, p=%.3g, E=%.3g erg,\nepse=%.3g, epsb=%.3g, nuobs=%.3g\n\nTime(yr)\t F_nu(microJy)\n",alpha,n,beta0,thetaobs,p,E,eps_e,eps_b,nu);
    s=0;
    for (s = 0; s <= smax; s++) {
        fprintf(filePtr, "%.5g\t\t%.5g\n", Tobs[s]/(86400*365), F_nutotal[s]);
    }
    fclose(filePtr);
    
    
 /*   char buf2[0x100];
    snprintf(buf2, sizeof(buf2), "/Users/Adi/Documents/intro/KNspec_nu_c_%g.txt",alpha);
    FILE *filePtr2= fopen(buf2,"w");
    //filePtr = fopen(buf,"w");
    fprintf(filePtr2, "alpha=%.3g,n=%.3g,beta0=%.3g, thetaobs=%.3g, p=%.3g, E=%.3g erg,\nepse=%.3g, epsb=%.3g\n\nTime(yr)\t nu_c(Hz)\n",alpha,n,beta0,thetaobs,p,E,eps_e,eps_b);
    s=0;
    for (s = 0; s <= smax; s++) {
        fprintf(filePtr2, "%.5g\t\t%.5g\n", Tobs[s]/(86400*365), nu_cobs[s]);
    }
    fclose(filePtr2);*/
     /* char buf2[0x100];
    snprintf(buf2, sizeof(buf2), "/Users/Adi/Documents/intro/E_tot_%g_%g.txt",alpha,beta0);
    FILE *filePtr2= fopen(buf2,"w");
    //filePtr = fopen(buf,"w");
    fprintf(filePtr2, "alpha=%.3g,n=%.3g,beta0=%.3g, thetaobs=%.3g, p=%.3g, E=%.3g erg,\nepse=%.3g, epsb=%.3g\n\nTime(yr)\t nu_c(Hz)\n",alpha,n,beta0,thetaobs,p,E,eps_e,eps_b);
    s=0;
    for (s = 0; s <= smax; s++) {
        fprintf(filePtr2, "%.5g\t\t%.5g\n", Tobs[s]/(86400*365), E_tot[s]);
    }
    fclose(filePtr2);*/
    //printf("(Fnutot[%d]=%.5e\n)",25,F_nutotal[25]);
    
    
    //printf("(l[%d,%d,%d]=%.5e\n)",1,1,1,T_array[6][1][1]);

    //printf("theta[%d]=%.5e\n,phi[%d]=%.5e\n,Tobs[%d]=%.5e\n,r_0[%d]=%.5e\n,T_array[%d,%d,%d]=%.5e\n,Nmax=%d,nmax=%d",itest,theta[itest],jtest,phi[jtest],stest,Tobs[stest],stest,r_0[stest],stest,itest,jtest,T_array[stest][itest][jtest],Nmax,nmax);
    
    //printf("gamma_0=%.5e\n, E_iso=%.5e\n, rdec[%g]=%.5e\n, betagamma[%g,%g]=%.5e\n,beta[%g]=%.5e\n,Gamma[%g]=%.5e\n,tlab[%g]=%.5e\n,tobs[%g,%g,%g,%g]=%.5e\n, theta[%d]=%.5e\n",g0(),E_iso(),n,R_dec(),r,alpha,betagamma(r),r,beta(r),r,Gamma(r),r,tlab(r),r,theta[itest],phi[jtest],thetaobs,tobs(r,itest,jtest,thetaobs),itest,theta[itest]);
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time);
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
    return 0;
 
              }




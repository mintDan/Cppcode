//PolarCoordsRK4
//Dan Krog
#include <iostream>
#include <array>
#include <fstream>
#include <math.h>
using namespace std;


double fx(double x0, double y0, double a)
{
    return -y0+a*x0*(x0*x0+y0*y0);
}

double fy(double x0,double y0,double a)
{
    return x0+a*y0*(x0*x0+y0*y0);
}

//Polar coordinates
double fr(double r0,double theta0, double a)
{
    return a*pow(r0,3);
}
double ftheta(double r0, double theta0, double a)
{
    return 1;
}

const double a = -1;
const double dt = 0.1;
const int nmax = 200;

int main()
{
    cout << "START PP HERE \n";

    //CARTESIAN COORDS
    //parameters

    array<double,nmax> t {0};
    array<double,nmax> x {0};
    array<double,nmax> y {2};

    double k1x,k1y, k2x, k2y, k3x, k3y, k4x, k4y;
    for (int n=0;n<y.size();++n)
    {
        t[n+1] = t[n] + dt;

        k1x = fx(x[n],y[n],a)*dt;
        k1y = fy(x[n],y[n],a)*dt;

        k2x = fx(x[n]+0.5*k1x,y[n]+0.5*k1y,a)*dt;
        k2y = fy(x[n]+0.5*k1x,y[n]+0.5*k1y,a)*dt;

        k3x = fx(x[n]+0.5*k2x,y[n]+0.5*k2y,a)*dt;
        k3y = fy(x[n]+0.5*k2x,y[n]+0.5*k2y,a)*dt;

        k4x = fx(x[n]+k3x,y[n]+k3y,a)*dt;
        k4y = fy(x[n]+k3x,y[n]+k3y,a)*dt;

        x[n+1] = x[n]+ (1.0/6.0)*(k1x+2*k2x+2*k3x+k4x);
        y[n+1] = y[n]+ (1.0/6.0)*(k1y+2*k2y+2*k3y+k4y);

    }
    // POLAR COORDS
    array<double,nmax> r {2};
    array<double,nmax> theta {3.141592/2.0}; //we start at approx x=0,y=3,

    double k1r,k1theta, k2r, k2theta, k3r, k3theta, k4r, k4theta;
    for (int n=0;n<y.size();++n)
    {
        k1r = fr(r[n],theta[n],a)*dt;
        k1theta = ftheta(r[n],theta[n],a)*dt;

        k2r = fr(r[n]+0.5*k1r,theta[n]+0.5*k1theta,a)*dt;
        k2theta = ftheta(r[n]+0.5*k1r,theta[n]+0.5*k1theta,a)*dt;

        k3r = fr(r[n]+0.5*k2r,theta[n]+0.5*k2theta,a)*dt;
        k3theta = ftheta(r[n]+0.5*k2r,theta[n]+0.5*k2theta,a)*dt;

        k4r = fr(r[n]+k3r,theta[n]+k3theta,a)*dt;
        k4theta = ftheta(r[n]+k3r,theta[n]+k3theta,a)*dt;

        r[n+1] = r[n]+ (1.0/6.0)*(k1r+2*k2r+2*k3r+k4r);
        theta[n+1] = theta[n]+ (1.0/6.0)*(k1theta+2*k2theta+2*k3theta+k4theta);

    }

    //Let's convert the polar coordinates to cartesian
    array<double,nmax> xconv {};
    array<double,nmax> yconv {}; //we start at approx x=0,y=3,
    for (int n=0;n<y.size();++n)
    {
        xconv[n] = r[n]*cos(theta[n]);
        yconv[n] = r[n]*sin(theta[n]);
    }


    ofstream myfile;
    myfile.open("PP.txt");
    for (int n=1;n<y.size();++n)
    {
        myfile << t[n] << " ";
        myfile << x[n] << " ";
        myfile << y[n] << " ";
        myfile << r[n] << " ";
        myfile << theta[n] << " ";
        myfile << xconv[n] << " ";
        myfile << yconv[n] << "\n";

    }
    myfile.close();

    //end program
    return 0;
}


//Malthusian model
//Dan Krog
#include <iostream>
#include <array>
#include <fstream>
#include <math.h>
using namespace std;


double UpdateL(double L0, double t0)
{
    double L1,s,A,b,dt,d,k;
    L1 = L0 + dt*(s*A*exp(k*t0)*pow(L0,b)-d*L0);
    return L1;
}

int main()
{
    //parameters
    double b = 0.5;
    double dt = 0.5;
    double A = 100;
    double s = 0.1;
    double d = 0.025;
    double k = 0.005;


    cout << "START MALTHUS HERE \n";
    const int nmax = 800;

    array<double,nmax> t {0}; //time, some unit
    array<double,nmax> L {1000000}; //Population, number of people
    array<double,nmax> Y {A*pow(L[0],b)}; //GDP
    array<double,nmax> y {Y[0]/L[0]}; //GDP/capita

    //Index numbers for easier plotting
    array<double,nmax> Lind {100};
    array<double,nmax> Yind {100};
    array<double,nmax> yind {100};


    //looping
    int n;
    int cutoff = 50;
    double Ymultiple = 0;

    for (n=0;n<L.size();++n)
    {
        t[n+1] = t[n]+ dt;
        L[n+1] = L[n] + dt*(s*A*exp(k*t[n])*pow(L[n],b)-d*L[n]);
        //L[n+1] = UpdateL(L[n],t[n]);
        Y[n+1] = A*exp(k*t[n])*pow(L[n],b);
        y[n+1] = Y[n+1]/L[n+1];

        //Index numbers

        Lind[n+1] = 100*L[n+1]/L[0];
        Yind[n+1] = 100*Y[n+1]/Y[0];
        yind[n+1] = 100*y[n+1]/y[0];
        //cout << t[n] << " " << L[n] << " " << Y[n] <<"\n";
    }


    ofstream myfile;
    myfile.open("Malth.txt");
    for (n=1;n<L.size();++n)
    {
        myfile << t[n] << " ";
        myfile << L[n] << " ";
        myfile << Y[n] << " ";
        myfile << y[n] << " ";
        myfile << Lind[n] << " ";
        myfile << Yind[n] << " ";
        myfile << yind[n] << "\n";
    }
    myfile.close();

    //end program
    return 0;
}


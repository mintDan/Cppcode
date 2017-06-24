/* iteration1.c */

//#include <stdio.h>
#include <math.h>
//#include <stdlib.h>
#include <array>
#include <fstream>
#include <iostream>
#include <tuple>
#include <cmath>
using namespace std;


double stuff(double);
double stuff(double a)
{
    printf("Hi from stuff");
    printf("a = %3d",a);
}

void ptest(int, double *);
void ptest(int n, double *fact)
{
    *fact = 1;
    printf("n = %3d",n);

}





double f1(double, double,double);
double f1(double t, double x, double y)
{
    double fx;
    fx = y;
    //cout << "f here" << fxy;
    //cout << fxy;
    return fx;

}


double f2(double,double,double);
double f2(double t, double x, double y)
{
    double fv,k;
    k = 1;
    fv = -k*x;
    //cout << "f here" << fxy;
    //cout << fxy;
    return fv;

}

//lidt usikker på arguments tilo runge kutta, men will find out
double RungeKutta(double,double,double,double);
double RungeKutta(double t, double x, double y, double h)
{
    double k1,k2,k3,k4,y1;
    double fakt = 0.1666666;
    k1 = h*f1(t,x,y);
    k2 = h*f1(t,x+h/2,y+k1/2);
    k3 = h*f1(t,x+h/2,y+k2/2);
    k4 = h*f1(t,x+h,y+k3);
    //cout << f(x,y) << " ";
    y1 = y + fakt*(k1+2*k2+2*k3+k4);
    //cout << k1 << " " << k2 << "\n" ;
    //cout << y1 << "\n";
    return y1;
}


std::tuple<double, double> RungeKutta1(double,double,double,double);
std::tuple<double, double> RungeKutta1(double t, double x, double y, double h)
{
    //jeg har her x+h/2....
    //MEN... h er jo dt... og x er position, så units passer ikke!
    //Det er noget med der er en k1x,k1y, right?
    //Hvis jeg husker rigtigt...
    //Så hmmm... måske lidt som partial derivatives
    //først step man x+k1x, og holder y constant?
    //så step man y+k1y og holder x constant? etc maybe?
    double k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,x1,y1;
    double fakt = 0.1666666;

    k1x = h*f1(t,x,y);
    k1y = h*f2(t,x,y);

    k2x = h*f1(t+h/2,x+0.5*k1x,y+0.5*k1y);
    k2y = h*f2(t+h/2,x+0.5*k1x,y+0.5*k1y);

    k3x = h*f1(t+h/2,x+0.5*k2x,y+0.5*k2y);
    k3y = h*f2(t+h/2,x+0.5*k2x,y+0.5*k2y);

    k4x = h*f1(t+h,x+k3x,y+k3y);
    k4y = h*f2(t+h,x+k3x,y+k3y);
    //cout << f(x,y) << " ";
    x1 = x + fakt*(k1x+2*k2x+2*k3x+k4x);
    y1 = y + fakt*(k1y+2*k2y+2*k3y+k4y);
    //cout << k1 << " " << k2 << "\n" ;
    //x1 = 4.2;
    //y1 = 5.4;
    //cout << x1 << " " << y1 << "\n";

    return std::make_tuple(x1, y1);
}


//Baseball functions
//Lige nu er BÅDE variables OG functions smame name...
//måske ikke så godt
double fx(double);
double fx(double vx){
    double fx;
    fx = vx;
    return fx;
}
double fy(double);
double fy(double vy){
    double fy;
    fy = vy;
    return fy;
}
double fz(double);
double fz(double vz){
    double fz;
    fz = vz;
    return fz;
}

double Fv(double, double, double);
double Fv(double v,double vd, double inc){
    double Fv;
    Fv = 0.0039 + 0.0058/(1+exp((v-vd)/inc));
    return Fv;
}

double fvx(double,double,double,double,double,double,double,double);
double fvx(double Fv,double v,double vx,double B, double w, double vz,double phi,double vy){
    double afvx;
    afvx = -Fv*v*vx + B*w*(vz*sin(phi)-vy*cos(phi));
    return afvx;
}
double fvy(double, double, double, double, double,double,double,double);
double fvy(double Fv,double v, double vx,double B, double w, double vz,double phi,double vy){
    double afvy;
    afvy = -Fv*v*vy + B*w*vx*cos(phi);
    return afvy;
}
double fvz(double,double,double,double,double,double,double,double,double);
double fvz(double Fv,double v, double vx,double B, double w, double vz,double phi,double vy,double g){
    double afvz;
    afvz = -g-Fv*v*vz - B*w*vx*sin(phi);
    return afvz;
}

//Baseball RungeKutta
std::tuple<double, double, double, double, double, double> RungeKutta2(double,double,double,double,double,
                                                                       double,double,double,double,double,
                                                                       double,double,double,double);
std::tuple<double, double, double, double, double, double> RungeKutta2(double t, double x, double y,
                                                                        double h,double z, double vx,
                                                                         double vy, double vz,double vd,double inc,
                                                                          double B,double phi,double g,
                                                                          double w){
    //jeg har her x+h/2....
    //MEN... h er jo dt... og x er position, så units passer ikke!
    //Det er noget med der er en k1x,k1y, right?
    //Hvis jeg husker rigtigt...
    //Så hmmm... måske lidt som partial derivatives
    //først step man x+k1x, og holder y constant?
    //så step man y+k1y og holder x constant? etc maybe?
    double k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y,x1,y1,z1;
    double vx1,vy1,vz1;
    double k1z,k2z,k3z,k4z;
    double k1vx,k2vx,k3vx,k4vx;
    double k1vy,k2vy,k3vy,k4vy;
    double k1vz,k2vz,k3vz,k4vz;
    double Fvcalc;
    double v;
    double fakt = 0.1666666;
    v = sqrt(vx*vx + vy*vy + vz*vz);
    Fvcalc = Fv(v,vd,inc);

    k1x = h*fx(vx);
    k1y = h*fy(vy);
    k1z = h*fz(vz);
    k1vx = h*fvx(Fvcalc,v,vx,B,w,vz,phi,vy);
    k1vy = h*fvy(Fvcalc,v,vx,B,w,vz,phi,vy);
    k1vz = h*fvz(Fvcalc,v,vx,B,w,vz,phi,vy,g);

    v = sqrt((vx+0.5*k1vx)*(vx+0.5*k1vx) + (vy+0.5*k1vy)*(vy+0.5*k1vy) + (vz+0.5*k1vz)*(vz+0.5*k1vz));
    Fvcalc = Fv(v,vd,inc);

    k2x = h*fx(vx+0.5*k1vx);
    k2y = h*fy(vy+0.5*k1vy);
    k2z = h*fz(vz+0.5*k1vz);
    k2vx = h*fvx(Fvcalc,v,vx+0.5*k1vx,B,w,vz+0.5*k1vz,phi,vy+0.5*k1vy);
    k2vy = h*fvy(Fvcalc,v,vx+0.5*k1vx,B,w,vz+0.5*k1vz,phi,vy+0.5*k1vy);
    k2vz = h*fvz(Fvcalc,v,vx+0.5*k1vx,B,w,vz+0.5*k1vz,phi,vy+0.5*k1vy,g);

    v = sqrt((vx+0.5*k2vx)*(vx+0.5*k2vx) + (vy+0.5*k2vy)*(vy+0.5*k2vy) + (vz+0.5*k2vz)*(vz+0.5*k2vz));
    Fvcalc = Fv(v,vd,inc);

    k3x = h*fx(vx+0.5*k2vx);
    k3y = h*fy(vy+0.5*k2vy);
    k3z = h*fz(vz+0.5*k2vz);
    k3vx = h*fvx(Fvcalc,v,vx+0.5*k2vx,B,w,vz+0.5*k2vz,phi,vy+0.5*k2vy);
    k3vy = h*fvy(Fvcalc,v,vx+0.5*k2vx,B,w,vz+0.5*k2vz,phi,vy+0.5*k2vy);
    k3vz = h*fvz(Fvcalc,v,vx+0.5*k2vx,B,w,vz+0.5*k2vz,phi,vy+0.5*k2vy,g);

    v = sqrt((vx+k3vx)*(vx+k3vx) + (vy+k3vy)*(vy+k3vy) + (vz+k3vz)*(vz+k3vz));
    Fvcalc = Fv(v,vd,inc);

    k4x = h*fx(vx+k3vx);
    k4y = h*fy(vy+k3vy);
    k4z = h*fz(vz+k3vz);
    k4vx = h*fvx(Fvcalc,v,vx+k1vx,B,w,vz+k1vz,phi,vy+k1vy);
    k4vy = h*fvy(Fvcalc,v,vx+k1vx,B,w,vz+k1vz,phi,vy+k1vy);
    k4vz = h*fvz(Fvcalc,v,vx+k1vx,B,w,vz+k1vz,phi,vy+k1vy,g);

    //cout << f(x,y) << " ";
    x1 = x + fakt*(k1x+2*k2x+2*k3x+k4x);
    y1 = y + fakt*(k1y+2*k2y+2*k3y+k4y);
    z1 = z + fakt*(k1z+2*k2z+2*k3z+k4z);
    vx1 = vx + fakt*(k1vx+2*k2vx+2*k3vx+k4vx);
    vy1 = vy + fakt*(k1vy+2*k2vy+2*k3vy+k4vy);
    vz1 = vz + fakt*(k1vz+2*k2vz+2*k3vz+k4vz);
    //cout << k1 << " " << k2 << "\n" ;
    //x1 = 4.2;
    //y1 = 5.4;
    //cout << x1 << " " << y1 << "\n";

    return std::make_tuple(x1, y1, z1, vx1, vy1, vz1);
}
//const int nmax = 200;


/* Set max. allowable no. of iterations */
//#define NITER 30
int main()
{
    //cout << "START RK LET'S GO HERE \n";
    //int count = 0;
    //do
    //    {++count;
    //    printf("lol");
    //    }
    //while (count < 10); // Test for convergence
    //stuff(3);

    //double fact;
    //ptest(4,&fact);
    cout << "\n";
    double k = 1,dt = 0.00001, dt2,dytest;
    double eps = 0.001;
    double dy = 0.01;
    double v0 = 37.9984;
    double theta = 0.0174532925;
    double phi = 0;
    double vd = 35;
    double inc = 5;
    double w = 188.495559;
    double g = 9.81;
    double B = 4.1/(10*10*10*10);
    double L = 18.44;
    double T = L/v0;
    const int nmax = 250;//T/dt;
    array<double,nmax> x {0};
    array<double,nmax> v {sqrt(k)};
    array<double,nmax> t {0};
    array<double,nmax> y {0};
    array<double,nmax> z {0};
    array<double,nmax> vx {v0*cos(theta)};
    array<double,nmax> vy {0};
    array<double,nmax> vz {v0*sin(theta)};
    //array<double,nmax> exps {1};

    double newvalx1,newvaly1,newvalz1,newvalx2,newvaly2,newvalz2;
    double newvx1,newvy1,newvz1,newvx2,newvy2,newvz2;


    for (int n=0;n<v.size();++n)
    {
        tie(newvalx2,newvaly2,newvalz2,newvx2,newvy2,newvz2) = RungeKutta2(t[n], x[n], y[n],dt,z[n],vx[n], vy[n], vz[n],vd,inc,B,phi,g,w);
        tie(newvalx1,newvaly1,newvalz1,newvx1,newvy1,newvz1) = RungeKutta2(t[n], x[n], y[n],0.5*dt,z[n],vx[n], vy[n], vz[n],vd,inc,B,phi,g,w);

        dy = std::abs(newvaly2-newvaly1);

        if (dy < eps)
        {


            t[n+1] = t[n]+dt;
            x[n+1] = newvalx2;
            y[n+1] = newvaly2;
            z[n+1] = newvalz2;
            vx[n+1] = newvx2;
            vy[n+1] = newvy2;
            vz[n+1] = newvz2;
            //hvis dy < eps, skal vi da også lige prøve
            //at se om vi kan gøre dt STØRRE.
            //Men, måske først EFTER dette step er lavet, så?
            while (dy < eps) {
                //skal måske lige lave en ny dt2 til testing
                //så jeg ikke ALTID overwrite den gamle dt...
                //ellers så overwrite vi den måske eng ang for meget
                dt2 = 2.0*dt;

                tie(newvalx2,newvaly2,newvalz2,newvx2,newvy2,newvz2) = RungeKutta2(t[n], x[n], y[n],dt2,z[n],vx[n], vy[n], vz[n],vd,inc,B,phi,g,w);
                tie(newvalx1,newvaly1,newvalz1,newvx1,newvy1,newvz1) = RungeKutta2(t[n], x[n], y[n],0.5*dt2,z[n],vx[n], vy[n], vz[n],vd,inc,B,phi,g,w);

                dytest = std::abs(newvaly2-newvaly1);
                if (dytest < eps){
                    dy = dytest;
                    dt = dt2;
                    cout << "Made dt bigger";
                    cout << "\n";
                }
                else{
                    break;
                }
            }
        }


        else if (dy > eps)
        {
            dt = 0.5*dt;
            cout << "Made dt smaller";
            cout << "\n";
            while (dy>eps)
                {
                    tie(newvalx2,newvaly2,newvalz2,newvx2,newvy2,newvz2) = RungeKutta2(t[n], x[n], y[n],dt,z[n],vx[n], vy[n], vz[n],vd,inc,B,phi,g,w);
                    tie(newvalx1,newvaly1,newvalz1,newvx1,newvy1,newvz1) = RungeKutta2(t[n], x[n], y[n],0.5*dt,z[n],vx[n], vy[n], vz[n],vd,inc,B,phi,g,w);

                    dy = std::abs(newvaly2-newvaly1);
                }
            if (dy < eps)
            {
                t[n+1] = t[n]+dt;
                x[n+1] = newvalx2;
                y[n+1] = newvaly2;
                z[n+1] = newvalz2;
                vx[n+1] = newvx2;
                vy[n+1] = newvy2;
                vz[n+1] = newvz2;
            }
        }

        //exps[n+1] = exp(x[n]);
        cout << newvalx2 << " " << newvaly2 << "\n";
        //cout << x[n+1] << " " << exps[n+1] << " " << y[n+1] << "\n";
        //cout << x[n+1] << " " << v[n+1] << "\n";
    }

    ofstream myfile;
    myfile.open("Baseball.txt");
    for (int n=1;n<v.size();++n)
    {
        myfile << t[n] << " ";
        myfile << x[n] << " ";
        myfile << y[n] << " ";
        myfile << z[n] << " ";
        myfile << vx[n] << " ";
        myfile << vy[n] << " ";
        myfile << vz[n] << " " << "\n";

    }
    myfile.close();

    return 0;
}

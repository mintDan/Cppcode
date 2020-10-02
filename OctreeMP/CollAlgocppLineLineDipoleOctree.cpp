#include<iostream>
#include <omp.h>
#include <math.h>
#include <cmath>        //Christoph: needed for std::abs
#include <random>

//#include <tuple>
 
//using namespace std;
/*int print_array(double * in_array, double * out_array, int size);*/
//double ** particle_list [200][15]


//g++ -shared -o CollAlgocppLineLineDipoleOctree.so -fopenmp -fPIC CollAlgocppLineLineDipoleOctree.cpp


//#include <iostream>
using namespace std;
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
//#include <math.h>
//#include <cmath>

//does both Node and Octree contain children?



//For CPU, memory usage, process id
#include "sys/types.h"
#include "sys/sysinfo.h"
#include <unistd.h> 
#include <fstream>



//For variable length array?
#include <memory>



struct sysinfo memInfo;





class Octree;

class Octree {

    

    //const int npmax = 300;


    
   // static Octree Qtsub1;
    
    
    double cx;
    double cy;
    double cz;
    double x1_boundary;
    double x2_boundary;
    double y1_boundary;
    double y2_boundary;
    double z1_boundary;
    double z2_boundary;
    double xhalflength;
    double yhalflength;
    double zhalflength;
    
      // Node *no1;

    Octree *ne1; 
    Octree *nw1; 
    Octree *sw1; 
    Octree *se1;
    Octree *ne2; 
    Octree *nw2; 
    Octree *sw2; 
    Octree *se2; 
    
    
    public:





        int npmax;// = 300;

        int max_level = 2;

        bool Leaf = 1;

        int level;

        int NodeNumber;
    
    
    //bool Leaf;
   // Node n;
    




    int pcounter = 0;





    //int points[300];
    //int n = 10;
     //double* points = new double[npmax];



     //https://stackoverflow.com/questions/25153153/initializing-array-of-variable-size-inside-a-class?noredirect=1&lq=1
     //One way here
     //int *pointsarr;


     //https://stackoverflow.com/questions/26198052/how-do-i-initialize-a-variable-size-array-in-a-c-class
     //one way here
     std::unique_ptr<int[]> points;



     
     Octree(int,double,double,double,double,double,double,double,double,double,int,int);



     //https://www.tutorialspoint.com/cplusplus/cpp_constructor_destructor.htm
     //her, destructor er sat under public:

     //https://www.geeksforgeeks.org/destructors-c/
     // ogsaa her, destructor er public


     /*~Octree() //destructor defined
       {
        //Destructor
       }*/

     //Destructor, to clear the tree recursiely, to remove it from memory after it has been used.
     ~Octree();





     //Octree methods
     //Insert points, split cells, collide particles
     
     //insert a point into the cell that it fits into
     void Insert (int,double particle_list[][15],int points_array[][64],int NN);
     
     // Helper function for the Insert function, due to duplicate coding sections
     void Insert_Helper (int, double particle_list[][15],int points_array[][64],int NN);


     //Split a node into eight new nodes, once the point capacity has been reached
     void split ();


     //collide particles in cells with the line-line algorithm
     int CollideParticles (double particle_list[][15],int n_S, int n_M,  double * Mass_adj_Precalc, double * R_adj_Precalc, int DimRange, double dr[][4],int ncollision, int points_array[][64],int NN);



};




Octree::Octree(int levelin, double cxin, double cyin, double czin, double x1_boundaryin , double x2_boundarin,
    double y1_boundaryin ,double y2_boundaryin, double z1_boundaryin, double z2_boundaryin, int npmaxin,int NN)
    : npmax{npmaxin}
    , points{new int[npmax]} //https://stackoverflow.com/questions/26198052/how-do-i-initialize-a-variable-size-array-in-a-c-class
{

    cx = cxin;
    cy = cyin;
    cz = czin;
    x1_boundary=x1_boundaryin;
    x2_boundary=x2_boundarin;
    y1_boundary=y1_boundaryin;
    y2_boundary=y2_boundaryin;
    z1_boundary=z1_boundaryin;
    z2_boundary=z2_boundaryin;
    xhalflength = (x2_boundary-x1_boundary)/2;
    yhalflength = (y2_boundary-y1_boundary)/2;
    zhalflength = (z2_boundary-z1_boundary)/2;
    Leaf = 1;



    ne1 = NULL; 
    nw1 = NULL; 
    sw1 = NULL; 
    se1 = NULL;
    ne2 = NULL; 
    nw2 = NULL; 
    sw2 = NULL; 
    se2 = NULL; 

/*


    ne1 = nullptr;
    nw1 = nullptr;
    sw1 = nullptr;
    se1 = nullptr;
    ne2 = nullptr;
    nw2 = nullptr;
    sw2 = nullptr;
    se2 = nullptr; */





    
    level =levelin;
    NodeNumber = NN;
    //npmax = npmaxin;

    //https://stackoverflow.com/questions/25153153/initializing-array-of-variable-size-inside-a-class?noredirect=1&lq=1
    //points = new pointsarr[npmax]; // <= allocate array
    //n = a;


    
    
}


//Desctructor
Octree::~Octree() {
    //Clear the octree recursively to free up memory used



    //if (ne1 != NULL)

    //ne1->~Octree();
    //nw1->~Octree(); 
    //sw1->~Octree(); 
    //se1->~Octree();
    //ne2->~Octree(); 
    //nw2->~Octree(); 
    //sw2->~Octree(); 
    //se2->~Octree(); 

/*
    if (ne1 != NULL) ne1->~Octree(); 
    if (nw1 != NULL) nw1->~Octree(); 
    if (sw1 != NULL) sw1->~Octree(); 
    if (se1 != NULL) se1->~Octree();
    if (ne2 != NULL) ne2->~Octree();
    if (nw2 != NULL) nw2->~Octree();
    if (sw2 != NULL) sw2->~Octree();  
    if (se2 != NULL) se2->~Octree(); */

    //for ( int i = 0; i < this->count; i++ ) {
    //    free ( this->points[i].data );
    //}

    //delete points;
    //delete &points;
    //delete []points;

/*    if (ne1 != nullptr) delete ne1; 
    if (nw1 != nullptr) delete nw1; 
    if (sw1 != nullptr) delete sw1; 
    if (se1 != nullptr) delete se1;
    if (ne2 != nullptr) delete ne2; 
    if (nw2 != nullptr) delete nw2; 
    if (sw2 != nullptr) delete sw2; 
    if (se2 != nullptr) delete se2; */




    if (ne1 != NULL) delete ne1; 
    if (nw1 != NULL) delete nw1; 
    if (sw1 != NULL) delete sw1; 
    if (se1 != NULL) delete se1;
    if (ne2 != NULL) delete ne2; 
    if (nw2 != NULL) delete nw2; 
    if (sw2 != NULL) delete sw2; 
    if (se2 != NULL) delete se2; 
    
   //cout << "Object is being deleted" << endl;
}




void Octree::Insert (int point_index, double particle_list[][15],int points_array[][64],int NN) { 
//Insert a point into a cell

//for max level, how do we want it?
//if node is leaf, then insert point to points list
//if not, call insert again...
//no splitting from this function...



//ahh jeg tror vi har et problem med points array... den er nok ikke shared mellem processes...
//det er ikke som particle_list, der er shared, maaske... hmmm...


//Maybe I would have to...maybe a multidimensional array??
//PointIndex[500][8] for example... that is passed to the insert function??
//and also to the collide function??
//So PointIndex[500][0] contains the list of indices for the first cell, etc...
//Maybe THEN it is shared??



//First we check the current cell is a leaf. If not leaf, we must check for one of its 8 children nodes that the point fits into.
if (Leaf == 1){
    
    //We are in leaf, and counter is less than npmax
    //Also, the point fits into this cell
    //points[pcounter] = point_index;



    //Insert point into outside defined multidimensional array
    points_array[pcounter][NN] = point_index;


    pcounter +=1;

// Node is not leaf, so we don't insert point_index to this node particle list, but to some subnode
} else {

    // find which existing child node the point fits into and insert there

    //cout << "Not leaf, insert into existing node. pcounter: " << pcounter << endl;

    Insert_Helper(point_index, particle_list,points_array,NN);


    return;

    
} 





return;
} //end Insert function




void Octree::Insert_Helper(int point_index, double particle_list[][15],int points_array[][64],int NN)
{

    if (particle_list[point_index][0]>= nw1->x1_boundary && particle_list[point_index][0]< nw1->x2_boundary && 
        particle_list[point_index][1]>= nw1->y1_boundary && particle_list[point_index][1]< nw1->y2_boundary &&
        particle_list[point_index][2]>= nw1->z1_boundary && particle_list[point_index][2]< nw1->z2_boundary)
    {
       //cout << "inserted " << point_index << endl;
       nw1->Insert(point_index,particle_list,points_array,0);
    }else if (particle_list[point_index][0]>= ne1->x1_boundary && particle_list[point_index][0]< ne1->x2_boundary && 
        particle_list[point_index][1]>= ne1->y1_boundary && particle_list[point_index][1]< ne1->y2_boundary &&
        particle_list[point_index][2]>= ne1->z1_boundary && particle_list[point_index][2]< ne1->z2_boundary)
    {
        //cout << "inserted " << point_index << endl;
        ne1->Insert(point_index,particle_list,points_array,1);
    }else if (particle_list[point_index][0]>= sw1->x1_boundary && particle_list[point_index][0]< sw1->x2_boundary && 
        particle_list[point_index][1]>= sw1->y1_boundary && particle_list[point_index][1]< sw1->y2_boundary &&
        particle_list[point_index][2]>= sw1->z1_boundary && particle_list[point_index][2]< sw1->z2_boundary)
    {
        //cout << "inserted " << point_index << endl;
        sw1->Insert(point_index,particle_list,points_array,2);
    }
    else if (particle_list[point_index][0]>= se1->x1_boundary && particle_list[point_index][0]< se1->x2_boundary && 
        particle_list[point_index][1]>= se1->y1_boundary && particle_list[point_index][1]< se1->y2_boundary &&
        particle_list[point_index][2]>= se1->z1_boundary && particle_list[point_index][2]< se1->z2_boundary)
    {
        //cout << "inserted " << point_index << endl;
        se1->Insert(point_index,particle_list,points_array,3);

    } else if(particle_list[point_index][0]>= nw2->x1_boundary && particle_list[point_index][0]< nw2->x2_boundary && 
        particle_list[point_index][1]>= nw2->y1_boundary && particle_list[point_index][1]< nw2->y2_boundary &&
        particle_list[point_index][2]>= nw2->z1_boundary && particle_list[point_index][2]< nw2->z2_boundary)
    {
       //cout << "inserted " << point_index << endl;
       nw2->Insert(point_index,particle_list,points_array,4);

    }else if (particle_list[point_index][0]>= ne2->x1_boundary && particle_list[point_index][0]< ne2->x2_boundary && 
        particle_list[point_index][1]>= ne2->y1_boundary && particle_list[point_index][1]< ne2->y2_boundary &&
        particle_list[point_index][2]>= ne2->z1_boundary && particle_list[point_index][2]< ne2->z2_boundary)
    {
        //cout << "inserted " << point_index << endl;
        ne2->Insert(point_index,particle_list,points_array,5);

    }else if (particle_list[point_index][0]>= sw2->x1_boundary && particle_list[point_index][0]< sw2->x2_boundary && 
        particle_list[point_index][1]>= sw2->y1_boundary && particle_list[point_index][1]< sw2->y2_boundary &&
        particle_list[point_index][2]>= sw2->z1_boundary && particle_list[point_index][2]< sw2->z2_boundary)
    {
        //cout << "inserted " << point_index << endl;
        sw2->Insert(point_index,particle_list,points_array,6);
    }
    else if (particle_list[point_index][0]>= se2->x1_boundary && particle_list[point_index][0]< se2->x2_boundary && 
        particle_list[point_index][1]>= se2->y1_boundary && particle_list[point_index][1]< se2->y2_boundary &&
        particle_list[point_index][2]>= se2->z1_boundary && particle_list[point_index][2]< se2->z2_boundary)
    {
        //cout << "inserted " << point_index << endl;
        se2->Insert(point_index,particle_list,points_array,7);
    }

    return;
}


void Octree::split ()
{   

    //Thew new Octree actually has to split on to max level on construction...
    //without any particles...
    //so like, pretty stand-alone stuff... it should be a max-level tree even before particles are inserted
    //and then the 100.000 particles can be inserted with OpenMP...
    //Should be faster tbh...
    //No reinsertion, and with multiprocessing... should be something tbh...!!
    //Ok, so MAYBE, the Construction itself doesn't have to call split?
    //I can just call the split() method right after construction... so
    //construct
    //split
    //insert
    //deconstruct
    //4 separate calls like that? why not??


    //hmmm... not so good with inside c++ array...
    //maybe with outside python array???




    
    
    int newlevel = level+1;

    double Newxhalflength = xhalflength/2;
    double Newyhalflength = yhalflength/2;
    double Newzhalflength = zhalflength/2;


    //w1 = new Octree(newlevel, 0.25,0.25,0.25,x1_boundary,x2_boundary,y1_boundary,y2_boundary,z1_boundary,z2_boundary);
    
    

    if (level < max_level){
    //Since level for this node is less than max_level, we split the node into 8 new nodes

    //No longer a leaf node
    //std::cout << level << " " << max_level << std::endl;

    Leaf = 0;

    //#Subquad1
    //#NE
    ne1 = new Octree(newlevel,
    cx+Newxhalflength,
    cy-Newyhalflength,
    cz+Newzhalflength,
    x1_boundary+xhalflength,
    x2_boundary,
    y1_boundary,
    y2_boundary-yhalflength,
    z1_boundary+zhalflength,
    z2_boundary,
    npmax,0);

   // #NW
  // # Subquad2
    nw1 = new Octree(newlevel,
    cx-Newxhalflength,
    cy-Newyhalflength,
    cz+Newzhalflength,
    x1_boundary,
    x2_boundary-xhalflength,
    y1_boundary,
    y2_boundary-yhalflength,
    z1_boundary+zhalflength,
    z2_boundary,
    npmax,1);


    //#subquad3 SW 
    sw1 = new Octree(newlevel,
    cx-Newxhalflength,
    cy-Newyhalflength,
    cz-Newzhalflength,
    x1_boundary,
    x2_boundary-xhalflength,
    y1_boundary,
    y2_boundary-yhalflength,
    z1_boundary,
    z2_boundary-zhalflength,
    npmax,2);


   // #subquad4 SE
    se1 = new Octree(newlevel,
    cx+Newxhalflength,
    cy-Newyhalflength,
    cz-Newzhalflength,
    x1_boundary+xhalflength,
    x2_boundary,
    y1_boundary,
    y2_boundary-yhalflength,
    z1_boundary,
    z2_boundary-zhalflength,
    npmax,3);

   //#subquad5 so this is like 1 but pushed along y axis
    ne2 = new Octree(newlevel,
    cx+Newxhalflength,
    cy+Newyhalflength,
    cz+Newzhalflength,
    x1_boundary+xhalflength,
    x2_boundary,
    y1_boundary+yhalflength,
    y2_boundary,
    z1_boundary+zhalflength,
    z2_boundary,
    npmax,4);


    //#subquad 6 so like 2 but pushed along y
    nw2 = new Octree(newlevel,
    cx-Newxhalflength,
    cy+Newyhalflength,
    cz+Newzhalflength,
    x1_boundary,
    x2_boundary-xhalflength,
    y1_boundary+yhalflength,
    y2_boundary,
    z1_boundary+zhalflength,
    z2_boundary,
    npmax,5);


    //#subquad 7 so like 3 but pushed along y
    sw2 = new Octree(newlevel,
    cx-Newxhalflength,
    cy+Newyhalflength,
    cz-Newzhalflength,
    x1_boundary,
    x2_boundary-xhalflength,
    y1_boundary+yhalflength,
    y2_boundary,
    z1_boundary,
    z2_boundary-zhalflength,
    npmax,6);


    //#subquad8 so like 4 but pushed in y direction
    se2 = new Octree(newlevel,
    cx+Newxhalflength,
    cy+Newyhalflength,
    cz-Newzhalflength,
    x1_boundary+xhalflength,
    x2_boundary,
    y1_boundary+yhalflength,
    y2_boundary,
    z1_boundary,
    z2_boundary-zhalflength,
    npmax,7);
    }




    return;
    
}






int Octree::CollideParticles(double particle_list[][15],int n_S, int n_M,  double * Mass_adj_Precalc, double * R_adj_Precalc, int DimRange, double dr[][4],int ncollision, int points_array[][64],int NN)
{

    double dxij;
    double dyij;
    double dzij;
    double dij;
    double rai;
    double raj;

    //initiate RNG
    //std::default_random_engine generator;
   // std::uniform_real_distribution<double> distribution(0.0,1.0);


    std::random_device rd;

    //
    // Engines 
    //
    std::mt19937 e2(rd());
    //std::knuth_b e2(rd());
    //std::default_random_engine e2(rd()) ;

    //
    // Distribtuions
    //
    std::uniform_real_distribution<> dist(0, 1);



    
    //std::cout << points_array[10][3] << std::endl;


    //int ncollision=0;

    if (Leaf == 1)
    {
        //We are in leaf, so this node might have points in it we can collide
        if (pcounter > 0)

        {

            //#pragma omp parallel private(npi)
            //{
            int npi;
            #pragma omp parallel for
            for (npi=0; npi<pcounter-1; npi++){
            //#pragma omp parallel for shared(points_array)
            //#pragma omp parallel for default(none) shared(ukp1, uk) schedule(static,8)
                //cout << points[npi] << endl;




            //Here, I have pointer i and pointer_i_n_M, they are both used
            int pointer_i;
            //pointer_i = points[npi];


            pointer_i = points_array[npi][NN];


            //should it be without -n_M??? hmm
            int pointer_i_n_M = pointer_i-n_M;
            //if i is alive
            //det her, eller maaske gøre if >0? bare check om det er større end 0? Maaske easier computationally
            if ((int)(particle_list[pointer_i][10]+0.5) == 1)
            {

            
            //rai = particle_list[pointer_i][11];

            double ri;
            double rj;

            //double coolparam;
            //double dij;

            double dri [3];
            double pi [3];
            double drj [3];
            double pj [3];

            double drixdrj [3];
            double drixdrjnorm;
            double drixdrjunit [3];

            double pjpi [3];

            double dijlineline;

            //for t factors
            double dridotdrj;
            double dridotdri;
            double drjdotdrj;

            double n1 [3];
            double n2 [3];

            double dridotn2;
            double drjdotn1;

            double pjpidotn2;
            double pjpidotn1;
            double ti;
            double tj;


            // //Actually should be without -n_M i think. Because upda is already given DimRange.
            //Yes, update_positioncpp gives dr beginning from n_M, so the first n_M values are just 0,0,0,0 etc...
            dri[0] = dr[pointer_i][0];
            dri[1] = dr[pointer_i][1];
            dri[2] = dr[pointer_i][2]; 
            //dri[0] = dr[pointer_i_n_M][0];
            //dri[1] = dr[pointer_i_n_M][1];
            //dri[2] = dr[pointer_i_n_M][2]; 


            //position particle i
            pi[0] = particle_list[pointer_i][0];
            pi[1] = particle_list[pointer_i][1];
            pi[2] = particle_list[pointer_i][2];
            
            //radius particle i
            ri = particle_list[pointer_i][11];

            //cout<<"There are "<<omp_get_num_threads()<< " threads" << endl;


            //Declaring j here makes it private for each loop apparantly...
            //Maybe i also need to declare stuff above in a pragma clause??
            //Like the very first declarations???
            //int j; //https://stackoverflow.com/questions/31489664/openmp-nested-for-loop-gives-unexpected-result?noredirect=1&lq=1
            int npj;
                //#pragma omp parallel for shared(points_array)
                for (npj=npi+1; npj<pcounter; npj++)
                {
                int pointer_j;
                pointer_j = points_array[npj][NN];
                //pointer_j = points[npj];


                int pointer_j_n_M = pointer_j-n_M;
               if ((int)(particle_list[pointer_j][10]+0.5) == 1)
                {
                // // //Actually should be without -n_M i think. Because upda is already given DimRange.
                // dri[0] = dr[0][i-n_M];
                // dri[1] = dr[1][i-n_M];
                // dri[2] = dr[2][i-n_M]; 



                // //position particle i
                // pi[0] = particle_list[i][0];
                // pi[1] = particle_list[i][1];
                // pi[2] = particle_list[i][2];
                
                // //radius particle i
                // ri = particle_list[i][11];



                //radius particle j
                rj = particle_list[pointer_j][11];

                pj[0] = particle_list[pointer_j][0];
                pj[1] = particle_list[pointer_j][1];
                pj[2] = particle_list[pointer_j][2];


                //Here, for dipoles, just like for outer loop i, do WITHOUT -n_M
                drj[0] = dr[pointer_j][0];
                drj[1] = dr[pointer_j][1];
                drj[2] = dr[pointer_j][2];
                //drj[0] = dr[pointer_j_n_M][0];
                //drj[1] = dr[pointer_j_n_M][1];
                //drj[2] = dr[pointer_j_n_M][2];



                //#==================================================================
                //#Manual version distance
                //drixdrj = {dri[1]*drj[2]-drj[1]*dri[2],drj[0]*dri[2]-dri[0]*drj[2],dri[0]*drj[1]-drj[0]*dri[1]};
                drixdrj[0] = dri[1]*drj[2]-drj[1]*dri[2];
                drixdrj[1] = drj[0]*dri[2]-dri[0]*drj[2];
                drixdrj[2] = dri[0]*drj[1]-drj[0]*dri[1];


                //Norm af cross product... |axb|=|a||b|sin(v)...
                //Maybe that formula on the right is faster to compute, if i can get the angle between vectors from spherical units??
                // drixdrjnorm = (drixdrj[0]**2+drixdrj[1]**2+drixdrj[2]**2)**0.5
                drixdrjnorm = sqrt(pow(drixdrj[0],2)+pow(drixdrj[1],2)+pow(drixdrj[2],2));
                drixdrjunit[0] = drixdrj[0]/drixdrjnorm;
                drixdrjunit[1] = drixdrj[1]/drixdrjnorm;
                drixdrjunit[2] = drixdrj[2]/drixdrjnorm;

                pjpi[0] = pi[0]-pj[0];
                pjpi[1] = pi[1]-pj[1];
                pjpi[2] = pi[2]-pj[2];
                
                //#http://mathworld.wolfram.com/Line-LineDistance.html
                //#Bruge her i stedet for
                //#a = dri
                //#b = drj
                //#c = pi-pj
                //pjpi = pi-pj
                //dijlineline = abs(drixdrjunit[0]*pjpi[0]+drixdrjunit[1]*pjpi[1]+drixdrjunit[2]*pjpi[2])
                //abs value or sqrt(k^2)???
                dijlineline = std::abs (drixdrjunit[0]*pjpi[0]+drixdrjunit[1]*pjpi[1]+drixdrjunit[2]*pjpi[2]);

                
               
                //rj = particle_list[pointer_j][11];

                //dxij = particle_list[pointer_j][0]-particle_list[pointer_i][0];
                //dyij = particle_list[pointer_j][1]-particle_list[pointer_i][1];
                //dzij = particle_list[pointer_j][2]-particle_list[pointer_i][2];


                //dij = sqrt(dxij*dxij+dyij*dyij+dzij*dzij);

                
                //raj = particle_list[pointer_j][11];
                    //cout << dxij << endl;
                if (dijlineline <=ri+rj)
                {
                    


                    dridotdrj = dri[0]*drj[0]+dri[1]*drj[1]+dri[2]*drj[2];
                    
                    //#This, is almost calculated in upda. I should just return dri.dri, but atm it does np.sqrt(dri.dri) in upda! So, can reduce speed use here
                    //dridotdri = dri[0]*dri[0]+dri[1]*dri[1]+dri[2]*dri[2];
                    //drjdotdrj = drj[0]*drj[0]+drj[1]*drj[1]+drj[2]*drj[2];

                    dridotdri = dr[pointer_i][3];
                    drjdotdrj = dr[pointer_j][3];
                    

                    
                    n1[0] = dridotdri*drj[0] - dridotdrj*dri[0];
                    n1[1] = dridotdri*drj[1] - dridotdrj*dri[1];
                    n1[2] = dridotdri*drj[2] - dridotdrj*dri[2];


                    n2[0] = drjdotdrj*dri[0] - dridotdrj*drj[0];
                    n2[1] = drjdotdrj*dri[1] - dridotdrj*drj[1];
                    n2[2] = drjdotdrj*dri[2] - dridotdrj*drj[2];
                    

                    dridotn2 = dri[0]*n2[0]+dri[1]*n2[1]+dri[2]*n2[2];
                    drjdotn1 = drj[0]*n1[0]+drj[1]*n1[1]+drj[2]*n1[2];
                    


                    pjpidotn2 = pjpi[0]*n2[0]+pjpi[1]*n2[1]+pjpi[2]*n2[2];
                    pjpidotn1 = pjpi[0]*n1[0]+pjpi[1]*n1[1]+pjpi[2]*n1[2];
                    ti = -pjpidotn2/dridotn2;
                    tj = pjpidotn1/drjdotn1;

                    if ((0<=ti) && (ti<=1) && (0<=tj) && (tj<=1))
                        {



                        //cout << "Collision" << endl;

                        //cout << sqrt(dridotdri)/(2.0*ri+2.0*rj) << " " << sqrt(drjdotdrj)/(2.0*rj+2.0*ri) << endl;

                        double cosanglelines;
                        double ksafe = 0.000000000001;
                        //cosanglelines = dridotdri.drjdotdrj/(sqrt(dridotdri)*sqrt(drjdotdrj))
                        //cosanglelines = std::abs(dridotdrj)/(sqrt(dridotdri+ksafe)*sqrt(drjdotdrj+ksafe)+ksafe);
                        cosanglelines = std::abs(dridotdrj)/(sqrt(dridotdri*drjdotdrj));
                        //cosanglelines = dridotdrj/(sqrt(dridotdri*drjdotdrj));
                        //Angle can be between 0 to pi, hence, cosangle can be between -1 or 1


                        double anglelines;
                        anglelines = acos(cosanglelines+0.0000001);
                        anglelines = std::abs(anglelines);

                        //double Lcoll;

                       //Lcoll = 2.0*(ri+rj)/(sin(anglelines/2.0+0.001)+0.00001);

                        //double Picoll = Lcoll/(sqrt(dridotdri)+0.00000000000001);
                        //double Pjcoll = Lcoll/(sqrt(drjdotdrj)+0.00000000000001);

                        


                        //double Angleprob;
                        //Angleprob = 1.0-anglelines/(3.141592);

                        //cout << Picoll << " " << Pjcoll << " " << Lcoll << " " << cosanglelines << " " << Angleprob << " " << anglelines << endl;


                        double CollRNG;
                        CollRNG = 0.0;
                        //CollRNG = distribution(generator);
                        //CollRNG = 1;//
                        CollRNG = dist(e2);

                        //CollRNG = 0.5; //maybe try seed, or create array with random numbers, or something...
                        //but for now, just set it at E[X] = 0.5, average value should be decent

                        //cout << cosanglelines << " " << Angleprob << " " << CollRNG << endl;
                        

                        //if (CollRNG <= Angleprob)
                        if (cosanglelines <= CollRNG )
                        {

                        ncollision = ncollision+1;
                        
                        //cout << "CollRNG " << endl;
                        //cout << "We have collision" << endl;
                
                        /*particle_list[1,1] += 100;*/
                        //kill particle j
                        particle_list[pointer_j][10] = 0;
                        // //KilledIndex = np.append(KilledIndex,j)
                            
                            
                        // //choose center of mass for new position
                        particle_list[pointer_i][0] = (particle_list[pointer_i][13]*particle_list[pointer_i][0]+particle_list[pointer_j][13]*particle_list[pointer_j][0])/(particle_list[pointer_i][13]+particle_list[pointer_j][13]);
                        particle_list[pointer_i][1] = (particle_list[pointer_i][13]*particle_list[pointer_i][1]+particle_list[pointer_j][13]*particle_list[pointer_j][1])/(particle_list[pointer_i][13]+particle_list[pointer_j][13]);
                        particle_list[pointer_i][2] = (particle_list[pointer_i][13]*particle_list[pointer_i][2]+particle_list[pointer_j][13]*particle_list[pointer_j][2])/(particle_list[pointer_i][13]+particle_list[pointer_j][13]);

                        // #inelastic collision: get new direction
                        particle_list[pointer_i][3] = (particle_list[pointer_i][13]*particle_list[pointer_i][3]+particle_list[pointer_j][13]*particle_list[pointer_j][3])/(particle_list[pointer_i][13]+particle_list[pointer_j][13]);
                        particle_list[pointer_i][4] = (particle_list[pointer_i][13]*particle_list[pointer_i][4]+particle_list[pointer_j][13]*particle_list[pointer_j][4])/(particle_list[pointer_i][13]+particle_list[pointer_j][13]);
                        particle_list[pointer_i][5] = (particle_list[pointer_i][13]*particle_list[pointer_i][5]+particle_list[pointer_j][13]*particle_list[pointer_j][5])/(particle_list[pointer_i][13]+particle_list[pointer_j][13]);
            
            
                        // #increase number of adatoms in cluster
                         particle_list[pointer_i][14]=particle_list[pointer_i][14]+particle_list[pointer_j][14];
                            
                            
                        // #Mass and radius
                        particle_list[pointer_i][13] = Mass_adj_Precalc[(int)particle_list[pointer_i][14]-1];
                        particle_list[pointer_i][11] = R_adj_Precalc[(int)particle_list[pointer_i][14]-1];

                        break;


                        } //end CollRNG if


                        } //end ti and tj if

                } //end dijlineline<=ri+rj if





                }//if pointer_j is alive
            }//end npj loop


            }//if pointer_i is alive

            
        } //end npi loop
 

        return ncollision;
        } else {

            //Well, we are in leaf, but doesn't contain points, so we just return out from this node

            return ncollision;
        }




    } else
    {
            //We are not in leaf!
    ncollision = ne1->CollideParticles(particle_list,n_S,n_M,Mass_adj_Precalc,R_adj_Precalc,DimRange,dr,ncollision,points_array,0); 
    ncollision = nw1->CollideParticles(particle_list,n_S,n_M,Mass_adj_Precalc,R_adj_Precalc,DimRange,dr,ncollision,points_array,1);
    ncollision = sw1->CollideParticles(particle_list,n_S,n_M,Mass_adj_Precalc,R_adj_Precalc,DimRange,dr,ncollision,points_array,2);
    ncollision = se1->CollideParticles(particle_list,n_S,n_M,Mass_adj_Precalc,R_adj_Precalc,DimRange,dr,ncollision,points_array,3);
    ncollision = ne2->CollideParticles(particle_list,n_S,n_M,Mass_adj_Precalc,R_adj_Precalc,DimRange,dr,ncollision,points_array,4);
    ncollision = nw2->CollideParticles(particle_list,n_S,n_M,Mass_adj_Precalc,R_adj_Precalc,DimRange,dr,ncollision,points_array,5);
    ncollision = sw2->CollideParticles(particle_list,n_S,n_M,Mass_adj_Precalc,R_adj_Precalc,DimRange,dr,ncollision,points_array,6);
    ncollision = se2->CollideParticles(particle_list,n_S,n_M,Mass_adj_Precalc,R_adj_Precalc,DimRange,dr,ncollision,points_array,7);


    return ncollision;
    }



    return ncollision;

}






//class Employee { //rest of the class body here 
//public: 
//Employee *manager; };

//and objects could have )))

int Tryf(double particle_list[][15],Octree *p,int n_S, int n_M,  double * Mass_adj_Precalc, double * R_adj_Precalc, int DimRange, double dr[][4],int points_array[][64])
{
    //cout << "func" <<endl;

    p->split();

    int point_index;
    #pragma omp parallel for
    for (point_index=n_M; point_index<n_M+n_S; point_index++){
    p->Insert(point_index,particle_list,points_array,0);
    }
    
    int ncollision = 0;
    ncollision = p->CollideParticles(particle_list,n_S,n_M,Mass_adj_Precalc,R_adj_Precalc,DimRange,dr,ncollision,points_array,0);
    return ncollision;
}

/*
int main()
{


srand (time(NULL));


//Employee salesManager;
//Employee salesWorker;
//salesWorker.manager = salesManager;

double particle_list [20000][3];

for (int i=0; i<20000; i++)
{   

    //double rand1 = rand() % 10;
    //double rand2 = rand1/10;
    double r1 = ((double) rand() / (RAND_MAX));
    double r2 = ((double) rand() / (RAND_MAX));
    double r3 = ((double) rand() / (RAND_MAX));
    particle_list[i][0]=r1;
    particle_list[i][1]=r2;
    particle_list[i][2]=r3;

    //cout << rand2 << endl;
}



Octree* NewTree = new Octree(0,0.5,0.5,0.5,0,1,0,1,0,1);


//Node node1;



//for (int point_index=0; point_index<20000; point_index++){
    //NewTree->Insert(point_index,particle_list);
//}




cout << "Now we gonna find the particles!!" << endl;

//NewTree->CollideParticles(particle_list);

Tryf(particle_list,NewTree);

    return 0;
}
*/




extern "C" int CollisionAlgocppOctree(double particle_list[][15], int n_S, int n_M,  double * Mass_adj_Precalc, double * R_adj_Precalc, int DimRange, double dr[][4], double lx,int npmax,int points_array[][64])
{
    
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    //int points_array[1000][8];

    //cout << "HEllo" << endl;

    int ncollision;
    //nt npmax = 300;
    //cout << npmax << endl;
    //Construct first tree
    Octree* NewTree = new Octree(0,lx/2,lx/2,lx/2,0,lx,0,lx,0,lx,npmax,0);
    ncollision = Tryf(particle_list,NewTree,n_S,n_M,Mass_adj_Precalc,R_adj_Precalc,DimRange,dr,points_array);
    //Destruct Octree recursively
    NewTree->~Octree();

    //cout << "coll from octree " << abereturn << endl;
    
    //cout << lx << endl;


/*    //Write memory usage to file
    sysinfo (&memInfo);
    long long physMemUsed = memInfo.totalram - memInfo.freeram;
    //Multiply in next statement to avoid int overflow on right hand side...
    physMemUsed *= memInfo.mem_unit;



    std::ofstream CppProcess;
    CppProcess.open("CppProcess.dat",std::ios::in | std::ios::app);
    CppProcess << getpid() << " " ;
    CppProcess << getppid() << " " ;
    CppProcess << physMemUsed;
    CppProcess << "\n";
    CppProcess.close();

*/

    return  ncollision;
    }

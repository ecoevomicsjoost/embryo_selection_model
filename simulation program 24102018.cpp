//// embryo selection program. Simulates de novo methylation after demethylation. 
//// Joost van den Heuvel, joost.vandenheuvel@wur.nl
//// last edits 24102018

#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using namespace std;

class TFBS
{
      public:
      int loci[10]; // every TFBS has 10 loci in the standard model
      };

class indi
{
      
      public:
             
             void firstback(TFBS& a){firstvec.push_back(a);}
             void secondback(TFBS& a){secondvec.push_back(a);}
                          
             float firsthow(int i, int j){return firstvec[i].loci[j];}
             float secondhow(int i, int j){return secondvec[i].loci[j];}
             
              float firstgive(int i, int j, int k){return firstvec[i].loci[j]=k;}
              float secondgive(int i, int j, int k){return secondvec[i].loci[j]=k;}

                void Switch(){ firstvec.swap(secondvec);}
          
        
           
             
                     vector<TFBS> firstvec;
                     vector<TFBS> secondvec;
             
             };



int main()
{      
     srand ((unsigned)time (NULL));
     
     indi embryo;
   
    float Q[10000];
    float epitot;

     int I = 50;    // number of cells
     int i;
     for (i=0; i<I; i++)
     {
     TFBS ins;
     embryo.firstback(ins);
     embryo.secondback(ins);
     }
    int j; 
    float help1, help2; 
     float sum;
     
     float prob;
  float rate=0.1;   

     
     ofstream outf ("methylation.txt");  // output file, will contain all relevant information for R script

     int ind;
float locis[10]={};

int n_l = 10;

float TF = 10;
float H = 25;
float mean;
float sumprob;
int f;
int m, w;
int p;

int r,q;
float a;     
 int T;       
   float qq;
   
   I = 50;
   n_l = 10;
   H = 10;
   rate =0.02;
   qq = 1.2;
   

w = 0;   

/// depending on how many simlations one wants to perform, q can be altered
   
   for (q=0; q<1; q++)
   {
     
    for (j=0; j<n_l; j++){

        Q[j] = qq;
}    
  


/// TF activity levels for the 75 TFBSs
    for (f=0; f<76; f++)
    {
  
    
  if (r==0){  TF = 5*(f+1);}

    
     for (i=0; i<I; i++)/// initially all values are zero // this makes sure it is so
     {
     for (j=0; j<n_l; j++)
     {
     embryo.firstgive(i,j,0);
     embryo.secondgive(i,j,0);
     }}

    cout << q<<'\t'<<m<<'\t'<<w<<'\t'<<  r<<'\t'<< f<<endl;
   
     
     for (ind=0; ind<100; ind++) //number of individuals
     {

 
     for (i=0; i<I; i++) // number of cells
     {
     for (j=0; j<n_l; j++)
     {
     embryo.firstgive(i,j,0);
    embryo.secondgive(i,j,0);
     }}
 
 
     for (T=0; T<1001; T++)
     {
     sumprob = 0;     
     
     for (j=0; j<n_l; j++){locis[j] = 0;}
     
     for (i=0; i<I; i++)
     {
     
     ///// first chromosome
     
     sum = 0;
     epitot =1;
     for (j=0; j<n_l; j++)
     {
     sum = sum + embryo.firsthow(i,j);
     locis[j] = locis[j]+embryo.firsthow(i,j);
     epitot = epitot*pow(Q[j],embryo.firsthow(i,j));
     }
     
     prob = pow(TF,2)/(pow(TF,2)+pow((H*epitot),2));
     sumprob = sumprob+prob;
 
     for (j=0; j<n_l; j++)
     {
     a = rand();
     a = a / (RAND_MAX+1);
     if (embryo.firsthow(i,j)==0){
     if (a<(1-prob)){embryo.firstgive(i,j,1);}
     }
     else 
     {
     if (a<rate){embryo.firstgive(i,j,0);}
     }
     }

    ////////second chromosome since the emrbyo is diploid
    
     sum = 0;
     epitot =1;
     for (j=0; j<n_l; j++)
     {
     sum = sum + embryo.secondhow(i,j);
     locis[j] = locis[j]+embryo.secondhow(i,j);
     epitot = epitot*pow(Q[j],embryo.secondhow(i,j));
     }
     
     prob = pow(TF,2)/(pow(TF,2)+pow((H*epitot),2));
     sumprob = sumprob+prob;
 
     for (j=0; j<n_l; j++)
     {
     a = rand();
     a = a / (RAND_MAX+1);
     if (embryo.secondhow(i,j)==0){
     if (a<(1-prob)){embryo.secondgive(i,j,1);}
     }
     else 
     {
     if (a<rate){embryo.secondgive(i,j,0);}
     }
     }
}  // i loop
if (T==999)
{
for (j=0; j<n_l; j++){outf<<H<<'\t'<< w<<'t'<<f<<'\t'<<r<<'\t'<< ind<<'\t'<< j<<'\t'<<locis[j]/(2*I)<<'\t'<< sumprob<<endl;}
}
} // T loop

  



} // ind
}} 
cout << "Done "<<endl;
cin.get();
cin.get();
     return 0;
     }

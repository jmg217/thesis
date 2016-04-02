#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>

//this function returns the payoff value
double payoff(double x, double k){
double h;
h=k-x;
	if(h<0){
	h=0;
	}

return h;
}

//this function returns the high bias mesh price
double MeshEstimator(double strike, double r, double delta_t, int b, double m,  std::vector< std::vector<double> >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V){

double H; //payoff variable 
double C; //continuation value variable
double sum, V_0;

// temp vector in Estimator loop
std::vector< double > tempvec;


//Mesh Estimator loop

for(int i=0; i<m; i++){
tempvec.clear();

	for(int j=0; j<b; j++){
		if(i==0){
		H=payoff(X[(m-1)-i][j], strike)*exp(-r*delta_t*(m-i));
		tempvec.push_back(H);
		}
	
		else{
//continue to develope continuation vale			
			sum=0;
			for(int k=0; k<b; k++){
			sum+=(W[(m-i)][k][j])*V[i-1][k]; //m-i when i=1 is 10-1=9.when i=9 m-i=1. we get V_0 separately by using W[0][k][j]	
		
			}

			C=(1/((double)b))*sum; //continuation value
			H=payoff(X[(m-1)-i][j], strike)*exp(-r*delta_t*(m-i));
		
			if(H>=C){
			tempvec.push_back(H);
			
			}

			else{
			tempvec.push_back(C);
			
			}	
		}
	}
V.push_back(tempvec); //high bias option value matrix
}


//sum over option values at the second time step
sum=0;
for(int k=0; k<b; k++){
sum+=V[m-1][k];        
}
//this is the high bias option value at time 0
V_0=(1/((double)b))*sum;


return V_0;
}

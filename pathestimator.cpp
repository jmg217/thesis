#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>


//declare transition density function. This is contained in meshgen.cpp
double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t);
//declare the box muller function. This is contained in meshgen.cpp
double boxmuller();
//this returns the payoff value
double Payoff(double x, double k){
double h;
h=k-x;
        if(h<0){
        h=0;
        }

return h;

}


//This function returns a low bias price of an option.
double PathEstimator(double strike, double r, double delta_t, int b, double m, double sigma, double delta, double X0, std::vector< std::vector<double> >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V){

double v_0, S_i, Z, C, H, sum, weight, w_s, sum_Z;

//Simulated path for sub optimal stoppage
std::vector<double > S;
//temp vector of S weights
std::vector<double > tempvec;
//weights matrix for S
std::vector< std::vector<double> > S_Weights; 

for(int i=0; i<m; i++){
	if(i==0){
	Z=boxmuller();
	S_i=X0 * (exp((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));//the second value in the simulated path
	}

	else{
	Z=boxmuller();
	S_i=S[i-1]*(exp((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));//the simulate path values
	}

S.push_back(S_i);//simulated path is stored in this vector 
}

//this for-loop generates the S weights 
for(int t=0; t<(m-1); t++){
tempvec.clear();
	for(int h=0; h<b; h++){   //h=k
	sum=0;
	w_s=density(S[t], X[t+1][h], sigma, r, delta, delta_t);
		for(int g=0; g<b; g++){   //g=l
		sum+=(1/((double)b))*density(X[t][g], X[t+1][h], sigma, r, delta, delta_t);
		}
	w_s=w_s/sum;	
	tempvec.push_back(w_s);
	}

S_Weights.push_back(tempvec); //vector storing S weights
}

double con_val; //continuation value variable
//sub optimal path loop
for(int i=0; i<m; i++){
	sum=0;
	if(i==m-1){
	C=0;//continuation value at the last time step
	}
	
	else{
		for(int k=0; k<b; k++){	

			weight=S_Weights[i][k];
			con_val=V[(m-1)-i-1][k];
		
			sum+=weight*con_val;			
		}
	C=(1/(double)b)*sum; //continuation value
	}	
	

H=Payoff(S[i], strike)*exp(-r*delta_t*((i+1)));

//check if continuation value is greater then the immediate payoff
	if(H>=C || i==m-1){
	v_0=H; 
	break;
	}

}

return v_0;

}

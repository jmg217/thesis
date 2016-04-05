#include<iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sstream>

int main(){

std::ifstream setting( "settings.txt" );
std::string line;
std::vector<std::string> settings;
int linenumber=0;
while(std::getline( setting, line))
    {
          if(linenumber%2==1)
              settings.push_back(line);
          linenumber++;
    }
setting.close();

double double_num;
int integer; 

std::vector < double > X0;
std::vector < double > delta;
std::vector <double> sigma;

std::istringstream ss(settings[0]);
std::string token;
while(std::getline(ss, token, ',')) 
    {
          X0.push_back(atof(token.c_str()));
    }

double T = atof(settings[1].c_str());
double m = atof(settings[2].c_str());
double delta_t=T/m;
int Rn;
double v_0, V_0, Z, Xi, Xj, w, wdenominator, v_sum, sum_Z=0, vtotal_sum=0, Vtotal_sum=0;
double r= atof(settings[3].c_str());

std::istringstream ss2(settings[4]);
while(std::getline(ss2, token, ','))
    {
        delta.push_back(atof(token.c_str()));
    }

std::istringstream ss3(settings[5]);
while(std::getline(ss3, token, ','))
    {
        sigma.push_back(atof(token.c_str()));
    }

double Path_estimator_iterations=atof(settings[6].c_str());
double strike=atof(settings[7].c_str());
int b=atoi(settings[8].c_str());
int N=atoi(settings[9].c_str());
double quantile=atof(settings[10].c_str());
int num_assets=atof(settings[11].c_str());

if(X0.size() != num_assets || sigma.size() != num_assets || delta.size() !=num_assets){
          std::cout<<"Either the starting price, volatility, number of assets or dividend yield was not specified for all assets"<<std::endl;
           exit (EXIT_FAILURE);
       }

std::cout<<"The parameters of this simulation are:"<<std::endl;
       //Print these values to screen 
for(integer=0; integer<X0.size(); integer++){
std::cout<<"Starting Price="<<X0[integer]<<std::endl;
}

std::cout<<"Time to expiry="<<T<<"\n"<<"Number of time steps="<<m<<"\n"<<"interest rate="<<r<<std::endl;

for(integer=0; integer<sigma.size(); integer++){
    std::cout<<"volatility="<<sigma[integer]<<std::endl;
}

for(integer=0; integer<delta.size(); integer++){
    std::cout<<"dividend yield="<<delta[integer]<<std::endl;
}
std::cout<<"number of iterations over path estimator="<<Path_estimator_iterations<<"\n"<<"strike  price="<<strike<<"\n"<<"number of nodes per time step="<<b<<"\n"<<"number mesh generations="<<N<<"\n"<<"Number of Assets="<<num_assets<<std::endl; 


}

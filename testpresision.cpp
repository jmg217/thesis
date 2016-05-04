#include <iostream>                                                                                                                                                              
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm> 

double UniRandom(double b) //This function returns uniformally distributed variables in the interval (0,b-1)
  {
  double rn;
  rn=((double)rand()/(double)RAND_MAX)*((double)(b-1));
  rn=round(rn);
  return rn;
  }



int main(){
srand((unsigned)time(NULL));
long double ss, total;
long double sss;
ss=0.0001000100010000001;
total=100;
//sss <<std::setprecision(30)<< total+ss<<std::endl;;

std::vector<double> m;

for(int i=0; i<10; i++){
m.push_back(UniRandom(100));
}

m.push_back(0.0000000000005);
m.push_back(0.0000000000001);
for(int j=0; j<10; j++){
std::cout<<m[j]<<std::endl;
}
std::sort(m.begin(), m.end());


for(int k=0; k<10; k++){
  std::cout<<"2="<<m[k]<<std::endl;
 }
double p=0;
for(int y=0; y<10; y++)
{
p+=m[y];
std::cout<<"sum="<<p<<std::endl;
}

//std::cout<< sss<<std::endl;
//ss >> total;

return 0;
}

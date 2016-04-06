#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>

void print_high_payoff(int b, double m, std::vector<std::vector< std::vector<double> > >& X, std::vector< std::vector<double> >& V, std::vector<double>& asset_amount ){
double x=0;
std::ofstream outFile("highpayoff.txt", std::ios_base::app | std::ios_base::out);

for (int t=0; t<m; t++){

	for (int i=0; i<b; i++){
	x=0;
		for(int tt=0; tt<asset_amount.size(); tt++){
		x+=asset_amount[tt]*X[m-1-t][i][tt];
		}
	outFile << m-t <<"\t"<< x <<"\t"<< V[t][i]<< std::endl;
	}
}

outFile.close();

}



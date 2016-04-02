#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>

#define PI 3.14159265358979323846

//This function is contained within meshestimator.cpp
double MeshEstimator(double strike, double r, double delta_t, int b, double m,  std::vector< std::vector<double> >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V);

//This function is contained within patestimator.cpp.
double PathEstimator(double strike, double r, double delta_t, int b, double m, double sigma, double delta, double X0, std::vector< std::vector<double> >& X, std::vector< std::vector< std::vector<double> > >& W, std::vector< std::vector<double> >& V);



double round( double value )//This function rounds the double sent to it to the nearest integer.
{
  return floor( value + 0.5 );
}



double boxmuller()//This function uses the Box Muller algorithm to return standard normally distributed varibles.
{
double U1, U2, R, Y, Z1, Z2;

	for(int i=0; i<100000; i++){

		U1=((double)rand()/(double)RAND_MAX);
		U2=((double)rand()/(double)RAND_MAX);

		U1=2*U1-1;
		U2=2*U2-1;

		R=U1*U1+U2*U2;

		if(R<=1){
			Y=sqrt((-2*log(R))/R);
			Z1=U1*Y;
			Z2=U2*Y;

		}
	}
return Z1;
}

int UniRandom(int b) //This function returns uniformally distributed variables in the interval (0,b-1)
{
double rn;
rn=((double)rand()/(double)RAND_MAX)*((double)(b-1));
rn=round(rn);
return rn;
}

double density(double Xold, double  Xnew, double sigma, double r, double delta, double delta_t)// This function returns the transition density between node values.
{
double f, x;
x=(1/(sigma*sqrt(delta_t)))*(log(Xnew)-log(Xold)-(r-delta-0.5*sigma*sigma)*delta_t);

f= (1/(sigma*sqrt(delta_t)*Xnew))*(1/(2*PI))*exp(-0.5*x*x); // this is the transition density

return f;
}


int main(){

srand((unsigned)time(NULL));

//Now we read in the parameters from settings.txt
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
double X0=atof(settings[0].c_str());
double T = atof(settings[1].c_str());
double m = atof(settings[2].c_str());
double delta_t=T/m;
int Rn;
double v_0, V_0, Z, Xi, Xj, w, wdenominator, v_sum, sum_Z=0, vtotal_sum=0, Vtotal_sum=0;
double r= atof(settings[3].c_str());
double delta=atof(settings[4].c_str());
double sigma=atof(settings[5].c_str());
double Path_estimator_iterations=atof(settings[6].c_str());
double strike=atof(settings[7].c_str());
int b=atoi(settings[8].c_str());
int N=atoi(settings[9].c_str());
double quantile=atof(settings[10].c_str());
int num_assets=atof(settings[11].c_str());
//Print these values to screen 
std::cout<<"The parameters of this simulation are:"<<"\n"<<"Starting Price="<<X0<<"\n"<<"Time to expiry="<<T<<"\n"<<"Number of time steps="<<m<<"\n"<<"interest rate="<<r<<"\n"<<"volatility="<<sigma<<"\n"<<"dividend yield="<<delta<<"\n"<<"number of iterations over path estimator="<<Path_estimator_iterations<<"\n"<<"strike price="<<strike<<"\n"<<"number of nodes per time step="<<b<<"\n"<<"number mesh generations="<<N<<"\n"<<"Number of Assets="<<num_assets<<std::endl;


//Mesh matrix
std::vector< std::vector<double> > X;
//WEIGHTS 3-dimensional matrix for step 1 and beyond
std::vector< std::vector< std::vector<double> > > W;
//temp vector in MeshGen for-loop
std::vector< double > myvector;
//2 d temp vector in WeightsGen for-loop
std::vector< std::vector<double> > dim2temp;
//1 d vector in Weightsgen for-loop
std:: vector<double> dim1temp;
//mesh estimator high bias 2-d matrix
std::vector< std::vector<double> > V;
//V values from each iteration over meshes
std::vector< double > Vvector;
//v values from each iteration over meshes
std::vector< double > vvector;

//for-loop over different meshes
for(int iterator=0; iterator<N; iterator++){
X.clear();
W.clear();
V.clear();
//for-loop to generate the mesh
for(int i=0; i<m; i++){


	myvector.clear();

	if(i==0){

	for(int l=0; l<b; l++){
		Z=boxmuller();//standard normally distributed variable 
		Xi=X0 * (exp ((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));//node value at the second time step
		myvector.push_back(Xi);	//store the value in a temp vector
	}
	}

	if(i>0){
	
	for(int j=0; j<b; j++){
		Z=boxmuller();
		Rn=UniRandom(b); 
		Xi=X[i-1][Rn];
		Xj=Xi * (exp ((r-delta-0.5*sigma*sigma)*delta_t + sigma*sqrt(delta_t)*Z));
		myvector.push_back(Xj);
	}
	}

X.push_back(myvector);

}


//Weights generation for-loop 
//NOTE: W^i_(j,k) IS REPRESENTED AT W[i][k][j] where k is at step i+1 and j is at step i.
for(int I=0; I<m; I++){

dim2temp.clear();//temporary vector
	
	if(I==0){
		for(int k=0; k<b; k++){
        	dim1temp.clear();
			w=1;// all weights from the starting node are equal to 1
		dim1temp.push_back(w);
		dim2temp.push_back(dim1temp);
		}
	}


	if(I>0){
	
		for(int k=0; k<b; k++){	
		dim1temp.clear();
		wdenominator=0;
	
			for(int j=0; j<b; j++){
			w=density(X[I-1][j], X[I][k], sigma, r, delta, delta_t);//step 1 in X is X[0] step 1 in W is W[1]	
			dim1temp.push_back(w);
			wdenominator+=w; //this generates the denominator value in the weights formula
			}
			
			//devide each element by the denominator
			for(int t=0; t<b; t++){
			dim1temp[t]=(((double)b)*dim1temp[t])/wdenominator;
			}		
		dim2temp.push_back(dim1temp); //dim1 is full therefore we add it onto dim2 vector
		}	
	
	}

W.push_back(dim2temp); //mesh weights matrix
}

//some weights checking
double check=0;
//check all the weights from X0 are 1
for(int e=0; e<b; e++){
if(W[0][e][0]!=1){
std::cout<<"there is an error with the weights. check that W[0][k][0]'s =1"<<std::endl;
}
}
//check that the weights going into a node sum to 1
for(int q=1; q<m; q++){ 
	for(int a=0; a<b; a++){
		check=0;
		for(int E=0; E<b; E++){
			check+=W[q][a][E];

		}
	}
}

V_0=MeshEstimator(strike, r, delta_t, b, m, X, W, V);//high bias option price
Vvector.push_back(V_0);//vector containing high bias option prices
Vtotal_sum+=V_0;

std::cout<<"High Bias price (V_0) for mesh iteration "<<iterator<<" is "<<V_0<<std::endl;

//average over path estimators
v_sum=0;
for(int f=0; f<Path_estimator_iterations; f++){
v_sum+=PathEstimator(strike, r, delta_t, b,  m, sigma, delta, X0, X, W, V);
}

v_0=(1/double(Path_estimator_iterations))*v_sum;
vvector.push_back(v_0);
vtotal_sum+=v_0;

std::cout<<"Low Bias price (v_0) for mesh iteration "<<iterator<<" is "<<v_0<<std::endl;
}//this is the end of the loop over the whole process.

//Calculate V(N) and v(N)
V_0=(1/double(N))*Vtotal_sum;
v_0=(1/double(N))*vtotal_sum;

//calculate errors
double std_div_V=0, std_div_v=0, squaresumV=0, squaresumv=0, Verror=0, verror=0;

for(int h=0; h<N; h++){
squaresumV+=(Vvector[h]-V_0)*(Vvector[h]-V_0);
squaresumv+=(vvector[h]-v_0)*(vvector[h]-v_0);
}
std_div_V=sqrt((1/double(N))*squaresumV); //standard deviation of V
std_div_v=sqrt((1/double(N))*squaresumv); //standard deviation of v

Verror=quantile*std_div_V*(1/sqrt(double(N)));
verror=quantile*std_div_v*(1/sqrt(double(N)));

std::cout<<"V(N)_0="<<V_0<<"\t"<<"V error="<<Verror<<std::endl;
std::cout<<"v(N)_0="<<v_0<<"\t"<<"v error="<<verror<<std::endl;


std::ofstream outFile("results.txt", std::ios_base::app | std::ios_base::out);

outFile << N <<"\t"<< b <<"\t"<< Path_estimator_iterations<<"\t"<< V_0 <<"\t"<< v_0 <<"\t"<< Verror+V_0 <<"\t"<< v_0-verror << std::endl;

outFile.close();


return 0;


}



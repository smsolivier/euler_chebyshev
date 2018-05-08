#include "DataObjects.H"
#include "Writer.H"

using namespace std; 

int main(int argc, char* argv[]) {
#ifndef ZERO 
	cout << "WARNING: unstable if ZERO not defined" << endl; 
#endif
	int N = 8; 
	if (argc > 1) N = atoi(argv[1]); 
	array<int,DIM> dims = {N,N,N}; 

	double T = 1; // end time 
	double K = 0.001; // time step 
	int Nt = T/K; // number of time steps 
	int NSAVES = 300; 
	int mod = Nt/NSAVES; 

	// setup writer 
	Writer writer; 
	writer.setFreq(mod); 

	// initialize variables 
	// vorticities 
	Vector omega0(dims); 
	Vector omega1(dims); 
	Vector omega(dims); 

	// velocities 
	Vector V0(dims, true); 
	Vector V1(dims); 
	Vector V(dims); 
	Vector Vhalf(dims); 

	// pressure 
	Scalar P0(dims); 
	Scalar P1(dims); 
	Scalar Pi(dims); 

	// divergence 
	Scalar div(dims); 

	// add variables to writer 
	writer.add(omega, "vorticity"); 
	writer.add(V, "velocity"); 
	writer.add(Pi, "Pi"); 
	writer.add(div, "div"); 

	// print out memory requirement 
	div.memory(); 

	for (int i=0; i<V0.size(); i++) {
		V0[0][i] = 1.; 
	}
	V0.forward(); 

	// store cross products 
	Vector cross0 = K*V0.cross(omega0); 

	// first time step with forward euler 
	omega0 = V0.curl(); 
	P0 = (V0 + cross0).divergence(); 
	P0.invert_laplacian(); 
	V1 = V0 + cross0 - K*P0.gradient(); 

	// setup for second time step 
	omega1 = V1.curl(); 

	// store K/2*(3 V1xomega1 - V0xomega0)
	Vector AB2_cross(dims); 

	// store K*V1xomega1 
	Vector cross1(dims); 

	// for (int t=1; t<Nt+1; t++) {

	// 	// cross products 
	// 	cross1 = V1.cross(omega1); 
	// 	AB2_cross = K/2*(3.*cross1 - cross0); 

	// 	// pressure 
	// 	Pi = 
	// }
}
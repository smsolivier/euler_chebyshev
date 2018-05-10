#include "DataObjects.H"
#include "Greens.H"
#include "Writer.H"
#include "Inverter.H"

using namespace std; 

double gaussian(double x, double y, double xc, double yc, double width) {
	return exp(-pow(x - xc,2)/width)*exp(-pow(y-yc,2)/width); 
}

int main(int argc, char* argv[]) {
#ifndef ZERO 
	WARNING("unstable if zero not defined"); 
#endif
#ifdef PC 
	cout << "paper cutter on" << endl; 
#endif 
#ifdef TAU 
	cout << "tau's on" << endl; 
#endif
	int N = 16; 
	if (argc > 1) N = atoi(argv[1]); 
	array<int,DIM> dims = {N,N,N}; 

	double T = 6; // end time 
	double K = 0.0001; // time step 
	int Nt = T/K; // number of time steps 
	int NSAVES = 300; 
	int mod = Nt/NSAVES; 

	// setup writer 
	Writer writer; 
	writer.setFreq(mod); 

	// initialize variables 
	// vorticities 
	Vector omega0(dims, true); 
	Vector omega1(dims); 
	Vector omega(dims); 

	// velocities 
	Vector V0(dims, false); 
	Vector V1(dims); 
	Vector V(dims); 
	Vector Vhalf(dims); 
	Vector V34(dims); 

	// pressure 
	Scalar P0(dims); 
	Scalar P1(dims); 
	Scalar Pi(dims); 
	Scalar Pih(dims); 

	// divergence 
	Scalar div(dims); 

	// store K/2*(3 V1xomega1 - V0xomega0)
	Vector AB2_cross(dims); 

	// store K*V_i x omega_i
	Vector cross1(dims); 
	Vector cross0(dims); 

	// setup green's functions 
	Greens greens(dims); 

	// streamline for initial velocity 
	Scalar psi(dims); 
	Vector grad(dims); 

	// add variables to writer 
	writer.add(omega, "vorticity"); 
	writer.add(V, "velocity"); 
	writer.add(V[0], "Vx"); 
	writer.add(V[1], "Vy"); 
	writer.add(V[2], "Vz"); 
	writer.add(Pih, "Pih"); 
	writer.add(div, "div"); 

	// print out memory requirement 
	div.memory(); 

	// set initial conditions 
	double width = M_PI/8; 
	double mag = 25; 
	double dist = M_PI/2; 
	for (int i=0; i<dims[0]; i++) {
		double x = 2*M_PI*i/dims[0]; 
		for (int j=0; j<dims[1]; j++) {
			double y = 2*M_PI*j/dims[1]; 
			for (int k=0; k<dims[2]; k++) {
				double z = cos(k*M_PI/(dims[2]-1)); 
				omega0[2](i,j,k) = mag*gaussian(x,y,M_PI-dist,M_PI,width) + 
					mag*gaussian(x,y,M_PI+dist, M_PI,width); 
			}
		}
	}
	omega0.forward(); 

	for (int i=0; i<dims[0]; i++) {
		double m = psi.freq(i,0); 
		for (int j=0; j<dims[1]; j++) {
			double n = psi.freq(j,1); 
			for (int k=0; k<dims[2]; k++) {
				if (n + m != 0) {
					psi(i,j,k) = -1./(m*m + n*n)*omega0[2](i,j,k); 
				}
			}
		}
	}
	grad = psi.gradient(); 
	V0[0] = -1.*grad[1]; 
	V0[1] = grad[0]; 
	
	// store cross products 
	cross0 = K*V0.cross(omega0); 

	// first time step with forward euler 
	omega0 = V0.curl(); 

	Vhalf = V0 + K*V0.cross(omega0); 
	Pih = Vhalf.divergence();
	Pih.invert_theta(Vhalf); 
	Vhalf = V0 + K*V0.cross(omega0) - Pih.gradient(); 
	greens.tau(Vhalf, V1); 

	div = V1.divergence(); 
	div.inverse(); 
	cout << "max div = " << div.max() << endl; 

	// setup for second time step 
	omega1 = V1.curl(); 

	for (int t=1; t<Nt+1; t++) {

		// // cross products 
		// cross1 = K*V1.cross(omega1); 
		// AB2_cross = .5*(3.*cross1 - cross0);

		// AB2 step 
		Vhalf = V1 + K/2*(3.*V1.cross(omega1) - V0.cross(omega0)); 

		// pressure term 
		Pih = Vhalf.divergence(); 
		Pih.invert_theta(Vhalf); 

		// green's function step 
		V34 = Vhalf - Pih.gradient(); 
		greens.tau(V34, V); 

		// setup for next time step 
		omega = V.curl(); 

		// check divergence is ok 
		div = V.divergence(); 
		div.inverse(); 

		// save histories 
		V0 = V1; 
		V1 = V; 
		omega0 = omega1; 
		omega1 = omega; 
		// cross0 = cross1; 

		// write to VTK 
		writer.write(); 

		double div_max = div.max(); 
		CHECK(div_max < 1e-3, "divergence greater than tolerance"); 

		// for (int k=0; k<dims[2]; k++) {
		// 	CHECK(V[2](0,0,k) == 0., "Vz must be zero"); 
		// }

		cout << t*K/T << "\r"; 
		cout.flush();
		// cout << div.max() << endl; 
	}
	cout << "final div = " << div.max() << endl; 
}
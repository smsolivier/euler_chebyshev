#include "DataObjects.H"
#include "Writer.H"
#include "Greens.H"
#include "Inverter.H"
#include "Sparse"

using namespace std; 

double gaussian(double x, double y, double xc, double yc, double width) {
	return exp(-pow(x - xc,2)/width)*exp(-pow(y-yc,2)/width); 
}

int main(int argc, char* argv[]) {
	int N = 8; 
	if (argc > 1) N = atoi(argv[1]); 
	array<int,DIM> dims = {N, N, N}; 

	double K = .001; 
	Writer writer; 

	Vector V0(dims, false); 
	Vector Vhalf(dims); 
	Vector V34(dims); 
	Vector V1(dims); 
	Vector omega0(dims, true); 
	Scalar div(dims); 
	Scalar Pih(dims); 

	Greens greens(dims); 

	// set initial conditions 
	double width = M_PI/8; 
	double mag = 1; 
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
	Scalar psi(dims); 
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
	Vector grad = psi.gradient(); 
	V0[0] = -1.*grad[1]; 
	V0[1] = grad[0]; 

	omega0 = V0.curl(); 

	Vhalf = V0 + K*V0.cross(omega0); 
	Pih = Vhalf.divergence(); 
	// Pih.invert_theta(Vhalf); 
	Inverter::instance().invert(Pih, Vhalf); 
	// for (int i=0; i<dims[0]; i++) {
	// 	for (int j=0; j<dims[1]; j++) {
	// 		Pih(i,j,dims[2]-1) = 0; 
	// 		Pih(i,j,dims[2]-2) = 0; 
	// 		for (int k=0; k<dims[2]; k++) {
	// 			Pih(i,j,dims[2]-1) += Vhalf[2](i,j,k); 
	// 			Pih(i,j,dims[2]-2) += Vhalf[2](i,j,k)*pow(-1,k); 
	// 		}
	// 	}
	// }
	// Inverter::instance().invert(Pih); 
	V34 = Vhalf - Pih.gradient(); 
	greens.tau(V34, V1); 

	div = V1.divergence(); 
	div.inverse(); 
	cout << "max div = " << div.max() << endl; 

	writer.add(div, "div"); 
	writer.add(Pih, "Pih"); 
	writer.add(V1, "V1"); 
	writer.add(V1[0], "V1_x"); 
	writer.add(V1[1], "V1_y"); 
	writer.add(V1[2], "V1_z"); 
	writer.write(); 
}
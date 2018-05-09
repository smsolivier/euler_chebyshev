#include "DataObjects.H"
#include "Writer.H"
#include "Inverter.H"

using namespace std; 

int main(int argc, char* argv[]) {
	int N = 8; 
	int BC = 1; 
	if (argc>1) N = atoi(argv[1]); 
	if (argc>2) BC = atoi(argv[2]); 
	array<int,DIM> dims = {N, N, N}; 

	Writer writer; 

	Scalar g(dims, true); 
	Scalar f(dims, false); 
	writer.add(g, "g"); 

	// set rhs 
	for (int i=0; i<dims[0]; i++) {
		double x = i*2*M_PI/dims[0]; 
		for (int j=0; j<dims[1]; j++) {
			double y = j*2*M_PI/dims[1]; 
			for (int k=0; k<dims[2]; k++) {
				double z = cos(k*M_PI/(dims[2]-1)); 
				g(i,j,k) = exp(-pow(x-M_PI,2))*exp(-pow(y-M_PI,2)); 
				// g(i,j,k) = 1.; 
			}
		}
	}
	g.forward();

	// set boundary conditions for inverter  
	for (int i=0; i<dims[0]; i++) {
		for (int j=0; j<dims[1]; j++) {
			// if (i==0 && j==0) 
			// 	g(i,j,dims[2]-1) = 1.; 
			// else 
			// 	g(i,j,dims[2]-1) = 0; 
			g(i,j,dims[2]-1) = 0; 
			g(i,j,dims[2]-2) = 0; 
		}
	}

	Inverter::instance().init(dims, BC); 
	Inverter::instance().invert(g); 

	writer.write(); 
}
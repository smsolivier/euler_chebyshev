#include "DataObjects.H"
#include "Inverter.H"
#include "Writer.H"

using namespace std; 

int main(int argc, char* argv[]) {
	int N = 8; 
	if (argc>1) N = atoi(argv[1]); 
	array<int,DIM> dims = {N, N, N}; 

	Inverter inverter(dims); 
	Writer writer; 

	Scalar g(dims, true); 
	Scalar f(dims, false); 
	writer.add(g, "g"); 

	for (int i=0; i<dims[0]; i++) {
		double x = i*2*M_PI/dims[0]; 
		for (int j=0; j<dims[1]; j++) {
			double y = j*2*M_PI/dims[1]; 
			for (int k=0; k<dims[2]; k++) {
				double z = cos(k*M_PI/(dims[2]-1)); 
				g(i,j,k) = exp(-pow(x-M_PI,2))*exp(-pow(y-M_PI,2)); 
			}
		}
	}
	g.forward(); 

	inverter(g); 

	writer.write(); 
}
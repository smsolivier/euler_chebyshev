#include "DataObjects.H"
#include "Writer.H"

using namespace std; 

int main() {
	int N = 16; 
	array<int,DIM> dims = {N, N, N}; 

	Writer writer; 

	Scalar u(dims, true); 
	for (int i=0; i<N; i++) {
		double x = i*2*M_PI/dims[0]; 
		for (int j=0; j<N; j++) {
			double y = j*2*M_PI/dims[1]; 
			for (int k=0; k<dims[2]; k++) {
				u(i,j,k) = exp(-pow(x-M_PI,2)/.01) * 
					exp(-pow(y-M_PI,2)/.01); 
			}
		}
	}
	u.forward(); 
	Scalar out = u; 
	writer.add(out, "temp"); 

	int Nt = 100; 
	double Tend = 5; 
	for (int t=1; t<Nt+1; t++) {
		double T = (double)t*Tend/Nt; 

		for (int i=0; i<dims[0]; i++) {
			double m = out.freq(i,0); 
			for (int j=0; j<dims[1]; j++) {
				double n = out.freq(j,1); 
				for (int k=0; k<dims[2]; k++) {
					out(i,j,k) = u(i,j,k) * exp(-(m*m + n*n)*T); 
				}
			}
		}
		writer.write(); 
	}
}
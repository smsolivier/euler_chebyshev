#include "Writer.H"

using namespace std; 

int main() {
	int N = 32; 
	array<int,DIM> dims = {N, N, N}; 

	Writer writer; 

	Scalar f(dims, true); 
	writer.add(f, "f"); 

	array<int,DIM> ind; 
	for (int i=0; i<N; i++) {
		double x = i*2*M_PI/N; 
		ind[0] = i; 
		for (int j=0; j<N; j++) {
			double y = i*2*M_PI/N; 
			ind[1] = j; 
			for (int k=0; k<N; k++) {
				ind[2] = k; 
				double z = cos(k*M_PI/(N-1)); 

				f[ind] = cos(x)*cos(y)*z*z; 
			}
		}
	}

	writer.write(); 

}
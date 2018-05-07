#include "DataObjects.H"
#include <iostream> 

using namespace std; 

#define TOL 1e-3

#define PRINT(name, wrong) \
	if (wrong) cout << name << ": WRONG!" << endl; \
	else cout << name << ": my man!" << endl; 

int main(int argc, char* argv[]) {
	int N = 32; 
	if (argc > 1) N = atoi(argv[1]); 
	array<int,DIM> dims = {N, N, N}; 

	Scalar f(dims, true); 
	array<int,DIM> ind; 
	for (ind[0]=0; ind[0]<dims[0]; ind[0]++) {
		double x = ind[0]*2*M_PI/dims[0]; 
		for (ind[1]=0; ind[1]<dims[1]; ind[1]++) {
			double y = ind[1]*2*M_PI/dims[1]; 
			for (ind[2]=0; ind[2]<dims[2]; ind[2]++) {
				double z = cos(ind[2]*M_PI/(dims[2]-1)); 

				f[ind] = sin(x)*sin(y)*z*z; 
			}
		}
	}
	f.forward(); 

	Vector grad = f.gradient(); 

	Vector test; 
	grad.inverse(test); 
	bool wrong = false; 
	for (ind[0]=0; ind[0]<dims[0]; ind[0]++) {
		double x = ind[0]*2*M_PI/dims[0]; 
		for (ind[1]=0; ind[1]<dims[1]; ind[1]++) {
			double y = ind[1]*2*M_PI/dims[1]; 
			for (ind[2]=0; ind[2]<dims[2]; ind[2]++) {
				double z = cos(ind[2]*M_PI/(dims[2]-1)); 
				if (abs(test[0][ind] - cos(x)*sin(y)*z*z) > TOL) wrong = true;  
				if (abs(test[1][ind] - sin(x)*cos(y)*z*z) > TOL) wrong = true; 
				if (abs(test[2][ind] - sin(x)*sin(y)*2*z) > TOL) wrong = true; 
			}
		}
	}
	PRINT("gradient", wrong); 

	// --- test laplacian --- 
	Scalar lap = f.laplacian(); 
	lap.inverse(); 
	wrong = false; 
	for (ind[0]=0; ind[0]<dims[0]; ind[0]++) {
		double x = ind[0]*2*M_PI/dims[0]; 
		for (ind[1]=0; ind[1]<dims[1]; ind[1]++) {
			double y = ind[1]*2*M_PI/dims[1]; 
			for (ind[2]=0; ind[2]<dims[2]; ind[2]++) {
				double z = cos(ind[2]*M_PI/(dims[2]-1)); 
				if (abs(lap[ind].real() - 
					(-2*sin(x)*sin(y)*z*z + sin(x)*sin(y)*2.)) > TOL) 
					wrong = true; 
			}
		}
	}
	PRINT("scalar laplacian", wrong); 
}
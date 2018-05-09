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

	// --- test cross product --- 
	Vector v1(dims, true); 
	Vector v2(dims, true); 
	for (int i=0; i<v1.size(); i++) {
		v1[2][i] = 1.; 
		v2[1][i] = 1.; 
	}
	v1.forward(); 
	v2.forward(); 
	Vector cross = v1.cross(v2); 
	cross.inverse(); 
	wrong = false; 
	for (int i=0; i<dims[0]; i++) {
		for (int j=0; j<dims[1]; j++) {
			for (int k=0; k<dims[2]; k++) {
				if (abs(cross[0](i,j,k).real() + 1) > TOL) wrong = true; 
				if (abs(cross[1](i,j,k).real()) > TOL) wrong = true; 
				if (abs(cross[2](i,j,k).real()) > TOL) wrong = true; 
			}
		}
	}
	PRINT("cross product", wrong); 

	// --- test curl --- 
	Vector curl = grad.curl(); 
	curl.inverse(); 
	wrong = false; 
	for (int i=0; i<curl.size(); i++) {
		for (int d=0; d<DIM; d++) {
			if (abs(curl[d][i].real()) > TOL) wrong = true; 
		}
	}
	PRINT("curl", wrong); 

	// --- divergence --- 
	Vector v(dims); 
	for (int d=0; d<DIM; d++) {
		v[d] = f; 
	}
	Scalar div = v.divergence(); 
	div.inverse(); 
	wrong = false; 
	for (int i=0; i<dims[0]; i++) {
		double x = i*2*M_PI/dims[0]; 
		for (int j=0; j<dims[1]; j++) {
			double y = j*2*M_PI/dims[1]; 
			for (int k=0; k<dims[2]; k++) {
				double z = cos(k*M_PI/(dims[2]-1)); 
				if (abs(div(i,j,k).real() - (
					cos(x)*sin(y)*z*z + sin(x)*cos(y)*z*z +
					sin(x)*sin(y)*2*z)) > TOL) wrong = true;  
			}
		}
	}
	PRINT("divergence", wrong); 

	// --- arithmetic --- 
	Vector a(dims, true); 
	Vector b(dims, true); 
	for (int i=0; i<a.size(); i++) {
		for (int d=0; d<DIM; d++) {
			a[d][i] = 1.; 
			b[d][i] = 1.; 
		}
	}

	Vector sum = a+b; 
	Vector sub = a - b; 
	Vector mult = 2.*a; 
	wrong = false; 
	for (int i=0; i<a.size(); i++) {
		for (int d=0; d<DIM; d++) {
			if (abs(sum[d][i] - 2.) > TOL) wrong = true; 
		}
	}
	PRINT("vector addition", wrong); 

	wrong = false; 
	for (int i=0; i<a.size(); i++) {
		for (int d=0; d<DIM; d++) {
			if (abs(sub[d][i]) > TOL) {
				wrong = true; 
				cout << sub[d][i] << endl; 
			}
		}
	}
	PRINT("vector subtraction", wrong); 

	wrong = false; 
	for (int i=0; i<a.size(); i++) {
		for (int d=0; d<DIM; d++) {
			if (abs(mult[d][i] - 2.) > TOL) {
				wrong = true; 
			}
		}
	}
	PRINT("double * Vector", wrong); 

	Scalar smult = 2.*a[0]; 
	wrong = false; 
	for (int i=0; i<smult.size(); i++) {
		if (abs(smult[i] - 2.) > TOL) wrong = true; 
	}
	PRINT("double * Scalar", wrong); 
}
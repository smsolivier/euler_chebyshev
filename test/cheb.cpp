#include "Cheb1D.H"
#include <vector> 
#include <complex> 
#include <iostream>

using namespace std; 

#define PRINT(wrong, message) if (wrong) cout << message << ": WRONG!" << endl; \
	else cout << message << ": pass" << endl; 

int main() {
	int N = 8; 

	vector<double> theta(N), z(N); 
	vector<complex<double>> f(N); 
	// vector<double> f(N); 
	Cheb1D cheb1d(N, 1, &f[0]); 

	for (int i=0; i<N; i++) {
		theta[i] = i*M_PI/(N-1); 
		z[i] = cos(theta[i]); 
		f[i] = 2*z[i]*z[i] - 1; 
		// f[i] = 1; 
	}

	cheb1d.transform(&f[0], -1); 
	bool wrong = false; 
	vector<double> ans(N); 
	ans[2] = 1; 
	for (int i=0; i<N; i++) {
		if (abs(ans[i] - f[i]) > 1e-5) wrong = true; 
	}
	PRINT(wrong, "1D forward chebyshev transform"); 

	cheb1d.transform(&f[0], 1); 
	wrong = false; 
	for (int i=0; i<N; i++) {
		if (abs(f[i] - (2*z[i]*z[i]-1)) > 1e-5) wrong = true; 
	}
	PRINT(wrong, "1D inverse chebyshev transform"); 

}
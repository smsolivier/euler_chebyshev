#include "Cheb1D.H"
#include <iostream> 

using namespace std; 

int main() {
	int N = 8; 
	vector<cdouble> f(N); 

	Cheb1D cheb(N, 1, &f[0]); 

	vector<double> theta(N), z(N); 

	for (int i=0; i<N; i++) {
		theta[i] = i*M_PI/(N-1); 
		z[i] = cos(theta[i]); 
		// f[i] = 2*z[i]*z[i] - 1; // T2 
		f[i] = z[i]*z[i]*z[i]; 
		// f[i] = z[i]*z[i] + z[i]*z[i]*z[i]; 
	}

	cheb.transform(&f[0], -1); 

	cheb.deriv(&f[0]); 

	cheb.transform(&f[0], 1); 

	for (int i=0; i<N; i++) {
		cout << f[i].real() << endl; 
	}
	cout << endl; 

	cheb.transform(&f[0], -1); 
	cheb.deriv(&f[0]); 
	cheb.transform(&f[0], 1); 

	for (int i=0; i<N; i++) {
		cout << f[i].real() << endl; 
	}
}
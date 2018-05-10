#include "Cheb1D.H"
#include <iostream> 

using namespace std; 

int main(int argc, char* argv[]) {
	int N = 8; 
	if (argc > 1) N = atoi(argv[1]); 
	vector<cdouble> f(N); 

	Cheb1D cheb(N, 1, &f[0]); 

	vector<cdouble> theta(N), z(N); 

	for (int i=0; i<N; i++) {
		theta[i] = i*M_PI/(N-1); 
		z[i] = cos(theta[i]); 
		f[i] = z[i]*z[i]*z[i]; 
	}

	vector<cdouble> d(N), d2(N); 

	cheb.transform(&f[0], -1); 
	cheb.deriv(&f[0], &d[0]); 
	cheb.transform(&d[0], 1); 

	bool wrong = false; 
	for (int i=0; i<N; i++) {
		if (abs(d[i] - 3.*z[i]*z[i]) > 1e-3) wrong = true; 
	}
	if (wrong) cout << "first derivative: WRONG!" << endl; 
	else cout << "first derivative: my man!" << endl; 

	cheb.transform(&d[0], -1); 
	cheb.deriv(&d[0], &d2[0]); 
	cheb.transform(&d2[0], 1); 

	wrong = false; 
	for (int i=0; i<N; i++) {
		if (abs(d2[i] - 6.*z[i]) > 1e-3) wrong = true; 
	}
	if (wrong) cout << "second derivative: WRONG!" << endl; 
	else cout << "second derivate: my man!" << endl; 
}
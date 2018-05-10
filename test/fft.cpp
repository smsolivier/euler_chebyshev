#include "FFT1D.H" 

using namespace std; 

int main() {
	int N = 32; 

	cdouble* f = new cdouble[N]; 
	double* ans = new double[N]; 

	FFT1D fft(N, 1, f); 

	for (int i=0; i<N; i++) {
		f[i] = (double)rand()/RAND_MAX; 
		ans[i] = f[i].real(); 
	}

	fft.transform(f, -1); 
	fft.transform(f, 1); 

	bool wrong = false; 
	for (int i=0; i<N; i++) {
		if (abs(f[i].real() - ans[i]) > 1e-3) wrong = true; 
	}
	if (wrong) cout << "1D FFT: WRONG!" << endl; 
	else cout << "1D FFT: pass" << endl; 

	delete f; 
	delete ans; 

}

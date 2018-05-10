#include "DataObjects.H"
#include <iostream> 

using namespace std;

int main(int argc, char* argv[]) {
#ifdef ZERO 
	ERROR("answers will be wrong with ZERO defined"); 
#endif
	int N = 8;  
	if (argc > 1) N = atoi(argv[1]); 
	array<int,DIM> dims = {N, N, N}; 

	Scalar s(dims, true); 
	Scalar ans(dims, true); 

	s.memory(); 

	for (int i=0; i<s.size(); i++) {
		s[i] = (double)rand()/RAND_MAX; 
		ans[i] = s[i].real(); 
	}

	s.forward(); 
	s.inverse(); 

	bool wrong = false; 
	for (int i=0; i<s.size(); i++) {
		if (abs(s[i].real() - ans[i].real()) > 1e-3) wrong = true; 
	}
	if (wrong) cout << "WRONG!" << endl; 
	else cout << "my man!" << endl; 

}
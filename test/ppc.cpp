#include "PaperCutter.H"
using namespace std; 

int main(int argc, char* argv[]) {
	int N = 32; 
	if (argc > 1) N = atoi(argv[1]); 
	Eigen::SparseMatrix<cdouble> theta(N, N); 
	double c; 
	for (int i=0; i<N; i++) {
		if (i==0) c = 2; 
		else c = 1; 
		for (int j=0; j<N; j++) {
			if ((i+j)%2==0 and j>i) {
				theta.insert(i,j) = 1/c*j*(j*j - i*i); 
			}
		}
	}
	for (int i=0; i<N; i++) {
		theta.coeffRef(i,i) -= 1.; 
		theta.coeffRef(N-2,i) = pow(i,2) * pow(-1,i+1); 
		theta.coeffRef(N-1,i) = pow(i,2);
		// theta.coeffRef(N-2,i) = pow(-1, i); 
		// theta.coeffRef(N-1,i) = 1;  
	}

	PaperCutter ppc; 
	ppc.init(&theta); 

	Eigen::VectorXcd rhs(N), sol(N); 
	for (int i=0; i<N; i++) {
		// rhs[i] = (double)rand()/RAND_MAX; 
		rhs[i] = i; 
	}

	ppc.solve(rhs, sol); 

	for (int i=0; i<1000; i++) {
		ppc.solve(rhs, sol); 
	}
	Eigen::VectorXcd mult = theta*sol; 
	bool wrong = false; 
	for (int i=0; i<N; i++) {
		if (abs(mult[i] - rhs[i]) > 1e-10) wrong = true; 
	}
	if (wrong) cout << "WRONG!" << endl; 
	else cout << "CORRECT! theta*x = rhs" << endl; 
}
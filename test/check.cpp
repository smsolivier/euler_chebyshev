#include "Scalar.H"

int main() {
	int N = 32; 
	array<int,DIM> dims = {N, N, N}; 
	Scalar phys(dims, true); 

	phys.inverse(); 
}
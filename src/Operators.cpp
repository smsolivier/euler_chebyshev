#include "DataObjects.H"

Vector operator+(const Vector& a, const Vector& b) {
	CH_TIMERS("Vector +"); 
	Vector ret(a.dims(), a.isPhysical()); 

	#pragma omp parallel for 
	for (int i=0; i<ret.size(); i++) {
		for (int d=0; d<DIM; d++) {
			ret[d][i] = a[d][i] + b[d][i]; 
		}
	}
	return ret; 
}

Vector operator-(const Vector& a, const Vector& b) {
	CH_TIMERS("Vector -"); 
	Vector ret(a.dims(), a.isPhysical()); 

	#pragma omp parallel for 
	for (int i=0; i<ret.size(); i++) {
		for (int d=0; d<DIM; d++) {
			ret[d][i] = a[d][i] - b[d][i]; 
		}
	}
	return ret; 
}

Vector operator*(const double alpha, const Vector& a) {
	CH_TIMERS("double * Vector"); 
	Vector ret(a.dims(), a.isPhysical()); 

	#pragma omp parallel for 
	for (int i=0; i<ret.size(); i++) {
		for (int d=0; d<DIM; d++) {
			ret[d][i] = alpha * a[d][i]; 
		}
	}
	return ret; 
}

Scalar operator*(const double alpha, const Scalar& a) {
	CH_TIMERS("double * Scalar"); 
	Scalar ret(a.dims(), a.isPhysical()); 

	#pragma omp parallel for 
	for (int i=0; i<ret.size(); i++) {
		ret[i] = alpha * a[i]; 
	}
	return ret; 
}
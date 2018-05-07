#ifndef __CHEB1D_H__
#define __CHEB1D_H__

#include <vector> 
#include <complex> 
#include <array> 
#ifdef FFTW 
#include "fftw3.h"
#endif

#include "Common.H"

using namespace std; 

/// Chebyshev transform in 1 dimensions 
class Cheb1D {
public:
	/// default constructor 
	Cheb1D(); 
	/// constructor 
	Cheb1D(int N, int stride, cdouble* input); 
	/// set stride and size to be able to reuse the plan 
	void init(int N, int stride, cdouble* input); 
	/// destructor: destroys plans 
	~Cheb1D(); 
	/// use the premade plan to transform 
	/** dir = -1 is forward \n
		dir = 1 is inverse 
	**/ 
	void transform(cdouble* input, int DIR); 
private:
	/// FFTW DCT plan 
	fftw_plan m_plan; 
	/// size of DCT 
	int m_N; 
	/// stride 
	int m_stride; 
}; 

#endif
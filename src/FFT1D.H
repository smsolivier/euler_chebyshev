#ifndef __FFT_1D_H__
#define __FFT_1D_H__

#include <vector> 
#include <complex> 
#include <array> 
#ifdef FFTW 
#include "fftw3.h" 
#endif

#include "Common.H"

using namespace std; 

/// do FFT in 1 dimension 
class FFT1D {
public:
	/// default constructor 
	FFT1D(); 
	/// constructor 
	FFT1D(int N, int stride, cdouble* input); 
	/// set stride and size to reuse plan 
	void init(int N, int stride, cdouble* input); 
	/// destructor 
	~FFT1D(); 
	/// use premade plans for transform 
	/** dir = -1 is forward transform \n
		dir = 1 is inverse transform 
	**/ 
	void transform(cdouble* input, int DIR); 
private: 
	/// fftw plan for many FFTs of the same size and stride 
	fftw_plan m_forward;
	/// fftw plan for inverse 
	fftw_plan m_backward;  
	/// size of FFT 
	int m_N; 
	/// stride 
	int m_stride; 
	/// FFT direction 
	int m_dir; 
	/// save alignment of original array (subsequent calls must have same alignment) 
	int m_alignment; 
}; 

#endif

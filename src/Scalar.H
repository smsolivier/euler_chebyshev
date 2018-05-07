#ifndef __SCALAR_H__
#define __SCALAR_H__
#include "FFT1D.H"
#include "Cheb1D.H"

/// represent a scalar variable 
class Scalar {
public: 
	/// default constructor 
	Scalar(); 
	/// set size in DIM dimensions. Defaults to creation in FFC space 
	Scalar(array<int,DIM> N, bool physical=false); 
	/// destructor 
	~Scalar(); 
	/// copy constructor. initializes and deeps copies values over 
	Scalar(const Scalar& scalar); 
	/// copy assignment. Deep copy values 
	void operator=(const Scalar& scalar); 

	/// initialize a Scalar 
	void init(array<int,DIM> N, bool physical=false); 
	/// index into data with 3D indexing 
	cdouble& operator[](array<int,DIM> ind); 
	/// const access from 3D indexing 
	const cdouble& operator[](array<int,DIM> ind) const; 
	/// explicit 3D indexing 
	cdouble& operator()(int a_i, int a_j, int a_k); 
	/// const explicit 3D indexing 
	const cdouble& operator()(int a_i, int a_j, int a_k) const; 
	/// index into data with flattened indexing 
	cdouble& operator[](int index); 
	/// const access with flattened indexing 
	const cdouble& operator[](int index) const; 
	/// get dimensions 
	array<int,DIM> dims() const; 
	/// get total size 
	int size() const; 

	/// in place transform to FFC space 
	void forward(); 
	/// in place transform to physical space 
	void inverse(); 
	/// out of place tnransform to physical space 
	void inverse(Scalar& scalar) const; 
	/// get fourier frequency number 
	double freq(int ind, int d) const; 

	/// gradient: \f$ \nabla f \f$ 
	Vector gradient() const; 
	/// laplacian: \f$ \nabla^2 f \f$ 
	Scalar laplacian() const; 
	/// invert \f$ (a - b\nabla^2) f = f \f$ 
	void invert_laplacian(double a=0, double b=1); 

	/// output total memory of initialized scalar 
	void memory() const; 

	/// query if (*this) is in physical space 
	bool isPhysical() const; 
	/// query if (*this) is in FFC space 
	bool isFFC() const; 
	/// set to physical space 
	void setPhysical(); 
	/// set to FFC space 
	void setFFC(); 

	/// return derivative in z direction 
	void ddz(int a_i, int a_j, cdouble* output) const; 
private:
	/// remove highest fourier modes 
	void zeroHighModes(); 
	/// 3D to 1D indexing 
	inline int index(array<int,DIM> ind) const; 

	/// size in DIM dimensions 
	array<int,DIM> m_N; 
	/// total size 
	int m_size; 
	/// store complex double data 
	vector<cdouble> m_data; 	
	/// fft in x direction 
	FFT1D m_fft_x; 
	/// fft in y direction 
	FFT1D m_fft_y; 
	/// cheb in z direction 
	Cheb1D m_cheb; 
	/// true if in FFC space, false if physical 
	bool m_ffc; 

	/// store total number of scalars allocated 
	static int m_nscalars;  
}; 

#endif
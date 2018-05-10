#ifndef __INVERTER_H__
#define __INVERTER_H__
#include <vector> 
#include "Common.H" 
#include "DataObjects.H"
#include "PaperCutter.H"

using namespace std; 

/// invert \f$ \nabla^2 f = \left( D^2 - \vec{k}^2 \right) \f$ 
class Inverter {
public: 
	/// singleton interface
	/** allow only one static instance to exist to prevent redoing expensive factorizations **/ 
	static Inverter& instance() {
		static Inverter inv; 
		return inv; 
	}
	/// destructor 
	~Inverter(); 
	/// initialize inverses 
	/** BC=0: double neumann boundaries \\BC=1: double dirichlet \\BC=2: mixed **/
	void init(array<int,DIM> N, int BC=0); 
	/// invert a scalar 
	void invert(Scalar& scalar, Vector& V); 
	/// invert a scalar 
	void invert(Scalar& scalar); 
private:
	/// default constructor 
	Inverter(); 
	/// constructor 
	Inverter(array<int,DIM> N);
	/// fourier frequencies 
	double freq(int ind, int d) const; 

	/// size in DIM dimensions 
	array<int,DIM> m_N; 
	/// second derivative matrix (in sparse Eigen format) 
	Eigen::SparseMatrix<cdouble> m_D2; 
	/// store all matrices 
	Eigen::SparseMatrix<cdouble>** m_theta; 
	/// pointers to PaperCutter's
	PaperCutter** m_ppc; 
	/// store if init has been called 
	bool m_init; 
}; 

#endif
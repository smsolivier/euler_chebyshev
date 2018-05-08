#ifndef __INVERTER_H__
#define __INVERTER_H__
#include <vector> 
#include "Common.H" 
#include "DataObjects.H"
#include "Sparse" 

using namespace std; 

/// invert \f$ \nabla^2 f = \left( D^2 - \vec{k}^2 \right) \f$ 
class Inverter {
public: 
	static Inverter& instance() {
		static Inverter inv; 
		return inv; 
	}
	/// destructor 
	~Inverter(); 
	/// initialize inverses 
	void init(array<int,DIM> N); 
	/// invert a scalar 
	void invert(Scalar& scalar); 
private:
	Inverter(); 
	/// constructor 
	Inverter(array<int,DIM> N); 
	/// size in DIM dimensions 
	array<int,DIM> m_N; 
	/// second derivative matrix (in sparse Eigen format) 
	Eigen::SparseMatrix<cdouble> m_D2; 
	/// inverse matrices 
	Eigen::SparseLU<Eigen::SparseMatrix<cdouble>>** m_solvers; 
	/// fourier frequencies 
	double freq(int ind, int d) const; 
	/// store if init has been called 
	bool m_init; 
}; 

#endif
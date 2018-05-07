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
	/// constructor 
	Inverter(array<int,DIM> N); 
	/// destructor 
	~Inverter(); 
	/// invert a scalar 
	void operator()(Scalar& scalar); 
private:
	/// size in DIM dimensions 
	array<int,DIM> m_N; 
	/// second derivative matrix (in sparse Eigen format) 
	Eigen::SparseMatrix<cdouble> m_D2; 
	/// inverse matrices 
	Eigen::SparseLU<Eigen::SparseMatrix<cdouble>>** m_solvers; 
	/// fourier frequencies 
	double freq(int ind, int d) const; 
}; 

#endif
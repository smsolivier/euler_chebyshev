#ifndef __GREENS_H__
#define	__GREENS_H__
#include "Common.H" 
#include "DataObjects.H" 

using namespace std; 

// number of green's functions 
#define NG 2

/// Represent the two Green's function when inverting the pressure 
class Greens {
public:
	/// constructor 
	Greens(array<int,DIM> N); 
	/// const access to scalars (green's functions should not need to be altered) 
	const Scalar& operator[](int a_i) const; 

	/// solve 2x2 system for \f$ \tau \f$ 
	void tau(Vector& V34, Vector& V); 
private:
	/// store NG scalars 
	array<Scalar,NG> m_scalars; 
	/// store NG laplacian's of m_scalars 
	array<Scalar,NG> m_lap; 
	/// store NG gradients of m_scalars 
	array<Vector,NG> m_grad; 
	/// size in DIM dimensions 
	array<int,DIM> m_N; 
	/// store M, M-1 coefficients 
	array<cdouble,NG*NG> m_coef; 
	/// store 2x2 components for all \f$k_x\,, k_y \f$ 
	vector<vector<cdouble>> m_a; 
	/// store 2x2 components for all \f$k_x\,, k_y \f$ 
	vector<vector<cdouble>> m_b; 
	/// store 2x2 components for all \f$k_x\,, k_y \f$ 
	vector<vector<cdouble>> m_c; 
	/// store 2x2 components for all \f$k_x\,, k_y \f$ 
	vector<vector<cdouble>> m_d; 
	/// store 2x2 determinant for all \f$k_x\,, k_y \f$ 
	vector<vector<cdouble>> m_det; 
	/// store \f$ T_M \hat{z} - \nabla G_1 \f$ 
	Vector m_v1; 
	/// store \f$ T_{M-1} \hat{z} - \nabla G_{-1} \f$ 
	Vector m_v2; 
}; 

#endif
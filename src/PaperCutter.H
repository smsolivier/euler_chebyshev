#ifndef __PAPERCUTTER_H__
#define __PAPERCUTTER_H__
#include "Sparse"
#include "Common.H"

/// solves a difficult matrix with the "paper cutter" method 
class PaperCutter {
public:
	/// initialize reused elements of paper cutter method 
	PaperCutter(Eigen::SparseMatrix<cdouble>& theta); 
	/// solve the system 
	void solve(Eigen::VectorXcd& rhs, Eigen::VectorXcd& sol); 
private:
	/// original matrix 
	Eigen::SparseMatrix<cdouble> m_theta; 
	/// SparseLU solver 
	Eigen::SparseLU<Eigen::SparseMatrix<cdouble>> m_solver; 
	/// reduced matrix 
	Eigen::SparseMatrix<cdouble> m_M; 
	/// left column \f$ M^{-1} L \f$ 
	Eigen::VectorXcd m_L; 
	/// right column \f$ M^{-1} R \f$ 
	Eigen::VectorXcd m_R; 
	/// upper removed row 
	Eigen::RowVectorXcd m_T; 
	/// lower removed row 
	Eigen::RowVectorXcd m_B; 
	/// values left over from paper cutter intersection 
	array<cdouble,4> m_vals;
	/// store coefficients of inverse 2x2 matrix  
	array<cdouble,4> m_inv; 
}; 

#endif
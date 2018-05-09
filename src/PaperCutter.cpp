#include "PaperCutter.H"

PaperCutter::PaperCutter() {}

void PaperCutter::init(Eigen::SparseMatrix<cdouble>& theta) {
	m_theta = theta; 

	m_M = m_theta.block(0,2, m_theta.rows()-2, m_theta.cols()-2); 
	m_solver.compute(m_M); 

	Eigen::VectorXcd L = m_theta.block(0,0, m_theta.rows()-2, 1); 
	m_L = m_solver.solve(L); 

	Eigen::VectorXcd R = m_theta.block(0,1, m_theta.rows()-2, 1); 
	m_R = m_solver.solve(R); 

	m_T = m_theta.block(m_theta.rows()-2,2, 1, m_theta.cols()-2); 
	m_B = m_theta.block(m_theta.rows()-1,2, 1, m_theta.cols()-2); 

	m_vals[0] = m_theta.coeff(m_theta.rows()-2,0); 
	m_vals[1] = m_theta.coeff(m_theta.rows()-2,1); 
	m_vals[2] = m_theta.coeff(m_theta.rows()-1,0); 
	m_vals[3] = m_theta.coeff(m_theta.rows()-1,1); 

	cdouble a = m_vals[0] - m_T.dot(m_L); 
	cdouble b = m_vals[1] - m_T.dot(m_R); 
	cdouble c = m_vals[2] - m_B.dot(m_L); 
	cdouble d = m_vals[3] - m_B.dot(m_R); 

	cdouble det = a*d - b*c; 
	m_skip = false; 
	if (det == 0. || m_solver.info() != 0) {
		cout << "skipping" << endl; 
		m_skip = true; 
	}

	// compute inverse 
	m_inv[0] = d/det; 
	m_inv[1] = -b/det; 
	m_inv[2] = -c/det; 
	m_inv[3] = a/det; 
}

void PaperCutter::solve(Eigen::VectorXcd& rhs, Eigen::VectorXcd& sol) {
	if (m_skip) return; 
	Eigen::VectorXcd g(rhs.rows()-2); 
	for (int i=0; i<rhs.rows()-2; i++) {
		g[i] = rhs[i]; 
	}

	Eigen::VectorXcd gt = m_solver.solve(g); 

	cdouble g0 = rhs[rhs.rows()-2] - m_T.dot(gt); 
	cdouble g1 = rhs[rhs.rows()-1] - m_B.dot(gt); 

	cdouble f0 = m_inv[0]*g0 + m_inv[1]*g1; 
	cdouble f1 = m_inv[2]*g0 + m_inv[3]*g1; 

	sol[0] = f0; 
	sol[1] = f1; 

	for (int i=0; i<rhs.rows()-2; i++) {
		sol[i+2] = gt[i] - f0*m_L[i] - f1*m_R[i]; 
	}

	#ifdef ErrorCheck 
	Eigen::VectorXcd check = m_theta*sol - rhs; 
	double max = -1; 
	for (int i=0; i<m_theta.rows(); i++) {
		if (abs(check[i]) > max) max = abs(check[i]); 
	}
	CHECK(max < 1e-6, "solver not exact" + to_string(m_skip)); 
	#endif
}
#include "Inverter.H"

#define DIRICHLET
// #define NEUMANN

Inverter::Inverter() {
	m_init = false; 
}

void Inverter::init(array<int,DIM> N, int BC) {
	if (m_init) return; 
	CH_TIMERS("initialize inverter"); 

	m_init = true; 
	m_N = N; 

	m_D2.resize(m_N[2], m_N[2]); 

	double c; 
	for (int i=0; i<m_N[2]; i++) {
		if (i==0) c = 2; 
		else c = 1; 
		for (int j=0; j<m_N[2]; j++) {
			if ((i+j)%2==0 and j>i) {
				m_D2.insert(i,j) = 1/c*j*(j*j - i*i); 
			}
		}
	}

	m_solvers = new Eigen::SparseLU<Eigen::SparseMatrix<cdouble>>*[m_N[0]*m_N[1]]; 
	m_theta = new Eigen::SparseMatrix<cdouble>*[m_N[0]*m_N[1]]; 

	Eigen::SparseMatrix<cdouble> theta(m_N[2], m_N[2]); 
	for (int i=0; i<m_N[0]; i++) {
		double m2 = pow(freq(i,0), 2); 
		for (int j=0; j<m_N[1]; j++) {
			int index = j+i*m_N[1]; 
			double n2 = pow(freq(j,1), 2); 
			m_theta[index] = new Eigen::SparseMatrix<cdouble>; 
			m_solvers[index] = new Eigen::SparseLU<Eigen::SparseMatrix<cdouble>>; 
			*m_theta[index] = m_D2; 
			for (int k=0; k<m_N[2]; k++) {
				m_theta[index]->coeffRef(k,k) -= (m2 + n2); 
				if (BC == 0) { // double neumann 
					m_theta[index]->coeffRef(m_N[2]-1,k) = pow(k,2); 
					m_theta[index]->coeffRef(m_N[2]-2,k) = pow(k,2)*pow(-1, k+1); 					
				} else if (BC == 1) { // double dirichlet 
					m_theta[index]->coeffRef(m_N[2]-2,k) = pow(-1,k); 
					m_theta[index]->coeffRef(m_N[2]-1,k) = 1; 					
				} else { // mixed 
					m_theta[index]->coeffRef(m_N[2]-1,k) = 1; 
					m_theta[index]->coeffRef(m_N[2]-2,k) = pow(k,2)*pow(-1, k+1);
					// m_theta[index]->coeffRef(m_N[2]-1,k) = pow(k,2);
					// m_theta[index]->coeffRef(m_N[2]-2,k) = pow(-1,k); 
				}
			}

			m_solvers[index]->compute(*m_theta[index]);
			// m_solvers[index]->analyzePattern(*m_theta[index]); 
			// m_solvers[index]->factorize(*m_theta[index]); 
		}
	}
}

Inverter::~Inverter() {
	if (m_init) {
		for (int i=0; i<m_N[0]*m_N[1]; i++) {
			delete m_solvers[i]; 
			delete m_theta[i]; 
		}
		delete m_solvers; 
		delete m_theta; 
	}

	m_init = false; 
}

void Inverter::invert(Scalar& scalar, Vector& V) {
	CH_TIMERS("invert theta"); 
	if (!m_init) init(scalar.dims()); 
	Scalar orig = scalar; // copy

	Eigen::VectorXcd sol(m_N[2]); 
	Eigen::VectorXcd rhs(m_N[2]); 
	Eigen::VectorXcd check(m_N[2]); 

	for (int j=0; j<m_N[1]; j++) {
		double n = freq(j,1); 
		for (int i=0; i<m_N[0]; i++) {
			double m = freq(i,0); 
			int index = j + i*m_N[1]; 
			for (int k=0; k<m_N[2]; k++) {
				rhs[k] = orig(i,j,k); 
			}
			rhs[m_N[2]-1] = 0; 
			rhs[m_N[2]-2] = 0; 
			for (int k=0; k<m_N[2]; k++) {
				rhs[m_N[2]-2] += V[2](i,j,k) * pow(-1., k); 
				rhs[m_N[2]-1] += V[2](i,j,k); 
			}
			if (m_solvers[index]->info() != 0) continue; 
			CHECK(m_solvers[index]->info() == 0, "factorization issue"); 
			sol = m_solvers[index]->solve(rhs); 
			for (int k=0; k<m_N[2]; k++) {
				scalar(i,j,k) = sol[k]; 
			}
			#ifdef ErrorCheck
			check = *m_theta[index]*sol - rhs; 
			double max = -1; 
			for (int k=0; k<m_N[2]; k++) {
				if (abs(check[k]) > max) max = abs(check[k]); 
			}
			CHECK(max < 1e-15, "solver not exact " + to_string(max)); 
			#endif
		}
	}
}

void Inverter::invert(Scalar& scalar) {
	CH_TIMERS("invert theta"); 
	if (!m_init) init(scalar.dims()); 
	Scalar orig = scalar; // copy

	Eigen::VectorXcd sol(m_N[2]); 
	Eigen::VectorXcd rhs(m_N[2]); 
	Eigen::VectorXcd check(m_N[2]); 

	for (int j=0; j<m_N[1]; j++) {
		double n = freq(j,1); 
		for (int i=0; i<m_N[0]; i++) {
			double m = freq(i,0); 
			int index = j + i*m_N[1]; 
			for (int k=0; k<m_N[2]; k++) {
				rhs[k] = orig(i,j,k); 
			}
			if (m_solvers[index]->info() != 0) continue; 
			CHECK((bool)(m_solvers[index]->info() == 0), "factorization issue"); 
			sol = m_solvers[index]->solve(rhs); 
			for (int k=0; k<m_N[2]; k++) {
				scalar(i,j,k) = sol[k]; 
			}
			#ifdef ErrorCheck
			check = *m_theta[index]*sol - rhs; 
			double max = -1; 
			for (int k=0; k<m_N[2]; k++) {
				if (abs(check[k]) > max) max = abs(check[k]); 
			}
			CHECK((bool)(max < 1e-6), "solver not exact " + to_string(max)); 
			#endif
		}
	}

	// Eigen::SparseLU<Eigen::SparseMatrix<cdouble>> solver; 
	// for (int i=0; i<m_N[0]; i++) {
	// 	double m = freq(i,0); 
	// 	for (int j=0; j<m_N[1]; j++) {
	// 		double n = freq(j,1); 
	// 		// if (n+m != 0) {
	// 		if (true) {
	// 			Eigen::SparseMatrix<cdouble> theta = m_D2; 
	// 			for (int k=0; k<m_N[2]; k++) {
	// 				theta.coeffRef(k,k) -= (m*m + n*n); 
	// 				theta.coeffRef(m_N[2]-2,k) = 1.; 
	// 				theta.coeffRef(m_N[2]-1,k) = pow(-1., k); 
	// 				rhs[k] = orig(i,j,k); 
	// 			}

	// 			rhs[m_N[2]-1] = 0; 
	// 			rhs[m_N[2]-2] = 0; 

	// 			solver.compute(theta); 
	// 			sol = solver.solve(rhs); 
	// 			for (int k=0; k<m_N[2]; k++) {
	// 				scalar(i,j,k) = sol[k]; 
	// 			}

	// 		}
	// 	}
	// }
}

double Inverter::freq(int ind, int d) const {
#ifndef NCHECK 
	if (d >= 2) ERROR("frequency not defined for z dimension"); 
#endif
	if (ind <= m_N[d]/2) return (double)ind; 
	else return -1.*(double)m_N[d] + (double)ind;  
}
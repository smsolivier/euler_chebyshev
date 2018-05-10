#include "Inverter.H"

Inverter::Inverter() {
	m_init = false; 
}

void Inverter::init(array<int,DIM> N, int BC) {
	if (m_init) return; 
	CH_TIMERS("initialize inverter"); 

	m_init = true; 
	m_N = N; 

	// compute second derivative matrix 
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

	// store solver object and matrices 
	m_ppc = new PaperCutter*[m_N[0]*m_N[1]]; 
	m_theta = new Eigen::SparseMatrix<cdouble>*[m_N[0]*m_N[1]]; 

	// run LU factorization on all matrices 
	#pragma omp parallel for 
	for (int i=0; i<m_N[0]; i++) {
		double m2 = pow(freq(i,0), 2); 
		for (int j=0; j<m_N[1]; j++) {
			int index = j+i*m_N[1]; 
			double n2 = pow(freq(j,1), 2); 

			// compute theta 
			m_theta[index] = new Eigen::SparseMatrix<cdouble>(m_N[2], m_N[2]); 
			*m_theta[index] = m_D2; // set to second derivative 
			for (int k=0; k<m_N[2]; k++) {
				m_theta[index]->coeffRef(k,k) -= (m2 + n2); // set to theta 

				// set boundary conditions 
				if (BC == 0) { // double neumann 
					m_theta[index]->coeffRef(m_N[2]-1,k) = pow(k,2); 
					m_theta[index]->coeffRef(m_N[2]-2,k) = pow(k,2)*pow(-1, k+1); 					
				} else if (BC == 1) { // double dirichlet 
					m_theta[index]->coeffRef(m_N[2]-2,k) = pow(-1,k); 
					m_theta[index]->coeffRef(m_N[2]-1,k) = 1; 					
				} else { // mixed 
					m_theta[index]->coeffRef(m_N[2]-1,k) = 1; 
					m_theta[index]->coeffRef(m_N[2]-2,k) = pow(k,2)*pow(-1, k+1);
					// theta.coeffRef(m_N[2]-1,k) = pow(k,2);
					// theta.coeffRef(m_N[2]-2,k) = pow(-1,k); 
				}
			}

			// initialize solver object 
			m_ppc[index] = new PaperCutter; 
			m_ppc[index]->init(m_theta[index]); 
		}
	}
}

Inverter::~Inverter() {
	if (m_init) {
		for (int i=0; i<m_N[1]*m_N[0]; i++) {
			delete m_ppc[i]; 
			delete m_theta[i]; 
		}
		delete m_ppc; 
		delete m_theta; 
	}

	m_init = false; 
}

void Inverter::invert(Scalar& scalar, Vector& V) {
	CH_TIMERS("invert theta"); 
	if (!m_init) init(scalar.dims()); 
	Scalar orig = scalar; // copy

	#pragma omp parallel 
	{
		Eigen::VectorXcd sol(m_N[2]); 
		Eigen::VectorXcd rhs(m_N[2]);
		#pragma omp for 
		for (int i=0; i<m_N[0]; i++) {
			for (int j=0; j<m_N[1]; j++) {
				int index = j + i*m_N[1]; // flatten 

				// copy to Eigen RHS vector 
				for (int k=0; k<m_N[2]; k++) {
					rhs[k] = orig(i,j,k); 
				}

				// compute V_z(z=+-1) (in physical space) 
				rhs[m_N[2]-1] = 0; 
				rhs[m_N[2]-2] = 0; 
				for (int k=0; k<m_N[2]; k++) {
					rhs[m_N[2]-2] += V[2](i,j,k) * pow(-1., k); 
					rhs[m_N[2]-1] += V[2](i,j,k); 
				}

				// call matrix solver 
				m_ppc[index]->solve(rhs, sol); 

				// copy from Eigen to Scalar 
				for (int k=0; k<m_N[2]; k++) {
					scalar(i,j,k) = sol[k]; 
				}
			}
		}
	}
}

void Inverter::invert(Scalar& scalar) {
	CH_TIMERS("invert theta"); 
	if (!m_init) init(scalar.dims()); 
	Scalar orig = scalar; // copy

	#pragma omp parallel 
	{
		Eigen::VectorXcd sol(m_N[2]); 
		Eigen::VectorXcd rhs(m_N[2]);
		#pragma omp for 
		for (int i=0; i<m_N[0]; i++) {
			for (int j=0; j<m_N[1]; j++) {
				int index = j + i*m_N[1]; 

				// copy to Eigen RHS vector 
				for (int k=0; k<m_N[2]; k++) {
					rhs[k] = orig(i,j,k); 
				}

				// call solver object 
				m_ppc[index]->solve(rhs, sol); 
				
				// copy back to Scalar 
				for (int k=0; k<m_N[2]; k++) {
					scalar(i,j,k) = sol[k]; 
				}
			}
		}		
	}
}

double Inverter::freq(int ind, int d) const {
#ifndef NCHECK 
	if (d >= 2) ERROR("frequency not defined for z dimension"); 
#endif
	if (ind <= m_N[d]/2) return (double)ind; 
	else return -1.*(double)m_N[d] + (double)ind;  
}
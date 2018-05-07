#include "Inverter.H"

Inverter::Inverter(array<int,DIM> N) {
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

	Eigen::SparseMatrix<cdouble> theta(m_N[2], m_N[2]); 
	for (int i=0; i<m_N[0]; i++) {
		double m2 = pow(freq(i,0), 2); 
		for (int j=0; j<m_N[1]; j++) {
			int index = j+i*m_N[1]; 
			double n2 = pow(freq(j,1), 2); 
			theta = m_D2; 
			for (int k=0; k<m_N[2]; k++) {
				theta.coeffRef(k,k) -= (m2 + n2); 
				theta.coeffRef(m_N[2]-2,k) = 1; 
				theta.coeffRef(m_N[2]-1,k) = pow(-1., k); 
			}

			m_solvers[index] = new Eigen::SparseLU<Eigen::SparseMatrix<cdouble>>; 
			m_solvers[index]->compute(theta); 
		}
	}
}

Inverter::~Inverter() {
	// delete m_solvers; 
}

void Inverter::operator()(Scalar& scalar) {
	Scalar orig(scalar); // copy 

	Eigen::VectorXcd sol(m_N[2]); 
	Eigen::VectorXcd rhs(m_N[2]); 

	Eigen::SparseMatrix<cdouble> theta(m_N[2], m_N[2]); 
	for (int j=0; j<m_N[1]; j++) {
		for (int i=0; i<m_N[0]; i++) {
			for (int k=0; k<m_N[2]; k++) {
				rhs[k] = orig(i,j,k); 
			}
			if (i==0 && j==0) {
				rhs[m_N[2]-1] = 1; 
				rhs[m_N[2]-2] = 1; 
			} else {
				rhs[m_N[2]-1] = 0; 
				rhs[m_N[2]-2] = 0; 
			}
			sol = m_solvers[j+i*m_N[1]]->solve(rhs); 
			for (int k=0; k<m_N[2]; k++) {
				scalar(i,j,k) = sol[k]; 
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
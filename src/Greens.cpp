#include "Greens.H" 
#include <vector> 
#include <cstring> 
#include <iostream> 
#include "Inverter.H"

#define ON 

Greens::Greens(array<int,DIM> N) {
#ifndef ON 
	WARNING("tau's turned off"); 
#endif
	m_N = N; 

	// initialize scalars in FFC space 
	for (int i=0; i<NG; i++) {
		m_scalars[i].init(m_N, false); 
	}

	m_v1.init(m_N); 
	m_v2.init(m_N); 

	// setup chebyshev transform object (for taking derivatives) 
	Cheb1D cheb(m_N[2], 1, &m_scalars[0][0]); 

	// store derivatives of T_M and and T_M-1
	vector<cdouble> Tm(m_N[2], 0), dTm(m_N[2], 0); 
	vector<cdouble> Tm1(m_N[2], 0), dTm1(m_N[2], 0); 

	// set to T_M 
	Tm[m_N[2]-1] = 1; 
	// set to T_M-1 
	Tm1[m_N[2]-2] = 1; 

	// take derivatives 
	cheb.deriv(&Tm[0], &dTm[0]); 
	cheb.deriv(&Tm1[0], &dTm1[0]); 

	// save dTm_M 
	m_coef[0] = dTm[m_N[2]-1]; 
	m_coef[1] = dTm1[m_N[2]-1]; 
	m_coef[2] = dTm[m_N[2]-2]; 
	m_coef[3] = dTm1[m_N[2]-2]; 

	// set boundary conditions 
	dTm[m_N[2]-1] = 1; 
	dTm[m_N[2]-2] = pow(-1, m_N[2]); 
	dTm1[m_N[2]-1] = 1; 
	dTm1[m_N[2]-2] = pow(-1, m_N[2]-1); 

	// copy data over to scalar to be able to use Inverter 
	for (int j=0; j<m_N[1]; j++) {
		for (int i=0; i<m_N[0]; i++) {
			memcpy(&m_scalars[0](i,j,0), &dTm[0], sizeof(cdouble)*m_N[2]); 
			memcpy(&m_scalars[1](i,j,0), &dTm1[0], sizeof(cdouble)*m_N[2]); 
		}
	}

	// call inverter 
	for (int i=0; i<NG; i++) {
		Inverter::instance().invert(m_scalars[i]); 
		m_lap[i] = m_scalars[i].laplacian(); 
		m_grad[i] = m_scalars[i].gradient(); 
	}

	// setup 2x2 
	m_a.resize(m_N[1]); 
	m_b.resize(m_N[1]); 
	m_c.resize(m_N[1]); 
	m_d.resize(m_N[1]); 
	m_det.resize(m_N[1]); 
	for (int j=0; j<m_N[1]; j++) {
		m_a[j].resize(m_N[0]); 
		m_b[j].resize(m_N[0]); 
		m_c[j].resize(m_N[0]); 
		m_d[j].resize(m_N[0]); 
		m_det[j].resize(m_N[0]); 
	}

	for (int j=0; j<m_N[1]; j++) {
		for (int i=0; i<m_N[0]; i++) {
			m_a[j][i] = m_lap[0](i,j,m_N[2]-1) - m_coef[0]; 
			m_b[j][i] = m_lap[1](i,j,m_N[2]-1) - m_coef[1]; 
			m_c[j][i] = m_lap[0](i,j,m_N[2]-2) - m_coef[2]; 
			m_d[j][i] = m_lap[1](i,j,m_N[2]-2) - m_coef[3]; 
			m_det[j][i] = m_a[j][i]*m_d[j][i] - m_b[j][i]*m_c[j][i]; 
		}
	}

	// compute T_m z - del G1 
	m_v1 = -1.*m_grad[0]; 
	for (int i=0; i<m_N[0]; i++) {
		for (int j=0; j<m_N[1]; j++) {
			m_v1[2](i,j,m_N[2]-1) += 1.; 
		}
	}

	// compute T_{M-1} z - del G_{-1} 
	m_v2 = -1.*m_grad[1]; 
	for (int i=0; i<m_N[0]; i++) {
		for (int j=0; j<m_N[1]; j++) {
			m_v2[2](i,j,m_N[2]-2) += 1.; 
		}
	}
}

const Scalar& Greens::operator[](int a_i) const {return m_scalars[a_i]; }

void Greens::tau(const Vector& V34, Vector& V) {
	Scalar div = V34.divergence(); 
	vector<vector<cdouble>> tau0, tau1; 
	if (tau0.size() == 0) {
		tau0.resize(m_N[1]); 
		tau1.resize(m_N[1]); 
		for (int j=0; j<m_N[1]; j++) {
			tau0[j].resize(m_N[0]); 
			tau1[j].resize(m_N[0]); 
		}
	}

#ifdef ON 
	for (int j=1; j<m_N[1]; j++) {
		for (int i=1; i<m_N[0]; i++) {
			CHECK((bool)(m_det[j][i] != 0.), 
				"determinant is zero " + to_string(i) + ", " + to_string(j)); 
			tau0[j][i] = 1./m_det[j][i]*(m_d[j][i]*div(i,j,m_N[2]-1) - 
				m_b[j][i]*div(i,j,m_N[2]-2)); 
			tau1[j][i] = 1./m_det[j][i]*(
				-m_c[j][i]*div(i,j,m_N[2]-1) + m_a[j][i]*div(i,j,m_N[2]-2)); 
			if (abs(tau0[j][i]) > 1e-5) cout << "tau0 = " << tau0[j][i] << endl; 
		}
	}
#endif

	for (int i=0; i<m_N[0]; i++) {
		for (int j=0; j<m_N[1]; j++) {
			for (int k=0; k<m_N[2]; k++) {
				for (int d=0; d<DIM; d++) {
					V[d](i,j,k) = V34[d](i,j,k) + tau0[j][i]*m_v1[d](i,j,k) + 
						tau1[j][i]*m_v2[d](i,j,k); 
				}
			}
		}
	}
}
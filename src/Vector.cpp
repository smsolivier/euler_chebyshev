#include "DataObjects.H" 

Vector::Vector() {} 

Vector::Vector(array<int,DIM> N, bool physical) {
	init(N, physical); 
}

void Vector::init(array<int,DIM> N, bool physical) {
	m_N = N; 
	for (int d=0; d<DIM; d++) {
		m_vector[d].init(m_N, physical); 
	}

	if (physical == true) setPhysical(); 
	else setFFC(); 
}

Scalar& Vector::operator[](int ind) {return m_vector[ind]; }
const Scalar& Vector::operator[](int ind) const {return m_vector[ind]; }

void Vector::forward() {
	CHECK(isPhysical(), "must start in physical space"); 
	for (int d=0; d<DIM; d++) {
		m_vector[d].forward(); 
	}
}

void Vector::inverse() {
	CHECK(isFFC(), "must start in FFC space"); 
	for (int d=0; d<DIM; d++) {
		m_vector[d].inverse(); 
	}
}

void Vector::inverse(Vector& vector) const {
	CHECK(isFFC(), "must start in FFC space"); 
	for (int d=0; d<DIM; d++) {
		vector[d] = (*this)[d]; 
		vector[d].inverse(); 
	}
}

Vector Vector::cross(const Vector& a_v) const {
	CH_TIMERS("cross product"); 
	CHECK(isFFC() && a_v.isFFC(), "must start in FFC space"); 

	Vector u; 
	Vector v; 
	inverse(u); 
	a_v.inverse(v); 

	#pragma omp parallel
	{
		array<cdouble,DIM> tmp; // store components of u_i to do in place 
		#pragma omp for 
		for (int i=0; i<u.size(); i++) {
			for (int d=0; d<DIM; d++) {
				tmp[d] = u[d][i]; 
			}
			u[0][i] = tmp[1]*v[2][i] - tmp[2]*v[1][i]; 
			u[1][i] = tmp[2]*v[0][i] - tmp[0]*v[2][i]; 
			u[2][i] = tmp[0]*v[1][i] - tmp[1]*v[0][i]; 
		}

	}
	u.forward(); 
	return u; 
}

Vector Vector::curl() const {
	CH_TIMERS("curl"); 
	CHECK(isFFC(), "must start in FFC space"); 

	Vector curl(dims(), false); 

	#pragma omp parallel 
	{
		cdouble imag(0,1); 
		cdouble* D[DIM]; 
		for	(int d=0; d<DIM; d++) {
			D[d] = new cdouble[m_N[2]]; 
		}

		#pragma omp for 
		for (int i=0; i<m_N[0]; i++) {
			double m = (*this)[0].freq(i, 0); 
			for (int j=0; j<m_N[1]; j++) {
				double n = (*this)[0].freq(j, 1); 
				for (int d=0; d<DIM; d++) {
					m_vector[d].ddz(i,j, D[d]); 
				}
				for (int k=0; k<m_N[2]; k++) {
					curl[0](i,j,k) = imag*n*(*this)[2](i,j,k) 
						- D[1][k]; 
					curl[1](i,j,k) = D[0][k] - imag*m*(*this)[2](i,j,k); 
					curl[2](i,j,k) = imag*m*(*this)[1](i,j,k) - 
						imag*n*(*this)[0](i,j,k); 
				}
			}
		}

		for (int d=0; d<DIM; d++) {
			delete D[d]; 
		}
	}

	return curl; 
}

Scalar Vector::divergence() const {
	CH_TIMERS("divergence"); 
	CHECK(isFFC(), "must start in FFC space"); 

	// return in FFC space 
	Scalar div(dims(), false);

	#pragma omp parallel 
	{
		cdouble imag(0,1); 
		cdouble* D = new cdouble[m_N[2]]; 
		#pragma omp for 
		for (int i=0; i<m_N[0]; i++) {
			double m = (*this)[0].freq(i,0); 
			for (int j=0; j<m_N[1]; j++) {
				double n = (*this)[0].freq(j,1); 
				m_vector[2].ddz(i,j, D); 
				for (int k=0; k<m_N[2]; k++) {
					div(i,j,k) = imag*m*(*this)[0](i,j,k) 
						+ imag*n*(*this)[1](i,j,k) 
						+ D[k]; 				
				}
			}
		} 
		delete D;
	} 

	return div; 
}

Vector Vector::laplacian() const {
	ERROR("not defined"); 
}

bool Vector::isFFC() const {
	bool ffc = true; 
	for (int d=0; d<DIM; d++) {
		if (!m_vector[d].isFFC()) ffc = false; 
	}
	return ffc; 
}

bool Vector::isPhysical() const {
	bool phys = true; 
	for (int d=0; d<DIM; d++) {
		if (!m_vector[d].isPhysical()) phys = false; 
	}
	return phys; 
}

int Vector::size() const {return m_vector[0].size(); }
array<int,DIM> Vector::dims() const {return m_vector[0].dims(); }

void Vector::setFFC() {
	for (int d=0; d<DIM; d++) {
		m_vector[d].setFFC(); 
	}
}

void Vector::setPhysical() {
	for (int d=0; d<DIM; d++) {
		m_vector[d].setPhysical(); 
	}
}
#include "DataObjects.H" 
#include <iostream> 

Scalar::Scalar() {} 

Scalar::Scalar(array<int,DIM> N, bool physical) {
	init(N, physical); 
}

Scalar::~Scalar() {
	// decrement number of scalars allocated 
	m_nscalars--; 
}

Scalar::Scalar(const Scalar& scalar) {
	init(scalar.dims(), scalar.isPhysical()); 

	// copy local data 
	for (int i=0; i<m_size; i++) {
		m_data[i] = scalar[i]; 
	}
}

void Scalar::operator=(const Scalar& scalar) {
	if (size() == 0) init(scalar.dims(), scalar.isPhysical()); 

	for (int i=0; i<scalar.size(); i++) {
		m_data[i] = scalar[i]; 
	}
}

void Scalar::init(array<int,DIM> N, bool physical) {
	// store dimensions 
	m_N = N; 
	m_size = 1; 
	for (int i=0; i<DIM; i++) {
		m_size *= m_N[i]; 
	}

	// setup FFC flag 
	if (physical) m_ffc = false; 
	else m_ffc = true; 

	// update number of allocated scalars 
	m_nscalars++; 

	// resize m_data 
	m_data.resize(m_size); 

	// setup transform plans 
	m_fft_x.init(m_N[0], m_N[2]*m_N[1], &m_data[0]); 
	m_fft_y.init(m_N[1], m_N[2], &m_data[0]); 
	m_cheb.init(m_N[2], 1, &m_data[0]);
}

cdouble& Scalar::operator[](array<int,DIM> ind) {
	return m_data[index(ind)]; 
}

const cdouble& Scalar::operator[](array<int,DIM> ind) const {
	return m_data[index(ind)]; 
}

cdouble& Scalar::operator[](int ind) {return m_data[ind]; }
const cdouble& Scalar::operator[](int ind) const {return m_data[ind]; }
cdouble& Scalar::operator()(int a_i, int a_j, int a_k) {
	return m_data[index({a_i, a_j, a_k})]; 
}
const cdouble& Scalar::operator()(int a_i, int a_j, int a_k) const {
	return m_data[index({a_i, a_j, a_k})]; 
}
array<int,DIM> Scalar::dims() const {return m_N; } 
int Scalar::size() const {return m_data.size(); }

void Scalar::forward() {
	CHECK(isPhysical(), "already in FFC space"); 

	// z 
	#pragma omp parallel for 
	for (int j=0; j<m_N[1]; j++) {
		for (int i=0; i<m_N[0]; i++) {
			m_cheb.transform(&(*this)(i,j,0), -1); 
		}
	}

	// y 
	#pragma omp parallel for 
	for (int k=0; k<m_N[2]; k++) {
		for (int i=0; i<m_N[0]; i++) {
			m_fft_y.transform(&(*this)(i,0,k), -1); 
		}
	}

	// x
	#pragma omp parallel for 
	for (int k=0; k<m_N[2]; k++) {
		for (int j=0; j<m_N[1]; j++) {
			m_fft_x.transform(&(*this)(0,j,k), -1); 
		}
	}
	zeroHighModes(); 
	setFFC(); 
}

void Scalar::inverse() {
	CHECK(isFFC(), "already in physical space"); 

	// x
	#pragma omp parallel for 
	for (int k=0; k<m_N[2]; k++) {
		for (int j=0; j<m_N[1]; j++) {
			m_fft_x.transform(&(*this)(0,j,k), 1);  
		}
	}

	// y 
	#pragma omp parallel for 
	for (int k=0; k<m_N[2]; k++) {
		for (int i=0; i<m_N[0]; i++) {
			m_fft_y.transform(&(*this)(i,0,k), 1); 
		}
	}

	// z 
	#pragma omp parallel for 
	for (int j=0; j<m_N[1]; j++) {
		for (int i=0; i<m_N[0]; i++) {
			m_cheb.transform(&(*this)(i,j,0), 1); 
		}
	}
	setPhysical(); 
}

void Scalar::inverse(Scalar& scalar) const {
	CHECK(isFFC(), "already in physical space"); 

	scalar = (*this); 
	scalar.inverse(); 
}

double Scalar::freq(int ind, int d) const {
#ifndef NCHECK 
	if (d >= 2) ERROR("frequency not defined for z dimension"); 
#endif
	if (ind <= m_N[d]/2) return (double)ind; 
	else return -1.*(double)m_N[d] + (double)ind;  
}

Vector Scalar::gradient() const {
	CHECK(isFFC(), "must start in FFC space"); 

	// return in FFC space 
	Vector grad(m_N, false); 

	cdouble imag(0, 1.); 

	#pragma omp parallel for 
	for (int i=0; i<m_N[0]; i++) {
		double m = freq(i, 0); 
		for (int j=0; j<m_N[1]; j++) {
			double n = freq(j, 1); 
			for (int k=0; k<m_N[2]; k++) {
				grad[0](i,j,k) = imag*m*(*this)(i,j,k); 
				grad[1](i,j,k) = imag*n*(*this)(i,j,k); 
			}
			m_cheb.deriv(&(*this)(i,j,0), &grad[2](i,j,0)); 
		}
	}

	return grad; 
}

Scalar Scalar::laplacian() const {
	CHECK(isFFC(), "must start in FFC space"); 

	// return in FFC space 
	Scalar lap(m_N, false);

	#pragma omp parallel 
	{
		cdouble* d = new cdouble[m_N[2]]; 
		cdouble* d2 = new cdouble[m_N[2]]; 

		#pragma omp for 
		for (int i=0; i<m_N[0]; i++) {
			double m2 = pow(freq(i,0), 2); 
			for (int j=0; j<m_N[1]; j++) {
				double n2 = pow(freq(j,1), 2); 

				m_cheb.deriv(&(*this)(i,j,0), d); 
				m_cheb.deriv(d, d2); 
				for (int k=0; k<m_N[2]; k++) {
					lap(i,j,k) = (-m2 - n2)*(*this)(i,j,k) + d2[k]; 
				}
			}
		}

		delete d; 
		delete d2; 
	}

	return lap; 
}

void Scalar::invert_laplacian(double a, double b) {
	ERROR("not defined"); 
}

void Scalar::memory() const {
	cout << "memory requirement = " << m_nscalars*m_size*sizeof(cdouble)/1e9 << 
		" GB" << endl; 
}

bool Scalar::isPhysical() const {return !m_ffc; }
bool Scalar::isFFC() const {return m_ffc; } 
void Scalar::setPhysical() {m_ffc = false; }
void Scalar::setFFC() {m_ffc = true; } 

void Scalar::ddz(int a_i, int a_j, cdouble* output) const {
	m_cheb.deriv(&(*this)(a_i, a_j, 0), output); 
}

void Scalar::zeroHighModes() {
#ifdef ZERO 
	array<int,DIM> half; 
	for (int i=0; i<DIM; i++) {
		half[i] = m_N[i]/2; 
	}

	// zero out x 
	array<int,DIM> ind = {half[0], 0, 0}; 
	for (ind[1]=0; ind[1]<m_N[1]; ind[1]++) {
		for (ind[2]=0; ind[2]<m_N[2]; ind[2]++) {
			(*this)[ind]=0; 
		}
	}

	// zero out y 
	ind = {0, half[1], 0}; 
	for (ind[2]=0; ind[2]<m_N[2]; ind[2]++) {
		for (ind[0]=0; ind[0]<m_N[0]; ind[0]++) {
			(*this)[ind] = 0; 
		}
	}
#endif
}

int Scalar::index(array<int,DIM> ind) const {
	int index = ind[2] + m_N[2]*ind[1] + m_N[2]*m_N[1]*ind[0]; 
	CHECK((bool)(index < size()) && (bool)(index >= 0), "index out of range (" + 
		to_string(index) + "/" + to_string(size()) + ")");  
	// return ind[0] + m_N[0]*ind[1] + m_N[0]*m_N[1]*ind[2]; 
	return index; 
}

int Scalar::m_nscalars=0; 
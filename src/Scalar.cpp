#include "Scalar.H" 
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
	m_fft_x.init(m_N[0], 1, &m_data[0]); 
	m_fft_y.init(m_N[1], m_N[0], &m_data[0]); 
	m_cheb.init(m_N[2], m_N[0]*m_N[1], &m_data[0]); 
}

cdouble& Scalar::operator[](array<int,DIM> ind) {
	return m_data[index(ind)]; 
}

cdouble Scalar::operator[](array<int,DIM> ind) const {
	return m_data[index(ind)]; 
}

cdouble& Scalar::operator[](int ind) {return m_data[ind]; }
cdouble Scalar::operator[](int ind) const {return m_data[ind]; }
array<int,DIM> Scalar::dims() const {return m_N; } 
int Scalar::size() const {return m_data.size(); }

double Scalar::freq(int ind, int d) const {
	if (ind <= m_N[d]/2) return (double)ind; 
	else return -1.*(double)m_N[d] + (double)ind;  
}

void Scalar::forward() {
	CHECK(isPhysical(), "already in FFC space"); 

	// z 
	#pragma omp parallel for 
	for (int j=0; j<m_N[1]; j++) {
		for (int i=0; i<m_N[0]; i++) {
			m_cheb.transform(&m_data[i + j*m_N[0]], -1); 
		}
	}

	// y 
	#pragma omp parallel for 
	for (int k=0; k<m_N[2]; k++) {
		for (int i=0; i<m_N[0]; i++) {
			m_fft_y.transform(&m_data[i + k*m_N[0]*m_N[1]], -1); 
		}
	}

	// x
	#pragma omp parallel for 
	for (int k=0; k<m_N[2]; k++) {
		for (int j=0; j<m_N[1]; j++) {
			m_fft_x.transform(&m_data[j*m_N[0] + k*m_N[0]*m_N[1]], -1); 
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
			m_fft_x.transform(&m_data[j*m_N[0] + k*m_N[0]*m_N[1]], 1); 
		}
	}

	// y 
	#pragma omp parallel for 
	for (int k=0; k<m_N[2]; k++) {
		for (int i=0; i<m_N[0]; i++) {
			m_fft_y.transform(&m_data[i + k*m_N[0]*m_N[1]], 1); 
		}
	}

	// z 
	#pragma omp parallel for 
	for (int j=0; j<m_N[1]; j++) {
		for (int i=0; i<m_N[0]; i++) {
			m_cheb.transform(&m_data[i + j*m_N[0]], 1); 
		}
	}
	setPhysical(); 
}

void Scalar::inverse(Scalar& scalar) const {
	scalar = (*this); 
	scalar.inverse(); 
}

// Vector Scalar::gradient() const {

// }

Scalar Scalar::laplacian() const {

}

void Scalar::invert_laplacian(double a, double b) {

}

void Scalar::memory() {
	cout << "memory requirement = " << m_nscalars*m_size*sizeof(cdouble)/1e9 << 
		" GB" << endl; 
}

bool Scalar::isPhysical() const {return !m_ffc; }
bool Scalar::isFFC() const {return m_ffc; } 
void Scalar::setPhysical() {m_ffc = false; }
void Scalar::setFFC() {m_ffc = true; } 

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
	return ind[0] + m_N[0]*ind[1] + m_N[0]*m_N[1]*ind[2]; 
}

int Scalar::m_nscalars=0; 
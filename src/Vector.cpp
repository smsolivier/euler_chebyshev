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
	for (int d=0; d<DIM; d++) {
		m_vector[d].forward(); 
	}
}

void Vector::inverse() {
	for (int d=0; d<DIM; d++) {
		m_vector[d].inverse(); 
	}
}

void Vector::inverse(Vector& vector) const {
	for (int d=0; d<DIM; d++) {
		vector[d] = (*this)[d]; 
		vector[d].inverse(); 
	}
}

Vector Vector::cross(const Vector& v) const {
	ERROR("not defined"); 
}

Vector Vector::curl() const {
	ERROR("not defined"); 
}

Scalar Vector::divergence() const {
	ERROR("not defined"); 
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
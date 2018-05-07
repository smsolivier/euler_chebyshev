#ifndef __VECTOR_H__ 
#define __VECTOR_H__

/// represent a vector variable in DIM dimensions 
class Vector {
public:
	/// default constructor 
	Vector(); 
	/// set dimensions in DIM dimensions. defaults to buidling in FFC space 
	Vector(array<int,DIM> N, bool physical=false); 
	/// initialize 
	void init(array<int,DIM> N, bool physical=false); 

	/// index into components of vector 
	Scalar& operator[](int ind); 
	const Scalar& operator[](int ind) const; 

	/// transform each component in place to FFC space 
	void forward(); 
	/// transform each component to physical space 
	void inverse();
	/// transform each component out of place to physical space 
	void inverse(Vector& vector) const; 

	/// cross product of two vectors 
	/** does cross product in physical space **/ 
	Vector cross(const Vector& v) const; 
	/// curl: \f$ \nabla \times \vec{f} \f$ 
	Vector curl() const; 
	/// divergence: \f$ \nabla \cdot \vec{f} \f$ 
	Scalar divergence() const; 
	/// vector laplacian: \f$ \nabla \cdot \nabla \vec{f} \f$ 
	Vector laplacian() const; 

	/// query if (*this) is in FFC space 
	bool isFFC() const; 
	/// query if (*this) is in physical space 
	bool isPhysical() const; 

	/// get size of each scalar 
	int size() const; 
	/// get dimensions of scalar 
	array<int,DIM> dims() const; 
private:
	/// flag (*this) as in FFC space 
	void setFFC(); 
	/// flag (*this) as in physical space 
	void setPhysical(); 

	/// size in DIM dimensions 
	array<int,DIM> m_N; 
	/// store DIM scalars 
	array<Scalar,DIM> m_vector; 
}; 

#endif
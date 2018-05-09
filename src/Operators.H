#ifndef __OPERATORS_H__
#define __OPERATORS_H__

Vector operator+(const Vector& a, const Vector& b); 
Vector operator-(const Vector& a, const Vector& b); 
Vector operator*(const double alpha, const Vector& a); 

Scalar operator*(const double alpha, const Scalar& a); 

#endif
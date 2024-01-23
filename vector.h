#pragma once
#ifndef VECTOR_H
#define VECTOR_H
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>

using namespace std;

class Matrix; // forward declaration for matrix-vector multiplication

class Vector {
private:
	double* vecValue;
	int vecElementCount;

public:
	// parameterized constructor
	explicit Vector(int size) : vecElementCount(size) {
		vecValue = new double[vecElementCount]; 
	}
	// destructor
	~Vector()
	{
		delete[] vecValue;
	}

	// user defined copy constructor
	Vector(const Vector& vec) : vecElementCount(vec.vecElementCount) {
		vecValue = new double[vecElementCount];
		copy(vec.vecValue, vec.vecValue + vecElementCount, vecValue);
	}

	// user defined move constructor
	Vector(Vector&& other) noexcept
		:vecValue{ other.vecValue },
		vecElementCount{ other.vecElementCount }
	{
		other.vecValue = nullptr;
	}

	// user defined move assignment
	Vector& operator= (Vector&& other) noexcept
	{
		if (&other == this)
		{
			return *this;
		}
		delete[] vecValue;
		 
		vecValue = other.vecValue;
		vecElementCount = other.vecElementCount;

		other.vecValue= nullptr;

		return *this;
	}

	// get the value at a certain index
	double operator[] (int index) const { // "const" enforces immutability to the compiler of the calling object
		return *(vecValue + index);
	}

	void set(int index, double value); // change elements at a given position in the vector
	void fill(double value); // populate the entire vector with value
	void zeros();
	int getSize() const; // returns number of elements in vector

	// Copy assignment operator
	Vector& operator= (const Vector& a) {
		if (this == &a) // parity check
		{
			return *this;
		}

		delete[] vecValue;

		vecElementCount = a.vecElementCount;
		vecValue = new double[vecElementCount]; // allocate on to heap

		copy(a.vecValue, a.vecValue + vecElementCount, vecValue); // copy over data

		

		return *this;
	}

	friend Vector operator* (const Matrix&, const Vector&); // matrix-vector multiplication
	friend Vector operator* (const Vector&, const double multiplier); // scalar-vector multiplication
	friend ostream& operator<< (ostream&, const Vector&); // stream-operator overload for vector display
	friend Vector operator+ (const Vector&,const Vector&); // vector addition
	friend Vector operator/ (const Vector&, const double divisor); // scalar-reciprocal vector multiplication
};

Vector operator/ (const Vector& A, const double divisor) {
	if (!(divisor == 0)){
		return A * (1 / divisor);
	} else
	{
		throw domain_error("Division by zero not allowed.");
	}
}

Vector operator* (const Vector& A, const double multiplier) {

	Vector product(A.getSize());

	for (int i{ 0 }; i < A.getSize(); i++) {
		product.set(i, A[i] * multiplier);
	}
	
	return product;
}

Vector operator+ (const Vector& A, const Vector& B) {

	Vector sum(A.getSize());
	for (int e{ 0 }; e < A.getSize(); e++) {
		sum.set(e, A[e] + B[e]);
	}
	return sum;
}

ostream& operator<< (ostream& out, const Vector& A) {
	for (int i{ 0 }; i < A.vecElementCount; i++) {
		out << left << setw(4) << A[i] << endl;
	}
	out << endl;
	return out;
}

int Vector::getSize() const {
	return vecElementCount;
}

void Vector::zeros() {
	fill(0);
}

void Vector::fill(double value)
{
	for (int i{ 0 }; i < vecElementCount; i++)
	{
		vecValue[i] = value;
	}
}

void Vector::set(int index, double value) {
	vecValue[index] = value;
}

#endif

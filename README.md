LinLib
======

JavaScript Linear Algebra library

Authored while studying for Linear Algebra class.

Goal
====

The goal of this project to to make a JavaScript based Linear Algebra library that is:

1. Immutable (All objects cannot change state)
2. Flexible (Linear Algebra primitives seperated from helper functions to solve problems)
3. Chainable (All object operations that "mutate" return objects that can further be operated on)
4. Encapsulated (All primitives must be constructed through string parsing, primitive internal data is encapsulated)

Since these goals inherently cause performance decrease, efficiency is not the #1 goal of this library. Ease of reasoning about the library is the priority.


Primitives
==========

Currently the library supports these primitives:

1. Fraction (rational number)
2. Vector
3. Matrix


Vectors and Matrices internally have Fraction objects as their data.

Future versions will have a RationalTerm, and RationalPolynomial primitive such that algebraic expressions can be operated on. RationalPolynomial will take the place of Fraction as the internal object of Vector and Matrix.

Interface Functions
===================

Functions that are exposed from the library under the namespace LinLib are as follows:

	{
	'parseMatrix' : parseMatrix,
	'parseVector' : parseVector,
	'parseVectors' : parseVectors,
	'parseFraction' : parseFraction,
	'changeOfBase' : changeOfBaseMatrix,
	'rowVecsToMatrix' : rowVecsToMatrix,
	'colVecsToMatrix' : colVecsToMatrix,
	'identity' : identity,
	'zero' : zero
	}

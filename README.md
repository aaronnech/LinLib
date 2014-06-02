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

Of course, the library primitives, once constructed, will have operations exposed to them. (Such as transpose() for Matrix).

Usage
=====

The library is non-mutable and chainable. So an example computation could be:

	LinLib.colVecsToMatrix(LinLib.parseMatrix('{{2,0},{0,2}}').invert().add(LinLib.identity(2)).transpose().columnVectors()).transpose().invert().invert().toString();
	
This computation converts the column vectors of:

	{2, 0}
	{0, 2}

Inverted, summed with the n = 2 identity matrix, and transposed.

It then converts those column vectors back to a Matrix ('cause why not), and takes the inverse of the inverse of the transpose of that Matrix. It then takes that result and converts it to a string for printing.

Obviously this example is not practical, but it serves as an example of how to use this library. To see all available commands for library primitives, check out linlib.js. It is fully documented with JavaDoc notation.

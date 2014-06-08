/**
 * Javascript Linear Algebra library
 *
 *
 * @author Aaron Nech
 * @license http://opensource.org/licenses/MIT MIT
 *
 * I am not liable for any use, distribution, or implications of this software.
 * It is provided free of charge for any use, modification, or redistribution.
 *
 * 
 */

var LinLib = (function() {

	//constants
	var ZERO = new Fraction(0, 1);
	var ONE = new Fraction(1, 1);
	var NEGATIVE_ONE = new Fraction(-1, 1);
	

	//Operation Recording
	var recordings = [];

	/**
	 * Starts a new recording.
	 * 
	 * @param  closure filter an optional parameter that will return
	 *         true or false on a operation signature. Returning true
	 *         indicates that this operation will be recorded, false
	 *         will not.
	 * @param  closure callback (if defined) will be called with each recording that
	 *         passes filter as it is called
	 * @return int id where id is the id of this recording
	 */
	var startRecording = function(filter, callback) {
		if(typeof filter == 'undefined')
			filter = function(x) {return true};
		if(typeof callback == 'undefined')
			callback = function(x) {};
		recordings.push({
			'filter' : filter,
			'callback' : callback,
			'ops' : []
		});
		return recordings.length - 1;
	};

	/**
	 * Stops a recording.
	 * 
	 * @param  int id the id of the recording to stop.
	 * @return array A where A is an array of operations
	 *         captured by the recording in the following format:
	 *         {
	 *         		sig : operation signature
	 *         		obj : object operated on
	 *         		args : arguments passed to operation
	 *         }
	 */
	var stopRecording = function(id) {
		return recordings.splice(id, 1)[0].ops;
	};

	/**
	 * Attempts to record an operation on all open
	 * recordings. (Sig must be accepted by a recording's
	 * filter)
	 * 
	 * @param  string sig the signature of the operation
	 * @param  Object? obj the object being operated on
	 * @param  array[Object?] args the arguments passed to the operation
	 * @return true if the operation was recorded by at least one recording,
	 *         false otherwise.
	 */
	var record = function(sig, obj, args) {
		if(recordings.length > 0) {
			for(var i = 0; i < recordings.length; i++) {
				if(recordings[i].filter(sig)) {
					var data = {
						'sig' : sig,
						'obj' : obj,
						'args' : args
					};
					recordings[i].ops.push(data);
					recordings[i].callback(data);
				}
			}
			return true;
		}
		return false;
	};

	/**
	 * Parses a Fraction from a string ignoring whitespace
	 *
	 * @requires string is well-formed: it is either an integer or
	 *           a/b format where a, b are integers.
	 * @param  string s string to parse
	 * @return Fraction f such that f represents s
	 */
	var parseFraction = function(s) {
		if(s == '0')
			return ZERO;
		if(s.indexOf('/') != -1) {
			var parts = s.split('/');
			return new Fraction(parseInt(parts[0]), parseInt(parts[1]));
		}
		return new Fraction(parseInt(s), 1);
	}

	/**
	 * Parses a Vector from a string ignoring whitespace
	 *
	 * @requires  string is well-formed: it is a comma seperated
	 *            list of Fractions able to be parsed by parseFraction
	 *            enclosed with brackets {}.
	 * @param  string s string to be parsed
	 * @return Vector v such that v represents s
	 */
	var parseVector = function(s) {
		s = s.replace(/\s+/g, '').substring(1, s.length-1);
		return new Vector(s.split(',').map(function(x) {
				return parseFraction(x);
		}));
	};

	/**
	 * Parses an array Vectors from an array of strings ignoring whitespace
	 *
	 * @requires  strings in arr are well-formed: they are a comma seperated
	 *            list of Fractions able to be parsed by parseFraction
	 *            enclosed with brackets {}.
	 * @param  array[string] arr the array of strings to be parsed
	 * @return array[Vector] a with each element v such that v represents
	 *         each corresponding element s of arr.
	 */
	var parseVectors = function(arr) {
		return arr.map(function(x) {
			return parseVector(x);
		});
	};

	/**
	 * Parses a Matrix from a string ignoring whitespace
	 *
	 * @requires string is well-formed: it is a comma seperated
	 *           list of row-vectors able to be parsed by parseVector
	 *           enclosed with brackets {}
	 * @param  string s string to parse
	 * @return Matrix m such that m represents s
	 */
	var parseMatrix = function(s) {
		s = s.replace(/\s+/g, '');
		if(s.charAt(0) != '{' || s.charAt(s.length - 1) != '}')
			throw "Illegal Array: Must start with '{' and end with '}'"
		var noBrackets = s.substring(1, s.length-1);
		return new Matrix(noBrackets.split('},{').map(function(x){
			return '{' + x.replace('}', '').replace('{','') + '}';
		}).map(parseVector).map(function(vec) {
			return vec.toArray();
		}));
	};

	/**
	 * Converts an array of row vectors to a matrix
	 *
	 * @requires vecs contains vectors of equal dimension
	 *           vecs.length > 0
	 * @param  array[Vector] vecs array of Vectors that are the
	 *         rows of the result
	 * @return Matrix A such that the rows of A are the vectors in vecs
	 */
	var rowVecsToMatrix = function(vecs) {
		return new Matrix(vecs.map(function(x) {
			return x.toArray();
		}));
	};

	/**
	 * Converts an array of column vectors to a matrix
	 *
	 * @requires vecs contains vectors of equal dimension
	 *           vecs.length > 0
	 * @param  array[Vector] vecs array of Vectors that are the
	 *         columns of the result
	 * @return Matrix A such that the columns of A are the vectors in vecs
	 */
	var colVecsToMatrix = function(vecs) {
		return rowVecsToMatrix(vecs).transpose();
	};

	/**
	 * Applies Least Squared Approximation to a 
	 * (sometimes) unsolvable system of equations representing
	 * a 2d line through all the points passed in vector <x,y> form.
	 *
	 * @requires  points is an array of 2 dimensional vectors
	 * @param  array[Vector] points The points to approximate
	 * @return Vector V such that the first element of V is the unknown m,
	 *          and the second value is the unknown b in the equation y = mx + b.
	 */
	var bestFitLine = function(points) {
		// Least Squared Approximation:
		// A^TA * Xaprox = A^Tb
		//
		// In this case, A is the matrix
		// formed by the coeffiences of the
		// line equation, and b is the 
		// vector of f(x) = y coordinates.

		var b = new Matrix([points.map(function(x) {
			return x.get(1);
		})]).transpose();

		// Get a matrix A such that
		// every row in A has two elements,
		// the first will be the coefficent x for m,
		// the second will be the coefficent 1 for b.
		var a = new Matrix(points.map(function(x) {
			return [x.get(0), ONE];
		}));
		
		var aT = a.transpose();

		// solve
		// A^TA * Xaprox = A^Tb via gauss elimination
		var result = aT.mult(a).augmentCol(aT.mult(b).columnVectors()[0]).gaussJordan();

		//return the resultant augmented column vector yielding {Maprox,Baprox}
		return result.columnVectors()[result.colCount() - 1];
	};

	/**
	 * Generates the identity matrix with dimension n.
	 * 
	 * @param  int n dimension of identity matrix
	 * @return Matrix m such that m is the identity matrix with
	 *         dimension n
	 */
	var identity = function(n) {
		var result = [];
		for(var i = 0; i < n; i++) {
			result.push([]);
			for(var j = 0; j < n; j++) {
				if(j == i) {
					result[i].push(ONE);
				} else {
					result[i].push(ZERO);
				}
			}
		}
		return new Matrix(result);
	}

	/**
	 * Generates the zero matrix that is nxm
	 * 
	 * @param  int n number of rows
	 * @param  int m number of columns
	 * @return Matrix m such that m is the zero matrix and is nxm
	 */
	var zero = function(n, m) {
		var result = [];
		for(var i = 0; i < n; i++) {
			result.push([]);
			for(var j = 0; j < m; j++) {
				result[i].push(ZERO);
			}
		}
		return new Matrix(result);
	}

	/**
	 * Finds the change of base matrix A such that
	 * A changes a vector with coordinates in respect to alpha
	 * to coordinates in respect to beta
	 * 
	 * @param  array[Vector] alpha array of Vectors that are a basis for alpha
	 * @param  array[Vector] beta array of Vectors that are a basis for beta
	 * @return Matrix A such that A is the change of basis matrix from alpha to beta
	 */
	var changeOfBaseMatrix = function(alpha, beta) {
		var baseMatrix = colVecsToMatrix(beta);
		results = [];
		for(var i = 0; i < alpha.length; i++) {
			var augmented = baseMatrix.augmentCol(alpha[i]);
			var solved = augmented.gaussJordan().columnVectors();
			results.push(solved[solved.length - 1]);
		}
		return colVecsToMatrix(results);
	};

	/**
	 * Returns true if a vector vec is orthogonal to a subspace formed
	 * by basis.
	 * 
	 * @param  Vector vec the vector to check for orthogonality
	 * @param  array[Vector] basis the basis of the subspace to check
	 * @return true if vec is orthogonal to the subspace formed by basis,
	 *         false otherwise.
	 */
	var orthogonalToSubspace = function(vec, basis) {
		for(var i = 0; i < basis.length; i++) {
			if(!basis[i].dot(vec).isZero())
				return false;
		}
		return true;
	};

	/**
	 * Returns a Vector V such that V is the projection
	 * of vec on to a subspace formed by an orthogonal basis
	 *
	 * @requires  basis is an orthogonal basis of a subspace
	 * @param  Vector vec the vector project
	 * @param  array[Vector] basis the orthogonal basis of the subspace to project on to
	 * @return Vector V such that V is the projection of vec onto basis
	 */
	var projectSubspace = function(vec, basis) {
		// Let S be a non-zero subspace of R^n with orthogonal
		// basis basis, and let vec be of R^n.
		// The projection of vec on to S is the sum of
		// the projections of vec on the elements of
		// basis.
		var result = vec.proj(basis[0]);
		for(var i = 1; i < basis.length; i++) {
			result = result.add(result.proj(basis[i]));
		}
		return result;
	};


	/**
	 * Zips each ith element in each array in an array
	 * together.
	 *
	 * @param arrays the arrays to zip together
	 * @return Array A such that A is a single array with each
	 *         element i corresponding to an array with
	 *         elements that are the ith elements of each array in arrays
	 */
	var zip = function(arrays) {
	    return arrays[0].map(function(_,i){
	        return arrays.map(function(array){return array[i]})
	    });
	}

	/**
	 * Fraction is a high level representation of the math primitive Fraction (rational number).
	 * 
	 * @requires den != 0
	 * @param int num, the numerator of the fraction
	 * @param int den, the denominator of the fraction
	 */
	function Fraction(num, den) {
		var self = this;
		var numerator = 0;
		var denominator = 0;
		if(num < 0 && den < 0) {
			numerator = -num;
			denominator = -den;
		} else if(num < 0 || den < 0) {
			numerator = -Math.abs(num);
			denominator = Math.abs(den);
		} else {
			numerator = num;
			denominator = den;
		}

		if(denominator == 0)
			throw "Divide by zero! num: " + numerator + " den: " + denominator;

		//private
		var gcd = function(a, b) {
			if(!b)
				return a;
			return gcd(b, a % b);
		};
		var normalize = function() {
			var gc = gcd(Math.abs(numerator), Math.abs(denominator));
			numerator /= gc;
			denominator /= gc;
		};
		normalize();


		//public
		
		/**
		 * Gets the numerator of the fraction
		 * 
		 * @return int n such that n is the numerator
		 */
		self.getNum = function() {
			record('Fraction.getNum', self);
			return numerator;
		};

		/**
		 * Gets the denominator of the fraction
		 * 
		 * @return int n such that n is the numerator
		 */
		self.getDen = function() {
			record('Fraction.getDen', self);
			return denominator;
		};

		/**
		 * Adds a fraction f to the current fraction returning the result
		 * 
		 * @param Fraction f the fraction to add
		 * @return Fraction a such that a = self + f
		 */
		self.add = function(f) {
			record('Fraction.add', self, [f]);
			if(f.isZero())
				return self;
			return new Fraction(self.getNum() * f.getDen() + f.getNum() * self.getDen(),
				f.getDen() * self.getDen());
		};

		/**
		 * Subtracts a fraction f from the current fraction returning the result
		 * 
		 * @param Fraction f the fraction to subtract
		 * @return Fraction a such that a = self - f
		 */
		self.sub = function(f) {
			record('Fraction.sub', self, [f]);
			if(f.isZero())
				return self;
			return new Fraction(self.getNum() * f.getDen() - f.getNum() * self.getDen(),
				f.getDen() * self.getDen());
		};

		/**
		 * Multiplies a fraction f with the current fraction returning the result
		 * 
		 * @param Fraction f the fraction to multiply
		 * @return Fraction a such that a = self * f
		 */
		self.mult = function(f) {
			record('Fraction.mult', self, [f]);
			if(f.isZero())
				return ZERO;
			return new Fraction(self.getNum() * f.getNum(),
				f.getDen() * self.getDen());
		};

		/**
		 * Divides self by the fraction f returning the result
		 * 
		 * @param Fraction f the fraction to divide by
		 * @return Fraction a such that a = self / f
		 */
		self.div = function(f) {
			record('Fraction.div', self, [f]);
			if(self.isZero())
				return ZERO;
			return new Fraction(self.getNum() * f.getDen(),
				f.getNum() * self.getDen());
		};

		/**
		 * Inverts self returning the result
		 * 
		 * @return Fraction a such that a = self^-1
		 */
		self.invert = function() {
			record('Fraction.invert', self);
			if(self.isZero())
				return ZERO;
			return new Fraction(self.getDen(), self.getNum());
		};

		/**
		 * Evaluates self, returning the decimal result
		 * 
		 * @return float a such that a = self (to highest precision possible)
		 */
		self.toDec = function() {
			record('Fraction.toDec', self);
			if(self.isZero())
				return 0;
			return self.getNum() / self.getDen();
		};

		/**
		 * Negates self, returning the result
		 * 
		 * @return float a such that a = -self
		 */
		self.negate = function() {
			record('Fraction.negate', self);
			if(self.isZero())
				return ZERO;
			return new Fraction(-self.getNum(), self.getDen());
		};

		/**
		 * Tests self == 0, returning the result
		 * 
		 * @return boolean b such that b = self == 0
		 */
		self.isZero = function() {
			record('Fraction.isZero', self);
			return self.getNum() == 0;
		};

		/**
		 * Tests if self == f
		 * (defined as num == num and den == den)
		 *
		 * @return boolean b such that b = self == f
		 */
		self.eq = function(f) {
			record('Fraction.eq', self, [f]);
			return self.getNum() == f.getNum() && self.getDen() == f.getDen();
		}

		/**
		 * Converts self to string representation
		 * 
		 * @return string s such that s represents self
		 */
		self.toString = function() {
			record('Fraction.toString', self);
			if(Math.abs(self.getNum()) == Math.abs(self.getDen()))
				return (self.getNum() < 0 ? "-" : "") + "1";
			return self.getNum() + (self.getDen() != 1 ? "/" + self.getDen() : "");
		};
	}

	/**
	 * Vector is a high level representation of the math primitive vector.
	 * 
	 * @requires arr is not empty and contains only Fraction elements
	 * @param (1 dimensional array) arr the array representation of the vector
	 */
	function Vector(arr) {
		var self = this;
		var data = arr;

		//private
		var sizesMatch = function(vec1, vec2) {
			record('Vector.sizesMatch', self, [vec1, vec2]);
			return vec1.dim() == vec2.dim();
		}

		//public
		
		/**
		 * Gets the dimension of the vector
		 * 
		 * @return int n where n is the number of elements in self
		 */
		self.dim = function() {
			record('Vector.dim', self);
			return data.length;
		};

		/**
		 * Gets the ith element in this matrix
		 *
		 * @requires i < self.dim() && i >= 0
		 * @param int i the index of the element to get
		 * @return Fraction f the element at i.
		 */
		self.get = function(i) {
			record('Vector.get', self, [i]);
			if(i >= self.dim() || i < 0)
				throw "Illegal index";
			return data[i];
		};

		/**
		 * Gets the array representation of this vector
		 * 
		 * @return array[Fraction] the array representation of self
		 */
		self.toArray = function() {
			record('Vector.toArray', self);
			return data.slice(0);
		}

		/**
		 * Adds a Vector vec to self returning the result
		 * 
		 * @param Vector vec the vector to add
		 * @return Vector a such that a = self + vec
		 */
		self.add = function(vec) {
			record('Vector.add', self, [vec]);
			if(sizesMatch(vec, self))
				return new Vector(zip([self.toArray(), vec.toArray()]).map(function(x) {
					return x[0].add(x[1]);
				}));
			return null;
		};

		/**
		 * Subtracts a Vector vec from self returning the result
		 * 
		 * @param Vector vec the vector to subtract from self
		 * @return Vector a such that a = self - vec
		 */
		self.sub = function(vec) {
			record('Vector.sub', self, [vec]);
			if(sizesMatch(vec, self))
				return new Vector(zip([self.toArray(), vec.toArray()]).map(function(x) {
					return x[0].sub(x[1]);
				}));
			return null;
		};

		/**
		 * Takes the dot product of self with vec returning the result
		 * 
		 * @param Vector vec the vector to dot self with
		 * @return Fraction a such that a = self (dot) vec
		 */
		self.dot = function(vec) {
			record('Vector.dot', self, [vec]);
			if(sizesMatch(vec, self))
				return zip([self.toArray(), vec.toArray()]).reduce(function(prev, x) {
					return prev.add(x[0].mult(x[1]));
				}, ZERO);
			return null;
		};

		/**
		 * Find the distance^2 between two N dimensional vectors
		 * 
		 * @param Vector vec the vector to find the distance from
		 * @return Fraction a such that a = dist(self, vec)^2
		 */
		self.distSq = function(vec) {
			record('Vector.dist', self, [vec]);
			var dif = self.sub(vec)
			return dif.dot(dif);
		};

		/**
		 * Find the length^2 of self
		 * 
		 * @return Fraction a such that a = ||self||^2
		 */
		self.lengthSq = function() {
			return self.dot(self);
		};

		/**
		 * Finds a vector V that is the projection
		 * of self on to vec.
		 * 
		 * @param  Vector vec the vector to project on to.
		 * @return Vector V such that V = proj(self, vec)
		 */
		self.proj = function(vec) {
			return vec.scale(self.dot(vec).div(vec.lengthSq()));
		};


		/**
		 * Scales self by scalar s
		 * 
		 * @param Fraction s the scalar to scale by
		 * @return Vector a such that a = self * s
		 */		
		self.scale = function(s) {
			record('Vector.scale', self, [s]);
			return new Vector(self.toArray().map(function(x) {
				return x.mult(s);
			}));
		};

		/**
		 * Tests if self == vec
		 * (defined as all elements being equal Fractions)
		 *
		 * @return boolean b such that b = self == vec
		 */
		self.eq = function(vec) {
			record('Vector.eq', self, [vec]);
			return zip([self.toArray(), vec.toArray()]).reduce(function(prev, x) {
				return prev && x[0].eq(x[1]);
			});
		}

		/**
		 * Tests if the current vector contains zero as an element.
		 * 
		 * @return boolean b such that b == true if self has zero as
		 *            an element, false otherwise
		 */
		self.hasZero = function() {
			record('Vector.hasZero', self);
			for(var i = 0; i < data.length; i++) {
				if(data[i].isZero())
					return true;
			}
			return false;
		}

		/**
		 * Converts self to string representation
		 * 
		 * @return string s such that s represents self
		 */
		self.toString = function() {
			record('Vector.toString', self);
			return '{' + data.map(function(x) {
				return x.toString();
			}).join(', ') + '}';
		};
	}

	/**
	 * Matrix is a high level representation of the math primitive matrix.
	 * @requires arr is well-formed with no empty arrays or sub arrays, and consists
	 *           of only Fraction elements
	 * @param (2 dimensional array) arr the array representation of the matrix,
	 *        in row-major order.
	 */
	function Matrix(arr) {
		var self = this;
		var data = arr;

		//private
		
		//public
		
		/**
		 * Tests if self == mat
		 * (defined as all row vectors being equal)
		 *
		 * @return boolean b such that b = self == mat
		 */
		self.eq = function(mat) {
			record('Matrix.eq', self, [mat]);
			return zip([self.rowVectors(), mat.rowVectors()]).reduce(function(prev, x) {
				return prev && x[0].eq(x[1]);
			});
		}
		
		/**
		 * Gets the number of columns of self
		 * 
		 * @return int a where a is the number of columns of self
		 */
		self.colCount = function() {
			record('Matrix.colCount', self);
			return data[0].length;
		};

		/**
		 * Gets the number of rows of self
		 * 
		 * @return int a where a is the number of rows of self
		 */
		self.rowCount = function() {
			record('Matrix.rowCount', self);
			return data.length;
		};

		/**
		 * Gets the array representation of self in row-major form
		 * 
		 * @return array[array[Fraction]] the array representation of self
		 */
		self.toArray = function() {
			record('Matrix.toArray', self);
			result = data.slice(0);
			for(var i = 0; i < result.length; i++) {
				result[i] = result[i].slice(0);
			}
			return result;
		};

		/**
		 * Gets the array of colomn vectors of self
		 * 
		 * @return array[Vector] the array representing columns as vectors
		 */
		self.columnVectors = function() {
			record('Matrix.columnVectors', self);
			result = [];
			for(var i = 0; i < self.colCount(); i++) {
				column = [];
				for(var j = 0; j < self.rowCount(); j++) {
					column.push(data[j][i]);
				}
				result.push(new Vector(column));
			}
			return result;
		};

		/**
		 * Gets the array of row vectors of self
		 * 
		 * @return array[Vector] the array representing rows as vectors
		 */
		self.rowVectors = function() {
			record('Matrix.rowVectors', self);
			result = [];
			for(var i = 0; i < self.rowCount(); i++) {
				result.push(new Vector(data[i].slice(0)));
			}
			return result;
		};

		/**
		 * Adds a column vector to the right side of the matrix as a column
		 * returning the result
		 *
		 * @requires vec.dim() == self.rowCount()
		 * @return Matrix m that augmented matrix, null if not possible to augment
		 */
		self.augmentCol = function(vec) {
			record('Matrix.augmentCol', self, [vec]);
			if(vec.dim() != self.rowCount())
				return null;
			var vecArray = vec.toArray();
			return new Matrix(self.toArray().map(function(x, index) {
				x.push(vecArray[index]);
				return x;
			}));
		};

		/**
		 * Adds a row vector to the bottom of the matrix as a row
		 * returning the result
		 *
		 * @requires vec.dim() == self.colCount()
		 * @return Matrix m that augmented matrix, null if not possible to augment
		 */
		self.augmentRow = function(vec) {
			record('Matrix.augmentRow', self, [vec]);
			if(vec.dim() != self.colCount())
				return null;
			var result = self.toArray();
			result.push(vec.toArray());
			return new Matrix(result);
		};

		/**
		 * Swaps row i with row j in self returning the result
		 *
		 * @requires i < self.rowCount() && i >= 0 && j < self.rowCount() && j >= 0
		 * @param int i the index of the first row to swap
		 * @param int j the index of the second row to swap
		 * @return Matrix m the altered matrix, null if not possible to swap rows
		 */
		self.swapRows = function(i, j) {
			record('Matrix.swapRows', self, [i, j]);
			if(i >= self.rowCount() || i < 0 || j >= self.rowCount() || j < 0)
				return null;
			var copy = self.toArray();
			var tmp = copy[i];
			copy[i] = copy[j];
			copy[j] = tmp;
			return new Matrix(copy);
		};

		/**
		 * Scales row i by s returning the result
		 *
		 * @requires i < self.rowCount() && i >= 0
		 * @param int i the index of row to scale
		 * @param Fraction s scalar to multiply row i by
		 * @return Matrix m the altered matrix, null if not possible to scale row
		 */
		self.scaleRow = function(i, s) {
			record('Matrix.scaleRow', self, [i, s]);
			if(i >= self.rowCount() || i < 0)
				return null;
			var copy = self.toArray();
			copy[i] = new Vector(copy[i]).scale(s).toArray();
			return new Matrix(copy);
		};

		/**
		 * Adds row i to row j returning the result
		 *
		 * @requires i < self.rowCount() && i >= 0 && j < self.rowCount() && j >= 0
		 * @param int i the index of the row to add
		 * @param int j the index of the row to store the resultant row
		 * @return Matrix m the altered matrix, null if not possible to add rows
		 */
		self.addRowToRow = function(i, j) {
			record('Matrix.addRowToRow', self, [i, j]);
			if(i >= self.rowCount() || i < 0 || j >= self.rowCount() || j < 0)
				return null;
			var copy = self.toArray();
			copy[j] = new Vector(copy[j]).add(new Vector(copy[i])).toArray();
			return new Matrix(copy);
		};

		/**
		 * Adds row i scaled by s to row j returning the result
		 *
		 * @requires i < self.rowCount() && i >= 0 && j < self.rowCount() && j >= 0
		 * @param int i the index of the row to add
		 * @param int j the index of the row to store the resultant row
		 * @param Fraction s scalar to multiply row i by
		 * @return Matrix m the altered matrix, null if not possible to add rows
		 */
		self.addScaledRowToRow = function(i, s, j) {
			record('Matrix.addScaledRowToRow', self, [i, s, j]);
			if(i >= self.rowCount() || i < 0 || j >= self.rowCount() || j < 0)
				return null;
			var copy = self.toArray();
			copy[j] = new Vector(copy[j]).add(new Vector(copy[i]).scale(s)).toArray();
			return new Matrix(copy);
		};

		/**
		 * Inverts self returning the result
		 *
		 * @requires self is NxN and non-singular
		 * @return Matrix m the inverted matrix
		 */
		self.invert = function() {
			record('Matrix.invert', self);
			if(self.rowCount() != self.colCount())
				throw "Inverse not defined for non NxN";
			if(self.determinant().isZero())
				throw "Inverse not defined for singular matrix";
			var n = self.rowCount();
			var id = identity(n).rowVectors();
			var result = self;
			for(var i = 0; i < n; i++) {
				result = result.augmentCol(id[i]);
			}
			var cols = result.gaussJordan().columnVectors();
			return colVecsToMatrix(cols.splice(n, cols.length - 1));
		};

		/**
		 * Applies gaussian elimination to self returning the result
		 *
		 * @requires self is non-singular
		 * @return Matrix m the altered matrix
		 */
		self.gauss = function() {
			record('Matrix.gauss', self);
			var help = function(matrix, i) {
				if(i >= matrix.rowCount() || i >= matrix.colCount())
					return matrix;
				if(matrix.toArray()[i][i].isZero() && i < Math.min(matrix.colCount() - 1, matrix.rowCount() - 1)) {
					var found = false;
					for(var j = i + 1; j < matrix.rowCount(); j++) {
						if(!matrix.toArray()[j][i].isZero()) {
							matrix = matrix.swapRows(i, j);
							found = true;
							break;
						}
					}
					if(!found)
						throw "Matrix is singular";
				}
				for(var j = i + 1; j < matrix.rowCount(); j++) {
					var factor = matrix.toArray()[j][i].div(matrix.toArray()[i][i]).negate();
					matrix = matrix.addScaledRowToRow(i, factor, j);
				}
				return help(matrix, i + 1);
			}
			return help(self, 0);
		};

		/**
		 * Applies gaussian-jordan elimination to self returning the result
		 *
		 * @requires self is non-singular
		 * @return Matrix m the altered matrix
		 */
		self.gaussJordan = function() {
			record('Matrix.gaussJordan', self);
			var help = function(matrix, i) {
				if(i < 0)
					return matrix;
				for(var j = i - 1; j >= 0; j--) {
					if(!matrix.toArray()[i][i].isZero()) {
						var factor = matrix.toArray()[j][i].div(matrix.toArray()[i][i]).negate();
						matrix = matrix.addScaledRowToRow(i, factor, j);
					}
				}
				var makeOne = matrix.toArray()[i][i].invert();
				if(!makeOne.isZero())
					matrix = matrix.scaleRow(i, makeOne);
				return help(matrix, i - 1);
			}
			return help(self.gauss(), Math.min(self.colCount() - 1, self.rowCount() - 1));
		};

		/**
		 * Removes column i, row j from self and returns the result
		 *
		 * @requires i < self.rowCount() && i >= 0 && j < self.rowCount() && j >= 0
		 * @param int i the index of the row to remove
		 * @param int j the index of the column to remove
		 * @return Matrix m the altered matrix, null if not possible to remove row and column
		 */
		self.removeCross = function(i, j) {
			record('Matrix.removeCross', self);
			if(i >= self.rowCount() || i < 0 || j >= self.rowCount() || j < 0)
				return null;

			var rows = self.rowVectors();
			rows.splice(i, 1);
			return new Matrix(rows.map(function(r) {
				var result = r.toArray();
				result.splice(j, 1);
				return result;
			}));
		};

		/**
		 * Gets the diagonal of self as a vector.
		 * 
		 * @return Vector V where V is the diagonal of this
		 *         matrix.
		 */
		self.getDiagonalVector = function() {
			record('Matrix.getDiagonalVector', self);
			var result = [];
			var bound = Math.min(self.colCount(), self.rowCount());
			for(var i = 0; i < bound; i++) {
				result.push(data[i][i]);
			}
			return new Vector(result);
		};

		/**
		 * Calculates the determinant of self
		 *
		 * @requires self is NxN
		 * @return Fraction a the determinant of self
		 */
		self.determinant = function() {
			record('Matrix.determinant', self);
			if(self.rowCount() != self.colCount())
				throw "Determinant not defined for non NxN";
			if(self.rowCount() == 1)
				return data[0][0];
			//convert to RRF and get the diagonal.
			var id = startRecording(function(sig){
				return sig == 'Matrix.swapRows' || sig == 'Matrix.addScaledRowToRow'
						|| sig == 'Matrix.scaleRow'
			});
			var rrf = self.gaussJordan();
			var steps = stopRecording(id);
			var diagonal = rrf.getDiagonalVector();
			if(diagonal.hasZero()) {
				return ZERO;
			} else {
				var multiplyRow = diagonal.toArray().reduce(function(prev, x) {
					return prev.mult(x);
				});
				var factor = ONE;
				for(var i = 0; i < steps.length; i++) {
					if(steps[i].sig == 'Matrix.swapRows') {
						factor = factor.mult(NEGATIVE_ONE);
					} else if(steps[i].sig == 'Matrix.scaleRow') {
						factor = factor.mult(steps[i].args[1].invert());
					}
				}
				return factor.mult(multiplyRow);
			}
		};

		/**
		 * Multiplies this matrix by mat
		 *
		 * @requires self.colCount() == mat.rowCount()
		 * @return Matrix m such that m = self * mat, null if not possible 
		 */
		self.mult = function(mat) {
			record('Matrix.mult', self, [mat]);
			if(self.colCount() == mat.rowCount()){
				var myRows = self.rowVectors();
				var otherCols = mat.columnVectors();
				result = [];
				for(var i = 0; i < myRows.length; i++) {
					//push a new row in
					result.push([]);
					for(var j = 0; j < otherCols.length; j++) {
						result[i].push(myRows[i].dot(otherCols[j]))
					}
				}
				return new Matrix(result);
			}
			return null;
		};

		/**
		 * Adds mat to self
		 *
		 * @requires self.colCount() == mat.rowCount() &&
		 *           self.rowCount() == mat.rowCount()
		 * @return Matrix m such that m = self + mat, null if not possible 
		 */
		self.add = function(mat) {
			record('Matrix.add', self, [mat]);
			if(self.rowCount() != mat.rowCount() ||
				self.colCount() != mat.colCount())
				return null;
			var result = mat.toArray();
			for(var i = 0; i < result.length; i++) {
				for(var j = 0; j < result[i].length; j++) {
					result[i][j] = result[i][j].add(data[i][j]);
				}
			}
			return new Matrix(result);
		};

		/**
		 * Subtracts mat from self
		 *
		 * @requires self.colCount() == mat.rowCount() &&
		 *           self.rowCount() == mat.rowCount()
		 * @return Matrix m such that m = self - mat, null if not possible 
		 */
		self.sub = function(mat) {
			record('Matrix.sub', self, [mat]);
			if(self.rowCount() != mat.rowCount() ||
				self.colCount() != mat.colCount())
				return null;
			var result = mat.toArray();
			for(var i = 0; i < result.length; i++) {
				for(var j = 0; j < result[i].length; j++) {
					result[i][j] = data[i][j].sub(result[i][j]);
				}
			}
			return new Matrix(result);
		};

		/**
		 * Multiplies self by self n times
		 *
		 * @requires self is NxN
		 * @return Matrix m such that m = self^n
		 */
		self.pow = function(n) {
			record('Matrix.pow', self, [n]);
			if(self.rowCount() != self.colCount())
				throw "Self-multiplication not defined for non NxN";
			var help = function(matrix, k) {
				if(k < 0)
					return matrix;
				return help(matrix.mult(matrix), k - 1);
			};

			return help(self, n);
		};

		/**
		 * Returns self transposed
		 *
		 * @return Matrix m such that m = self^T
		 */
		self.transpose = function() {
			record('Matrix.transpose', self);
			var result = [];
			for(var i = 0; i < self.colCount(); i++) {
				result[i] = [];
				for(var j = 0; j < self.rowCount(); j++) {
					result[i].push(data[j][i]);
				}
			}
			return new Matrix(result);
		};

		/**
		 * Scales self by Fraction s
		 *
		 * @param Fraction s scalar fraction
		 * @return Matrix m such that m = self * s
		 */
		self.scale = function(s) {
			record('Matrix.scale', self, [s]);
			return new Matrix(self.rowVectors().map(function(x) {
				return x.scale(s).toArray();
			}));
		};

		/**
		 * Converts self to string representation
		 * 
		 * @return string s such that s represents self
		 */
		self.toString = function() {
			record('Matrix.toString', self);
			return self.rowVectors().map(function(x) {
				return x.toString();
			}).join("\n");
		};

	}


	return {'parseMatrix' : parseMatrix,
			'parseVector' : parseVector,
			'parseVectors' : parseVectors,
			'parseFraction' : parseFraction,
			'changeOfBase' : changeOfBaseMatrix,
			'rowVecsToMatrix' : rowVecsToMatrix,
			'colVecsToMatrix' : colVecsToMatrix,
			'orthogonalToSubspace' : orthogonalToSubspace,
			'projectSubspace' : projectSubspace,
			'identity' : identity,
			'bestFitLine' : bestFitLine,
			'debug' : {
				'startRecording' : startRecording,
				'stopRecording' : stopRecording
			},
			'zero' : zero};
})()

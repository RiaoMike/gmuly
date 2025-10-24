use ark_ff::Field;
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial, Polynomial,
};
use ark_std::vec::Vec;
use crate::SparseMatrix;
use std::collections::HashMap;

/// Multilinear Extension (MLE) representation
/// 
/// An MLE is a multilinear polynomial that extends a function from {0,1}^n to F^n.
/// The polynomial is represented in coefficient form with respect to the tensor product basis.
/// 
/// The basis is constructed as the tensor product (1, x_n) ⊗ ... ⊗ (1, x_2) ⊗ (1, x_1),
/// where x_1 is the first variable and x_n is the last variable.
/// 
/// For 2 variables: basis = (1, x_2) ⊗ (1, x_1) = (1, x_1, x_2, x_1*x_2)
/// For 3 variables: basis = (1, x_3) ⊗ (1, x_2) ⊗ (1, x_1) = (1, x_1, x_2, x_1*x_2, x_3, x_1*x_3, x_2*x_3, x_1*x_2*x_3)
/// 
/// The coefficient at index i corresponds to the basis element formed by the product
/// of variables selected by the binary representation of i:
/// - i = 0b00: coefficient for constant term 1
/// - i = 0b01: coefficient for x_1 (first variable, bit 0)
/// - i = 0b10: coefficient for x_2 (second variable, bit 1)
/// - i = 0b11: coefficient for x_1*x_2
/// 
/// For implementation convenience, the evaluate function accepts points in the natural order [x_1, x_2, ..., x_n]
/// 
#[derive(Debug, Clone, PartialEq)]
pub struct MLE<F: Field> {
    /// Number of variables
    pub num_vars: usize,
    /// Coefficients of the multilinear polynomial in tensor product basis
    /// Length must be 2^num_vars
    pub coefficients: Vec<F>,
}

impl<F: Field> MLE<F> {
    /// Create an MLE from coefficients
    /// 
    /// # Arguments
    /// * `coefficients` - Coefficients in tensor product basis
    /// 
    /// # Panics
    /// Panics if the length of coefficients is not a power of 2
    pub fn new(coefficients: Vec<F>) -> Self {
        let len = coefficients.len();
        assert!(len.is_power_of_two(), "Length must be a power of 2");
        
        let num_vars = len.trailing_zeros() as usize;
        
        Self {
            num_vars,
            coefficients,
        }
    }

    /// Create an MLE from evaluations at boolean hypercube points
    /// Converts evaluations to coefficient form
    pub fn from_evaluations(evaluations: &[F]) -> Self {
        // Pad to the nearest power of 2 if necessary
        let len = evaluations.len();
        let padded_len = len.next_power_of_two();
        
        let mut evals = evaluations.to_vec();
        evals.resize(padded_len, F::ZERO);
        
        // Convert evaluations to coefficients using Walsh-Hadamard transform
        let coefficients = Self::evaluations_to_coefficients(&evals);
        Self::new(coefficients)
    }
    
    /// Create an MLE from a vector (1D array)
    /// The vector is treated as evaluations at boolean points
    pub fn from_vector(vec: &[F]) -> Self {
        Self::from_evaluations(vec)
    }

    /// Create an MLE from a sparse matrix
    /// 
    /// Converts the sparse matrix to a dense format and then creates an MLE.
    /// 
    /// # Arguments
    /// * `matrix` - The sparse matrix to convert
    /// 
    /// # Returns
    /// An MLE with m1 + m2 variables where:
    /// - m1 = log2(rows.next_power_of_two())
    /// - m2 = log2(cols.next_power_of_two())
    pub fn from_sparse_matrix(matrix: &SparseMatrix<F>) -> Self {
        let rows = matrix.rows();
        let cols = matrix.cols();
        
        // Convert sparse matrix to dense format
        let mut dense_matrix = Vec::with_capacity(rows);
        for i in 0..rows {
            let mut row = Vec::with_capacity(cols);
            for j in 0..cols {
                row.push(matrix.get(i, j));
            }
            dense_matrix.push(row);
        }
        
        // Use from_matrix to create the MLE
        Self::from_matrix(&dense_matrix, rows, cols)
    }

    /// Create an MLE from a matrix (2D array)
    /// The matrix is treated as a 2D function where the first m1 variables
    /// encode the row index and the next m2 variables encode the column index.
    /// 
    /// # Arguments
    /// * `matrix` - 2D matrix to convert to MLE
    /// * `rows` - Number of rows
    /// * `cols` - Number of columns
    /// 
    /// # Returns
    /// An MLE with m1 + m2 variables where:
    /// - m1 = log2(rows.next_power_of_two())
    /// - m2 = log2(cols.next_power_of_two())
    pub fn from_matrix(matrix: &[Vec<F>], rows: usize, cols: usize) -> Self {
        assert_eq!(matrix.len(), rows, "Matrix rows mismatch");
        
        // Calculate padded dimensions
        let padded_rows = rows.next_power_of_two();
        let padded_cols = cols.next_power_of_two();
        let _m1 = padded_rows.trailing_zeros() as usize;
        let _m2 = padded_cols.trailing_zeros() as usize;
        
        // Create padded evaluations array
        let mut evaluations = vec![F::ZERO; padded_rows * padded_cols];
        
        // Copy matrix values into padded array in row-major order
        for (i, row) in matrix.iter().enumerate() {
            assert_eq!(row.len(), cols, "Matrix cols mismatch");
            for (j, value) in row.iter().enumerate() {
                evaluations[i * padded_cols + j] = *value;
            }
        }
        
        // Convert evaluations to coefficients
        let coefficients = Self::evaluations_to_coefficients(&evaluations);
        Self::new(coefficients)
    }


    /// Get the number of variables
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Get the coefficients
    pub fn coefficients(&self) -> &[F] {
        &self.coefficients
    }
    
    /// Get evaluations at boolean hypercube points
    pub fn evaluations(&self) -> Vec<F> {
        Self::coefficients_to_evaluations(&self.coefficients)
    }

    /// Partially evaluate the MLE by binding the first r variables to a point
    /// 
    /// # Arguments
    /// * `point` - Values for the first r variables [x_1, x_2, ..., x_r]
    ///   where r = point.len() <= num_vars
    /// 
    /// # Returns
    /// A new MLE with (num_vars - r) variables, representing f(point[0], point[1], ..., point[r-1], x_{r+1}, ..., x_n)
    /// 
    /// # Example
    /// For a 3-variable MLE f(x_1, x_2, x_3), calling partial_evaluate(&[a, b]) returns:
    /// g(x_3) = f(a, b, x_3)
    /// 
    /// # Implementation
    /// This function processes all coefficients in a single pass for optimal performance.
    /// For each coefficient in the original polynomial, we:
    /// 1. Extract which of the first r variables appear in the term
    /// 2. Compute the monomial value for those variables
    /// 3. Add the contribution to the corresponding term in the new polynomial
    /// 
    /// Complexity: O(2^n) where n is the number of variables
    pub fn partial_evaluate(&self, point: &[F]) -> Self {
        assert!(
            point.len() <= self.num_vars,
            "Point dimension must not exceed number of variables"
        );
        
        if point.is_empty() {
            return self.clone();
        }
        
        let r = point.len();
        let remaining_vars = self.num_vars - r;
        let new_size = 1 << remaining_vars;
        let mut new_coefficients = vec![F::ZERO; new_size];
        
        // Process all coefficients in a single pass
        for i in 0..self.coefficients.len() {
            // Extract which of the first r variables are in this term (lower r bits)
            let bound_vars_mask = i & ((1 << r) - 1);
            
            // Compute the monomial value for the bound variables
            // monomial = product of point[j] for all j where bit j of bound_vars_mask is 1
            let mut monomial = F::ONE;
            for j in 0..r {
                if (bound_vars_mask >> j) & 1 == 1 {
                    monomial *= point[j];
                }
            }
            
            // The new index is the upper (n-r) bits (variables that remain)
            let new_idx = i >> r;
            
            // Add the contribution: coefficient * monomial
            new_coefficients[new_idx] += self.coefficients[i] * monomial;
        }
        
        Self::new(new_coefficients)
    }

    /// Evaluate the MLE at a point
    /// 
    /// # Arguments
    /// * `point` - Point at which to evaluate the MLE (length must equal num_vars)
    ///   where point = [x_1, x_2, ..., x_n] for convenience, even though
    ///   the variable ordering in the tensor product is [x_n, ..., x_2, x_1]
    /// 
    /// # Returns
    /// The value of the MLE at the given point
    pub fn evaluate(&self, point: &[F]) -> F {
        assert_eq!(
            point.len(),
            self.num_vars,
            "Point dimension must match number of variables"
        );
        
        // Use partial_evaluate to bind all variables
        let result_mle = self.partial_evaluate(point);
        
        // The result should be a 0-variable MLE (constant polynomial)
        assert_eq!(result_mle.num_vars, 0, "Result should be a constant");
        assert_eq!(result_mle.coefficients.len(), 1, "Result should have exactly one coefficient");
        
        result_mle.coefficients[0]
    }

    /// Compute the Lagrange basis polynomial χ_b at point x
    /// χ_b(x) = Π_i (b_i * x_i + (1 - b_i) * (1 - x_i))
    #[allow(dead_code)]
    fn chi(b: usize, point: &[F]) -> F {
        let mut result = F::ONE;

        for (i, x_i) in point.iter().enumerate() {
            // b's binary representation: lowest bit corresponds to x_1, highest bit to x_n
            // point array order: [x_1, x_2, ..., x_n]
            let b_i = (b >> i) & 1;
            if b_i == 1 {
                result *= x_i;
            } else {
                result *= F::ONE - x_i;
            }
        }
        
        result
    }

    /// Bind a variable to a value and return a new MLE with one less variable
    /// 
    /// # Arguments
    /// * `var_index` - Index of the variable to bind (0-indexed, where 0 corresponds to x_1, 1 to x_2, etc.)
    /// * `value` - Value to bind the variable to
    /// 
    /// # Note
    /// The variable indexing follows the binary representation where bit i corresponds to x_{i+1}
    pub fn bind(&self, var_index: usize, value: F) -> Self {
        assert!(var_index < self.num_vars, "Variable index out of bounds");
        
        let new_num_vars = self.num_vars - 1;
        let new_size = 1 << new_num_vars;
        let mut new_coefficients = vec![F::ZERO; new_size];
        
        // When binding x_j to value v, the new polynomial is obtained by:
        // - For terms not containing x_j: coefficient remains the same
        // - For terms containing x_j: multiply coefficient by v and add to corresponding term without x_j
        
        for i in 0..self.coefficients.len() {
            let coeff = self.coefficients[i];
            let new_idx = (i & ((1 << var_index) - 1)) 
                        | ((i >> (var_index + 1)) << var_index);
            
            if (i >> var_index) & 1 == 0 {
                // This term doesn't contain x_{var_index}
                new_coefficients[new_idx] += coeff;
            } else {
                // This term contains x_{var_index}
                new_coefficients[new_idx] += coeff * value;
            }
        }
        
        Self::new(new_coefficients)
    }

    /// Get the multilinear polynomial representation as coefficients
    /// The coefficients follow the tensor product basis (1, x_n) ⊗ ... ⊗ (1, x_1),
    /// resulting in the order (1, x_1, x_2, x_1*x_2, x_3, x_1*x_3, x_2*x_3, x_1*x_2*x_3, ...)
    /// where coefficient at index i has variables determined by i's binary representation
    pub fn to_coefficients(&self) -> Vec<F> {
        self.coefficients.clone()
    }

    /// Create an MLE that represents the zero polynomial
    pub fn zero(num_vars: usize) -> Self {
        let size = 1 << num_vars;
        Self::new(vec![F::ZERO; size])
    }

    /// Create an MLE that represents the constant polynomial 1
    pub fn one(num_vars: usize) -> Self {
        let size = 1 << num_vars;
        let mut coefficients = vec![F::ZERO; size];
        coefficients[0] = F::ONE;
        Self::new(coefficients)
    }

    /// Add two MLEs (in-place modification)
    pub fn add(&mut self, other: &Self) {
        assert_eq!(self.num_vars, other.num_vars, "MLEs must have same number of variables");
        
        for i in 0..self.coefficients.len() {
            self.coefficients[i] += other.coefficients[i];
        }
    }

    /// Scale MLE by a constant (in-place modification)
    pub fn scale(&mut self, scalar: F) {
        for i in 0..self.coefficients.len() {
            self.coefficients[i] *= scalar;
        }
    }

    /// Multiply two MLEs to produce a general multivariate polynomial
    /// 
    /// Computes h(x) = f(x) * g(x) as a polynomial product in the field F.
    /// The result is NOT necessarily multilinear - it can have degree 2 terms.
    /// 
    /// This version uses a HashMap to efficiently merge terms with the same monomial,
    /// reducing redundant work when many products result in the same term.
    /// 
    /// # Arguments
    /// * `other` - The other MLE to multiply with
    /// 
    /// # Returns
    /// A SparsePolynomial representing the product f(x) * g(x)
    /// 
    /// # Example
    /// If f(x,y) = 1 + x and g(x,y) = 1 + y, then f*g = 1 + x + y + xy
    /// 
    /// # Complexity
    /// O(T1 * T2) where T1, T2 are the number of non-zero terms in each MLE.
    /// In the worst case (dense MLEs), this is O(2^n * 2^n) = O(4^n), but in practice
    /// MLEs are often sparse, making this much faster than the worst case.
    pub fn polynomial_multiply(&self, other: &Self) -> SparsePolynomial<F, SparseTerm> {
        assert_eq!(
            self.num_vars, other.num_vars,
            "MLEs must have same number of variables for multiplication"
        );
        
        // Convert both MLEs to SparsePolynomial
        let poly1 = self.to_sparse_polynomial();
        let poly2 = other.to_sparse_polynomial();
        
        // Use HashMap to merge terms with same monomials efficiently
        let mut term_map: HashMap<Vec<(usize, usize)>, F> = HashMap::new();
        
        for (coeff1, term1) in poly1.terms.iter() {
            for (coeff2, term2) in poly2.terms.iter() {
                // Multiply coefficients
                let new_coeff = *coeff1 * *coeff2;
                
                // Skip zero terms
                if new_coeff == F::ZERO {
                    continue;
                }
                
                // Multiply terms by adding exponents for each variable
                let mut var_powers: HashMap<usize, usize> = HashMap::new();
                
                // Add powers from term1
                for (var, pow) in term1.iter() {
                    *var_powers.entry(*var).or_insert(0) += pow;
                }
                
                // Add powers from term2
                for (var, pow) in term2.iter() {
                    *var_powers.entry(*var).or_insert(0) += pow;
                }
                
                // Convert to sorted vector for consistent key
                let mut var_vec: Vec<(usize, usize)> = var_powers.into_iter().collect();
                var_vec.sort_by_key(|(var, _)| *var);
                
                // Add to term_map, merging with existing terms
                *term_map.entry(var_vec).or_insert(F::ZERO) += new_coeff;
            }
        }
        
        // Convert HashMap to Vec, filtering out zero coefficients
        let result_terms: Vec<(F, SparseTerm)> = term_map
            .into_iter()
            .filter(|(_, coeff)| *coeff != F::ZERO)
            .map(|(var_vec, coeff)| (coeff, SparseTerm::new(var_vec)))
            .collect();
        
        SparsePolynomial::from_coefficients_vec(self.num_vars, result_terms)
    }

    /// Compute partial sum over the first r variables
    /// Returns sum_{x_1,...,x_r ∈ {0,1}} f(x_1,...,x_r,x_{r+1},...,x_n)
    /// 
    /// This is useful in the sumcheck protocol where we need to sum over
    /// a subset of variables while keeping others symbolic.
    /// 
    /// # Arguments
    /// * `num_sum_vars` - Number of variables to sum over (r), starting from x_1
    /// 
    /// # Returns
    /// A new MLE with (num_vars - num_sum_vars) variables representing the partial sum
    /// 
    /// # Example
    /// For a 3-variable MLE f(x_1, x_2, x_3), calling mle_partial_sum(2) returns:
    /// g(x_3) = f(0,0,x_3) + f(1,0,x_3) + f(0,1,x_3) + f(1,1,x_3)
    /// 
    /// # Implementation
    /// This function operates directly on coefficients without converting to evaluations.
    /// For a term c_S * ∏_{i∈S} x_i, when summing over variables {x_1,...,x_r}:
    /// sum_{x_1,...,x_r} c_S * ∏_{i∈S} x_i = c_S * 2^(r - |S∩{1,...,r}|) * ∏_{i∈S\{1,...,r}} x_i
    /// 
    /// Complexity: O(2^n) where n is the number of variables
    pub fn mle_partial_sum(&self, num_sum_vars: usize) -> Self {
        assert!(num_sum_vars <= self.num_vars, "Number of sum variables exceeds total variables");
        
        if num_sum_vars == 0 {
            return self.clone();
        }
        
        if num_sum_vars == self.num_vars {
            // Sum over all variables gives a constant polynomial
            // For each coefficient c_S, its contribution is c_S * 2^(n - |S|)
            let mut sum = F::ZERO;
            for i in 0..self.coefficients.len() {
                let num_vars_in_term = i.count_ones() as usize;
                let factor_exp = num_sum_vars - num_vars_in_term;
                let factor = F::from((1u64 << factor_exp) as u64);
                sum += self.coefficients[i] * factor;
            }
            return Self::new(vec![sum]);
        }
        
        let remaining_vars = self.num_vars - num_sum_vars;
        let new_size = 1 << remaining_vars;
        let mut new_coefficients = vec![F::ZERO; new_size];
        
        // For each coefficient in the original polynomial
        for i in 0..self.coefficients.len() {
            // Extract which of the first num_sum_vars variables are in this term (lower bits)
            let sum_vars_mask = i & ((1 << num_sum_vars) - 1);
            let num_sum_vars_in_term = sum_vars_mask.count_ones() as usize;
            
            // The new index is the upper bits (variables that are not summed over)
            let new_idx = i >> num_sum_vars;
            
            // The contribution factor is 2^(num_sum_vars - num_sum_vars_in_term)
            // This accounts for the fact that summing over a variable not in the term doubles it,
            // while summing over a variable in the term evaluates it at 0 and 1
            let factor_exp = num_sum_vars - num_sum_vars_in_term;
            let factor = F::from((1u64 << factor_exp) as u64);
            
            new_coefficients[new_idx] += self.coefficients[i] * factor;
        }
        
        Self::new(new_coefficients)
    }
    
    /// Get the split count
    pub fn get_split_count(&self) -> usize {
        self.split
    }
    
    /// Convert MLE to a SparseMatrix
    /// 
    /// Directly maps MLE coefficients to a matrix representation.
    /// The dimensions (rows, cols) must be provided and their product
    /// must equal 2^num_vars. Both rows and cols must be powers of 2.
    /// 
    /// # Arguments
    /// * `rows` - Number of rows in the resulting matrix
    /// * `cols` - Number of columns in the resulting matrix
    /// 
    /// # Returns
    /// A SparseMatrix with the specified dimensions
    /// 
    /// # Panics
    /// Panics if rows * cols != 2^num_vars or if dimensions are not powers of 2
    pub fn to_sparse_matrix(&self) -> SparseMatrix<F> {
        let cols = 1 << self.split;
        let rows = 1 << (self.num_vars - self.split);
        
        // Build dense matrix representation directly from coefficients
        // Directly convert coefficients to sparse matrix without intermediate dense representation
        let mut values = Vec::new();
        let mut col_indices = Vec::new();
        let mut row_ptr = Vec::with_capacity(rows + 1);
        row_ptr.push(0);
        
        for row in 0..rows {
            for col in 0..cols {
                let idx = row * cols + col;
                let val = self.coefficients[idx];
                if !val.is_zero() {
                    values.push(val);
                    col_indices.push(col);
                }
            }
            row_ptr.push(values.len());
        }
        
        SparseMatrix::from_csr(cols, values, col_indices, row_ptr).unwrap()
    }
    
    /// Convert MLE to SparsePolynomial
    /// 
    /// Converts the MLE representation to a sparse multivariate polynomial.
    /// Each coefficient in the MLE corresponds to a term in the polynomial.
    /// 
    /// # Returns
    /// A SparsePolynomial with the same number of variables
    /// 
    /// # Example
    /// ```rust,ignore
    /// // MLE with 2 variables: coefficients [c0, c1, c2, c3]
    /// // represents: c0 + c1*x_0 + c2*x_1 + c3*x_0*x_1
    /// let mle = MLE::new(vec![1, 2, 3, 4]);
    /// let sparse_poly = mle.to_sparse_polynomial();
    /// ```
    pub fn to_sparse_polynomial(&self) -> SparsePolynomial<F, SparseTerm> {
        let mut terms = Vec::new();
        
        // For each coefficient in the MLE
        for (index, coeff) in self.coefficients.iter().enumerate() {
            if coeff.is_zero() {
                continue;
            }
            
            // Convert index to term
            // The binary representation of index determines which variables appear
            let mut term_vec = Vec::new();
            for var_idx in 0..self.num_vars {
                if (index >> var_idx) & 1 == 1 {
                    // Variable x_{var_idx} appears in this term with power 1
                    term_vec.push((var_idx, 1));
                }
            }
            
            let term = SparseTerm::new(term_vec);
            terms.push((*coeff, term));
        }
        
        SparsePolynomial::from_coefficients_vec(self.num_vars, terms)
    }
    
    /// Create MLE from SparsePolynomial
    /// 
    /// Converts a multivariate sparse polynomial to MLE representation.
    /// The polynomial must be multilinear (all variables have degree at most 1).
    /// 
    /// # Arguments
    /// * `poly` - The sparse polynomial to convert
    /// 
    /// # Returns
    /// An MLE with the same number of variables
    /// 
    /// # Panics
    /// Panics if the polynomial is not multilinear (any variable has degree > 1)
    /// 
    /// # Example
    /// ```rust,ignore
    /// // Create a polynomial: 1 + 2*x_0 + 3*x_1 + 4*x_0*x_1
    /// let poly = SparsePolynomial::from_coefficients_vec(...);
    /// let mle = MLE::from_sparse_polynomial(&poly);
    /// ```
    pub fn from_sparse_polynomial(poly: &SparsePolynomial<F, SparseTerm>) -> Self {
        let num_vars = poly.num_vars;
        let size = 1 << num_vars;
        let mut coefficients = vec![F::ZERO; size];
        
        // For each term in the polynomial
        for (coeff, term) in poly.terms.iter() {
            // Convert term to index in MLE coefficient array
            let mut index = 0;
            
            for (var, power) in term.iter() {
                // Check that the polynomial is multilinear
                assert!(
                    *power <= 1,
                    "Polynomial must be multilinear (all powers <= 1), but found variable {} with power {}",
                    var, power
                );
                
                if *power == 1 {
                    index |= 1 << var;
                }
            }
            
            coefficients[index] += *coeff;
        }
        
        Self::new(coefficients)
    }

    /// Convert evaluations at boolean hypercube points to coefficient representation
    /// Uses the inverse Walsh-Hadamard transform
    fn evaluations_to_coefficients(evaluations: &[F]) -> Vec<F> {
        let n = evaluations.len();
        assert!(n.is_power_of_two(), "Length must be a power of 2");
        
        let mut coefficients = evaluations.to_vec();
        let num_vars = n.trailing_zeros() as usize;
        
        // Apply inverse Walsh-Hadamard transform
        // For each variable, we perform a butterfly operation
        for j in 0..num_vars {
            let stride = 1 << j;
            let pairs = n / (2 * stride);
            
            for k in 0..pairs {
                for i in 0..stride {
                    let idx1 = 2 * k * stride + i;
                    let idx2 = (2 * k + 1) * stride + i;
                    
                    let a = coefficients[idx1];
                    let b = coefficients[idx2];
                    
                    coefficients[idx1] = a;
                    coefficients[idx2] = b - a;
                }
            }
        }
        
        coefficients
    }
    
    /// Convert coefficient representation to evaluations at boolean hypercube points
    /// Uses the Walsh-Hadamard transform
    fn coefficients_to_evaluations(coefficients: &[F]) -> Vec<F> {
        let n = coefficients.len();
        assert!(n.is_power_of_two(), "Length must be a power of 2");
        
        let mut evaluations = coefficients.to_vec();
        let num_vars = n.trailing_zeros() as usize;
        
        // Apply Walsh-Hadamard transform
        // For each variable, we perform a butterfly operation in reverse order
        for j in (0..num_vars).rev() {
            let stride = 1 << j;
            let pairs = n / (2 * stride);
            
            for k in 0..pairs {
                for i in 0..stride {
                    let idx1 = 2 * k * stride + i;
                    let idx2 = (2 * k + 1) * stride + i;
                    
                    let a = evaluations[idx1];
                    let b = evaluations[idx2];
                    
                    evaluations[idx1] = a;
                    evaluations[idx2] = a + b;
                }
            }
        }
        
        evaluations
    }
}

/// Helper function to convert index to binary vector
/// binary = 1001 -> [true, false, false, true]
/// for variables [x_1, x_2, x_3, x_4], this corresponds to x_1*x_4
/// this is the reverse of the natural order of variables
#[allow(dead_code)]
fn index_to_binary_vec(index: usize, num_vars: usize) -> Vec<bool> {
    let mut result = vec![false; num_vars];
    for i in 0..num_vars {
        result[i] = ((index >> i) & 1) == 1;
    }
    result
}

/// Helper function to convert binary vector to index
#[allow(dead_code)]
fn binary_vec_to_index(binary: &[bool]) -> usize {
    let mut result = 0;
    for (i, &b) in binary.iter().enumerate() {
        if b {
            result |= 1 << i;
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_mle_from_vector() {
        // Test with evaluation vector [3, 5, 7, 9]
        let vec = vec![
            Fr::from(3u64),
            Fr::from(5u64),
            Fr::from(7u64),
            Fr::from(9u64),
        ];
        
        let mle = MLE::from_vector(&vec);
        assert_eq!(mle.num_vars(), 2);
        assert_eq!(mle.coefficients().len(), 4);
        
        // Check that coefficients are correct [3, 2, 4, 0]
        assert_eq!(mle.coefficients()[0], Fr::from(3u64));
        assert_eq!(mle.coefficients()[1], Fr::from(2u64));
        assert_eq!(mle.coefficients()[2], Fr::from(4u64));
        assert_eq!(mle.coefficients()[3], Fr::from(0u64));
    }

    #[test]
    fn test_mle_from_matrix() {
        // Test with 2x2 matrix [[1, 2], [3, 4]]
        let matrix = vec![
            vec![Fr::from(1u64), Fr::from(2u64)],
            vec![Fr::from(3u64), Fr::from(4u64)],
        ];
        
        let mle = MLE::from_matrix(&matrix, 2, 2);
        assert_eq!(mle.num_vars(), 2);
        
        // The matrix is flattened to evaluations [1, 2, 3, 4]
        // These should convert to coefficients [1, 1, 2, 0]
        assert_eq!(mle.coefficients()[0], Fr::from(1u64));
        assert_eq!(mle.coefficients()[1], Fr::from(1u64));
        assert_eq!(mle.coefficients()[2], Fr::from(2u64));
        assert_eq!(mle.coefficients()[3], Fr::from(0u64));
    }

    #[test]
    fn test_mle_evaluate() {
        // Create MLE with coefficients for polynomial:
        // f(x_1, x_2) = 3 + 2*x_1 + 4*x_2 + 0*x_1*x_2
        // This gives evaluations: f(0,0)=3, f(1,0)=5, f(0,1)=7, f(1,1)=9
        // Note: basis ordering is [1, x_1, x_2, x_1*x_2]
        let coefficients = vec![
            Fr::from(3u64),  // constant term
            Fr::from(2u64),  // x_1 coefficient
            Fr::from(4u64),  // x_2 coefficient
            Fr::from(0u64),  // x_1*x_2 coefficient
        ];
        
        let mle = MLE::new(coefficients);
        
        // Test evaluation at boolean points
        assert_eq!(mle.evaluate(&[Fr::from(0u64), Fr::from(0u64)]), Fr::from(3u64)); // f(0,0) = 3
        assert_eq!(mle.evaluate(&[Fr::from(1u64), Fr::from(0u64)]), Fr::from(5u64)); // f(1,0) = 3 + 2 = 5
        assert_eq!(mle.evaluate(&[Fr::from(0u64), Fr::from(1u64)]), Fr::from(7u64)); // f(0,1) = 3 + 4 = 7
        assert_eq!(mle.evaluate(&[Fr::from(1u64), Fr::from(1u64)]), Fr::from(9u64)); // f(1,1) = 3 + 2 + 4 = 9
        
        // Test evaluation at intermediate point (0.5, 0.5)
        let half = Fr::from(2u64).inverse().unwrap();
        let result = mle.evaluate(&[half, half]);
        // Expected: 3 + 2*0.5 + 4*0.5 + 0*0.5*0.5 = 3 + 1 + 2 = 6
        assert_eq!(result, Fr::from(6u64));
    }

    #[test]
    fn test_mle_bind() {
        // Create MLE from evaluations [1, 2, 3, 4, 5, 6, 7, 8]
        // This is a 3-variable MLE
        let evaluations: Vec<Fr> = (1..=8).map(|i| Fr::from(i as u64)).collect();
        let mle = MLE::from_evaluations(&evaluations);
        
        // Bind x_1 (index 0) to 0
        let mle_bound = mle.bind(0, Fr::from(0u64));
        assert_eq!(mle_bound.num_vars(), 2);
        // Should get evaluations at (0, x_2, x_3)
        let bound_evals = mle_bound.evaluations();
        assert_eq!(bound_evals, vec![
            Fr::from(1u64), Fr::from(3u64),
            Fr::from(5u64), Fr::from(7u64),
        ]);
        
        // Bind x_1 to 1
        let mle_bound = mle.bind(0, Fr::from(1u64));
        let bound_evals = mle_bound.evaluations();
        assert_eq!(bound_evals, vec![
            Fr::from(2u64), Fr::from(4u64),
            Fr::from(6u64), Fr::from(8u64),
        ]);
    }

    #[test]
    fn test_mle_add() {
        // Create MLEs from evaluations
        let evals1 = vec![
            Fr::from(1u64), Fr::from(2u64),
            Fr::from(3u64), Fr::from(4u64),
        ];
        let mle1 = MLE::from_evaluations(&evals1);
        
        let evals2 = vec![
            Fr::from(5u64), Fr::from(6u64),
            Fr::from(7u64), Fr::from(8u64),
        ];
        let mle2 = MLE::from_evaluations(&evals2);
        
        // Test addition (in-place)
        let mut sum = mle1.clone();
        sum.add(&mle2);
        let sum_evals = sum.evaluations();
        assert_eq!(sum_evals, vec![
            Fr::from(6u64), Fr::from(8u64),
            Fr::from(10u64), Fr::from(12u64),
        ]);
    }

    #[test]
    fn test_coefficient_ordering() {
        // For 2 variables, coefficients are ordered as:
        // (c_0, c_1, c_2, c_3) corresponding to basis (1, x_1, x_2, x_1*x_2)
        // where bit 0 corresponds to x_1 and bit 1 corresponds to x_2
        
        // Create MLE where f(x_1, x_2) = 1 + 2*x_1 + 3*x_2 + 4*x_1*x_2
        let coefficients = vec![
            Fr::from(1u64),  // coefficient for 1
            Fr::from(2u64),  // coefficient for x_1
            Fr::from(3u64),  // coefficient for x_2
            Fr::from(4u64),  // coefficient for x_1*x_2
        ];
        
        let mle = MLE::new(coefficients);
        
        // Check evaluations at boolean points
        // Now evaluate takes point = [x_1, x_2]
        assert_eq!(mle.evaluate(&[Fr::from(0u64), Fr::from(0u64)]), Fr::from(1u64)); // f(0,0) = 1
        assert_eq!(mle.evaluate(&[Fr::from(1u64), Fr::from(0u64)]), Fr::from(3u64)); // f(1,0) = 1 + 2 = 3
        assert_eq!(mle.evaluate(&[Fr::from(0u64), Fr::from(1u64)]), Fr::from(4u64)); // f(0,1) = 1 + 3 = 4
        assert_eq!(mle.evaluate(&[Fr::from(1u64), Fr::from(1u64)]), Fr::from(10u64)); // f(1,1) = 1 + 2 + 3 + 4 = 10
    }
    
    #[test]
    fn test_evaluations_coefficients_conversion() {
        // Test conversion between evaluations and coefficients
        let evaluations = vec![
            Fr::from(3u64),  // f(0,0)
            Fr::from(5u64),  // f(1,0)
            Fr::from(7u64),  // f(0,1)
            Fr::from(9u64),  // f(1,1)
        ];
        
        // Convert evaluations to coefficients
        let coefficients = MLE::<Fr>::evaluations_to_coefficients(&evaluations);
        
        // Expected coefficients: [3, 2, 4, 0]
        assert_eq!(coefficients[0], Fr::from(3u64));
        assert_eq!(coefficients[1], Fr::from(2u64));
        assert_eq!(coefficients[2], Fr::from(4u64));
        assert_eq!(coefficients[3], Fr::from(0u64));
        
        // Convert back to evaluations
        let evals_back = MLE::<Fr>::coefficients_to_evaluations(&coefficients);
        assert_eq!(evals_back, evaluations);
    }
    
    #[test]
    fn test_mle_from_matrix_with_padding() {
        // Test with 3x5 matrix to verify m1+m2 variables behavior
        let matrix = vec![
            vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64), Fr::from(5u64)],
            vec![Fr::from(6u64), Fr::from(7u64), Fr::from(8u64), Fr::from(9u64), Fr::from(10u64)],
            vec![Fr::from(11u64), Fr::from(12u64), Fr::from(13u64), Fr::from(14u64), Fr::from(15u64)],
        ];
        
        let mle = MLE::from_matrix(&matrix, 3, 5);
        
        // 3 rows -> next_power_of_two = 4 = 2^2, so m1 = 2
        // 5 cols -> next_power_of_two = 8 = 2^3, so m2 = 3
        // Total variables = m1 + m2 = 2 + 3 = 5
        assert_eq!(mle.num_vars(), 5);
        
        // The MLE should represent a 4x8 = 32 element padded matrix
        let evaluations = mle.evaluations();
        assert_eq!(evaluations.len(), 32);
        
        // Verify that the first 15 values match the original matrix values
        // padded matrix layout (row-major, 8 cols per row):
        // Row 0: 1, 2, 3, 4, 5, 0, 0, 0
        // Row 1: 6, 7, 8, 9, 10, 0, 0, 0
        // Row 2: 11, 12, 13, 14, 15, 0, 0, 0
        // Row 3: 0, 0, 0, 0, 0, 0, 0, 0
        assert_eq!(evaluations[0], Fr::from(1u64));   // (0,0)
        assert_eq!(evaluations[1], Fr::from(2u64));   // (0,1)
        assert_eq!(evaluations[2], Fr::from(3u64));   // (0,2)
        assert_eq!(evaluations[3], Fr::from(4u64));   // (0,3)
        assert_eq!(evaluations[4], Fr::from(5u64));   // (0,4)
        assert_eq!(evaluations[5], Fr::from(0u64));   // (0,5) - padding
        assert_eq!(evaluations[8], Fr::from(6u64));   // (1,0)
        assert_eq!(evaluations[9], Fr::from(7u64));   // (1,1)
        assert_eq!(evaluations[16], Fr::from(11u64)); // (2,0)
        assert_eq!(evaluations[24], Fr::from(0u64));  // (3,0) - padding
    }
    
    #[test]
    fn test_bind_function() {
        // Test bind function with a simple 2-variable polynomial
        // The basis ordering follows (1, x_2) ⊗ (1, x_1) = (1, x_1, x_2, x_1*x_2)
        // 
        // Binary index mapping:
        // - index 0 (0b00): coefficient for 1
        // - index 1 (0b01): coefficient for x_1 (bit 0 is set)
        // - index 2 (0b10): coefficient for x_2 (bit 1 is set)
        // - index 3 (0b11): coefficient for x_1*x_2
        //
        // So for f(x_1, x_2) = 1 + 2*x_1 + 3*x_2 + 4*x_1*x_2:
        let coefficients = vec![
            Fr::from(1u64),  // index 0: coefficient for 1
            Fr::from(2u64),  // index 1: coefficient for x_1
            Fr::from(3u64),  // index 2: coefficient for x_2
            Fr::from(4u64),  // index 3: coefficient for x_1*x_2
        ];
        
        let mle = MLE::new(coefficients);
        
        // Test binding x_1 to value 5
        // Expected result: f(5, x_2) = 1 + 2*5 + 3*x_2 + 4*5*x_2 = 11 + 23*x_2
        let bound_x1 = mle.bind(0, Fr::from(5u64)); // var_index=0 for x_1
        assert_eq!(bound_x1.num_vars(), 1);
        assert_eq!(bound_x1.coefficients[0], Fr::from(11u64)); // constant term
        assert_eq!(bound_x1.coefficients[1], Fr::from(23u64)); // coefficient of x_2
        
        // Test binding x_2 to value 3
        // Expected result: f(x_1, 3) = 1 + 2*x_1 + 3*3 + 4*x_1*3 = 10 + 14*x_1
        let bound_x2 = mle.bind(1, Fr::from(3u64)); // var_index=1 for x_2
        assert_eq!(bound_x2.num_vars(), 1);
        assert_eq!(bound_x2.coefficients[0], Fr::from(10u64)); // constant term
        assert_eq!(bound_x2.coefficients[1], Fr::from(14u64)); // coefficient of x_1
    }
    
    #[test]
    fn test_evaluate_correctness() {
        // Test that evaluation produces correct results
        // MLE from evaluations represents polynomial:
        // f(x_1, x_2) = 3 + 2*x_1 + 4*x_2 + 0*x_1*x_2
        // Evaluation ordering: evaluations[i] where i's binary representation 
        // determines the boolean values with bit j corresponding to x_{j+1}
        // So: [f(0,0), f(1,0), f(0,1), f(1,1)]
        let evaluations = vec![
            Fr::from(3u64),  // index 0 = 0b00 => f(0,0) = 3
            Fr::from(5u64),  // index 1 = 0b01 => f(1,0) = 3 + 2 = 5
            Fr::from(7u64),  // index 2 = 0b10 => f(0,1) = 3 + 4 = 7
            Fr::from(9u64),  // index 3 = 0b11 => f(1,1) = 3 + 2 + 4 = 9
        ];
        
        let mle = MLE::from_evaluations(&evaluations);
        
        // Test at boolean points
        assert_eq!(mle.evaluate(&[Fr::from(0u64), Fr::from(0u64)]), Fr::from(3u64)); // f(0,0) = 3
        assert_eq!(mle.evaluate(&[Fr::from(1u64), Fr::from(0u64)]), Fr::from(5u64)); // f(1,0) = 3 + 2 = 5
        assert_eq!(mle.evaluate(&[Fr::from(0u64), Fr::from(1u64)]), Fr::from(7u64)); // f(0,1) = 3 + 4 = 7
        assert_eq!(mle.evaluate(&[Fr::from(1u64), Fr::from(1u64)]), Fr::from(9u64)); // f(1,1) = 3 + 2 + 4 = 9
        
        // Test at (0.5, 0.5)
        // f(0.5, 0.5) = 3 + 2*0.5 + 4*0.5 + 0*0.5*0.5 = 3 + 1 + 2 = 6
        let half = Fr::from(2u64).inverse().unwrap();
        assert_eq!(mle.evaluate(&[half, half]), Fr::from(6u64));
        
        // Test with larger MLE - 3 variables
        // Evaluation ordering: evaluations[i] where i = x_1 + 2*x_2 + 4*x_3
        // So: [f(0,0,0), f(1,0,0), f(0,1,0), f(1,1,0), f(0,0,1), f(1,0,1), f(0,1,1), f(1,1,1)]
        let large_evals = vec![
            Fr::from(1u64), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64),
            Fr::from(5u64), Fr::from(6u64), Fr::from(7u64), Fr::from(8u64),
        ];
        let large_mle = MLE::from_evaluations(&large_evals);
        
        // Test at boolean points
        assert_eq!(large_mle.evaluate(&[Fr::from(0u64), Fr::from(0u64), Fr::from(0u64)]), Fr::from(1u64)); // f(0,0,0) = 1
        assert_eq!(large_mle.evaluate(&[Fr::from(1u64), Fr::from(0u64), Fr::from(0u64)]), Fr::from(2u64)); // f(1,0,0) = 2
        assert_eq!(large_mle.evaluate(&[Fr::from(0u64), Fr::from(1u64), Fr::from(0u64)]), Fr::from(3u64)); // f(0,1,0) = 3
        assert_eq!(large_mle.evaluate(&[Fr::from(0u64), Fr::from(0u64), Fr::from(1u64)]), Fr::from(5u64)); // f(0,0,1) = 5
        assert_eq!(large_mle.evaluate(&[Fr::from(1u64), Fr::from(1u64), Fr::from(1u64)]), Fr::from(8u64)); // f(1,1,1) = 8
    }

    #[test]
    fn test_partialsum_basic() {
        // Test with a simple 2-variable MLE
        // Evaluations: [1, 2, 3, 4] for f(x_1, x_2) at points (0,0), (1,0), (0,1), (1,1)
        let evaluations = vec![
            Fr::from(1u64), // f(0,0) = 1
            Fr::from(2u64), // f(1,0) = 2
            Fr::from(3u64), // f(0,1) = 3
            Fr::from(4u64), // f(1,1) = 4
        ];
        let mle = MLE::from_evaluations(&evaluations);
        
        // Sum over x_1: g(x_2) = f(0,x_2) + f(1,x_2)
        // g(0) = f(0,0) + f(1,0) = 1 + 2 = 3
        // g(1) = f(0,1) + f(1,1) = 3 + 4 = 7
        let sum_x1 = mle.mle_partial_sum(1);
        assert_eq!(sum_x1.num_vars(), 1);
        let sum_evals = sum_x1.evaluations();
        assert_eq!(sum_evals[0], Fr::from(3u64)); // g(0) = 3
        assert_eq!(sum_evals[1], Fr::from(7u64)); // g(1) = 7
        
        // Sum over both variables: constant = f(0,0) + f(1,0) + f(0,1) + f(1,1)
        // = 1 + 2 + 3 + 4 = 10
        let sum_all = mle.mle_partial_sum(2);
        assert_eq!(sum_all.num_vars(), 0);
        assert_eq!(sum_all.coefficients()[0], Fr::from(10u64));
    }

    #[test]
    fn test_partialsum_three_vars() {
        // Test with a 3-variable MLE
        // Evaluations at points: (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1,1,1)
        let evaluations = vec![
            Fr::from(1u64),  // f(0,0,0)
            Fr::from(2u64),  // f(1,0,0)
            Fr::from(3u64),  // f(0,1,0)
            Fr::from(4u64),  // f(1,1,0)
            Fr::from(5u64),  // f(0,0,1)
            Fr::from(6u64),  // f(1,0,1)
            Fr::from(7u64),  // f(0,1,1)
            Fr::from(8u64),  // f(1,1,1)
        ];
        let mle = MLE::from_evaluations(&evaluations);
        
        // Sum over x_1: g(x_2, x_3) = f(0,x_2,x_3) + f(1,x_2,x_3)
        // g(0,0) = f(0,0,0) + f(1,0,0) = 1 + 2 = 3
        // g(1,0) = f(0,1,0) + f(1,1,0) = 3 + 4 = 7
        // g(0,1) = f(0,0,1) + f(1,0,1) = 5 + 6 = 11
        // g(1,1) = f(0,1,1) + f(1,1,1) = 7 + 8 = 15
        let sum_x1 = mle.mle_partial_sum(1);
        assert_eq!(sum_x1.num_vars(), 2);
        let sum_evals = sum_x1.evaluations();
        assert_eq!(sum_evals[0], Fr::from(3u64));  // g(0,0)
        assert_eq!(sum_evals[1], Fr::from(7u64));  // g(1,0)
        assert_eq!(sum_evals[2], Fr::from(11u64)); // g(0,1)
        assert_eq!(sum_evals[3], Fr::from(15u64)); // g(1,1)
        
        // Sum over x_1 and x_2: h(x_3) = sum_{x_1,x_2} f(x_1,x_2,x_3)
        // h(0) = f(0,0,0) + f(1,0,0) + f(0,1,0) + f(1,1,0) = 1 + 2 + 3 + 4 = 10
        // h(1) = f(0,0,1) + f(1,0,1) + f(0,1,1) + f(1,1,1) = 5 + 6 + 7 + 8 = 26
        let sum_x1_x2 = mle.mle_partial_sum(2);
        assert_eq!(sum_x1_x2.num_vars(), 1);
        let sum_evals = sum_x1_x2.evaluations();
        assert_eq!(sum_evals[0], Fr::from(10u64)); // h(0)
        assert_eq!(sum_evals[1], Fr::from(26u64)); // h(1)
        
        // Sum over all variables
        let sum_all = mle.mle_partial_sum(3);
        assert_eq!(sum_all.num_vars(), 0);
        assert_eq!(sum_all.coefficients()[0], Fr::from(36u64)); // 1+2+3+4+5+6+7+8
    }

    #[test]
    fn test_partialsum_zero_vars() {
        // Test that mle_partial_sum(0) returns a clone of the original MLE
        let evaluations = vec![
            Fr::from(1u64), Fr::from(2u64),
            Fr::from(3u64), Fr::from(4u64),
        ];
        let mle = MLE::from_evaluations(&evaluations);
        
        let sum_zero = mle.mle_partial_sum(0);
        assert_eq!(sum_zero.num_vars(), mle.num_vars());
        assert_eq!(sum_zero.coefficients(), mle.coefficients());
    }

    #[test]
    fn test_partialsum_for_sumcheck() {
        // Test a typical sumcheck protocol scenario
        // For a polynomial f(x_1, x_2, x_3), in round 1 of sumcheck,
        // we compute g_1(x_1) = sum_{x_2,x_3 ∈ {0,1}} f(x_1, x_2, x_3)
        
        let evaluations = vec![
            Fr::from(1u64),  // f(0,0,0)
            Fr::from(2u64),  // f(1,0,0)
            Fr::from(3u64),  // f(0,1,0)
            Fr::from(4u64),  // f(1,1,0)
            Fr::from(5u64),  // f(0,0,1)
            Fr::from(6u64),  // f(1,0,1)
            Fr::from(7u64),  // f(0,1,1)
            Fr::from(8u64),  // f(1,1,1)
        ];
        let mle = MLE::from_evaluations(&evaluations);
        
        // To compute g_1(x_1), we first need to realize that mle_partial_sum
        // currently sums over the FIRST r variables. For sumcheck, we want
        // to sum over the LAST (n-1) variables to get g_1(x_1).
        // 
        // However, we can use the current implementation by understanding the
        // index structure. Let's verify the behavior.
        
        // Sum over the last 2 variables: This doesn't match the typical sumcheck,
        // so let's think of it differently.
        // 
        // Actually, for this test, let's just verify that the function works correctly.
        // In a real sumcheck implementation, we might need additional helper functions.
        
        // Sum over x_1, x_2: we get h(x_3)
        let h = mle.mle_partial_sum(2);
        assert_eq!(h.num_vars(), 1);
        
        // Verify h(0) and h(1)
        assert_eq!(h.evaluate(&[Fr::from(0u64)]), Fr::from(10u64)); // sum of first 4 values
        assert_eq!(h.evaluate(&[Fr::from(1u64)]), Fr::from(26u64)); // sum of last 4 values
    }

    #[test]
    fn test_sparse_matrix_conversion() {
        // Create a dense matrix (must be power-of-2 dimensions)
        let dense = vec![
            vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)],
            vec![Fr::from(5u64), Fr::from(6u64), Fr::from(7u64), Fr::from(8u64)],
        ];
        
        // Convert to SparseMatrix
        let sparse = SparseMatrix::from_dense(&dense);
        
        // Convert to MLE (treating matrix values as evaluations)
        let mle = MLE::from_sparse_matrix(&sparse);
        
        // The MLE should have 2x4 = 8 elements
        assert_eq!(mle.num_vars(), 3); // log2(8) = 3
        assert_eq!(mle.get_split_count(), 2); // log2(4) = 2 (columns)
        
        // Verify the MLE encodes the matrix correctly by checking evaluations
        // M[row, col] -> mle.evaluate([col_bits..., row_bits...])
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(0u64), Fr::from(0u64)]), Fr::from(1u64)); // M[0,0]
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(0u64), Fr::from(0u64)]), Fr::from(2u64)); // M[0,1]
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(1u64), Fr::from(0u64)]), Fr::from(3u64)); // M[0,2]
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(1u64), Fr::from(0u64)]), Fr::from(4u64)); // M[0,3]
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(0u64), Fr::from(1u64)]), Fr::from(5u64)); // M[1,0]
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(0u64), Fr::from(1u64)]), Fr::from(6u64)); // M[1,1]
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(1u64), Fr::from(1u64)]), Fr::from(7u64)); // M[1,2]
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(1u64), Fr::from(1u64)]), Fr::from(8u64)); // M[1,3]
    }
    
    #[test]
    fn test_square_sparse_matrix() {
        // Create a 4x4 matrix
        let dense = vec![
            vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)],
            vec![Fr::from(5u64), Fr::from(6u64), Fr::from(7u64), Fr::from(8u64)],
            vec![Fr::from(9u64), Fr::from(10u64), Fr::from(11u64), Fr::from(12u64)],
            vec![Fr::from(13u64), Fr::from(14u64), Fr::from(15u64), Fr::from(16u64)],
        ];
        
        let sparse = SparseMatrix::from_dense(&dense);
        let mle = MLE::from_sparse_matrix(&sparse);
        
        // Should have 4x4 = 16 elements, num_vars = 4 (log2(16))
        assert_eq!(mle.num_vars(), 4);
        assert_eq!(mle.get_split_count(), 2); // log2(4) = 2 for 4 columns
        
        // Verify the MLE encodes the matrix correctly
        // For 4x4 matrix: [col_bit0, col_bit1, row_bit0, row_bit1]
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(0u64), Fr::from(0u64), Fr::from(0u64)]), Fr::from(1u64)); // M[0,0]
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(1u64), Fr::from(1u64), Fr::from(1u64)]), Fr::from(16u64)); // M[3,3]
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(0u64), Fr::from(1u64), Fr::from(0u64)]), Fr::from(5u64)); // M[1,0]
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(0u64), Fr::from(0u64), Fr::from(1u64)]), Fr::from(10u64)); // M[2,1]
    }
    
    #[test]
    fn test_sparse_matrix_zero() {
        // Test with a zero matrix (must be power-of-2 dimensions)
        let sparse = SparseMatrix::<Fr>::new(2, 4);
        let mle = MLE::from_sparse_matrix(&sparse);
        
        // Should create an MLE with all zero coefficients
        assert_eq!(mle.num_vars(), 3); // log2(2*4) = 3
        assert_eq!(mle.get_split_count(), 2); // log2(4) = 2
        
        for coeff in &mle.coefficients {
            assert_eq!(*coeff, Fr::ZERO);
        }
    }

    #[test]
    fn test_mle_to_sparse_polynomial() {
        // Create MLE: 1 + 2*x_0 + 3*x_1 + 4*x_0*x_1
        // Coefficients: [1, 2, 3, 4]
        let coefficients = vec![
            Fr::from(1u64),  // constant term
            Fr::from(2u64),  // x_0
            Fr::from(3u64),  // x_1
            Fr::from(4u64),  // x_0*x_1
        ];
        let mle = MLE::new(coefficients);
        
        // Convert to SparsePolynomial
        let poly = mle.to_sparse_polynomial();
        
        // Verify the polynomial has correct number of variables
        assert_eq!(poly.num_vars, 2);
        
        // Verify evaluation at several points
        // At (0, 0): should be 1
        assert_eq!(poly.evaluate(&vec![Fr::ZERO, Fr::ZERO]), Fr::from(1u64));
        
        // At (1, 0): should be 1 + 2 = 3
        assert_eq!(poly.evaluate(&vec![Fr::ONE, Fr::ZERO]), Fr::from(3u64));
        
        // At (0, 1): should be 1 + 3 = 4
        assert_eq!(poly.evaluate(&vec![Fr::ZERO, Fr::ONE]), Fr::from(4u64));
        
        // At (1, 1): should be 1 + 2 + 3 + 4 = 10
        assert_eq!(poly.evaluate(&vec![Fr::ONE, Fr::ONE]), Fr::from(10u64));
    }
    
    #[test]
    fn test_sparse_polynomial_to_mle() {
        // Create polynomial: 5 + 7*x_0 + 11*x_1 + 13*x_0*x_1
        let mut terms = Vec::new();
        terms.push((Fr::from(5u64), SparseTerm::new(vec![])));  // constant
        terms.push((Fr::from(7u64), SparseTerm::new(vec![(0, 1)])));  // x_0
        terms.push((Fr::from(11u64), SparseTerm::new(vec![(1, 1)])));  // x_1
        terms.push((Fr::from(13u64), SparseTerm::new(vec![(0, 1), (1, 1)])));  // x_0*x_1
        
        let poly = SparsePolynomial::from_coefficients_vec(2, terms);
        
        // Convert to MLE
        let mle = MLE::from_sparse_polynomial(&poly);
        
        // Verify MLE structure
        assert_eq!(mle.num_vars, 2);
        assert_eq!(mle.coefficients.len(), 4);
        
        // Verify coefficients
        assert_eq!(mle.coefficients[0], Fr::from(5u64));   // constant
        assert_eq!(mle.coefficients[1], Fr::from(7u64));   // x_0
        assert_eq!(mle.coefficients[2], Fr::from(11u64));  // x_1
        assert_eq!(mle.coefficients[3], Fr::from(13u64));  // x_0*x_1
        
        // Verify evaluation
        assert_eq!(mle.evaluate(&vec![Fr::ZERO, Fr::ZERO]), Fr::from(5u64));
        assert_eq!(mle.evaluate(&vec![Fr::ONE, Fr::ZERO]), Fr::from(12u64));  // 5+7
        assert_eq!(mle.evaluate(&vec![Fr::ZERO, Fr::ONE]), Fr::from(16u64));  // 5+11
        assert_eq!(mle.evaluate(&vec![Fr::ONE, Fr::ONE]), Fr::from(36u64));   // 5+7+11+13
    }
    
    #[test]
    fn test_mle_sparse_polynomial_round_trip() {
        // Create MLE with 3 variables
        let coefficients = vec![
            Fr::from(1u64),
            Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64),
            Fr::from(5u64),
            Fr::from(6u64),
            Fr::from(7u64),
            Fr::from(8u64),
        ];
        let original_mle = MLE::new(coefficients.clone());
        
        // Convert to SparsePolynomial and back
        let poly = original_mle.to_sparse_polynomial();
        let recovered_mle = MLE::from_sparse_polynomial(&poly);
        
        // Verify the coefficients match
        assert_eq!(recovered_mle.num_vars, original_mle.num_vars);
        assert_eq!(recovered_mle.coefficients, original_mle.coefficients);
        
        // Verify evaluations match at several points
        let test_points = vec![
            vec![Fr::ZERO, Fr::ZERO, Fr::ZERO],
            vec![Fr::ONE, Fr::ZERO, Fr::ZERO],
            vec![Fr::ZERO, Fr::ONE, Fr::ZERO],
            vec![Fr::ONE, Fr::ONE, Fr::ZERO],
            vec![Fr::ZERO, Fr::ZERO, Fr::ONE],
            vec![Fr::ONE, Fr::ZERO, Fr::ONE],
            vec![Fr::ZERO, Fr::ONE, Fr::ONE],
            vec![Fr::ONE, Fr::ONE, Fr::ONE],
        ];
        
        for point in test_points {
            let mle_eval = original_mle.evaluate(&point);
            let poly_eval = poly.evaluate(&point);
            let recovered_eval = recovered_mle.evaluate(&point);
            assert_eq!(mle_eval, poly_eval);
            assert_eq!(mle_eval, recovered_eval);
        }
    }
    
    #[test]
    #[should_panic(expected = "Polynomial must be multilinear")]
    fn test_non_multilinear_polynomial_panics() {
        // Create polynomial with x_0^2 (not multilinear)
        let terms = vec![
            (Fr::from(1u64), SparseTerm::new(vec![(0, 2)])),  // x_0^2
        ];
        
        let poly = SparsePolynomial::from_coefficients_vec(2, terms);
        
        // This should panic
        let _mle = MLE::from_sparse_polynomial(&poly);
    }

    #[test]
    fn test_partial_evaluate_basic() {
        // Test with a 3-variable MLE f(x_1, x_2, x_3)
        // Evaluations: [1, 2, 3, 4, 5, 6, 7, 8] at (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1), (1,1,1)
        let evaluations: Vec<Fr> = (1..=8).map(|i| Fr::from(i as u64)).collect();
        let mle = MLE::from_evaluations(&evaluations);
        
        // Partial evaluate binding x_1 to 0
        let partial = mle.partial_evaluate(&[Fr::from(0u64)]);
        assert_eq!(partial.num_vars(), 2);
        
        // Should get f(0, x_2, x_3) which has evaluations at (0,0), (1,0), (0,1), (1,1)
        let partial_evals = partial.evaluations();
        assert_eq!(partial_evals[0], Fr::from(1u64)); // f(0,0,0)
        assert_eq!(partial_evals[1], Fr::from(3u64)); // f(0,1,0)
        assert_eq!(partial_evals[2], Fr::from(5u64)); // f(0,0,1)
        assert_eq!(partial_evals[3], Fr::from(7u64)); // f(0,1,1)
        
        // Partial evaluate binding x_1 to 1
        let partial = mle.partial_evaluate(&[Fr::from(1u64)]);
        assert_eq!(partial.num_vars(), 2);
        let partial_evals = partial.evaluations();
        assert_eq!(partial_evals[0], Fr::from(2u64)); // f(1,0,0)
        assert_eq!(partial_evals[1], Fr::from(4u64)); // f(1,1,0)
        assert_eq!(partial_evals[2], Fr::from(6u64)); // f(1,0,1)
        assert_eq!(partial_evals[3], Fr::from(8u64)); // f(1,1,1)
    }

    #[test]
    fn test_partial_evaluate_multiple_vars() {
        // Test with a 3-variable MLE, binding the first 2 variables
        let evaluations: Vec<Fr> = (1..=8).map(|i| Fr::from(i as u64)).collect();
        let mle = MLE::from_evaluations(&evaluations);
        
        // Partial evaluate binding x_1 and x_2 to (0, 0)
        let partial = mle.partial_evaluate(&[Fr::from(0u64), Fr::from(0u64)]);
        assert_eq!(partial.num_vars(), 1);
        let partial_evals = partial.evaluations();
        assert_eq!(partial_evals[0], Fr::from(1u64)); // f(0,0,0)
        assert_eq!(partial_evals[1], Fr::from(5u64)); // f(0,0,1)
        
        // Partial evaluate binding x_1 and x_2 to (1, 1)
        let partial = mle.partial_evaluate(&[Fr::from(1u64), Fr::from(1u64)]);
        assert_eq!(partial.num_vars(), 1);
        let partial_evals = partial.evaluations();
        assert_eq!(partial_evals[0], Fr::from(4u64)); // f(1,1,0)
        assert_eq!(partial_evals[1], Fr::from(8u64)); // f(1,1,1)
    }

    #[test]
    fn test_partial_evaluate_all_vars() {
        // Test that partial_evaluate with all variables equals evaluate
        let evaluations = vec![
            Fr::from(3u64),
            Fr::from(5u64),
            Fr::from(7u64),
            Fr::from(9u64),
        ];
        let mle = MLE::from_evaluations(&evaluations);
        
        // Test several points
        let test_points = vec![
            vec![Fr::from(0u64), Fr::from(0u64)],
            vec![Fr::from(1u64), Fr::from(0u64)],
            vec![Fr::from(0u64), Fr::from(1u64)],
            vec![Fr::from(1u64), Fr::from(1u64)],
        ];
        
        for point in test_points {
            let partial_result = mle.partial_evaluate(&point);
            assert_eq!(partial_result.num_vars(), 0);
            assert_eq!(partial_result.coefficients[0], mle.evaluate(&point));
        }
    }

    #[test]
    fn test_partial_evaluate_empty() {
        // Test that partial_evaluate with empty point returns the same MLE
        let evaluations = vec![
            Fr::from(1u64), Fr::from(2u64),
            Fr::from(3u64), Fr::from(4u64),
        ];
        let mle = MLE::from_evaluations(&evaluations);
        
        let partial = mle.partial_evaluate(&[]);
        assert_eq!(partial.num_vars(), mle.num_vars());
        assert_eq!(partial.coefficients(), mle.coefficients());
    }

    #[test]
    fn test_partial_evaluate_with_field_values() {
        // Test partial_evaluate with non-boolean field values
        // f(x_1, x_2) = 1 + 2*x_1 + 3*x_2 + 4*x_1*x_2
        let coefficients = vec![
            Fr::from(1u64),  // constant
            Fr::from(2u64),  // x_1
            Fr::from(3u64),  // x_2
            Fr::from(4u64),  // x_1*x_2
        ];
        let mle = MLE::new(coefficients);
        
        // Partial evaluate at x_1 = 5
        // Expected: f(5, x_2) = 1 + 2*5 + 3*x_2 + 4*5*x_2 = 11 + 23*x_2
        let partial = mle.partial_evaluate(&[Fr::from(5u64)]);
        assert_eq!(partial.num_vars(), 1);
        assert_eq!(partial.coefficients[0], Fr::from(11u64)); // constant term
        assert_eq!(partial.coefficients[1], Fr::from(23u64)); // coefficient of x_2
        
        // Verify by evaluating the partial result
        assert_eq!(partial.evaluate(&[Fr::from(0u64)]), Fr::from(11u64)); // f(5, 0) = 11
        assert_eq!(partial.evaluate(&[Fr::from(1u64)]), Fr::from(34u64)); // f(5, 1) = 11 + 23 = 34
        assert_eq!(partial.evaluate(&[Fr::from(2u64)]), Fr::from(57u64)); // f(5, 2) = 11 + 46 = 57
    }

    #[test]
    fn test_evaluate_uses_partial_evaluate() {
        // Verify that the refactored evaluate function produces correct results
        let coefficients = vec![
            Fr::from(3u64),  // constant
            Fr::from(2u64),  // x_1
            Fr::from(4u64),  // x_2
            Fr::from(0u64),  // x_1*x_2
        ];
        let mle = MLE::new(coefficients);
        
        // Test evaluation at boolean points
        assert_eq!(mle.evaluate(&[Fr::from(0u64), Fr::from(0u64)]), Fr::from(3u64)); // f(0,0) = 3
        assert_eq!(mle.evaluate(&[Fr::from(1u64), Fr::from(0u64)]), Fr::from(5u64)); // f(1,0) = 5
        assert_eq!(mle.evaluate(&[Fr::from(0u64), Fr::from(1u64)]), Fr::from(7u64)); // f(0,1) = 7
        assert_eq!(mle.evaluate(&[Fr::from(1u64), Fr::from(1u64)]), Fr::from(9u64)); // f(1,1) = 9
        
        // Test at intermediate point
        let half = Fr::from(2u64).inverse().unwrap();
        let result = mle.evaluate(&[half, half]);
        assert_eq!(result, Fr::from(6u64)); // f(0.5, 0.5) = 3 + 1 + 2 = 6
    }

    #[test]
    fn test_polynomial_multiply_simple() {
        // Test: f(x) = 1 + x, g(x) = 1 + x
        // Result: f*g = (1+x)(1+x) = 1 + 2x + x^2
        let mle1 = MLE::from_evaluations(&[Fr::from(1u64), Fr::from(2u64)]);
        let mle2 = MLE::from_evaluations(&[Fr::from(1u64), Fr::from(2u64)]);
        
        let product = mle1.polynomial_multiply(&mle2);
        
        // Verify at several points
        // At x=0: (1+0)(1+0) = 1
        assert_eq!(product.evaluate(&vec![Fr::from(0u64)]), Fr::from(1u64));
        // At x=1: (1+1)(1+1) = 4
        assert_eq!(product.evaluate(&vec![Fr::from(1u64)]), Fr::from(4u64));
        // At x=2: (1+2)(1+2) = 9
        assert_eq!(product.evaluate(&vec![Fr::from(2u64)]), Fr::from(9u64));
    }

    #[test]
    fn test_polynomial_multiply_two_vars() {
        // Test: f(x,y) = 1 + x, g(x,y) = 1 + y
        // Evaluations for f: [1, 2, 1, 2] at (0,0), (1,0), (0,1), (1,1)
        // Evaluations for g: [1, 1, 2, 2]
        let mle1 = MLE::from_evaluations(&[
            Fr::from(1u64), Fr::from(2u64), 
            Fr::from(1u64), Fr::from(2u64)
        ]);
        let mle2 = MLE::from_evaluations(&[
            Fr::from(1u64), Fr::from(1u64),
            Fr::from(2u64), Fr::from(2u64)
        ]);
        
        let product = mle1.polynomial_multiply(&mle2);
        
        // (1+x)(1+y) = 1 + x + y + xy
        // At (0,0): 1*1 = 1
        assert_eq!(product.evaluate(&vec![Fr::from(0u64), Fr::from(0u64)]), Fr::from(1u64));
        // At (1,0): 2*1 = 2
        assert_eq!(product.evaluate(&vec![Fr::from(1u64), Fr::from(0u64)]), Fr::from(2u64));
        // At (0,1): 1*2 = 2
        assert_eq!(product.evaluate(&vec![Fr::from(0u64), Fr::from(1u64)]), Fr::from(2u64));
        // At (1,1): 2*2 = 4
        assert_eq!(product.evaluate(&vec![Fr::from(1u64), Fr::from(1u64)]), Fr::from(4u64));
        // At (2,3): (1+2)(1+3) = 3*4 = 12
        assert_eq!(product.evaluate(&vec![Fr::from(2u64), Fr::from(3u64)]), Fr::from(12u64));
    }

    #[test]
    fn test_polynomial_multiply_produces_quadratic() {
        // Test that multiplying f(x) = x with g(x) = x produces x^2
        // f has evaluations [0, 1], g has evaluations [0, 1]
        let mle1 = MLE::from_evaluations(&[Fr::from(0u64), Fr::from(1u64)]);
        let mle2 = MLE::from_evaluations(&[Fr::from(0u64), Fr::from(1u64)]);
        
        let product = mle1.polynomial_multiply(&mle2);
        
        // x * x = x^2
        // At x=0: 0
        assert_eq!(product.evaluate(&vec![Fr::from(0u64)]), Fr::from(0u64));
        // At x=1: 1
        assert_eq!(product.evaluate(&vec![Fr::from(1u64)]), Fr::from(1u64));
        // At x=2: 4
        assert_eq!(product.evaluate(&vec![Fr::from(2u64)]), Fr::from(4u64));
        // At x=3: 9
        assert_eq!(product.evaluate(&vec![Fr::from(3u64)]), Fr::from(9u64));
    }

    #[test]
    fn test_polynomial_multiply_constant() {
        // Test multiplying by a constant
        // f(x) = 1 + x (evals [1, 2]), g(x) = 3 (evals [3, 3])
        let mle1 = MLE::from_evaluations(&[Fr::from(1u64), Fr::from(2u64)]);
        let mle2 = MLE::from_evaluations(&[Fr::from(3u64), Fr::from(3u64)]);
        
        let product = mle1.polynomial_multiply(&mle2);
        
        // (1+x)*3 = 3 + 3x
        // At x=0: 3
        assert_eq!(product.evaluate(&vec![Fr::from(0u64)]), Fr::from(3u64));
        // At x=1: 6
        assert_eq!(product.evaluate(&vec![Fr::from(1u64)]), Fr::from(6u64));
        // At x=5: 18
        assert_eq!(product.evaluate(&vec![Fr::from(5u64)]), Fr::from(18u64));
    }

    #[test]
    fn test_from_sparse_matrix() {
        // Create a 2x2 sparse matrix
        // M = [[1, 2],
        //      [3, 4]]
        let matrix_data = vec![
            vec![Fr::from(1u64), Fr::from(2u64)],
            vec![Fr::from(3u64), Fr::from(4u64)],
        ];
        let sparse_matrix = SparseMatrix::from_dense(&matrix_data);
        
        // Create MLE from sparse matrix
        let mle = MLE::from_sparse_matrix(&sparse_matrix);
        
        // Verify the MLE has correct number of variables
        // For 2x2 matrix: log2(2) + log2(2) = 1 + 1 = 2
        assert_eq!(mle.num_vars(), 2);
        
        // Note: The variable encoding is column-major for the first m2 vars, then row for m1 vars
        // For a 2x2 matrix with vars [x0, x1]:
        // mle.evaluate([x0, x1]) = M[x1, x0]
        // So: x0 encodes column index, x1 encodes row index
        
        // Verify evaluations match matrix values  
        // M[0,0] = 1 -> mle.evaluate([col=0, row=0]) = mle.evaluate([0,0])
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(0u64)]), Fr::from(1u64));
        // M[0,1] = 2 -> mle.evaluate([col=1, row=0]) = mle.evaluate([1,0])
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(0u64)]), Fr::from(2u64));
        // M[1,0] = 3 -> mle.evaluate([col=0, row=1]) = mle.evaluate([0,1])
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(1u64)]), Fr::from(3u64));
        // M[1,1] = 4 -> mle.evaluate([col=1, row=1]) = mle.evaluate([1,1])
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(1u64)]), Fr::from(4u64));
    }

    #[test]
    fn test_from_sparse_matrix_non_power_of_two() {
        // Create a 2x3 sparse matrix (non-power-of-two columns)
        // M = [[1, 2, 3],
        //      [4, 5, 6]]
        let matrix_data = vec![
            vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)],
            vec![Fr::from(4u64), Fr::from(5u64), Fr::from(6u64)],
        ];
        let sparse_matrix = SparseMatrix::from_dense(&matrix_data);
        
        // Create MLE from sparse matrix
        let mle = MLE::from_sparse_matrix(&sparse_matrix);
        
        // For 2x3 matrix, padded to 2x4
        // num_vars = log2(2) + log2(4) = 1 + 2 = 3
        assert_eq!(mle.num_vars(), 3);
        
        // Variable encoding: [x0, x1, x2] where x0,x1 encode column (log2(4)=2 bits), x2 encodes row (log2(2)=1 bit)
        // M[row, col] -> mle.evaluate([col_bit0, col_bit1, row_bit0])
        // M[0,0] -> mle.evaluate([0,0,0])
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(0u64), Fr::from(0u64)]), Fr::from(1u64));
        // M[0,1] -> mle.evaluate([1,0,0])
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(0u64), Fr::from(0u64)]), Fr::from(2u64));
        // M[0,2] -> mle.evaluate([0,1,0])
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(1u64), Fr::from(0u64)]), Fr::from(3u64));
        // M[1,0] -> mle.evaluate([0,0,1])
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(0u64), Fr::from(1u64)]), Fr::from(4u64));
        // M[1,1] -> mle.evaluate([1,0,1])
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(0u64), Fr::from(1u64)]), Fr::from(5u64));
        // M[1,2] -> mle.evaluate([0,1,1])
        assert_eq!(mle.evaluate(&vec![Fr::from(0u64), Fr::from(1u64), Fr::from(1u64)]), Fr::from(6u64));
        
        // Padded positions should be zero  
        // M[0,3] -> mle.evaluate([1,1,0])
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(1u64), Fr::from(0u64)]), Fr::from(0u64));
        // M[1,3] -> mle.evaluate([1,1,1])
        assert_eq!(mle.evaluate(&vec![Fr::from(1u64), Fr::from(1u64), Fr::from(1u64)]), Fr::from(0u64));
    }
}

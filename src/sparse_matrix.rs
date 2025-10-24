use ark_ff::Field;
use ark_std::vec::Vec;
use core::ops::Mul;

/// Compressed Sparse Row (CSR) format sparse matrix
/// 
/// CSR format stores a sparse matrix using three arrays:
/// - `values`: non-zero values
/// - `col_indices`: column index for each non-zero value
/// - `row_ptr`: pointers to the start of each row in `values` and `col_indices`
/// 
/// For a matrix with `m` rows and `n` columns:
/// - `row_ptr` has length `m + 1`
/// - `values` and `col_indices` have length equal to the number of non-zero elements
/// 
/// The elements of row `i` are stored in `values[row_ptr[i]..row_ptr[i+1]]`
/// 
/// Note: The number of rows can be derived from row_ptr.len() - 1
#[derive(Debug, Clone, PartialEq)]
pub struct SparseMatrix<F: Field> {
    /// Number of columns
    cols: usize,
    /// Non-zero values in row-major order
    values: Vec<F>,
    /// Column indices corresponding to each value
    col_indices: Vec<usize>,
    /// Row pointers: row_ptr[i] is the index in values/col_indices where row i starts
    /// row_ptr has length rows + 1, where row_ptr[rows] = values.len()
    row_ptr: Vec<usize>,
}

impl<F: Field> SparseMatrix<F> {
    /// Create a new empty sparse matrix with given dimensions
    pub fn new(rows: usize, cols: usize) -> Self {
        Self {
            cols,
            values: Vec::new(),
            col_indices: Vec::new(),
            row_ptr: vec![0; rows + 1],
        }
    }

    /// Get the number of rows
    pub fn rows(&self) -> usize {
        self.row_ptr.len().saturating_sub(1)
    }

    /// Create a sparse matrix from dense data
    /// 
    /// # Arguments
    /// * `data` - Dense matrix in row-major order (rows x cols)
    pub fn from_dense(data: &[Vec<F>]) -> Self {
        if data.is_empty() {
            return Self::new(0, 0);
        }

        let cols = data[0].len();
        let mut values = Vec::new();
        let mut col_indices = Vec::new();
        let mut row_ptr = vec![0];

        for row in data {
            assert_eq!(row.len(), cols, "All rows must have the same length");
            for (col_idx, val) in row.iter().enumerate() {
                if !val.is_zero() {
                    values.push(*val);
                    col_indices.push(col_idx);
                }
            }
            row_ptr.push(values.len());
        }

        Self {
            cols,
            values,
            col_indices,
            row_ptr,
        }
    }

    /// Create a sparse matrix from CSR format components
    pub fn from_csr(
        cols: usize,
        values: Vec<F>,
        col_indices: Vec<usize>,
        row_ptr: Vec<usize>,
    ) -> Result<Self, &'static str> {
        // Validate the CSR format
        if row_ptr.is_empty() {
            return Err("row_ptr cannot be empty");
        }
        if values.len() != col_indices.len() {
            return Err("values and col_indices must have the same length");
        }
        if row_ptr[0] != 0 {
            return Err("row_ptr must start with 0");
        }
        if *row_ptr.last().unwrap() != values.len() {
            return Err("last element of row_ptr must equal values.len()");
        }
        let rows = row_ptr.len() - 1;
        for i in 0..rows {
            if row_ptr[i] > row_ptr[i + 1] {
                return Err("row_ptr must be non-decreasing");
            }
        }
        for &col in &col_indices {
            if col >= cols {
                return Err("column index out of bounds");
            }
        }

        Ok(Self {
            cols,
            values,
            col_indices,
            row_ptr,
        })
    }

    /// Get the number of columns
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Get the number of non-zero elements
    pub fn nnz(&self) -> usize {
        self.values.len()
    }

    /// Get element at (row, col)
    pub fn get(&self, row: usize, col: usize) -> F {
        assert!(row < self.rows(), "Row index out of bounds");
        assert!(col < self.cols, "Column index out of bounds");

        let start = self.row_ptr[row];
        let end = self.row_ptr[row + 1];

        for i in start..end {
            if self.col_indices[i] == col {
                return self.values[i];
            }
        }

        F::ZERO
    }

    /// Set element at (row, col)
    /// Note: This is inefficient for CSR format. Use from_dense or from_csr for bulk construction.
    pub fn set(&mut self, row: usize, col: usize, value: F) {
        let rows = self.rows();
        assert!(row < rows, "Row index out of bounds");
        assert!(col < self.cols, "Column index out of bounds");

        let start = self.row_ptr[row];
        let end = self.row_ptr[row + 1];

        // Check if the element already exists
        for i in start..end {
            if self.col_indices[i] == col {
                if value.is_zero() {
                    // Remove the element
                    self.values.remove(i);
                    self.col_indices.remove(i);
                    for r in (row + 1)..=rows {
                        self.row_ptr[r] -= 1;
                    }
                } else {
                    // Update the value
                    self.values[i] = value;
                }
                return;
            } else if self.col_indices[i] > col {
                // Insert at position i
                if !value.is_zero() {
                    self.values.insert(i, value);
                    self.col_indices.insert(i, col);
                    for r in (row + 1)..=rows {
                        self.row_ptr[r] += 1;
                    }
                }
                return;
            }
        }

        // Element doesn't exist, insert at end of row
        if !value.is_zero() {
            self.values.insert(end, value);
            self.col_indices.insert(end, col);
            for r in (row + 1)..=rows {
                self.row_ptr[r] += 1;
            }
        }
    }

    /// Get a reference to the values array
    pub fn values(&self) -> &[F] {
        &self.values
    }

    /// Get a reference to the column indices array
    pub fn col_indices(&self) -> &[usize] {
        &self.col_indices
    }

    /// Get a reference to the row pointers array
    pub fn row_ptr(&self) -> &[usize] {
        &self.row_ptr
    }

    /// Matrix-vector multiplication: result = self * vec
    /// 
    /// Deprecated: Use the `Mul` trait instead: `matrix * &vec`
    #[deprecated(since = "0.1.0", note = "Use the `Mul` trait instead: `matrix.mul_vec(&vec)`")]
    pub fn mul_vec(&self, vec: &[F]) -> Vec<F> {
        self * vec
    }

    /// Create an identity matrix
    pub fn identity(size: usize) -> Self {
        let mut values = Vec::with_capacity(size);
        let mut col_indices = Vec::with_capacity(size);
        let mut row_ptr = Vec::with_capacity(size + 1);

        for i in 0..size {
            values.push(F::ONE);
            col_indices.push(i);
            row_ptr.push(i);
        }
        row_ptr.push(size);

        Self {
            cols: size,
            values,
            col_indices,
            row_ptr,
        }
    }

    /// Create a zero matrix
    pub fn zero(rows: usize, cols: usize) -> Self {
        Self::new(rows, cols)
    }

    /// Transpose the matrix
    pub fn transpose(&self) -> Self {
        let rows = self.rows();
        let cols = self.cols;
        
        // Count non-zeros in each column
        let mut col_counts = vec![0; cols];
        for &col in &self.col_indices {
            col_counts[col] += 1;
        }

        // Build row_ptr for transposed matrix
        let mut new_row_ptr = vec![0];
        for count in col_counts {
            new_row_ptr.push(new_row_ptr.last().unwrap() + count);
        }

        // Fill values and col_indices
        let nnz = self.nnz();
        let mut new_values = vec![F::ZERO; nnz];
        let mut new_col_indices = vec![0; nnz];
        let mut insert_positions = new_row_ptr[..cols].to_vec();

        for row in 0..rows {
            let start = self.row_ptr[row];
            let end = self.row_ptr[row + 1];

            for i in start..end {
                let col = self.col_indices[i];
                let pos = insert_positions[col];
                new_values[pos] = self.values[i];
                new_col_indices[pos] = row;
                insert_positions[col] += 1;
            }
        }

        Self {
            cols: rows,
            values: new_values,
            col_indices: new_col_indices,
            row_ptr: new_row_ptr,
        }
    }
}

// Implement Mul trait for matrix-vector multiplication (SparseMatrix * &[F])
impl<F: Field> Mul<&[F]> for &SparseMatrix<F> {
    type Output = Vec<F>;

    fn mul(self, vec: &[F]) -> Self::Output {
        assert_eq!(
            vec.len(),
            self.cols,
            "Vector length must match number of columns"
        );

        let rows = self.rows();
        let mut result = vec![F::ZERO; rows];

        for row in 0..rows {
            let start = self.row_ptr[row];
            let end = self.row_ptr[row + 1];

            for i in start..end {
                result[row] += self.values[i] * vec[self.col_indices[i]];
            }
        }

        result
    }
}

// Implement Mul trait for matrix-vector multiplication (SparseMatrix * Vec<F>)
impl<F: Field> Mul<Vec<F>> for &SparseMatrix<F> {
    type Output = Vec<F>;

    fn mul(self, vec: Vec<F>) -> Self::Output {
        self * &vec[..]
    }
}

// Implement Mul trait for matrix-vector multiplication (SparseMatrix * &Vec<F>)
impl<F: Field> Mul<&Vec<F>> for &SparseMatrix<F> {
    type Output = Vec<F>;

    fn mul(self, vec: &Vec<F>) -> Self::Output {
        self * &vec[..]
    }
}

// Implement Mul trait for matrix-matrix multiplication (SparseMatrix * SparseMatrix)
impl<F: Field> Mul<&SparseMatrix<F>> for &SparseMatrix<F> {
    type Output = SparseMatrix<F>;

    fn mul(self, other: &SparseMatrix<F>) -> Self::Output {
        assert_eq!(
            self.cols, other.rows(),
            "Matrix dimensions incompatible for multiplication"
        );

        let m = self.rows();
        let n = other.cols;

        // Transpose the second matrix for efficient column access
        let other_t = other.transpose();

        let mut result_data = vec![vec![F::ZERO; n]; m];

        // Compute C[i][j] = A[i] · B^T[j]
        for i in 0..m {
            let a_start = self.row_ptr[i];
            let a_end = self.row_ptr[i + 1];

            for j in 0..n {
                let b_start = other_t.row_ptr[j];
                let b_end = other_t.row_ptr[j + 1];

                let mut sum = F::ZERO;
                let mut a_idx = a_start;
                let mut b_idx = b_start;

                // Merge two sorted arrays to compute dot product
                while a_idx < a_end && b_idx < b_end {
                    let a_col = self.col_indices[a_idx];
                    let b_col = other_t.col_indices[b_idx];

                    if a_col == b_col {
                        sum += self.values[a_idx] * other_t.values[b_idx];
                        a_idx += 1;
                        b_idx += 1;
                    } else if a_col < b_col {
                        a_idx += 1;
                    } else {
                        b_idx += 1;
                    }
                }

                result_data[i][j] = sum;
            }
        }

        SparseMatrix::from_dense(&result_data)
    }
}

// Implement Mul trait for owned matrix-matrix multiplication
impl<F: Field> Mul<SparseMatrix<F>> for &SparseMatrix<F> {
    type Output = SparseMatrix<F>;

    fn mul(self, other: SparseMatrix<F>) -> Self::Output {
        self * &other
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_empty_matrix() {
        let matrix = SparseMatrix::<Fr>::new(3, 4);
        assert_eq!(matrix.rows(), 3);
        assert_eq!(matrix.cols(), 4);
        assert_eq!(matrix.nnz(), 0);
        assert_eq!(matrix.get(0, 0), Fr::from(0u64));
    }

    #[test]
    fn test_from_dense() {
        let dense = vec![
            vec![Fr::from(1u64), Fr::from(0u64), Fr::from(2u64)],
            vec![Fr::from(0u64), Fr::from(3u64), Fr::from(0u64)],
            vec![Fr::from(4u64), Fr::from(0u64), Fr::from(5u64)],
        ];

        let sparse = SparseMatrix::from_dense(&dense);
        assert_eq!(sparse.rows(), 3);
        assert_eq!(sparse.cols(), 3);
        assert_eq!(sparse.nnz(), 5);
        assert_eq!(sparse.get(0, 0), Fr::from(1u64));
        assert_eq!(sparse.get(0, 2), Fr::from(2u64));
        assert_eq!(sparse.get(1, 1), Fr::from(3u64));
        assert_eq!(sparse.get(2, 0), Fr::from(4u64));
        assert_eq!(sparse.get(2, 2), Fr::from(5u64));
        assert_eq!(sparse.get(1, 0), Fr::from(0u64));
    }

    #[test]
    fn test_set_and_get() {
        let mut matrix = SparseMatrix::<Fr>::new(3, 3);
        
        matrix.set(0, 0, Fr::from(1u64));
        matrix.set(1, 2, Fr::from(5u64));
        matrix.set(2, 1, Fr::from(3u64));

        assert_eq!(matrix.get(0, 0), Fr::from(1u64));
        assert_eq!(matrix.get(1, 2), Fr::from(5u64));
        assert_eq!(matrix.get(2, 1), Fr::from(3u64));
        assert_eq!(matrix.get(0, 1), Fr::from(0u64));
        assert_eq!(matrix.nnz(), 3);
    }

    #[test]
    fn test_mul_vec() {
        let dense = vec![
            vec![Fr::from(1u64), Fr::from(2u64), Fr::from(0u64)],
            vec![Fr::from(0u64), Fr::from(3u64), Fr::from(4u64)],
            vec![Fr::from(5u64), Fr::from(0u64), Fr::from(6u64)],
        ];

        let sparse = SparseMatrix::from_dense(&dense);
        let vec = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        let result = sparse.mul_vec(&vec);

        // [1, 2, 0] * [1, 2, 3]^T = 1*1 + 2*2 + 0*3 = 5
        // [0, 3, 4] * [1, 2, 3]^T = 0*1 + 3*2 + 4*3 = 18
        // [5, 0, 6] * [1, 2, 3]^T = 5*1 + 0*2 + 6*3 = 23
        assert_eq!(result[0], Fr::from(5u64));
        assert_eq!(result[1], Fr::from(18u64));
        assert_eq!(result[2], Fr::from(23u64));
    }

    #[test]
    fn test_identity() {
        let identity = SparseMatrix::<Fr>::identity(3);
        assert_eq!(identity.rows(), 3);
        assert_eq!(identity.cols(), 3);
        assert_eq!(identity.nnz(), 3);
        assert_eq!(identity.get(0, 0), Fr::from(1u64));
        assert_eq!(identity.get(1, 1), Fr::from(1u64));
        assert_eq!(identity.get(2, 2), Fr::from(1u64));
        assert_eq!(identity.get(0, 1), Fr::from(0u64));
    }

    #[test]
    fn test_from_csr() {
        // Matrix:
        // [1, 0, 2]
        // [0, 3, 0]
        // [4, 0, 5]
        let values = vec![
            Fr::from(1u64), Fr::from(2u64),
            Fr::from(3u64),
            Fr::from(4u64), Fr::from(5u64),
        ];
        let col_indices = vec![0, 2, 1, 0, 2];
        let row_ptr = vec![0, 2, 3, 5];

        let sparse = SparseMatrix::from_csr(3, values, col_indices, row_ptr).unwrap();
        assert_eq!(sparse.rows(), 3);
        assert_eq!(sparse.cols(), 3);
        assert_eq!(sparse.nnz(), 5);
        assert_eq!(sparse.get(0, 0), Fr::from(1u64));
        assert_eq!(sparse.get(0, 2), Fr::from(2u64));
        assert_eq!(sparse.get(1, 1), Fr::from(3u64));
    }

    #[test]
    fn test_transpose() {
        let dense = vec![
            vec![Fr::from(1u64), Fr::from(2u64), Fr::from(0u64)],
            vec![Fr::from(0u64), Fr::from(3u64), Fr::from(4u64)],
        ];

        let sparse = SparseMatrix::from_dense(&dense);
        let transposed = sparse.transpose();

        // Transposed should be 3x2:
        // [1, 0]
        // [2, 3]
        // [0, 4]
        assert_eq!(transposed.rows(), 3);
        assert_eq!(transposed.cols(), 2);
        assert_eq!(transposed.get(0, 0), Fr::from(1u64));
        assert_eq!(transposed.get(0, 1), Fr::from(0u64));
        assert_eq!(transposed.get(1, 0), Fr::from(2u64));
        assert_eq!(transposed.get(1, 1), Fr::from(3u64));
        assert_eq!(transposed.get(2, 0), Fr::from(0u64));
        assert_eq!(transposed.get(2, 1), Fr::from(4u64));
    }

    #[test]
    fn test_mul_trait_vector() {
        let dense = vec![
            vec![Fr::from(1u64), Fr::from(2u64), Fr::from(0u64)],
            vec![Fr::from(0u64), Fr::from(3u64), Fr::from(4u64)],
            vec![Fr::from(5u64), Fr::from(0u64), Fr::from(6u64)],
        ];

        let sparse = SparseMatrix::from_dense(&dense);
        let vec = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
        
        // Using Mul trait
        let result = &sparse * &vec;

        // [1, 2, 0] * [1, 2, 3]^T = 1*1 + 2*2 + 0*3 = 5
        // [0, 3, 4] * [1, 2, 3]^T = 0*1 + 3*2 + 4*3 = 18
        // [5, 0, 6] * [1, 2, 3]^T = 5*1 + 0*2 + 6*3 = 23
        assert_eq!(result[0], Fr::from(5u64));
        assert_eq!(result[1], Fr::from(18u64));
        assert_eq!(result[2], Fr::from(23u64));

        // Test with owned vector
        let result2 = &sparse * vec.clone();
        assert_eq!(result2, result);
    }

    #[test]
    fn test_mul_trait_matrix() {
        // Matrix A (2x3):
        // [1, 2, 0]
        // [0, 3, 4]
        let dense_a = vec![
            vec![Fr::from(1u64), Fr::from(2u64), Fr::from(0u64)],
            vec![Fr::from(0u64), Fr::from(3u64), Fr::from(4u64)],
        ];

        // Matrix B (3x2):
        // [5, 6]
        // [0, 7]
        // [8, 0]
        let dense_b = vec![
            vec![Fr::from(5u64), Fr::from(6u64)],
            vec![Fr::from(0u64), Fr::from(7u64)],
            vec![Fr::from(8u64), Fr::from(0u64)],
        ];

        let matrix_a = SparseMatrix::from_dense(&dense_a);
        let matrix_b = SparseMatrix::from_dense(&dense_b);

        // C = A * B (2x2)
        let result = &matrix_a * &matrix_b;

        // C[0][0] = 1*5 + 2*0 + 0*8 = 5
        // C[0][1] = 1*6 + 2*7 + 0*0 = 20
        // C[1][0] = 0*5 + 3*0 + 4*8 = 32
        // C[1][1] = 0*6 + 3*7 + 4*0 = 21
        assert_eq!(result.rows(), 2);
        assert_eq!(result.cols(), 2);
        assert_eq!(result.get(0, 0), Fr::from(5u64));
        assert_eq!(result.get(0, 1), Fr::from(20u64));
        assert_eq!(result.get(1, 0), Fr::from(32u64));
        assert_eq!(result.get(1, 1), Fr::from(21u64));
    }

    #[test]
    fn test_identity_multiplication() {
        let dense = vec![
            vec![Fr::from(1u64), Fr::from(2u64)],
            vec![Fr::from(3u64), Fr::from(4u64)],
        ];

        let matrix = SparseMatrix::from_dense(&dense);
        let identity = SparseMatrix::identity(2);

        // I * A = A
        let result1 = &identity * &matrix;
        assert_eq!(result1.get(0, 0), Fr::from(1u64));
        assert_eq!(result1.get(0, 1), Fr::from(2u64));
        assert_eq!(result1.get(1, 0), Fr::from(3u64));
        assert_eq!(result1.get(1, 1), Fr::from(4u64));

        // A * I = A
        let result2 = &matrix * &identity;
        assert_eq!(result2.get(0, 0), Fr::from(1u64));
        assert_eq!(result2.get(0, 1), Fr::from(2u64));
        assert_eq!(result2.get(1, 0), Fr::from(3u64));
        assert_eq!(result2.get(1, 1), Fr::from(4u64));
    }
}


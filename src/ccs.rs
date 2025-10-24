use ark_ff::Field;
use ark_std::vec::Vec;
use crate::sparse_matrix::SparseMatrix;

/// Customizable Constraint System (CCS)
/// 
/// CCS is a generalization of constraint systems like R1CS, Plonkish, and AIR.
/// It represents constraints of the form:
/// 
/// ∑_i c_i · ⊙_{j ∈ S_i} M_j · z = 0
/// 
/// where:
/// - c_i are constants
/// - S_i are multisets of matrix indices
/// - M_j are constraint matrices
/// - z is the witness vector
/// - ⊙ denotes the Hadamard (element-wise) product
/// 
/// Each constraint is formed by:
/// 1. Multiplying matrices M_j by the witness vector z for each j in multiset S_i
/// 2. Taking the Hadamard product of these results
/// 3. Multiplying by constant c_i
/// 4. Summing over all i to get zero
#[derive(Debug, Clone)]
pub struct CCS<F: Field> {
    /// Constants for each constraint term
    pub constants: Vec<F>,
    /// Sets of matrix indices for Hadamard products
    /// Each multiset represents which matrices to multiply and hadamard product
    pub multisets: Vec<Vec<usize>>,
    /// Constraint matrices in CSR format
    pub matrices: Vec<SparseMatrix<F>>,
}

impl<F: Field> CCS<F> {
    /// Create a new empty CCS
    pub fn new() -> Self {
        Self {
            constants: Vec::new(),
            multisets: Vec::new(),
            matrices: Vec::new(),
        }
    }

    /// Create a CCS with given components
    pub fn new_with(
        constants: Vec<F>,
        multisets: Vec<Vec<usize>>,
        matrices: Vec<SparseMatrix<F>>,
    ) -> Result<Self, &'static str> {
        // Validate the structure
        if constants.len() != multisets.len() {
            return Err("constants and multisets must have the same length");
        }

        // Check that all matrix indices in multisets are valid
        let num_matrices = matrices.len();
        for multiset in &multisets {
            for &idx in multiset {
                if idx >= num_matrices {
                    return Err("matrix index in multiset out of bounds");
                }
            }
        }

        // Check that all matrices have compatible dimensions
        if !matrices.is_empty() {
            let rows = matrices[0].rows();
            let cols = matrices[0].cols();
            for matrix in &matrices {
                if matrix.rows() != rows || matrix.cols() != cols {
                    return Err("all matrices must have the same dimensions");
                }
            }
        }

        Ok(Self {
            constants,
            multisets,
            matrices,
        })
    }

    /// Get the number of constraint terms
    pub fn num_terms(&self) -> usize {
        self.constants.len()
    }

    /// Get the number of matrices
    pub fn num_matrices(&self) -> usize {
        self.matrices.len()
    }

    /// Get the number of constraints (rows in matrices)
    pub fn num_constraints(&self) -> usize {
        if self.matrices.is_empty() {
            0
        } else {
            self.matrices[0].rows()
        }
    }

    /// Get the number of variables (columns in matrices)
    pub fn num_variables(&self) -> usize {
        if self.matrices.is_empty() {
            0
        } else {
            self.matrices[0].cols()
        }
    }

    /// Check if a witness satisfies the constraints
    /// 
    /// Evaluates: ∑_i c_i · ⊙_{j ∈ S_i} M_j · z = 0
    pub fn is_satisfied(&self, witness: &[F]) -> bool {
        if self.matrices.is_empty() {
            return true;
        }

        assert_eq!(
            witness.len(),
            self.num_variables(),
            "Witness length must match number of variables"
        );

        let num_constraints = self.num_constraints();
        let mut result = vec![F::ZERO; num_constraints];

        // For each term
        for i in 0..self.num_terms() {
            let constant = self.constants[i];
            let multiset = &self.multisets[i];

            if multiset.is_empty() {
                continue;
            }

            // Compute Hadamard product of (M_j · z) for all j in multiset
            let mut hadamard_product = vec![F::ONE; num_constraints];

            for &matrix_idx in multiset {
                let matrix = &self.matrices[matrix_idx];
                let mat_vec_product = matrix.mul_vec(witness);

                // Hadamard product (element-wise multiplication)
                for k in 0..num_constraints {
                    hadamard_product[k] *= mat_vec_product[k];
                }
            }

            // Add c_i * hadamard_product to result
            for k in 0..num_constraints {
                result[k] += constant * hadamard_product[k];
            }
        }

        // Check if all constraints are satisfied (all elements are zero)
        result.iter().all(|&x| x.is_zero())
    }

    /// Add a new constraint term
    pub fn add_term(&mut self, constant: F, multiset: Vec<usize>) -> Result<(), &'static str> {
        // Validate matrix indices
        let num_matrices = self.matrices.len();
        for &idx in &multiset {
            if idx >= num_matrices {
                return Err("matrix index in multiset out of bounds");
            }
        }

        self.constants.push(constant);
        self.multisets.push(multiset);
        Ok(())
    }

    /// Add a new matrix
    pub fn add_matrix(&mut self, matrix: SparseMatrix<F>) -> Result<(), &'static str> {
        // Check dimension compatibility
        if !self.matrices.is_empty() {
            let rows = self.matrices[0].rows();
            let cols = self.matrices[0].cols();
            if matrix.rows() != rows || matrix.cols() != cols {
                return Err("matrix dimensions must match existing matrices");
            }
        }

        self.matrices.push(matrix);
        Ok(())
    }
}

impl<F: Field> Default for CCS<F> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_ccs_creation() {
        let ccs = CCS::<Fr>::new();
        assert_eq!(ccs.num_terms(), 0);
        assert_eq!(ccs.num_matrices(), 0);
    }

    #[test]
    fn test_ccs_new_with() {
        let matrix1 = SparseMatrix::<Fr>::identity(3);
        let matrix2 = SparseMatrix::<Fr>::zero(3, 3);

        let constants = vec![Fr::from(1u64), Fr::from(2u64)];
        let multisets = vec![vec![0], vec![1]];
        let matrices = vec![matrix1, matrix2];

        let ccs = CCS::new_with(constants, multisets, matrices).unwrap();
        assert_eq!(ccs.num_terms(), 2);
        assert_eq!(ccs.num_matrices(), 2);
        assert_eq!(ccs.num_constraints(), 3);
        assert_eq!(ccs.num_variables(), 3);
    }

    #[test]
    fn test_ccs_validation() {
        let matrix = SparseMatrix::<Fr>::identity(3);

        // Mismatched lengths
        let result = CCS::new_with(
            vec![Fr::from(1u64)],
            vec![vec![0], vec![0]],
            vec![matrix.clone()],
        );
        assert!(result.is_err());

        // Invalid matrix index
        let result = CCS::new_with(
            vec![Fr::from(1u64)],
            vec![vec![1]], // No matrix at index 1
            vec![matrix.clone()],
        );
        assert!(result.is_err());
    }

    #[test]
    fn test_simple_r1cs_as_ccs() {
        // Simple R1CS constraint: x * y = z
        // Can be expressed as: A·z ⊙ B·z - C·z = 0
        // where z = [1, x, y, z] (1 is the constant, rest are variables)

        // A matrix: [0, 1, 0, 0] (selects x)
        let mut matrix_a = SparseMatrix::<Fr>::new(1, 4);
        matrix_a.set(0, 1, Fr::from(1u64));

        // B matrix: [0, 0, 1, 0] (selects y)
        let mut matrix_b = SparseMatrix::<Fr>::new(1, 4);
        matrix_b.set(0, 2, Fr::from(1u64));

        // C matrix: [0, 0, 0, 1] (selects z)
        let mut matrix_c = SparseMatrix::<Fr>::new(1, 4);
        matrix_c.set(0, 3, Fr::from(1u64));

        let constants = vec![Fr::from(1u64), -Fr::from(1u64)];
        let multisets = vec![
            vec![0, 1], // A ⊙ B
            vec![2],    // C
        ];
        let matrices = vec![matrix_a, matrix_b, matrix_c];

        let ccs = CCS::new_with(constants, multisets, matrices).unwrap();

        // Test with x=3, y=4, z=12 (should be satisfied)
        let witness_satisfied = vec![
            Fr::from(1u64),  // constant 1
            Fr::from(3u64),  // x
            Fr::from(4u64),  // y
            Fr::from(12u64), // z
        ];
        assert!(ccs.is_satisfied(&witness_satisfied));

        // Test with x=3, y=4, z=10 (should NOT be satisfied)
        let witness_unsatisfied = vec![
            Fr::from(1u64),  // constant 1
            Fr::from(3u64),  // x
            Fr::from(4u64),  // y
            Fr::from(10u64), // z (wrong!)
        ];
        assert!(!ccs.is_satisfied(&witness_unsatisfied));
    }

    #[test]
    fn test_add_term_and_matrix() {
        let mut ccs = CCS::<Fr>::new();

        let matrix = SparseMatrix::<Fr>::identity(3);
        ccs.add_matrix(matrix).unwrap();

        ccs.add_term(Fr::from(1u64), vec![0]).unwrap();

        assert_eq!(ccs.num_terms(), 1);
        assert_eq!(ccs.num_matrices(), 1);
    }

    #[test]
    fn test_ccs_with_multiple_hadamard() {
        // Test with three-way Hadamard product
        let mut matrix1 = SparseMatrix::<Fr>::new(1, 3);
        matrix1.set(0, 0, Fr::from(2u64));

        let mut matrix2 = SparseMatrix::<Fr>::new(1, 3);
        matrix2.set(0, 1, Fr::from(3u64));

        let mut matrix3 = SparseMatrix::<Fr>::new(1, 3);
        matrix3.set(0, 2, Fr::from(4u64));

        let constants = vec![Fr::from(1u64)];
        let multisets = vec![vec![0, 1, 2]]; // Three-way Hadamard
        let matrices = vec![matrix1, matrix2, matrix3];

        let ccs = CCS::new_with(constants, multisets, matrices).unwrap();

        // M1·z = [2*1, 0, 0] = [2, 0, 0]
        // M2·z = [0, 3*5, 0] = [15, 0, 0] (wrong - let me recalculate)
        // Actually for a 1x3 matrix with one element at (0,1)=3, times vector [1,5,7]:
        // Result is 1-element vector: [3*5] = [15]

        let witness = vec![Fr::from(1u64), Fr::from(5u64), Fr::from(7u64)];
        
        // M1·[1,5,7] = [2*1] = [2]
        // M2·[1,5,7] = [3*5] = [15]
        // M3·[1,5,7] = [4*7] = [28]
        // Hadamard: [2] ⊙ [15] ⊙ [28] = [2*15*28] = [840]
        // 1 * [840] = [840], which is not zero
        
        assert!(!ccs.is_satisfied(&witness));
    }
}

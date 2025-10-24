pub mod sparse_matrix;
pub mod ccs;
pub mod mle;
pub mod gmuly_protocol;
pub mod utils;

pub use sparse_matrix::SparseMatrix;
pub use ccs::CCS;
pub use mle::MLE;
pub use gmuly_protocol::{GmulyProof, GmulyProver, GmulyVerifier, GmulyStatement, GmulyWitness};
pub use utils::{partialsum, all_sum, bind_variable, univariate_poly_to_coeffs, coeffs_to_univariate_poly};

#[cfg(test)]
mod tests {
    use super::*;
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_sparse_matrix_creation() {
        let matrix = SparseMatrix::<Fr>::new(3, 3);
        assert_eq!(matrix.rows(), 3);
        assert_eq!(matrix.cols(), 3);
    }

    #[test]
    fn test_ccs_creation() {
        let ccs = CCS::<Fr> {
            constants: vec![Fr::from(1u64), Fr::from(2u64)],
            multisets: vec![vec![0, 1], vec![1, 2]],
            matrices: vec![SparseMatrix::new(3, 3)],
        };
        assert_eq!(ccs.constants.len(), 2);
        assert_eq!(ccs.multisets.len(), 2);
        assert_eq!(ccs.matrices.len(), 1);
    }
}

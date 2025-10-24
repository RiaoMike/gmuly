use ark_ff::Field;
use ark_std::vec::Vec;
use crate::SparseMatrix;

/// Statement of Linear CCS
/// 包含矩阵、随机向量r以及求和值u
#[derive(Debug, Clone)]
pub struct GmulyStatement<F: Field> {
    /// 矩阵M
    pub matrix: SparseMatrix<F>,
    /// 随机向量r
    pub random_vector: Vec<F>,
    /// 求和值u
    pub sum_value: F,
}

impl<F: Field> GmulyStatement<F> {
    /// 创建新的Statement
    ///
    /// # Arguments
    /// * `matrix` - 矩阵M
    /// * `random_vector` - 随机向量r
    /// * `sum_value` - 求和值u
    pub fn new(matrix: SparseMatrix<F>, random_vector: Vec<F>, sum_value: F) -> Self {
        assert_eq!(
            random_vector.len(),
            matrix.rows().next_power_of_two().trailing_zeros() as usize,
        );

        Self {
            matrix,
            random_vector,
            sum_value,
        }
    }
}

/// Witness of Linear CCS
/// 包含一个向量
#[derive(Debug, Clone)]
pub struct GmulyWitness<F: Field> {
    /// 见证向量w
    pub vector: Vec<F>,
}

impl<F: Field> GmulyWitness<F> {
    /// 创建新的Witness
    ///
    /// # Arguments
    /// * `vector` - 见证向量w
    pub fn new(vector: Vec<F>) -> Self {
        Self { vector }
    }
}

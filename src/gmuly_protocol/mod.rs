//! Gmuly协议实现
//! 
//! 该模块包含了Gmuly协议的完整实现，包括:
//! - Statement和Witness结构 (lccs.rs)
//! - Prover实现 (prover.rs)
//! - Verifier实现 (verifier.rs)

pub mod lccs;
pub mod prover;
pub mod verifier;

// Re-export主要类型
pub use lccs::{GmulyStatement, GmulyWitness};
pub use prover::{GmulyProof, GmulyProver};
pub use verifier::GmulyVerifier;

#[cfg(test)]
mod tests {
    use super::*;
    use ark_test_curves::bls12_381::Fr;
    use ark_poly::{
        multivariate::{SparsePolynomial, SparseTerm, Term},
        DenseMVPolynomial, Polynomial,
    };
    use ark_ff::Field;
    use crate::{utils::partialsum, all_sum, MLE};

    #[test]
    fn test_gmuly_protocol_simple() {
        // 创建一个简单的2x2矩阵测试
        // M = [[1, 2],
        //      [3, 4]]
        let matrix_data = vec![
            vec![Fr::from(1u64), Fr::from(2u64)],
            vec![Fr::from(3u64), Fr::from(4u64)],
        ];
        let matrix = crate::SparseMatrix::from_dense(&matrix_data);
        
        // 创建witness向量 w = [1, 1]
        let witness_vector = vec![Fr::from(1u64), Fr::from(1u64)];
        
        // 计算 M * w = [3, 7]
        // 随机向量 r = [r0] (1个元素，因为结果向量长度是2，log2(2)=1)
        let random_vector = vec![Fr::from(5u64)]; // 使用固定值便于测试
        
        // 计算期望的和值
        // MLE(M*w) 在 r=[5] 处的值
        // (M*w)[0] = 3, (M*w)[1] = 7
        // MLE插值: f(0) = 3, f(1) = 7
        // f(x) = 3 + 4x
        // f(5) = 3 + 4*5 = 23
        let sum_value = Fr::from(23u64);
        
        // 创建Statement和Witness
        let statement = GmulyStatement::new(matrix, random_vector, sum_value);
        let witness = GmulyWitness::new(witness_vector);
        
        // 折叠次数设为1（对应witness向量的比特长度 log2(2) = 1）
        let fold = 1;
        
        // 创建Prover并生成证明
        let prover = GmulyProver::new(statement.clone(), witness, fold);
        let proof = prover.prove();
        
        // 创建Verifier并验证证明
        let verifier = GmulyVerifier::new(statement, fold);
        let is_valid = verifier.verify(&proof);
        
        assert!(is_valid, "Proof should be valid");
        println!("✓ Gmuly protocol test passed!");
    }

    #[test]
    fn test_gmuly_protocol_identity_matrix() {
        // 使用单位矩阵测试
        // I = [[1, 0],
        //      [0, 1]]
        let matrix_data = vec![
            vec![Fr::from(2u64), Fr::from(0u64), Fr::from(5u64)],
            vec![Fr::from(3u64), Fr::from(1u64), Fr::from(0u64)],
            vec![Fr::from(3u64), Fr::from(0u64), Fr::from(5u64)],
        ];
        let matrix = crate::SparseMatrix::from_dense(&matrix_data);
        
        let random_vector = vec![Fr::from(39u64), Fr::from(3329u64)];
        // witness向量 w = [2, 3, 1]
        let witness_vector = vec![Fr::from(2u64), Fr::from(3u64), Fr::from(1u64)];
        
        let poly = MLE::from_vector(&witness_vector);
        let poly_matrix = MLE::from_sparse_matrix(&matrix);
        let bound_poly_matrix = poly_matrix.partial_evaluate(&random_vector);
        let g_poly = poly.polynomial_multiply(&bound_poly_matrix);
        
        let sum_value = all_sum(&g_poly);
        
        let statement = GmulyStatement::new(matrix, random_vector, sum_value);
        let witness = GmulyWitness::new(witness_vector);
        let fold = 1;
        
        let prover = GmulyProver::new(statement.clone(), witness, fold);
        let proof = prover.prove();
        
        let verifier = GmulyVerifier::new(statement, fold);
        let is_valid = verifier.verify(&proof);
        
        assert!(is_valid, "Proof should be valid for identity matrix");
        println!("✓ Identity matrix test passed!");
    }

    #[test]
    fn test_partialsum_simple() {
        // 测试 f(x0, x1, x2) = x0 + x1 + x2
        // 对 x1, x2 求和 (start_index=1)
        // 结果应该是: 4*x0 + 4
        let mut terms = Vec::new();
        terms.push((Fr::from(1u64), SparseTerm::new(vec![(0, 1)]))); // x0
        terms.push((Fr::from(1u64), SparseTerm::new(vec![(1, 1)]))); // x1
        terms.push((Fr::from(1u64), SparseTerm::new(vec![(2, 1)]))); // x2
        let poly = SparsePolynomial::from_coefficients_vec(3, terms);

        let result = partialsum(&poly, 1);

        // 验证结果多项式的变量数
        assert_eq!(result.num_vars, 1, "Result should have 1 variable");

        // 验证在不同点的评估值
        // result(0) = 4*0 + 4 = 4
        assert_eq!(result.evaluate(&vec![Fr::ZERO]), Fr::from(4u64));
        // result(1) = 4*1 + 4 = 8
        assert_eq!(result.evaluate(&vec![Fr::ONE]), Fr::from(8u64));
    }

    #[test]
    fn test_partialsum_product() {
        // 测试 f(x0, x1) = x0 * x1
        // 对 x1 求和 (start_index=1)
        // 结果应该是: x0 (因为 sum_{x1∈{0,1}} x0*x1 = x0*0 + x0*1 = x0)
        let terms = vec![(Fr::from(1u64), SparseTerm::new(vec![(0, 1), (1, 1)]))];
        let poly = SparsePolynomial::from_coefficients_vec(2, terms);

        let result = partialsum(&poly, 1);

        // 验证结果
        assert_eq!(result.num_vars, 1);
        assert_eq!(result.evaluate(&vec![Fr::ZERO]), Fr::ZERO);
        assert_eq!(result.evaluate(&vec![Fr::ONE]), Fr::ONE);
    }

    #[test]
    fn test_partialsum_constant() {
        // 测试常数多项式 f(x0, x1) = 5
        // 对 x1 求和 (start_index=1)
        // 结果应该是: 10 (因为 sum_{x1∈{0,1}} 5 = 5 + 5 = 10)
        let terms = vec![(Fr::from(5u64), SparseTerm::new(vec![]))];
        let poly = SparsePolynomial::from_coefficients_vec(2, terms);

        let result = partialsum(&poly, 1);

        // 验证结果
        assert_eq!(result.num_vars, 1);
        // 结果应该是常数 10
        assert_eq!(result.evaluate(&vec![Fr::ZERO]), Fr::from(10u64));
        assert_eq!(result.evaluate(&vec![Fr::ONE]), Fr::from(10u64));
    }

    #[test]
    fn test_partialsum_all_vars() {
        // 测试 f(x0, x1) = x0 + x1 + 1
        // 对所有变量求和 (start_index=0)
        // f(0,0)=1, f(0,1)=2, f(1,0)=2, f(1,1)=3
        // 结果应该是: sum_{x0,x1∈{0,1}} (x0 + x1 + 1) = 1 + 2 + 2 + 3 = 8
        let mut terms = Vec::new();
        terms.push((Fr::from(1u64), SparseTerm::new(vec![(0, 1)]))); // x0
        terms.push((Fr::from(1u64), SparseTerm::new(vec![(1, 1)]))); // x1
        terms.push((Fr::from(1u64), SparseTerm::new(vec![]))); // 常数1
        let poly = SparsePolynomial::from_coefficients_vec(2, terms);

        let result = partialsum(&poly, 0);

        // 结果应该是常数多项式 8
        assert_eq!(result.num_vars, 0);
        assert_eq!(result.evaluate(&vec![]), Fr::from(8u64));
    }
}

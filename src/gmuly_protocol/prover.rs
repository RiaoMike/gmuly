use ark_ff::Field;
use ark_std::rand::{rngs::StdRng, SeedableRng};
use crate::{bind_variable, partialsum, univariate_poly_to_coeffs, GmulyStatement, GmulyWitness, MLE};
use ark_std::vec::Vec;

/// Gmuly协议的证明
#[derive(Debug, Clone)]
pub struct GmulyProof<F: Field> {
    /// 每一轮的单变量多项式系数
    /// round_polynomials[i]表示第i轮的单变量多项式
    pub round_polynomials: Vec<Vec<F>>,
}

/// Gmuly协议的Prover
pub struct GmulyProver<F: Field> {
    /// statement
    statement: GmulyStatement<F>,

    /// witness
    witness: GmulyWitness<F>,

    /// 折叠次数
    fold: usize,
}

impl<F: Field> GmulyProver<F> {
    /// 创建新的Gmuly Prover
    pub fn new(
        statement: GmulyStatement<F>,
        witness: GmulyWitness<F>,
        fold: usize,
    ) -> Self {
        Self {
            statement,
            witness,
            fold,
        }
    }

    pub fn prove(&self) -> GmulyProof<F> {
        let mut round_polynomials: Vec<Vec<F>> = Vec::new();
        
        // 将statement中的矩阵转换为MLE形式
        let bound_matrix_mle = MLE::from_sparse_matrix(&self.statement.matrix)
            .partial_evaluate(&self.statement.random_vector);
        // 将witness中的向量转换为MLE形式
        let witness_mle = MLE::from_vector(&self.witness.vector);

        // 计算两个 MLE 的多项式乘积，得到一个二次多项式
        let mut current_poly = bound_matrix_mle.polynomial_multiply(&witness_mle);
        current_poly = partialsum(&current_poly, self.fold);

        // 使用固定种子创建随机数生成器，验证方可以用相同种子重现随机挑战
        // 在实际应用中，这个种子应该由Fiat-Shamir变换从transcript中派生
        let seed = 42u64; // TODO: 从Fiat-Shamir transcript获取
        let mut rng = StdRng::seed_from_u64(seed);
        
        for i in 0..self.fold {
            let sum_g = partialsum(&current_poly, 1);
            round_polynomials.push(univariate_poly_to_coeffs(&sum_g));

            let random_element = F::rand(&mut rng);
            current_poly = bind_variable(&current_poly, 0, random_element);
        }

        // todo
        // 添加greyhound子证明系统
        GmulyProof {
            round_polynomials,
        }
    }
}

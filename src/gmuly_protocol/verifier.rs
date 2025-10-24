use ark_ff::Field;
use ark_poly::Polynomial;
use ark_std::rand::{rngs::StdRng, SeedableRng};
use ark_std::UniformRand;
use crate::{bind_variable, coeffs_to_univariate_poly, GmulyStatement};

use super::prover::GmulyProof;

/// Gmuly协议的Verifier
pub struct GmulyVerifier<F: Field> {
    /// statement
    statement: GmulyStatement<F>,

    /// fold times
    fold: usize,
}

impl<F: Field> GmulyVerifier<F> {
    /// 创建新的Gmuly Verifier
    pub fn new(statement: GmulyStatement<F>, fold: usize) -> Self {
        Self {
            statement,
            fold,
        }
    }

    /// 验证gmuly证明
    /// 
    /// # Arguments
    /// * `proof` - gmuly证明
    /// 
    /// # Returns
    /// 如果证明有效返回true，否则返回false
    pub fn verify(
        &self,
        proof: &GmulyProof<F>,
    ) -> bool {
        // 多项式个数应该等于折叠次数
        if proof.round_polynomials.len() != self.fold {
            return false;
        }

        let mut target = self.statement.sum_value;
        let seed = 42u64;
        let mut rng = StdRng::seed_from_u64(seed);

        // sum check
        // sum_value = g_0(0) + g_0(1)
        let mut current_poly = coeffs_to_univariate_poly(&proof.round_polynomials[0]);
        let mut sum_check = current_poly.evaluate(&vec![F::ZERO]) + current_poly.evaluate(&vec![F::ONE]);
        if sum_check != target {
            return false;
        }

        for i in 1..self.fold - 1 {
            let last_poly = current_poly;
            current_poly = coeffs_to_univariate_poly(&proof.round_polynomials[i]);

            let rand_element = F::rand(&mut rng);
            // i round check the validation of g_{i-1}
            // 检查g_i(r_i) = g_{i+1}(0) + g_{i+1}(1)
            sum_check = current_poly.evaluate(&vec![F::ZERO]) + current_poly.evaluate(&vec![F::ONE]);
            target = last_poly.evaluate(&vec![rand_element]);
            if sum_check != target {
                return false;
            }
        }

        // 检查g_fold
        true
    }
}

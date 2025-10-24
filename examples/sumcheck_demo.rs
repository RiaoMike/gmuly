use gmuly::{SumcheckProver, SumcheckVerifier};
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial,
};
use ark_test_curves::bls12_381::Fr;

fn main() {
    println!("=== Sumcheck Protocol Demo ===\n");

    // 示例1: 简单的线性多项式
    example1_linear_polynomial();
    
    // 示例2: 乘积多项式
    example2_product_polynomial();
    
    // 示例3: 更复杂的多项式
    example3_complex_polynomial();
}

/// 示例1: 线性多项式 f(x1, x2) = x1 + x2
fn example1_linear_polynomial() {
    println!("📝 Example 1: Linear Polynomial f(x1, x2) = x1 + x2");
    println!("---------------------------------------------------");
    
    // 构造多项式 f(x1, x2) = x1 + x2
    let mut terms = Vec::new();
    terms.push((Fr::from(1u64), SparseTerm::new(vec![(0, 1)]))); // x1 (变量0的1次方)
    terms.push((Fr::from(1u64), SparseTerm::new(vec![(1, 1)]))); // x2 (变量1的1次方)
    
    let poly = SparsePolynomial::from_coefficients_vec(2, terms);
    
    // 计算期望的和
    // f(0,0) + f(0,1) + f(1,0) + f(1,1) = 0 + 1 + 1 + 2 = 4
    let claimed_sum = Fr::from(4u64);
    
    println!("多项式: f(x1, x2) = x1 + x2");
    println!("声称的和: {}", claimed_sum);
    println!("实际计算:");
    println!("  f(0,0) = 0");
    println!("  f(0,1) = 1");
    println!("  f(1,0) = 1");
    println!("  f(1,1) = 2");
    println!("  总和 = 4\n");
    
    // Prover生成证明
    let prover = SumcheckProver::new(poly, 2, claimed_sum);
    let (proof, challenges, final_eval) = prover.prove();
    
    println!("✅ Prover生成了证明");
    println!("随机挑战: {:?}", challenges);
    println!("最终评估值: {}\n", final_eval);
    
    // Verifier验证证明
    let verifier = SumcheckVerifier::new(2, claimed_sum);
    let is_valid = verifier.verify(&proof, &challenges, final_eval);
    
    if is_valid {
        println!("✓ 验证成功! 证明有效。\n");
    } else {
        println!("✗ 验证失败! 证明无效。\n");
    }
}

/// 示例2: 乘积多项式 f(x1, x2) = x1 * x2
fn example2_product_polynomial() {
    println!("📝 Example 2: Product Polynomial f(x1, x2) = x1 * x2");
    println!("----------------------------------------------------");
    
    // 构造多项式 f(x1, x2) = x1 * x2
    let terms = vec![(Fr::from(1u64), SparseTerm::new(vec![(0, 1), (1, 1)]))];
    let poly = SparsePolynomial::from_coefficients_vec(2, terms);
    
    // 计算期望的和
    // f(0,0) + f(0,1) + f(1,0) + f(1,1) = 0 + 0 + 0 + 1 = 1
    let claimed_sum = Fr::from(1u64);
    
    println!("多项式: f(x1, x2) = x1 * x2");
    println!("声称的和: {}", claimed_sum);
    println!("实际计算:");
    println!("  f(0,0) = 0");
    println!("  f(0,1) = 0");
    println!("  f(1,0) = 0");
    println!("  f(1,1) = 1");
    println!("  总和 = 1\n");
    
    let prover = SumcheckProver::new(poly, 2, claimed_sum);
    let (proof, challenges, final_eval) = prover.prove();
    
    println!("✅ Prover生成了证明");
    println!("随机挑战: {:?}", challenges);
    println!("最终评估值: {}\n", final_eval);
    
    let verifier = SumcheckVerifier::new(2, claimed_sum);
    let is_valid = verifier.verify(&proof, &challenges, final_eval);
    
    if is_valid {
        println!("✓ 验证成功! 证明有效。\n");
    } else {
        println!("✗ 验证失败! 证明无效。\n");
    }
}

/// 示例3: 更复杂的多项式 f(x1, x2, x3) = 2*x1 + 3*x2 + x1*x2 + x2*x3
fn example3_complex_polynomial() {
    println!("📝 Example 3: Complex Polynomial f(x1, x2, x3) = 2*x1 + 3*x2 + x1*x2 + x2*x3");
    println!("-------------------------------------------------------------------------------");
    
    // 构造多项式
    let mut terms = Vec::new();
    terms.push((Fr::from(2u64), SparseTerm::new(vec![(0, 1)]))); // 2*x1
    terms.push((Fr::from(3u64), SparseTerm::new(vec![(1, 1)]))); // 3*x2
    terms.push((Fr::from(1u64), SparseTerm::new(vec![(0, 1), (1, 1)]))); // x1*x2
    terms.push((Fr::from(1u64), SparseTerm::new(vec![(1, 1), (2, 1)]))); // x2*x3
    
    let poly = SparsePolynomial::from_coefficients_vec(3, terms);
    
    // 计算期望的和（需要遍历所有8个布尔值组合）
    // f(0,0,0)=0, f(1,0,0)=2, f(0,1,0)=3, f(1,1,0)=6
    // f(0,0,1)=0, f(1,0,1)=2, f(0,1,1)=4, f(1,1,1)=7
    // 总和 = 0+2+3+6+0+2+4+7 = 24
    let claimed_sum = Fr::from(24u64);
    
    println!("多项式: f(x1, x2, x3) = 2*x1 + 3*x2 + x1*x2 + x2*x3");
    println!("声称的和: {}", claimed_sum);
    println!("实际计算:");
    println!("  f(0,0,0) = 0");
    println!("  f(1,0,0) = 2");
    println!("  f(0,1,0) = 3");
    println!("  f(1,1,0) = 6");
    println!("  f(0,0,1) = 0");
    println!("  f(1,0,1) = 2");
    println!("  f(0,1,1) = 4");
    println!("  f(1,1,1) = 7");
    println!("  总和 = 24\n");
    
    let prover = SumcheckProver::new(poly, 3, claimed_sum);
    let (proof, challenges, final_eval) = prover.prove();
    
    println!("✅ Prover生成了证明");
    println!("随机挑战: {:?}", challenges);
    println!("最终评估值: {}\n", final_eval);
    
    let verifier = SumcheckVerifier::new(3, claimed_sum);
    let is_valid = verifier.verify(&proof, &challenges, final_eval);
    
    if is_valid {
        println!("✓ 验证成功! 证明有效。\n");
    } else {
        println!("✗ 验证失败! 证明无效。\n");
    }
}

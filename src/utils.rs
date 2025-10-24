/// 通用工具函数模块
/// 
/// 提供稀疏多项式操作的通用函数，用于GMulY协议
use ark_ff::Field;
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial, Polynomial,
};

/// 对变量 x_k, x_{k+1}, ..., x_n 在布尔域 {0,1} 上进行部分求和
/// 
/// 返回一个新的稀疏多项式，只包含变量 x_0, ..., x_{k-1}
/// 
/// # 优化版本
/// 利用稀疏多项式特性和布尔域性质
/// 时间复杂度：O(项数) 而不是 O(2^n)
/// 
/// # 数学原理
/// 对于布尔变量 x ∈ {0,1}:
/// - sum_{x∈{0,1}} x^0 = 2 (不参与该项)
/// - sum_{x∈{0,1}} x^p = 1 (当 p > 0，参与该项)
/// 
/// 对于每一项 coeff * x_0^{p_0} * ... * x_n^{p_n}:
/// ```text
/// sum_{x_k,...,x_n ∈ {0,1}} coeff * x_0^{p_0} * ... * x_n^{p_n}
/// = coeff * x_0^{p_0} * ... * x_{k-1}^{p_{k-1}} * 2^(索引>=k且p=0)
/// ```
/// 
/// # Arguments
/// * `poly` - 输入的多变量多项式
/// * `start_index` - 求和起始索引 k，对 x_k, x_{k+1}, ..., x_n 求和
/// 
/// # Returns
/// 新的稀疏多项式，只包含 x_0, ..., x_{k-1}
/// 
/// # Example
/// ```rust,ignore
/// use gmuly::utils::partialsum;
/// 
/// // 对于多项式 f(x0, x1, x2) = x0 + x1 + x2
/// // partialsum(&poly, 1) 对 x1, x2 求和
/// // 结果: 一个包含 x0 的多项式
/// ```
pub fn partialsum<F: Field>(
    poly: &SparsePolynomial<F, SparseTerm>,
    start_index: usize,
) -> SparsePolynomial<F, SparseTerm> {
    let num_vars = poly.num_vars;

    if start_index >= num_vars {
        // 没有变量需要求和，返回原多项式的克隆
        return poly.clone();
    }

    let sum_vars_count = num_vars - start_index;
    let mut new_terms = Vec::new();

    // 直接从稀疏多项式的项计算，复杂度 O(项数)
    for (coeff, term) in poly.terms.iter() {
        // 构建新项：只保留索引 < start_index 的变量
        let mut new_term_vec = Vec::new();
        let mut participating_vars = 0;

        for (var, power) in term.iter() {
            if *var < start_index {
                // 保留这个变量到新项中
                new_term_vec.push((*var, *power));
            } else if *power > 0 {
                // 索引 >= start_index 且幂次 > 0，这是参与求和的变量
                participating_vars += 1;
            }
        }

        // 不参与的变量数量
        let non_participating = sum_vars_count - participating_vars;

        // 新系数 = 原系数 * 2^non_participating
        let new_coeff = *coeff * F::from((1u64 << non_participating) as u64);

        if !new_coeff.is_zero() {
            let new_term = SparseTerm::new(new_term_vec);
            new_terms.push((new_coeff, new_term));
        }
    }

    // 新多项式的变量数是 start_index
    SparsePolynomial::from_coefficients_vec(start_index, new_terms)
}

/// 对多项式的所有变量在布尔域 {0,1} 上求和，返回标量值
/// 
/// 这是 `partialsum` 的特殊情况，对所有变量求和并返回结果的常数值
/// 
/// # Arguments
/// * `poly` - 输入的多变量多项式
/// 
/// # Returns
/// 对所有变量求和后的标量值
/// 
/// # Example
/// ```rust,ignore
/// use gmuly::utils::all_sum;
/// 
/// // 对于多项式 f(x0, x1) = x0 + x1 + 1
/// // f(0,0)=1, f(0,1)=2, f(1,0)=2, f(1,1)=3
/// // all_sum(&poly) = 1 + 2 + 2 + 3 = 8
/// ```
pub fn all_sum<F: Field>(poly: &SparsePolynomial<F, SparseTerm>) -> F {
    // 对所有变量求和，得到一个常数多项式（num_vars = 0）
    let constant_poly = partialsum(poly, 0);
    
    // 对常数多项式求值（空向量），得到标量值
    constant_poly.evaluate(&vec![])
}

/// 用值绑定多项式的某个变量
/// # Arguments
/// * `poly` - 当前的多变量多项式
/// * `var_index` - 要绑定的变量索引
/// * `value` - 绑定的值
/// 
/// # Returns
/// 新的多项式，变量数为 `poly.num_vars - 1`，索引大于 `var_index` 的变量索引会减 1
/// 
/// # Example
/// ```rust,ignore
/// use gmuly::utils::bind_variable;
/// 
/// // 对于多项式 f(x0, x1) = x0 * x1
/// // bind_variable(&poly, 0, value) 固定 x0 = value
/// // 结果: g(x1) = value * x1
/// ```
pub fn bind_variable<F: Field>(
    poly: &SparsePolynomial<F, SparseTerm>,
    var_index: usize,
    value: F,
) -> SparsePolynomial<F, SparseTerm> {
    let num_vars = poly.num_vars;
    assert!(
        var_index < num_vars,
        "var_index must be less than number of variables"
    );

    let mut new_terms = Vec::new();

    for (coeff, term) in poly.terms.iter() {
        // 获取当前项中var_index变量的次数
        let var_power = term
            .iter()
            .find(|(var, _)| *var == var_index)
            .map(|(_, power)| *power)
            .unwrap_or(0);

        // 计算value^var_power
        let mut value_power = F::ONE;
        for _ in 0..var_power {
            value_power *= value;
        }

        // 创建新项，去掉var_index变量
        let new_term_vec: Vec<(usize, usize)> = term
            .iter()
            .filter(|(var, _)| *var != var_index)
            .map(|(var, power)| {
                // 重新映射变量索引
                let new_var = if *var > var_index { *var - 1 } else { *var };
                (new_var, *power)
            })
            .collect();

        let new_coeff = *coeff * value_power;
        if !new_coeff.is_zero() {
            let new_term = SparseTerm::new(new_term_vec);
            new_terms.push((new_coeff, new_term));
        }
    }

    SparsePolynomial::from_coefficients_vec(num_vars - 1, new_terms)
}

/// 将单变量稀疏多项式转换为系数向量
/// 
/// 将一个只包含一个变量的 SparsePolynomial 转换为 Vec<F>，
/// 其中 Vec[i] 表示 x^i 的系数，系数按照幂次从低到高排列。
/// 
/// # Arguments
/// * `poly` - 单变量稀疏多项式
/// 
/// # Returns
/// 系数向量，Vec[i] 对应 x^i 的系数
/// 
/// # Panics
/// 如果多项式包含多个变量或者变量索引不为 0
/// 
/// # Example
/// 对于多项式 3 + 2x + 5x^2，返回 vec![3, 2, 5]
pub fn univariate_poly_to_coeffs<F: Field>(
    poly: &SparsePolynomial<F, SparseTerm>,
) -> Vec<F> {
    assert_eq!(
        poly.num_vars, 1,
        "Polynomial must be univariate (have exactly 1 variable)"
    );

    // 如果没有项，返回空向量
    if poly.terms.is_empty() {
        return vec![F::ZERO];
    }

    // 辅助函数：获取项的次数
    let get_degree = |term: &SparseTerm| -> usize {
        if term.is_constant() {
            0
        } else {
            term.iter()
                .find(|(var, _)| *var == 0)
                .map(|(_, power)| *power)
                .unwrap_or(0)
        }
    };

    // 一次遍历：同时找到最高次数并收集所有 (次数, 系数) 对
    let mut max_degree = 0;
    let mut degree_coeff_pairs = Vec::with_capacity(poly.terms.len());
    
    for (coeff, term) in poly.terms.iter() {
        let degree = get_degree(term);
        max_degree = std::cmp::max(max_degree, degree);
        degree_coeff_pairs.push((degree, *coeff));
    }
    
    // 初始化系数向量（全为零）
    let mut coeffs = vec![F::ZERO; max_degree + 1];
    
    // 填充系数
    for (degree, coeff) in degree_coeff_pairs {
        coeffs[degree] += coeff;
    }
    
    coeffs
}

/// 将系数向量转换为单变量稀疏多项式
/// 
/// 将一个系数向量 Vec<F> 转换为只包含一个变量的 SparsePolynomial，
/// 其中 coeffs[i] 表示 x^i 的系数。
/// 
/// # Arguments
/// * `coeffs` - 系数向量，coeffs[i] 对应 x^i 的系数
/// 
/// # Returns
/// 单变量稀疏多项式（num_vars = 1）
/// 
/// # Example
/// 对于 vec![3, 2, 5]，返回多项式 3 + 2x + 5x^2
pub fn coeffs_to_univariate_poly<F: Field>(
    coeffs: &[F],
) -> SparsePolynomial<F, SparseTerm> {
    let mut terms = Vec::new();
    
    for (degree, coeff) in coeffs.iter().enumerate() {
        if *coeff != F::ZERO {
            let term = if degree == 0 {
                SparseTerm::new(vec![]) // 常数项
            } else {
                SparseTerm::new(vec![(0, degree)]) // x^degree
            };
            terms.push((*coeff, term));
        }
    }
    
    SparsePolynomial::from_coefficients_vec(1, terms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Field;
    use ark_poly::Polynomial;
    use ark_test_curves::bls12_381::Fr;

    #[test]
    fn test_partialsum_linear() {
        // f(x0, x1, x2) = x0 + x1 + x2
        let mut terms = Vec::new();
        terms.push((Fr::from(1u64), SparseTerm::new(vec![(0, 1)]))); // x0
        terms.push((Fr::from(1u64), SparseTerm::new(vec![(1, 1)]))); // x1
        terms.push((Fr::from(1u64), SparseTerm::new(vec![(2, 1)]))); // x2
        let poly = SparsePolynomial::from_coefficients_vec(3, terms);

        // 对 x1, x2 求和
        let result = partialsum(&poly, 1);

        assert_eq!(result.num_vars, 1);
        assert_eq!(result.evaluate(&vec![Fr::ZERO]), Fr::from(4u64));
        assert_eq!(result.evaluate(&vec![Fr::ONE]), Fr::from(8u64));
    }

    #[test]
    fn test_bind_variable_simple() {
        // f(x0, x1) = x0 * x1
        let terms = vec![(Fr::from(1u64), SparseTerm::new(vec![(0, 1), (1, 1)]))];
        let poly = SparsePolynomial::from_coefficients_vec(2, terms);

        // 绑定 x0 = 2
        let result = bind_variable(&poly, 0, Fr::from(2u64));

        assert_eq!(result.num_vars, 1);
        // g(x1) = 2 * x1
        assert_eq!(result.evaluate(&vec![Fr::ZERO]), Fr::ZERO);
        assert_eq!(result.evaluate(&vec![Fr::ONE]), Fr::from(2u64));
    }

    #[test]
    fn test_univariate_poly_to_coeffs_simple() {
        // 测试多项式 3 + 2x + 5x^2
        let mut terms = Vec::new();
        terms.push((Fr::from(3u64), SparseTerm::new(vec![]))); // 常数项 3
        terms.push((Fr::from(2u64), SparseTerm::new(vec![(0, 1)]))); // 2x
        terms.push((Fr::from(5u64), SparseTerm::new(vec![(0, 2)]))); // 5x^2
        let poly = SparsePolynomial::from_coefficients_vec(1, terms);

        let coeffs = univariate_poly_to_coeffs(&poly);

        assert_eq!(coeffs.len(), 3);
        assert_eq!(coeffs[0], Fr::from(3u64)); // x^0 系数
        assert_eq!(coeffs[1], Fr::from(2u64)); // x^1 系数
        assert_eq!(coeffs[2], Fr::from(5u64)); // x^2 系数
    }

    #[test]
    fn test_univariate_poly_to_coeffs_with_gaps() {
        // 测试多项式 1 + 0x + 0x^2 + 4x^3 (有跳过的次数)
        let mut terms = Vec::new();
        terms.push((Fr::from(1u64), SparseTerm::new(vec![]))); // 常数项 1
        terms.push((Fr::from(4u64), SparseTerm::new(vec![(0, 3)]))); // 4x^3
        let poly = SparsePolynomial::from_coefficients_vec(1, terms);

        let coeffs = univariate_poly_to_coeffs(&poly);

        assert_eq!(coeffs.len(), 4);
        assert_eq!(coeffs[0], Fr::from(1u64)); // x^0 系数
        assert_eq!(coeffs[1], Fr::ZERO); // x^1 系数
        assert_eq!(coeffs[2], Fr::ZERO); // x^2 系数
        assert_eq!(coeffs[3], Fr::from(4u64)); // x^3 系数
    }

    #[test]
    fn test_univariate_poly_to_coeffs_empty() {
        // 测试空多项式
        let poly: SparsePolynomial<Fr, SparseTerm> = SparsePolynomial::from_coefficients_vec(1, vec![]);
        let coeffs = univariate_poly_to_coeffs(&poly);
        
        assert_eq!(coeffs.len(), 1);
        assert_eq!(coeffs[0], Fr::ZERO);
    }

    #[test]
    fn test_coeffs_to_univariate_poly_simple() {
        // 测试 vec![3, 2, 5] -> 3 + 2x + 5x^2
        let coeffs = vec![Fr::from(3u64), Fr::from(2u64), Fr::from(5u64)];
        let poly = coeffs_to_univariate_poly(&coeffs);

        assert_eq!(poly.num_vars, 1);
        // 验证在不同点的求值
        assert_eq!(poly.evaluate(&vec![Fr::ZERO]), Fr::from(3u64)); // 3
        assert_eq!(poly.evaluate(&vec![Fr::ONE]), Fr::from(10u64)); // 3 + 2 + 5 = 10
        assert_eq!(poly.evaluate(&vec![Fr::from(2u64)]), Fr::from(27u64)); // 3 + 4 + 20 = 27
    }

    #[test]
    fn test_coeffs_to_univariate_poly_with_zeros() {
        // 测试 vec![1, 0, 0, 4] -> 1 + 4x^3
        let coeffs = vec![Fr::from(1u64), Fr::ZERO, Fr::ZERO, Fr::from(4u64)];
        let poly = coeffs_to_univariate_poly(&coeffs);

        assert_eq!(poly.num_vars, 1);
        // 只应该有 2 项（非零系数）
        assert_eq!(poly.terms.len(), 2);
        
        // 验证在不同点的求值
        assert_eq!(poly.evaluate(&vec![Fr::ZERO]), Fr::from(1u64)); // 1
        assert_eq!(poly.evaluate(&vec![Fr::ONE]), Fr::from(5u64)); // 1 + 4 = 5
        assert_eq!(poly.evaluate(&vec![Fr::from(2u64)]), Fr::from(33u64)); // 1 + 4*8 = 33
    }

    #[test]
    fn test_round_trip_conversion() {
        // 测试往返转换：poly -> coeffs -> poly
        let original_coeffs = vec![
            Fr::from(7u64),
            Fr::from(3u64),
            Fr::ZERO,
            Fr::from(11u64),
        ];
        
        let poly = coeffs_to_univariate_poly(&original_coeffs);
        let coeffs = univariate_poly_to_coeffs(&poly);
        
        // 系数应该相同
        assert_eq!(coeffs.len(), original_coeffs.len());
        for (a, b) in coeffs.iter().zip(original_coeffs.iter()) {
            assert_eq!(*a, *b);
        }
        
        // 再次转换回多项式并验证求值
        let poly2 = coeffs_to_univariate_poly(&coeffs);
        for i in 0..5 {
            let x = Fr::from(i);
            assert_eq!(
                poly.evaluate(&vec![x]),
                poly2.evaluate(&vec![x])
            );
        }
    }

    #[test]
    fn test_all_sum() {
        // 测试 f(x0, x1) = x0 + x1 + 1
        // f(0,0)=1, f(0,1)=2, f(1,0)=2, f(1,1)=3
        // all_sum = 1 + 2 + 2 + 3 = 8
        let mut terms = Vec::new();
        terms.push((Fr::from(1u64), SparseTerm::new(vec![(0, 1)]))); // x0
        terms.push((Fr::from(1u64), SparseTerm::new(vec![(1, 1)]))); // x1
        terms.push((Fr::from(1u64), SparseTerm::new(vec![]))); // 常数1
        let poly = SparsePolynomial::from_coefficients_vec(2, terms);

        let result = all_sum(&poly);
        assert_eq!(result, Fr::from(8u64));
    }

    #[test]
    fn test_all_sum_constant() {
        // 测试常数多项式 f(x0, x1) = 5
        // f(0,0)=5, f(0,1)=5, f(1,0)=5, f(1,1)=5
        // all_sum = 5 + 5 + 5 + 5 = 20
        let terms = vec![(Fr::from(5u64), SparseTerm::new(vec![]))];
        let poly = SparsePolynomial::from_coefficients_vec(2, terms);

        let result = all_sum(&poly);
        assert_eq!(result, Fr::from(20u64));
    }

    #[test]
    fn test_all_sum_product() {
        // 测试 f(x0, x1) = x0 * x1
        // f(0,0)=0, f(0,1)=0, f(1,0)=0, f(1,1)=1
        // all_sum = 0 + 0 + 0 + 1 = 1
        let terms = vec![(Fr::from(1u64), SparseTerm::new(vec![(0, 1), (1, 1)]))];
        let poly = SparsePolynomial::from_coefficients_vec(2, terms);

        let result = all_sum(&poly);
        assert_eq!(result, Fr::from(1u64));
    }
}

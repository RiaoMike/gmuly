# Partialsum 函数使用说明

## 概述

`partialsum` 函数对稀疏多项式的部分变量在布尔域 {0,1} 上进行高效求和。

## 函数签名

```rust
pub fn partialsum(
    &self,
    poly: &SparsePolynomial<F, SparseTerm>,
    start_index: usize,
) -> SparsePolynomial<F, SparseTerm>
```

## 参数说明

- `poly`: 输入的多变量稀疏多项式
- `start_index`: 求和起始索引 k，对变量 x_k, x_{k+1}, ..., x_n 进行求和

## 返回值

返回一个新的稀疏多项式，只包含变量 x_0, x_1, ..., x_{k-1}

## 时间复杂度

**O(项数)** - 显著优于暴力枚举的 O(2^n)

## 数学原理

对于布尔变量 x ∈ {0,1}:
- sum_{x∈{0,1}} x^0 = 2 (变量不参与该项)
- sum_{x∈{0,1}} x^p = 1 (当 p > 0，变量参与该项)

因此对于每一项 `coeff * x_0^{p_0} * ... * x_n^{p_n}`:
```
sum_{x_k,...,x_n ∈ {0,1}} coeff * x_0^{p_0} * ... * x_n^{p_n}
= coeff * x_0^{p_0} * ... * x_{k-1}^{p_{k-1}} * 2^(不参与的变量数)
```

## 使用示例

### 示例 1: 线性多项式

```rust
use ark_test_curves::bls12_381::Fr;
use ark_poly::multivariate::{SparsePolynomial, SparseTerm};

// f(x0, x1, x2) = x0 + x1 + x2
let mut terms = Vec::new();
terms.push((Fr::from(1u64), SparseTerm::new(vec![(0, 1)]))); // x0
terms.push((Fr::from(1u64), SparseTerm::new(vec![(1, 1)]))); // x1  
terms.push((Fr::from(1u64), SparseTerm::new(vec![(2, 1)]))); // x2
let poly = SparsePolynomial::from_coefficients_vec(3, terms);

let prover = GmulyProver::new(poly.clone(), 3, Fr::ZERO);

// 对 x1, x2 求和
let result = prover.partialsum(&poly, 1);

// result(x0) = sum_{x1,x2∈{0,1}} (x0 + x1 + x2)
//            = 4*x0 + 4
// 因为:
// - x0 项: 系数 1 * 2^2 = 4
// - x1 项: 不包含 x0，求和后消失但贡献常数
// - x2 项: 不包含 x0，求和后消失但贡献常数
```

### 示例 2: 乘积多项式

```rust
// f(x0, x1) = x0 * x1
let terms = vec![(Fr::from(1u64), SparseTerm::new(vec![(0, 1), (1, 1)]))];
let poly = SparsePolynomial::from_coefficients_vec(2, terms);

let prover = GmulyProver::new(poly.clone(), 2, Fr::ZERO);

// 对 x1 求和
let result = prover.partialsum(&poly, 1);

// result(x0) = sum_{x1∈{0,1}} x0*x1 = x0*0 + x0*1 = x0
```

### 示例 3: 常数多项式

```rust
// f(x0, x1) = 5
let terms = vec![(Fr::from(5u64), SparseTerm::new(vec![]))];
let poly = SparsePolynomial::from_coefficients_vec(2, terms);

let prover = GmulyProver::new(poly.clone(), 2, Fr::ZERO);

// 对 x1 求和  
let result = prover.partialsum(&poly, 1);

// result(x0) = sum_{x1∈{0,1}} 5 = 5 + 5 = 10
```

### 示例 4: 对所有变量求和

```rust
// f(x0, x1) = x0 + x1 + 1
let mut terms = Vec::new();
terms.push((Fr::from(1u64), SparseTerm::new(vec![(0, 1)]))); // x0
terms.push((Fr::from(1u64), SparseTerm::new(vec![(1, 1)]))); // x1
terms.push((Fr::from(1u64), SparseTerm::new(vec![]))); // 常数1
let poly = SparsePolynomial::from_coefficients_vec(2, terms);

let prover = GmulyProver::new(poly.clone(), 2, Fr::ZERO);

// 对所有变量求和
let result = prover.partialsum(&poly, 0);

// result = sum_{x0,x1∈{0,1}} (x0 + x1 + 1)
//        = f(0,0) + f(0,1) + f(1,0) + f(1,1)
//        = 1 + 2 + 2 + 3 = 8
```

## 性能优势

对于 n 变量稀疏多项式，具有 m 项：
- **暴力枚举**: O(m * 2^n) - 遍历所有 2^n 种组合
- **Partialsum**: O(m) - 只遍历 m 项

当 `m << 2^n` 时（稀疏多项式的常见情况），性能提升非常显著。

例如，对于 n=20 的多项式：
- 暴力枚举: ~1,048,576 次操作
- Partialsum (假设 m=100): ~100 次操作
- **加速比**: ~10,000x

## 注意事项

1. 函数只适用于布尔域 {0,1} 上的求和
2. 返回的多项式变量数为 `start_index`
3. 如果 `start_index >= poly.num_vars`，返回原多项式的克隆
4. 复杂度与多项式的项数线性相关，与变量数无关

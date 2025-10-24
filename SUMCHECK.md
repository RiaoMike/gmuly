# Sumcheck协议模块

本模块实现了sumcheck协议，这是一个交互式证明系统，用于证明多变量多项式在布尔超立方体上的求和。

## 概述

Sumcheck协议允许Prover向Verifier证明：对于一个n变量多项式g(x₁, x₂, ..., xₙ)，其在布尔超立方体{0,1}ⁿ上的求和等于某个声称的值H：

```
H = Σ_{x∈{0,1}ⁿ} g(x)
```

协议的通信复杂度为O(n·d)，其中n是变量数量，d是多项式的最大次数。

## 主要组件

### 1. `SumcheckProof<F>`
存储sumcheck协议的证明数据
- `round_polynomials`: 每一轮的单变量多项式系数

### 2. `SumcheckProver<F>`
Prover端实现
- 创建: `SumcheckProver::new(polynomial, num_vars, claimed_sum)`
- 生成证明: `prove()` → `(SumcheckProof, Vec<F>, F)`

### 3. `SumcheckVerifier<F>`
Verifier端实现
- 创建: `SumcheckVerifier::new(num_vars, claimed_sum)`
- 验证证明: `verify(&proof, &challenges, final_eval)` → `bool`

## 协议流程

1. **初始化阶段**
   - Prover和Verifier都知道多项式g和声称的和H
   - Verifier知道n（变量数量）

2. **交互轮次**（共n轮）
   - 每一轮i：
     - Prover发送单变量多项式gᵢ(X)
     - Verifier检查gᵢ(0) + gᵢ(1) = 上一轮的评估值
     - Verifier发送随机挑战rᵢ
     - Prover用rᵢ绑定第i个变量

3. **最终验证**
   - Prover发送g(r₁, r₂, ..., rₙ)
   - Verifier可以独立计算或使用oracle验证

## 使用示例

### 基础用法

```rust
use gmuly::{SumcheckProver, SumcheckVerifier};
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial,
};
use ark_test_curves::bls12_381::Fr;

// 创建多项式 f(x1, x2) = x1 + x2
let mut terms = Vec::new();
terms.push((Fr::from(1u64), SparseTerm::new(vec![(0, 1)]))); // x1
terms.push((Fr::from(1u64), SparseTerm::new(vec![(1, 1)]))); // x2

let poly = SparsePolynomial::from_coefficients_vec(2, terms);

// 计算期望的和: f(0,0) + f(0,1) + f(1,0) + f(1,1) = 0+1+1+2 = 4
let claimed_sum = Fr::from(4u64);

// Prover生成证明
let prover = SumcheckProver::new(poly, 2, claimed_sum);
let (proof, challenges, final_eval) = prover.prove();

// Verifier验证证明
let verifier = SumcheckVerifier::new(2, claimed_sum);
let is_valid = verifier.verify(&proof, &challenges, final_eval);

assert!(is_valid);
```

### 多项式定义

使用arkworks的`SparsePolynomial`和`SparseTerm`：

```rust
// 常数多项式: f = 5
let terms = vec![(Fr::from(5u64), SparseTerm::new(vec![]))];

// 单项式: f = x1
let terms = vec![(Fr::from(1u64), SparseTerm::new(vec![(0, 1)]))];

// 乘积: f = x1 * x2
let terms = vec![(Fr::from(1u64), SparseTerm::new(vec![(0, 1), (1, 1)]))];

// 复杂多项式: f = 2*x1 + 3*x2 + x1*x2
let terms = vec![
    (Fr::from(2u64), SparseTerm::new(vec![(0, 1)])),
    (Fr::from(3u64), SparseTerm::new(vec![(1, 1)])),
    (Fr::from(1u64), SparseTerm::new(vec![(0, 1), (1, 1)])),
];
```

**注意**: `SparseTerm`使用`(变量索引, 次数)`的元组格式。例如：
- `(0, 1)` 表示变量0的1次方
- `(1, 2)` 表示变量1的2次方
- `[(0, 1), (1, 1)]` 表示 x₀ * x₁

## 运行示例

```bash
# 运行sumcheck示例
cargo run --example sumcheck_demo --features examples

# 运行测试
cargo test --lib sumcheck
```

## 技术细节

### 拉格朗日插值

本实现使用拉格朗日插值将多项式在点{0, 1, ..., d}的评估值转换为系数形式。对于度数为d的多项式p(x)，给定评估值p(0), p(1), ..., p(d)，可以恢复系数。

### 变量绑定

在每一轮中，Prover需要将当前变量绑定为Verifier提供的随机值。实现中总是绑定第0个变量，并重新映射剩余变量的索引。

### 安全性考虑

本实现中，随机挑战使用确定性生成（用于演示目的）。在实际应用中：
- 随机挑战应由Verifier生成并发送给Prover
- 或使用Fiat-Shamir变换将交互协议转换为非交互式

## 依赖

- `ark-ff`: 有限域运算
- `ark-poly`: 多变量多项式
- `ark-std`: 标准库抽象

## 参考资料

- [Sumcheck Protocol - Proofs, Arguments, and Zero-Knowledge (PAZK)](https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf)
- [The Unreasonable Power of the Sumcheck Protocol](https://people.cs.georgetown.edu/jthaler/blogpost.pdf)

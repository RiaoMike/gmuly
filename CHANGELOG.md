# Changelog

## [Latest] - 2025-01-20

### Added - Sumcheck Protocol Module

新增了完整的sumcheck协议实现，使用arkworks的多变量多项式模块。

#### 主要组件

1. **SumcheckProof<F>**
   - 存储每一轮的单变量多项式系数
   - 支持序列化和克隆

2. **SumcheckProver<F>**
   - 实现Prover端逻辑
   - 自动生成每轮的单变量多项式
   - 使用拉格朗日插值将评估值转换为系数
   - 支持任意次数的多变量多项式

3. **SumcheckVerifier<F>**
   - 实现Verifier端验证逻辑
   - 检查每轮多项式的一致性
   - 验证最终评估值

#### 技术特性

- **拉格朗日插值**: 实现了从评估值恢复多项式系数的算法
- **变量绑定**: 每轮自动绑定第0个变量并重新映射索引
- **多项式表示**: 使用arkworks的`SparsePolynomial`和`SparseTerm`
- **通信复杂度**: O(n·d)，其中n是变量数，d是最大次数

#### 测试覆盖

新增4个测试用例：
- `test_sumcheck_simple`: 线性多项式 f(x1, x2) = x1 + x2
- `test_sumcheck_constant`: 常数多项式 f(x1, x2) = 5
- `test_sumcheck_product`: 乘积多项式 f(x1, x2) = x1 * x2
- 通过MLE模块的`test_partialsum_for_sumcheck`测试集成

#### 示例

新增`sumcheck_demo`示例，展示：
1. 线性多项式的sumcheck证明
2. 乘积多项式的sumcheck证明
3. 三变量复杂多项式的sumcheck证明

#### 文档

- 新增`SUMCHECK.md`详细文档
- 更新主`README.md`添加sumcheck特性和示例
- 包含完整的API文档和使用指南

#### 依赖更新

- 新增`ark-poly = "0.4"`用于多变量多项式支持

### 当前状态

- ✅ 35个测试全部通过
- ✅ 3个工作示例（matrix_operations, mle_demo, sumcheck_demo）
- ✅ 完整的文档覆盖
- ⚠️ 3个编译警告（使用了已弃用的方法，计划在未来版本修复）

### 模块结构

```
gmuly/
├── src/
│   ├── lib.rs           - 模块导出
│   ├── sparse_matrix.rs - CSR稀疏矩阵
│   ├── ccs.rs          - CCS约束系统
│   ├── mle.rs          - 多线性扩展
│   └── sumcheck.rs     - Sumcheck协议 (NEW)
├── examples/
│   ├── matrix_operations.rs
│   ├── mle_demo.rs
│   └── sumcheck_demo.rs    (NEW)
├── README.md
├── SUMCHECK.md             (NEW)
└── CHANGELOG.md            (NEW)
```

### 下一步计划

可能的改进方向：
1. 实现Fiat-Shamir变换使协议非交互式
2. 优化拉格朗日插值性能
3. 添加更多多项式类型的测试
4. 与MLE模块的更深度集成
5. 实现基于sumcheck的应用（如GKR协议）

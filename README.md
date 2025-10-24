# Gmuly - CCS Constraint System with CSR Sparse Matrices

A Rust implementation of Customizable Constraint Systems (CCS) using Compressed Sparse Row (CSR) format sparse matrices, along with Multilinear Extensions (MLE) and the Sumcheck protocol.

## Features

- **CSR Sparse Matrix**: Efficient sparse matrix implementation using the standard CSR format
  - Memory-efficient storage with only 3 arrays: `values`, `col_indices`, `row_ptr`
  - Automatic row count derivation from `row_ptr` length (no redundant storage)
  
- **Operator Overloading**: Natural matrix operations using Rust's `Mul` trait
  - Matrix-vector multiplication: `&matrix * &vec`
  - Matrix-matrix multiplication: `&matrix_a * &matrix_b`
  
- **CCS Constraint System**: Generalized constraint system supporting R1CS, Plonkish, AIR, and more
  - Customizable constraint terms with constants and matrix multisets
  - Hadamard (element-wise) product support
  - Witness satisfaction checking

- **Multilinear Extension (MLE)**: Multilinear polynomials with coefficient representation
  - Stores coefficients in tensor product basis: `(1, x_1, x_2, x_1*x_2)` for 2 variables
  - Efficient conversion between coefficient and evaluation representations
  - Support for evaluation, variable binding, and arithmetic operations
  - Automatic padding to power-of-2 dimensions

- **Sumcheck Protocol**: Interactive proof system for polynomial evaluations
  - Prover and Verifier implementations
  - Support for multivariate polynomials using arkworks
  - Lagrange interpolation for coefficient recovery
  - O(n·d) communication complexity

## Quick Start

### Matrix-Vector Multiplication

```rust
use gmuly::SparseMatrix;
use ark_test_curves::bls12_381::Fr;

// Create a sparse matrix
let dense = vec![
    vec![Fr::from(1u64), Fr::from(2u64), Fr::from(0u64)],
    vec![Fr::from(0u64), Fr::from(3u64), Fr::from(4u64)],
];
let matrix = SparseMatrix::from_dense(&dense);

// Multiply with a vector using the Mul trait
let vec = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
let result = &matrix * &vec;  // [5, 18]
```

### Matrix-Matrix Multiplication

```rust
// Create two matrices
let matrix_a = SparseMatrix::from_dense(&dense_a);
let matrix_b = SparseMatrix::from_dense(&dense_b);

// Multiply matrices using the Mul trait
let result = &matrix_a * &matrix_b;
```

### CCS Constraint System

```rust
use gmuly::{CCS, SparseMatrix};
use ark_test_curves::bls12_381::Fr;

// Example: R1CS constraint x * y = z
// Expressed as: A·z ⊙ B·z - C·z = 0

// Define matrices (witness format: [1, x, y, z])
let mut matrix_a = SparseMatrix::new(1, 4);
matrix_a.set(0, 1, Fr::from(1u64)); // Select x

let mut matrix_b = SparseMatrix::new(1, 4);
matrix_b.set(0, 2, Fr::from(1u64)); // Select y

let mut matrix_c = SparseMatrix::new(1, 4);
matrix_c.set(0, 3, Fr::from(1u64)); // Select z

// Build CCS
let ccs = CCS::new_with(
    vec![Fr::from(1u64), -Fr::from(1u64)], // Constants
    vec![vec![0, 1], vec![2]],              // Multisets: [A⊙B, C]
    vec![matrix_a, matrix_b, matrix_c],
).unwrap();

// Check witness
let witness = vec![
    Fr::from(1u64),  // constant
    Fr::from(3u64),  // x
    Fr::from(4u64),  // y
    Fr::from(12u64), // z = x*y
];
assert!(ccs.is_satisfied(&witness));
```

## Architecture

### CSR Sparse Matrix Format

The CSR format stores a sparse matrix using three arrays:

- `values`: Non-zero values in row-major order
- `col_indices`: Column index for each non-zero value
- `row_ptr`: Row pointers (length = rows + 1)

**Optimization**: The number of rows is derived from `row_ptr.len() - 1`, eliminating redundant storage.

### CCS Structure

```rust
pub struct CCS<F: Field> {
    /// Constants for each constraint term
    pub constants: Vec<F>,
    /// Sets of matrix indices for Hadamard products
    pub multisets: Vec<Vec<usize>>,
    /// Constraint matrices in CSR format
    pub matrices: Vec<SparseMatrix<F>>,
}
```

Represents constraints: **∑ᵢ cᵢ · ⊙_{j ∈ Sᵢ} Mⱼ · z = 0**

## API Highlights

### SparseMatrix Methods

- `new(rows, cols)` - Create empty matrix
- `from_dense(data)` - Convert from dense format
- `from_csr(cols, values, col_indices, row_ptr)` - Build from CSR components
- `identity(size)` - Create identity matrix
- `transpose()` - Transpose the matrix
- `rows()` - Get number of rows (derived from `row_ptr`)
- `cols()` - Get number of columns
- `nnz()` - Get number of non-zero elements
- `get(row, col)` - Get element value
- `set(row, col, value)` - Set element value

### Mul Trait Implementations

- `&SparseMatrix<F> * &[F]` → `Vec<F>` (matrix-vector)
- `&SparseMatrix<F> * &Vec<F>` → `Vec<F>` (matrix-vector)
- `&SparseMatrix<F> * Vec<F>` → `Vec<F>` (matrix-vector, owned)
- `&SparseMatrix<F> * &SparseMatrix<F>` → `SparseMatrix<F>` (matrix-matrix)
- `&SparseMatrix<F> * SparseMatrix<F>` → `SparseMatrix<F>` (matrix-matrix, owned)

### CCS Methods

- `new()` - Create empty CCS
- `new_with(constants, multisets, matrices)` - Create with validation
- `is_satisfied(witness)` - Check if witness satisfies constraints
- `add_term(constant, multiset)` - Add constraint term
- `add_matrix(matrix)` - Add constraint matrix
- `num_terms()`, `num_matrices()`, `num_constraints()`, `num_variables()` - Query dimensions

### MLE Methods

- `new(coefficients)` - Create from coefficient vector in tensor product basis
- `from_evaluations(evaluations)` - Create from evaluations at boolean hypercube
- `from_vector(vec)` - Convert vector to MLE (treats as evaluations)
- `from_matrix(matrix, rows, cols)` - Convert matrix to MLE (treats as evaluations)
- `evaluate(point)` - Evaluate MLE at a point
- `bind(var_index, value)` - Bind variable to value
- `add(other)`, `multiply(other)`, `scale(scalar)` - Arithmetic operations
- `coefficients()` - Get coefficient vector
- `evaluations()` - Get evaluations at boolean hypercube
- `num_vars()` - Get number of variables

### Multilinear Extension (MLE)

```rust
use gmuly::MLE;
use ark_test_curves::bls12_381::Fr;

// Create MLE from coefficients
// For polynomial f(x_1, x_2) = 1 + 2*x_1 + 3*x_2 + 4*x_1*x_2
let coefficients = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64), Fr::from(4u64)];
let mle = MLE::new(coefficients);

// Create MLE from evaluations
let evaluations = vec![Fr::from(3u64), Fr::from(5u64), Fr::from(7u64), Fr::from(9u64)];
let mle2 = MLE::from_evaluations(&evaluations);

// Evaluate at a point
let point = vec![Fr::from(1u64), Fr::from(0u64)];
let value = mle2.evaluate(&point); // Returns 5

// Create MLE from matrix
let matrix = vec![
    vec![Fr::from(1u64), Fr::from(2u64)],
    vec![Fr::from(3u64), Fr::from(4u64)],
];
let mle_matrix = MLE::from_matrix(&matrix, 2, 2);

// Variable binding
let mle_bound = mle.bind(0, Fr::from(0u64)); // Bind first variable to 0

// MLE operations
let mle2 = MLE::from_vector(&vec![Fr::from(1u64); 4]);
let sum = mle.add(&mle2);
let product = mle.multiply(&mle2);
```

### Sumcheck Protocol

```rust
use gmuly::{SumcheckProver, SumcheckVerifier};
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial,
};
use ark_test_curves::bls12_381::Fr;

// Create polynomial f(x1, x2) = x1 + x2
let mut terms = Vec::new();
terms.push((Fr::from(1u64), SparseTerm::new(vec![(0, 1)]))); // x1
terms.push((Fr::from(1u64), SparseTerm::new(vec![(1, 1)]))); // x2

let poly = SparsePolynomial::from_coefficients_vec(2, terms);
let claimed_sum = Fr::from(4u64); // sum over {0,1}^2

// Prover generates proof
let prover = SumcheckProver::new(poly, 2, claimed_sum);
let (proof, challenges, final_eval) = prover.prove();

// Verifier verifies proof
let verifier = SumcheckVerifier::new(2, claimed_sum);
let is_valid = verifier.verify(&proof, &challenges, final_eval);
assert!(is_valid);
```

For detailed sumcheck documentation, see [SUMCHECK.md](SUMCHECK.md).

## Examples

Run the examples:

```bash
# Matrix operations demo
cargo run --example matrix_operations --features ark-test-curves

# MLE demo
cargo run --example mle_demo --features ark-test-curves

# Sumcheck protocol demo
cargo run --example sumcheck_demo --features examples
```

## Testing

Run all tests:

```bash
cargo test
```

Current test coverage:
- 27+ tests covering:
  - Sparse matrix operations, CSR format, matrix multiplication, transpose
  - CCS constraint satisfaction and validation
  - MLE creation from vectors/matrices, evaluation, binding, and arithmetic
  - Sumcheck protocol: prover/verifier with linear, product, and constant polynomials
- All tests passing ✓

## Dependencies

- `ark-ff` - Finite field arithmetic
- `ark-poly` - Multivariate polynomial support
- `ark-std` - Standard library abstractions
- `ark-test-curves` - For testing (optional)

## License

This project uses the same license as your Rust installation.


use gmuly::{MLE, SparseMatrix};
use ark_test_curves::bls12_381::Fr;
use ark_ff::Zero;

fn main() {
    println!("=== MLE and SparseMatrix Conversion Demo ===\n");

    // Example 1: Convert SparseMatrix to MLE and back
    example_sparse_to_mle();
    
    // Example 2: Square matrix conversion
    example_square_matrix();
}

fn example_sparse_to_mle() {
    println!("📊 Example 1: SparseMatrix → MLE → SparseMatrix");
    println!("------------------------------------------------");
    
    // Create a sparse matrix with some non-zero values (must be power-of-2 dimensions)
    let dense = vec![
        vec![Fr::from(1u64), Fr::from(0u64), Fr::from(3u64), Fr::from(0u64)],
        vec![Fr::from(0u64), Fr::from(5u64), Fr::from(0u64), Fr::from(7u64)],
        vec![Fr::from(9u64), Fr::from(0u64), Fr::from(11u64), Fr::from(0u64)],
        vec![Fr::from(0u64), Fr::from(13u64), Fr::from(0u64), Fr::from(15u64)],
    ];
    
    println!("Original matrix (4x4):");
    for (i, row) in dense.iter().enumerate() {
        print!("  Row {}: ", i);
        for val in row {
            if val.is_zero() {
                print!("0 ");
            } else {
                print!("{:?} ", val);
            }
        }
        println!();
    }
    
    // Convert to SparseMatrix
    let sparse = SparseMatrix::from_dense(&dense);
    println!("\nSparse matrix created:");
    println!("  Non-zero elements: {}", sparse.nnz());
    println!("  Dimensions: {}x{}", sparse.rows(), sparse.cols());
    
    // Convert to MLE (values are treated as coefficients, not evaluations)
    let mle = MLE::from_sparse_matrix(&sparse);
    println!("\nConverted to MLE:");
    println!("  Number of variables: {}", mle.num_vars());
    println!("  Split count: {}", mle.get_split_count());
    println!("  Coefficients length: {}", mle.coefficients().len());
    
    // Matrix values are directly mapped as coefficients
    println!("\nCoefficient mapping (row-major order):");
    println!("  First coefficient: {:?}", mle.coefficients()[0]);
    println!("  Last coefficient: {:?}", mle.coefficients()[15]);
    
    // Convert back to sparse matrix
    let sparse_back = mle.to_sparse_matrix();
    println!("\nConverted back to SparseMatrix ({}x{}):", sparse_back.rows(), sparse_back.cols());
    
    // Verify original values are preserved
    for i in 0..4 {
        print!("  Row {}: ", i);
        for j in 0..4 {
            let val = sparse_back.get(i, j);
            if val.is_zero() {
                print!("0 ");
            } else {
                print!("{:?} ", val);
            }
        }
        println!();
    }
    println!();
}

fn example_square_matrix() {
    println!("📊 Example 2: Square Matrix Conversion");
    println!("---------------------------------------");
    
    // Create a 2x2 matrix
    let dense = vec![
        vec![Fr::from(1u64), Fr::from(2u64)],
        vec![Fr::from(3u64), Fr::from(4u64)],
    ];
    
    println!("Original 2x2 matrix:");
    for (i, row) in dense.iter().enumerate() {
        print!("  Row {}: ", i);
        for val in row {
            print!("{:?} ", val);
        }
        println!();
    }
    
    // Convert to sparse and then to MLE
    let sparse = SparseMatrix::from_dense(&dense);
    let mle = MLE::from_sparse_matrix(&sparse);
    
    println!("\nMLE properties:");
    println!("  Number of variables: {}", mle.num_vars());
    println!("  Split count: {}", mle.get_split_count());
    
    // Convert to sparse matrix
    let square_sparse = mle.to_sparse_matrix();
    
    println!("\nSquare matrix conversion:");
    println!("  Dimensions: {}x{}", square_sparse.rows(), square_sparse.cols());
    
    // Display the square matrix (values should match original)
    for i in 0..square_sparse.rows() {
        print!("  Row {}: ", i);
        for j in 0..square_sparse.cols() {
            let val = square_sparse.get(i, j);
            if val.is_zero() {
                print!("0 ");
            } else {
                print!("{:?} ", val);
            }
        }
        println!();
    }
    
    // The MLE coefficients are directly the matrix values
    println!("\nMLE coefficients (in row-major order):");
    println!("  Coefficient 0 (matrix[0][0]): {:?}", mle.coefficients()[0]);
    println!("  Coefficient 1 (matrix[0][1]): {:?}", mle.coefficients()[1]);
    println!("  Coefficient 2 (matrix[1][0]): {:?}", mle.coefficients()[2]);
    println!("  Coefficient 3 (matrix[1][1]): {:?}", mle.coefficients()[3]);
    
    println!();
}

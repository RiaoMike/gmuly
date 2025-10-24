use gmuly::SparseMatrix;
use ark_test_curves::bls12_381::Fr;

fn main() {
    println!("=== Sparse Matrix Operations Demo ===\n");

    // Create a sparse matrix from dense data
    let dense = vec![
        vec![Fr::from(1u64), Fr::from(2u64), Fr::from(0u64)],
        vec![Fr::from(0u64), Fr::from(3u64), Fr::from(4u64)],
        vec![Fr::from(5u64), Fr::from(0u64), Fr::from(6u64)],
    ];
    let matrix = SparseMatrix::from_dense(&dense);
    
    println!("Matrix A (3x3):");
    println!("[1, 2, 0]");
    println!("[0, 3, 4]");
    println!("[5, 0, 6]\n");

    // Matrix-Vector Multiplication using Mul trait
    println!("--- Matrix-Vector Multiplication ---");
    let vec = vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)];
    println!("Vector: [1, 2, 3]");
    
    // Using the new Mul trait
    let result = &matrix * &vec;
    println!("A * v = [{}, {}, {}]", result[0], result[1], result[2]);
    println!("Expected: [5, 18, 23]\n");

    // Matrix-Matrix Multiplication using Mul trait
    println!("--- Matrix-Matrix Multiplication ---");
    let dense_b = vec![
        vec![Fr::from(1u64), Fr::from(0u64)],
        vec![Fr::from(0u64), Fr::from(1u64)],
        vec![Fr::from(1u64), Fr::from(1u64)],
    ];
    let matrix_b = SparseMatrix::from_dense(&dense_b);
    
    println!("Matrix B (3x2):");
    println!("[1, 0]");
    println!("[0, 1]");
    println!("[1, 1]\n");
    
    // Using the new Mul trait for matrix multiplication
    let result_matrix = &matrix * &matrix_b;
    println!("Matrix C = A * B (3x2):");
    for i in 0..result_matrix.rows() {
        print!("[");
        for j in 0..result_matrix.cols() {
            print!("{}", result_matrix.get(i, j));
            if j < result_matrix.cols() - 1 {
                print!(", ");
            }
        }
        println!("]");
    }
    println!();

    // Transpose
    println!("--- Matrix Transpose ---");
    let transposed = matrix.transpose();
    println!("A^T (3x3):");
    for i in 0..transposed.rows() {
        print!("[");
        for j in 0..transposed.cols() {
            print!("{}", transposed.get(i, j));
            if j < transposed.cols() - 1 {
                print!(", ");
            }
        }
        println!("]");
    }
    println!();

    // Identity matrix multiplication
    println!("--- Identity Matrix ---");
    let identity = SparseMatrix::<Fr>::identity(3);
    let result_identity = &matrix * &identity;
    println!("A * I = A (verified)");
    assert_eq!(result_identity.get(0, 0), Fr::from(1u64));
    assert_eq!(result_identity.get(1, 1), Fr::from(3u64));
    assert_eq!(result_identity.get(2, 2), Fr::from(6u64));
    println!("All identity tests passed!\n");

    println!("=== Demo Complete ===");
}


use gmuly::MLE;
use ark_test_curves::bls12_381::Fr;
use ark_ff::Field;

fn main() {
    println!("=== Multilinear Extension (MLE) Demo ===\n");

    // Example 1: MLE from Vector
    println!("--- 1. MLE from Vector ---");
    let vec = vec![
        Fr::from(3u64),
        Fr::from(5u64),
        Fr::from(7u64),
        Fr::from(9u64),
    ];
    
    println!("Original vector: [3, 5, 7, 9]");
    let mle_vec = MLE::from_vector(&vec);
    println!("MLE has {} variables", mle_vec.num_vars());
    println!("Coefficient ordering: (1, x_2, x_1, x_1*x_2) = (1, x_1) ⊗ (1, x_2)\n");
    
    // Verify at boolean points
    println!("Evaluation at boolean points:");
    println!("  f(0,0) = {}", mle_vec.evaluate(&[Fr::from(0u64), Fr::from(0u64)]));
    println!("  f(1,0) = {}", mle_vec.evaluate(&[Fr::from(1u64), Fr::from(0u64)]));
    println!("  f(0,1) = {}", mle_vec.evaluate(&[Fr::from(0u64), Fr::from(1u64)]));
    println!("  f(1,1) = {}", mle_vec.evaluate(&[Fr::from(1u64), Fr::from(1u64)]));
    
    // Evaluate at intermediate point
    let half = Fr::from(2u64).inverse().unwrap();
    println!("\nEvaluation at (1/2, 1/2) = {}", 
            mle_vec.evaluate(&[half, half]));
    println!();

    // Example 2: MLE from Matrix
    println!("--- 2. MLE from Matrix ---");
    let matrix = vec![
        vec![Fr::from(1u64), Fr::from(2u64), Fr::from(3u64)],
        vec![Fr::from(4u64), Fr::from(5u64), Fr::from(6u64)],
    ];
    
    println!("Original 2x3 matrix:");
    println!("  [1, 2, 3]");
    println!("  [4, 5, 6]");
    
    let mle_matrix = MLE::from_matrix(&matrix, 2, 3);
    println!("\nFlattened (row-major): [1, 2, 3, 4, 5, 6]");
    println!("Padded to power of 2: [1, 2, 3, 4, 5, 6, 0, 0]");
    println!("MLE has {} variables (2^3 = 8 points)", mle_matrix.num_vars());
    println!();

    // Example 3: Variable Binding
    println!("--- 3. Variable Binding ---");
    let evaluations = vec![
        Fr::from(1u64), Fr::from(2u64),
        Fr::from(3u64), Fr::from(4u64),
        Fr::from(5u64), Fr::from(6u64),
        Fr::from(7u64), Fr::from(8u64),
    ];
    
    let mle_3var = MLE::new(evaluations);
    println!("3-variable MLE with evaluations [1, 2, 3, 4, 5, 6, 7, 8]");
    
    // Bind first variable to 0
    let mle_bound_0 = mle_3var.bind(0, Fr::from(0u64));
    println!("\nAfter binding x_1 = 0:");
    println!("  Resulting 2-variable MLE: {:?}", mle_bound_0.evaluations());
    println!("  (These are evaluations at (0, x_2, x_3))");
    
    // Bind first variable to 1
    let mle_bound_1 = mle_3var.bind(0, Fr::from(1u64));
    println!("\nAfter binding x_1 = 1:");
    println!("  Resulting 2-variable MLE: {:?}", mle_bound_1.evaluations());
    println!("  (These are evaluations at (1, x_2, x_3))");
    
    // Bind to intermediate value
    let mle_bound_half = mle_3var.bind(0, half);
    println!("\nAfter binding x_1 = 1/2:");
    println!("  Resulting 2-variable MLE: {:?}", mle_bound_half.evaluations());
    println!();

    // Example 4: MLE Operations
    println!("--- 4. MLE Arithmetic Operations ---");
    let mle1 = MLE::new(vec![
        Fr::from(1u64), Fr::from(2u64),
        Fr::from(3u64), Fr::from(4u64),
    ]);
    
    let mle2 = MLE::new(vec![
        Fr::from(5u64), Fr::from(6u64),
        Fr::from(7u64), Fr::from(8u64),
    ]);
    
    println!("MLE 1: [1, 2, 3, 4]");
    println!("MLE 2: [5, 6, 7, 8]");
    
    let sum = mle1.add(&mle2);
    println!("\nAddition: {:?}", sum.evaluations());
    
    let product = mle1.multiply(&mle2);
    println!("Multiplication: {:?}", product.evaluations());
    
    let scaled = mle1.scale(Fr::from(3u64));
    println!("MLE 1 scaled by 3: {:?}", scaled.evaluations());
    println!();

    // Example 5: Understanding the coefficient ordering
    println!("--- 5. Coefficient Ordering Explanation ---");
    println!("For n variables, coefficients follow binary index pattern:");
    println!("2 vars: index binary -> basis");
    println!("  0:    00    -> 1");
    println!("  1:    01    -> x_1");
    println!("  2:    10    -> x_2");
    println!("  3:    11    -> x_1*x_2");
    println!("\nThis follows the tensor product: (1, x_1) ⊗ (1, x_2)");
    println!();

    println!("3 vars: index binary -> basis");
    println!("  0:    000   -> 1");
    println!("  1:    001   -> x_1");
    println!("  2:    010   -> x_2");
    println!("  3:    011   -> x_1*x_2");
    println!("  4:    100   -> x_3");
    println!("  5:    101   -> x_1*x_3");
    println!("  6:    110   -> x_2*x_3");
    println!("  7:    111   -> x_1*x_2*x_3");
    println!("\nThis follows: (1, x_1) ⊗ (1, x_2) ⊗ (1, x_3)");
    
    println!("\n=== Demo Complete ===");
}




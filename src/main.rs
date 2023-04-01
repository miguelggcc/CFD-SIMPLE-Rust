mod plotter;

use std::time::Instant;

use cfd_rust::Cases;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let now = Instant::now();

    let problem = //"backward_facing_step";
    "pipe_flow";
    //"lid_driven_cavity";

    let mut case = Cases::new(problem, 400, 80, 10.0);

    let iterations = 1500;
    let mut iter = 0;

    for i in 1..iterations + 1 {
        iter = i;
        if !case.iterate() {
            println!("Solution diverged!");
            break;
        }

        println!("Iteration {} / {}", i, iterations);
    }

    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);

    case.postprocessing(iter);

    Ok(())
}

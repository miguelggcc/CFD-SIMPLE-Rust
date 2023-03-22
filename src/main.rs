mod plotter;

use std::time::Instant;

use CFDplotlib::Cases;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let now = Instant::now();

    let problem = "lid_driven_cavity";

    let mut case = Cases::new(problem, 80, 80, 100.0);

    let iterations = 800;
    let mut iter = 0;

    for i in 1..iterations + 1 {
        iter = i;
        if !case.iterate() {
            println!("Solution diverged!");
            break;
        }

        println!("iteration {} / {}", i, iterations);
    }

    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);

    //case.postprocessing(iter);

    Ok(())
}

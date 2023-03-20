mod maths;
pub mod physics;

use std::{fs, time::Instant};

use crate::physics::{Enviroment, Physics};
use CFDplotlib::{Env, Plot};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let now = Instant::now();

    let env = Env::new();
    let mut plot = Plot::new(&env, 25);
    let title = "Cavity Flow Pressure+velocity";

    let enviroment_data =
        fs::read_to_string("Enviroments/cavity_flow.json").expect("Unable to read JSON file");

    let enviroment: Enviroment =
        serde_json::from_str(&enviroment_data).expect("JSON does not have correct format.");

    let mut physics = Physics::new(enviroment);

    let mut iter = 0;

    plot.pcolormesh(&physics.x, &physics.y, &physics.u, "viridis", title);
    //plot.quiver(&physics.x, &physics.y, &physics.u, &physics.v);
    plot.setup_animation();
    plot.update_frame(&physics.u, &physics.u, &physics.v, physics.nx, physics.ny);
    let mut var = physics.u.iter().fold(0.0, |acc,x|acc+x);

    while iter < 1000 &&!var.is_nan(){
        physics.iterate();
        var = physics.u.iter().fold(0.0, |acc,x|acc+x);
        plot.update_frame(&physics.p, &physics.u, &physics.v, physics.nx,physics.ny);

        iter += 1;
        dbg!(iter);
    }
    plot.update_frame(&physics.u, &physics.u, &physics.v,physics.nx, physics.ny);
    plot.finish_animation();

    plot.clf();
    plot.streamplot(&physics.x, &physics.y, &physics.u, &physics.v);
    plot.save("plot.png");
    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
    plot.clf();

    let plot_r = Plot::new(&env, 0);
    let x_axis: Vec<f32> = (0..iter).map(|i| i as f32).collect();
    plot_r.semilogy(&x_axis, &physics.residuals.u);
    plot_r.semilogy(&x_axis, &physics.residuals.v);
    plot_r.semilogy(&x_axis, &physics.residuals.pressure);
    let epsilon = vec![f32::EPSILON;iter];
    plot_r.semilogy(&x_axis, &epsilon);
    plot_r.save("residuals.png");

    Ok(())
}

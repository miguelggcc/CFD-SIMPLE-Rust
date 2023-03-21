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
    let mut var = physics.u.iter().fold(0.0, |acc, x| acc + x);

    while iter < 50 && !var.is_nan() {
        physics.iterate();
        var = physics.u.iter().fold(0.0, |acc, x| acc + x);
        //plot.update_frame(&physics.p, &physics.u, &physics.v, physics.nx,physics.ny);

        iter += 1;
        dbg!(iter);
    }
    plot.update_frame(&physics.u, &physics.u, &physics.v, physics.nx, physics.ny);
    plot.finish_animation();

    plot.clf();
    plot.streamplot(&physics.x, &physics.y, &physics.u, &physics.v);
    plot.save("plot.png");
    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
    plot.clf();

    plot.contourf(&physics.x, &physics.y, &physics.u, "jet", "u velocity");
    plot.save("u.png");
    plot.clf();

    plot.contourf(&physics.x, &physics.y, &physics.v, "jet", "v velocity");
    plot.save("v.png");
    plot.clf();

    plot.contourf(&physics.x, &physics.y, &physics.p, "jet", "p");
    plot.save("p.png");
    plot.clf();

    let mut div = vec![0.0; physics.nx * physics.ny];
    div.iter_mut().zip(physics.faces).for_each(|(div, f)| {
        *div = ((f.u_e - f.u_w) + (f.v_n - f.v_s)).abs();
    });
    plot.contourf(&physics.x, &physics.y, &div, "jet", "div");
    plot.save("div.png");
    plot.clf();

    let ghia_u = vec![
        0.0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, -0.20581,
        -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1.00000,
    ];
    let ghia_y = vec![
        0.0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5000, 0.6172, 0.7344,
        0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0,
    ];

    let ghia_v = vec![0.0,0.09233 ,0.10091 ,0.10890,0.12317 ,0.16077 ,0.17507  ,0.17527 ,0.05454 ,-0.24533,-0.22445 ,-0.16914 ,- 0.10313 ,-0.08864 ,-0.07391,-0.05906,0.0 ];
    let ghia_x = vec![0.0,0.0625,0.0703  ,0.0781  ,0.0938  ,0.1563  ,0.2266  ,0.2344  ,0.5000  ,0.8047 ,0.8594 ,0.9063 , 0.9453 ,0.9531 ,0.9609 ,0.9688 ,1.0];

    let mut ghia_u_simple = vec![];
    for j in 0..physics.ny {
        ghia_u_simple.push(physics.u[physics.nx / 2 + j * physics.nx]);
    }

    let mut ghia_v_simple = vec![];
    for i in 0..physics.nx {
        ghia_v_simple.push(physics.v[i + physics.ny/2 * physics.nx]);
    }

    plot.plot(&physics.y, &ghia_u_simple);
    plot.scatter(&ghia_y, &ghia_u);
    plot.plot(&physics.x, &ghia_v_simple);
    plot.scatter(&ghia_x, &ghia_v);

    plot.legend(&["u solution","u real","v solution","v real"]);

    plot.save("ghia.png");
    plot.clf();

    let plot_r = Plot::new(&env, 0);
    let x_axis: Vec<f64> = (0..iter).map(|i| i as f64).collect();
    plot_r.semilogy(&x_axis, &physics.residuals.u);
    plot_r.semilogy(&x_axis, &physics.residuals.v);
    plot_r.semilogy(&x_axis, &physics.residuals.pressure);
    plot_r.save("residuals.png");

    Ok(())
}

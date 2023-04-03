use std::path::Path;

use crate::{plotter::Plot, tools::vector_delta};

use super::PipeFlow;

impl PipeFlow {
    pub fn plot(&self, plot: &mut Plot, iter: u32) {
        let export_path = Path::new("./results_pipe_flow/");
        std::fs::create_dir_all(export_path).expect("Error creating path");

        let x_centers = vector_delta(self.dx * 0.5, self.dx, self.nx);
        let y_centers = vector_delta(self.dy * 0.5, self.dy, self.ny);
        let x = vector_delta(0.0, self.width / (self.nx - 1) as f64, self.nx);
        let y = vector_delta(0.0, self.height / (self.ny - 1) as f64, self.ny);

        plot.clf();
        let vel: Vec<f64> = self
            .u
            .iter()
            .zip(&self.v)
            .map(|(u, v)| (u * u + v * v).sqrt())
            .collect();
        plot.contourf(
            &x,
            &y,
            &vel,
            "plasma",
            &format!("Velocity magnitude (Re = {:.0})", self.re),
        );
        plot.streamplot(&x, &y, &self.u, &self.v, 1.0);
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("velocity_m.svg"));
        plot.clf();

        plot.pcolormesh(
            &x,
            &y,
            &self.u,
            "jet",
            &format!("u velocity (Re = {:.0})", self.re),
        );
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("u.svg"));
        plot.clf();

        plot.contourf(
            &x,
            &y,
            &self.v,
            "jet",
            &format!("v velocity (Re = {:.0})", self.re),
        );
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("v.png"));
        plot.clf();

        plot.contourf(
            &x,
            &y,
            &self.p,
            "jet",
            &format!("Pressure (Re = {:.0})", self.re),
        );
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("p.png"));
        plot.clf();

        let h = self.height*0.5;

        //Get u profile at nx-10 and compare it to analytical solution
        let mut u_profile = vec![];
        for j in 0..self.ny {
            u_profile.push(self.u[self.nx - 10 + j * self.nx]);
        }

        let u_profile_analytical: Vec<f64> = y
            .iter()
            .map(|y| 1.5 * self.u_in * (1.0 - (y - h).powi(2) / (h * h)))
            .collect();

        let x_profile = x_centers[self.nx - 10];

        plot.dpi(150.0);
        plot.plot(&u_profile, &y_centers);
        plot.plot_color(&u_profile_analytical, &y, "m--");
        plot.grid(true);
        plot.xlabel("u");
        plot.ylabel("y");
        plot.legend(&[
            &format!("u profile at x = {:.2}", x_profile),
            "u profile, analytical solution",
        ]);
        plot.save(export_path.join("u_profile.png"));
        plot.clf();

        //Get pressure drop at pipe center and compare it to  analytical solution
        let mut p_drop = vec![];
        for i in 0..self.nx {
            p_drop.push(self.p[i + self.ny / 2 * self.nx]);
        }

        let p_drop_analytical: Vec<f64> = x
            .iter()
            .map(|x| (-3.0 * self.nu * self.rho / (h * h) * self.u_in) * (x - self.width))
            .collect();

        plot.plot(&x_centers, &p_drop);
        plot.plot_color(&x, &p_drop_analytical, "m--");
        plot.grid(true);
        plot.xlabel("x");
        plot.ylabel("p");
        plot.legend(&["Pressure drop", "Pressure drop, analytical solution"]);
        plot.save(export_path.join("p_drop.png"));
        plot.clf();

        let iter_axis: Vec<f64> = (1..iter).map(|i| i as f64).collect();
        plot.semilogy(&iter_axis, &self.residuals.u[1..iter as usize]);
        plot.semilogy(&iter_axis, &self.residuals.v[1..iter as usize]);
        plot.semilogy(&iter_axis, &self.residuals.pressure[1..iter as usize]);
        plot.ylabel("Residual");
        plot.xlabel("Iterations");
        plot.grid(true);
        plot.legend(&[
            "u velocity residual",
            "v velocity residual",
            "Pressure correction residual",
        ]);

        plot.save(export_path.join("residuals.png"));
    }
}

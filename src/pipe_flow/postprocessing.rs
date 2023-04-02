use std::path::Path;

use crate::plotter::Plot;

use super::PipeFlow;

impl PipeFlow {
    pub fn plot(&self, plot: &mut Plot, iter: u32) {
        let export_path = Path::new("./results_pipe_flow/");
        std::fs::create_dir_all(export_path).expect("Error creating path");

        plot.clf();
        let vel: Vec<f64> = self
            .u
            .iter()
            .zip(&self.v)
            .map(|(u, v)| (u * u + v * v).sqrt())
            .collect();
        plot.contourf(&self.x, &self.y, &vel, "plasma", &format!("Velocity magnitude (Re = {:.0})", self.re));
        plot.streamplot(&self.x, &self.y, &self.u, &self.v, 1.0);
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("velocity_m.png"));
        plot.clf();

        plot.contourf(&self.x, &self.y, &self.u, "jet",&format!("u velocity (Re = {:.0})", self.re));
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("u.png"));
        plot.clf();

        plot.contourf(&self.x, &self.y, &self.v, "jet", &format!("v velocity (Re = {:.0})", self.re));
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("v.png"));
        plot.clf();

        plot.contourf(&self.x, &self.y, &self.p, "jet", &format!("Pressure (Re = {:.0})", self.re));
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("p.png"));
        plot.clf();

        let h = 0.5;
        let x_end = 5.0;

        //Get u profile at nx-10 and compare it to  analytical solution
        let mut u_profile = vec![];
        for j in 0..self.ny {
            u_profile.push(self.u[self.nx - 10 + j * self.nx]);
        }

        let u_profile_analytical: Vec<f64> = self
            .y
            .iter()
            .map(|y| 1.5 * self.u_in * (1.0 - (y - h).powi(2) / (h * h)))
            .collect();

        let x_profile = self.x[self.nx - 10];

        plot.plot(&u_profile, &self.y);
        plot.plot(&u_profile_analytical, &self.y);
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

        let p_drop_analytical: Vec<f64> = self
            .x
            .iter()
            .map(|x| (-3.0 * self.nu * self.rho / (h * h) * self.u_in) * (x - x_end))
            .collect();

        plot.plot(&self.x, &p_drop);
        plot.plot(&self.x, &p_drop_analytical);
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

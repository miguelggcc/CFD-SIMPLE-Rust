use std::path::Path;

use crate::plotter::Plot;

use super::BackwardFacingStep;

impl BackwardFacingStep {
    pub fn plot(&self, plot: &mut Plot, iter: u32) {
        let export_path = Path::new("./results_backward_facing_step/");
        std::fs::create_dir_all(export_path).expect("Error creating path");

        plot.clf();

        plot.contourf(&self.x, &self.y, &self.u, "jet", &format!("u velocity (Re = {:.0})", self.re));
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

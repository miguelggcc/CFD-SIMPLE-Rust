use std::path::Path;

use crate::plotter::Plot;

use super::BackwardFacingStep;

impl BackwardFacingStep {
    pub fn plot(&self, plot: &mut Plot, iter: usize) {
        let export_path = Path::new("./results_backward_facing_step/");
        plot.clf();
        plot.streamplot(&self.x, &self.y, &self.u, &self.v);
        plot.axis("scaled");
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("streamlines.png"));

        plot.clf();

        plot.pcolormesh(&self.x, &self.y, &self.u, "jet", "u velocity");
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("u.png"));
        plot.clf();

        plot.contourf(&self.x, &self.y, &self.v, "jet", "v velocity");
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("v.png"));
        plot.clf();

        plot.contourf(&self.x, &self.y, &self.p, "jet", "p");
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
        plot.contourf(&self.x, &self.y, &vel, "plasma", "velocity magnitude");
        plot.streamplot(&self.x, &self.y, &self.u, &self.v);
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("velocity.png"));
        plot.clf();

        let iter_axis: Vec<f64> = (0..iter).map(|i| i as f64).collect();
        plot.semilogy(&iter_axis, &self.residuals.u);
        plot.semilogy(&iter_axis, &self.residuals.v);
        plot.semilogy(&iter_axis, &self.residuals.pressure);
        plot.xlabel("Residual");
        plot.ylabel("iter");

        plot.save(export_path.join("residuals.png"));
    }
}

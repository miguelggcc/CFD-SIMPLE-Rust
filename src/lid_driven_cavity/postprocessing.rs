use std::path::Path;

use crate::plotter::Plot;

use super::LidDrivenCavity;

impl LidDrivenCavity {
    pub fn plot(&self, plot: &mut Plot, iter: usize) {
        let export_path = Path::new("./results_lid_driven_cavity/");
        plot.clf();
        plot.streamplot(&self.x, &self.y, &self.u, &self.v);
        plot.axis("scaled");
        plot.xlabel("x");
        plot.ylabel("y");
        plot.save(export_path.join("streamlines.png"));

        plot.clf();

        plot.contourf(&self.x, &self.y, &self.u, "jet", "u velocity");
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

        let ghia_u = vec![
            0.0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662, -0.21090, -0.20581,
            -0.13641, 0.00332, 0.23151, 0.68717, 0.73722, 0.78871, 0.84123, 1.00000,
        ];
        let ghia_y = vec![
            0.0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5000, 0.6172, 0.7344,
            0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0,
        ];

        let ghia_v = vec![
            0.0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507, 0.17527, 0.05454, -0.24533,
            -0.22445, -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0.0,
        ];
        let ghia_x = vec![
            0.0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5000, 0.8047, 0.8594,
            0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1.0,
        ];

        let mut ghia_u_solution = vec![];
        for j in 0..self.ny {
            ghia_u_solution.push(self.u[self.nx / 2 + j * self.nx]);
        }

        let mut ghia_v_solution = vec![];
        for i in 0..self.nx {
            ghia_v_solution.push(self.v[i + self.ny / 2 * self.nx]);
        }

        plot.plot(&self.y, &ghia_u_solution);
        plot.scatter(&ghia_y, &ghia_u);
        plot.plot(&self.x, &ghia_v_solution);
        plot.scatter(&ghia_x, &ghia_v);
        plot.grid(true);
        plot.legend(&[
            "u solution",
            "u from Ghia et al.",
            "v solution",
            "v from Ghia et al.",
        ]);
        plot.xlabel("x, y");
        plot.ylabel("u, v");
        plot.save(export_path.join("ghia.png"));
        plot.clf();

        let iter_axis: Vec<f32> = (0..iter).map(|i| i as f32).collect();
        plot.semilogy(&iter_axis, &self.residuals.u);
        plot.semilogy(&iter_axis, &self.residuals.v);
        plot.semilogy(&iter_axis, &self.residuals.pressure);
        plot.xlabel("Residual");
        plot.ylabel("iter");

        plot.save(export_path.join("residuals.png"));
    }
}

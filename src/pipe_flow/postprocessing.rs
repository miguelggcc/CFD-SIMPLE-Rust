use std::path::Path;

use crate::plotter::Plot;

use super::PipeFlow;

impl PipeFlow {
    pub fn plot(&self, plot: &mut Plot, iter: usize) {
        let export_path = Path::new("./results_pipe_flow/");

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

        let h = 0.5;

        //Get u profile at nx-10 and compare it to  analytical solution
        let mut u_profile = vec![];
        for j in 0..self.ny {
            u_profile.push(self.u[self.nx -10 + j * self.nx]);
        }

        let u_profile_analytical: Vec<f64> = self.y.iter().map(|y|1.5*self.u_in*(1.0-(y-h).powi(2)/(h*h))).collect();

        let x_profile = self.x[self.nx-10];

        plot.plot( &u_profile, &self.y);
        plot.plot( &u_profile_analytical, &self.y);
        plot.grid(true);
        plot.xlabel("y");
        plot.ylabel("u");
        plot.legend(&[&format!("u profile at x = {:.2}", x_profile), "u profile, analytical solution"]);
        plot.save(export_path.join("u_profile.png"));
        plot.clf();


        //Get pressure drop at pipe center and compare it to  analytical solution
        let mut p_drop = vec![];
        for i in 0..self.nx {
            p_drop.push(self.p[i + self.ny/2 * self.nx]);
        }

        let p_drop_analytical: Vec<f64> = self.x.iter().map(|x| (-3.0*self.nu*self.rho/(h*h)*self.u_in)*(x-5.0)).collect();


        plot.plot( &self.x, &p_drop);
        plot.plot( &self.x, &p_drop_analytical);
        plot.grid(true);
        plot.xlabel("x");
        plot.ylabel("p");
        plot.legend(&["pressure drop" , "pressure drop, analytical solution"]);
        plot.save(export_path.join("p_drop.png"));
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

mod correct_parameters;
mod face_velocity;
mod get_links_momentum;
mod get_links_pressure_correction;
mod postprocessing;
mod residuals;
mod solver;
mod solver_correction;

use itertools_num::linspace;

use crate::{tools::ix, Case};

use self::{face_velocity::Faces, residuals::Residuals};

pub struct BackwardFacingStep {
    pub nx: usize,
    pub ny: usize,
    pub nx_in: usize,
    pub ny_in: usize,
    pub u_in: f64,
    pub p_out: f64,
    pub re: f64, //Reynolds number
    pub dx: f64,
    pub dy: f64,
    pub nu: f64,
    pub rho: f64,
    pub relax_uv: f64,
    pub relax_p: f64,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub links: Vec<Links>,
    pub plinks: Vec<Links>,
    pub source_x: Vec<f64>,
    pub source_y: Vec<f64>,
    pub source_p: Vec<f64>,
    pub a_0: Vec<f64>,
    pub a_p0: Vec<f64>,
    pub faces: Vec<Faces>,
    pub u: Vec<f64>,
    pub v: Vec<f64>,
    pub p: Vec<f64>,
    pub pc: Vec<f64>,
    pub residuals: Residuals,
}

impl BackwardFacingStep {
    pub fn new(nx: usize, ny: usize, re: f64) -> Self {
        let dx = 0.1 / (nx as f64);
        let dy = 0.02 / (ny as f64);

        let u_in = 0.1;
        let p_out = 0.0;

        let nx_in = nx * 2 / 10 + 1;
        let ny_in = ny / 2;
        let x = linspace::<f64>(0.0, 0.1, nx).collect();
        let y = linspace::<f64>(0.0, 0.02, ny).collect();
        let mut u = vec![u_in; nx * ny];

        let v = vec![0.0; nx * ny];
        let p = vec![0.0; nx * ny];
        let pc = vec![0.0; nx * ny];
        let links = vec![Links::default(); nx * ny];
        let plinks = links.clone();
        let source_x = vec![0.0; nx * ny];
        let source_y = source_x.clone();
        let source_p = source_x.clone();
        let a_0 = vec![0.0; nx * ny];
        let a_p0 = a_0.clone();
        let mut faces = vec![Faces::new(u_in, u_in, 0.0, 0.0); ny * nx];

        let n = nx;
        for j in 0..ny_in {
            for i in 0..nx_in {
                u[ix(i, j, n)] = 0.0;
                faces[ix(i, j, n)].u_e = 0.0;
                faces[ix(i, j, n)].u_w = 0.0;
            }
        }
        let residuals = Residuals::default();

        let relax_uv = 0.8;
        let relax_p = 0.2;

        Self {
            nx,
            ny,
            nx_in,
            ny_in,
            u_in,
            p_out,
            re,
            nu: 2e-5,
            dx,
            dy,
            rho: 1.0,
            u,
            v,
            relax_uv,
            relax_p,
            x,
            y,
            links,
            plinks,
            source_x,
            source_y,
            source_p,
            a_0,
            a_p0,
            faces,
            p,
            pc,
            residuals,
        }
    }
}
impl Case for BackwardFacingStep {
    fn iterate(&mut self) -> bool {
        self.get_links_momentum();

        self.save_u_residual();
        let mut u = std::mem::take(&mut self.u);
        self.solver_correction(&mut u, &self.a_0, &self.links, &self.source_x, 2, 0.2);
        self.u = u;

        dbg!(self.u[ix(self.nx - 10, self.ny - 10, self.nx)]);

        self.save_v_residual();
        let mut v = std::mem::take(&mut self.v);
        self.solver_correction(&mut v, &self.a_0, &self.links, &self.source_y, 2, 0.2);
        self.v = v;
        dbg!(self.v[ix(self.nx_in, self.ny - 10, self.nx)]);
        self.get_face_velocities();
        self.get_links_pressure_correction();

        self.save_pressure_residual();
        self.pc = vec![0.0; self.nx * self.ny];
        let mut pc = std::mem::take(&mut self.pc);
        self.solver_correction(&mut pc, &self.a_p0, &self.plinks, &self.source_p, 20, 0.0);
        self.pc = pc;
        dbg!(self.pc.iter().fold(0.0, |acc, x| acc + x.abs()));

        self.correct_cell_velocities();
        self.correct_face_velocities();
        self.correct_pressure();
        !self.u.iter().fold(0.0, |acc, x| acc + x.abs()).is_nan()
    }

    fn postprocessing(&self, plot: &mut crate::plotter::Plot<'_>, iter: usize) {
        self.plot(plot, iter);
    }
}

#[derive(Clone, Default)]
pub struct Links {
    pub a_e: f64,
    pub a_w: f64,
    pub a_n: f64,
    pub a_s: f64,
}

impl Links {
    pub fn set_links(&mut self, a_e: f64, a_w: f64, a_n: f64, a_s: f64) {
        self.a_e = a_e;
        self.a_w = a_w;
        self.a_n = a_n;
        self.a_s = a_s;
    }
}

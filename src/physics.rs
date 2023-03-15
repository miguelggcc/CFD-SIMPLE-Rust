use itertools_num::linspace;
use serde::{Deserialize, Serialize};

pub mod face_velocity;
pub mod get_links_momentum;
pub mod get_links_pressure_correction;
mod solver;

pub use face_velocity::*;
pub use get_links_momentum::*;
pub use get_links_pressure_correction::*;

pub struct Physics {
    pub nx: usize,
    pub ny: usize,
    pub nt: usize,
    pub dx: f32,
    pub dy: f32,
    pub nu: f32,
    pub rho: f32,
    pub re: f32, //Reynolds number
    pub dt: f32,
    pub f: [f32; 2],
    relax_uv: f32,
    relax_p: f32,
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub links: Vec<Links>,
    pub plinks: Vec<Links>,
    pub source_x: Vec<f32>,
    pub source_y: Vec<f32>,
    pub source_p: Vec<f32>,
    pub a_0: Vec<f32>,
    pub a_p0: Vec<f32>,
    pub faces: Vec<Faces>,
    pub u: Vec<f32>,
    pub v: Vec<f32>,
    pub p: Vec<f32>,
    pub pc: Vec<f32>,
    pub bc_n: BC,
    pub bc_w: BC,
    pub bc_s: BC,
    pub bc_e: BC,
}

impl Physics {
    pub fn new(e: Enviroment) -> Self {
        let dx = 1.0 / (e.nx as f32);
        let dy = 1.0 / (e.ny as f32);
        let dt = e.dt;
        let re = e.rho * 1.0 / e.nu;
        let x = linspace::<f32>(0.0, 1.0, e.nx).collect();
        let y = linspace::<f32>(0.0, 1.0, e.ny).collect();
        let mut u = vec![0.0; e.nx * e.ny];
        for i in 0..e.nx{
            u[ix(i,e.ny-1,e.nx)]=1.0;
        }
        let v = vec![0.0; e.nx * e.ny];
        let p = vec![0.0; e.nx * e.ny];
        let pc = vec![0.0; e.nx * e.ny];
        let links = vec![Links::default(); e.nx * e.ny];
        let plinks = links.clone();
        let source_x = vec![0.0; e.nx * e.ny];
        let source_y = source_x.clone();
        let source_p = source_x.clone();
        let a_0 = vec![0.0; e.nx * e.ny];
        let a_p0 = a_0.clone();
        let mut faces = vec![Faces::default(); e.ny * e.nx];
       
        Self {
            nx: e.nx,
            ny: e.ny,
            nt: e.nt,
            nu: e.nu,
            dx,
            dy,
            rho: e.rho,
            u,
            v,
            re,
            dt,
            relax_uv: 0.8,
            relax_p: 0.2,
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
            f: e.f,
            bc_n: e.bc_n,
            bc_w: e.bc_w,
            bc_s: e.bc_s,
            bc_e: e.bc_e,
        }
    }

    /*pub fn solve_x_momentum(&mut self) {
        let n = self.nx;
        for _ in 0..3 {
            for j in 1..self.ny - 1 {
                for i in 1..self.nx - 1 {
                    let l = &self.links[ix(i, j, n)];
                    self.u[ix(i, j, n)] = -(l.a_e * self.u[ix(i - 1, j, n)]
                        + l.a_w * self.u[ix(i + 1, j, n)]
                        + l.a_n * self.u[ix(i, j + 1, n)]
                        + l.a_s * self.u[ix(i, j - 1, n)]
                        - l.s_x)
                        / self.a_0[ix(i, j, n)];
                }
            }
        }
    }*/

    /*pub fn solve_y_momentum(&mut self) {
        let n = self.nx;
        for _ in 0..3 {
            for j in 1..self.ny - 1 {
                for i in 1..self.nx - 1 {
                    let l = &self.links[ix(i, j, n)];
                    self.v[ix(i, j, n)] = -(l.a_e * self.v[ix(i + 1, j, n)]
                        + l.a_w * self.v[ix(i - 1, j, n)]
                        + l.a_n * self.v[ix(i, j + 1, n)]
                        + l.a_s * self.v[ix(i, j - 1, n)]
                        - l.s_y)
                        / self.a_0[ix(i, j, n)];
                }
            }
        }
    }*/

    /*pub fn solve_pressure(&mut self) {
        let n = self.nx;
        for _ in 0..40 {
            for j in 1..self.ny - 1 {
                for i in 1..self.nx - 1 {
                    let pl = &self.plinks[ix(i, j, n)];
                    self.pc[ix(i, j, n)] = -(pl.a_e * self.pc[ix(i + 1, j, n)]
                        + pl.a_w * self.pc[ix(i - 1, j, n)]
                        + pl.a_n * self.pc[ix(i, j + 1, n)]
                        + pl.a_s * self.pc[ix(i, j - 1, n)]
                        - pl.s_x)
                        / self.a_p0[ix(i, j, n)];
                    if pl.s_x != 0.0 {
                        dbg!(pl.s_x);
                    }
                }
            }
        }
    }*/

    pub fn correct_velocity(&mut self) {
        let n = self.nx;

        //Interior cells
        for j in 1..self.ny - 1 {
            for i in 1..self.nx - 1 {
                let a_0 = self.a_0[ix(i, j, n)];
                let a_0_e = self.a_0[ix(i + 1, j, n)];
                let a_0_w = self.a_0[ix(i - 1, j, n)];
                let a_0_n = self.a_0[ix(i, j + 1, n)];
                let a_0_s = self.a_0[ix(i, j - 1, n)];

                self.u[ix(i, j, n)] += self.relax_uv
                    * 0.5
                    * (self.pc[ix(i - 1, j, n)] + self.pc[ix(i + 1, j, n)])
                    * self.dy
                    / a_0;
                self.v[ix(i, j, n)] += self.relax_uv
                    * 0.5
                    * (self.pc[ix(i, j - 1, n)] + self.pc[ix(i, j + 1, n)])
                    * self.dx
                    / a_0;

                let dpcdx_eface = self.pc[ix(i + 1, j, n)] - self.pc[ix(i, j, n)];
                let dpcdx_wface = self.pc[ix(i, j, n)] - self.pc[ix(i - 1, j, n)];
                let dpcdx_nface = self.pc[ix(i, j + 1, n)] - self.pc[ix(i, j, n)];
                let dpcdx_sface = self.pc[ix(i, j, n)] - self.pc[ix(i, j - 1, n)];

                let f = self.faces.get_mut(ix(i, j, n)).unwrap();
                f.u_e += self.relax_uv * 0.5 * self.dy * (1.0 / a_0 + 1.0 / a_0_e) * (dpcdx_eface);
                f.u_w += self.relax_uv * 0.5 * self.dy * (1.0 / a_0 + 1.0 / a_0_w) * (dpcdx_wface);
                f.v_n += self.relax_uv * 0.5 * self.dx * (1.0 / a_0 + 1.0 / a_0_n) * (dpcdx_nface);
                f.v_s += self.relax_uv * 0.5 * self.dx * (1.0 / a_0 + 1.0 / a_0_s) * (dpcdx_sface);
            }
        }

        //Bottom wall
        let j = 0;
        for i in 1..self.nx - 1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_e = self.a_0[ix(i + 1, j, n)];
            let a_0_w = self.a_0[ix(i - 1, j, n)];
            let a_0_n = self.a_0[ix(i, j + 1, n)];

            self.u[ix(i, j, n)] += self.relax_uv
                * 0.5
                * (self.pc[ix(i - 1, j, n)] + self.pc[ix(i + 1, j, n)])
                * self.dy
                / a_0;
            self.v[ix(i, j, n)] +=
                self.relax_uv * 0.5 * (self.pc[ix(i, j, n)] + self.pc[ix(i, j + 1, n)]) * self.dx
                    / a_0;

            let dpcdx_eface = self.pc[ix(i + 1, j, n)] - self.pc[ix(i, j, n)];
            let dpcdx_wface = self.pc[ix(i, j, n)] - self.pc[ix(i - 1, j, n)];
            let dpcdx_nface = self.pc[ix(i, j + 1, n)] - self.pc[ix(i, j, n)];

            let f = self.faces.get_mut(ix(i, j, n)).unwrap();
            f.u_e += self.relax_uv * 0.5 * self.dy * (1.0 / a_0 + 1.0 / a_0_e) * (dpcdx_eface);
            f.u_w += self.relax_uv * 0.5 * self.dy * (1.0 / a_0 + 1.0 / a_0_w) * (dpcdx_wface);
            f.v_n += self.relax_uv * 0.5 * self.dx * (1.0 / a_0 + 1.0 / a_0_n) * (dpcdx_nface);
        }

        //Left  wall
        let i = 0;
        for j in 1..self.ny - 1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_w = self.a_0[ix(i - 1, j, n)];
            let a_0_n = self.a_0[ix(i, j + 1, n)];
            let a_0_s = self.a_0[ix(i, j - 1, n)];

            self.u[ix(i, j, n)] +=
                self.relax_uv * 0.5 * (self.pc[ix(i - 1, j, n)] + self.pc[ix(i, j, n)]) * self.dy
                    / a_0;
            self.v[ix(i, j, n)] += self.relax_uv
                * 0.5
                * (self.pc[ix(i, j - 1, n)] + self.pc[ix(i, j + 1, n)])
                * self.dx
                / a_0;

            let dpcdx_wface = self.pc[ix(i, j, n)] - self.pc[ix(i - 1, j, n)];
            let dpcdx_nface = self.pc[ix(i, j + 1, n)] - self.pc[ix(i, j, n)];
            let dpcdx_sface = self.pc[ix(i, j, n)] - self.pc[ix(i, j - 1, n)];

            let f = self.faces.get_mut(ix(i, j, n)).unwrap();
            f.u_w += self.relax_uv * 0.5 * self.dy * (1.0 / a_0 + 1.0 / a_0_w) * (dpcdx_wface);
            f.v_n += self.relax_uv * 0.5 * self.dx * (1.0 / a_0 + 1.0 / a_0_n) * (dpcdx_nface);
            f.v_s += self.relax_uv * 0.5 * self.dx * (1.0 / a_0 + 1.0 / a_0_s) * (dpcdx_sface);
        }

        //Right wall
        let i = self.nx - 1;
        for j in 1..self.ny - 1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_w = self.a_0[ix(i - 1, j, n)];
            let a_0_n = self.a_0[ix(i, j + 1, n)];
            let a_0_s = self.a_0[ix(i, j - 1, n)];

            self.u[ix(i, j, n)] +=
                self.relax_uv * 0.5 * (self.pc[ix(i - 1, j, n)] + self.pc[ix(i, j, n)]) * self.dy
                    / a_0;
            self.v[ix(i, j, n)] += self.relax_uv
                * 0.5
                * (self.pc[ix(i, j - 1, n)] + self.pc[ix(i, j + 1, n)])
                * self.dx
                / a_0;

            let dpcdx_wface = self.pc[ix(i, j, n)] - self.pc[ix(i - 1, j, n)];
            let dpcdx_nface = self.pc[ix(i, j + 1, n)] - self.pc[ix(i, j, n)];
            let dpcdx_sface = self.pc[ix(i, j, n)] - self.pc[ix(i, j - 1, n)];

            let f = self.faces.get_mut(ix(i, j, n)).unwrap();
            f.u_w += self.relax_uv * 0.5 * self.dy * (1.0 / a_0 + 1.0 / a_0_w) * (dpcdx_wface);
            f.v_n += self.relax_uv * 0.5 * self.dx * (1.0 / a_0 + 1.0 / a_0_n) * (dpcdx_nface);
            f.v_s += self.relax_uv * 0.5 * self.dx * (1.0 / a_0 + 1.0 / a_0_s) * (dpcdx_sface);
        }

        //Top wall
        let j = self.ny - 1;
        for i in 1..self.nx - 1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_e = self.a_0[ix(i + 1, j, n)];
            let a_0_w = self.a_0[ix(i - 1, j, n)];
            let a_0_s = self.a_0[ix(i, j - 1, n)];

            self.u[ix(i, j, n)] += self.relax_uv
                * 0.5
                * (self.pc[ix(i - 1, j, n)] + self.pc[ix(i + 1, j, n)])
                * self.dy
                / a_0;
            self.v[ix(i, j, n)] +=
                self.relax_uv * 0.5 * (self.pc[ix(i, j - 1, n)] + self.pc[ix(i, j, n)]) * self.dx
                    / a_0;

            let dpcdx_eface = self.pc[ix(i + 1, j, n)] - self.pc[ix(i, j, n)];
            let dpcdx_wface = self.pc[ix(i, j, n)] - self.pc[ix(i - 1, j, n)];
            let dpcdx_sface = self.pc[ix(i, j, n)] - self.pc[ix(i, j - 1, n)];

            let f = self.faces.get_mut(ix(i, j, n)).unwrap();
            f.u_e += self.relax_uv * 0.5 * self.dy * (1.0 / a_0 + 1.0 / a_0_e) * (dpcdx_eface);
            f.u_w += self.relax_uv * 0.5 * self.dy * (1.0 / a_0 + 1.0 / a_0_w) * (dpcdx_wface);
            f.v_s += self.relax_uv * 0.5 * self.dx * (1.0 / a_0 + 1.0 / a_0_s) * (dpcdx_sface);
        }
    }
    pub fn correct_pressure(&mut self) {
        let n = self.nx;

        self.p
            .iter_mut()
            .zip(&self.pc)
            .for_each(|(p, pc)| *p += self.relax_p * pc);

        //bottom wall
        for i in 1..self.nx - 1 {
            self.p[ix(i, 0, n)] = self.p[ix(i, 1, n)];
        }

        //left wall
        for j in 1..self.ny - 1 {
            self.p[ix(0, j, n)] = self.p[ix(1, j, n)];
        }

        //right wall
        for j in 1..self.ny - 1 {
            self.p[ix(self.nx - 1, j, n)] = self.p[ix(self.nx - 2, j, n)];
        }

        //top wall
        for i in 1..self.nx - 1 {
            self.p[ix(i, self.ny - 1, n)] = self.p[ix(i, self.ny - 2, n)];
        }

        //bottom left corner
        self.p[ix(0, 0, n)] =
            (self.p[ix(1, 1, n)] + self.p[ix(0, 1, n)] + self.p[ix(1, 0, n)]) / 3.0;

        //bottom right corner
        self.p[ix(self.nx - 1, 0, n)] = (self.p[ix(self.nx - 2, 1, n)]
            + self.p[ix(self.nx - 1, 1, n)]
            + self.p[ix(self.nx - 2, 0, n)])
            / 3.0;

        //top left corner
        self.p[ix(0, self.ny - 1, n)] =
            (self.p[ix(1, self.ny - 2, n)] + self.p[ix(0, 1, n)] + self.p[ix(1, self.ny - 1, n)])
                / 3.0;

        //top right corner
        self.p[ix(self.nx - 1, self.ny - 1, n)] = (self.p[ix(self.nx - 2, self.ny - 2, n)]
            + self.p[ix(self.nx - 1, self.ny - 2, n)]
            + self.p[ix(self.nx - 2, self.ny - 1, n)])
            / 3.0;
    }

    pub fn iterate(&mut self) {
        self.get_links_momentum();
        //dbg!(self.v[ix(self.nx - 20, self.ny - 2, self.nx)]);

        let mut u = std::mem::take(&mut self.u);
        self.solver(&mut u, &self.a_0, &self.links, &self.source_x, 2);
        self.u = u;
        let mut v = std::mem::take(&mut self.v);
        self.solver(&mut v, &self.a_0, &self.links, &self.source_y, 2);
        self.v = v;

        self.face_velocity();
        self.get_links_pressure_correction();
        //dbg!(self.p.iter().fold(0.0, |acc, x| acc + x.abs()));
        let mut pc = std::mem::take(&mut self.pc);
        self.solver(&mut pc, &self.a_p0, &self.plinks, &self.source_p, 40);
        self.pc = pc;
        dbg!(self.pc[ix(self.nx - 20, self.ny - 2, self.nx)]);

        self.correct_velocity();
        dbg!(self.u[ix(self.nx - 20, self.ny - 2, self.nx)]);
        self.correct_pressure();
        dbg!(self.p[ix(self.nx - 20, self.ny - 2, self.nx)]);
    }
}

#[inline(always)]
pub fn ix(i: usize, j: usize, width: usize) -> usize {
    i + j * width
}

#[derive(Serialize, Deserialize)]
pub struct BC {
    u: BCType,
    v: BCType,
    p: BCType,
}

#[derive(Serialize, Deserialize)]
pub enum BCType {
    Nothing,
    Dirichlet(f32),
    DirichletVec(Vec<f32>),
    Neumann(f32),
    NeumannVec(Vec<f32>),
    Periodic,
}

#[derive(Serialize, Deserialize)]
pub struct Enviroment {
    pub nx: usize,
    pub ny: usize,
    pub nt: usize,
    pub nu: f32,
    pub rho: f32,
    pub dt: f32,
    pub f: [f32; 2],
    pub bc_n: BC,
    pub bc_w: BC,
    pub bc_s: BC,
    pub bc_e: BC,
}

#[derive(Clone, Default)]
pub struct Links {
    pub a_e: f32,
    pub a_w: f32,
    pub a_n: f32,
    pub a_s: f32,
}

impl Links {
    pub fn new(a_e: f32, a_w: f32, a_n: f32, a_s: f32) -> Self {
        Self { a_e, a_w, a_n, a_s }
    }

    pub fn set_links(&mut self, a_e: f32, a_w: f32, a_n: f32, a_s: f32) {
        self.a_e = a_e;
        self.a_w = a_w;
        self.a_n = a_n;
        self.a_s = a_s;
    }
}

impl Faces {
    pub fn new(u_e: f32, u_w: f32, v_n: f32, v_s: f32) -> Self {
        Self { u_e, u_w, v_n, v_s }
    }
}

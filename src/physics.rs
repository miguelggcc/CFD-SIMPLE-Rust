use itertools_num::linspace;
use serde::{Deserialize, Serialize};

pub struct Physics {
    pub nx: usize,
    pub ny: usize,
    pub nt: usize,
    pub dx: f32,
    pub dy: f32,
    pub nu: f32,
    pub rho: f32,
    pub dt: f32,
    pub f: [f32; 2],
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub u: Vec<f32>,
    pub v: Vec<f32>,
    pub p: Vec<f32>,
    pub bc_n: BC,
    pub bc_w: BC,
    pub bc_s: BC,
    pub bc_e: BC,
}

impl Physics {
    pub fn new(e: Enviroment) -> Self {
        let dx = 2.0 / (e.nx as f32 - 1.0);
        let dy = 2.0 / (e.ny as f32 - 1.0);
        let dt = e.dt;
        let x = linspace::<f32>(0.0, 2.0, e.nx).collect();
        let y = linspace::<f32>(0.0, 2.0, e.ny).collect();
        let u = Vec::default();
        let v = Vec::default();
        let p = Vec::default();
        Self {
            nx: e.nx,
            ny: e.ny,
            nt: e.nt,
            nu: e.nu,
            dx,
            dy,
            rho: e.rho,
            dt,
            x,
            y,
            u,
            v,
            p,
            f: e.f,
            bc_n: e.bc_n,
            bc_w: e.bc_w,
            bc_s: e.bc_s,
            bc_e: e.bc_e,
        }
    }

    pub fn initialize(&mut self) {
        let u = vec![0.0; self.nx * self.ny];
        let v = vec![0.0; self.nx * self.ny];
        let p = vec![0.0; self.nx * self.ny];

        /*for i in 1..self.nx - 1 {
            u[ix(i, self.ny-1, self.nx)] = 1.0;
        }*/
        self.u = u;
        self.v = v;
        self.p = p;
    }

    pub fn calculate_pressure(&mut self) {
        let n = self.nx;

        for _ in 0..50 {
            let pn = self.p.clone();

            for j in 1..self.ny - 1 {
                for i in 1..self.nx - 1 {
                    self.p[ix(i, j, n)] =
                        (self.dy * self.dy * (pn[ix(i + 1, j, n)] + pn[ix(i - 1, j, n)])
                            + self.dx * self.dx * (pn[ix(i, j + 1, n)] + pn[ix(i, j - 1, n)]))
                            / (2.0 * (self.dx * self.dx + self.dy * self.dy))
                            - self.rho * self.dx * self.dx * self.dy * self.dy
                                / (2.0 * (self.dx * self.dx + self.dy * self.dy))
                                * self.get_b(i, j, &self.u, &self.v);
                }
            }

            for j in 1..self.ny - 1 {
                // x=nx
                match self.bc_e.p {
                    BCType::Dirichlet(p0) => {
                        self.p[ix(self.nx - 1, j, n)] = p0;
                    }
                    BCType::Neumann(p0) => {
                        self.p[ix(self.nx - 1, j, n)] =
                            p0 * self.p[ix(self.nx - 1, j, n)] + self.p[ix(self.nx - 2, j, n)];
                    }
                    BCType::Periodic => {
                        self.p[ix(self.nx - 1, j, n)] = (self.dy
                            * self.dy
                            * (pn[ix(0, j, n)] + pn[ix(self.nx - 2, j, n)])
                            + self.dx
                                * self.dx
                                * (pn[ix(self.nx - 1, j + 1, n)] + pn[ix(self.nx - 1, j - 1, n)]))
                            / (2.0 * (self.dx * self.dx + self.dy * self.dy))
                            - self.rho * self.dx * self.dx * self.dy * self.dy
                                / (2.0 * (self.dx * self.dx + self.dy * self.dy))
                                * self.get_b(self.nx - 1, j, &self.u, &self.v);
                    }
                    _ => {}
                }

                // x=0
                match self.bc_w.p {
                    BCType::Dirichlet(p0) => {
                        self.p[ix(0, j, n)] = p0;
                    }
                    BCType::Neumann(p0) => {
                        self.p[ix(0, j, n)] = p0 * self.p[ix(0, j, n)] + self.p[ix(1, j, n)];
                    }
                    BCType::Periodic => {
                        self.p[ix(0, j, n)] =
                            (self.dy * self.dy * (pn[ix(1, j, n)] + pn[ix(self.nx - 1, j, n)])
                                + self.dx * self.dx * (pn[ix(0, j + 1, n)] + pn[ix(0, j - 1, n)]))
                                / (2.0 * (self.dx * self.dx + self.dy * self.dy))
                                - self.rho * self.dx * self.dx * self.dy * self.dy
                                    / (2.0 * (self.dx * self.dx + self.dy * self.dy))
                                    * self.get_b(0, j, &self.u, &self.v);
                    }
                    _ => {}
                }
            }

            for i in 0..self.nx {
                // y= 0
                match self.bc_s.p {
                    BCType::Dirichlet(p0) => {
                        self.p[ix(i, 0, n)] = p0;
                    }
                    BCType::Neumann(p0) => {
                        self.p[ix(i, 0, n)] = p0 * self.p[ix(i, 0, n)] + self.p[ix(i, 1, n)];
                    }
                    _ => {}
                }
                // y= ny
                match self.bc_n.p {
                    BCType::Dirichlet(p0) => {
                        self.p[ix(i, self.ny - 1, n)] = p0;
                    }
                    BCType::Neumann(p0) => {
                        self.p[ix(i, self.ny - 1, n)] =
                            p0 * self.p[ix(i, self.ny - 1, n)] + self.p[ix(i, self.ny - 2, n)];
                    }
                    _ => {}
                }
            }
        }
    }
    pub fn iterate(&mut self) {
        let n = self.nx;
        let un = self.u.clone();
        let vn = self.v.clone();

        self.calculate_pressure();

        for j in 1..self.ny - 1 {
            for i in 1..self.nx - 1 {
                self.u[ix(i, j, n)] = un[ix(i, j, n)]
                    - (self.dt / self.dx
                        * un[ix(i, j, n)]
                        * (un[ix(i, j, n)] - un[ix(i - 1, j, n)]))
                    - (self.dt / self.dy
                        * vn[ix(i, j, n)]
                        * (un[ix(i, j, n)] - un[ix(i, j - 1, n)]))
                    - (self.dt / (self.rho * 2.0 * self.dx)
                        * (self.p[ix(i + 1, j, n)] - self.p[ix(i - 1, j, n)]))
                    + (self.nu * self.dt / (self.dx * self.dx)
                        * (un[ix(i + 1, j, n)] - 2.0 * un[ix(i, j, n)] + un[ix(i - 1, j, n)]))
                    + (self.nu * self.dt / (self.dy * self.dy)
                        * (un[ix(i, j + 1, n)] - 2.0 * un[ix(i, j, n)] + un[ix(i, j - 1, n)]))
                    + self.dt * self.f[0];

                self.v[ix(i, j, n)] = vn[ix(i, j, n)]
                    - (self.dt / self.dx
                        * un[ix(i, j, n)]
                        * (vn[ix(i, j, n)] - vn[ix(i - 1, j, n)]))
                    - (self.dt / self.dy
                        * vn[ix(i, j, n)]
                        * (vn[ix(i, j, n)] - vn[ix(i, j - 1, n)]))
                    - self.dt / (self.rho * 2.0 * self.dy)
                        * (self.p[ix(i, j + 1, n)] - self.p[ix(i, j - 1, n)])
                    + (self.nu * self.dt / (self.dx * self.dx)
                        * (vn[ix(i + 1, j, n)] - 2.0 * vn[ix(i, j, n)] + vn[ix(i - 1, j, n)]))
                    + (self.nu * self.dt / (self.dy * self.dy)
                        * (vn[ix(i, j + 1, n)] - 2.0 * vn[ix(i, j, n)] + vn[ix(i, j - 1, n)]))
                    + self.dt * self.f[1];
            }
        }

        for j in 1..self.ny - 1 {
            // u at x=xn
            match self.bc_e.u {
                BCType::Dirichlet(u0) => {
                    self.u[ix(self.nx - 1, j, n)] = u0;
                }
                BCType::Neumann(u0) => {
                    self.u[ix(self.nx - 1, j, n)] =
                        u0 * self.u[ix(self.nx - 1, j, n)] + self.u[ix(self.nx - 2, j, n)];
                }
                BCType::Periodic => {
                    self.u[ix(self.nx - 1, j, n)] = un[ix(self.nx - 1, j, n)]
                        - (self.dt / self.dx
                            * un[ix(self.nx - 1, j, n)]
                            * (un[ix(self.nx - 1, j, n)] - un[ix(self.nx - 2, j, n)]))
                        - (self.dt / self.dy
                            * vn[ix(self.nx - 1, j, n)]
                            * (un[ix(self.nx - 1, j, n)] - un[ix(self.nx - 1, j - 1, n)]))
                        - (self.dt / (self.rho * 2.0 * self.dx)
                            * (self.p[ix(0, j, n)] - self.p[ix(self.nx - 2, j - 1, n)]))
                        + (self.nu * self.dt / (self.dx * self.dx)
                            * (un[ix(0, j, n)] - 2.0 * un[ix(self.nx - 1, j, n)]
                                + un[ix(self.nx - 2, j, n)]))
                        + (self.nu * self.dt / (self.dy * self.dy)
                            * (un[ix(self.nx - 1, j + 1, n)] - 2.0 * un[ix(self.nx - 1, j, n)]
                                + un[ix(self.nx - 1, j - 1, n)]))
                        + self.dt * self.f[0];
                }
                _ => {}
            }

            // u at x = 0
            match self.bc_w.u {
                BCType::Dirichlet(u0) => {
                    self.u[ix(0, j, n)] = u0;
                }
                BCType::Neumann(u0) => {
                    self.u[ix(0, j, n)] = u0 * self.p[ix(0, j, n)] + self.u[ix(1, j, n)];
                }
                BCType::Periodic => {
                    self.u[ix(0, j, n)] = un[ix(0, j, n)]
                        - (self.dt / self.dx
                            * un[ix(0, j, n)]
                            * (un[ix(0, j, n)] - un[ix(self.nx - 1, j, n)]))
                        - (self.dt / self.dy
                            * vn[ix(0, j, n)]
                            * (un[ix(0, j, n)] - un[ix(0, j - 1, n)]))
                        - (self.dt / (self.rho * 2.0 * self.dx)
                            * (self.p[ix(0, j, n)] - self.p[ix(self.nx - 2, j, n)]))
                        + (self.nu * self.dt / (self.dx * self.dx)
                            * (un[ix(1, j, n)] - 2.0 * un[ix(0, j, n)]
                                + un[ix(self.nx - 1, j, n)]))
                        + (self.nu * self.dt / (self.dy * self.dy)
                            * (un[ix(0, j + 1, n)] - 2.0 * un[ix(0, j, n)] + un[ix(0, j - 1, n)]))
                        + self.dt * self.f[0];
                }
                _ => {}
            }

            // v at x = xn
            match self.bc_e.v {
                BCType::Dirichlet(v0) => {
                    self.v[ix(self.nx - 1, j, n)] = v0;
                }
                BCType::Neumann(v0) => {
                    self.v[ix(self.nx - 1, j, n)] =
                        v0 * self.v[ix(self.nx - 1, j, n)] + self.v[ix(self.nx - 2, j, n)];
                }
                BCType::Periodic => {
                    self.v[ix(self.nx - 1, j, n)] = vn[ix(self.nx - 1, j, n)]
                        - (self.dt / self.dx
                            * un[ix(self.nx - 1, j, n)]
                            * (vn[ix(self.nx - 1, j, n)] - vn[ix(self.nx - 2, j, n)]))
                        - (self.dt / self.dy
                            * vn[ix(self.nx - 1, j, n)]
                            * (vn[ix(self.nx - 1, j, n)] - vn[ix(self.nx - 1, j - 1, n)]))
                        - self.dt / (self.rho * 2.0 * self.dy)
                            * (self.p[ix(self.nx - 1, j + 1, n)]
                                - self.p[ix(self.nx - 1, j - 1, n)])
                        + (self.nu * self.dt / (self.dx * self.dx)
                            * (vn[ix(0, j, n)] - 2.0 * vn[ix(self.nx - 1, j, n)]
                                + vn[ix(self.nx - 2, j, n)]))
                        + (self.nu * self.dt / (self.dy * self.dy)
                            * (vn[ix(self.nx - 1, j + 1, n)] - 2.0 * vn[ix(self.nx - 1, j, n)]
                                + vn[ix(self.nx - 1, j - 1, n)]))
                        + self.dt * self.f[1];
                }
                _ => {}
            }

            //v at x = 0

            match self.bc_w.v {
                BCType::Dirichlet(v0) => {
                    self.v[ix(0, j, n)] = v0;
                }
                BCType::Neumann(v0) => {
                    self.v[ix(0, j, n)] = v0 * self.p[ix(0, j, n)] + self.v[ix(1, j, n)];
                }
                BCType::Periodic => {
                    self.v[ix(0, j, n)] = vn[ix(0, j, n)]
                        - (self.dt / self.dx
                            * un[ix(0, j, n)]
                            * (vn[ix(0, j, n)] - vn[ix(self.nx - 1, j, n)]))
                        - (self.dt / self.dy
                            * vn[ix(0, j, n)]
                            * (vn[ix(0, j, n)] - vn[ix(0, j - 1, n)]))
                        - self.dt / (self.rho * 2.0 * self.dy)
                            * (self.p[ix(0, j + 1, n)] - self.p[ix(0, j - 1, n)])
                        + (self.nu * self.dt / (self.dx * self.dx)
                            * (vn[ix(1, j, n)] - 2.0 * vn[ix(0, j, n)]
                                + vn[ix(self.nx - 1, j, n)]))
                        + (self.nu * self.dt / (self.dy * self.dy)
                            * (vn[ix(0, j + 1, n)] - 2.0 * vn[ix(0, j, n)] + vn[ix(0, j - 1, n)]))
                        + self.dt * self.f[1];
                }
                _ => {}
            }
        }

        for i in 1..self.nx - 1 {
            // u at y = 0
            match self.bc_s.u {
                BCType::Dirichlet(u0) => {
                    self.u[ix(i, 0, n)] = u0;
                }
                BCType::Neumann(u0) => {
                    self.u[ix(i, 0, n)] = u0 * self.u[ix(i, 0, n)] + self.u[ix(i, 1, n)];
                }
                _ => {}
            }
            // u at y = yn
            match self.bc_n.u {
                BCType::Dirichlet(u0) => {
                    self.u[ix(i, self.ny - 1, n)] = u0;
                }
                BCType::Neumann(u0) => {
                    self.u[ix(i, self.ny - 1, n)] =
                        u0 * self.u[ix(i, self.ny - 1, n)] + self.u[ix(i, self.ny - 2, n)];
                }
                _ => {}
            }

            //v at y = 0
            match self.bc_s.v {
                BCType::Dirichlet(v0) => {
                    self.v[ix(i, 0, n)] = v0;
                }
                BCType::Neumann(v0) => {
                    self.v[ix(i, 0, n)] = v0 * self.v[ix(i, 0, n)] + self.v[ix(i, 1, n)];
                }
                _ => {}
            }
            //v at y = yn
            match self.bc_n.v {
                BCType::Dirichlet(v0) => {
                    self.v[ix(i, self.ny - 1, n)] = v0;
                }
                BCType::Neumann(v0) => {
                    self.v[ix(i, self.ny - 1, n)] =
                        v0 * self.v[ix(i, self.ny - 1, n)] + self.v[ix(i, self.ny - 2, n)];
                }
                _ => {}
            }
        }
    }

    #[inline(always)]
    fn get_b(&self, i: usize, j: usize, un: &[f32], vn: &[f32]) -> f32 {
        let n = self.nx;
        let edge = self.nx - 1; 

        if i == edge {
            1.0 / self.dt
                * ((un[ix(edge, j, n)] - un[ix(edge - 1, j, n)]) / self.dx / 2.0
                    + (vn[ix(edge, j + 1, n)] - vn[ix(edge, j - 1, n)]) / self.dy / 2.0)
                - ((un[ix(0, j, n)] - un[ix(edge - 1, j, n)]) / self.dx / 2.0).powi(2)
                - ((vn[ix(edge, j + 1, n)] - vn[ix(edge, j - 1, n)]) / self.dy / 2.0).powi(2)
                - 2.0 * (un[ix(edge, j + 1, n)] - un[ix(edge, j - 1, n)]) / self.dy / 2.0
                    * (vn[ix(edge, j, n)] - vn[ix(edge, j - 1, n)])
                    / self.dx
                    / 2.0
        } else if i == 0 {
            1.0 / self.dt
                * ((un[ix(1, j, n)] - un[ix(edge, j, n)]) / self.dx / 2.0
                    + (vn[ix(0, j + 1, n)] - vn[ix(0, j - 1, n)]) / self.dy / 2.0)
                - ((un[ix(1, j, n)] - un[ix(edge, j, n)]) / self.dx / 2.0).powi(2)
                - ((vn[ix(0, j + 1, n)] - vn[ix(0, j - 1, n)]) / self.dy / 2.0).powi(2)
                - 2.0 * (un[ix(0, j + 1, n)] - un[ix(0, j - 1, n)]) / self.dy / 2.0
                    * (vn[ix(1, j, n)] - vn[ix(0, j - 1, n)])
                    / self.dx
                    / 2.0
        } else {
            1.0 / self.dt
                * ((un[ix(i + 1, j, n)] - un[ix(i - 1, j, n)]) / self.dx / 2.0
                    + (vn[ix(i, j + 1, n)] - vn[ix(i, j - 1, n)]) / self.dy / 2.0)
                - ((un[ix(i + 1, j, n)] - un[ix(i - 1, j, n)]) / self.dx / 2.0).powi(2)
                - ((vn[ix(i, j + 1, n)] - vn[ix(i, j - 1, n)]) / self.dy / 2.0).powi(2)
                - 2.0 * (un[ix(i, j + 1, n)] - un[ix(i, j - 1, n)]) / self.dy / 2.0
                    * (vn[ix(i + 1, j, n)] - vn[ix(i, j - 1, n)])
                    / self.dx
                    / 2.0
        }
    }

    pub fn get_u(&self, i: usize, j: usize) -> f32 {
        self.u[ix(i, j, self.nx)]
    }

    pub fn get_v(&self, i: usize, j: usize) -> f32 {
        self.v[ix(i, j, self.nx)]
    }

    pub fn get_p(&self, i: usize, j: usize) -> f32 {
        self.p[ix(i, j, self.nx)]
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

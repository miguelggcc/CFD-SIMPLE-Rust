use super::Links;
pub use super::Physics;
use crate::{maths::tridiagonal_solver, physics::ix};

impl Physics {
    pub fn solver(
        &self,
        x: &mut [f32],
        a_0: &[f32],
        links: &[Links],
        sources: &[f32],
        iter: usize,
    ) {
        let n = self.nx;
        let mut diagonal = vec![0.0; self.nx * self.ny];
        let mut ax = vec![0.0; self.nx * self.ny - 1];
        let mut cx = vec![0.0; self.nx * self.ny - 1];
        let mut rhs = vec![0.0; self.nx * self.ny];

        for _ in 0..iter {
            //Interior nodes

            for j in 1..self.ny - 1 {
                for i in 1..self.nx - 1 {
                    let l = &links[ix(i, j, n)];
                    diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
                    ax[ix(i - 1, j, n)] = l.a_w;
                    cx[ix(i, j, n)] = l.a_e;
                    rhs[ix(i, j, n)] = -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)]
                        + sources[ix(i, j, n)];
                    rhs[ix(i, j, n)] += -l.a_e * x[ix(i + 1, j, n)]
                        - l.a_w * x[ix(i - 1, j, n)]
                        - a_0[ix(i, j, n)] * x[ix(i, j, n)];
                }
            }
            //Bottom left corner
            let j = 0;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
            cx[ix(i, j, n)] = l.a_e;
            rhs[ix(i, j, n)] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
            rhs[ix(i, j, n)] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Bottom wall
            let j = 0;

            for i in 1..self.nx - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
                ax[ix(i - 1, j, n)] = l.a_w;
                cx[ix(i, j, n)] = l.a_e;
                rhs[ix(i, j, n)] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
                rhs[ix(i, j, n)] += -l.a_e * x[ix(i + 1, j, n)]
                    - l.a_w * x[ix(i - 1, j, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Bottom right corner
            let j = 0;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
            ax[ix(i - 1, j, n)] = l.a_w;
            cx[ix(i, j, n)] = 0.0;
            rhs[ix(i, j, n)] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
            rhs[ix(i, j, n)] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Left wall
            let i = 0;

            for j in 1..self.ny - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
                ax[ix(i - 1, j, n)] = 0.0;
                cx[ix(i, j, n)] = l.a_e;
                rhs[ix(i, j, n)] =
                    -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[ix(i, j, n)] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Right wall
            let i = self.nx - 1;

            for j in 1..self.ny - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
                ax[ix(i - 1, j, n)] = l.a_w;
                cx[ix(i, j, n)] = 0.0;
                rhs[ix(i, j, n)] =
                    -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[ix(i, j, n)] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Top left corner
            let j = self.ny - 1;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
            ax[ix(i - 1, j, n)] = 0.0;
            cx[ix(i, j, n)] = l.a_e;
            rhs[ix(i, j, n)] = -l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
            rhs[ix(i, j, n)] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Top wall
            let j = self.ny - 1;

            for i in 1..self.nx - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
                ax[ix(i - 1, j, n)] = l.a_w;
                cx[ix(i, j, n)] = l.a_e;
                rhs[ix(i, j, n)] = -l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[ix(i, j, n)] += -l.a_e * x[ix(i + 1, j, n)]
                    - l.a_w * x[ix(i - 1, j, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
                //dbg!(l.a_e,l.a_w,l.a_s, a_0[ix(i, j, n)], rhs[ix(i,j,n)]);
            }
            //dbg!(links[ix(self.nx-20,self.ny-1,self.nx)].a_e,self.links[ix(self.nx-20,self.ny-1,self.nx)].a_w);

            //Top right corner
            let j = self.ny - 1;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
            ax[ix(i - 1, j, n)] = l.a_w;
            rhs[ix(i, j, n)] = -l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
            rhs[ix(i, j, n)] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //--------------------------Solve--------------------------------------------------------------------------

            tridiagonal_solver(&ax, &diagonal, &cx, &mut rhs);

            x.iter_mut().zip(&rhs).for_each(|(x, rhs)| *x += rhs);
        }

        //------------------------------------------------------------------------------------------------------------------------------------------------------

        for _ in 0..iter {
            //Interior nodes

            for j in 1..self.ny - 1 {
                for i in 1..self.nx - 1 {
                    let l = &links[ix(i, j, n)];
                    diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
                    ax[ix(i - 1, j, n)] = l.a_s;
                    cx[ix(i, j, n)] = l.a_n;
                    rhs[ix(i, j, n)] = -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)]
                        + sources[ix(i, j, n)];
                    rhs[ix(i, j, n)] += -l.a_n * x[ix(i, j + 1, n)]
                        - l.a_w * x[ix(i, j - 1, n)]
                        - a_0[ix(i, j, n)] * x[ix(i, j, n)];
                }
            }
            //Bottom left corner
            let j = 0;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
            cx[ix(i, j, n)] = l.a_n;
            rhs[ix(i, j, n)] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
            rhs[ix(i, j, n)] += -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Bottom wall
            let j = 0;

            for i in 1..self.nx - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
                ax[ix(i - 1, j, n)] = 0.0;
                cx[ix(i, j, n)] = l.a_n;
                rhs[ix(i, j, n)] =
                    -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[ix(i, j, n)] += -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Bottom right corner
            let j = 0;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
            ax[ix(i - 1, j, n)] = 0.0;
            cx[ix(i, j, n)] = l.a_n;
            rhs[ix(i, j, n)] = -l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
            rhs[ix(i, j, n)] += -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Left wall
            let i = 0;

            for j in 1..self.ny - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
                ax[ix(i - 1, j, n)] = l.a_s;
                cx[ix(i, j, n)] = l.a_n;
                rhs[ix(i, j, n)] =
                    -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[ix(i, j, n)] += -l.a_n * x[ix(i, j + 1, n)]
                    - -l.a_s * x[ix(i, j - 1, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Right wall
            let i = self.nx - 1;

            for j in 1..self.ny - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
                ax[ix(i - 1, j, n)] = l.a_s;
                cx[ix(i, j, n)] = l.a_n;
                rhs[ix(i, j, n)] =
                    -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[ix(i, j, n)] += -l.a_n * x[ix(i, j + 1, n)]
                    - l.a_s * x[ix(i, j - 1, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Top left corner
            let j = self.ny - 1;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
            ax[ix(i - 1, j, n)] = l.a_s;
            cx[ix(i, j, n)] = 0.0;
            rhs[ix(i, j, n)] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
            rhs[ix(i, j, n)] += -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Top wall
            let j = self.ny - 1;

            for i in 1..self.nx - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
                ax[ix(i - 1, j, n)] = l.a_s;
                cx[ix(i, j, n)] = 0.0;
                rhs[ix(i, j, n)] =
                    -l.a_e * x[ix(i - 1, j, n)] - l.a_w * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
                rhs[ix(i, j, n)] += -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];
                //dbg!(l.a_e,l.a_w,l.a_s, a_0[ix(i, j, n)], rhs[ix(i,j,n)]);
            }
            //dbg!(links[ix(self.nx-20,self.ny-1,self.nx)].a_e,self.links[ix(self.nx-20,self.ny-1,self.nx)].a_w);

            //Top right corner
            let j = self.ny - 1;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[ix(i, j, n)] = a_0[ix(i, j, n)];
            ax[ix(i - 1, j, n)] = l.a_s;
            rhs[ix(i, j, n)] = -l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
            rhs[ix(i, j, n)] += -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //--------------------------Solve--------------------------------------------------------------------------

            tridiagonal_solver(&ax, &diagonal, &cx, &mut rhs);

            x.iter_mut().zip(&rhs).for_each(|(x, rhs)| *x += rhs);
        }
    }
}

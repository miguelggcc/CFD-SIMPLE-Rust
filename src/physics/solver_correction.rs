use super::Links;
pub use super::Physics;
use crate::{maths::tridiagonal_solver, physics::ix};

impl Physics {
    pub fn solver_correction(
        &self,
        x: &mut [f32],
        a_0: &[f32],
        links: &[Links],
        sources: &[f32],
        iter: usize,
        dumping: f32,
    ) {
        let n = self.nx;

        for _ in 0..iter {
            let mut diagonal = vec![0.0; self.nx];
            let mut ax = vec![0.0; self.nx];
            let mut cx = vec![0.0; self.nx - 1];
            let mut rhs = vec![0.0; self.nx];

            //Bottom left corner
            let j = 0;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
            cx[i] = l.a_e;
            rhs[i] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
            rhs[i] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Bottom wall
            let j = 0;

            for i in 1..self.nx - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[i] = l.a_w;
                cx[i] = l.a_e;
                rhs[i] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
                rhs[i] += -l.a_e * x[ix(i + 1, j, n)]
                    - l.a_w * x[ix(i - 1, j, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Bottom right corner
            let j = 0;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
            ax[i] = l.a_w;
            rhs[i] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
            rhs[i] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            tridiagonal_solver(&ax, &diagonal, &cx, &mut rhs);

            replace_row(0, x, &rhs, self.nx);

            let mut diagonal = vec![0.0; self.nx];
            let mut ax = vec![0.0; self.nx];
            let mut cx = vec![0.0; self.nx - 1];
            let mut rhs = vec![0.0; self.nx];

            for j in 1..self.ny - 1 {
                //Left wall
                let i = 0;

                let l = &links[ix(i, j, n)];
                diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
                cx[i] = l.a_e;
                rhs[i] =
                    -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[i] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                //Interior nodes

                for i in 1..self.nx - 1 {
                    let l = &links[ix(i, j, n)];
                    diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
                    ax[i] = l.a_w;
                    cx[i] = l.a_e;
                    rhs[i] = -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)]
                        + sources[ix(i, j, n)];
                    rhs[i] += -l.a_e * x[ix(i + 1, j, n)]
                        - l.a_w * x[ix(i - 1, j, n)]
                        - a_0[ix(i, j, n)] * x[ix(i, j, n)];
                }

                //Right wall
                let i = self.nx - 1;

                let l = &links[ix(i, j, n)];
                diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[i] = l.a_w;
                rhs[i] =
                    -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[i] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                tridiagonal_solver(&ax, &diagonal, &cx, &mut rhs);

                replace_row(j, x, &rhs, self.nx);
            }

            let mut diagonal = vec![0.0; self.nx];
            let mut ax = vec![0.0; self.nx];
            let mut cx = vec![0.0; self.nx - 1];
            let mut rhs = vec![0.0; self.nx];

            //Top left corner
            let j = self.ny - 1;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
            cx[i] = l.a_e;
            rhs[i] = -l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
            rhs[i] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Top wall
            let j = self.ny - 1;

            for i in 1..self.nx - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[i] = l.a_w;
                cx[i] = l.a_e;
                rhs[i] = -l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[i] += -l.a_e * x[ix(i + 1, j, n)]
                    - l.a_w * x[ix(i - 1, j, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Top right corner
            let j = self.ny - 1;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
            ax[i] = l.a_w;
            rhs[i] = -l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
            rhs[i] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            tridiagonal_solver(&ax, &diagonal, &cx, &mut rhs);

            replace_row(self.nx - 1, x, &rhs, self.nx);

            //------------------------------------------------------------------------------------------------------------------------------------------------------

            let mut diagonal = vec![0.0; self.ny];
            let mut ax = vec![0.0; self.ny];
            let mut cx = vec![0.0; self.ny - 1];
            let mut rhs = vec![0.0; self.ny];

            //Bottom left corner
            let j = 0;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[j] = (1.0 + dumping) * a_0[ix(i, j, n)];
            cx[j] = l.a_n;
            rhs[j] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
            rhs[j] += -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Left wall
            let i = 0;

            for j in 1..self.ny - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[j] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[j] = l.a_s;
                cx[j] = l.a_n;
                rhs[j] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
                rhs[j] += -l.a_n * x[ix(i, j + 1, n)]
                    - l.a_s * x[ix(i, j - 1, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Top left corner
            let j = self.ny - 1;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[j] = (1.0 + dumping) * a_0[ix(i, j, n)];
            ax[j] = l.a_s;
            rhs[j] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
            rhs[j] += -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            tridiagonal_solver(&ax, &diagonal, &cx, &mut rhs);

            replace_column(0, x, &rhs, self.nx, self.ny);

            let mut diagonal = vec![0.0; self.ny];
            let mut ax = vec![0.0; self.ny];
            let mut cx = vec![0.0; self.ny - 1];
            let mut rhs = vec![0.0; self.ny];

            for i in 1..self.nx - 1 {
                //Bottom wall

                let j = 0;

                let l = &links[ix(i, j, n)];
                diagonal[j] = (1.0 + dumping) * a_0[ix(i, j, n)];
                cx[j] = l.a_n;
                rhs[j] =
                    -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[j] += -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                //Interior nodes

                for j in 1..self.ny - 1 {
                    let l = &links[ix(i, j, n)];
                    diagonal[j] = (1.0 + dumping) * a_0[ix(i, j, n)];
                    ax[j] = l.a_s;
                    cx[j] = l.a_n;
                    rhs[j] = -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)]
                        + sources[ix(i, j, n)];
                    rhs[j] += -l.a_n * x[ix(i, j + 1, n)]
                        - l.a_s * x[ix(i, j - 1, n)]
                        - a_0[ix(i, j, n)] * x[ix(i, j, n)];
                }

                //Top wall
                let j = self.ny - 1;
                let l = &links[ix(i, j, n)];
                diagonal[j] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[j] = l.a_s;
                rhs[j] =
                    -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[j] += -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                tridiagonal_solver(&ax, &diagonal, &cx, &mut rhs);

                replace_column(i, x, &rhs, self.nx, self.ny);
            }

            let mut diagonal = vec![0.0; self.ny];
            let mut ax = vec![0.0; self.ny];
            let mut cx = vec![0.0; self.ny - 1];
            let mut rhs = vec![0.0; self.ny];

            //Bottom right corner
            let j = 0;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[j] = (1.0 + dumping) * a_0[ix(i, j, n)];
            cx[j] = l.a_n;
            rhs[j] = -l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
            rhs[j] += -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Right wall
            let i = self.nx - 1;

            for j in 1..self.ny - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[j] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[j] = l.a_s;
                cx[j] = l.a_n;
                rhs[j] = -l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[j] += -l.a_n * x[ix(i, j + 1, n)]
                    - l.a_s * x[ix(i, j - 1, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Top right corner
            let j = self.ny - 1;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[j] = (1.0 + dumping) * a_0[ix(i, j, n)];
            ax[j] = l.a_s;
            rhs[j] = -l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
            rhs[j] += -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            tridiagonal_solver(&ax, &diagonal, &cx, &mut rhs);

            replace_column(self.nx - 1, x, &rhs, self.nx, self.ny);
        }
    }
}

#[inline(always)]
fn replace_row(row: usize, x: &mut [f32], solx: &[f32], nx: usize) {
    for i in 0..nx {
        x[ix(i, row, nx)] += solx[i];
    }
}

#[inline(always)]
fn replace_column(column: usize, x: &mut [f32], soly: &[f32], nx: usize, ny: usize) {
    for j in 0..ny {
        x[ix(column, j, nx)] += soly[j];
    }
}

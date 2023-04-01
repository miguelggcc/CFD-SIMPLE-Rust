use crate::tools::{ix, tridiagonal_solver};

use super::{BackwardFacingStep, Links};

impl BackwardFacingStep {
    pub fn solver_correction(
        &self,
        x: &mut [f64],
        a_0: &[f64],
        links: &[Links],
        sources: &[f64],
        iter: usize,
        dumping: f64,
    ) {
        let n = self.nx;
        let mut diagonal = vec![0.0; self.nx];
        let mut ax = vec![0.0; self.nx];
        let mut cx = vec![0.0; self.nx - 1];
        let mut rhs = vec![0.0; self.nx];

        for _ in 0..iter {
            //Bottom left corner
            let j = 0;
            let i = self.nx_in;

            let l = &links[ix(i, j, n)];
            diagonal[i - self.nx_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
            cx[i - self.nx_in] = l.a_e;
            rhs[i - self.nx_in] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
            rhs[i - self.nx_in] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Bottom wall
            let j = 0;

            for i in self.nx_in + 1..self.nx - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[i - self.nx_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[i - self.nx_in] = l.a_w;
                cx[i - self.nx_in] = l.a_e;
                rhs[i - self.nx_in] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
                rhs[i - self.nx_in] += -l.a_e * x[ix(i + 1, j, n)]
                    - l.a_w * x[ix(i - 1, j, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Bottom right corner
            let j = 0;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[i - self.nx_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
            ax[i - self.nx_in] = l.a_w;
            rhs[i - self.nx_in] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
            rhs[i - self.nx_in] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.nx - self.nx_in);

            replace_row(0, x, &rhs, (self.nx_in, self.nx), n);

            for j in 1..self.ny_in {
                //Down left wall
                let i = self.nx_in;

                let l = &links[ix(i, j, n)];
                diagonal[i - self.nx_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
                cx[i - self.nx_in] = l.a_e;
                rhs[i - self.nx_in] =
                    -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[i - self.nx_in] +=
                    -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                //Down interior nodes

                for i in self.nx_in + 1..self.nx - 1 {
                    let l = &links[ix(i, j, n)];
                    diagonal[i - self.nx_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
                    ax[i - self.nx_in] = l.a_w;
                    cx[i - self.nx_in] = l.a_e;
                    rhs[i - self.nx_in] = -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)]
                        + sources[ix(i, j, n)];
                    rhs[i - self.nx_in] += -l.a_e * x[ix(i + 1, j, n)]
                        - l.a_w * x[ix(i - 1, j, n)]
                        - a_0[ix(i, j, n)] * x[ix(i, j, n)];
                }

                //Down right wall
                let i = self.nx - 1;

                let l = &links[ix(i, j, n)];
                diagonal[i - self.nx_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[i - self.nx_in] = l.a_w;
                rhs[i - self.nx_in] =
                    -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[i - self.nx_in] +=
                    -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.nx - self.nx_in);

                replace_row(j, x, &rhs, (self.nx_in, self.nx), n);
            }

            for j in self.ny_in..self.ny - 1 {
                //Inlet left wall
                let i = 0;

                let l = &links[ix(i, j, n)];
                diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
                cx[i] = l.a_e;
                rhs[i] =
                    -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[i] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                //Up interior nodes

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

                //Up right wall
                let i = self.nx - 1;

                let l = &links[ix(i, j, n)];
                diagonal[i] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[i] = l.a_w;
                rhs[i] =
                    -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[i] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.nx);

                replace_row(j, x, &rhs, (0, self.nx), n);
            }

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

            tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.nx);

            replace_row(self.ny - 1, x, &rhs, (0, self.nx), n);

            //------------------------------------------------------------------------------------------------------------------------------------------------------
            let mut diagonal = vec![0.0; self.nx];
            let mut ax = vec![0.0; self.nx];
            let mut cx = vec![0.0; self.nx - 1];
            let mut rhs = vec![0.0; self.nx];

            //Bottom left corner
            let j = self.ny_in;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[j - self.ny_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
            cx[j - self.ny_in] = l.a_n;
            rhs[j - self.ny_in] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
            rhs[j - self.ny_in] += -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Inlet left wall
            let i = 0;

            for j in self.ny_in + 1..self.ny - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[j - self.ny_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[j - self.ny_in] = l.a_s;
                cx[j - self.ny_in] = l.a_n;
                rhs[j - self.ny_in] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
                rhs[j - self.ny_in] += -l.a_n * x[ix(i, j + 1, n)]
                    - l.a_s * x[ix(i, j - 1, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Top left corner
            let j = self.ny - 1;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[j - self.ny_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
            ax[j - self.ny_in] = l.a_s;
            rhs[j - self.ny_in] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
            rhs[j - self.ny_in] += -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.ny - self.ny_in);

            replace_column(0, x, &rhs, (self.ny_in, self.ny), n);

            for i in 1..self.nx_in {
                //Inlet Bottom wall

                let j = self.ny_in;

                let l = &links[ix(i, j, n)];
                diagonal[j - self.ny_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
                cx[j - self.ny_in] = l.a_n;
                rhs[j - self.ny_in] =
                    -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[j - self.ny_in] +=
                    -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                //Interior nodes

                for j in self.ny_in + 1..self.ny - 1 {
                    let l = &links[ix(i, j, n)];
                    diagonal[j - self.ny_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
                    ax[j - self.ny_in] = l.a_s;
                    cx[j - self.ny_in] = l.a_n;
                    rhs[j - self.ny_in] = -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)]
                        + sources[ix(i, j, n)];
                    rhs[j - self.ny_in] += -l.a_n * x[ix(i, j + 1, n)]
                        - l.a_s * x[ix(i, j - 1, n)]
                        - a_0[ix(i, j, n)] * x[ix(i, j, n)];
                }

                //Top wall
                let j = self.ny - 1;
                let l = &links[ix(i, j, n)];
                diagonal[j - self.ny_in] = (1.0 + dumping) * a_0[ix(i, j, n)];
                ax[j - self.ny_in] = l.a_s;
                rhs[j - self.ny_in] =
                    -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[j - self.ny_in] +=
                    -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.ny - self.ny_in);

                replace_column(i, x, &rhs, (self.ny_in, self.ny), n);
            }

            for i in self.nx_in..self.nx - 1 {
                //Outlet Bottom wall

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

                tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.ny);

                replace_column(i, x, &rhs, (0, self.ny), n);
            }

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

            tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.ny);

            replace_column(self.nx - 1, x, &rhs, (0, self.ny), n);
        }
    }
}

#[inline(always)]
fn replace_row(row: usize, x: &mut [f64], solx: &[f64], domain_x: (usize, usize), n: usize) {
    for i in domain_x.0..domain_x.1 {
        x[ix(i, row, n)] += solx[i - domain_x.0];
    }
}

#[inline(always)]
fn replace_column(column: usize, x: &mut [f64], soly: &[f64], domain_y: (usize, usize), n: usize) {
    for j in domain_y.0..domain_y.1 {
        x[ix(column, j, n)] += soly[j - domain_y.0];
    }
}

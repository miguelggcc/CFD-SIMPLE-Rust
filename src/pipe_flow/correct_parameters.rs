use crate::tools::ix;

use super::PipeFlow;

impl PipeFlow {
    pub fn correct_cell_velocities(&mut self) {
        let n = self.nx;
        //------u velocity------------------
        //Interior cells
        for j in 0..self.ny {
            for i in 1..self.nx - 1 {
                let a_0 = self.a_0[ix(i, j, n)];
                self.u[ix(i, j, n)] += self.relax_uv
                    * 0.5
                    * (self.pc[ix(i - 1, j, n)] - self.pc[ix(i + 1, j, n)])
                    * self.dy
                    / a_0;
            }
        }
        //Left wall
        let i = 0;
        for j in 0..self.ny {
            let a_0 = self.a_0[ix(i, j, n)];
            self.u[ix(i, j, n)] +=
                self.relax_uv * 0.5 * (self.pc[ix(i, j, n)] - self.pc[ix(i + 1, j, n)]) * self.dy
                    / a_0; //p0 - (pe+p0)/2
        }
        //Right wall
        let i = self.nx - 1;
        for j in 0..self.ny {
            let a_0 = self.a_0[ix(i, j, n)];
            self.u[ix(i, j, n)] +=
                self.relax_uv * 0.5 * (self.pc[ix(i - 1, j, n)] + self.pc[ix(i, j, n)]) * self.dy
                    / a_0;
        }

        //------v velocity------------------
        //Interior cells
        for j in 1..self.ny - 1 {
            for i in 0..self.nx {
                let a_0 = self.a_0[ix(i, j, n)];
                self.v[ix(i, j, n)] += self.relax_uv
                    * 0.5
                    * (self.pc[ix(i, j - 1, n)] - self.pc[ix(i, j + 1, n)])
                    * self.dx
                    / a_0;
            }
        }
        //Bottom wall
        let j = 0;
        for i in 0..self.nx {
            let a_0 = self.a_0[ix(i, j, n)];
            self.v[ix(i, j, n)] +=
                self.relax_uv * 0.5 * (self.pc[ix(i, j, n)] - self.pc[ix(i, j + 1, n)]) * self.dx
                    / a_0; //p0 - (pn+p0)/2
        }
        //Top wall
        let j = self.ny - 1;
        for i in 0..self.nx {
            let a_0 = self.a_0[ix(i, j, n)];
            self.v[ix(i, j, n)] +=
                self.relax_uv * 0.5 * (self.pc[ix(i, j - 1, n)] - self.pc[ix(i, j, n)]) * self.dx
                    / a_0; //(ps+p0)/2 - p0
        }
    }

    pub fn correct_face_velocities(&mut self) {
        let n = self.nx;
        //------u velocity------------------
        //Interior cells
        for j in 0..self.ny {
            for i in 1..self.nx - 1 {
                let a_0 = self.a_0[ix(i, j, n)];
                let a_0e = self.a_0[ix(i + 1, j, n)];
                let a_0w = self.a_0[ix(i - 1, j, n)];

                let f = self.faces.get_mut(ix(i, j, n)).unwrap();

                f.u_e += self.relax_uv
                    * 0.5
                    * (self.pc[ix(i, j, n)] - self.pc[ix(i + 1, j, n)])
                    * self.dy
                    / (1.0 / a_0 + 1.0 / a_0e);
                f.u_w += self.relax_uv
                    * 0.5
                    * (self.pc[ix(i - 1, j, n)] - self.pc[ix(i, j, n)])
                    * self.dy
                    / (1.0 / a_0 + 1.0 / a_0w);
            }
        }
        //Left wall
        let i = 0;
        for j in 0..self.ny {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0e = self.a_0[ix(i + 1, j, n)];

            let f = self.faces.get_mut(ix(i, j, n)).unwrap();

            f.u_e +=
                self.relax_uv * 0.5 * (self.pc[ix(i, j, n)] - self.pc[ix(i + 1, j, n)]) * self.dy
                    / (1.0 / a_0 + 1.0 / a_0e);
        }
        //Right wall
        let i = self.nx - 1;
        for j in 0..self.ny {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0w = self.a_0[ix(i - 1, j, n)];

            let f = self.faces.get_mut(ix(i, j, n)).unwrap();

            f.u_e = self.u[ix(i, j, n)];
            f.u_w +=
                self.relax_uv * 0.5 * (self.pc[ix(i - 1, j, n)] - self.pc[ix(i, j, n)]) * self.dy
                    / (1.0 / a_0 + 1.0 / a_0w);
        }

        //------v velocity------------------
        //Interior cells
        for j in 1..self.ny - 1 {
            for i in 0..self.nx {
                let a_0 = self.a_0[ix(i, j, n)];
                let a_0n = self.a_0[ix(i, j + 1, n)];
                let a_0s = self.a_0[ix(i, j - 1, n)];

                let f = self.faces.get_mut(ix(i, j, n)).unwrap();

                f.v_n += self.relax_uv
                    * 0.5
                    * (self.pc[ix(i, j, n)] - self.pc[ix(i, j + 1, n)])
                    * self.dx
                    / (1.0 / a_0 + 1.0 / a_0n);
                f.v_s += self.relax_uv
                    * 0.5
                    * (self.pc[ix(i, j - 1, n)] - self.pc[ix(i, j, n)])
                    * self.dx
                    / (1.0 / a_0 + 1.0 / a_0s);
            }
        }
        //Bottom wall
        let j = 0;
        for i in 0..self.nx {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0n = self.a_0[ix(i, j + 1, n)];

            let f = self.faces.get_mut(ix(i, j, n)).unwrap();

            f.v_n +=
                self.relax_uv * 0.5 * (self.pc[ix(i, j, n)] - self.pc[ix(i, j + 1, n)]) * self.dx
                    / (1.0 / a_0 + 1.0 / a_0n);
        }
        //Top wall
        let j = self.ny - 1;
        for i in 0..self.nx {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0s = self.a_0[ix(i, j - 1, n)];

            let f = self.faces.get_mut(ix(i, j, n)).unwrap();

            f.v_s +=
                self.relax_uv * 0.5 * (self.pc[ix(i, j - 1, n)] - self.pc[ix(i, j, n)]) * self.dx
                    / (1.0 / a_0 + 1.0 / a_0s);
        }
    }

    pub fn correct_pressure(&mut self) {
        self.p
            .iter_mut()
            .zip(&self.pc)
            .for_each(|(p, pc)| *p += self.relax_p * pc);
    }
}

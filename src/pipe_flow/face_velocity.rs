use crate::tools::ix;

use super::PipeFlow;

impl PipeFlow {
    pub fn get_face_velocities(&mut self) {
        let n = self.nx;

        #[inline(always)]
        fn pwim(
            vel: &[f64],
            a_0: &[f64],
            dpdx_0: f64,
            dpdx_other: f64,
            dpdx_face: f64,
            index: usize,
            index_other: usize,
            d: f64,
        ) -> f64 {
            0.5 * (vel[index] + vel[index_other])
                + 0.5
                    * (dpdx_0 / a_0[index] + dpdx_other / a_0[index_other]
                        - (1.0 / a_0[index] + 1.0 / a_0[index_other]) * dpdx_face)
                    * d
        }

        //-----------------uface-----------------------------------------------------------------------------------------------------

        //-------------------------------inner cells---------------------------------------------------
        for j in 0..self.ny {
            for i in 2..self.nx - 2 {
                let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

                let dpdx_0 = 0.5 * (self.p[ix(i + 1, j, n)] - self.p[ix(i - 1, j, n)]);
                let dpdx_e = 0.5 * (self.p[ix(i + 2, j, n)] - self.p[ix(i, j, n)]);
                let dpdx_w = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i - 2, j, n)]);
                let dpdx_eface = self.p[ix(i + 1, j, n)] - self.p[ix(i, j, n)];
                let dpdx_wface = self.p[ix(i, j, n)] - self.p[ix(i - 1, j, n)];

                f.u_e = pwim(
                    &self.u,
                    &self.a_0,
                    dpdx_0,
                    dpdx_e,
                    dpdx_eface,
                    ix(i, j, n),
                    ix(i + 1, j, n),
                    self.dy,
                );

                f.u_w = pwim(
                    &self.u,
                    &self.a_0,
                    dpdx_0,
                    dpdx_w,
                    dpdx_wface,
                    ix(i, j, n),
                    ix(i - 1, j, n),
                    self.dy,
                );
            }
        }

        //-------------------------------Left wall---------------------------------------------------
        let i = 1;
        for j in 0..self.ny {
            let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

            let dpdx_0 = 0.5 * (self.p[ix(i + 1, j, n)] - self.p[ix(i - 1, j, n)]);
            let dpdx_e = 0.5 * (self.p[ix(i + 2, j, n)] - self.p[ix(i, j, n)]);
            let dpdx_w = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i - 1, j, n)]); //(p0+pw)/2-pw
            let dpdx_eface = self.p[ix(i + 1, j, n)] - self.p[ix(i, j, n)];
            let dpdx_wface = self.p[ix(i, j, n)] - self.p[ix(i - 1, j, n)];

            f.u_e = pwim(
                &self.u,
                &self.a_0,
                dpdx_0,
                dpdx_e,
                dpdx_eface,
                ix(i, j, n),
                ix(i + 1, j, n),
                self.dy,
            );

            f.u_w = pwim(
                &self.u,
                &self.a_0,
                dpdx_0,
                dpdx_w,
                dpdx_wface,
                ix(i, j, n),
                ix(i - 1, j, n),
                self.dy,
            );
        }

        let i = 0;
        for j in 0..self.ny {
            let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

            let dpdx_0 = 0.5 * (self.p[ix(i + 1, j, n)] - self.p[ix(i, j, n)]); //(pe+p0)/2-p0
            let dpdx_e = 0.5 * (self.p[ix(i + 2, j, n)] - self.p[ix(i, j, n)]);
            let dpdx_eface = self.p[ix(i + 1, j, n)] - self.p[ix(i, j, n)];

            f.u_e = pwim(
                &self.u,
                &self.a_0,
                dpdx_0,
                dpdx_e,
                dpdx_eface,
                ix(i, j, n),
                ix(i + 1, j, n),
                self.dy,
            );

            f.u_w = self.u_in;
        }

        //-------------------------------Right wall---------------------------------------------------
        let i = self.nx - 2;
        for j in 0..self.ny {
            let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

            let dpdx_0 = 0.5 * (self.p[ix(i + 1, j, n)] - self.p[ix(i - 1, j, n)]);
            let dpdx_e = self.p_out - 0.5 * (self.p[ix(i + 1, j, n)] + self.p[ix(i, j, n)]);
            let dpdx_w = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i - 2, j, n)]);
            let dpdx_eface = self.p[ix(i + 1, j, n)] - self.p[ix(i, j, n)];
            let dpdx_wface = self.p[ix(i, j, n)] - self.p[ix(i - 1, j, n)];

            f.u_e = pwim(
                &self.u,
                &self.a_0,
                dpdx_0,
                dpdx_e,
                dpdx_eface,
                ix(i, j, n),
                ix(i + 1, j, n),
                self.dy,
            );

            f.u_w = pwim(
                &self.u,
                &self.a_0,
                dpdx_0,
                dpdx_w,
                dpdx_wface,
                ix(i, j, n),
                ix(i - 1, j, n),
                self.dy,
            );
        }

        let i = self.nx - 1;
        for j in 0..self.ny {
            let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

            let dpdx_0 = self.p_out - 0.5 * (self.p[ix(i, j, n)] + self.p[ix(i - 1, j, n)]);
            let dpdx_w = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i - 2, j, n)]);
            let dpdx_wface = self.p[ix(i, j, n)] - self.p[ix(i - 1, j, n)];

            f.u_e = self.u[ix(i, j, n)];

            f.u_w = pwim(
                &self.u,
                &self.a_0,
                dpdx_0,
                dpdx_w,
                dpdx_wface,
                ix(i, j, n),
                ix(i - 1, j, n),
                self.dy,
            );
        }

        //-----------------v face-----------------------------------------------------------------------------------------------------

        //-------------------------------inner cells---------------------------------------------------
        for j in 2..self.ny - 2 {
            for i in 0..self.nx {
                let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

                let dpdy_0 = 0.5 * (self.p[ix(i, j + 1, n)] - self.p[ix(i, j - 1, n)]);
                let dpdy_n = 0.5 * (self.p[ix(i, j + 2, n)] - self.p[ix(i, j, n)]);
                let dpdy_s = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i, j - 2, n)]);
                let dpdy_nface = self.p[ix(i, j + 1, n)] - self.p[ix(i, j, n)];
                let dpdy_sface = self.p[ix(i, j, n)] - self.p[ix(i, j - 1, n)];

                f.v_n = pwim(
                    &self.v,
                    &self.a_0,
                    dpdy_0,
                    dpdy_n,
                    dpdy_nface,
                    ix(i, j, n),
                    ix(i, j + 1, n),
                    self.dx,
                );

                f.v_s = pwim(
                    &self.v,
                    &self.a_0,
                    dpdy_0,
                    dpdy_s,
                    dpdy_sface,
                    ix(i, j, n),
                    ix(i, j - 1, n),
                    self.dx,
                );
            }
        }

        //-------------------------------Top wall---------------------------------------------------
        let j = self.ny - 2;
        for i in 0..self.nx {
            let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

            let dpdy_0 = 0.5 * (self.p[ix(i, j + 1, n)] - self.p[ix(i, j - 1, n)]);
            let dpdy_n = 0.5 * (self.p[ix(i, j + 1, n)] - self.p[ix(i, j, n)]); //pn-(p0+pn)/2
            let dpdy_s = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i, j - 2, n)]);
            let dpdy_nface = self.p[ix(i, j + 1, n)] - self.p[ix(i, j, n)];
            let dpdy_sface = self.p[ix(i, j, n)] - self.p[ix(i, j - 1, n)];

            f.v_n = pwim(
                &self.v,
                &self.a_0,
                dpdy_0,
                dpdy_n,
                dpdy_nface,
                ix(i, j, n),
                ix(i, j + 1, n),
                self.dx,
            );

            f.v_s = pwim(
                &self.v,
                &self.a_0,
                dpdy_0,
                dpdy_s,
                dpdy_sface,
                ix(i, j, n),
                ix(i, j - 1, n),
                self.dx,
            );
        }

        let j = self.ny - 1;
        for i in 0..self.nx {
            let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

            let dpdy_0 = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i, j - 1, n)]); //(p0+pn)/2-pn
            let dpdy_s = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i, j - 2, n)]);
            let dpdy_sface = self.p[ix(i, j, n)] - self.p[ix(i, j - 1, n)];

            f.v_n = 0.0;

            f.v_s = pwim(
                &self.v,
                &self.a_0,
                dpdy_0,
                dpdy_s,
                dpdy_sface,
                ix(i, j, n),
                ix(i, j - 1, n),
                self.dx,
            );
        }

        //-------------------------------Bottom wall---------------------------------------------------
        let j = 1;
        for i in 0..self.nx {
            let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

            let dpdy_0 = 0.5 * (self.p[ix(i, j + 1, n)] - self.p[ix(i, j - 1, n)]);
            let dpdy_n = 0.5 * (self.p[ix(i, j + 2, n)] - self.p[ix(i, j, n)]);
            let dpdy_s = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i, j - 1, n)]); //(p0+ps)/2-ps
            let dpdy_nface = self.p[ix(i, j + 1, n)] - self.p[ix(i, j, n)];
            let dpdy_sface = self.p[ix(i, j, n)] - self.p[ix(i, j - 1, n)];

            f.v_n = pwim(
                &self.v,
                &self.a_0,
                dpdy_0,
                dpdy_n,
                dpdy_nface,
                ix(i, j, n),
                ix(i, j + 1, n),
                self.dx,
            );

            f.v_s = pwim(
                &self.v,
                &self.a_0,
                dpdy_0,
                dpdy_s,
                dpdy_sface,
                ix(i, j, n),
                ix(i, j - 1, n),
                self.dx,
            );
        }

        let j = 0;
        for i in 0..self.nx {
            let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

            let dpdy_0 = 0.5 * (self.p[ix(i, j + 1, n)] - self.p[ix(i, j, n)]); //(p0+ps)/2-ps
            let dpdy_n = 0.5 * (self.p[ix(i, j + 2, n)] - self.p[ix(i, j, n)]);
            let dpdy_nface = self.p[ix(i, j + 1, n)] - self.p[ix(i, j, n)];

            f.v_n = pwim(
                &self.v,
                &self.a_0,
                dpdy_0,
                dpdy_n,
                dpdy_nface,
                ix(i, j, n),
                ix(i, j + 1, n),
                self.dx,
            );

            f.v_s = 0.0;
        }
    }
}

#[derive(Clone, Default)]
pub struct Faces {
    pub u_e: f64,
    pub u_w: f64,
    pub v_n: f64,
    pub v_s: f64,
}

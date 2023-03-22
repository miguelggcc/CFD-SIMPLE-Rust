use crate::tools::ix;

use super::LidDrivenCavity;

impl LidDrivenCavity {
    pub fn get_face_velocities(&mut self) {
        let n = self.nx;

        #[inline(always)]
        fn pwim(
            vel: &[f32],
            a_0: &[f32],
            dpdx_0: f32,
            dpdx_other: f32,
            dpdx_face: f32,
            i: usize,
            j: usize,
            i_other: usize,
            j_other: usize,
            d: f32,
            n: usize,
        ) -> f32 {
            0.5 * (vel[ix(i, j, n)] + vel[ix(i_other, j_other, n)])
                + 0.5
                    * (dpdx_0 / a_0[ix(i, j, n)] + dpdx_other / a_0[ix(i_other, j_other, n)]
                        - (1.0 / a_0[ix(i, j, n)] + 1.0 / a_0[ix(i_other, j_other, n)]) * dpdx_face)
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
                    i,
                    j,
                    i + 1,
                    j,
                    self.dy,
                    n,
                );

                f.u_w = pwim(
                    &self.u,
                    &self.a_0,
                    dpdx_0,
                    dpdx_w,
                    dpdx_wface,
                    i,
                    j,
                    i - 1,
                    j,
                    self.dy,
                    n,
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
                i,
                j,
                i + 1,
                j,
                self.dy,
                n,
            );

            f.u_w = pwim(
                &self.u,
                &self.a_0,
                dpdx_0,
                dpdx_w,
                dpdx_wface,
                i,
                j,
                i - 1,
                j,
                self.dy,
                n,
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
                i,
                j,
                i + 1,
                j,
                self.dy,
                n,
            );

            f.u_w = 0.0;
        }

        //-------------------------------Right wall---------------------------------------------------
        let i = self.nx - 2;
        for j in 0..self.ny {
            let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

            let dpdx_0 = 0.5 * (self.p[ix(i + 1, j, n)] - self.p[ix(i - 1, j, n)]);
            let dpdx_e = 0.5 * (self.p[ix(i + 1, j, n)] - self.p[ix(i, j, n)]); //(pe+p0)/2-pe
            let dpdx_w = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i - 2, j, n)]);
            let dpdx_eface = self.p[ix(i + 1, j, n)] - self.p[ix(i, j, n)];
            let dpdx_wface = self.p[ix(i, j, n)] - self.p[ix(i - 1, j, n)];

            f.u_e = pwim(
                &self.u,
                &self.a_0,
                dpdx_0,
                dpdx_e,
                dpdx_eface,
                i,
                j,
                i + 1,
                j,
                self.dy,
                n,
            );

            f.u_w = pwim(
                &self.u,
                &self.a_0,
                dpdx_0,
                dpdx_w,
                dpdx_wface,
                i,
                j,
                i - 1,
                j,
                self.dy,
                n,
            );
        }

        let i = self.nx - 1;
        for j in 0..self.ny {
            let f: &mut Faces = self.faces.get_mut(ix(i, j, n)).unwrap();

            let dpdx_0 = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i - 1, j, n)]); //(pe+p0/2)-pe
            let dpdx_w = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i - 2, j, n)]);
            let dpdx_wface = self.p[ix(i, j, n)] - self.p[ix(i - 1, j, n)];

            f.u_e = 0.0;

            f.u_w = pwim(
                &self.u,
                &self.a_0,
                dpdx_0,
                dpdx_w,
                dpdx_wface,
                i,
                j,
                i - 1,
                j,
                self.dy,
                n,
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
                    i,
                    j,
                    i,
                    j + 1,
                    self.dx,
                    n,
                );

                f.v_s = pwim(
                    &self.v,
                    &self.a_0,
                    dpdy_0,
                    dpdy_s,
                    dpdy_sface,
                    i,
                    j,
                    i,
                    j - 1,
                    self.dx,
                    n,
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
                i,
                j,
                i,
                j + 1,
                self.dx,
                n,
            );

            f.v_s = pwim(
                &self.v,
                &self.a_0,
                dpdy_0,
                dpdy_s,
                dpdy_sface,
                i,
                j,
                i,
                j - 1,
                self.dx,
                n,
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
                i,
                j,
                i,
                j - 1,
                self.dx,
                n,
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
                i,
                j,
                i,
                j + 1,
                self.dx,
                n,
            );

            f.v_s = pwim(
                &self.v,
                &self.a_0,
                dpdy_0,
                dpdy_s,
                dpdy_sface,
                i,
                j,
                i,
                j - 1,
                self.dx,
                n,
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
                i,
                j,
                i,
                j + 1,
                self.dx,
                n,
            );

            f.v_s = 0.0;
        }
    }
}

#[derive(Clone, Default)]
pub struct Faces {
    pub u_e: f32,
    pub u_w: f32,
    pub v_n: f32,
    pub v_s: f32,
}

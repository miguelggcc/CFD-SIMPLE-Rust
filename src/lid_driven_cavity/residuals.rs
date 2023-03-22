//use std::io::Write;
//use std::{fs::File, io};

use crate::tools::ix;

use super::LidDrivenCavity;

impl LidDrivenCavity {
    pub fn save_pressure_residual(&mut self) {
        let n = self.nx;
        let mut res = 0.0;
        //let file = File::create("out.txt").unwrap();
        //let mut buffer = io::BufWriter::new(file);
        for j in 1..self.ny - 1 {
            //let mut line = String::new();
            for i in 1..self.nx - 1 {
                let f = &self.faces[ix(i, j, n)];
                let l = &self.plinks[ix(i, j, n)];

                let mdot = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;
                res += (l.a_e * self.pc[ix(i + 1, j, n)]
                    + l.a_w * self.pc[ix(i - 1, j, n)]
                    + l.a_n * self.pc[ix(i, j + 1, n)]
                    + l.a_s * self.pc[ix(i, j - 1, n)]
                    + self.a_p0[ix(i, j, n)] * self.pc[ix(i, j, n)]
                    - mdot)
                    .powi(2);

                //line.push_str(format!["{},{} ", f.u_w - f.u_e, f.v_n - f.v_s].as_str());
            }
            //writeln!(buffer, "{}", line).unwrap();
        }
        res = res.sqrt();
        self.residuals.pressure.push(res);
    }

    pub fn save_u_residual(&mut self) {
        let n = self.nx;
        let mut res = 0.0;
        //let file = File::create("out.txt").unwrap();
        //let mut buffer = io::BufWriter::new(file);
        for j in 1..self.ny - 1 {
            //let mut line = String::new();
            for i in 1..self.nx - 1 {
                let l = &self.links[ix(i, j, n)];

                let source_x = self.source_x[ix(i, j, n)];
                res += (l.a_e * self.u[ix(i + 1, j, n)]
                    + l.a_w * self.u[ix(i - 1, j, n)]
                    + l.a_n * self.u[ix(i, j + 1, n)]
                    + l.a_s * self.u[ix(i, j - 1, n)]
                    + self.a_0[ix(i, j, n)] * self.u[ix(i, j, n)]
                    - source_x)
                    .powi(2);

                /*line.push_str(
                    format![
                        "{} ",
                        (l.a_e * self.u[ix(i + 1, j, n)]
                + l.a_w * self.u[ix(i - 1, j, n)]
                + l.a_n * self.u[ix(i, j + 1, n)]
                + l.a_s * self.u[ix(i, j - 1, n)]
                + self.a_0[ix(i, j, n)] * self.u[ix(i, j, n)]
                - source_x)
                .powi(2)
                    ]
                    .as_str(),
                );*/
            }
            //writeln!(buffer, "{}", line).unwrap();
        }
        res = res.sqrt();
        self.residuals.u.push(res);
    }

    pub fn save_v_residual(&mut self) {
        let n = self.nx;
        let mut res = 0.0;
        for j in 1..self.ny - 1 {
            for i in 1..self.nx - 1 {
                let l = &self.links[ix(i, j, n)];
                let source_y = self.source_y[ix(i, j, n)];
                res += (l.a_e * self.v[ix(i + 1, j, n)]
                    + l.a_w * self.v[ix(i - 1, j, n)]
                    + l.a_n * self.v[ix(i, j + 1, n)]
                    + l.a_s * self.v[ix(i, j - 1, n)]
                    + self.a_0[ix(i, j, n)] * self.v[ix(i, j, n)]
                    - source_y)
                    .powi(2);
            }
        }
        res = res.sqrt();
        self.residuals.v.push(res);
    }
}

#[derive(Default)]
pub struct Residuals {
    pub pressure: Vec<f64>,
    pub u: Vec<f64>,
    pub v: Vec<f64>,
}

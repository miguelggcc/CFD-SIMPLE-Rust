pub fn tridiagonal_solver(a: &[f64], b: &[f64], c: &[f64], d: &mut [f64]) {
    let n = d.len();
    /*let mut w = vec![0.0; n];
        let mut g = b[0];
        w[0] = c[0] / g;
        d[0] = d[0] / g;
        for i in 1..n {
            g = 1.0 / (b[i] - a[i] * w[i - 1]);
            if i < n - 1 {
                w[i] = c[i] * g;
            }
            d[i] = (d[i] - a[i] * d[i - 1]) * g;
        }
        for i in (0..n - 1).rev() {
            d[i] -= w[i + 1] * d[i + 1];
        }
    }*/
    let mut c_prime = vec![0.0; n - 1];
    c_prime[0] = c[0] / b[0];
    d[0] /= b[0];

    for i in 1..n - 1 {
        let m = 1.0 / (b[i] - a[i] * c_prime[i - 1]);
        c_prime[i] = c[i] * m;
        d[i] = (d[i] - a[i] * d[i - 1]) * m;
    }
    d[n - 1] = (d[n - 1] - a[n - 1] * d[n - 2]) / (b[n - 1] - a[n - 1] * c_prime[n - 2]);

    for i in (0..n - 1).rev() {
        d[i] -= c_prime[i] * d[i + 1];
    }
}

pub fn tridiagonal_solver(a: &[f32], b: &[f32], c: &[f32], d: &mut [f32]) {
    let n = d.len();
    let mut w = vec![0.0; n];
    let mut g = b[0];
    d[0] = d[0] / g;
    for i in 1..n {
            w[i] = c[i - 1] / g;
        g = b[i] - a[i - 1] * w[i];
        d[i] = (d[i] - a[i - 1] * d[i - 1]) / g;
    }
    for i in (1..n - 1).rev() {
        d[i] -= w[i + 1] * d[i + 1];
    }
}

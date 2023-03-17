pub fn tridiagonal_solver(a: &[f32], b: &[f32], c: &[f32], d: &mut [f32]) {
    let n = d.len();
    let mut w = vec![0.0; n];
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
}

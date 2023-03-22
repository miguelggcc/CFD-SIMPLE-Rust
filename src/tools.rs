#[inline(always)]
pub fn ix(i: usize, j: usize, width: usize) -> usize {
    i + j * width
}

#[inline(always)]
pub fn tridiagonal_solver(a: &[f32], b: &[f32], c: &mut [f32], d: &mut [f32]) {

    //a: Vector with range 1..n (index 0 is never read)
    //b: Vector with range 0..n
    //c: Vector with range 0..n-1
    //d: Vector with range 0..n


    let n = d.len();

    c[0] /= b[0];
    d[0] /= b[0];

    assert!(c.len() + 1 == b.len());
    for i in 1..n - 1 {
        let m = 1.0 / (b[i] - a[i] * c[i - 1]);
        c[i] = c[i] * m;
        d[i] = (d[i] - a[i] * d[i - 1]) * m;
    }
    d[n - 1] = (d[n - 1] - a[n - 1] * d[n - 2]) / (b[n - 1] - a[n - 1] * c[n - 2]);

    for i in (0..n - 1).rev() {
        d[i] -= c[i] * d[i + 1];
    }
}

#[inline(always)]
pub fn ix(i: usize, j: usize, width: usize) -> usize {
    i + j * width
}

#[inline(always)]
pub fn tridiagonal_solver(a: &[f64], b: &[f64], c: &mut [f64], d: &mut [f64], n_eq: usize) {
    //a: Vector with range 1..n (index 0 is never read)
    //b: Vector with range 0..n
    //c: Vector with range 0..n-1
    //d: Vector with range 0..n

    c[0] /= b[0];
    d[0] /= b[0];

    assert!(c.len() + 1 == b.len());
    for i in 1..n_eq - 1 {
        let m = 1.0 / (b[i] - a[i] * c[i - 1]);
        c[i] *= m;
        d[i] = (d[i] - a[i] * d[i - 1]) * m;
    }
    d[n_eq - 1] =
        (d[n_eq - 1] - a[n_eq - 1] * d[n_eq - 2]) / (b[n_eq - 1] - a[n_eq - 1] * c[n_eq - 2]);

    assert!(c.len() + 1 == d.len());
    for i in (0..n_eq - 1).rev() {
        d[i] -= c[i] * d[i + 1];
    }
}

#[inline(always)]
pub fn scientific_notation_f64(num: &f64) -> String {
    let mut num = format!("{:.precision$e}", num, precision = 2);

    let e_loc = match num.find('e') {
        Some(num) => num,
        None => return String::from("NaN"),
    };

    let exp = num.split_off(e_loc);

    let (sign, exp) = if let Some(stripped) = exp.strip_prefix("e-") {
        ('-', stripped)
    } else {
        ('+', &exp[1..])
    };
    num.push_str(&format!("e{}{:0>pad$}", sign, exp, pad = 2));
    num
}

#[inline(always)]
pub fn vector_delta(start: f64, delta: f64, n: usize) -> Vec<f64> {
    let mut v = Vec::with_capacity(n);
    for i in 0..n {
        v.push(start + delta * i as f64);
    }
    v
}

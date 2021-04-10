//! Riemann and Hurwitz zeta function.
use num_traits::ToPrimitive;

use crate::complex::FromComplex;
use crate::consts::*;
use crate::utils::*;

/// Calculate the Riemann zeta function, which is defined as
/// 
/// $$\zeta(s)=\sum_{k=1}^{\infty}\frac{1}{k^{s}}$$
/// 
/// for $\mathfrak{R}[s]>1$. 
/// 
/// And it has a unique analytic continuation to entire complex plane, 
/// excluding the point $s=1$.
/// 
/// Then, a globally convergent series for the Riemann zeta function is given by
/// 
/// $$\zeta(s)=\frac{1}{1-2^{1-s}}\sum_{n=0}^{\infty}\frac{1}{2^{n+1}}\sum_{k=0}^{n}(-1)^{k}\binom{n}{k}(k+1)^{1-s}$$
pub fn zeta<T: FromComplex>(s: T) -> T {
    let one = T::one();

    if s == one {
        return T::from_f64(INFINITY).unwrap();
    }

    let two = one + one;

    let mut akn = vec![one];
    let mut two_pow = one / two;
    let head = one / (one - pow_os(two, s));
    let mut tail_prev = T::zero();
    let mut tail = two_pow * akn[0];

    while (tail - tail_prev).abs().to_f64().unwrap() >= EPSILON {
        update_akn(&mut akn, s);
        two_pow /= two;
        tail_prev = tail;
        tail += two_pow * akn.iter().map(|&e| e).sum();
    }

    head * tail
}

fn update_akn<T: FromComplex>(akn: &mut Vec<T>, s: T) {
    let n = akn.len() - 1;

    let n1 = T::from_usize(n + 1).unwrap();
    let one = T::one();

    let _ = akn.iter_mut()
                .enumerate()
                .map(|(k, a)| {
                    let num = n1;
                    let den = T::from_usize(n + 1 - k).unwrap();
                    *a *= num / den;
                }).collect::<()>();
    let p1 = -T::one() / n1;
    let p2 = if s.im().to_f64().unwrap() < EPSILON {
        (n1 / (n1 + one)).powf(s.re())
    } else {
        T::from_complex((n1 / (n1 + one)).powc(T::complex(s.re(), s.im())))
    };
    akn.push(p1 * p2 * akn[n]);
}

/// Calculate the Hurwitz zeta function, which is defined as
/// 
/// $$\zeta(s, q)=\sum_{k=0}^{\infty}\frac{1}{(q+k)^{s}}$$
/// 
/// for $\mathfrak{R}[s]>1$, where any term with $q+k=0$ is excluded. 
/// 
/// And it has a unique analytic continuation to entire complex plane, 
/// excluding the point $s=1$, where any term with $q+k=0$ is excluded. 
/// 
/// Then, a globally convergent series for the Hurwitz zeta function is given by
/// 
/// $$\zeta(s, q)=\frac{1}{s-1}\sum_{n=0}^{\infty}\frac{1}{n+1}\sum_{k=0}^{n}(-1)^{k}\binom{n}{k}(q+k)^{1-s}\tag{1}$$
/// 
/// This function calculate formula (1) directly.
pub fn zetah_raw<T: FromComplex>(s: T, q: T) -> T {
    let one = T::one();
    let real = s.re().to_f64().unwrap();
    let imag = s.im().to_f64().unwrap();
    
    if s == one {
        return T::from_f64(INFINITY).unwrap();
    }
    if imag < EPSILON && real < 0.0 && real.floor() == real {
        return T::from_f64(INFINITY).unwrap();
    }

    let mut qkos = vec![pow_os(q, s)];
    let mut fact: Vec<T> = vec![one];
    let mut bk_inv = one;
    let head = one / (s - one);
    let mut tail_prev = T::zero();
    let mut tail = qkos[0] * fact[0] / bk_inv;
    let mut dt_prev = T::from_f64(INFINITY).unwrap();
    
    while (tail - tail_prev).abs().to_f64().unwrap() >= EPSILON {
        let n1 = T::from_usize(qkos.len()).unwrap();
        let _ = fact.iter_mut()
                    .enumerate()
                    .skip(1)
                    .map(|(k, a)| {
                        let n1k = n1 - T::from_usize(k).unwrap();
                        *a = *a / n1k * n1;
                    })
                    .collect::<()>();
        fact.push(one);
        if qkos.len() % 2 == 0 {
            qkos.push(pow_os(q + n1, s));
        } else {
            qkos.push(-pow_os(q + n1, s));
        }

        bk_inv += one;
        
        let dt = fact.iter().zip(qkos.iter()).map(|(&f, &q)| f * q).sum::<T>() / bk_inv;
        if dt_prev.abs() <= dt.abs() {
            break;
        }

        tail_prev = tail;
        dt_prev = dt;
        tail += dt;
    }

    head * tail
}

/// The Gregory's coefficients for calculating the Hurwitz zeta function.
const G_N: [f64; 19] = [
    1.0 / 2.0,
    -1.0 / 12.0,
    1.0 / 24.0,
    -19.0 / 720.0,
    3.0 / 160.0,
    -863.0 / 60480.0,
    275.0 / 24192.0,
    -33953.0 / 3628800.0,
    8183.0 / 1036800.0,
    -3250433.0 / 479001600.0,
    4671.0 / 788480.0,
    -13695779093.0 / 2615348736000.0,
    2224234463.0 / 475517952000.0,
    -132282840127.0 / 3138418432000.0,
    2639651053.0 / 689762304000.0,
    -111956703448001.0 / 32011868528640000.0,
    50188465.0 / 15613165568.0,
    -2334028946344463.0 / 786014494949376000.0,
    301124035185049.0 / 109285437800448000.0,
];

fn update_gn<T: FromComplex>(gn: &mut Vec<T>) {
    let one = T::one();
    let n = gn.len() + 1;
    let nt = T::from_usize(n).unwrap();
    let head = if n % 2 == 1 {
        one / (nt + one)
    } else {
        -one / (nt + one)
    };
    let mut tail = T::zero();
    let mut on = if n % 2 == 1 {
        one
    } else {
        -one
    };
    let mut den = nt + one;
    for &gk in gn.iter() {
        den -= one;
        on = -on;
        tail += on * gk / den;
    }
    
    gn.push(head + tail)
}

/// Calculate the Hurwitz zeta function using the Gregory's coefficient.
/// 
/// The Gregory's coefficient is defined either via their generating function
/// 
/// $$\frac{z}{\ln{(1+z)}}=1+\sum_{n=1}^{\infty}G_{n}z^{n},\quad|z|<1,$$
/// 
/// or explicitly via the recurrence relation
/// 
/// $$G_{n}=\frac{(-1)^{n-1}}{n+1}+\sum_{k=1}^{n-1}\frac{(-1)^{n+1-k}G_{k}}{n+1-k},\quad G_{1}=\frac{1}{2},\quad n=2,3,4,\ldots$$
/// 
/// Using this, the Hurwitz zeta function can be represented as
/// 
/// $$\zeta(s, q)=\frac{q^{1-s}}{s-1}+\sum_{n=0}^{\infty}|G_{n+1}|\sum_{k=0}^{n}(-1)^{k}\binom{n}{k}(q+k)^{-s}$$
/// 
/// This globally series can be calculated faster than formula (1) and has higher accuracy.
/// 
/// (Please refer to Blagouchine 2018 "Three notes on Ser's and Hasse's representations for the zeta-functions".)
/// 
/// When $q=1$, the Hurwitz zeta function is same as the Riemann zeta function.
pub fn zetah<T: FromComplex>(s: T, q: T) -> T {
    let one = T::one();
    let real = s.re().to_f64().unwrap();
    let imag = s.im().to_f64().unwrap();
    
    if s == one {
        return T::from_f64(INFINITY).unwrap();
    }
    if q == one {
        return zeta(s);
    }
    if imag < EPSILON && real < 0.0 && real.floor() == real {
        return T::from_f64(INFINITY).unwrap();
    }

    let mut qkms = vec![pow_ms(q, s)];
    let mut fact: Vec<T> = vec![one];
    let mut gn = G_N.iter().map(|&e| T::from_f64(e).unwrap()).collect::<Vec<T>>();
    let n_calc = gn.len();
    let head = pow_os(q, s) / (s - one);
    let mut tail_prev = T::zero();
    let mut tail = qkms[0] * fact[0] * gn[0];
    let mut dt_prev = T::from_f64(INFINITY).unwrap();
    
    let mut n = 0;
    while (tail - tail_prev).abs().to_f64().unwrap() >= EPSILON {
        n += 1;
        let nt = T::from_usize(n).unwrap();
        let _ = fact.iter_mut()
                    .enumerate()
                    .skip(1)
                    .map(|(k, a)| {
                        let n1k = nt - T::from_usize(k).unwrap();
                        *a = *a / n1k * nt;
                    })
                    .collect::<()>();
        fact.push(one);
        if n >= n_calc {
            update_gn(&mut gn);
        }
        if qkms.len() % 2 == 0 {
            qkms.push(pow_ms(q + nt, s));
        } else {
            qkms.push(-pow_ms(q + nt, s));
        }
        
        let dt = fact.iter().zip(qkms.iter()).map(|(&f, &q)| f * q).sum::<T>() * T::from_real(gn[n].abs());
        if dt_prev.abs() <= dt.abs() {
            break;
        }

        tail_prev = tail;
        dt_prev = dt;
        tail += dt;
    }

    head + tail
}
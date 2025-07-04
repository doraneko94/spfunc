//! Riemann and Hurwitz zeta function.
use num_complex::Complex as NumComplex;
use num_traits::Float;
use num_traits::ToPrimitive;
use std::f64::consts::PI;

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
    let mut two_pow = one / two;

    if (s.re() - two_pow.re()).to_f64().unwrap() < EPSILON && s.im().to_f64().unwrap() > 0.0 {
        let c: NumComplex<T::Real> = zeta_riemann_siegel(s.im());
        return T::from_complex(T::complex(c.re, c.im));
    }

    let pow = pow_os(two, s);
    let denom = one - pow;

    let denom_abs = denom.abs().to_f64().unwrap_or(f64::INFINITY);
    if denom_abs < 1e-12 {
        return T::from_f64(f64::NAN).unwrap();
    }

    let mut akn = vec![one];
    let head = one / denom;

    let mut tail_prev = T::zero();
    let mut tail = two_pow * akn[0];

    let mut iter = 0;
    const MAX_ITER: usize = 500;

    while (tail - tail_prev).abs().to_f64().unwrap_or(f64::INFINITY) >= EPSILON {
        if iter >= MAX_ITER {
            break;
        }

        update_akn(&mut akn, s);
        two_pow = two_pow / two;
        tail_prev = tail;

        let term_sum: T = akn.iter().copied().sum();
        if term_sum.re().to_f64().unwrap().is_nan() || term_sum.re().to_f64().unwrap().is_infinite() {
            break;
        }

        tail += two_pow * term_sum;
        iter += 1;
    }

    let result = head * tail;
    if result.re().to_f64().unwrap_or(f64::INFINITY).is_nan() || result.re().to_f64().unwrap_or(0.0).is_infinite() {
        T::from_f64(INFINITY).unwrap()
    } else {
        result
    }
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
        let p = n1 / (n1 + one);
        if p.to_f64().unwrap() <= 0.0 {
            return;
        } else {
            p.powf(s.re())
        }
    } else {
        T::from_complex((n1 / (n1 + one)).powc(T::complex(s.re(), s.im())))
    };
    akn.push(p1 * p2 * akn[n]);
}

/// Riemann–Siegel theta function
fn theta<T: Float>(t: T) -> T {
    let one = T::one();
    let two = one + one;
    let eight = two + two + two + two;
    let pi = T::from(PI).unwrap();
    t / two * (t / (two * pi)).ln() - t / two - pi / eight
}

/// Riemann–Siegel Z(t) approximation
pub fn zeta_riemann_siegel<T: Float>(t: T) -> NumComplex<T> {
    if t <= T::zero() {
        panic!("Riemann–Siegel approximation is only valid for t > 0");
    }
    let one = T::one();
    let two = one + one;
    let pi = T::from(PI).unwrap();

    let n_max = (t / (two * pi)).sqrt().floor().to_usize().unwrap();
    let theta_t = theta(t);
    let mut z = T::zero();

    for n in 1..=n_max {
        let n_f = T::from(n).unwrap();
        z = z + (theta_t - t * n_f.ln()).cos() / n_f.sqrt();
    }

    // Z(t) is real, so return ζ(0.5 + it) = Z(t) * e^{-iθ(t)}
    let phase = NumComplex::from_polar(one, -theta_t);
    NumComplex::new(z, T::zero()) * phase
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

    let mut iter = 0;
    const MAX_ITER: usize = 500;
    
    while (tail - tail_prev).abs().to_f64().unwrap() >= EPSILON {
        if iter >= MAX_ITER {
            break;
        }

        let n1 = T::from_usize(qkos.len()).unwrap();

        for (k, a) in fact.iter_mut().enumerate().skip(1) {
            let n1k = n1 - T::from_usize(k).unwrap();
            *a = *a / n1k * n1;
        }
        fact.push(one);

        let term = if qkos.len() % 2 == 0 {
            pow_os(q + n1, s)
        } else {
            -pow_os(q + n1, s)
        };
        qkos.push(term);
        bk_inv += one;
        
        let dt = fact.iter().zip(qkos.iter()).map(|(&f, &q)| f * q).sum::<T>() / bk_inv;
        if dt.to_f64().unwrap().is_nan() || dt.to_f64().unwrap().is_infinite() {
            break;
        }
        if dt_prev.abs() <= dt.abs() {
            break;
        }

        tail_prev = tail;
        dt_prev = dt;
        tail += dt;
        iter += 1;
    }

    let result = head * tail;
    if result.re().to_f64().unwrap().is_nan() || result.re().to_f64().unwrap().is_infinite() {
        T::from_f64(INFINITY).unwrap()
    } else {
        result
    }
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
    let mut tail = qkms[0] * fact[0] * T::from_real(gn[0].abs());
    let mut dt_prev = T::from_f64(INFINITY).unwrap();
    
    let mut n = 0;
    const MAX_ITER: usize = 500;

    while (tail - tail_prev).abs().to_f64().unwrap() >= EPSILON {
        if n >= MAX_ITER {
            break;
        }

        n += 1;
        let nt = T::from_usize(n).unwrap();
        for (k, a) in fact.iter_mut().enumerate().skip(1) {
            let n1k = nt - T::from_usize(k).unwrap();
            *a = *a / n1k * nt;
        }
        fact.push(one);

        if n >= n_calc {
            update_gn(&mut gn);
        }

        let term = if qkms.len() % 2 == 0 {
            pow_os(q + nt, s)
        } else {
            -pow_os(q + nt, s)
        };
        qkms.push(term);
        
        let dt = fact.iter().zip(qkms.iter()).map(|(&f, &q)| f * q).sum::<T>() * T::from_real(gn[n].abs());
        if dt_prev.abs() <= dt.abs() {
            break;
        }

        tail_prev = tail;
        dt_prev = dt;
        tail += dt;
    }

    let result = head + tail;
    if result.re().to_f64().unwrap().is_nan() || result.re().to_f64().unwrap().is_infinite() {
        T::from_f64(INFINITY).unwrap()
    } else {
        result
    }
}
//! The family of the gamma function.
//! 
//! $$\Gamma(z)\equiv\int_{0}^{\infty}t^{z-1}e^{-t}dt\quad(\mathfrak{R}[z]>0)$$
use cauchy::Scalar;
use num_traits::ToPrimitive;

use crate::consts::SQRT_2_PI;
use crate::complex::FromComplex;
use crate::zeta::*;

/// Coefficients for calculating the gamma function.
const G_COF: [f64; 7] = [
    1.000000000190015,
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5,
];

/// Calculate $\ln{\Gamma(z)}$.
pub fn gamma_ln<T: Scalar>(z: T) -> T {
    if z.re().to_f64().unwrap() < 0.5 {
        let pi = T::from_f64(std::f64::consts::PI).unwrap();
        return pi.ln() - (pi * z).sin().ln() - gamma_ln(T::one() - z);
    }

    let g_cof = G_COF.iter().map(|&c| T::from_f64(c).unwrap()).collect::<Vec<T>>();
    let sqrt_2_pi = T::from_f64(SQRT_2_PI).unwrap();
    let tmp = z + T::from(5.5).unwrap();
    let mut ser = g_cof[0];
    for i in 1..7 {
        ser += g_cof[i] / (z + T::from(i).unwrap());
    }
    -tmp + (z + T::from(0.5).unwrap()) * tmp.ln() + (sqrt_2_pi * ser / z).ln()
}

/// Calculate $\Gamma(z)$.
/// 
/// The result is given as $\exp(\ln{\Gamma(z)})$.
pub fn gamma<T: Scalar>(z: T) -> T {
    gamma_ln(z).exp()
}

/// Coefficients for calculating the digamma function.
const D_COF: [f64; 9] = [
    12.0,
    -0.57721566490153286,
    1.6449340668482264365,
    1e-6,
    1.0 / 12.0,
    1.0 / 120.0,
    1.0 / 252.0,
    1.0 / 240.0,
    1.0 / 132.0,
];

/// Calculate the digamma function $\Psi(z)$, which is defined as the derivative of 
/// the natural logarithm of the gamma function.
/// 
/// $$\Psi(z)\equiv\frac{d}{dz}\ln{\Gamma(z)}$$
/// 
/// This implementation is based on statrs crate (statrs::function::gamma::digamma).
pub fn digamma<T: Scalar>(z: T) -> T {
    let one = T::one();
    let real = z.re().to_f64().unwrap();
    let imag = z.im().to_f64().unwrap();

    if imag < D_COF[3] {
        if real == std::f64::NEG_INFINITY || real.is_nan() {
            return T::from_f64(std::f64::NAN).unwrap()
        }
        if real <= 0.0 && real.floor() == real {
            return T::from_f64(std::f64::NEG_INFINITY).unwrap()
        }
        if real < 0.0 {
            let pi = T::from_f64(std::f64::consts::PI).unwrap();
            return digamma(one - z) + pi / (-pi * z).tan()
        }
        if real <= D_COF[3] {
            let d1 = T::from_f64(D_COF[1]).unwrap();
            let d2 = T::from_f64(D_COF[2]).unwrap();
            return d1 - one / z + d2 * z;
        }
    }

    let mut result = T::zero();
    let mut zr = z;
    while zr.re().to_f64().unwrap() < D_COF[0] {
        result -= one / zr;
        zr += one;
    }

    let zpf = T::from_f64(0.5).unwrap();
    let d_cof = D_COF.iter().map(|&e| T::from_f64(e).unwrap()).collect::<Vec<T>>();
    if zr.re().to_f64().unwrap() >= D_COF[0] {
        let mut r = one / zr;
        result += zr.ln() - zpf * r;
        r *= r;

        result -= r * (d_cof[4] - (r * (d_cof[5] - (r * (d_cof[6] - (r * (d_cof[7] - (r * d_cof[8]))))))));
    }

    result
}

/// Calculate the polygamma function $\Psi^{(n)}(z)$, which is defined as the n-th derivative of 
/// the digamma function.
/// 
/// $$\Psi^{(n)}(z)\equiv\frac{d^{n+1}}{dz^{n+1}}\ln{\Gamma(z)}=\frac{d^{n}}{dz^{n}}\Psi(z)$$
/// 
/// When n = 0, it is same as the digamma function $\Psi(z)$.
pub fn polygamma<T: FromComplex>(z: T, n: usize) -> T {
    if n == 0 {
        digamma(z)
    } else {
        let n1 = T::from(n + 1).unwrap();
        if (n + 1) % 2 == 0 {
            gamma(n1) * zetah(n1, z)
        } else {
            -gamma(n1) * zetah(n1, z)
        }
    }
}
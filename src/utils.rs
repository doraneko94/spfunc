//! Utility functions.
use num_traits::ToPrimitive;

use crate::complex::FromComplex;
use crate::consts::*;

/// Calculate $z^{-s}$.
pub fn pow_ms<T: FromComplex>(z: T, s: T) -> T {
    if z == T::zero() {
        return T::zero()
    }
    let imag = s.im().to_f64().unwrap();
    if imag < EPSILON {
        z.powf(-s.re())
    } else {
        T::from_complex(z.powc(-T::complex(s.re(), s.im())))
    }
}

/// Calculate $z^{1-s}$.
pub fn pow_os<T: FromComplex>(z: T, s: T) -> T {
    if z == T::zero() {
        return T::zero()
    }
    let imag = s.im().to_f64().unwrap();
    if imag < EPSILON {
        z.powf(T::real(1.0) - s.re())
    } else {
        T::from_complex(z.powc(T::complex(1.0, 0.0) - T::complex(s.re(), s.im())))
    }
}
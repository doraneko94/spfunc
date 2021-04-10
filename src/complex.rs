//! How to deal with complex numbers.
use cauchy::{Scalar, c32, c64};

/// Convert Scalar::Complex to Scalar.
pub trait FromComplex: Scalar {
    fn from_complex(c: Self::Complex) -> Self;
}

impl FromComplex for f32 {
    fn from_complex(c: Self::Complex) -> Self {
        c.re()
    }
}

impl FromComplex for f64 {
    fn from_complex(c: Self::Complex) -> Self {
        c.re()
    }
}

impl FromComplex for c32 {
    fn from_complex(c: Self::Complex) -> Self {
        Self::new(c.re(), c.im())
    }
}

impl FromComplex for c64 {
    fn from_complex(c: Self::Complex) -> Self {
        Self::new(c.re(), c.im())
    }
}
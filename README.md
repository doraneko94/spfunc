# spfunc
Rust crate for numerical calculation of special functions.

[![Crate](http://meritbadge.herokuapp.com/spfunc)](https://crates.io/crates/spfunc)
[![docs.rs](https://docs.rs/spfunc/badge.svg)](https://docs.rs/spfunc)

This crate can calculate each special function for f32, f64, Complex32, Complex64 (from num_complex crate).

# Note
This crate is still in the development stage and the numerical calculations are not so accurate (especially Hurwitz zeta function).

# Functions
## The Gamma Function
- The gamma function
- The digamma function
- The polygamma function

## The Zeta Function
- The Riemann zeta function
- The Hurwitz zeta function

# How to use
```rust
use spfunc::gamma::*;
use cauchy::{c32, c64};

fn main() {
    println!("Gamma(1.0) = {}", gamma(1.0));
    println!("ln(Gamma(1.0)) = {}", gamma_ln(1.0));

    println!("Gamma(1.2+3.4i) = {}", gamma(c32::new(1.2, 3.4)));
    println!("ln(Gamma(1.2+3.4i)) = {}", gamma_ln(c32::new(1.2, 3.4)));
    
    println!("DiGamma(1.2+3.4i) = {}", digamma(c64::new(1.2, 3.4)));
    println!("TriGamma(1.2+3.4i) = {}", polygamma(c64::new(1.2, 3.4), 1))
}
```
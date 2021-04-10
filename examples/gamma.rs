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
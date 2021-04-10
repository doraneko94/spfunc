use spfunc::zeta::*;
use cauchy::c64;

fn main() {
    println!("zeta(1.0) = {}", zeta(1.0));
    println!("zeta(2.0) = {}", zeta(2.0));
    println!("zeta(2.0, 1.0) = {}", zetah(2.0, 1.0));
    println!("zeta(1.2+3.4i) = {}", zeta(c64::new(1.2, 3.4)));
    println!("zeta(1.2+3.4i, 1.0+0.0i) = {}", zetah(c64::new(1.2, 3.4), c64::new(1.0, 0.0)));
    println!("zeta(1.2+3.4i, 5.6+7.8i) = {}", zetah(c64::new(1.2, 3.4), c64::new(5.6, 7.8)));
}
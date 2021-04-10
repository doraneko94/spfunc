#![feature(test)]

use spfunc::zeta::*;
extern crate test;

#[bench]
fn zetah_raw20(b: &mut test::Bencher) {
    b.iter(|| {
        for i in 2..20 {
            let a = i as f64;
            zetah_raw(a, a);
        }
    });
}

#[bench]
fn zetah_20(b: &mut test::Bencher) {
    b.iter(|| {
        for i in 2..20 {
            let a = i as f64;
            zetah(a, a);
        }
    });
}
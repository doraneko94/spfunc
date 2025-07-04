use criterion::{criterion_group, criterion_main, Criterion};
use spfunc::zeta::*;

fn benchmark_zetah_raw(c: &mut Criterion) {
    c.bench_function("zetah_raw (2..20)", |b| {
        b.iter(|| {
            for i in 2..20 {
                let a = i as f64;
                zetah_raw(a, a);
            }
        })
    });
}

fn benchmark_zetah(c: &mut Criterion) {
    c.bench_function("zetah (2..20)", |b| {
        b.iter(|| {
            for i in 2..20 {
                let a = i as f64;
                zetah(a, a);
            }
        })
    });
}

criterion_group!(benches, benchmark_zetah_raw, benchmark_zetah);
criterion_main!(benches);
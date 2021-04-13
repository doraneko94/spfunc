//! Example of obtaining the parameters of the gamma distribution 
//! by Maximum Likelihood Estimation (MLE).
//! 
//! For a description of the algorithm used in this code, see
//! 
//! https://tminka.github.io/papers/minka-gamma.pdf

use spfunc::gamma::{gamma, digamma, polygamma};
use rand_distr::{Distribution, Gamma};

fn main() {
    let a_true = 1.2; // shape: Set any value > 0.
    let b_true = 3.4; // scale: Set any value > 0.

    let gd = Gamma::new(a_true, b_true).unwrap();
    let x: Vec<f64> = gd.sample_iter(&mut rand::thread_rng()).take(10000).collect();

    let x_mean = x.iter().fold(0.0, |m, &e| m + e) / 10000.0;
    let xln_mean = x.iter().fold(0.0, |m, &e| m + e.ln()) / 10000.0;
    let x_mean_ln = x_mean.ln();

    let mut a = 0.5 / (x_mean_ln - xln_mean);
    let mut b = x_mean / a;

    let mut log_like_prev = std::f64::NEG_INFINITY;
    let mut log_like = calc_log_like(&x, a, b);
    
    while (log_like - log_like_prev).abs() >= 1e-7 {
        let num = xln_mean - digamma(a) - x_mean_ln + a.ln();
        let den = a * (1.0 - a * polygamma(a, 1));
        a = 1.0 / (1.0 / a + num / den);
        b = x_mean / a;

        log_like_prev = log_like;
        log_like = calc_log_like(&x, a, b);
    }

    println!("== True parameters ==");
    println!("a: shape = {}", a_true);
    println!("b: scale = {}\n", b_true);

    println!("== Estimated parameters ==");
    println!("a: shape = {}", a);
    println!("b: scale = {}", b);
}

/// The probability density function of the gamma distribution (f(x|a,b)).
/// a: shape
/// b: scale
fn gamma_pdf(x: f64, a: f64, b: f64) -> f64 {
    x.powf(a - 1.0) * (-x / b).exp() / (gamma(a) * b.powf(a))
}

/// Calculating log-likeliood of the dataset X in the estimated f(x|a,b).
fn calc_log_like(x: &Vec<f64>, a: f64, b: f64) -> f64 {
    x.iter().fold(0.0, |m, &e| m + gamma_pdf(e, a, b).ln())
}


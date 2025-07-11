use plotters::prelude::*;
use spfunc::zeta::zeta;
use cauchy::c64;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let t_start = 765.0;
    let t_end = 770.0;
    let t_step = 0.01;
    let epsilon = 1e-4;

    let root = BitMapBackend::new("zeta_critical_line.png", (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    // Scale plotting area
    let mut chart = ChartBuilder::on(&root)
        .caption("Rust: spfunc::zeta::zeta", ("sans-serif", 25))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(t_start..t_end, -1.5f64..3.0f64)?;

    chart.configure_mesh()
        .x_desc("0.5 + xi")
        .y_desc("y")
        .draw()?;

    let mut samples_re = Vec::new();
    let mut samples_im = Vec::new();
    let mut zero_spikes = Vec::new();

    for i in 0..=((t_end - t_start) / t_step) as usize {
        let t = t_start + i as f64 * t_step;
        let s = c64::new(0.5, t);
        let z = zeta(s);

        samples_re.push((t, z.re));
        samples_im.push((t, z.im));

        // Spike |ζ(s)| is nearly zero.
        let norm_sq = z.re * z.re + z.im * z.im;
        if norm_sq < epsilon {
            // zero_spikes.push(vec![(t, -30.0), (t, 30.0)]);
            zero_spikes.push((t, 1.0))
        } else {
            zero_spikes.push((t, 0.0))
        }
    }

    // Plot
   chart.draw_series(LineSeries::new(samples_re, &BLUE))?
        .label("real")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], &BLUE));

    chart.draw_series(LineSeries::new(samples_im, &RED))?
        .label("image")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], &RED));

    //for spike in zero_spikes {
    //    chart.draw_series(LineSeries::new(spike, &GREEN))?;
    //}
    chart.draw_series(LineSeries::new(zero_spikes, &GREEN))?
        .label("detector")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], &GREEN));


    chart.configure_series_labels().border_style(&BLACK).draw()?;

    println!("Complete: zeta_critical_line.png");
    Ok(())
}
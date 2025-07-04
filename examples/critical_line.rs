use plotters::prelude::*;
use spfunc::zeta::zeta;
use cauchy::c64;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let t_start = 765.0;
    let t_end = 770.0;
    let t_step = 0.01;
    let epsilon = 1e-6;

    let root = BitMapBackend::new("zeta_critical_line.png", (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    // Scale plotting area
    let mut chart = ChartBuilder::on(&root)
        .caption("Riemann Zeta on Critical Line", ("sans-serif", 25))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(t_start..t_end, -1.5f64..3.0f64)?;

    chart.configure_mesh()
        .x_desc("t in s = 0.5 + i·t")
        .y_desc("Re / Im / zero detector")
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
        if norm_sq < epsilon * epsilon {
            zero_spikes.push(vec![(t, -30.0), (t, 30.0)]);
        }
    }

    // Plot
   chart.draw_series(LineSeries::new(samples_re, &BLUE))?
        .label("Re(ζ(0.5 + i·t))")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], &BLUE));

    chart.draw_series(LineSeries::new(samples_im, &RED))?
        .label("Im(ζ(0.5 + i·t))")
        .legend(|(x, y)| PathElement::new([(x, y), (x + 20, y)], &RED));

    for spike in zero_spikes {
        chart.draw_series(LineSeries::new(spike, &GREEN))?;
    }

    chart.configure_series_labels().border_style(&BLACK).draw()?;

    println!("Complete: zeta_critical_line.png");
    Ok(())
}
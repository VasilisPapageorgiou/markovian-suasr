# Markovian SuASR Epidemic Model

This repository contains R implementations of the stochastic SuASR epidemic framework with external factor states.

- `suasr_variant1.R`: joint distribution of symptomatic and asymptomatic infections for Variant I.
- `suasr_variant1_maximum.R`: maximum infectious burden distribution for Variant I.
- `suasr_variant2_h1n1_joint.R`: joint infection-count distribution for the H1N1 Variant II example.
- `suasr_variant2_h1n1_maximum.R`: maximum infectious burden distribution for the H1N1 Variant II example.

The scripts use exact matrix-based computations and Gillespie simulation for numerical validation. Required R packages include `Matrix` and `ggplot2`.

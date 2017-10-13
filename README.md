# iMKT
Integrative McDonald and Kreitman test

R package to compute the iMKT. This is, to quantify the different selection regimes acting on any given region from polymorphism and divergence genomic data.

In summary, the diverse selection regimes are:

   - Strongly deleterious (d)
   - Weakly deleterious (b)
   - Neutral (f)
   - Adaptive fixations (α)

-----------------------------------------------------------

Functions included in the package:

   - fitMKmodel: two-step nls2() model fit at a given level of precision (res)
   - predictNLS: get a CI using Monte Carlo simulation based upon a fitted model
   - watchdog_asymptoticMK: check input data for asymptoticMK error handling
   - asymptoticMK: compute alpha asymptotic (uses fitMKmodel & predictNLS)
   - DGRP: compute alpha following DGRP methodology at 2 cutoffs: 0.05 & 0.1
   - iMK: compute d, b and f (uses asymptoticMK & DGRP)
   - plotIMK: horizontal plot of d, b and f fractions (uses previous)
   - plotALPHA: scatter plot of alpha values + shaded b fraction (uses previous)
   - plotDAF: double scatter plot of DAF distributions

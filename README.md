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
   - asymptoticMK: compute alpha asymptotic (uses fitMKmodel & predictNLS)
   - integrativeMK: compute d, b and f (uses asymptoticMK)

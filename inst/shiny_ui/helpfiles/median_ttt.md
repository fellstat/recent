# Median Time from Diagnosis to Treatment

For the RITA2 screening algorithm we require a curve for the probability that an individual 
is treated given they were diagnosed t days ago. Assuming that treatment cessation
during the first tau years is rare, the curve is identical to the time from
diagnosis to treatment curve. This application supports using an
exponential distribution to represent this.

Median time to treatment may be obtainable through administrative records.

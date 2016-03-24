Small script to compute the Genome Repeat Index (GRI) from raw sequencing reads.

The GRI is the percentage of all k-mer words which are predicted to be
repetitive. We predict that genomes with a lower GRI will be easier to assemble
than those with a higher one.

The k-mers counted for this statistic need to be reasonably large (>= 31bp) for this
metric to be meaningful, as we have found that shorter k-mer words can lead to
inaccuracies in the k-mer spectra produced.

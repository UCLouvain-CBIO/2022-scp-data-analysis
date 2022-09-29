This repository contains the scripts allowing the reproduce the
figures of *The current state of single-cell proteomics data analysis*
by Christophe Vanderaa and Laurent Gatto.

The main messages of the paper are:

- Data analysis plays a central role in science. The field of
  single-cell proteomics (SCP) has experienced an explosion of
  technical workflows, but its data analysis practices are lagging
  behind.
- Too little attention is given to computational analyses and the
  current efforts are not focused on quantitative data processing.
- Quantitative data processing workflows are a sequence of processing
  steps, and this sequence influences the analysis.
- Every SCP publication comes with a different computational pipeline
  that are fundamentally different.
- There is a lack of consensus, artefacts are propagated to the
  protein data, and current methods are taken from bulk proteomics
  (see figure 1).
- The solution is benchmarking, but it requires standardized
  computational tools and well-designed data.
- `scp` and `sceptre` are the only two computational environments that
  are flexible and designed for SCP.
- `scpdata` offers access to published data, but rigorous benchmark
  data is still missing.
- Replication of current SCP studies increases the trust in the
  results, highlights issues and remaining challenges, and provides
  didactic material.
- Different workflows lead to different results, but controled designs
  are needed for measuring the performance.
- Lessons learned from replication of SCP data analysis: (1) A good
  analysis requires good tools. (2) Consistent input formats
  facilitate data analysis. (3) Beware of confounding effects -
  because SCP requires large datasets, confounding effects become
  inevitable and need to be correctly modeled.
- We need to harmonize the workflows now before more analyses
  propagate current bad practices.

The results are based on the the SCP replications documented in
[SCP.Replication](https://uclouvain-cbio.github.io/SCP.replication/)
using the [scp](https://uclouvain-cbio.github.io/scp/) R/Bioconductor
package and the data in
[scpdata](https://www.bioconductor.org/packages/release/data/experiment/html/scpdata.html).

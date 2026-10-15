# CAMDAC 0.2.0

* Integrated RRBS and WGBS analysis under a single call to the `pipeline()` function.
* `CamConfig()` now resolves `min_mapq` per protocol: 0 for RRBS, 1 for WGBS.
  RRBS alignments are validated on read start site rather than mapping quality,
  so filtering on MAPQ would bias against reads carrying the alternate allele.
* `CamConfig()` default `min_cov` is 3, matching the RRBS and WGBS pipelines
  that preceded this release.
* `CamConfig()` gained `het_baf`, the normal BAF window for heterozygous SNP
  selection, defaulting to 0.15/0.85 for RRBS and 0.08/0.92 for WGBS.
* Added option for RRBS paired end analysis
* Upgraded ASCAT to version 3.2.0 (commit: 44ddd3080723a2c3640d1cfead13437a093c21d1) 

# CAMDAC 0.1.0

* Minor documentation updates.

# Cancer Hotspots BED File

This repository contains a BED file derived from [cancerhotspots.org](https://www.cancerhotspots.org/) data, reformatted to define genomic regions relevant for NGS analysis.

## Overview

I used the cancerhotspots.org API to retrieve hotspot data and converted it into BED format. Since many hotspot entries are single base positions, I expanded these to cover:

- Full codons for single base substitutions  
- Full region spans for insertions/deletions  

Splice variants have remained unchanged from the API output.

The aim was to represent regions in a way that would be more useful for a pipeline to identify variants in hotspots, to aid S-VIG UK code automation. This resource may be helpful for NGS labs wanting to incorporate these hotspot regions into their pipelines.

The transcripts provided by cancerhotspots.org API were not always accurate, so in the generation of this bed file, a check of the residue against the transcript's protein sequence was performed, in the case of mismatch, the mane_select transcript was tried and residue check also performed. If neither the transcript provided or mane_select transcript were a correct match, then then hotspot was attempted to be resolved manually.

## Use at Your Own Risk

Please validate the regions yourself if you're using this in production or clinical pipelines.
- Not all entries have been manually reviewed—some regions may be inaccurate or imprecise.
- BED coordinates are based on the best interpretation of the API output, but edge cases may need review.
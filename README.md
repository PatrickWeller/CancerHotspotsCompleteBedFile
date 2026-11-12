# Cancer Hotspots BED File

This repository contains a BED file derived from [cancerhotspots.org](https://www.cancerhotspots.org/) data, reformatted to define genomic regions relevant for NGS analysis.

## Overview

I used the cancerhotspots.org API to retrieve hotspot data and converted it into BED format. Since many hotspot entries are single base positions, I expanded these to cover:

- Full codons for single base substitutions  
- Full region spans for insertions/deletions  

Splice site regions are taken unchanged from the 2022 [sigven/cancerHotspots](https://github.com/sigven/cancerHotspots) BED file, since the API gives no genomic position for them.

The aim was to represent regions in a way that would be more useful for a pipeline to identify variants in hotspots, to aid S-VIG UK code automation. This resource may be helpful for NGS labs wanting to incorporate these hotspot regions into their pipelines.

The transcripts provided by the cancerhotspots.org API were not always accurate. So for single-residue hotspots, the residue was checked against the transcript's protein sequence. If it didn't match, or no transcript was given, the MANE Select transcript was tried with the same check. If neither matched, the hotspot was resolved manually in IGV. Indel hotspots were only mapped on the transcript provided, and resolved manually when that failed.

Regions were built on GRCh38 and lifted over to GRCh37 with CrossMap.

## Use at Your Own Risk

Please validate the regions yourself if you're using this in production or clinical pipelines.
- Not all entries have been manually reviewed—some regions may be inaccurate or imprecise.
- BED coordinates are based on the best interpretation of the API output, but edge cases may need review.

## Reproducing the BED files

The scripts and inputs used to build these files are in [`pipeline/`](pipeline/). See [`pipeline/README.md`](pipeline/README.md) for how to run them.

# Pipeline

`build.py` produces `output/cancer_hotspots_GRCh38.bed` and `output/cancer_hotspots_GRCh37.bed`.

```bash
conda activate crossMap        # has python + requests + CrossMap
python build.py                # download fresh data, then build (~2–3 min, mostly waiting on Ensembl)
python build.py --cached       # reuse the last run's downloads (~2 s), e.g. after editing manual/*.bed
```

## Files

| File | What |
|------|------|
| `build.py` | Downloads the data and runs steps 1–9 in order |
| `ensembl.py` | Ensembl downloads (batch requests, with retries) and protein position → genome coordinates |
| `snv_hotspots.py` | Step 1 |
| `indel_hotspots.py` | Step 5 |
| `inputs/mane_select_cache.csv` | NCBI MANE summary: each gene's MANE Select transcript, used as the fallback in step 1 |
| `inputs/20220125_hotspots_all_38.bed` | 2022-01-25 cancerHotspots BED (sigven/cancerHotspots); source of the splice regions |
| `inputs/hg38ToHg19.over.chain.gz` | UCSC chain file for step 9 |
| `manual/snvs_manual.bed` | SNV (and splice) hotspots that failed automatically, resolved by hand in IGV |
| `manual/indels_manual.bed` | Indel hotspots that failed automatically, resolved by hand in IGV |

## Downloads

1. All hotspots from `cancerhotspots.org/api/hotspots/single?version=v3`: 1274 `single residue` and 55 `in-frame indel` (as of November 2026). Without `?version=v3` the API returns 1165 (1110 + 55), which is what the 2025 BEDs were built from. Most of the extra v3 hotspots come with no transcript, so step 1 uses MANE Select for them.
2. From Ensembl (GRCh38), every transcript those hotspots might use: each hotspot's own transcript plus its gene's MANE Select transcript. This includes their exons and protein sequences, fetched with Ensembl's batch endpoints, so it's a handful of requests.

Both are saved in `output/` (`0_cancerhotspots_api.json`, `0_ensembl_transcripts.json`), as a record of what the build used and so `--cached` can rebuild without downloading. Everything after this runs locally: a codon's position in the genome is worked out from the transcript's exons.

If Ensembl or cancerhotspots.org is down, the download retries for about 2 minutes, then stops with an error. Try again later.

## Steps

| # | What | Output |
|---|------|--------|
| 1 | **SNV hotspots → codons.** Try the hotspot's transcript (if it has one), then MANE Select; use the first whose protein has the hotspot's amino acid at that position | `1_snvs_auto.bed`, `1_snv_results.tsv` |
| 2 | **Drop splice rows** from step 1 (residue `X`, e.g. `PTEN_X70`); a codon is the wrong region for a splice site | `2_snvs_no_splice.bed` |
| 3 | **Hand-curated SNVs** from `manual/snvs_manual.bed` | |
| 4 | **Splice hotspots**: the `_sice_` rows of the 2022 cancerHotspots BED; coordinates unchanged, `chr` prefix and extra columns removed | `4_splice.bed` |
| 5 | **Indel hotspots → genomic spans**, on the hotspot's own transcript only. Indels have no amino acid to check, so switching to MANE wouldn't be safe | `5_indels_auto.bed`, `5_indel_results.tsv` |
| 6 | **Hand-curated indels** from `manual/indels_manual.bed` | |
| 7 | **Check every hotspot has a region** (below). If any don't, the build **stops** | `7_needs_manual_curation.tsv` |
| 8 | **Merge** 2 + 3 + 4 + 5 + 6 and sort by chromosome, then start | `cancer_hotspots_GRCh38.bed` |
| 9 | **Lift over to GRCh37** with CrossMap | `cancer_hotspots_GRCh37.bed` |

`1_snv_results.tsv` and `5_indel_results.tsv` have one row per hotspot: whether it mapped, which transcript was used, and why it failed or why MANE Select was used.

## Step 7: manual curation

The manual files are hotspots that failed automatically and were checked by hand in IGV. They're only ever read, and every line in them goes into the BED. Each line is matched to its hotspot by its name, `GENE_RESIDUE`. Notes in brackets are ignored, so `MYD88_L265(L252)` matches MYD88 L265 and `NOTCH1_356-359` matches that indel range.

- **Hotspot with no automatic or manual region**: listed in `output/7_needs_manual_curation.tsv` with the reason, and the build stops before merging. Check it in IGV, add it to the right manual file, then run `python build.py --cached`.
- **Hotspot with a manual entry that now also maps automatically**: your entry is kept and the automatic lines are dropped. It prints `NOTE` if the coordinates are the same, or `CHECK` if they differ, so you can look again.
- **Manual entry that matches no hotspot**: kept, with a `NOTE`.

`--allow-incomplete` builds the BED anyway, leaving the unresolved hotspots out.

### Splice hotspots

Splice hotspots (residue `X`, e.g. `TP53_X187`) aren't mapped by this code. Their regions come from the 2022 sigven/cancerHotspots BED: single intronic bases (the canonical ±1/±2 splice positions, plus any other mutated bases) next to the exon edge at that codon. A splice hotspot missing from that BED is flagged in step 7. Look at the hotspot's codon in IGV to find the splice site, and add it to `manual/snvs_manual.bed` named `GENE_X###`.

`RECQL4_X766` is a special case. Its transcript, ENST00000428558 (GRCh37 only), had a fake 1-bp intron to work around an extra C in the GRCh37 reference, so mutations at GRCh37 chr8:145738769 were called "splice". In GRCh38 that base is in codon 766 of the MANE transcript, so the manual entry is that codon.

## Before replacing the committed BEDs

Diff `output/cancer_hotspots_GRCh3*.bed` against `../cancer_hotspots_GRCh3*.bed` and check that every difference makes sense.

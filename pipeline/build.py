"""Build cancer_hotspots_GRCh38.bed and cancer_hotspots_GRCh37.bed from cancerhotspots.org.

First, download the hotspots (cancerhotspots.org) and every transcript they need (Ensembl).
Then (each step's output is saved in output/ with its step number):
  1. Map SNV (single residue) hotspots to their codons
  2. Remove the splice rows (residue X) from step 1; the codon mapping gets these wrong
  3. Load hand-curated SNV regions                         manual/snvs_manual.bed
  4. Load splice hotspots from the 2022 cancerHotspots BED  inputs/20220125_hotspots_all_38.bed
  5. Map in-frame indel hotspots to genomic spans
  6. Load hand-curated indel regions                       manual/indels_manual.bed
  7. Check every hotspot has a region, from steps 1-6. Any that don't are listed in
     output/7_needs_manual_curation.tsv and the build STOPS, so you can check them in IGV and
     add them to manual/*.bed. Manual entries always win over automatic ones.
  8. Merge and sort                                        output/cancer_hotspots_GRCh38.bed
  9. Lift over to GRCh37 with CrossMap                     output/cancer_hotspots_GRCh37.bed

Usage (from this folder, in the crossMap conda env, which has requests and CrossMap):
  python build.py            # download fresh data, then build
  python build.py --cached   # reuse the downloads from the last run (e.g. after editing manual/*.bed)
"""
import argparse
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import requests

import ensembl
from indel_hotspots import map_indel_hotspots
from snv_hotspots import load_mane_select, map_snv_hotspots

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
MANUAL = HERE / "manual"
HOTSPOTS_API = "https://www.cancerhotspots.org/api/hotspots/single?version=v3"


# ---------------------------------------------------------------- file helpers

def read_bed(path):
	"""BED lines, skipping blank lines."""
	with open(path) as f:
		return [line.rstrip("\n") for line in f if line.strip()]


def write_lines(path, lines, header=None):
	with open(path, "w") as f:
		if header:
			f.write(header + "\n")
		for line in lines:
			f.write(line + "\n")
	print(f"  -> {path.name} ({len(lines)} lines)")


def sort_bed(lines):
	"""Sort by chromosome (as text: 1, 10, 11, ..., 2, ..., X), then start position."""
	def key(line):
		chrom, start = line.split("\t")[:2]
		return (chrom, int(start), line)
	return sorted(lines, key=key)


def bed_name_key(bed_line):
	"""BED name -> (gene, residue), ignoring notes in brackets.
	'MYD88_L265(L252)' -> ('MYD88', 'L265');  'NOTCH1_356-359' -> ('NOTCH1', '356-359')"""
	name = bed_line.split("\t")[3]
	gene, residue = re.sub(r"\(.*?\)", "", name).rsplit("_", 1)
	return gene, residue


def group_by_hotspot(bed_lines):
	"""{(gene, residue): [bed lines]}. A codon or indel split by an intron has several lines."""
	groups = defaultdict(list)
	for line in bed_lines:
		groups[bed_name_key(line)].append(line)
	return groups


def coordinates(bed_lines):
	return sorted(tuple(line.split("\t")[:3]) for line in bed_lines)


# ---------------------------------------------------------------- downloading

def download(out, use_cached):
	"""Returns (hotspots, transcripts), downloading them or reusing the copies saved in out/."""
	hotspots_file = out / "0_cancerhotspots_api.json"
	transcripts_file = out / "0_ensembl_transcripts.json"
	mane = load_mane_select(INPUTS / "mane_select_cache.csv")

	if use_cached:
		if not (hotspots_file.exists() and transcripts_file.exists()):
			sys.exit(f"--cached needs {hotspots_file.name} and {transcripts_file.name} in {out} from an earlier run")
		print(f"Using the downloads saved in {out}")
		hotspots = json.loads(hotspots_file.read_text())
		records = json.loads(transcripts_file.read_text())
	else:
		try:
			print("Downloading hotspots from cancerhotspots.org")
			hotspots = ensembl.get_json(HOTSPOTS_API)
			# Every transcript any hotspot might use: its own, plus its gene's MANE Select
			wanted = {h["transcriptId"] for h in hotspots} | {mane[h["hugoSymbol"]] for h in hotspots if h["hugoSymbol"] in mane}
			print(f"Downloading {len(wanted)} transcripts and their proteins from Ensembl (a few batch requests)")
			records = ensembl.download_transcripts(wanted)
		except requests.RequestException as e:
			sys.exit(f"\nDownload failed, even after retrying: {e}\nEnsembl or cancerhotspots.org is probably busy; try again later.")
		# Saved so that --cached can rebuild offline, and as a record of what this build used
		hotspots_file.write_text(json.dumps(hotspots, indent=1))
		transcripts_file.write_text(json.dumps(records, indent=1))

	print(f"  {len(hotspots)} hotspots, {len(records)} transcripts")
	transcripts = {tid: ensembl.Transcript(tid, rec) for tid, rec in records.items()}
	return hotspots, transcripts, mane


# ---------------------------------------------------------------- step 7

def check_curation(hotspots, auto_lines, splice_lines, manual_lines, failed):
	"""Compare every hotspot from the API with what steps 1-6 produced.

	failed: {(gene, residue): (transcript, reason)} for hotspots steps 1 and 5 couldn't map.
	Returns (todo, replaced_by_manual):
	  todo                rows (type, gene, residue, transcript, reason) for hotspots with no region
	  replaced_by_manual  hotspots with a manual entry that now also map automatically; the manual
	                      entry wins (you checked it in IGV) and the automatic lines are dropped
	"""
	auto = group_by_hotspot(auto_lines)
	manual = group_by_hotspot(manual_lines)
	splice_names = {line.split("\t")[3] for line in splice_lines}

	todo = []
	api_keys = set()
	for hotspot in hotspots:
		gene, residue, transcript = hotspot["hugoSymbol"], hotspot["residue"], hotspot["transcriptId"] or "-"

		if hotspot["type"] == "in-frame indel":
			kind = "indel"
			first, _, last = residue.partition("-")
			residue = f"{first}-{last or first}"
		elif residue.startswith("X"):
			kind = "splice"   # covered by the 2022 splice BED (step 4), not by step 1
		else:
			kind = "SNV"

		key = (gene, residue)
		api_keys.add(key)
		if key in manual:
			continue
		if kind == "splice":
			if not any(name.startswith(f"{gene}_aa{residue}_sice_") for name in splice_names):
				todo.append((kind, gene, residue, transcript, "not in the 2022 splice BED"))
		elif key not in auto:
			reason = failed.get(key, (transcript, "no automatic result"))[1]
			todo.append((kind, gene, residue, transcript, reason))

	replaced_by_manual = set()
	for key, manual_entry in manual.items():
		label = f"{key[0]} {key[1]}"
		if key not in api_keys:
			print(f"  NOTE: manual entry {label} doesn't match any hotspot from the API (kept anyway)")
		elif key in auto:
			replaced_by_manual.add(key)
			if coordinates(auto[key]) == coordinates(manual_entry):
				print(f"  NOTE: {label} now maps automatically, the same as your manual entry (duplicate dropped)")
			else:
				print(f"  CHECK: {label} now maps automatically, but differently from your manual entry (manual kept)")
				print(f"         automatic: {coordinates(auto[key])}")
				print(f"         manual:    {coordinates(manual_entry)}")

	return todo, replaced_by_manual


# ---------------------------------------------------------------- main

def step(n, title):
	print(f"\n== Step {n}: {title}")


def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("--cached", action="store_true",
		help="reuse the downloads saved by the last run instead of downloading again")
	parser.add_argument("--allow-incomplete", action="store_true",
		help="build the BED even if some hotspots still need curating (they'll be missing from it)")
	parser.add_argument("--output", type=Path, default=HERE / "output", help="output folder (default: pipeline/output)")
	parser.add_argument("--crossmap", default=shutil.which("CrossMap"),
		help="path to the CrossMap program (default: CrossMap on PATH)")
	args = parser.parse_args()

	if not args.crossmap:
		sys.exit("CrossMap not found. Activate the crossMap conda env, or pass --crossmap /path/to/CrossMap")
	out = args.output
	out.mkdir(parents=True, exist_ok=True)
	results_header = "gene\tresidue\tresult\ttranscript\tdetail"

	hotspots, transcripts, mane = download(out, args.cached)

	step(1, "map SNV hotspots to codons")
	snvs_auto, snv_results = map_snv_hotspots(hotspots, transcripts, mane)
	write_lines(out / "1_snvs_auto.bed", snvs_auto)
	write_lines(out / "1_snv_results.tsv", ["\t".join(r) for r in snv_results], results_header)

	step(2, "remove splice (residue X) rows from step 1")
	snvs_no_splice = [line for line in snvs_auto if not re.search(r"_X\d+$", line.split("\t")[3])]
	print(f"  removed {len(snvs_auto) - len(snvs_no_splice)} rows")
	write_lines(out / "2_snvs_no_splice.bed", snvs_no_splice)

	step(3, "hand-curated SNV regions")
	snvs_manual = read_bed(MANUAL / "snvs_manual.bed")
	print(f"  {len(snvs_manual)} regions from manual/snvs_manual.bed")

	step(4, "splice hotspots from the 2022 cancerHotspots BED")
	splice = []
	for line in read_bed(INPUTS / "20220125_hotspots_all_38.bed"):
		if "_sice_" in line:
			chrom, start, end, name = line.split("\t")[:4]
			splice.append(f"{chrom.replace('chr', '')}\t{start}\t{end}\t{name}")
	write_lines(out / "4_splice.bed", splice)

	step(5, "map in-frame indel hotspots to genomic spans")
	indels_auto, indel_results = map_indel_hotspots(hotspots, transcripts)
	write_lines(out / "5_indels_auto.bed", indels_auto)
	write_lines(out / "5_indel_results.tsv", ["\t".join(r) for r in indel_results], results_header)

	step(6, "hand-curated indel regions")
	indels_manual = read_bed(MANUAL / "indels_manual.bed")
	print(f"  {len(indels_manual)} regions from manual/indels_manual.bed")

	step(7, "check every hotspot has a region")
	failed = {(gene, residue): (transcript, detail)
		for gene, residue, result, transcript, detail in snv_results + indel_results if result == "failed"}
	todo, replaced_by_manual = check_curation(
		hotspots, snvs_no_splice + indels_auto, splice, snvs_manual + indels_manual, failed)

	todo_file = out / "7_needs_manual_curation.tsv"
	if todo_file.exists():
		todo_file.unlink()   # from an earlier run
	if not todo:
		print("  all hotspots covered")
	else:
		write_lines(todo_file, ["\t".join(row) for row in todo], "type\tgene\tresidue\ttranscript\treason")
		print(f"\n  {len(todo)} hotspots have no region and need curating by hand (check in IGV):")
		for kind, gene, residue, transcript, reason in todo:
			print(f"    {kind:6} {gene:10} {residue:10} {transcript:16} {reason}")
		print("  Add them to manual/snvs_manual.bed or manual/indels_manual.bed, named GENE_RESIDUE")
		print("  (notes in brackets are fine, e.g. MYD88_L265(L252)), then run: python build.py --cached")
		if not args.allow_incomplete:
			sys.exit("\nStopped before merging. Use --allow-incomplete to build the BED without them.")
		print("  --allow-incomplete: carrying on; these hotspots will be MISSING from the BED")

	step(8, "merge and sort")
	auto_kept = [line for line in snvs_no_splice + indels_auto if bed_name_key(line) not in replaced_by_manual]
	grch38 = sort_bed(auto_kept + snvs_manual + splice + indels_manual)
	grch38_file = out / "cancer_hotspots_GRCh38.bed"
	write_lines(grch38_file, grch38)

	step(9, "lift over to GRCh37 (CrossMap)")
	lifted = out / "9_lifted.tmp.bed"
	unmapped = out / "9_lifted.tmp.bed.unmap"   # CrossMap writes regions it can't lift here
	crossmap = subprocess.run(
		[args.crossmap, "bed", str(INPUTS / "hg38ToHg19.over.chain.gz"), str(grch38_file), str(lifted)],
		stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	if crossmap.returncode != 0:
		sys.exit(f"CrossMap failed:\n{crossmap.stderr}")
	grch37 = sort_bed([line.replace("chr", "", 1) if line.startswith("chr") else line for line in read_bed(lifted)])
	lifted.unlink()
	write_lines(out / "cancer_hotspots_GRCh37.bed", grch37)
	if unmapped.exists() and read_bed(unmapped):
		unmapped.rename(out / "9_liftover_failed.bed")
		print("  WARNING: some regions couldn't be lifted over; see 9_liftover_failed.bed")
	elif unmapped.exists():
		unmapped.unlink()

	print("\nDone. Diff output/cancer_hotspots_GRCh3*.bed against the committed files before replacing them.")


if __name__ == "__main__":
	main()

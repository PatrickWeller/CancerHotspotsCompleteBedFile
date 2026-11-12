"""Step 1: map single-residue (SNV) hotspots to their genomic codons.

For each hotspot, try its cancerhotspots.org transcript, then the gene's MANE Select transcript.
The first transcript whose protein has the hotspot's amino acid at that position is used, and the
hotspot's codon is located from that transcript's exons. If neither works, the hotspot fails and
step 7 lists it for curating by hand.
"""


def load_mane_select(mane_file):
	"""{gene symbol: MANE Select Ensembl transcript ID (without version)} from the NCBI MANE summary."""
	mane = {}
	with open(mane_file) as f:
		header = f.readline().lstrip("#").rstrip("\n").split("\t")
		for line in f:
			row = dict(zip(header, line.rstrip("\n").split("\t")))
			if row["MANE_status"] == "MANE Select":
				mane[row["symbol"]] = row["Ensembl_nuc"].split(".")[0]
	return mane


def residue_problem(transcript, amino_acid, position):
	"""Why this transcript can't be used for the hotspot, or None if it can.
	Splice hotspots (amino acid X) can't be checked, so only their position is."""
	if transcript is None:
		return "not in Ensembl GRCh38"
	if not transcript.protein:
		return "no protein sequence"
	if position > len(transcript.protein):
		return f"position {position} is beyond the protein's length ({len(transcript.protein)})"
	found = transcript.protein[position - 1]
	if amino_acid != "X" and found != amino_acid:
		return f"residue mismatch: {found} at {position}, not {amino_acid}"
	return None


def map_snv_hotspots(hotspots, transcripts, mane):
	"""Returns (bed_lines, results), where results has one row per hotspot:
	(gene, residue, "mapped" or "failed", transcript, detail)."""
	bed_lines, results = [], []

	for hotspot in hotspots:
		if hotspot["type"] != "single residue":
			continue
		gene, residue = hotspot["hugoSymbol"], hotspot["residue"]   # e.g. "NRAS", "Q61"
		amino_acid, position = residue[0], int(residue[1:])
		given = hotspot["transcriptId"]

		candidates = [given] if given else []   # cancerhotspots.org gives no transcript for a few (PAK7 M173)
		if mane.get(gene) and mane[gene] != given:
			candidates.append(mane[gene])

		problems = [] if given else ["no transcript given by cancerhotspots.org"]
		for transcript_id in candidates:
			transcript = transcripts.get(transcript_id)
			problem = residue_problem(transcript, amino_acid, position)
			if problem:
				problems.append(f"{transcript_id}: {problem}")
				continue
			for chrom, start, end in transcript.codon_regions(position):
				bed_lines.append(f"{chrom}\t{start}\t{end}\t{gene}_{residue}")
			detail = f"MANE Select, because {problems[0]}" if problems else ""
			results.append((gene, residue, "mapped", transcript_id, detail))
			break
		else:
			results.append((gene, residue, "failed", given or "-", " | ".join(problems)))

	mapped = sum(1 for r in results if r[2] == "mapped")
	print(f"  {mapped} of {len(results)} SNV hotspots mapped; {len(results) - mapped} failed")
	return bed_lines, results

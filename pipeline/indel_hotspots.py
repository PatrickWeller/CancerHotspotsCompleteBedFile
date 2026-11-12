"""Step 5: map in-frame indel hotspots to genomic spans.

An indel hotspot's residue is a position ("546") or a range ("65-77"). It's located on the
hotspot's own transcript only: indels have no reference amino acid to check, so there's no safe
way to switch to the MANE Select transcript. If the transcript isn't usable, the hotspot fails and
step 7 lists it for curating by hand.
"""


def map_indel_hotspots(hotspots, transcripts):
	"""Returns (bed_lines, results), where results has one row per hotspot:
	(gene, residue, "mapped" or "failed", transcript, detail)."""
	bed_lines, results = [], []

	for hotspot in sorted(hotspots, key=lambda h: h["hugoSymbol"]):
		if hotspot["type"] != "in-frame indel":
			continue
		gene, transcript_id = hotspot["hugoSymbol"], hotspot["transcriptId"]
		first, _, last = hotspot["residue"].partition("-")
		residue = f"{first}-{last or first}"   # named like the BED: "546" -> "546-546"

		transcript = transcripts.get(transcript_id)
		if transcript is None:
			problem = "not in Ensembl GRCh38"
		elif not transcript.protein:
			problem = "no protein sequence"
		elif int(last or first) > len(transcript.protein):
			problem = f"residues beyond the protein's length ({len(transcript.protein)})"
		else:
			problem = None

		if problem:
			results.append((gene, residue, "failed", transcript_id, problem))
			continue
		for chrom, start, end in transcript.codon_regions(int(first), int(last or first)):
			bed_lines.append(f"{chrom}\t{start}\t{end}\t{gene}_{residue}")
		results.append((gene, residue, "mapped", transcript_id, ""))

	mapped = sum(1 for r in results if r[2] == "mapped")
	print(f"  {mapped} of {len(results)} indel hotspots mapped; {len(results) - mapped} failed")
	return bed_lines, results

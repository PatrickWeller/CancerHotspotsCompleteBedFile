"""Download transcripts and protein sequences from Ensembl, and turn protein positions into
genome coordinates.

All the transcripts a run needs come down in a few batch requests (instead of several requests
per hotspot), and codon positions are then worked out locally from each transcript's exons.
"""
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ENSEMBL = "https://rest.ensembl.org"   # GRCh38
TIMEOUT = 120         # seconds per request; Ensembl can be very slow
LOOKUP_BATCH = 200    # transcripts per request (Ensembl allows up to 1000)
SEQUENCE_BATCH = 50   # protein sequences per request (Ensembl's maximum)

# Retry when Ensembl is busy or down; give up (raise) after about 2 minutes of trying.
_retry = Retry(
	total=6,
	backoff_factor=2,
	status_forcelist=[429, 500, 502, 503, 504],
	respect_retry_after_header=True,
	allowed_methods=["GET", "POST"],
)
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=_retry))
session.headers.update({"Content-Type": "application/json", "Accept": "application/json"})


def get_json(url):
	response = session.get(url, timeout=TIMEOUT)
	response.raise_for_status()
	return response.json()


def post_json(url, body):
	response = session.post(url, json=body, timeout=TIMEOUT)
	response.raise_for_status()
	return response.json()


def in_chunks(items, size):
	for i in range(0, len(items), size):
		yield items[i:i + size]


# ---------------------------------------------------------------- downloading

def download_transcripts(transcript_ids):
	"""{transcript_id: Ensembl lookup record (with exons and translation)}.
	IDs Ensembl doesn't know (e.g. retired transcripts) are simply missing from the result."""
	transcript_ids = sorted({tid for tid in transcript_ids if tid})
	records = {}
	for batch in in_chunks(transcript_ids, LOOKUP_BATCH):
		found = post_json(f"{ENSEMBL}/lookup/id?expand=1", {"ids": batch})
		records.update({tid: rec for tid, rec in found.items() if rec})

	protein_ids = [rec["Translation"]["id"] for rec in records.values() if "Translation" in rec]
	sequences = {}
	for batch in in_chunks(protein_ids, SEQUENCE_BATCH):
		for entry in post_json(f"{ENSEMBL}/sequence/id", {"ids": batch}):
			sequences[entry["id"]] = entry["seq"]

	# Keep just what's needed, so the saved download stays small and readable
	slim = {}
	for tid, rec in records.items():
		translation = rec.get("Translation")
		slim[tid] = {
			"gene": rec.get("display_name"),
			"chrom": rec["seq_region_name"],
			"strand": rec["strand"],
			"exons": sorted((e["start"], e["end"]) for e in rec["Exon"]),
			"cds_start": translation["start"] if translation else None,
			"cds_end": translation["end"] if translation else None,
			"protein": sequences.get(translation["id"]) if translation else None,
		}
	return slim


# ---------------------------------------------------------------- protein -> genome

class Transcript:
	"""One transcript from download_transcripts(), able to locate codons in the genome."""

	def __init__(self, transcript_id, record):
		self.id = transcript_id
		self.chrom = record["chrom"]
		self.strand = record["strand"]
		self.protein = record["protein"]   # amino acid sequence, or None if non-coding
		self.cds_pieces = self._cds_pieces(record)

	@staticmethod
	def _cds_pieces(record):
		"""Coding parts of each exon as (start, end) 1-based genome positions, in the order they're
		read (5' to 3' of the transcript, so descending on the minus strand)."""
		if record["cds_start"] is None:
			return []
		pieces = []
		for exon_start, exon_end in record["exons"]:
			start = max(exon_start, record["cds_start"])
			end = min(exon_end, record["cds_end"])
			if start <= end:
				pieces.append((start, end))
		return pieces if record["strand"] == 1 else pieces[::-1]

	def codon_regions(self, first_residue, last_residue=None):
		"""Genome regions covering codons first_residue..last_residue (1-based amino acid positions),
		as BED tuples (chrom, start, end). A codon split by an intron gives one region per exon."""
		last_residue = last_residue or first_residue
		want_from = 3 * (first_residue - 1)   # 0-based offsets into the coding sequence
		want_to = 3 * last_residue            # (half-open)

		regions = []
		offset = 0   # coding bases before this exon piece
		for start, end in self.cds_pieces:
			length = end - start + 1
			lo = max(want_from, offset) - offset   # overlap, relative to the start of this piece
			hi = min(want_to, offset + length) - offset
			if lo < hi:
				if self.strand == 1:
					regions.append((self.chrom, start - 1 + lo, start - 1 + hi))
				else:
					regions.append((self.chrom, end - hi, end - lo))
			offset += length
		return regions

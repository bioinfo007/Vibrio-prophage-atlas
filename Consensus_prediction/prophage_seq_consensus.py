#!/usr/bin/env python3
"""
prophage_seq_consensus.py
============================
Consensus prophage pipeline.

ALGORITHM----AQIB JAVAID
----------
For each contig:

  1. Write ALL geNomad + VIBRANT fragments for this contig into a single
     multi-FASTA query file.  Run ONE minimap2 call (batch mode).
     Demultiplex hits by query ID.

  2. Filter each fragment's hits:
       identity  = 1 - NM/aln_len >= --min_identity   (default 0.99)
       coverage  = query_aligned / frag_len >= --min_self_coverage
       Deduplicate placements: two hits are the "same" only if they
       overlap by >= 90% of the shorter interval (tightened from 50%).
       0 passing distinct placements → failed
       >1 passing distinct placements → ambiguous (→ ambiguous.fasta)
       1 passing → anchored; record contig [start, end)

  3. Pre-merge same-tool anchored intervals with gap <= --merge_gap.
     This collapses split predictions from one tool before cross-tool
     comparison.

  4. Bipartite matching:
       Left nodes  = merged geNomad intervals
       Right nodes = merged VIBRANT intervals
       Edge weight = min(overlap_frac_g, overlap_frac_v) using TRUE union
       Edge eligible: overlap_bp >= --min_overlap_bp AND both fractions
                      >= --min_reciprocal_overlap
       Find maximum-weight bipartite matching (scipy linear_sum_assignment).
       Each geNomad interval → at most one VIBRANT interval.

  5. Matched pairs → consensus groups.
     Unmatched intervals → genomad_only / vibrant_only groups.

  6. For each consensus group:
       - outer sequence: contig_seq[min_start : max_end]  (union envelope)
       - inner coords:   [max(g_start, v_start) : min(g_end, v_end)]
         both reported; inner coords in merge_report.tsv for stricter use.

OUTPUT FILES
-------------
<outdir>/consensus.fasta
<outdir>/genomad_only.fasta
<outdir>/vibrant_only.fasta
<outdir>/ambiguous.fasta
<outdir>/anchor_report.tsv    — per-fragment; includes identity + candidate_coords
<outdir>/merge_report.tsv     — per-group; includes overlap metrics, boundary
                                 disagreement, inner/outer coords, support_score
<outdir>/contig_summary.tsv

USAGE
------
python prophage_seq_consensusV5.py \\
    --genomad_fasta   Genomad_predicted_phages.fasta \\
    --vibrant_fasta   VIBRANT_predicted_phages.fasta \\
    --assembly_fasta  combined_vibrio_genomes.fasta \\
    --outdir          consensus_output/ \\
    [--min_self_coverage      0.99]  \\
    [--min_identity           0.99]  \\
    [--min_overlap_bp         1000]  \\
    [--min_reciprocal_overlap 0.50]  \\
    [--merge_gap              0]     \\
    [--no_strict_consensus]          \\
    [--workers                4]     \\
    [--minimap_threads        2]     \\
    [--keep_tmp]

DEPENDENCIES
-------------
- minimap2      (must be in PATH)
- BioPython
- scipy         (for linear_sum_assignment — bipartite matching)
- pyfaidx       (optional; strongly recommended for large assemblies)
"""

import argparse
import logging
import multiprocessing
import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ImportError:
    sys.exit("BioPython is required: pip install biopython")

try:
    from scipy.optimize import linear_sum_assignment
    import numpy as np
except ImportError:
    sys.exit("scipy and numpy are required: pip install scipy numpy")

try:
    from pyfaidx import Fasta
except ImportError:
    try:
        from pyfaidx import Faidx as Fasta
    except ImportError:
        Fasta = None


logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Header parsing
# ─────────────────────────────────────────────────────────────────────────────

def parse_genomad_header(header: str) -> dict:
    """
    SAMN02604040_1|provirus_1114776_1137807  -> provirus on contig SAMN02604040_1
    SAMN02604040_3                            -> whole contig
    Hint coords stored for sanity-check only; never used in final output.
    """
    header = header.lstrip(">").strip()
    info = {
        "original_header": header,
        "tool":            "genomad",
        "contig":          None,
        "hint_start":      None,
        "hint_end":        None,
    }
    if "|provirus_" in header:
        contig_part, coord_part = header.split("|provirus_", 1)
        info["contig"] = contig_part
        parts = coord_part.split("_")
        if len(parts) == 2:
            try:
                info["hint_start"] = int(parts[0])
                info["hint_end"]   = int(parts[1])
            except ValueError:
                pass
    else:
        info["contig"] = header
    return info


def parse_vibrant_header(header: str) -> dict:
    """
    SAMN02604040_1_fragment_3  -> contig SAMN02604040_1
    SAMN02604040_3             -> whole contig
    """
    header = header.lstrip(">").strip()
    info = {
        "original_header": header,
        "tool":            "vibrant",
        "contig":          None,
        "hint_start":      None,
        "hint_end":        None,
    }
    frag_match = re.search(r"_fragment_(\d+)$", header)
    if frag_match:
        info["contig"] = header[: frag_match.start()]
    else:
        info["contig"] = header
    return info


# ─────────────────────────────────────────────────────────────────────────────
# FASTA loading  (metadata in parallel dict — no monkey-patching)
# ─────────────────────────────────────────────────────────────────────────────

def load_fasta_by_contig(fasta_path: str, tool: str):
    """
    Returns:
      by_contig : { contig_id: [SeqRecord, ...] }
      meta      : { record_id: info_dict }
    """
    by_contig = defaultdict(list)
    meta      = {}
    parse_fn  = parse_genomad_header if tool == "genomad" else parse_vibrant_header

    for rec in SeqIO.parse(fasta_path, "fasta"):
        info = parse_fn(rec.id)
        meta[rec.id] = info
        by_contig[info["contig"]].append(rec)

    n_frags = sum(len(v) for v in by_contig.values())
    log.info("Loaded %d %s fragments across %d contigs",
             n_frags, tool, len(by_contig))
    return by_contig, meta


def check_contig_id_overlap(g_by_contig: dict, v_by_contig: dict):
    g_set  = set(g_by_contig.keys())
    v_set  = set(v_by_contig.keys())
    shared = g_set & v_set
    if not shared:
        log.warning(
            "SANITY CHECK FAILED: geNomad and VIBRANT share ZERO contig IDs. "
            "This will produce zero consensus groups. "
            "Check that both FASTAs use identical contig naming. "
            "geNomad example: '%s'  |  VIBRANT example: '%s'",
            next(iter(g_set), "N/A"), next(iter(v_set), "N/A"),
        )
    else:
        log.info("Contig ID overlap: %d shared  (geNomad=%d, VIBRANT=%d)",
                 len(shared), len(g_set), len(v_set))


# ─────────────────────────────────────────────────────────────────────────────
# Assembly accessor
# ─────────────────────────────────────────────────────────────────────────────

class AssemblyAccessor:
    """pyfaidx indexed when available; in-memory fallback with size warning."""

    def __init__(self, path: str):
        self.path  = path
        self._mode = None

        if Fasta is not None:
            try:
                self._fa   = Fasta(path)
                self._mode = "indexed"
                log.info("Assembly accessor: pyfaidx indexed mode (%s)", path)
                return
            except Exception as exc:
                log.warning("pyfaidx failed (%s) — loading into memory", exc)

        log.info("Loading combined assembly into memory ...")
        self._mem  = {}
        total_bp   = 0
        for rec in SeqIO.parse(path, "fasta"):
            seq = str(rec.seq)
            self._mem[rec.id] = seq
            total_bp += len(seq)
        self._mode = "memory"
        log.info("Loaded %d contigs (~%.0f Mb). Install pyfaidx to avoid this.",
                 len(self._mem), total_bp / 1e6)

    def get_contig_seq(self, contig_id: str):
        if self._mode == "indexed":
            return str(self._fa[contig_id][:].seq) \
                   if contig_id in self._fa else None
        return self._mem.get(contig_id)


# ─────────────────────────────────────────────────────────────────────────────
# Batch minimap2  (one subprocess call per contig)
# ─────────────────────────────────────────────────────────────────────────────

def _write_multi_fasta(records, path: str):
    with open(path, "w") as fh:
        for seq_id, seq in records:
            fh.write(f">{seq_id}\n{seq}\n")


def _write_single_fasta(seq_id: str, seq: str, path: str):
    with open(path, "w") as fh:
        fh.write(f">{seq_id}\n{seq}\n")


def _run_minimap2(query: str, target: str, paf: str, threads: int) -> bool:
    cmd = [
        "minimap2", "-c", "--cs", "-x", "asm5",
        "-t", str(threads), "--secondary=yes", "-N", "10",
        target, query,
    ]
    try:
        with open(paf, "w") as fh:
            subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, check=True)
        return True
    except subprocess.CalledProcessError as exc:
        log.error("minimap2 failed: %s", exc.stderr.decode())
        return False
    except FileNotFoundError:
        sys.exit("minimap2 not found in PATH.")


def _parse_paf(paf_path: str) -> dict:
    """Return { query_id: [hit_dict, ...] } without any filtering."""
    hits = defaultdict(list)
    if not os.path.exists(paf_path):
        return hits
    with open(paf_path) as fh:
        for line in fh:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue
            nm = None
            for tag in cols[12:]:
                if tag.startswith("NM:i:"):
                    nm = int(tag.split(":")[2])
                    break
            hits[cols[0]].append({
                "query":        cols[0],
                "query_len":    int(cols[1]),
                "query_start":  int(cols[2]),
                "query_end":    int(cols[3]),
                "strand":       cols[4],
                "target_start": int(cols[7]),
                "target_end":   int(cols[8]),
                "aln_len":      int(cols[10]),
                "nm":           nm,
            })
    return hits


def _same_placement(h1: dict, h2: dict, threshold: float = 0.90) -> bool:
    """
    Two hits are the 'same placement' only if they overlap by >= threshold
    of the SHORTER interval.  Raised from 0.50 (V4) to 0.90 (V5) to avoid
    calling clearly different loci the same placement.
    """
    s1, e1  = h1["target_start"], h1["target_end"]
    s2, e2  = h2["target_start"], h2["target_end"]
    overlap = max(0, min(e1, e2) - max(s1, s2))
    shorter = min(e1 - s1, e2 - s2)
    return shorter > 0 and (overlap / shorter) >= threshold


def _anchor_from_hits(
    frag_id:      str,
    frag_len:     int,
    hits:         list,
    min_cov:      float,
    min_identity: float,
) -> dict:
    """
    Apply identity + coverage filters, deduplicate placements, classify.
    Identity uses  1 - NM/aln_len  (tolerates minor boundary trimming).
    """
    passing = []
    for h in hits:
        if h["nm"] is None or h["aln_len"] == 0:
            continue
        identity = 1.0 - h["nm"] / h["aln_len"]
        if identity < min_identity:
            continue
        cov = (h["query_end"] - h["query_start"]) / frag_len
        if cov < min_cov:
            continue
        passing.append((h, cov, identity))

    if not passing:
        return {"status": "failed", "reason": "no_qualifying_hit",
                "n_hits_passing": 0, "candidate_coords": ""}

    # Deduplicate at 90% overlap threshold
    distinct = []
    for h, cov, ident in passing:
        if not any(_same_placement(h, ex[0]) for ex in distinct):
            distinct.append((h, cov, ident))

    cand_str = ";".join(
        f"{h['target_start']}-{h['target_end']}" for h, _, _ in distinct
    )

    if len(distinct) > 1:
        return {"status": "ambiguous", "n_hits_passing": len(distinct),
                "candidate_coords": cand_str}

    best, cov, ident = distinct[0]
    return {
        "status":           "anchored",
        "start":            best["target_start"],
        "end":              best["target_end"],
        "strand":           best["strand"],
        "coverage":         round(cov,   4),
        "identity":         round(ident, 5),
        "n_hits_passing":   1,
        "candidate_coords": cand_str,
    }


def batch_anchor_contig(
    contig_id:    str,
    contig_seq:   str,
    all_recs:     list,         # [(SeqRecord, tool_str), ...]
    tmpdir:       str,
    min_cov:      float,
    min_identity: float,
    threads:      int,
) -> dict:
    """
    ONE minimap2 call for all fragments on a contig.
    Returns { frag_id: anchor_result_dict }.
    """
    safe = re.sub(r"[^\w.-]", "_", contig_id)
    q    = os.path.join(tmpdir, f"{safe}__query.fna")
    t    = os.path.join(tmpdir, f"{safe}__target.fna")
    paf  = os.path.join(tmpdir, f"{safe}__self.paf")

    _write_multi_fasta([(r.id, str(r.seq)) for r, _ in all_recs], q)
    _write_single_fasta(contig_id, contig_seq, t)

    ok  = _run_minimap2(q, t, paf, threads)
    for p in (q, t):
        try: os.remove(p)
        except OSError: pass

    if not ok:
        return {r.id: {"status": "failed", "reason": "minimap2_error",
                       "n_hits_passing": 0, "candidate_coords": ""}
                for r, _ in all_recs}

    hits_by_query = _parse_paf(paf)
    try: os.remove(paf)
    except OSError: pass

    return {
        r.id: _anchor_from_hits(r.id, len(r.seq),
                                hits_by_query.get(r.id, []),
                                min_cov, min_identity)
        for r, _ in all_recs
    }


# ─────────────────────────────────────────────────────────────────────────────
# Interval utilities
# ─────────────────────────────────────────────────────────────────────────────

def _merge_intervals_gap(intervals: list, gap: int) -> list:
    """
    Merge a list of (start, end, meta_dict) tuples with gap tolerance.
    Returns list of merged dicts:
      { start, end, members: [meta_dict, ...] }
    Used for same-tool pre-merging before bipartite matching.
    """
    if not intervals:
        return []
    ordered = sorted(intervals, key=lambda x: x[0])
    groups  = [{"start": ordered[0][0], "end": ordered[0][1],
                "members": [ordered[0][2]]}]
    for s, e, meta in ordered[1:]:
        if s <= groups[-1]["end"] + gap:
            groups[-1]["end"] = max(groups[-1]["end"], e)
            groups[-1]["members"].append(meta)
        else:
            groups.append({"start": s, "end": e, "members": [meta]})
    return groups


def union_length(intervals) -> int:
    """True covered bp after merging a list of (start, end) pairs."""
    if not intervals:
        return 0
    merged = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append([s, e])
    return sum(e - s for s, e in merged)


def union_overlap_bp(ivals_a: list, ivals_b: list) -> int:
    """Overlap bp between the UNION of interval list A and list B."""
    if not ivals_a or not ivals_b:
        return 0

    def _merge(iv):
        out = []
        for s, e in sorted(iv):
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append([s, e])
        return out

    ma, mb = _merge(ivals_a), _merge(ivals_b)
    total = 0
    i = j = 0
    while i < len(ma) and j < len(mb):
        s1, e1 = ma[i]
        s2, e2 = mb[j]
        total += max(0, min(e1, e2) - max(s1, s2))
        if e1 < e2: i += 1
        else:        j += 1
    return total


# ─────────────────────────────────────────────────────────────────────────────
# Bipartite reciprocal-best matching
# ─────────────────────────────────────────────────────────────────────────────

def bipartite_match(
    g_groups:               list,   # list of merged geNomad groups
    v_groups:               list,   # list of merged VIBRANT groups
    min_overlap_bp:         int,
    min_reciprocal_overlap: float,
    strict_consensus:       bool,
) -> tuple:
    """
    Maximum-weight bipartite matching between geNomad and VIBRANT
    interval groups.

    Each geNomad group matches at most ONE VIBRANT group and vice versa.
    This prevents transitive chaining (A—B—C becoming one component).

    Edge weight = min(overlap_frac_g, overlap_frac_v) using TRUE union bp.
    An edge is eligible only when:
      overlap_bp >= min_overlap_bp  AND
      overlap_frac_g >= min_reciprocal_overlap  AND
      overlap_frac_v >= min_reciprocal_overlap
    (when strict_consensus=True; otherwise any positive overlap qualifies)

    Uses scipy.optimize.linear_sum_assignment (Hungarian algorithm) on the
    negated weight matrix to find the maximum-weight matching.

    Returns:
      matched_pairs : list of (g_group, v_group, overlap_metrics_dict)
      unmatched_g   : list of geNomad groups with no valid match
      unmatched_v   : list of VIBRANT groups with no valid match
    """
    ng = len(g_groups)
    nv = len(v_groups)

    if ng == 0 or nv == 0:
        return [], g_groups[:], v_groups[:]

    # Build weight matrix  (ng × nv)
    weight = np.zeros((ng, nv), dtype=float)
    edge_metrics = {}   # (gi, vi) -> metrics dict

    for gi, gg in enumerate(g_groups):
        g_ivals = [(m["start"], m["end"]) for m in gg["members"]]
        g_union = union_length(g_ivals)

        for vi, vg in enumerate(v_groups):
            v_ivals = [(m["start"], m["end"]) for m in vg["members"]]
            v_union = union_length(v_ivals)

            ovlp_bp = union_overlap_bp(g_ivals, v_ivals)
            if ovlp_bp == 0:
                continue

            frac_g = ovlp_bp / g_union if g_union > 0 else 0.0
            frac_v = ovlp_bp / v_union if v_union > 0 else 0.0

            if strict_consensus:
                if ovlp_bp < min_overlap_bp:
                    continue
                if frac_g < min_reciprocal_overlap:
                    continue
                if frac_v < min_reciprocal_overlap:
                    continue

            score = min(frac_g, frac_v)
            weight[gi, vi] = score
            edge_metrics[(gi, vi)] = {
                "overlap_bp":           ovlp_bp,
                "genomad_union_bp":     g_union,
                "vibrant_union_bp":     v_union,
                "overlap_frac_genomad": round(frac_g, 4),
                "overlap_frac_vibrant": round(frac_v, 4),
            }

    # Hungarian algorithm on negated weights (minimises → maximises)
    row_ind, col_ind = linear_sum_assignment(-weight)

    matched_pairs = []
    matched_g     = set()
    matched_v     = set()

    for gi, vi in zip(row_ind, col_ind):
        if weight[gi, vi] > 0:           # zero weight = no valid edge
            metrics = edge_metrics[(gi, vi)].copy()
            matched_pairs.append((g_groups[gi], v_groups[vi], metrics))
            matched_g.add(gi)
            matched_v.add(vi)

    unmatched_g = [g_groups[i] for i in range(ng) if i not in matched_g]
    unmatched_v = [v_groups[i] for i in range(nv) if i not in matched_v]

    return matched_pairs, unmatched_g, unmatched_v


# ─────────────────────────────────────────────────────────────────────────────
# Overlap metrics and support score
# ─────────────────────────────────────────────────────────────────────────────

def compute_boundary_metrics(g_group: dict, v_group: dict) -> dict:
    """
    Compute inner (intersection) and outer (union envelope) boundaries,
    and the boundary disagreement, for one matched consensus pair.
    """
    g_starts = [m["start"] for m in g_group["members"]]
    g_ends   = [m["end"]   for m in g_group["members"]]
    v_starts = [m["start"] for m in v_group["members"]]
    v_ends   = [m["end"]   for m in v_group["members"]]

    outer_start = min(min(g_starts), min(v_starts))
    outer_end   = max(max(g_ends),   max(v_ends))
    inner_start = max(min(g_starts), min(v_starts))   # intersection start
    inner_end   = min(max(g_ends),   max(v_ends))     # intersection end

    ext_left  = max(0, inner_start - outer_start)
    ext_right = max(0, outer_end   - inner_end)

    return {
        "outer_start":              outer_start,
        "outer_end":                outer_end,
        "consensus_inner_start":    inner_start,
        "consensus_inner_end":      inner_end,
        "boundary_extension_left":  ext_left,
        "boundary_extension_right": ext_right,
        "boundary_disagreement_bp": ext_left + ext_right,
    }


def compute_support_score(
    overlap_metrics: dict,
    boundary_metrics: dict,
    merged_len: int,
    label: str,
) -> float:
    """
    Composite QC ranking score (0–1).
    For consensus groups:
      0.4 × min(overlap_frac_g, overlap_frac_v)     — reciprocal support
      0.4 × max(0, 1 - boundary_disagreement / len)  — boundary agreement
      0.2 × 1.0                                      — dual-tool confirmed
    For single-tool groups: 0.5 (no cross-tool validation possible).
    Note: treat this as a QC ranking metric, not a probabilistic confidence.
    """
    if label == "consensus":
        recip = min(
            overlap_metrics.get("overlap_frac_genomad", 0.0),
            overlap_metrics.get("overlap_frac_vibrant",  0.0),
        )
        disagree = boundary_metrics.get("boundary_disagreement_bp", 0)
        bnd_sim  = max(0.0, 1.0 - disagree / merged_len) if merged_len > 0 else 0.0
        return round(0.4 * recip + 0.4 * bnd_sim + 0.2, 4)
    return 0.5


# ─────────────────────────────────────────────────────────────────────────────
# Per-contig processing
# ─────────────────────────────────────────────────────────────────────────────

_EMPTY_METRICS = {k: "" for k in [
    "overlap_bp", "genomad_union_bp", "vibrant_union_bp",
    "overlap_frac_genomad", "overlap_frac_vibrant",
    "boundary_extension_left", "boundary_extension_right",
    "boundary_disagreement_bp", "consensus_inner_start", "consensus_inner_end",
]}


def process_contig(
    contig_id:              str,
    g_recs:                 list,
    v_recs:                 list,
    g_meta:                 dict,
    v_meta:                 dict,
    assembly:               AssemblyAccessor,
    tmpdir:                 str,
    min_self_cov:           float,
    min_identity:           float,
    merge_gap:              int,
    threads:                int,
    strict_consensus:       bool,
    min_overlap_bp:         int,
    min_reciprocal_overlap: float,
) -> dict:

    contig_seq = assembly.get_contig_seq(contig_id)
    if contig_seq is None:
        log.warning("Contig %s not in assembly — skipping "
                    "(%d geNomad, %d VIBRANT lost)",
                    contig_id, len(g_recs), len(v_recs))
        return {
            "consensus": [], "genomad_only": [], "vibrant_only": [],
            "ambiguous": [], "anchor_rows": [], "merge_rows": [],
            "counts": {
                "contig": contig_id, "status": "missing_from_assembly",
                "genomad_total": len(g_recs), "vibrant_total": len(v_recs),
            },
        }

    # ── 1. Batch anchor all fragments (one minimap2 call) ──────────────────
    all_recs = [(r, "genomad") for r in g_recs] + [(r, "vibrant") for r in v_recs]
    anchor_results = batch_anchor_contig(
        contig_id    = contig_id,
        contig_seq   = contig_seq,
        all_recs     = all_recs,
        tmpdir       = tmpdir,
        min_cov      = min_self_cov,
        min_identity = min_identity,
        threads      = threads,
    )

    anchor_rows    = []
    g_anchored     = []    # (start, end, meta_dict)
    v_anchored     = []
    ambiguous_recs = []

    for rec, tool in all_recs:
        result   = anchor_results.get(rec.id, {
            "status": "failed", "reason": "missing_from_batch",
            "n_hits_passing": 0, "candidate_coords": "",
        })
        rec_meta = (g_meta if tool == "genomad" else v_meta).get(rec.id, {})

        anchor_rows.append({
            "contig":           contig_id,
            "tool":             tool,
            "frag_id":          rec.id,
            "frag_len":         len(rec.seq),
            "status":           result["status"],
            "start":            result.get("start",           ""),
            "end":              result.get("end",             ""),
            "strand":           result.get("strand",          ""),
            "coverage":         result.get("coverage",        ""),
            "identity":         result.get("identity",        ""),
            "n_hits_passing":   result.get("n_hits_passing",  ""),
            "candidate_coords": result.get("candidate_coords",""),
            "reason":           result.get("reason",          ""),
            "hint_start":       rec_meta.get("hint_start",   ""),
            "hint_end":         rec_meta.get("hint_end",     ""),
        })

        if result["status"] == "anchored":
            iv = (result["start"], result["end"],
                  {"frag_id": rec.id, "start": result["start"],
                   "end": result["end"], "tool": tool})
            if tool == "genomad":
                g_anchored.append(iv)
            else:
                v_anchored.append(iv)
        elif result["status"] == "ambiguous":
            ambiguous_recs.append(rec)

    # ── 2. Pre-merge same-tool intervals (applies merge_gap correctly) ─────
    g_groups = _merge_intervals_gap(g_anchored, merge_gap)
    v_groups = _merge_intervals_gap(v_anchored, merge_gap)

    # ── 3. Bipartite matching ──────────────────────────────────────────────
    matched_pairs, unmatched_g, unmatched_v = bipartite_match(
        g_groups               = g_groups,
        v_groups               = v_groups,
        min_overlap_bp         = min_overlap_bp,
        min_reciprocal_overlap = min_reciprocal_overlap,
        strict_consensus       = strict_consensus,
    )

    consensus_recs    = []
    genomad_only_recs = []
    vibrant_only_recs = []
    merge_rows        = []
    label_counters    = defaultdict(int)

    def _emit(label, start, end, members, ovlp_metrics, bnd_metrics):
        label_counters[label] += 1
        group_id   = f"{label}_{label_counters[label]}"
        merged_len = end - start
        frag_ids   = [m["frag_id"] for m in members]
        tools_used = sorted({m["tool"] for m in members})
        score      = compute_support_score(ovlp_metrics, bnd_metrics,
                                           merged_len, label)
        seq        = contig_seq[start:end]

        # Use member_ids as FASTA header
        member_ids_header = "|".join(frag_ids)

        rec = SeqRecord(
            Seq(seq), id=member_ids_header,
            description=(
                f"group_id={group_id} contig={contig_id} start={start} end={end} "
                f"tools={'|'.join(tools_used)} support_score={score}"
            ),
        )

        row = {
            "group_id":                 group_id,
            "contig":                   contig_id,
            "start":                    start,
            "end":                      end,
            "length":                   merged_len,
            "label":                    label,
            "tools":                    "|".join(tools_used),
            "n_members":                len(members),
            "member_ids":               "|".join(frag_ids)[:500],
            "support_score":            score,
            # overlap metrics
            "overlap_bp":               ovlp_metrics.get("overlap_bp",           ""),
            "genomad_union_bp":         ovlp_metrics.get("genomad_union_bp",     ""),
            "vibrant_union_bp":         ovlp_metrics.get("vibrant_union_bp",     ""),
            "overlap_frac_genomad":     ovlp_metrics.get("overlap_frac_genomad", ""),
            "overlap_frac_vibrant":     ovlp_metrics.get("overlap_frac_vibrant", ""),
            # boundary metrics
            "boundary_extension_left":  bnd_metrics.get("boundary_extension_left",  ""),
            "boundary_extension_right": bnd_metrics.get("boundary_extension_right", ""),
            "boundary_disagreement_bp": bnd_metrics.get("boundary_disagreement_bp", ""),
            # dual output: outer (union envelope) and inner (intersection)
            "consensus_inner_start":    bnd_metrics.get("consensus_inner_start", ""),
            "consensus_inner_end":      bnd_metrics.get("consensus_inner_end",   ""),
        }
        return rec, row

    # Consensus pairs
    for gg, vg, ovlp_m in matched_pairs:
        bnd_m    = compute_boundary_metrics(gg, vg)
        start    = bnd_m["outer_start"]
        end      = bnd_m["outer_end"]
        members  = gg["members"] + vg["members"]
        rec, row = _emit("consensus", start, end, members, ovlp_m, bnd_m)
        consensus_recs.append(rec)
        merge_rows.append(row)

    # Unmatched geNomad
    for gg in unmatched_g:
        start    = gg["start"]
        end      = gg["end"]
        rec, row = _emit("genomad_only", start, end,
                         gg["members"], dict(_EMPTY_METRICS), dict(_EMPTY_METRICS))
        genomad_only_recs.append(rec)
        merge_rows.append(row)

    # Unmatched VIBRANT
    for vg in unmatched_v:
        start    = vg["start"]
        end      = vg["end"]
        rec, row = _emit("vibrant_only", start, end,
                         vg["members"], dict(_EMPTY_METRICS), dict(_EMPTY_METRICS))
        vibrant_only_recs.append(rec)
        merge_rows.append(row)

    counts = {
        "contig":              contig_id,
        "status":              "ok",
        "genomad_total":       len(g_recs),
        "vibrant_total":       len(v_recs),
        "anchored_ok":         sum(1 for r in anchor_rows if r["status"] == "anchored"),
        "ambiguous":           sum(1 for r in anchor_rows if r["status"] == "ambiguous"),
        "anchor_failed":       sum(1 for r in anchor_rows if r["status"] == "failed"),
        "consensus_groups":    len(consensus_recs),
        "genomad_only_groups": len(genomad_only_recs),
        "vibrant_only_groups": len(vibrant_only_recs),
    }

    return {
        "consensus":    consensus_recs,
        "genomad_only": genomad_only_recs,
        "vibrant_only": vibrant_only_recs,
        "ambiguous":    ambiguous_recs,
        "anchor_rows":  anchor_rows,
        "merge_rows":   merge_rows,
        "counts":       counts,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Multiprocessing
# ─────────────────────────────────────────────────────────────────────────────

_ASSEMBLY = None


def _init_worker(assembly_path: str):
    global _ASSEMBLY
    _ASSEMBLY = AssemblyAccessor(assembly_path)


def _worker(args_tuple):
    (contig_id, g_recs, v_recs, g_meta, v_meta,
     tmpdir, min_self_cov, min_identity, merge_gap, threads,
     strict_consensus, min_overlap_bp, min_reciprocal_overlap) = args_tuple

    return process_contig(
        contig_id              = contig_id,
        g_recs                 = g_recs,
        v_recs                 = v_recs,
        g_meta                 = g_meta,
        v_meta                 = v_meta,
        assembly               = _ASSEMBLY,
        tmpdir                 = tmpdir,
        min_self_cov           = min_self_cov,
        min_identity           = min_identity,
        merge_gap              = merge_gap,
        threads                = threads,
        strict_consensus       = strict_consensus,
        min_overlap_bp         = min_overlap_bp,
        min_reciprocal_overlap = min_reciprocal_overlap,
    ), contig_id


# ─────────────────────────────────────────────────────────────────────────────
# Output
# ─────────────────────────────────────────────────────────────────────────────

def write_fasta_append(records: list, path: str):
    if not records:
        return
    with open(path, "a") as fh:
        SeqIO.write(records, fh, "fasta")


ANCHOR_FIELDS = [
    "contig", "tool", "frag_id", "frag_len", "status",
    "start", "end", "strand", "coverage", "identity",
    "n_hits_passing", "candidate_coords", "reason",
    "hint_start", "hint_end",
]

MERGE_FIELDS = [
    "group_id", "contig", "start", "end", "length",
    "label", "tools", "n_members", "member_ids", "support_score",
    # overlap
    "overlap_bp", "genomad_union_bp", "vibrant_union_bp",
    "overlap_frac_genomad", "overlap_frac_vibrant",
    # boundary
    "boundary_extension_left", "boundary_extension_right",
    "boundary_disagreement_bp",
    # dual coordinates
    "consensus_inner_start", "consensus_inner_end",
]

COUNT_FIELDS = [
    "contig", "status", "genomad_total", "vibrant_total",
    "anchored_ok", "ambiguous", "anchor_failed",
    "consensus_groups", "genomad_only_groups", "vibrant_only_groups",
]


def write_tsv(rows: list, fields: list, path: str):
    with open(path, "w") as fh:
        fh.write("\t".join(fields) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(f, "")) for f in fields) + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Prophage consensus pipeline:\n"
            "batch minimap2 + bipartite reciprocal-best matching + "
            "true union coverage + dual boundary output."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--genomad_fasta",   required=True)
    parser.add_argument("--vibrant_fasta",   required=True)
    parser.add_argument("--assembly_fasta",  required=True)
    parser.add_argument("--outdir",          required=True)

    parser.add_argument("--min_self_coverage", type=float, default=0.99,
                        help="Min fraction of fragment length aligned (default 0.99)")
    parser.add_argument("--min_identity", type=float, default=0.99,
                        help="Min alignment identity 1-NM/aln_len (default 0.99). "
                             "Tolerates minor boundary trimming by prediction tools.")
    parser.add_argument("--min_overlap_bp", type=int, default=1000,
                        help="Min raw overlap bp for a consensus match (default 1000)")
    parser.add_argument("--min_reciprocal_overlap", type=float, default=0.50,
                        help="Min overlap fraction required from BOTH tools (default 0.50)")
    parser.add_argument("--no_strict_consensus", action="store_true",
                        help="Disable strict reciprocal overlap. Not recommended.")
    parser.add_argument("--merge_gap", type=int, default=0,
                        help="Gap (bp) for pre-merging same-tool intervals on one contig "
                             "before bipartite matching. Default 0 (strict overlap only). "
                             "Use cautiously — large values can merge distinct insertions.")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--minimap_threads", type=int, default=2,
                        help="minimap2 threads per contig (default 2; "
                             "batch mode makes extra threads worthwhile)")
    parser.add_argument("--keep_tmp", action="store_true")

    args             = parser.parse_args()
    strict_consensus = not args.no_strict_consensus

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("prophage_seq_consensusV5  —  starting")
    log.info("min_self_coverage        : %.3f", args.min_self_coverage)
    log.info("min_identity             : %.3f", args.min_identity)
    log.info("strict_consensus         : %s",   strict_consensus)
    if strict_consensus:
        log.info("  min_overlap_bp         : %d",   args.min_overlap_bp)
        log.info("  min_reciprocal_overlap : %.2f", args.min_reciprocal_overlap)
    log.info("merge_gap (same-tool)    : %d bp", args.merge_gap)
    log.info("=" * 60)

    log.info("Loading geNomad fragments ...")
    g_by_contig, g_meta = load_fasta_by_contig(args.genomad_fasta, "genomad")
    log.info("Loading VIBRANT fragments ...")
    v_by_contig, v_meta = load_fasta_by_contig(args.vibrant_fasta, "vibrant")

    check_contig_id_overlap(g_by_contig, v_by_contig)

    all_contigs = sorted(set(g_by_contig) | set(v_by_contig))
    log.info("Total contigs with >= 1 prediction: %d", len(all_contigs))

    if args.keep_tmp:
        tmpdir = str(outdir / "tmp_anchor")
        Path(tmpdir).mkdir(exist_ok=True)
    else:
        _td    = tempfile.TemporaryDirectory(dir=str(outdir), prefix="tmp_anchor_")
        tmpdir = _td.name

    fa_consensus = outdir / "consensus.fasta"
    fa_g_only    = outdir / "genomad_only.fasta"
    fa_v_only    = outdir / "vibrant_only.fasta"
    fa_ambig     = outdir / "ambiguous.fasta"
    tsv_anchor   = outdir / "anchor_report.tsv"
    tsv_merge    = outdir / "merge_report.tsv"
    tsv_counts   = outdir / "contig_summary.tsv"

    for f in [fa_consensus, fa_g_only, fa_v_only, fa_ambig]:
        f.unlink(missing_ok=True)

    worker_args = [
        (cid,
         g_by_contig.get(cid, []), v_by_contig.get(cid, []),
         g_meta, v_meta, tmpdir,
         args.min_self_coverage, args.min_identity, args.merge_gap,
         args.minimap_threads, strict_consensus,
         args.min_overlap_bp, args.min_reciprocal_overlap)
        for cid in all_contigs
    ]

    all_anchor = []
    all_merge  = []
    all_counts = []

    log.info("Processing %d contigs (workers=%d) ...",
             len(all_contigs), args.workers)
    with multiprocessing.Pool(
        processes   = args.workers,
        initializer = _init_worker,
        initargs    = (args.assembly_fasta,),
    ) as pool:
        for result, _ in pool.imap_unordered(_worker, worker_args):
            write_fasta_append(result["consensus"],    str(fa_consensus))
            write_fasta_append(result["genomad_only"], str(fa_g_only))
            write_fasta_append(result["vibrant_only"], str(fa_v_only))
            write_fasta_append(result["ambiguous"],    str(fa_ambig))
            all_anchor.extend(result["anchor_rows"])
            all_merge.extend(result["merge_rows"])
            all_counts.append(result["counts"])

    if not args.keep_tmp:
        _td.cleanup()

    write_tsv(all_anchor, ANCHOR_FIELDS, str(tsv_anchor))
    write_tsv(all_merge,  MERGE_FIELDS,  str(tsv_merge))
    write_tsv(all_counts, COUNT_FIELDS,  str(tsv_counts))

    tc  = sum(c.get("consensus_groups",    0) for c in all_counts)
    tg  = sum(c.get("genomad_only_groups", 0) for c in all_counts)
    tv  = sum(c.get("vibrant_only_groups", 0) for c in all_counts)
    ta  = sum(c.get("ambiguous",           0) for c in all_counts)
    tf  = sum(c.get("anchor_failed",       0) for c in all_counts)
    tm  = sum(1 for c in all_counts if c.get("status") == "missing_from_assembly")

    log.info("-" * 60)
    log.info("DONE")
    log.info("Contigs processed:       %d", len(all_contigs))
    if tm:
        log.warning("Contigs missing from assembly: %d (skipped)", tm)
    log.info("Consensus groups:        %d  ->  %s", tc, fa_consensus)
    log.info("geNomad-only groups:     %d  ->  %s", tg, fa_g_only)
    log.info("VIBRANT-only groups:     %d  ->  %s", tv, fa_v_only)
    log.info("Ambiguous fragments:     %d  ->  %s", ta, fa_ambig)
    log.info("Self-anchor failures:    %d  (see anchor_report.tsv)", tf)
    log.info("Anchor report:           %s", tsv_anchor)
    log.info("Merge report:            %s", tsv_merge)
    log.info("Contig summary:          %s", tsv_counts)
    log.info("-" * 60)


if __name__ == "__main__":
    main()

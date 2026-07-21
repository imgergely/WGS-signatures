#!/usr/bin/env python2
""" IsoMut configuration for the 42-clone joint call in this protocol.

Adapted from the distributed IsoMut example script:
https://github.com/riblidezso/isomut/blob/master/isomut_example_script.py

Run with Python 2.7 from the package directory after editing the path constants
and, when needed, EXCLUDED_SAMPLE_IDS below. Table_S1_clone_key.csv defines
the exact analysed panel and records which libraries are omitted.
"""

import csv
import glob
import os
import sys

# ---------------------------- settings to edit ----------------------------
ISOMUT_DIR = "."
REF_FASTA = "reference/hg38.analysisSet.fa"
INPUT_DIR = "bams"
OUTPUT_DIR = "isomut_output"
CLONE_KEY = "Table_S1_clone_key.csv"

# Original clone IDs intentionally excluded from the joint call after QC.
# Every ID listed here must be marked "(omitted)" in Table S1. Excluded BAMs
# may remain in INPUT_DIR for QC; only the 42 non-omitted BAMs are passed to
# IsoMut through params["bam_filenames"].
EXCLUDED_SAMPLE_IDS = set(["BAP2KO3"])
# -------------------------------------------------------------------------

# Import the IsoMut wrapper and expose the IsoMut executable on PATH.
sys.path.insert(0, os.path.join(ISOMUT_DIR, "src"))
os.environ["PATH"] += os.pathsep + os.path.join(ISOMUT_DIR, "src")
from isomut_wrappers import run_isomut


def fail(message):
    raise RuntimeError(message)


def load_panel(clone_key_path):
    """Return analysed IDs and omitted IDs from Table S1."""
    if not os.path.isfile(clone_key_path):
        fail("Clone key not found: %s" % clone_key_path)

    handle = open(clone_key_path, "rb")
    try:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            fail("Table S1 is empty: %s" % clone_key_path)
        required = set(["Protocol_name", "Original_clone_ID"])
        if not required.issubset(set(reader.fieldnames)):
            fail("Table S1 must contain Protocol_name and Original_clone_ID columns.")

        analysed_ids = []
        omitted_ids = []
        for row in reader:
            protocol_name = row["Protocol_name"].strip()
            original_id = row["Original_clone_ID"].strip()
            if not original_id:
                fail("Table S1 contains an empty Original_clone_ID value.")
            if protocol_name.lower() in set(["omitted", "(omitted)"]):
                omitted_ids.append(original_id)
            else:
                analysed_ids.append(original_id)
    finally:
        handle.close()

    if len(analysed_ids) != 42 or len(set(analysed_ids)) != 42:
        fail(
            "Expected 42 unique analysed clone IDs in Table S1, found %d."
            % len(set(analysed_ids))
        )
    if set(omitted_ids) != EXCLUDED_SAMPLE_IDS:
        fail(
            "EXCLUDED_SAMPLE_IDS must exactly match the clone IDs marked omitted "
            "in Table S1. Config: %s; Table S1: %s"
            % (sorted(EXCLUDED_SAMPLE_IDS), sorted(set(omitted_ids)))
        )
    return sorted(analysed_ids), sorted(omitted_ids)


def has_bam_index(bam_path):
    """Accept either sample.bam.bai or sample.bai index naming."""
    root, extension = os.path.splitext(bam_path)
    return os.path.isfile(bam_path + ".bai") or os.path.isfile(root + ".bai")


if not os.path.isfile(REF_FASTA):
    fail("Reference FASTA not found: %s" % REF_FASTA)
if not os.path.isfile(REF_FASTA + ".fai"):
    fail(
        "Reference FASTA index not found: %s.fai. Run `samtools faidx %s`."
        % (REF_FASTA, REF_FASTA)
    )
if not os.path.isdir(INPUT_DIR):
    fail("BAM directory not found: %s" % INPUT_DIR)
if not os.path.isfile(os.path.join(ISOMUT_DIR, "src", "isomut_wrappers.py")):
    fail("IsoMut wrapper not found under %s/src." % ISOMUT_DIR)
if not os.path.isfile(os.path.join(ISOMUT_DIR, "src", "isomut")):
    fail("Compiled IsoMut executable not found under %s/src." % ISOMUT_DIR)

analysed_ids, omitted_ids = load_panel(CLONE_KEY)
expected_bams = [sample_id + ".final.bam" for sample_id in analysed_ids]
excluded_bams = set([sample_id + ".final.bam" for sample_id in omitted_ids])
observed_bams = sorted(
    [os.path.basename(path) for path in glob.glob(os.path.join(INPUT_DIR, "*.final.bam"))]
)

missing = sorted(set(expected_bams) - set(observed_bams))
unexpected = sorted(
    set(observed_bams) - set(expected_bams) - excluded_bams
)
if missing or unexpected:
    fail(
        "BAM panel mismatch. Missing analysed BAMs: %s; unexpected BAMs: %s"
        % (missing or "none", unexpected or "none")
    )

excluded_present = sorted(set(observed_bams) & excluded_bams)
if excluded_present:
    print "Excluded from the IsoMut joint call: %s" % ", ".join(excluded_present)

missing_indexes = []
for bam_name in expected_bams:
    bam_path = os.path.join(INPUT_DIR, bam_name)
    if not has_bam_index(bam_path):
        missing_indexes.append(bam_name)
if missing_indexes:
    fail(
        "Missing BAM indexes for: %s. Run `samtools index <bam>` before IsoMut."
        % ", ".join(missing_indexes)
    )

if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# mutation-calling parameters used in this study.
params = {}
params["n_min_block"] = 200
params["n_conc_blocks"] = 8
params["ref_fasta"] = REF_FASTA
params["input_dir"] = INPUT_DIR + os.sep
params["output_dir"] = OUTPUT_DIR + os.sep
params["bam_filenames"] = expected_bams
params["min_sample_freq"] = 0.21
params["min_other_ref_freq"] = 0.93
params["cov_limit"] = 5
params["base_quality_limit"] = 30
params["min_gap_dist_snv"] = 0
params["min_gap_dist_indel"] = 20

print "Running IsoMut on %d analysed BAMs." % len(expected_bams)
run_isomut(params)

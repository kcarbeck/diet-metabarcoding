#!/usr/bin/env python3
import sys, os
import csv

def error(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

manifest, metadata = sys.argv[1:3]

# Validate manifest
with open(manifest) as f:
    reader = csv.reader(f, delimiter='\t')
    rows = list(reader)
    if rows[0] != ['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath']:
        error("Manifest header must be: sample-id, forward-absolute-filepath, reverse-absolute-filepath")
    sample_ids = set()
    for row in rows[1:]:
        if len(row) != 3:
            error(f"Manifest row has wrong number of columns: {row}")
        sid, fwd, rev = row
        if sid in sample_ids:
            error(f"Duplicate sample-id in manifest: {sid}")
        sample_ids.add(sid)
        if not os.path.exists(fwd):
            error(f"Forward file does not exist: {fwd}")
        if not os.path.exists(rev):
            error(f"Reverse file does not exist: {rev}")

# Validate metadata
with open(metadata) as f:
    reader = csv.reader(f, delimiter='\t')
    rows = list(reader)
    if not rows[0][0].startswith('#SampleID'):
        error("Metadata first column must be #SampleID")
    if not rows[1][0].startswith('#q2:types'):
        error("Metadata second row must start with #q2:types")
    meta_ids = set(row[0] for row in rows[2:])
    if sample_ids != meta_ids:
        error(f"Sample IDs in manifest and metadata do not match.\nManifest: {sample_ids}\nMetadata: {meta_ids}")

print("Validation successful.") 
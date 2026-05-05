# RAD known issues

## `scan-wl -w <kit-key>` — bitlist whitelist match returns 0 perfect matches

**Status:** open as of 2026-04-28. Discovered while benchmarking against the seqsim
empirical-rageseq simulation.

### Symptom

When running `rad scan-wl` with the `-w` flag pointing to a bundled kit key
(e.g. `-w 10x_3v1`, `-w 10x_3v2`, or `-w 10x_3v3`), the run loads the
correct number of whitelist barcodes (737,280 for v1/v2 or 6,794,880 for
v3), reports a non-zero number of extracted barcodes from the FASTQ, but
then reports `Total perfect matches: 0` and `Match rate: 0.00%`. As a
result, no barcodes pass the floor and `Final barcodes (TRUE): 0`.

### Reproducer

Using the empirical-rageseq sim from seqsim's manuscript artifacts:

```bash
# Fails — bitlist path
rad scan-wl \
  -i sim_rageseq_emp.fq.gz \
  -p CTACACGACGCTCTTCCGATCT \
  -n 16 -e 0.3 \
  -w 10x_3v2 \
  -t 4 \
  -o sw_bitlist
# → "Total perfect matches: 0"
# → "Final barcodes (TRUE): 0"

# Works — passing the raw text whitelist as a path
rad scan-wl \
  -i sim_rageseq_emp.fq.gz \
  -p CTACACGACGCTCTTCCGATCT \
  -n 16 -e 0.3 \
  -w /path/to/737K-august-2016.txt \
  -t 4 \
  -o sw_text
# → "Match rate: 33.78%"
# → "Final barcodes (TRUE): 2270" (or 3458 with --max-error tuned)
```

### Verification that the data is fine

The same FASTQ on the same RAD binary, with `-w 10x_3v2` vs `-w <text-path>`,
produces wildly different match rates with the same extracted-barcode set
upstream. Specifically: with `-w 10x_3v2`, the top extracted barcode
`CGTAGGCTCGCGTTTC` does **not** match (count=0), even though the same string
is provably present in the 737K-august-2016 text whitelist that 10x_3v2
references in `rad_config`'s key→file map. So extraction is right and the
text WL contents are right; the bug is in the bitlist→string comparison
path used by `scan-wl` when `-w` is a kit key.

### Affected component

`rad/include/rad/whitelist_generator.hpp` (or the bitlist loader it calls).
Likely the bitlist-decode-and-compare logic produces a different canonical
form for the extracted strings than for the whitelist entries.

### Workaround for now

Pass the raw whitelist text file path directly via `-w /path/to/whitelist.txt`
rather than a kit key. Resolve the kit key→file path manually by inspecting
`rad list` output and using the `Files:` column.

### Cross-reference

This was diagnosed in the course of seqsim manuscript benchmarking (see
`~/Desktop/manuscripts/seqsim/data/demux_bench_v2/rad_scanwl/`). seqsim's
CHANGELOG 0.2.0 references this as a known external dependency issue.

# Discovering read structure and reference whitelists

`rad discover --wl` identifies a nonspatial reference whitelist from a known
flanking primer and barcode length. It extracts raw barcodes using `scan-wl`,
counts exact matching reads independently against each available local catalog,
and returns the catalog with the largest count:

```sh
rad discover -i reads.fastq.gz --wl \
  -p CTACACGACGCTCTTCCGATCT -n 16 -o catalog -t 4
```

This mode requires `-p` and either `-n` or a masked pattern containing uppercase
`N` positions. It bypasses layout discovery. The command does not modify the
layout or whitelist registry.

## Whitelist identification

Raw extraction happens before reference filtering and cell calling. An extracted
barcode is compared independently with every eligible local catalog of its
length. If two catalogs contain that barcode, both receive its read count; a
nonmatch to one catalog remains available to the others.

A unique positive maximum is reported as `matched`. Equal maximum read counts
are reported as `tie`, retaining all tied catalogs. If every count is zero, the
result is `no_hits`. Distinct barcode counts and catalog sizes are descriptive;
they do not affect the winner or break a tie. Repeated observations of one
barcode each add to its catalog's read count, and dense catalogs receive no
penalty. Selection is performed separately for each extraction recipe.

Masked patterns use the same convention as `scan-wl`: uppercase `N` captures a
barcode base; lowercase `n` is a noncapturing wildcard. The capture length is
inferred from uppercase `N` positions, including noncontiguous captures:

```sh
rad discover -i reads.fastq.gz --wl \
  -p NNNNNNNNTCTTCAGCGTTCCCGAGANNNNNNNN -o catalog
```

Catalogs come from the RAD whitelist registry, including local user overrides.
Missing resources cannot contribute evidence. Aliases resolving to the same
catalog file share one row and one candidate. Different catalog files remain
separate candidates even when their contents are identical. Thus a winner can
identify a shared barcode family without distinguishing the chemistry aliases
that use it. It also does not identify which observed barcodes are cells; use
`scan-wl` for subsequent whitelist calling.

Spatial whitelist identification is unsupported. Known spatial registry aliases,
including `10x_Vis*` and `visium_hd*`, and their resolved catalog paths are excluded
from the comparison. Custom catalog names do not reliably identify a platform:
an otherwise eligible custom DNA catalog is compared as supplied, without
inferring whether it represents a spatial assay. This command does not recover
spatial coordinates or solve joint spatial-barcode boundaries.

## Layout discovery

Layout discovery remains experimental; the known-primer whitelist-only mode
bypasses it.

`--layout` samples an input file, finds known static sequences, and uses their
order and spacing to propose a read layout:

```sh
rad discover -i reads.fastq.gz --layout -o structure -t 4
```

For compatibility, omitting both mode flags requests layout discovery and
whitelist comparison using inferred extraction recipes. Supplying both
`--layout --wl` does the same. An explicit `-p` is accepted only with `--wl`
alone.

Discovery uses the static sequences in registered RAD layouts as its vocabulary.
Additional known primers, adapters or linkers can be supplied in a CSV:

```csv
id,seq
primer_a,CTACACGACGCTCTTCCGATCT
linker_b,CCACAGTCTCAAGCACGTGGAT
```

```sh
rad discover -i reads.fastq.gz --layout --parts parts.csv -o structure
```

The supplied rows identify parts to search for; their file order does not
prescribe their order in the reads. Discovery examines both orientations,
combines evidence across sampled reads, and emits the supported order with
intervening variable regions. Known flanking-part contracts supply initial
barcode roles and lengths; a stable bounded gap can adjust the remaining UMI
length. When whitelist comparison is requested, raw catalog matches supply
additional barcode evidence. A usable result is a native RAD layout CSV,
suitable for inspection and subsequent `rad prep` or `rad demux`:

```sh
rad prep --read-layout -l structure_layout.csv -o prepared
```

This is discovery from a vocabulary of known parts. It does not assemble novel
adapter sequences from arbitrary reads, and does not use the concatemer HMM to
learn a new topology. A missing primer or an unresolved variable region can
prevent a complete layout even when some parts are well supported. A random
sequence's diversity alone does not establish that it is a UMI or locate its
boundary against an unknown insert.

At least two supported literal anchors are required for a resolved scaffold;
poly-A and poly-T runs alone do not establish one. Repeated occurrences of the
same anchor in a substantial fraction of reads also leave structure discovery
partial. This version retains that evidence instead of treating repeated
capture units as a single molecule. Whitelist comparison on a capture recipe
can still produce a unique maximum while the layout remains partial. The
existing concatemer HMM is not invoked to split units during discovery.

## Sampling and outputs

Defaults are the first **10,000 records**, subject to a **50,000,000-base** cap.
These are prefix samples. Change the limits explicitly when the input is ordered
or the prefix is unrepresentative:

```sh
rad discover -i reads.fastq.gz --wl \
  -p CTACACGACGCTCTTCCGATCT -n 16 -o catalog \
  --sample-reads 25000 --sample-bases 100000000
```

The base cap is applied to complete records; a record that would exceed the cap
is not included. Increasing the sample can add evidence but cannot resolve
structural distinctions absent from the sequenced reads.

| File | Content |
| --- | --- |
| `<prefix>_whitelist_match.csv` | Per-recipe `matched`, `tie` or `no_hits` result, with the winning or tied catalog rows |
| `<prefix>_whitelists.csv` | Raw per-recipe catalog counts, sizes, aliases and extraction denominators, including losing and zero-hit candidates |
| `<prefix>_layout.csv` | Complete proposed native layout, when requested and resolved |
| `<prefix>_layout.partial.csv` | Partial structure when a requested complete layout cannot be inferred |
| `<prefix>_parts.csv` | Static-part support, including orientation counts, when structure inference runs |
| `<prefix>_discover.log` | Sampling totals, inference notes and outcome |

The match table columns are `status`, `recipe_id`, `catalog`, `aliases`,
`matched_reads`, `extracted_reads` and `matched_distinct`. A tie has one row per
tied catalog. A `no_hits` row has no catalog. The raw count table retains every
eligible catalog independently.

Whitelist-only mode with an explicit pattern writes the two whitelist tables
and log. Layout-only mode writes layout, parts and log. Existing output files
are not overwritten. Exit status is `0` when every requested whitelist recipe
has a unique positive maximum and any requested layout is complete; `2` means
a tie, no hits or incomplete layout; `1` means invalid arguments or an execution
error. A partial layout records evidence for review; it is not a resolved
demultiplexing recipe.

Mixed structures, unsequenced boundaries and low support can leave the layout
partial. Catalog selection reports the greatest exact-match read count; it does
not establish a unique chemistry or estimate mixture proportions.

# Concatemer boundaries and selective whitelist loading

The concatemer HMM identifies ordered observable elements. Its first and last
observable elements are not necessarily physical molecule ends. For example,
reverse Curio reads have an insert before poly(A), and a barcode after the
reverse-complement linker. A midpoint between the linker of one copy and the
poly(A) of the next copy can therefore lie inside the next insert.

Model construction precomputes physical head and tail ranges for every observable
state from the complete layout, including variable inserts and barcode/spacer
lengths. Per-read boundary placement uses these cached ranges: both facing edges
bounded, one bounded, or neither bounded. It does not walk the layout at each
junction. Decoder opening/closing flags still describe observable state topology.
They do not, by themselves, authorize clipping. This simplification preserves
the same cuts, uncertainty windows and unresolved-boundary policy.

- When both facing layout edges are bounded, the split and extraction windows
  preserve their uncertainty and barcode blocks.
- When only one edge is bounded, that edge supplies a conservative allocation.
  The insert and unknown adjacent scaffold remain available for extraction.
  Such records have `HB:Z:layout_edge` in normal and debug FASTQ/A output. This
  tag does not claim an exactly located circle junction or scaffold removal.
- When neither side locates a safe junction, existing fold/strand evidence may
  resolve it. Otherwise, adjacent units receive no forced barcode assignment.
  The original parent sequence and quality are saved once in
  `<prefix>_boundary_unresolved.fq.gz` or `.fa.gz`, with
  `HU:Z:boundary_unresolved`. Other independently resolved units of that same
  parent may still appear in the primary output. These raw parent copies are
  not additional demultiplexed molecules.

Explicitly unresolved coordinates take precedence over
`--concat-hmm-abstain=legacy`, including when a construct cap truncates the
returned segment list. The ordinary HMM abstention policy otherwise remains
unchanged. `--concat-hmm` remains opt-in.

`rad demux --whitelist-load auto` is the default whitelist loading policy.
For replayable regular input with one non-joint barcode class and a single
catalog file of at least 32 MiB on disk, RAD first collects native barcode
queries using the same prepared layout, read limit and extraction route as the
final run. Forward/reverse aliases count as one barcode class. RAD then streams
the source catalog, retaining only targets reachable by native global lookup:

- Raw barcode identities and their reverse complements.
- Exact same-length barcode sequences in the expanded correction window.
- Substitution neighbors within the runtime `-m` correction radius.

The relevance index uses packed, disjoint seeds followed by full Hamming
distance confirmation. A shared seed alone never retains a catalog entry.
Source rows are streamed rather than stored, and index/retained-entry memory
scales with the queries and relevant targets. The entire catalog must still be
read once; this is not an on-disk random-access index.

The original source cardinality determines whether a catalog is global or
true. A small retained subset of a large global source remains global, so
reducing memory does not switch correction algorithms. Small true catalogs
keep their complete native indel/infix search. Dual catalogs, joint barcodes,
multiple barcode classes and unsupported long queries retain full loading.
This also avoids changing downstream barcode queries when correction of one
barcode moves another barcode's extraction coordinates.

Use `--whitelist-load observed` to request query discovery for an eligible
smaller file (including compressed files below the automatic size threshold).
Unsupported input/layouts are rejected in this explicit mode. A small source
still remains a complete true whitelist. `--whitelist-load full` retains the
original loader and is useful for controlled comparisons. The existing
distinction between `-m` runtime mutation lookup and `-M` whitelist setup is
unchanged.

Streaming logs report source rows, original role, retained entries, observed
keys, correction radius and index bytes. The demux log records the requested
loading mode and whether query discovery ran. Relevance is a lookup
optimization, not a tissue or spatial-density mask.

Run the focused regression checks with:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target rad concat_hmm_boundary_test concat_hmm_router_test \
  whitelist_relevance_test whitelist_loader_test
ctest --test-dir build --output-on-failure
```

The tests cover known insert boundaries in both orientations, mixed and
unlocated junctions, retention across construct caps, sequence/quality
preservation, indexed-versus-exhaustive target selection, and native full versus
selective lookup behavior. They distinguish boundary correctness from the
separate question of whether an internal adapter-like sequence was correctly
classified as another capture unit.

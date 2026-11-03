# Concatemer boundaries and whitelist loading

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
  preserve their uncertainty and barcode blocks. An overlap of up to five bases
  is accepted, matching the decoder's existing microhomology allowance. Both
  children retain the shared interval; larger conflicts and cuts outside the
  observed-marker interval remain unresolved.
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

`rad demux` loads the supplied whitelist in full. Single-catalog imports stream
parsed rows into resident storage after bounded role classification. Dual-catalog
imports consume their staging sets during population, avoiding a second complete
retained copy. These loading changes preserve barcode identities, true/global
whitelist roles and the existing correction policy.

The distinction between `-m` runtime mutation lookup and `-M` whitelist setup is
unchanged.

Run the focused regression checks with:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target rad concat_hmm_boundary_test concat_hmm_router_test \
  whitelist_loader_test packed_sequence_test barcode_storage_test
ctest --test-dir build --output-on-failure
```

The tests cover known insert boundaries in both orientations, mixed and
unlocated junctions, retention across construct caps, sequence/quality
preservation, whitelist formats and role classification, shared catalog loading,
and barcode packing and matching across widths. They distinguish boundary
correctness from the separate question of whether an internal adapter-like
sequence was correctly classified as another capture unit.

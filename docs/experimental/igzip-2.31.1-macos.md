# ISA-L igzip 2.31.1 on Apple Silicon

Test dates: 2026-07-26 and 2026-07-27

Status: concluded experiment; not adopted on `dev`

## Question

Could ISA-L's `igzip` CLI replace pigz as RAD's gzip subprocess on Apple
Silicon, and why did the Linux build work while the initial macOS build
corrupted output or hung?

## Verdict

The Linux/macOS difference was real, but it was not a different compression
algorithm or a missing optimization flag. ISA-L 2.31.1's threaded `igzip` CLI
declares a static pthread pool without initializing its mutex or condition
variable, then publishes completed jobs without synchronization. Zero-filled
pthread objects happen to match glibc's static initializers on Linux. Darwin's
pthread objects require nonzero signatures, so the same invalid assumption
fails on macOS. The remaining job-publication races make the upstream behavior
undefined even where the zero-filled objects appear to work.

A Darwin probe made the distinction explicit: lock/unlock and
signal/broadcast operations on zero-filled pthread storage returned `EINVAL`
(22), while objects created with `PTHREAD_MUTEX_INITIALIZER` and
`PTHREAD_COND_INITIALIZER` succeeded.

The private validation fork initializes those objects, publishes completed
jobs under the pool mutex, uses the condition variable while draining jobs in
order, and synchronizes worker shutdown. It preserves the ISA-L library API
and changes only `programs/igzip_cli.c`.

The patched compressor was substantially faster in isolation: median level-1
compression fell from 2.94 seconds with pigz-4 to 1.59 seconds with igzip-8, a
45.9% reduction. That did not fully translate to an all-igzip RAD run because
the `igzip` CLI's decompressor is single-threaded and slower on this Mac. The
best observed configuration was pigz for input plus igzip-8 for output: 22.87
seconds median versus 23.32 seconds for pigz in both directions, a tentative
1.9% improvement.

That end-to-end margin is small relative to host noise. pigz therefore remains
RAD's general default; the hybrid is the configuration worth revisiting on
larger production inputs.

## Provenance and test system

- Apple M2 Mac mini (`Mac14,3`), 8 cores, 8 GiB RAM
- macOS 26.2
- Homebrew LLVM 22.1.2
- pigz 2.8
- ISA-L tag `v2.31.1`, commit
  `c387163fcbca62701d646149564c550c78a4f985`
- RAD `dev` base commit
  `8426b9f72d696d3c7c2108cecf52f5e4a018af52`
- RAD experiment commits `4968062` and `4788982`
- Private ISA-L tag `v2.31.1-rad-pthread-sync1`, commit
  `89e894cb697e0dbff2256afe274cbdab3d35ee7f`
- Patch-only commit `3b7fa48ba6958bf0bb0887668c224aac7e15b98a`
- Patched CLI marker
  `igzip command line interface 2.31.1-rad-pthread-sync1`

The build forced ISA-L's `aarch64` target so `Makefile.unx` would not silently
select portable C when `uname -m` reported Apple's `arm64` spelling. Homebrew
LLVM was used because Apple's clang in this environment could not assemble all
of ISA-L 2.31.1's AArch64 CRC sources. Building the same unpatched source with
Autotools did not fix the corruption; synchronizing the pthread pool did.

## Correctness validation

The final macOS arm64 binary SHA-256 was:

```text
4f8e88e3f7c1869f1aeba04a53cab9f124b78ed766e772fb3b1f2fe0e444e1c0
```

Focused validation covered 381 cases with no corruption, hangs, or unexpected
exit status:

- levels 0 through 3 at `-T1`, `-T2`, `-T3`, `-T4`, and `-T8`;
- empty, one-byte, 1 MiB, 15 MiB, and job-ring boundary inputs;
- file and standard-input modes, including slow chunked input;
- sequential runs, concurrent processes, and multi-file pool reuse; and
- truncated streams and incorrect CRC/ISIZE controls.

Every valid stream passed Apple gzip, patched igzip, and pigz 2.8, followed by
decompressed SHA checks. A separate 257 MiB incompressible test passed at every
thread count, and a `-T8` ThreadSanitizer build reported no race.

A final 768 MiB matrix added 60 large runs: levels 0 through 3 crossed with
`-T1`, `-T2`, `-T3`, `-T4`, and `-T8`, repeated three times per cell. All 60
passed all three decoders and the complete decompressed SHA check, with zero
timeouts, repeat-SHA mismatches, or size mismatches.

After preserving the patch in the private fork, its tagged source was rebuilt
from a clean archive and passed an additional 40-case matrix spanning two
32 MiB inputs, levels 0–3, and thread counts 1, 2, 3, 4, and 8, plus an
eight-thread standard-input round trip. The rebuilt binary reproduced the
same version marker and SHA-256.

## Compression microbenchmark

The codec test used the decompressed 500K scTagger fixture (1.6 GiB), one warm
run per arm, and six timed runs in balanced permutation order. Every output
passed Apple `gzip -t`.

| Backend | Threads | Timed range | Median | Output bytes |
| --- | ---: | ---: | ---: | ---: |
| pigz 2.8, level 1 | 4 | 2.79–3.03 s | 2.94 s | 307,964,008 |
| patched igzip, level 1 | 4 | 1.45–2.06 s | 1.62 s | 304,367,056 |
| patched igzip, level 1 | 8 | 1.43–2.22 s | 1.59 s | 304,367,056 |

Patched igzip-8 was 45.9% faster than pigz-4 and produced output 1.17% smaller.
Four and eight igzip threads were close, suggesting the fixture reached a
storage or serial-overhead limit before eight workers were fully useful.

## RAD end-to-end benchmark

The end-to-end series used the 98,873-input-read scTagger fixture,
`rad demux -l sctagger -t 4`, one warm run per arm, and four timed repetitions.
The timed order was the balanced rotation `ABCD / BCDA / CDAB / DABC`.

| Read / write backend | External threads | Timed range | Median | Median output bytes |
| --- | ---: | ---: | ---: | ---: |
| pigz / pigz | 4 | 21.60–24.82 s | 23.32 s | 33,419,564 |
| igzip / igzip | 4 | 23.71–26.26 s | 24.86 s | 32,992,049 |
| igzip / igzip | 8 | 23.51–25.38 s | 24.02 s | 32,991,820 |
| pigz / patched igzip | 8 | 22.83–23.09 s | 22.87 s | 32,990,440 |

All-igzip at eight threads was 3.0% slower than pigz end to end despite its
much faster writer. The hybrid was 1.9% faster and produced output 1.28%
smaller. The hybrid runs were tightly grouped, but host load averaged around
7.5 from `cloudd`, GlobalProtect, WindowServer, and an unrelated rsync. The
small end-to-end win is promising rather than conclusive.

Each run reported 101,045 processed records after concatenated inputs were
split, with 34,032–34,037 passing `SigString` objects. The FASTQ files contain
76 additional records because each remaining `read_type == "concatenate"`
object intentionally writes one forward and one reverse record while the
summary counter increments once. That backend-independent difference is not a
compression error. Every gzip stream passed both Apple gzip and patched igzip.

## Historical unpatched baseline

| Operation | Backend | Threads | Median | Result |
| --- | --- | ---: | ---: | --- |
| Compress, level 1 | pigz 2.8 | 4 | 3.16 s | valid |
| Compress, level 1 | upstream igzip 2.31.1 | 1 | 3.80 s | valid |
| Compress, level 1 | upstream igzip 2.31.1 | 2 | 2.13 s | valid |
| Compress, level 1 | upstream igzip 2.31.1 | 3 or 4 | — | corrupt/hang |

The apparently safe two-thread result was not a portability boundary; it was
another outcome of the same races. The earlier all-igzip RAD median was 27.03
seconds versus 23.37 seconds for pigz.

## Reproduction

The source is retained in the private `indianewok/isa-l` fork at
`v2.31.1-rad-pthread-sync1`. Build that tag with:

```bash
make -f Makefile.unx -j8 \
  CC=/opt/homebrew/opt/llvm/bin/clang \
  DEBUG= \
  host_cpu=aarch64 \
  arch=aarch64 \
  version=2.31.1-rad-pthread-sync1 \
  programs/igzip
```

The RAD integration remains isolated on the local `experimental/igzip` branch.
It launches igzip as a subprocess rather than linking ISA-L and uses separate
read/write backend controls. The fastest observed combination was:

```bash
RAD_GZIP_READ_BACKEND=pigz \
RAD_GZIP_WRITE_BACKEND=igzip \
RAD_GZIP_THREADS=8 \
build-igzip/rad demux ...
```

## Decision and follow-up

- Keep pigz as RAD's default.
- Retain the synchronized ISA-L source and immutable tag privately.
- Do not propose the patch upstream until the private fork's build and
  regression coverage are ready for external review.
- Revisit the pigz-read/igzip-write hybrid on larger production data where
  output compression is a larger share of total runtime.

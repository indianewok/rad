# RAD (Read-structure Agnostic Demultiplexer)

RAD is a read-structure agnostic demultiplexer for dealing with long-read sequencing. TL;DR: all you should have to do is define the read structure if you've got a super-wonky custom sequencing format, pass it to RAD, see if it preps nicely, and let it demultiplex! I've used it with a whole bunch of stuff--weirdest so far has been long-read targeted enrichment of BCR/TCR from Visium HD data, so if you've got weirder than that I'd love to see whether RAD works for you!

## Docs map

| Need | File |
| --- | --- |
| Install + first run (`prep -> demux -> reformat`) | [`docs/installation.md`](docs/installation.md) |
| Command flags + output files | [`docs/cli-reference.md`](docs/cli-reference.md) |
| Layout + whitelist details (origins, sizes, pairings) | [`docs/layouts-and-whitelists.md`](docs/layouts-and-whitelists.md) |
| What RAD does under the hood | [`docs/architecture.md`](docs/architecture.md) |
| End-to-end smoke test on simulated data | [`test_data/README.md`](test_data/README.md) |

## Pipeline at a glance

```mermaid
flowchart LR
    A["rad prep (layout + optional position map)"] --> B["rad demux (layout-aware extraction + filtering)"]
    B --> C["rad reformat (header rewrite and/or CB split)"]
```

## Minimal run

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

build/rad prep -l five_prime --position-map -q reads.fq.gz -o run/demo
build/rad demux -l five_prime -q reads.fq.gz -o demo -d run -t 8
build/rad reformat -q run/demo.fq.gz --split-bc -o run/by_barcode -t 8
```

## Repo layout

```text
.
├── src/                  # CLI entrypoints
├── include/rad/          # core pipeline + algorithms
├── resources/
│   ├── read_layout/      # bundled layout templates
│   └── wl/               # bundled whitelist resources
├── docs/                 # user + methods docs
└── CMakeLists.txt
```

## Operational notes

- RAD looks for `resources/` relative to the executable location, then `./resources`.
- `pigz` is optional, but it usually improves gzip throughput a lot.
- `rad demux --bc_split` is shown in help, but split output is handled by `rad reformat --split-bc` in the current build.
- `rad_config set/rm` is process-local in the current build, so those updates won't persist across separate invocations.
- After big source/header edits, do a clean rebuild:

```bash
cmake --build build --clean-first
```

## License

[`LICENSE`](LICENSE)

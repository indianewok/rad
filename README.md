# RAD (Read-structure Agnostic Demultiplexer)

RAD is a read-structure agnostic demultiplexer for layout-aware sequencing demultiplexing.

Short version: define the read structure, map positions, run demux, then reformat/split as needed.

## Docs map

| Need | File |
| --- | --- |
| Install/build | [`docs/installation.md`](docs/installation.md) |
| First run (`prep -> demux -> reformat`) | [`docs/quickstart.md`](docs/quickstart.md) |
| Exact command/flag behavior | [`docs/cli-reference.md`](docs/cli-reference.md) |
| Layout + whitelist details (origins, sizes, pairings) | [`docs/layouts-and-whitelists.md`](docs/layouts-and-whitelists.md) |
| What RAD is doing under the hood | [`docs/architecture.md`](docs/architecture.md) |
| Output file contract | [`docs/output-files.md`](docs/output-files.md) |
| Failure diagnosis | [`docs/troubleshooting.md`](docs/troubleshooting.md) |

## Pipeline at a glance

```mermaid
flowchart LR
    A["rad prep\n(layout + optional position map)"] --> B["rad demux\n(layout-aware extraction + filtering)"]
    B --> C["rad reformat\n(header rewrite and/or CB split)"]
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
├── docs/                 # user/methods docs
└── CMakeLists.txt
```

## Operational notes

- RAD resolves `resources/` relative to executable location, then `./resources`.
- `pigz` is optional, but usually worth it for gzip throughput.
- After big source/header edits, do a clean rebuild:

```bash
cmake --build build --clean-first
```

## Current caveats

- `rad demux --bc_split` is shown in help, but split execution there is stubbed; use `rad reformat --split-bc`.
- `rad_config set/rm` is process-local right now (not persistent across independent invocations).

## License

[`LICENSE`](LICENSE)

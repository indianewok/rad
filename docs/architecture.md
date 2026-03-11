# Architecture

This is the short methods-level view of what RAD does internally.

## Core pieces

| Component | Job |
| --- | --- |
| `src/main.cpp` | CLI dispatch (`prep`, `demux`, `reformat`) |
| `include/rad/read_layout.hpp` | layout parsing/normalization + position-map logic |
| `include/rad/barcode_correction.hpp` | whitelist import and barcode set structures |
| `include/rad/sigstring.hpp` | per-read extraction/alignment/filtering (`sigalign`) |
| `include/rad/io_streaming.hpp` | chunk streaming, pigz/gzip I/O, write queues |
| `src/config_tools.cpp` | `rad_config` alias interface |

## Command-level flow

```mermaid
flowchart TD
    A["CLI parse"] --> B{"command"}
    B --> C["prep"]
    B --> D["demux"]
    B --> E["reformat"]

    C --> C1["load + normalize layout"]
    C1 --> C2["optional misalignment sampling"]
    C2 --> C3["write layout + position map"]

    D --> D1["resolve/import layout + map"]
    D1 --> D2["load whitelist maps"]
    D2 --> D3["chunked sigalign"]
    D3 --> D4["write FASTQ + optional debug outputs"]

    E --> E1["stream reads"]
    E1 --> E2["header collapse and/or CB split"]
    E2 --> E3["write rewritten/split outputs"]
```

## Layout processing (`prep_new_layout`)

What happens:

1. parse CSV rows
2. normalize fields (`class`, `class_id`, lengths)
3. inject sentinels (`seq_start`, `seq_stop`)
4. auto-generate reverse-complement elements unless rows are single-sided
5. build ordered indexes that `demux` uses later

If `--position-map` is set, RAD samples reads, computes misalignment stats, then writes `_layout.csv` and `_position_map.csv`.

## Demux pipeline (`sigalign`)

Runtime flow:

1. load or build layout/map
2. resolve/load whitelist data (`true_bcs` + `global_bcs`)
3. stream reads in chunks
4. process each read with:
   - `sigalign_static`
   - `sigalign_variable`
   - `sigalign_filter`
5. emit passing reads and optional debug channels
6. write whitelist summaries

Per-read logic:

```mermaid
flowchart LR
    A["input read"] --> B["static anchor alignment"]
    B --> C["variable segment extraction"]
    C --> D["whitelist/correction checks"]
    D --> E{"pass?"}
    E -->|yes| F["write output read"]
    E -->|no| G["filtered (debug only if enabled)"]
```

## Whitelist model

- `true_bcs`: stricter accepted set
- `global_bcs`: broader correction candidate set

Import behavior:
- one source: RAD assigns based on set size policy
- two sources: smaller usually maps to `true`, larger to `global`

## Parallelism and I/O

- read path: chunked streaming, `pigz` if available
- compute path: OpenMP across chunks/reads
- write path: buffered/asynchronous writer path

Main performance knobs:
- `--threads`
- `--chunk_size`
- `pigz` availability/config

## Practical implementation notes

- `demux --bc_split` appears in help, but split execution is handled in `reformat --split-bc` in the current build.
- `rad_config set/rm` is process-local in the current build, so updates won't persist across independent invocations.

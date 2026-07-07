# Installation and first run

This page takes you from a fresh machine to your first demultiplexed file. If you've never built a C++ tool before, you can copy-paste straight down the page — nothing here assumes prior experience.

## 1) What RAD needs

RAD is a small C++ program. Building it needs a handful of standard, free tools — that's the whole list:

| Dependency | What it's for |
| --- | --- |
| **CMake ≥ 3.18** | runs the build |
| **A C++17 compiler** | compiles RAD (clang or gcc, recent enough for C++17) |
| **OpenMP runtime** | multithreading (the `-t` flag) |
| **Boost** (`filesystem`, `iostreams`) | file paths + compressed I/O |
| **zlib** | reading/writing `.gz` files |
| **pigz** | parallel gzip — RAD uses it for all `.gz` I/O and it's dramatically faster than plain zlib |

Everything else RAD needs (edlib, ssw, kseq, the CSV parser…) is bundled **inside the repo**, so there's nothing else to chase down.

## 2) Install the dependencies

Pick the **one** line that matches your machine.

**On a Mac (Homebrew):**

```bash
brew install cmake boost libomp pigz llvm
```

(`llvm` gives you a clang with clean OpenMP support, which is the smoothest setup on macOS.)

**On Linux, an HPC cluster, or anywhere you don't have admin/`sudo` rights (conda/mamba):**

```bash
mamba create -n rad -c conda-forge cxx-compiler cmake boost zlib pigz
mamba activate rad
```

This puts the compiler **and** every dependency into one self-contained environment — no `sudo`, which is the usual situation on a shared cluster.

> **Don't have `mamba`/`conda` yet?** Install [Miniforge](https://github.com/conda-forge/miniforge) first — it's the standard, no-admin-rights way to get `mamba`. Then come back and run the line above.

## 3) Build RAD

```bash
git clone https://github.com/indianewok/rad.git
cd rad

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

On a Mac, point CMake at Homebrew's clang for the cleanest OpenMP linkage:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++"
cmake --build build -j
```

When it finishes you'll have a `build/rad` binary. Confirm it works:

```bash
build/rad --help     # lists every command (prep, scan-wl, demux, reformat, list, modify, help)
build/rad list       # shows the bundled layouts + whitelists RAD can find
```

If `rad list` prints layout and whitelist paths, your build found its bundled `resources/` and you're ready to demultiplex.

> **Typical build time:** compiling RAD takes **~15–20 seconds** from a clean checkout on a modern desktop (measured on an Apple M2 Max with `-j`). Installing the dependencies above adds a few minutes, depending on your platform and network.

> **On Bioconda:** RAD is available on Bioconda — `mamba install -c bioconda rad` installs the binary and all dependencies (pigz included), so you can skip the manual dependency setup and build above.

## 4) Your first run

Let's run RAD end-to-end on a small simulated dataset so you can see it work.

**Get the test reads** (~99k simulated long reads, 50 known cells):

```bash
curl -L -o test.fq.gz \
  https://github.com/indianewok/rad/releases/download/test-data-v1/shuffled_S_lr_synth.fq.gz
```

**(Optional) Eyeball the read structure first.** RAD ships layout templates by key — here we use `sctagger`:

```bash
build/rad prep -l sctagger --read-layout
```

This just prints how RAD will interpret the reads; nothing is written. (For a custom format, pass a CSV path instead of a key: `-l /path/to/layout.csv`.)

**Demultiplex:**

```bash
build/rad demux -l sctagger -q test.fq.gz -d run -o demo -t 4
```

This writes `run/demo.fq.gz` with corrected `CB:Z` (cell barcode) and `UB:Z` (UMI) tags in each header. On this dataset you should see **~33,700 reads pass** the filter — if you do, everything is working. The run takes **~20 seconds** on a modern desktop (Apple M2 Max, `-t 4`).

The `sctagger` layout carries a default whitelist, so this corrects against the full 10x list. On your own data you usually want to **discover the barcodes actually present first** — run `rad scan-wl` (see [`test_data/README.md`](../test_data/README.md)), or add `-A/--auto-wl` to the `demux` line to fold that scan into a single command. (By default `UB:Z` is reverse-complemented on reverse reads so it matches `CB:Z`'s orientation; `--no-umi-rc` keeps UMIs exactly as extracted.)

**(Optional) Tidy up the output.** Rewrite headers, or split into one file per cell barcode:

```bash
build/rad reformat -q run/demo.fq.gz --reformat-header -t 4
build/rad reformat -q run/demo.fq.gz --split-bc -o run/by_barcode -t 4
```

That's the pipeline: **discover barcodes with `scan-wl`** (or `demux -A/--auto-wl`), **`demux`** against that whitelist, then optional **`reformat`**. (`prep` is a preliminary layout/position-map step; it's not part of the per-run barcode path.) For the full smoke test — including de-novo barcode discovery with `scan-wl` and exact expected numbers — see [`test_data/README.md`](../test_data/README.md).

## 5) Runtime notes

- RAD finds its `resources/` folder by climbing up from the `rad` binary, then `./resources`. If you move the binary, move `resources/` with it (or just run from the repo) — or point `RAD_RESOURCES` at it.
- Handy environment knobs:
  - `RAD_RESOURCES=/path/to/resources` — override where RAD looks for its bundled resources (highest precedence)
  - `RAD_PIGZ=/path/to/pigz` — point at a specific pigz
  - `RAD_NO_PIGZ=1` — disable pigz, fall back to zlib
  - `RAD_PIGZ_THREADS=<N>` — cap pigz threads
- After big source/header edits, do a clean rebuild:

```bash
cmake --build build --clean-first
```

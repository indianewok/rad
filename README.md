# Quickstart: Building **rad** from Source

Follow these steps in your terminal to install dependencies and build **rad** in Release mode.

---

## 1. Prerequisites

* **macOS** or **Linux**
* A working **bash** (or similar) shell
* **Homebrew** (macOS) or **conda** (Linux/HPC)

---

## 2. Install Dependencies

### macOS (Homebrew)

```bash
brew install llvm libomp boost cmake pigz zlib
```

### Linux / HPC (conda)

```bash
conda install -c conda-forge boost-cpp cmake pigz zlib libomp-dev
```

> **Note:** `pigz` is optional but recommended for parallel gzip decompression. Set `RAD_NO_PIGZ=1` to disable at runtime.

---

## 3. Clone the Repository

```bash
git clone https://github.com/indianewok/rad.git
```

---

## 4. Create and Enter the Build Directory

```bash
cd rad
mkdir build && cd build
```

---

## 5. Configure & Build in Release Mode

```bash
# Generate Makefiles for a Release build
cmake -DCMAKE_BUILD_TYPE=Release ..

# Compile
cmake --build .
```

---

## 6. Run

Your executables will now be in `build/`:

```bash
./rad --help
./rad prep --help
./rad demux --help
./rad reformat --help
```

You're all set!

---

## 7. Test Data

A simulated long-read scRNA-seq dataset (scTagger benchmark, 10x 3' v3 chemistry) is available as a GitHub release asset for end-to-end testing.

```bash
# Download (~102 MB)
curl -L -o shuffled_S_lr_synth.fq.gz \
  https://github.com/indianewok/rad/releases/download/test-data-v1/shuffled_S_lr_synth.fq.gz

# Run rad demux with the pre-registered sctagger layout
./rad demux \
  --layout sctagger_sim \
  -i shuffled_S_lr_synth.fq.gz \
  -o test_out/ \
  -t 1
```

The dataset contains 98,873 reads across 50 simulated cells. True cell barcodes are encoded in each read header as `@<16bp_barcode>-<counter>`, so you can validate the `CB:Z` tag in the demultiplexed output directly against the header — no separate truth file required.

See the [release page](https://github.com/indianewok/rad/releases/tag/test-data-v1) for full details.

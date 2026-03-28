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

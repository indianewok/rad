# Quickstart: Building **rad** from Source

Follow these steps in your terminal to install dependencies and build **rad** in Release mode.

---

## 1. Prerequisites

* **macOS** or **Linux**
* A working **bash** (or similar) shell
* **Homebrew** (Package manager for macOS/Linux)

---

## 2. Install Homebrew (if needed)

```bash
which brew
```

* If you see a path (e.g. `/usr/local/bin/brew`), you already have it.
* Otherwise, install from [https://brew.sh/](https://brew.sh/) by following their instructions.

---

## 3. Install Dependencies

```bash
# Install LLVM (for up-to-date clang), OpenMP and Boost
brew install llvm libomp boost
```

> **Tip:** If you already have these, `brew install` will skip or upgrade them.

---

## 4. Clone the Repository

```bash
git clone --recurse-submodules https://github.com/indianewok/rad.git
```

Make a note of the clone path (e.g. `~/projects/rad`).

---

## 5. Create and Enter the Build Directory

```bash
cd path/to/rad          # e.g. cd ~/projects/rad
mkdir build && cd build
```

---

## 6. Configure & Build in Release Mode

```bash
# 6A) Generate Makefiles for a Release build
cmake -DCMAKE_BUILD_TYPE=Release ..

# 6B) Compile!
cmake --build .
```

---

## 7. Run

Your `rad` executable will now be in `build/`. For example:

```bash
./rad --help
```

You’re all set!

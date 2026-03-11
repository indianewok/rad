# Installation

RAD is a CMake/C++ build. Once deps are installed, build steps are the same across platforms.

## Required deps

- CMake `>= 3.18`
- C++17 compiler
- OpenMP runtime
- Boost (`filesystem`, `system`, `iostreams`)
- zlib

Optional:

- `pigz` for faster `.gz` I/O

## Platform-agnostic build

```bash
git clone --recurse-submodules https://github.com/indianewok/rad.git
cd rad

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Quick verification:

```bash
build/rad --help
build/rad_config --help
```

## macOS (Homebrew)

```bash
brew install cmake boost libomp llvm pigz
```

Use Homebrew LLVM for cleaner OpenMP linkage:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$(brew --prefix llvm)/bin/clang" \
  -DCMAKE_CXX_COMPILER="$(brew --prefix llvm)/bin/clang++"
cmake --build build -j
```

## Ubuntu / Debian

```bash
sudo apt update
sudo apt install -y \
  build-essential cmake \
  libboost-all-dev zlib1g-dev \
  libomp-dev pigz
```

## RHEL / Rocky / Alma

```bash
sudo dnf install -y \
  gcc gcc-c++ cmake \
  boost-devel zlib-devel \
  libgomp pigz
```

## Conda/Mamba user-space toolchain

```bash
mamba create -n rad-build -y \
  cmake cxx-compiler compilers \
  boost zlib libgomp pigz
mamba activate rad-build

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

## Runtime path note

RAD expects `resources/` to be reachable relative to executable/CWD search logic.
If binaries are moved, move `resources/` with them (or run from repo tree).

## Runtime env knobs

- `RAD_PIGZ=/path/to/pigz`
- `RAD_NO_PIGZ=1`
- `RAD_PIGZ_THREADS=<N>`

## Clean rebuild

```bash
cmake --build build --clean-first
```

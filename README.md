# ERI-engine

[![CI](https://github.com/ifilot/eri-engine/actions/workflows/build.yml/badge.svg)](https://github.com/ifilot/eri-engine/actions/workflows/build.yml)

## Compilation

Make sure the following packages are installed

```bash
sudo apt-get update
sudo apt-get install -y \
    meson \
    ninja-build \
    g++ \
    python3 \
    libbenchmark-dev \
    libjsoncpp-dev
```

```bash
meson setup build
ninja -C build
```

and test your compilation via

```bash
meson test -C build
```

or using verbose output

```bash
meson test -C build --verbose
```

## Benchmarking

```bash
meson compile -C build bench-eri
```
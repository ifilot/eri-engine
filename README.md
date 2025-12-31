# ERI-engine

[![CI](https://github.com/ifilot/eri-engine/actions/workflows/build.yml/badge.svg)](https://github.com/ifilot/eri-engine/actions/workflows/build.yml)

## Compilation

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
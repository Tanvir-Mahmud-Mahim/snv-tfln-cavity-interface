# Building Meep for the loaded-cavity FDTD

The loaded one-sided cavity (`sim/cavity_fdtd.py`, revision 1) is
simulated with the FDTD code [Meep](https://github.com/NanoComp/meep).
Meep is not on PyPI. The usual route is conda-forge
(`conda install -c conda-forge pymeep`). If that works for you, use it and
skip this file.

Our computer could not reach the conda channels, so we built Meep and its
dependencies from source. This file records exactly what we did, so that
the FDTD numbers can be reproduced with the same versions.

## Versions used

| Package | Version | Source |
|---|---|---|
| Meep | 1.36.0-beta, git commit `f9cc3474bc95365f8e12858747563926017dc83d` | github.com/NanoComp/meep |
| MPB | 1.11.1 (release tarball) | github.com/NanoComp/mpb/releases |
| libctl | 4.5.1 (release tarball) | github.com/NanoComp/libctl/releases |
| Harminv | 1.4.3 (git tag `v1.4.3`) | github.com/NanoComp/harminv |
| Python | 3.13 (virtual environment), numpy 2.4.4, h5py 3.16.0 | |
| System (Ubuntu) | libfftw3-dev 3.3.10, libgsl-dev 2.7.1, libhdf5-dev 1.10.10, liblapack-dev 3.12.0 | apt |

Two things did not work and are worth knowing:

- libctl from the git master branch failed to build (`geom-ctl-io.c`
  missing), so we used the 4.5.1 release tarball.
- Without MPB, the Python module `_meep.so` failed to import with an
  undefined `populate_maxwell_dielectric` symbol. MPB must be built
  before Meep. MPB is also what Meep uses for the eigenmode
  decomposition into the TE<sub>0</sub> mode of the output waveguide.

## Steps

```
export PREFIX=$HOME/meep-install
export PYTHON=/path/to/venv/bin/python     # numpy and h5py installed

# 1. Harminv
cd harminv && sh autogen.sh --prefix=$PREFIX --enable-shared && make -j2 && make install

# 2. libctl 4.5.1 (release tarball)
cd libctl-4.5.1 && ./configure --prefix=$PREFIX --enable-shared --without-guile \
  && make -j2 && make install

# 3. MPB 1.11.1 (release tarball)
cd mpb-1.11.1 && ./configure --prefix=$PREFIX --enable-shared --without-libctl \
  --with-hermitian-eps CPPFLAGS="-I$PREFIX/include" LDFLAGS="-L$PREFIX/lib" \
  && make -j2 && make install

# 4. Meep (git commit above)
cd meep && sh autogen.sh --prefix=$PREFIX --enable-shared --with-python \
  --without-mpi --without-scheme \
  CPPFLAGS="-I$PREFIX/include -I/usr/include/hdf5/serial" \
  LDFLAGS="-L$PREFIX/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial -Wl,-rpath,$PREFIX/lib" \
  PYTHON=$PYTHON && make -j2 && make install
```

## Running

```
export PYTHONPATH=$PREFIX/lib/python3.13/site-packages
export LD_LIBRARY_PATH=$PREFIX/lib
export OMP_NUM_THREADS=1
cd sim
python cavity_fdtd.py bulk 40        # dipole in bulk diamond (LDOS reference)
python cavity_fdtd.py 0 40 1000      # N_out = 0, run until the field has decayed
python cavity_fdtd.py 1 40 300       # N_out = 1 to 3: ring-down of 300 time units
python cavity_fdtd.py 2 40 300
python cavity_fdtd.py 3 40 300
python cavity_fdtd.py 0 50 300       # resolution check
python fdtd_analysis.py              # -> results/cavity_fdtd_summary.json, cavity_loaded.json
```

The arguments are the number of output mirror periods N<sub>out</sub>, the
resolution in pixels per micrometer, and the ring-down time after the
source in Meep time units (micrometers / c). `sim/fdtd_jobs.sh` runs the
same list in sequence. On one core of our shared computer one run at
resolution 40 took about 25 to 70 minutes, and the resolution-50 run
about an hour.

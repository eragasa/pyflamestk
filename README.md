# pyflamestk
pyFlamesTk - Python FLAMES Tookit

The python FLAMES (Florida Laboratory for Materials Engineering Simulation) Toolkit is a collection of functions, and libraries that I wrote while at the University of Florida to study materials science.  By it's very nature, it is an incomplete toolkit and a work in progress.

pyflamestk is written in python3 using packages found in the Anaconda3 distribution of python.  The choice of using Anaconda is driven by the desire to have pyflamestk portable across multiple platforms.  However, not all code is platform agnostic particularly when it comes to dealing with filenames, due to differences between how Windows and Unix-style operating systems behave.

I hope to fully document the system, including providing appropriate references for particular techniques.  However, this is fairly time-consuming so references will normally be provided inline.

Release-readiness work is documented in
[`docs/release-readiness.md`](docs/release-readiness.md). Generated calculator
artifacts staged for operator-managed offline storage are described in
[`OFFLINE_ARTIFACTS.md`](OFFLINE_ARTIFACTS.md).

## Installation boundary

PyFlamestk is an alpha research toolkit. The maintained wheel contains only the
`pyflamestk` package; `dev/`, `examples/`, `tests/`, offline artifacts, and
calculator output are not package content.

Core numerical modules require NumPy and SciPy. Plotting and DAKOTA YAML support
are explicit extras:

```bash
python -m pip install 'pyflamestk[plot]'
python -m pip install 'pyflamestk[dakota]'
```

Installing or importing the package does not install, discover, configure, or
run VASP, LAMMPS, GULP, DAKOTA, MPI, or a scheduler. Examples and development
scripts are historical research material, not portable supported workflows;
review their paths and process effects before any manual use. A source release
makes no claim of calculation reproducibility, numerical verification, or
scientific validation.

Interfaces to Major Packages:
- VASP (calculator for DFT simulations)
- LAMMPS (calculator for Molecular Dynamics Simulations)
- DAKOTA (calculator for VVUQ Simulations)

Primary Functionality:
- Fitting of empirical interatomic potentials via pareto surface
- Simulation of Infrared Spectra, via VASP
- Simulation of Raman Spectra, via VASP

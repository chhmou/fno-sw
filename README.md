# Learning the solution operator of standing waves with a Fourier neural operator

Code for the paper

> C. Mou, S. Ammons, D. Qi, and V. M. Hur, *Learning the solution operator of standing waves with a Fourier neural operator*, arXiv:2610.36747, 2026. https://arxiv.org/abs/2610.36747

A Fourier neural operator (FNO) maps the rest surface profile and a normalized time to the surface elevation, velocity potential, and dipole density of a standing water wave in infinite depth over a quarter period. Two experiments test extrapolation toward the free-fall limit (Experiment I) and the linear limit (Experiment II). A residual of the governing equations is evaluated to assess predictions without reference solutions.

## Contents

- FNO model, training, and evaluation for Experiments I and II
- Error and residual diagnostics
- Exploratory predictions beyond the available data
- Trained model checkpoints

## Requirements

Python 3.10+, PyTorch, NumPy, SciPy, Matplotlib.

```bash
pip install torch numpy scipy matplotlib
```

## Data

The reference solutions are not included and are available upon reasonable request. Place them in `data/` or set

```bash
export STANDING_WAVE_DATA=/path/to/data
```

## Authors

- Changhong Mou
- Savana Ammons
- Di Qi
- Vera Mikyoung Hur

## Citation

If you use this code, please cite:

```bibtex
@article{mou2026standing,
  title   = {Learning the solution operator of standing waves with a {F}ourier neural operator},
  author  = {Mou, Changhong and Ammons, Savana and Qi, Di and Hur, Vera Mikyoung},
  journal = {arXiv preprint arXiv:2610.36747},
  year    = {2026}
}
```

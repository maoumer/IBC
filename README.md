# Interpolation-Inspired Barrier Certificates (IBCs)
SOS and SMT solver code to search for IBCs with examples

This repo contains the code for our paper:

1. M. A. Oumer, V. Murali, A. Trivedi and M. Zamani, "Safety Verification of Discrete-Time Systems via Interpolation-Inspired Barrier Certificates," in IEEE Control Systems Letters, 2024. ([link](https://doi.org/10.1109/LCSYS.2024.3521356)).

## Problem Statement
We address the problem of safety verification of discrete-time dynamical systems by introducing the notion of interpolation for barrier certificates (BCs). The search for BCs is typically automated by first fixing a template of functions and then using sum-of-squares (SOS) programming or satisfiability modulo theory (SMT) solvers to find them. Unfortunately, it may not be possible to find a single function in this fixed template. To tackle this challenge, we propose the notion of interpolation-inspired barrier certificate (IBC). Instead of a single function, an interpolation-inspired barrier certificate consists of a set of functions such that the union of their sublevel sets over-approximate the reachable set of states. An additional practical benefit of IBCs is in addressing the complexity of SOS programs by allowing multiple lower degree functions to act as proofs of safety. We use SOS in our paper but sample SMT code is included here as well.

## Files
A. `sos_ibc.jl` contains the main SOS Julia code to search for both standard BCs and IBCs for any system defined under the `systems` folder.

B. `z3_smt_barrier.py` contains the SMT Z3 Python code to search IBCs for a 1D system given inside the file, and plot them once found.

C. `check_plot_barrier.py` contains the SMT Z3 Python code to verify the results from `sos_ibc.jl` and plot them. This is needed sometimes when the SOS solver does not terminate with an "OPTIMAL" status.

D. `systems` folder contains the systems of interest and their specifications in Julia files.

E. `media` folder contains the plots generated from running the Python files in ".svg" format.

## Running
1. Download [Julia](https://julialang.org/downloads/) and install all the packages at the top of `sos_ibc.jl` file.
2. Install [MOSEK](https://docs.mosek.com/11.0/install/index.html) for Julia and set it up according to the instructions. Otherwise, feel free to change the SOS solver inside the Julia file to work with other solvers.
3. Download [Python](https://www.python.org/downloads/) and install the [Z3](https://github.com/Z3Prover/z3) library.


Bibtex of our paper:
```
@ARTICLE{10811976,
  author={Oumer, Mohammed Adib and Murali, Vishnu and Trivedi, Ashutosh and Zamani, Majid},
  journal={IEEE Control Systems Letters}, 
  title={Safety Verification of Discrete-Time Systems via Interpolation-Inspired Barrier Certificates}, 
  year={2024},
  volume={8},
  number={},
  pages={3183-3188},
  doi={10.1109/LCSYS.2024.3521356}}
```

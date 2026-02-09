# Installation

AxiSEM3D is a compiled code. This means that it is not enough just to
download it into one place, you must also build and compile it. This will produce a binary, and the binary is what you actually ‘run’ to produce
simulations.

The installation of **AxiSEM3D** consists of three parts:

- the mesher  
- the solver  
- a few tools for pre- and post-processing  

## System Requirements

- Unix-like operating system (Linux, macOS)  
- A C++ compiler supporting the **C++17** standard  
- `conda` (Anaconda or Miniconda)  
- `cmake` **4.0 or above** (available from conda-forge)  
- MPI implementation (serial build possible but not recommended)

---

```{toctree}
---
maxdepth: 1
---
mesher.md
solver/index.md
```

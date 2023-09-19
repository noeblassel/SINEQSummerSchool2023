# ANR SINEQ Summer School 2023
## Sampling high-dimension probability measures: applications in (non)equilibrium molecular dynamics and statistics

This repository contains the various files for the hands-on sessions of the SINEQ summer school @ CERMICS, September 2023. The schedule can be found [here](https://sites.google.com/view/aleiac/anr-sineq/summer-school-mol-dyn-on-julia).
The files for each sessions are given below.
* Monday, Sept. 25: [Discretization of Langevin dynamics and its Metropolization](notebooks/langevin.ipynb)
* Tuesday, Sept. 26: [Nonequilibrium systems and coupling methods]()
* Wednesday, Sept. 27, morning: [Accelerated molecular dynamics methods]()
* Wednesday, Sept. 27, afternoon: [Adaptive multilevel splitting methods]()
* Thursday, Sept. 28: [Computing average properties in molecular dynamics with Molly.jl](notebooks/molly_average.ipynb)
* Friday, Sept. 29: [Computation of transport coefficients with Molly.jl]()

## Various installations
We provide Julia notebooks in this repository. The three main ingredients required are Python (>=3.3 or 2.7, needed for Jupyter Notebook), Jupyter Notebook and Julia.

### Installing Python, Jupyter and scientific computing libraries
* If Python is not installed on your system, install it using Anaconda. It conveniently installs Python, the Jupyter Notebook, and other commonly used packages for scientific computing and data sciences. You can download Anaconda [here](https://www.anaconda.com/download) (install it with the latest Python version). Then install the version of Anaconda which you downloaded, following the instructions on the download page. Once installed, you can simply launch Jupyter Notebook by running in a shell:
    ```bash
    jupyter-notebook
    ```

* If you already have Python installed (>= 3.3 or 2.7), run the following commands in a shell (use pip instead of pip3 if you use 2.7):
    ```bash
    pip3 install --upgrade pip
    pip3 install jupyter
    ```

### Installing Julia
Please refer to https://github.com/JuliaLang/juliaup where you will find all the needed resource to install Julia on your OS.
Once installed, open a Julia REPL by typing `julia` in a shell. Then type:
```julia
using Pkg
Pkg.add("IJulia")
```
so that Julia can be run inside Jupyter Notebooks.

### Verifying your installations and enable multithreading in Jupyter
Run the notebook `installation.ipynb` using the following in a shell:
```bash
jupyter-notebook installation.ipynb
```

### Installing PyMol
For the Molly session, the molecular visualization software PyMol is recommended. You can install it by following the instructions provided [here](https://pymol.org/2/support.html?#installation).

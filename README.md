# ANR SINEQ Summer School 2023
## Sampling high-dimension probability measures: applications in (non)equilibrium molecular dynamics and statistics

This repository contains the various files for the hands-on sessions of the SINEQ summer school @ CERMICS, September 2023. The schedule can be found [here](https://sites.google.com/view/aleiac/anr-sineq/summer-school-mol-dyn-on-julia).
The files for each sessions are given below.
* Monday, Sept. 25: [Discretization of Langevin dynamics and its Metropolization](notebooks/langevin.ipynb)
* Tuesday, Sept. 26: [Nonequilibrium systems and coupling methods]()
* Wednesday, Sept. 27, morning: [Accelerated molecular dynamics methods]()
* Wednesday, Sept. 27, afternoon: [Adaptive multilevel splitting methods]()
* Thursday, Sept. 28: [Computing average properties in molecular dynamics with Molly.jl]()
* Friday, Sept. 29: [Computation of transport coefficients with Molly.jl]()

## Various installations
We provide Julia notebooks in this repository. The three main ingredients required are Python (>=3.3 or 2.7, needed for Jupyter Notebook), Jupyter Notebook and Julia.

### Installing Julia
Please refer to https://github.com/JuliaLang/juliaup where you will find all the needed resource to install Julia on your OS.
Once installed, open a Julia REPL by typing `julia` in a shell. Then type:
```julia
using Pkg
Pkg.add("IJulia")
```
so that Julia can be run inside Jupyter Notebooks.

### Installing Jupyter
* If you already have Python installed (>= 3.3 or 2.7), run the following commands in a shell (use pip instead of pip3 if you use 2.7):
```bash
pip3 install --upgrade pip
pip3 install jupyter
```
* If Python is not installed on your system, you have two options:
    * if you plan to use Python often in the future, we recommend installing it with Anaconda. It convieniently installs Python, the Jupyter Notebook, and other commonly used packages for scientific computing and data sciences. You can download Anaconda [here](https://www.anaconda.com/download) (install it with the latest Python version). Then install the version of Anaconda which you downloaded, following the instructions on the download page. Once installed, you can simply launch Jupyter Notebook by running in a shell:
    ```bash
    jupyter notebook
    ```
    * if you do not plan to use Python that much, then after installing Julia and adding IJulia, simply run:
    ```julia
    using IJulia
    notebook()
    ```
    The first time you run `notebook()`, it will prompt you whether it should install Jupyter. You can hit enter to have it use the `Conda.jl` package to install a minimal Python+Jupyter distribution using Miniconda (instead of Anaconda) that is private to Julia (it won't appear in your `PATH`). You will then be able to launch notebooks with this method. Alternatively, you may write directly in a shell
    ```bash
    jupyter notebook
    ```
    however this will only work if `jupyter` is in you `PATH`. To know where Conda installed `jupyter`, run in a Julia REPL:
    ```julia
    import Conda
    Conda.SCRIPTDIR
    ```
    Then add this location to the `PATH` environment variable. 
    
    To add a folder to `PATH` on Windows, open a shell as administrator, and run:
    ```bash
    setx Path "%Path%;C:\path\to\folder"
    ```
    If you wish to modify the `PATH` of every users, add `/M` between `setx` and `Path`. 
    
    On Mac OS X, modify your `.bash_profile` file (located in your home directory) by adding to a new line:
    ```bash
    PATH="$PATH:the/absolute/path/to/your/folder"
    ```
    Save the file, then close the text editor, then run `source ~/.bash_profile` (or restart your terminal). 
    
    On Linux, modify the `.bashrc` file in the same manner as for Mac OS X.


### Installing Julia packages for the hands-on sessions
In the folder containing this `README.md`, execute the script:
```bash
julia setup.jl
```
which will add the necessary packages to your Julia environment. As this is a lengthy process, be sure to do it **before** the start of the summer school.
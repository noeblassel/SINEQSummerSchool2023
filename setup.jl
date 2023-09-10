import Pkg

for pak in ["IJulia","Plots","LaTeXStrings","Cubature","StatsBase", "Molly","StaticArrays","Unitful","Bio3DView","KernelDensity","Measurements","Zygote"]
    Pkg.add(pak),
end
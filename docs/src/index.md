```@meta
CurrentModule = StirredReactor
```

# StirredReactor
The Cstr model simulates a continuous stirred tank reactor. This type of reactor model may be used to simulate a packed bed reactor, when the length to diameter ratio is small, or when a packed bed reactor is operating under differential conditions. The code integrates the following governing equation 


Documentation for [StirredReactor](https://github.com/vinodjanardhanan/StirredReactor.jl).


# Installation
To install the package, use the following commands in the julia REPL
```julia
julia> using Pkg
julia> Pkg.add("StirredReactor")
```

# General interfaces

```@index
```

```@autodocs
Modules = [StirredReactor]
```

# Governing equations
```math
V \frac{dC_k}{dt} = \dot{q}_{in} C_{k,in} - \dot{q}_{out} C_k + \dot{s_k} A_s + \dot{\omega_k} V
```
Here $C_k$ is the concentration of species $k$, $V$ is the reactor volume (m$^3$), $\dot{q}$ is the volumetric flow rate (m$^3$/s), $\dot{s}_k$ is the molar production rate (mol/m$^2$-s) of species $k$ due to surface reactions, $\dot{\omega}_k$ is the molar production rate (mol/m$^3$-s) of species $k$ due to gas phase reactions, and $A_s$ is the surface are in $m^2$ The subscripts *in* and *out* respectively refers to the inlet and outlet conditions. The outlet volumetric flow rate is related to the inlet volumetric flow rate according to 

```math
\dot{q}_{out} = \frac{q_{in} \bar{M_{in}} }{\bar{M}}
```

Here $\bar{M}$ refers the average molecular weights. The above equation for outlet volumetric flow rate is valid only when there is no mass accumulation within the reactor. 

# Executing the code 
## Surface chemistry
For solving a surface chemistry problem: On the Julia REPL 
```julia
julia>using StirredReactor
julia>plug("cstr.xml","lib/", surfchem=true)
```
## Gasphase chemistry
For solving a gasphase chemistry problem: On the Julia REPL 
```julia
julia>using StirredReactor
julia>plug("cstr.xml", "lib/", gaschem=true)
```

In the above calls, it is assumed that the input file *cstr.xml* is present in the working directory and *../lib/* is the path to the *lib* directory relative to the current working directory. The function can be called from any directory and in that case the first argument must point to the *cstr.xml* file relative to the current working directory. The output files will be generated in the directory where *plug.xml* is present. In general the function takes three optional keyword arguments *gaschem*, *surfchem*, and *sens*. *gaschem* must be true to simulate gasphase chemistry, *surfchem* must be true for surface chemistry, and *sens* must be true whenever sensitivity analysis is performed. 

## User defined chemistry
For solving the model with user defined chemistry: On the Julia REPL 
```julia
julia>using StirredReactor, ReactionCommons
julia>plug("cstr.xml", "lib/", udf)
```
*udf* is a function having the following signature
```
function udf(state::UserDefinedState)
```
where state is a structure defined as follows
```
struct UserDefinedState
    T::Float64
    p::Float64
    molefracs::Array{Float64,1}
    molwt::Array{Float64,1}
    species::Array{String,1}
    source::Array{Float64,1}
end
```
The program expects the species source terms in *source* mols/m3-s depending on whether you are solving surface chemistry problem or gasphase chemistry problem. The example call provided in the *runtests.jl* returns zero molar production and consumption rates. Within the code the source terms are multiplied with the molecular weight. The order ing of species source terms must be same as the order in wich the species appears in UserState.species.

## Input file for surface chemistry
The method takes *file\_name* as the argument. The file_name points to the input XML file. The structure of the XML input file is shown below.
```
<?xml version="1.0" encoding="ISO-8859-1"?>
<cstr>
    <gasphase>CH4 H2O H2 CO CO2 O2 N2</gasphase>
    <molefractions>CH4=0.25,H2O=0.25,N2=0.5</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
    <volume>1e-5</volume>
    <flow-rate>1.66e-6</flow-rate>
    <AsV>5e-3</AsV>
    <time>10</time>
    <surface_mech>data/ch4ni.xml</surface_mech>
</cstr>
```

## Input file for user defined chemistry
```
<?xml version="1.0" encoding="ISO-8859-1"?>
<cstr>
    <gasphase>CH4 H2O H2 CO CO2 O2 N2</gasphase>
    <molefractions>CH4=0.25,H2O=0.25,N2=0.5</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
    <volume>1e-5</volume>
    <flow-rate>1.66e-6</flow-rate>
    <AsV>9.8e6</AsV>
    <time>10</time>
</cstr>
```


## Input file for gasphase chemistry problems
```
<?xml version="1.0" encoding="ISO-8859-1"?>
<cstr>
    <molefractions>CH4=0.25,H2O=0.25,N2=0.5</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
    <volume>1e-5</volume>
    <flow-rate>1.66e-6</flow-rate>
    <time>10</time>
    <gas_mech>h2o2.dat</gas_mech>
</cstr>
```


The major difference between the input file for surface chemistry problem and gasphase chemistry problem is the <gasphase> tag of xml input. In the case of gasphase chemistry problem, the participating species are read from the mechanism input file, which is specified using the <gas_mech> tag

The meaning of different tags is specified below.

- <cstr> : The root XML tag for Cstr
- <gasphase> : list of gas-phase species. The species names must be separated by white spaces or tab (Required only in the case of surface mechanism)
- <molefractions> : inlet mole fraction of the gas-phase species. Instead of mole fractions, mass fractions may also be specified. In that case, the tag must be <massfractions>. You must ensure that the sum of mass or mole fractions specified is unity. There are no internal checks to ascertain this. 
- <T>: operating temperature in K
- <p>: initial pressure in Pa
- <volume> : reactor volume in m$^3$
- <flow-rate> : volumetric flow rate in m$^3$/s
- <AsV> : surface area per unit volume (1/m)
- <time> : integration time in s
- <surface_mech> : name of the surface reaction mechanism. 
- <gas_mech> : name of the gasphase reaction mechanism. 


## Input file download
The xml input file and the *lib* directory containig other required input files may be downloaded from [here](https://github.com/vinodjanardhanan/Cstr.jl/tree/main/test).


## Output
The code generates two output files in the same directory where the input file **`cstr.xml`** is present. 
The file **`gas_profile.dat`** contains the mole fraction of the gas phase species as a function of time.
The file **`surf_profile.dat`** contains the surface coverages of adsorbed species as a function of time. 
In addition to these two files, the code also generates terminal output, which shows integration progress.
The terminal output is nothing by the integration time. 

An example terminal output is shown below

```
julia> str("cstr_surf/cstr.xml", "lib/", surfchem=true)
0.0000e+00
6.6941e-16
3.5140e-15
8.7886e-15
1.4063e-14

..
..
..
9.7634e+00
9.9916e+00
1.0000e+01
:Success
```

A sample output of **`gas_profile.dat`** is shown below 
```
         t           T           p         rho         CH4         H2O          H2          CO         CO2          O2          N2
0.0000e+00      1.0732e+03      1.0000e+05      2.5240e-01      2.5000e-01      2.5000e-01      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      5.0000e-01
6.6941e-16      1.0732e+03      1.0000e+05      2.5240e-01      2.5000e-01      2.5000e-01      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      5.0000e-01
3.5140e-15      1.0732e+03      1.0000e+05      2.5240e-01      2.5000e-01      2.5000e-01      7.7918e-25      1.1878e-55      0.0000e+00      9.8832e-49      5.0000e-01
...
...
9.6792e+00      1.0732e+03      1.0000e+05      2.2750e-01      1.7600e-01      1.5516e-01      1.6890e-01      2.8364e-02      2.0917e-02      3.9562e-19      4.5066e-01
1.0000e+01      1.0732e+03      1.0000e+05      2.2745e-01      1.7584e-01      1.5499e-01      1.6922e-01      2.8459e-02      2.0929e-02      3.9756e-19      4.5056e-01
```

A sample output of **`surf_covg.dat`** is shown below 
```
         t           T        (NI)       H(NI)       O(NI)     CH4(NI)     H2O(NI)     CO2(NI)      CO(NI)      OH(NI)       C(NI)     HCO(NI)      CH(NI)     CH3(NI)     CH2(NI)
.0000e+00      1.0732e+03      6.0000e-01      0.0000e+00      0.0000e+00      0.0000e+00      4.0000e-01      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00
6.6941e-16      1.0732e+03      6.0000e-01      3.3170e-08      0.0000e+00      1.0070e-10      4.0000e-01      0.0000e+00      0.0000e+00      3.3170e-08      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      4.6583e-20
3.5140e-15      1.0732e+03      6.0001e-01      1.7415e-07      2.5760e-11      3.8046e-10      3.9999e-01      0.0000e+00      2.5729e-44      1.7410e-07      1.5650e-24      1.0506e-42      7.6562e-22      1.0501e-14      2.4240e-18
...
...
...
9.6792e+00      1.0732e+03      6.2128e-01      1.7478e-01      3.9515e-03      8.4703e-10      2.7792e-04      3.7446e-06      1.9966e-01      3.0604e-05      2.3685e-05      6.3767e-12      2.0463e-11      2.6785e-10      1.1690e-10
1.0000e+01      1.0732e+03      6.2095e-01      1.7485e-01      3.9375e-03      8.4584e-10      2.7747e-04      3.7438e-06      1.9992e-01      3.0525e-05      2.3713e-05      6.4007e-12      2.0464e-11      2.6771e-10      1.1682e-10
```
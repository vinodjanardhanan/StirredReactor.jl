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
\rho V \frac{dY_k}{dt} = \dot{q}_{in} \rho_{in} Y_{k,in} - \dot{q}_{out} \rho Y_k + \dot{s_k} M_k A_s + \dot{\omega_k}M_k V
```
Here $Y_k$ is the mass fraction of species $k$, $\rho$ is the density (kg/m$^3$) of the mixture, $V$ is the reactor volume (m$^3$), $M_k$ is the molecular weight of species $k$ (Kg/mol), $\dot{q}$ is the volumetric flow rate (m$^3$/s), $\dot{s}_k$ is the molar production rate (mol/m$^2$-s) of species $k$ due to surface reactions, $\dot{\omega}_k$ is the molar production rate (mol/m$^3$-s) of species $k$ due to gas phase reactions, and $A_s$ is the surface are in $m^2$ The subscripts *in* and *out* respectively refers to the inlet and outlet conditions. The outlet volumetric flow rate is related to the inlet volumetric flow rate according to 

```math
\dot{q}_{out} = \frac{q_{in} \bar{M_in} }{\bar{M}}
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
    <As>5e-3</As>
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
    <As>5e-3</As>
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
    <As>5e-3</As>
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
- <As> : surface area available for surface reactions m$^2$
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
0.0000e+00  1.0732e+03  1.0000e+05  2.5240e-01  2.5000e-01  2.5000e-01  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  5.0000e-01
5.4134e-12  1.0732e+03  1.0000e+05  2.5240e-01  2.5000e-01  2.5001e-01  1.4392e-15  1.4653e-29  1.9600e-41  1.8072e-32  4.9999e-01
7.3599e-12  1.0732e+03  1.0000e+05  2.5240e-01  2.5000e-01  2.5001e-01  3.6770e-15  2.9837e-29  3.0298e-40  5.6564e-32  4.9999e-01
...
...
9.9232e+00  1.0732e+03  1.0000e+05  2.2739e-01  1.7569e-01  1.5482e-01  1.6954e-01  2.8553e-02  2.0940e-02  3.9303e-19  4.5046e-01
1.0000e+01  1.0732e+03  1.0000e+05  2.2738e-01  1.7566e-01  1.5479e-01  1.6961e-01  2.8572e-02  2.0942e-02  4.0255e-19  4.5044e-01
```

A sample output of **`surf_covg.dat`** is shown below 
```
         t           T        (NI)       H(NI)       O(NI)     CH4(NI)     H2O(NI)     CO2(NI)      CO(NI)      OH(NI)       C(NI)     HCO(NI)      CH(NI)     CH3(NI)     CH2(NI)
0.0000e+00  1.0732e+03  6.0000e-01  0.0000e+00  0.0000e+00  0.0000e+00  4.0000e-01  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00  0.0000e+00
5.4134e-12  1.0732e+03  6.0847e-01  3.0035e-04  3.3271e-05  1.1713e-09  3.9096e-01  2.8741e-33  6.1498e-21  2.3380e-04  1.1015e-12  2.3785e-23  8.9706e-13  4.7230e-11  8.9960e-12
7.3599e-12  1.0732e+03  6.1146e-01  4.2230e-04  5.9785e-05  1.1839e-09  3.8776e-01  6.0780e-32  1.8970e-20  3.0273e-04  2.9147e-12  8.4330e-23  1.7009e-12  6.0154e-11  1.5013e-11
...
...
...
9.9232e+00  1.0732e+03  6.2063e-01  1.7493e-01  3.9236e-03  8.4466e-10  2.7702e-04  3.7430e-06  2.0019e-01  3.0447e-05  2.3741e-05  6.4247e-12  2.0465e-11  2.6758e-10  1.1674e-10
1.0000e+01  1.0732e+03  6.2056e-01  1.7494e-01  3.9209e-03  8.4442e-10  2.7693e-04  3.7428e-06  2.0024e-01  3.0431e-05  2.3746e-05  6.4295e-12  2.0465e-11  2.6755e-10  1.1672e-10
```
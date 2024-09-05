var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = StirredReactor","category":"page"},{"location":"#StirredReactor","page":"Home","title":"StirredReactor","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The Cstr model simulates a continuous stirred tank reactor. This type of reactor model may be used to simulate a packed bed reactor, when the length to diameter ratio is small, or when a packed bed reactor is operating under differential conditions. The code integrates the following governing equation ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Documentation for StirredReactor.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the package, use the following commands in the julia REPL","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg\njulia> Pkg.add(\"StirredReactor\")","category":"page"},{"location":"#General-interfaces","page":"Home","title":"General interfaces","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [StirredReactor]","category":"page"},{"location":"#StirredReactor.cstr-NTuple{6, Any}","page":"Home","title":"StirredReactor.cstr","text":"A function for call from other packages, mainly intended for reactor network modeling\n\nUsage\n\ncstr(inletcomp, T, p, time; Asv=1.0, chem, thermoobj, md)\n\ninlet_comp : A dictionary of species and its mole fractions at the reactor inlet \nT : operating temperature (K)\np : operating pressure (Pa) \nq : volumetric flow rate (m3/s)\ntime : integration time (s)   \nAsV : Surface area to volume ratio (important in the case of surface reactions. 1 in the case of gasphase chemistry)\nthermo_obj : Species thermo object (Please refer IdealGas package documentation)\nmd : MechanismDefinition (Please refer to ReactionCommons documentation)\n\n\n\n\n\n","category":"method"},{"location":"#StirredReactor.cstr-Tuple{AbstractString, AbstractString, Function}","page":"Home","title":"StirredReactor.cstr","text":"This is the calling function for executing the cstr reactor with userdefined rate calculation\n\nUsage\n\ncstr(inputfile::AbstractString, libdir::AbstractString; sens= false, surfchem=false, gaschem=false)\n\ninput_file: the xml input file for batch reactor\nlib_dir: the direcrtory in which the data files are present. It must be the relative path\nuser_defined: A function which calculates the species source terms, to be supplied by the user\n\n\n\n\n\n","category":"method"},{"location":"#StirredReactor.cstr-Tuple{AbstractString, AbstractString}","page":"Home","title":"StirredReactor.cstr","text":"This is the calling function for executing the cstr reactor with chemistry input files \n\nUsage\n\ncstr(inputfile::AbstractString, libdir::AbstractString; sens= false, surfchem=false, gaschem=false)\n\ninput_file: the xml input file for batch reactor\nlib_dir: the direcrtory in which the data files are present. It must be the relative path\nsens: sensitivity calculation \nsurfchem : surface chemistry. false means no surface chemistry calculation and \n\n\n\n\n\n","category":"method"},{"location":"#Governing-equations","page":"Home","title":"Governing equations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"V fracdC_kdt = dotq_in C_kin - dotq_out C_k + dots_k A_s + dotomega_k V","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here C_k is the concentration of species k, V is the reactor volume (m^3), dotq is the volumetric flow rate (m^3/s), dots_k is the molar production rate (mol/m^2-s) of species k due to surface reactions, dotomega_k is the molar production rate (mol/m^3-s) of species k due to gas phase reactions, and A_s is the surface are in m^2 The subscripts in and out respectively refers to the inlet and outlet conditions. The outlet volumetric flow rate is related to the inlet volumetric flow rate according to ","category":"page"},{"location":"","page":"Home","title":"Home","text":"dotq_out = fracq_in barM_in barM","category":"page"},{"location":"","page":"Home","title":"Home","text":"Here barM refers the average molecular weights. The above equation for outlet volumetric flow rate is valid only when there is no mass accumulation within the reactor. ","category":"page"},{"location":"#Executing-the-code","page":"Home","title":"Executing the code","text":"","category":"section"},{"location":"#Surface-chemistry","page":"Home","title":"Surface chemistry","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For solving a surface chemistry problem: On the Julia REPL ","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia>using StirredReactor\njulia>plug(\"cstr.xml\",\"lib/\", surfchem=true)","category":"page"},{"location":"#Gasphase-chemistry","page":"Home","title":"Gasphase chemistry","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For solving a gasphase chemistry problem: On the Julia REPL ","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia>using StirredReactor\njulia>plug(\"cstr.xml\", \"lib/\", gaschem=true)","category":"page"},{"location":"","page":"Home","title":"Home","text":"In the above calls, it is assumed that the input file cstr.xml is present in the working directory and ../lib/ is the path to the lib directory relative to the current working directory. The function can be called from any directory and in that case the first argument must point to the cstr.xml file relative to the current working directory. The output files will be generated in the directory where plug.xml is present. In general the function takes three optional keyword arguments gaschem, surfchem, and sens. gaschem must be true to simulate gasphase chemistry, surfchem must be true for surface chemistry, and sens must be true whenever sensitivity analysis is performed. ","category":"page"},{"location":"#User-defined-chemistry","page":"Home","title":"User defined chemistry","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For solving the model with user defined chemistry: On the Julia REPL ","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia>using StirredReactor, ReactionCommons\njulia>plug(\"cstr.xml\", \"lib/\", udf)","category":"page"},{"location":"","page":"Home","title":"Home","text":"udf is a function having the following signature","category":"page"},{"location":"","page":"Home","title":"Home","text":"function udf(state::UserDefinedState)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where state is a structure defined as follows","category":"page"},{"location":"","page":"Home","title":"Home","text":"struct UserDefinedState\n    T::Float64\n    p::Float64\n    molefracs::Array{Float64,1}\n    molwt::Array{Float64,1}\n    species::Array{String,1}\n    source::Array{Float64,1}\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"The program expects the species source terms in source mols/m3-s depending on whether you are solving surface chemistry problem or gasphase chemistry problem. The example call provided in the runtests.jl returns zero molar production and consumption rates. Within the code the source terms are multiplied with the molecular weight. The order ing of species source terms must be same as the order in wich the species appears in UserState.species.","category":"page"},{"location":"#Input-file-for-surface-chemistry","page":"Home","title":"Input file for surface chemistry","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The method takes file_name as the argument. The file_name points to the input XML file. The structure of the XML input file is shown below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<cstr>\n    <gasphase>CH4 H2O H2 CO CO2 O2 N2</gasphase>\n    <molefractions>CH4=0.25,H2O=0.25,N2=0.5</molefractions>\n    <T>1073.15</T>\n    <p>1e5</p>\n    <volume>1e-5</volume>\n    <flow-rate>1.66e-6</flow-rate>\n    <AsV>5e-3</AsV>\n    <time>10</time>\n    <surface_mech>data/ch4ni.xml</surface_mech>\n</cstr>","category":"page"},{"location":"#Input-file-for-user-defined-chemistry","page":"Home","title":"Input file for user defined chemistry","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<cstr>\n    <gasphase>CH4 H2O H2 CO CO2 O2 N2</gasphase>\n    <molefractions>CH4=0.25,H2O=0.25,N2=0.5</molefractions>\n    <T>1073.15</T>\n    <p>1e5</p>\n    <volume>1e-5</volume>\n    <flow-rate>1.66e-6</flow-rate>\n    <AsV>9.8e6</AsV>\n    <time>10</time>\n</cstr>","category":"page"},{"location":"#Input-file-for-gasphase-chemistry-problems","page":"Home","title":"Input file for gasphase chemistry problems","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<cstr>\n    <molefractions>CH4=0.25,H2O=0.25,N2=0.5</molefractions>\n    <T>1073.15</T>\n    <p>1e5</p>\n    <volume>1e-5</volume>\n    <flow-rate>1.66e-6</flow-rate>\n    <time>10</time>\n    <gas_mech>h2o2.dat</gas_mech>\n</cstr>","category":"page"},{"location":"","page":"Home","title":"Home","text":"The major difference between the input file for surface chemistry problem and gasphase chemistry problem is the <gasphase> tag of xml input. In the case of gasphase chemistry problem, the participating species are read from the mechanism input file, which is specified using the <gas_mech> tag","category":"page"},{"location":"","page":"Home","title":"Home","text":"The meaning of different tags is specified below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<cstr> : The root XML tag for Cstr\n<gasphase> : list of gas-phase species. The species names must be separated by white spaces or tab (Required only in the case of surface mechanism)\n<molefractions> : inlet mole fraction of the gas-phase species. Instead of mole fractions, mass fractions may also be specified. In that case, the tag must be <massfractions>. You must ensure that the sum of mass or mole fractions specified is unity. There are no internal checks to ascertain this. \n<T>: operating temperature in K\n<p>: initial pressure in Pa\n<volume> : reactor volume in m^3\n<flow-rate> : volumetric flow rate in m^3/s\n<AsV> : surface area per unit volume (1/m)\n<time> : integration time in s\n<surface_mech> : name of the surface reaction mechanism. \n<gas_mech> : name of the gasphase reaction mechanism. ","category":"page"},{"location":"#Input-file-download","page":"Home","title":"Input file download","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The xml input file and the lib directory containig other required input files may be downloaded from here.","category":"page"},{"location":"#Output","page":"Home","title":"Output","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The code generates four output files in the same directory where the input file cstr.xml is present.  The file gas_profile.dat, and gas_profile.csv contains the mole fraction of the gas phase species as a function of time. One is a tab separated file and other one is a comma separated one. The file surf_profile.dat. surf_profile.csv contains the surface coverages of adsorbed species as a function of time. One is a tab separated file and other one is a comma separated one. In addition to these files, the code also generates terminal output, which shows integration progress. The terminal output is nothing by the integration time. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"An example terminal output is shown below","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> str(\"cstr_surf/cstr.xml\", \"lib/\", surfchem=true)\n0.0000e+00\n6.6941e-16\n3.5140e-15\n8.7886e-15\n1.4063e-14\n\n..\n..\n..\n9.7634e+00\n9.9916e+00\n1.0000e+01\n:Success","category":"page"},{"location":"","page":"Home","title":"Home","text":"A sample output of gas_profile.dat is shown below ","category":"page"},{"location":"","page":"Home","title":"Home","text":"         t           T           p         rho         CH4         H2O          H2          CO         CO2          O2          N2\n0.0000e+00      1.0732e+03      1.0000e+05      2.5240e-01      2.5000e-01      2.5000e-01      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      5.0000e-01\n6.6941e-16      1.0732e+03      1.0000e+05      2.5240e-01      2.5000e-01      2.5000e-01      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      5.0000e-01\n3.5140e-15      1.0732e+03      1.0000e+05      2.5240e-01      2.5000e-01      2.5000e-01      7.7918e-25      1.1878e-55      0.0000e+00      9.8832e-49      5.0000e-01\n...\n...\n9.6792e+00      1.0732e+03      1.0000e+05      2.2750e-01      1.7600e-01      1.5516e-01      1.6890e-01      2.8364e-02      2.0917e-02      3.9562e-19      4.5066e-01\n1.0000e+01      1.0732e+03      1.0000e+05      2.2745e-01      1.7584e-01      1.5499e-01      1.6922e-01      2.8459e-02      2.0929e-02      3.9756e-19      4.5056e-01","category":"page"},{"location":"","page":"Home","title":"Home","text":"A sample output of surf_covg.dat is shown below ","category":"page"},{"location":"","page":"Home","title":"Home","text":"         t           T        (NI)       H(NI)       O(NI)     CH4(NI)     H2O(NI)     CO2(NI)      CO(NI)      OH(NI)       C(NI)     HCO(NI)      CH(NI)     CH3(NI)     CH2(NI)\n.0000e+00      1.0732e+03      6.0000e-01      0.0000e+00      0.0000e+00      0.0000e+00      4.0000e-01      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00\n6.6941e-16      1.0732e+03      6.0000e-01      3.3170e-08      0.0000e+00      1.0070e-10      4.0000e-01      0.0000e+00      0.0000e+00      3.3170e-08      0.0000e+00      0.0000e+00      0.0000e+00      0.0000e+00      4.6583e-20\n3.5140e-15      1.0732e+03      6.0001e-01      1.7415e-07      2.5760e-11      3.8046e-10      3.9999e-01      0.0000e+00      2.5729e-44      1.7410e-07      1.5650e-24      1.0506e-42      7.6562e-22      1.0501e-14      2.4240e-18\n...\n...\n...\n9.6792e+00      1.0732e+03      6.2128e-01      1.7478e-01      3.9515e-03      8.4703e-10      2.7792e-04      3.7446e-06      1.9966e-01      3.0604e-05      2.3685e-05      6.3767e-12      2.0463e-11      2.6785e-10      1.1690e-10\n1.0000e+01      1.0732e+03      6.2095e-01      1.7485e-01      3.9375e-03      8.4584e-10      2.7747e-04      3.7438e-06      1.9992e-01      3.0525e-05      2.3713e-05      6.4007e-12      2.0464e-11      2.6771e-10      1.1682e-10","category":"page"}]
}

module StirredReactor


using LightXML, Printf
using DifferentialEquations, Sundials
using SurfaceReactions, GasphaseReactions, ReactionCommons, IdealGas, RxnHelperUtils

export cstr

struct UsrMech <: MechanismDefinition
end

global o_streams 

struct ConstParams{T1}
    #ρ_in # inlet density
    q_in    #inlet flow rate 
    avg_molwt_in   #average molecular weight at the inlet
    AsV  #surface area
    V   #reactor volume
    #mass_fracs_in::Array{T1,1}  #inlet mass fractions
    c_in::Array{T1,1}
    T   # temperature 
    p   # pressure
    #ρq_yk::Array{T1,1} #inlet mass flow rate ρ x q x y_k
end


"""
This is the calling function for executing the cstr reactor with userdefined rate calculation
#   Usage
cstr(input_file::AbstractString, lib_dir::AbstractString; sens= false, surfchem=false, gaschem=false)
-   input_file: the xml input file for batch reactor
-   lib_dir: the direcrtory in which the data files are present. It must be the relative path
-   user_defined: A function which calculates the species source terms, to be supplied by the user
"""
function cstr(input_file::AbstractString, lib_dir::AbstractString, user_defined::Function; sens= false)    
    chem = Chemistry(false, false, true, user_defined)
    cstr(input_file, lib_dir, sens, chem)
end

"""
This is the calling function for executing the cstr reactor with chemistry input files 
#   Usage
cstr(input_file::AbstractString, lib_dir::AbstractString; sens= false, surfchem=false, gaschem=false)
-   input_file: the xml input file for batch reactor
-   lib_dir: the direcrtory in which the data files are present. It must be the relative path
-   sens: sensitivity calculation 
-   surfchem : surface chemistry. false means no surface chemistry calculation and 
"""
function cstr(input_file::AbstractString, lib_dir::AbstractString; sens= false, surfchem=false, gaschem=false)
    chem = Chemistry(surfchem, gaschem, false, f->())
    cstr(input_file, lib_dir, sens, chem)
end

"""
A function for call from other packages, mainly intended for reactor network modeling
# Usage
cstr(inlet_comp, T, p, time; Asv=1.0, chem, thermo_obj, md)
-   inlet_comp : A dictionary of species and its mole fractions at the reactor inlet 
-   T : operating temperature (K)
-   p : operating pressure (Pa) 
-   q : volumetric flow rate (m3/s)
-   time : integration time (s)   
-   AsV : Surface area to volume ratio (important in the case of surface reactions. 1 in the case of gasphase chemistry)
-   thermo_obj : Species thermo object (Please refer IdealGas package documentation)
-   md : MechanismDefinition (Please refer to ReactionCommons documentation)
"""
function  cstr(inlet_comp, T, p, q, V, time; AsV=1.0, chem, thermo_obj, md)

    function  get_mole_fracs(species, inlet_comp)
        mole_fracs = zeros(length(species))
        for k in eachindex(species)
            if in(species[k], collect(keys(inlet_comp)))
                mole_fracs[k] = inlet_comp[species[k]]
            end
        end
        mole_fracs
    end

    if chem.surfchem
        species = collect(keys(inlet_comp))        
        n_species = length(species) + length(md.sm.si.ini_covg)
    end
    if chem.gaschem
        species = md.gm.species
        n_species = length(species)
    end
    
    mole_fracs = get_mole_fracs(species, inlet_comp)    

    geom = (V, AsV)
    conditions = (T, p, q, mole_fracs)

    end_time, u = cstr_common( conditions, geom, chem, md, n_species, time, thermo_obj, false, interface_call = true)    
    mole_fracs = u[1:length(species)]/sum(u[1:length(species)])
    return end_time, Dict(species .=> mole_fracs)    
end

#=
General interface for the solution of CSTR
=#
function cstr(input_file::AbstractString, lib_dir::AbstractString, sens::Bool, chem::Chemistry)

    xmldoc = parse_file(input_file)    
    xmlroot = root(xmldoc)

    local gasphase = Array{String,1}

    """ In case the gasphase mechanism is present then the gasphase species are read from the 
     gasphase mechanism file and not from the xml input. So this must be the first action before any
     other input parameters are read 
     """
    if chem.gaschem
        mech_file = get_text_from_xml(xmlroot,"gas_mech")
        # mech_file = lib_dir*"/"*mech_file        
        mech_file = get_path(lib_dir, mech_file)
        gmd = compile_gaschemistry(mech_file)        
        gasphase = gmd.gm.species
    else
        #get the gasphase from xml input
        gasphase = get_collection_from_xml(xmlroot,"gasphase")
    end
    

    #create the thermo object
    thermo_file = "therm.dat"    
    thermo_file = get_path(lib_dir, thermo_file)
    thermo_obj = SurfaceReactions.IdealGas.create_thermo(gasphase,thermo_file)

    #get the molefractions from xml input
    mole_fracs = get_molefraction_from_xml(xmlroot,thermo_obj.molwt,gasphase)
    
    #get the operating temperature
    local T = get_value_from_xml(xmlroot,"T")

    #get the initial pressure 
    local p = get_value_from_xml(xmlroot,"p")

    #get the reactor volume 
    local V = get_value_from_xml(xmlroot,"volume")

    #get the flow rate. The flow rate is expected in m3/s
    local q = get_value_from_xml(xmlroot,"flow-rate")

    #get the integration time
    it = get_value_from_xml(xmlroot,"time")    
    

    #toal number of species
    n_species = length(gasphase) 


    #file output stream for saving the data
    g_stream = open(output_file(input_file, "gas_profile.dat"),"w")
    s_stream = open(output_file(input_file, "surf_covg.dat"),"w")
    global o_streams = (g_stream, s_stream)
    create_header(g_stream,["t","T","p","rho"],gasphase)    

    # As = 1.0
    #get the surface area
    AsV = get_value_from_xml(xmlroot,"AsV")
    # surface chemistry related 
    if chem.surfchem        
        #get the mechanism input file
        local mech_file = get_text_from_xml(xmlroot,"surface_mech")
        mech_file = lib_dir*"/"*mech_file
        smd = SurfaceReactions.compile_mech(mech_file,thermo_obj,gasphase)
        n_species += length(smd.sm.species)
        create_header(s_stream,"t","T", smd.sm.species)
    end

    geom = (V, AsV)
    conditions = (T, p, q, mole_fracs, gasphase)
    if chem.surfchem && !chem.gaschem
        mech_def = smd
    elseif chem.gaschem && !chem.surfchem
        mech_def = gmd
    elseif chem.gaschem && chem.surfchem
        mech_def = (gmd, smd)
    elseif chem.userchem 
        mech_def = UsrMech()
    end
    
    cstr_common( conditions, geom, chem, mech_def, n_species, it, thermo_obj, sens)
    

end

#=
This is a common interface that is called by all other methods 
=#
function cstr_common(conditions, geom, chem, mech_def, n_species,  time, thermo_obj, sens; interface_call = false)
    T = conditions[1]
    p = conditions[2]
    q = conditions[3]
    mole_fracs = conditions[4]
    V = geom[1]
    AsV = geom[2]
    
    # species source terms 
    source = zeros(n_species)
    # all concentrations
    conc_all = zeros(n_species)
    mass_fracs_in = zeros(length(mole_fracs))
    molefrac_to_massfrac!(mass_fracs_in, mole_fracs, thermo_obj.molwt)
    avg_molwt = average_molwt(mole_fracs,thermo_obj.molwt)
    #density at the inlet 
    ρ = density(mole_fracs,thermo_obj.molwt,T,p)
    #solution vector
    c_in = (p/RxnHelperUtils.R/T)*mole_fracs
    sol = deepcopy(c_in)

    
    if chem.gaschem
        gmd = mech_def
    end
    if chem.surfchem && !chem.gaschem
        smd = mech_def
    elseif chem.surfchem && chem.gaschem
        gmd = mech_def[1]
        smd = mech_def[2]
    end
 
    if chem.surfchem        
        #initial coverage as specified in the mechanism
        covg = smd.sm.si.ini_covg
        #append the coverage to the solution vector
        append!(sol,smd.sm.si.ini_covg)
        surf_conc = similar(smd.sm.si.ini_covg)
        rxn_rate = zeros(length(smd.sm.reactions))        
        #create state object 
        sr_state = SurfaceRxnState(T, p, mole_fracs, covg, surf_conc, rxn_rate, source, conc_all)
    end

    if chem.gaschem        
        g_rxn_rate = zeros(length(gmd.gm.reactions))
        Kp = similar(g_rxn_rate)
        source = similar(mole_fracs)
        g_all = similar(mole_fracs)
        conc = similar(mole_fracs)
        gr_state = GasphaseState(T, p, mole_fracs, conc , g_rxn_rate, source, g_all, Kp)
    end


    if chem.userchem && !chem.surfchem && !chem.gaschem
        ud_state = UserDefinedState(T, p, mole_fracs, thermo_obj.molwt, conditions[end], source)
    end

    
    #inlet conditions    
    t_span = (0,time)    
    #the following is the first term in the right hand side of gasphase species residual
    # ρq_mass_fracs = mass_fracs_in * (ρ*q)

    
    
    if chem.surfchem
        cp = ConstParams(q,avg_molwt,AsV,V,c_in, T, p)
        params = (sr_state, thermo_obj, smd, cp, chem)
    elseif chem.gaschem
        # As is not required in the case of gas chemistry
        cp = ConstParams(q,avg_molwt,1.0, V, c_in, T, p)
        params = (gr_state, thermo_obj, gmd, cp, chem)
    elseif chem.gaschem && chem.surfchem
        cp = ConstParams(q,avg_molwt,AsV,V,c_in, T, p)
        params = (gr_state, sr_state, thermo_obj, smd, gmd, cp, chem)
    else
        cp = ConstParams(q,avg_molwt,AsV, V, c_in, T, p)
        params = (ud_state, thermo_obj, UsrMech(), cp, chem)
    end
        
    prob = ODEProblem(residual!,sol,t_span,params)
    if sens == true
        return (params, prob, t_span)
    end
    # cb = FunctionCallingCallback(save_data)
    # soln = solve(prob,alg_hints=[:stiff] , reltol=1e-6, abstol = 1e-8, save_everystep=false, callback=cb)
    # soln = solve(prob, CVODE_BDF(), reltol=1e-6, abstol=1e-10, save_everystep=false,callback=cb);   
    # println("Solver Integration: ", soln.retcode, "\t", typeof(soln.retcode))

    if !interface_call
        cb = FunctionCallingCallback(save_data)
        soln = solve(prob, CVODE_BDF(), reltol=1e-10, abstol=1e-15, save_everystep=false,callback=cb);   
        # soln = solve(SteadyStateProblem(prob), DynamicSS(CVODE_BDF()), dt=1e-4, reltol=1e-10, abstol=1e-15, save_everystep=false,callback=cb)
        #close the files
        close(o_streams[1])
        close(o_streams[2])    
        return soln.retcode
    else
        soln = solve(prob, CVODE_BDF(), reltol=1e-10, abstol=1e-15, save_everystep=false);   
        soln = solve(SteadyStateProblem(prob), DynamicSS(CVODE_BDF()), dt=1e-10, reltol=1e-10, abstol=1e-15, save_everystep=false)
        return 0, soln.u        
    end

end




function residual!(du,u,p,t)
    # state, thermo_obj, md, cp, o_streams = p
    state = 1
    thermo_obj = 2
    md = 3
    cp = p[4]   # ConstParams
    chem = 5
    

    ng = length(p[state].mole_frac)
    
    conc = u[1:ng]
    #update the state with latest mole fractions
    p[state].mole_frac = conc/sum(conc)



    #outlet volumetric flow rate 
    q_out = cp.q_in*cp.avg_molwt_in/average_molwt(p[state].mole_frac,p[thermo_obj].molwt)


    #update the state with latest coverages
    if p[chem].surfchem
        ns = length(p[md].sm.species)
        p[state].covg = u[ng+1:ng+ns]
        #calculate the molar production rates        
        SurfaceReactions.calculate_molar_production_rates!(p[state],p[thermo_obj],p[md])
        # rgVec = (p[state].source[1:ng] .* p[thermo_obj].molwt)*cp.As/cp.V    
        rgVec = p[state].source[1:ng] * cp.AsV
        #surface species residual 
        du[ng+1:ng+ns] = (p[state].source[ng+1:ng+ns] .* p[md].sm.si.site_coordination)/(p[md].sm.si.density*1e4)            
    end
    
    if p[chem].gaschem
        GasphaseReactions.calculate_molar_production_rates!(p[state], p[md], p[thermo_obj])
        rgVec = p[state].source[1:ng] 
    end


    # call to get user defined rates 
    if p[chem].userchem
        p[chem].udf(p[state])
        rgVec = p[state].source[1:ng] * cp.AsV
    end
    # Gas species residual 
    du[1:ng] = ( (cp.q_in/cp.V)*cp.c_in - (q_out/cp.V) * u[1:ng]) + rgVec 

    
end
    

function save_data(u,t,integrator)
    state = integrator.p[1]
    thermo_obj = integrator.p[2]
    g_stream, s_stream = o_streams
    d = density(state.mole_frac,thermo_obj.molwt,state.T,state.p)
    write_to_file(g_stream,t,state.T,state.p,d,state.mole_frac)
    if integrator.p[5].surfchem
        write_to_file(s_stream,t,state.T,state.covg)    
    end
    @printf("%.4e\n",t)
end


end
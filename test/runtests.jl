using StirredReactor
using Test
using IdealGas, RxnHelperUtils, SurfaceReactions, GasphaseReactions, ReactionCommons

@testset "StirredReactor.jl" begin
    
    if Sys.isapple()  || Sys.islinux()
        lib_dir = "lib/"
    elseif Sys.iswindows()
        lib_dir = "lib\\"
    end

    # @testset "Testing surface chemistry" begin
    #     input_file = joinpath("cstr_surf", "cstr.xml")
    #     retcode = cstr(input_file, lib_dir, surfchem=true)
    #     @test retcode == Symbol("Success")        
    # end

    # @testset "Testing GRI mech chemistry" begin
    #     input_file = joinpath("cstr_ch4", "cstr.xml")
    #     retcode = cstr(input_file, lib_dir, gaschem=true)
    #     @test retcode == Symbol("Success")        
    # end

    # @testset "Testing h2o2 mech chemistry" begin
    #     input_file = joinpath("cstr_h2o2", "cstr.xml")
    #     retcode = cstr(input_file, lib_dir, gaschem=true)
    #     @test retcode == Symbol("Success")        
    # end

    # @testset "Testing surface chemistry with interface call " begin
    #     inlet_comp = Dict("CH4"=>0.25,"H2O"=>0.25, "H2"=>0.0, "CO"=>0.0, "CO2"=>0.0, "O2"=>0.0, "N2"=>0.5)
    #     T = 1073.15
    #     p = 1e5
    #     t = 10
    #     q = 1.66e-6
    #     V = 1e-5
    #     gasphase = collect(keys(inlet_comp))
    #     thermo_obj = IdealGas.create_thermo(gasphase, get_path(lib_dir,"therm.dat"))
    #     md = SurfaceReactions.compile_mech(get_path(lib_dir,"ch4ni.xml"),thermo_obj,gasphase)
    #     chem = Chemistry(true, false, false, f->())        
    #     retcodes = cstr(inlet_comp, T, p, q, V, t; AsV = 5e2, chem=chem, thermo_obj=thermo_obj, md=md)                      
    #     @test retcodes[1][end] == t
    # end


    # @testset "Testing gasphase chemistry with interface call " begin
        
    #     mech_file = get_path(lib_dir, "h2o2.dat")
    #     gmd = compile_gaschemistry(mech_file)        
    #     gasphase = gmd.gm.species
    #     thermo_obj = IdealGas.create_thermo(gasphase, get_path(lib_dir,"therm.dat"))        
    #     inlet_comp = Dict("O2" => 0.25,"N2" => 0.5, "H2" => 0.25)
        
    #     T = 1073.15
    #     p = 1e5
    #     t = 10
    #     q = 1.66e-6
    #     V = 1e-5        
    #     chem = Chemistry(false, true, false, f->())              
    #     retcodes = cstr(inlet_comp, T, p, q, V, t; AsV =1.0, chem=chem, thermo_obj=thermo_obj, md=gmd)                      
        
    #     @test retcodes[1][end] == t
    # end


    @testset "Testing user defined chemistry " begin        
        function udf(state)
            # println(state.species)
            # println(state.mole_frac)
            # println(state.T)
            # println(state.p)
            pp = zeros(length(state.mole_frac))
            ch4 = 1
            h2o = 2
            h2 = 3
            co = 4
            co2 = 5
            o2 = 6
            n2 = 7
            pp = (state.p/GasphaseReactions.R/state.T) * state.mole_frac
            #CH4 + 0.5 O2 -> 2H2 + CO
            k10 = 1e6
            E1 = 40e3            
            k20 = 1e15
            E2 = 30e3
            k30 = 1e15
            E3 = 30e3

            k1 = k10*exp(-E1/GasphaseReactions.R/state.T)
            k2 = k20 * exp(-E2/GasphaseReactions.R/state.T)
            k3 = k30 * exp(-E3/GasphaseReactions.R/state.T)
            if pp[o2] < 0
                r1 = 0.0
                r2 = 0.0
                r3 = 0.0
            else
                r1 = k1 * pp[ch4] * pp[o2]^0.5
                r2 = k2 * pp[co] * pp[o2]^0.5
                r3 = k3 * pp[h2] * pp[o2]^0.5
            end

            state.source[1:end] .= 0.0   
            state.source[ch4]  = -r1
            state.source[o2] = -0.5(r1+r2+r3)
            state.source[co] = r1 - r2
            state.source[h2] = 2r1 - r3
            state.source[co2] = r2
            state.source[h2o] = r3

            
        end
        input_file = joinpath("cstr_udf", "cstr.xml")
        retcode = cstr(input_file, lib_dir, udf)
        @test retcode == Symbol("Success")        
    end

end

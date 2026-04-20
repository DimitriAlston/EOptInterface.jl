# reference
# [1] Gujer, W., Henze, M., Mino, T., & Van Loosdrecht, M. (1999).
# Activated sludge model No. 3. Water Science and Technology, 39(1), 183-193.
#
# [2] Henze, M., Gujer, W., Mino, T., & Van Loosedrecht, M. (2006).
# Activated Sludge Models ASM1, ASM2, ASM2d and ASM3. IWA Publishing.

# This file follows the local ASM3_MPC case:
# influent + splitter recycle -> mixer1 -> clarifier recycle mixer -> five serial
# CSTRs -> splitter -> clarifier, with reactor2/reactor4 KLa as the manipulated inputs.

using ModelingToolkit, OrdinaryDiffEq, SciMLBase, Unitful, IfElse, NativeNaNMath
using ModelingToolkit: t_nounits as t, D_nounits as D

const biokinetic_model = "ASM3"

include("Initial_Condition.jl")
include("ASM3_Connector.jl")
include("ASM3_Constant_Inflow.jl")
include("ASM3_Mixer.jl")
include("ASM3_Splitter.jl")
include("ASM3_Clarifier.jl")
include("ASM3_CSTR_w_Aerator.jl")

const ASM3_REACTOR_VOL_M3 = 1000.0
const ASM3_TSPAN_INIT = (0.0, 1.0e-4)

@mtkmodel UConn_BNR_Series_CSTRs_Emu begin
    @components begin
        Inflow = InletStream(comp = Ini2vec, flow_rate = Ini2vecflow)
        reactor1 = CSTR(vol = ASM3_REACTOR_VOL_M3, switch = 0, Temp = Ini2vecT)
        reactor2 = CSTR(vol = ASM3_REACTOR_VOL_M3, switch = 1, Temp = Ini2vecT)
        reactor3 = CSTR(vol = ASM3_REACTOR_VOL_M3, switch = 0, Temp = Ini2vecT)
        reactor4 = CSTR(vol = ASM3_REACTOR_VOL_M3, switch = 1, Temp = Ini2vecT)
        reactor5 = CSTR(vol = ASM3_REACTOR_VOL_M3, switch = 0, Temp = Ini2vecT)
        splitter1 = Splitter_1in_2out(out_factor_1 = 0.5)
        clarifier = Clarifier(R_1 = 0.4, w_1 = 0.03)
        mixer1 = Mixer_2in_1out()
        mixer3 = Mixer_2in_1out()
    end

    @equations begin
        connect(Inflow.port, mixer1.In1)
        connect(splitter1.Out2, mixer1.In2)
        connect(mixer1.Out1, mixer3.In2)
        connect(clarifier.recycle_stream, mixer3.In1)
        connect(mixer3.Out1, reactor1.In)
        connect(reactor1.Out, reactor2.In)
        connect(reactor2.Out, reactor3.In)
        connect(reactor3.Out, reactor4.In)
        connect(reactor4.Out, reactor5.In)
        connect(reactor5.Out, splitter1.In)
        connect(splitter1.Out1, clarifier.inlet_stream)
    end
end

@mtkcompile sys = UConn_BNR_Series_CSTRs_Emu()

function initialize_asm3_series_with_clarifier(sys)
    op = Dict{Any, Any}()
    for i in 1:13
        op[sys.Inflow.comp[i]] = Ini2vec[i]
        for r in (sys.reactor1, sys.reactor2, sys.reactor3, sys.reactor4, sys.reactor5)
            op[r.x[i]] = Ini1vec[i]
        end
    end
    op[sys.Inflow.flow_rate] = Ini2vecflow

    guesses = Dict(
        sys.splitter1.In.flow_rate => Ini2vecflow,
        sys.mixer1.Out1.flow_rate => Ini2vecflow * 1.5,
        sys.mixer3.Out1.flow_rate => Ini2vecflow * 1.75,
        sys.clarifier.inlet_stream.flow_rate => Ini2vecflow,
    )

    prob = ODEProblem(sys, op, ASM3_TSPAN_INIT; guesses = guesses, warn_initialize_determined = false)
    return solve(
        prob,
        Rodas5P();
        initializealg = OrdinaryDiffEqNonlinearSolve.BrownFullBasicInit(),
        maxiters = 1.0e6,
    )
end

sol = initialize_asm3_series_with_clarifier(sys)

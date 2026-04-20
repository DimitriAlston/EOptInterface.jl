using ModelingToolkit

include("ASM3_Connector.jl")

@mtkmodel Clarifier begin
    @components begin
        inlet_stream = MaterialStream()
        outlet_stream = MaterialStream()
        waste_stream = MaterialStream()
        recycle_stream = MaterialStream()
    end

    @parameters begin
        R_1
        w_1
        b_X = 1 / (R_1 + w_1 * (1 - R_1))
    end

    @equations begin
        outlet_stream.flow_rate ~ inlet_stream.flow_rate * (1 - R_1) * (1 - w_1)
        outlet_stream.S_O ~ inlet_stream.S_O
        outlet_stream.S_I ~ inlet_stream.S_I
        outlet_stream.S_S ~ inlet_stream.S_S
        outlet_stream.S_NH ~ inlet_stream.S_NH
        outlet_stream.S_N2 ~ inlet_stream.S_N2
        outlet_stream.S_NO ~ inlet_stream.S_NO
        outlet_stream.S_ALK ~ inlet_stream.S_ALK
        outlet_stream.X_I ~ 0.0
        outlet_stream.X_S ~ 0.0
        outlet_stream.X_H ~ 0.0
        outlet_stream.X_STO ~ 0.0
        outlet_stream.X_A ~ 0.0
        outlet_stream.X_TS ~ 0.0

        recycle_stream.flow_rate ~ inlet_stream.flow_rate * R_1
        recycle_stream.S_O ~ inlet_stream.S_O
        recycle_stream.S_I ~ inlet_stream.S_I
        recycle_stream.S_S ~ inlet_stream.S_S
        recycle_stream.S_NH ~ inlet_stream.S_NH
        recycle_stream.S_N2 ~ inlet_stream.S_N2
        recycle_stream.S_NO ~ inlet_stream.S_NO
        recycle_stream.S_ALK ~ inlet_stream.S_ALK
        recycle_stream.X_I ~ b_X * inlet_stream.X_I
        recycle_stream.X_S ~ b_X * inlet_stream.X_S
        recycle_stream.X_H ~ b_X * inlet_stream.X_H
        recycle_stream.X_STO ~ b_X * inlet_stream.X_STO
        recycle_stream.X_A ~ b_X * inlet_stream.X_A
        recycle_stream.X_TS ~ b_X * inlet_stream.X_TS

        waste_stream.flow_rate ~ inlet_stream.flow_rate * w_1 * (1 - R_1)
        waste_stream.S_O ~ inlet_stream.S_O
        waste_stream.S_I ~ inlet_stream.S_I
        waste_stream.S_S ~ inlet_stream.S_S
        waste_stream.S_NH ~ inlet_stream.S_NH
        waste_stream.S_N2 ~ inlet_stream.S_N2
        waste_stream.S_NO ~ inlet_stream.S_NO
        waste_stream.S_ALK ~ inlet_stream.S_ALK
        waste_stream.X_I ~ b_X * inlet_stream.X_I
        waste_stream.X_S ~ b_X * inlet_stream.X_S
        waste_stream.X_H ~ b_X * inlet_stream.X_H
        waste_stream.X_STO ~ b_X * inlet_stream.X_STO
        waste_stream.X_A ~ b_X * inlet_stream.X_A
        waste_stream.X_TS ~ b_X * inlet_stream.X_TS
    end
end

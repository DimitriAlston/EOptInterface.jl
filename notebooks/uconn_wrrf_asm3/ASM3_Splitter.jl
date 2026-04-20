using ModelingToolkit

include("ASM3_Connector.jl")

@mtkmodel Splitter_1in_2out begin
    @components begin
        In = MaterialStream()
        Out1 = MaterialStream()
        Out2 = MaterialStream()
    end

    @parameters begin
        out_factor_1
    end

    @equations begin
        Out1.S_O   ~ In.S_O
        Out2.S_O   ~ In.S_O
        Out1.S_I   ~ In.S_I
        Out2.S_I   ~ In.S_I
        Out1.S_S   ~ In.S_S
        Out2.S_S   ~ In.S_S
        Out1.S_NH  ~ In.S_NH
        Out2.S_NH  ~ In.S_NH
        Out1.S_N2  ~ In.S_N2
        Out2.S_N2  ~ In.S_N2
        Out1.S_NO  ~ In.S_NO
        Out2.S_NO  ~ In.S_NO
        Out1.S_ALK ~ In.S_ALK
        Out2.S_ALK ~ In.S_ALK
        Out1.X_I   ~ In.X_I
        Out2.X_I   ~ In.X_I
        Out1.X_S   ~ In.X_S
        Out2.X_S   ~ In.X_S
        Out1.X_H   ~ In.X_H
        Out2.X_H   ~ In.X_H
        Out1.X_STO ~ In.X_STO
        Out2.X_STO ~ In.X_STO
        Out1.X_A   ~ In.X_A
        Out2.X_A   ~ In.X_A
        Out1.X_TS  ~ In.X_TS
        Out2.X_TS  ~ In.X_TS

        Out1.flow_rate ~ In.flow_rate * out_factor_1
        Out2.flow_rate ~ In.flow_rate - Out1.flow_rate
    end
end

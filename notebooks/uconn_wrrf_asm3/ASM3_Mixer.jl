using ModelingToolkit

include("ASM3_Connector.jl")

@mtkmodel Mixer_2in_1out begin
    @components begin
        In1 = MaterialStream()
        In2 = MaterialStream()
        Out1 = MaterialStream()
    end

    @equations begin
        Out1.flow_rate ~ In1.flow_rate + In2.flow_rate
        Out1.S_O   ~ (In1.S_O   * In1.flow_rate + In2.S_O   * In2.flow_rate) / Out1.flow_rate
        Out1.S_I   ~ (In1.S_I   * In1.flow_rate + In2.S_I   * In2.flow_rate) / Out1.flow_rate
        Out1.S_S   ~ (In1.S_S   * In1.flow_rate + In2.S_S   * In2.flow_rate) / Out1.flow_rate
        Out1.S_NH  ~ (In1.S_NH  * In1.flow_rate + In2.S_NH  * In2.flow_rate) / Out1.flow_rate
        Out1.S_N2  ~ (In1.S_N2  * In1.flow_rate + In2.S_N2  * In2.flow_rate) / Out1.flow_rate
        Out1.S_NO  ~ (In1.S_NO  * In1.flow_rate + In2.S_NO  * In2.flow_rate) / Out1.flow_rate
        Out1.S_ALK ~ (In1.S_ALK * In1.flow_rate + In2.S_ALK * In2.flow_rate) / Out1.flow_rate
        Out1.X_I   ~ (In1.X_I   * In1.flow_rate + In2.X_I   * In2.flow_rate) / Out1.flow_rate
        Out1.X_S   ~ (In1.X_S   * In1.flow_rate + In2.X_S   * In2.flow_rate) / Out1.flow_rate
        Out1.X_H   ~ (In1.X_H   * In1.flow_rate + In2.X_H   * In2.flow_rate) / Out1.flow_rate
        Out1.X_STO ~ (In1.X_STO * In1.flow_rate + In2.X_STO * In2.flow_rate) / Out1.flow_rate
        Out1.X_A   ~ (In1.X_A   * In1.flow_rate + In2.X_A   * In2.flow_rate) / Out1.flow_rate
        Out1.X_TS  ~ (In1.X_TS  * In1.flow_rate + In2.X_TS  * In2.flow_rate) / Out1.flow_rate
    end
end

using ModelingToolkit

include("ASM3_Connector.jl")

@mtkmodel InletStream begin
    @components begin
        port = MaterialStream()
    end

    @parameters begin
        comp[1:13]
        flow_rate
    end

    @equations begin
        port.S_O ~ comp[1]
        port.S_I ~ comp[2]
        port.S_S ~ comp[3]
        port.S_NH ~ comp[4]
        port.S_N2 ~ comp[5]
        port.S_NO ~ comp[6]
        port.S_ALK ~ comp[7]
        port.X_I ~ comp[8]
        port.X_S ~ comp[9]
        port.X_H ~ comp[10]
        port.X_STO ~ comp[11]
        port.X_A ~ comp[12]
        port.X_TS ~ comp[13]
        port.flow_rate ~ flow_rate
    end
end

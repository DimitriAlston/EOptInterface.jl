using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

include("ASM3.jl")
include("ASM3_Connector.jl")

@mtkmodel CSTR begin
    @extend ASM3()

    @components begin
        In = MaterialStream()
        Out = MaterialStream()
    end

    @parameters begin
        vol
        KLa(t) = 240
        switch
        Temp
    end

    @variables begin
        SO_sat_temp(t)
        KLa_temp(t)
    end

    @equations begin
        T_ref ~ Temp

        D(x[1])  ~ 1.0 / vol * (In.flow_rate * (In.S_O - x[1])) + reac1  + switch * KLa_temp * (SO_sat_temp - x[1])
        D(x[2])  ~ 1.0 / vol * (In.flow_rate * (In.S_I - x[2])) + reac2
        D(x[3])  ~ 1.0 / vol * (In.flow_rate * (In.S_S - x[3])) + reac3
        D(x[4])  ~ 1.0 / vol * (In.flow_rate * (In.S_NH - x[4])) + reac4
        D(x[5])  ~ 1.0 / vol * (In.flow_rate * (In.S_N2 - x[5])) + reac5
        D(x[6])  ~ 1.0 / vol * (In.flow_rate * (In.S_NO - x[6])) + reac6
        D(x[7])  ~ 1.0 / vol * (In.flow_rate * (In.S_ALK - x[7])) + reac7
        D(x[8])  ~ 1.0 / vol * (In.flow_rate * (In.X_I - x[8])) + reac8
        D(x[9])  ~ 1.0 / vol * (In.flow_rate * (In.X_S - x[9])) + reac9
        D(x[10]) ~ 1.0 / vol * (In.flow_rate * (In.X_H - x[10])) + reac10
        D(x[11]) ~ 1.0 / vol * (In.flow_rate * (In.X_STO - x[11])) + reac11
        D(x[12]) ~ 1.0 / vol * (In.flow_rate * (In.X_A - x[12])) + reac12
        D(x[13]) ~ 1.0 / vol * (In.flow_rate * (In.X_TS - x[13])) + reac13

        Out.S_O ~ x[1]
        Out.S_I ~ x[2]
        Out.S_S ~ x[3]
        Out.S_NH ~ x[4]
        Out.S_N2 ~ x[5]
        Out.S_NO ~ x[6]
        Out.S_ALK ~ x[7]
        Out.X_I ~ x[8]
        Out.X_S ~ x[9]
        Out.X_H ~ x[10]
        Out.X_STO ~ x[11]
        Out.X_A ~ x[12]
        Out.X_TS ~ x[13]
        Out.flow_rate ~ In.flow_rate

        SO_sat_temp ~ 0.9997743214 * 8.0 / 10.5 * (56.12 * 6791.5 * exp(-66.7354 + 87.4755 / ((Temp + 273.15) / 100.0) + 24.4526 * log((Temp + 273.15) / 100.0)))
        KLa_temp ~ KLa * 1.024^(Temp - 15.0)
    end
end

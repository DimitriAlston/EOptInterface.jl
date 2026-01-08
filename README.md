![CSTRfig](https://github.com/user-attachments/assets/d4fe014c-9ca3-4c65-bd47-f3ae5bb01549)

# EOptInterface.jl

EOptInterface.jl is an abstraction layer for automatically formulating JuMP mathematical programming models from ModelingToolkit equation-oriented/acausal models.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PSORLab.github.io/EOptInterface.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PSORLab.github.io/EOptInterface.jl/dev/)
[![Build Status](https://github.com/PSORLab/EOptInterface.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/PSORLab/EOptInterface.jl/actions/workflows/CI.yml?query=branch%3Amaster)

## Feature Summary

```julia
decision_vars(::System)
```
Displays the optimization problem decision variables.

```julia
register_nlsystem(::Model, ::System, obj::Num, ineqs::Vector{Num})
```
Registers algebraic JuMP constraints and objective from ModelingToolkit algebraic `System`s built using `@mtkbuild`.

```julia
full_solutions(::Model, ::System)
```
Returns a dictionary of optimal solution values for all eliminated variables from ModelingToolkit's structural simplification step.

```julia
register_odesystem(::Model, ::System, tspan::Tuple{Number,Number}, tstep::Number, solver::String)
```
Registers algebraic JuMP constraints from ModelingToolkit ODE `System`s built using `@mtkbuild`. Available integration schemes: `"EE", "IE"`




## Example Usage: Algebraic Models

An optimal reactor-separator-recycle process design problem is used to demonstrate the use of `register_nlsystem` to formulate and solve a reduced-space `Model` using the deterministic global optimizer `EAGO`.

![Uploadi<svg width="3465" height="1815" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" xml:space="preserve" overflow="hidden"><defs><linearGradient x1="2582.09" y1="1659.38" x2="2582.66" y2="1790.36" gradientUnits="userSpaceOnUse" spreadMethod="reflect" id="fill0"><stop offset="0" stop-color="#7E7E7E"/><stop offset="0.5" stop-color="#B6B6B6"/><stop offset="1" stop-color="#D9D9D9"/></linearGradient><linearGradient x1="2733.82" y1="1724.22" x2="2734.39" y2="1855.19" gradientUnits="userSpaceOnUse" spreadMethod="reflect" id="fill1"><stop offset="0" stop-color="#7E7E7E"/><stop offset="0.5" stop-color="#B6B6B6"/><stop offset="1" stop-color="#D9D9D9"/></linearGradient></defs><g transform="translate(-1171 -815)"><rect x="3981" y="1145" width="330" height="1155" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#E2F0D9"/><text font-family="Arial,Arial_MSFontService,sans-serif" font-weight="400" font-size="110" transform="matrix(1 0 0 1 4026.74 1764)">sep2</text><path d="M0-6.875 165-6.875C168.797-6.875 171.875-3.79696 171.875 0L171.875 316.313 158.125 316.313 158.125 0 165 6.875 0 6.875ZM195.001 278.527 165 329.958 134.999 278.527C133.085 275.248 134.193 271.038 137.473 269.125 140.753 267.211 144.962 268.319 146.876 271.599L170.938 312.849 159.061 312.849 183.124 271.599C185.037 268.319 189.247 267.211 192.526 269.125 195.806 271.038 196.914 275.248 195.001 278.527Z" transform="matrix(1.83697e-16 1 1 -1.83697e-16 4146 2300)"/><path d="M4139.12 1145.52 4139.12 977C4139.12 973.203 4142.2 970.125 4146 970.125L4462.31 970.125 4462.31 983.875 4146 983.875 4152.88 977 4152.88 1145.52ZM4424.53 946.999 4475.96 977 4424.53 1007C4421.25 1008.91 4417.04 1007.81 4415.12 1004.53 4413.21 1001.25 4414.32 997.037 4417.6 995.124L4458.85 971.062 4458.85 982.939 4417.6 958.876C4414.32 956.963 4413.21 952.753 4415.12 949.473 4417.04 946.194 4421.25 945.086 4424.53 946.999Z"/><path d="M4097 1145 4146.5 1056 4196 1145Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#FF1A1A" fill-rule="evenodd"/><path d="M4196 2300 4146.5 2389 4097 2300Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#FF1A1A" fill-rule="evenodd"/><path d="M3881 1673 3970 1722.5 3881 1772Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#1A1AFF" fill-rule="evenodd"/><rect x="3320" y="1145" width="330" height="1155" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#E2F0D9"/><text font-family="Arial,Arial_MSFontService,sans-serif" font-weight="400" font-size="110" transform="matrix(1 0 0 1 3365.08 1764)">sep1</text><path d="M3534 2300 3484.5 2389 3435 2300Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#FF1A1A" fill-rule="evenodd"/><path d="M3218 1673 3307 1722.5 3218 1772Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#1A1AFF" fill-rule="evenodd"/><path d="M3435 1145 3484.5 1056 3534 1145Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#FF1A1A" fill-rule="evenodd"/><path d="M1667 1102 1931 1250.5 1667 1399Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#E2F0D9" fill-rule="evenodd"/><text font-family="Arial,Arial_MSFontService,sans-serif" font-weight="400" font-size="110" transform="matrix(1 0 0 1 1652.63 1083)">mixer</text><path d="M1931 1201 2020 1250.5 1931 1300Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#FF1A1A" fill-rule="evenodd"/><path d="M1566 1300 1655 1349.5 1566 1399Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#1A1AFF" fill-rule="evenodd"/><path d="M1566 1104 1655 1153.5 1566 1203Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#1A1AFF" fill-rule="evenodd"/><path d="M2020 1251 2258.45 1251 2258.45 1310.4" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="none" fill-rule="evenodd"/><path d="M0 0 72.1037 0 72.1037 414.221 135.889 414.221" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="none" fill-rule="evenodd" transform="matrix(1 0 0 -1 3082 2137.22)"/><rect x="2096" y="1409" width="891" height="891" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#E2F0D9"/><text font-family="Arial,Arial_MSFontService,sans-serif" font-weight="400" font-size="110" transform="matrix(1 0 0 1 2453.29 2094)">cstr</text><path d="M0 0 0.000360892 624.324" stroke="#000000" stroke-width="6.875" stroke-miterlimit="8" fill="none" fill-rule="evenodd" transform="matrix(-0.919567 -0.392934 -0.392934 0.919567 2903.56 1183.18)"/><path d="M2506.51 1692.45C2514.25 1674.34 2554.49 1674.16 2596.39 1692.06 2638.29 1709.97 2665.98 1739.17 2658.24 1757.29 2650.5 1775.41 2610.25 1775.58 2568.36 1757.68 2526.46 1739.77 2498.77 1710.57 2506.51 1692.45Z" stroke="#172C51" stroke-width="4.58333" stroke-miterlimit="8" fill="url(#fill0)" fill-rule="evenodd"/><path d="M2658.24 1757.29C2665.98 1739.17 2706.22 1739 2748.12 1756.9 2790.02 1774.8 2817.71 1804 2809.97 1822.12 2802.22 1840.24 2761.98 1840.41 2720.08 1822.51 2678.19 1804.61 2650.5 1775.41 2658.24 1757.29Z" stroke="#172C51" stroke-width="4.58333" stroke-miterlimit="8" fill="url(#fill1)" fill-rule="evenodd"/><path d="M36.7088 112.372 29.0004 108.28 21.8736 103.662 15.5921 98.678 10.1695 93.3634 5.62179 87.7549 1.96877 81.8901-0.766203 75.8103C-0.833481 75.6589-0.889471 75.5035-0.934264 75.3439L-2.55714 69.563C-2.60415 69.3955-2.6382 69.2258-2.65942 69.0531L-3.3782 63.2035C-3.39993 63.0266-3.40787 62.8489-3.40198 62.6708L-3.20798 56.7957C-3.20211 56.6179-3.18244 56.4407-3.14917 56.266L-2.03371 50.4085C-2.00131 50.2383-1.95611 50.0709-1.89848 49.9076L0.147133 44.1107C0.201998 43.9553 0.267649 43.8046 0.3442 43.6586L3.32864 37.9656 7.50087 32.0274 12.6542 26.3423 18.7824 20.9501 25.7187 16.0001 33.2513 11.6125 41.3052 7.7856 49.8096 4.51766 58.6949 1.80814 67.8925-0.342155 77.3343-1.93126 86.9525-2.9562 96.6793-3.4131 106.447-3.29737 116.187-2.60367 125.832-1.32603 135.315 0.542281 144.566 3.00871 153.517 6.08157 162.099 9.77002 170.059 13.9834 177.185 18.6008 183.467 23.5849 188.889 28.8994 193.437 34.508 197.09 40.3728 199.825 46.4526C199.892 46.6038 199.948 46.7597 199.993 46.9189L201.616 52.6999C201.663 52.8671 201.697 53.0374 201.718 53.2098L202.437 59.0594C202.459 59.2361 202.467 59.4141 202.461 59.5921L202.267 65.4672C202.261 65.6448 202.241 65.8223 202.208 65.9968L201.093 71.8544C201.06 72.0244 201.015 72.1921 200.957 72.3553L198.912 78.1521C198.857 78.3076 198.791 78.4583 198.715 78.6043L195.73 84.2973 191.558 90.2355 186.405 95.9206 180.277 101.313 170.692 107.894 165.774 110.464 162.618 104.356 166.812 102.218 175.749 96.1394 181.327 91.2861 185.95 86.2592 189.641 81.1053 192.626 75.4122 192.429 75.8644 194.474 70.0675 194.339 70.5683 195.454 64.7108 195.396 65.2404 195.59 59.3652 195.613 59.8979 194.895 54.0482 194.997 54.5581 193.374 48.7771 193.542 49.2435 191.241 43.9849 188.083 38.8213 184.065 33.7975 179.183 28.962 173.438 24.3647 166.835 20.0556 159.377 16.0835 151.279 12.5821 142.789 9.65027 133.981 7.28666 124.925 5.48885 115.695 4.25365 106.361 3.5771 96.9976 3.45452 87.6766 3.88057 78.4706 4.84919 69.4522 6.35357 60.6941 8.38603 52.2686 10.9379 44.2474 13.9992 36.7017 17.5589 29.7008 21.6045 23.3102 26.1235 17.7323 30.9768 13.1094 36.0037 9.41769 41.1576 6.43325 46.8507 6.63031 46.3985 4.5847 52.1954 4.71993 51.6946 3.60446 57.5521 3.66327 57.0225 3.46927 62.8977 3.44548 62.365 4.16426 68.2147 4.06198 67.7048 5.68486 73.4858 5.5168 73.0194 7.81841 78.278 10.9754 83.4416 14.9939 88.4654 19.8761 93.3009 25.6209 97.8982 32.2243 102.207 39.9328 106.3ZM195.816 112.47 158.135 110.542 178.361 78.6915C179.378 77.0889 181.503 76.6146 183.105 77.6323 184.708 78.65 185.182 80.7742 184.164 82.3769L167.098 109.253 164.372 103.977 196.167 105.604C198.063 105.701 199.521 107.316 199.424 109.212 199.327 111.108 197.712 112.567 195.816 112.47Z" transform="matrix(0.919567 0.392934 0.392934 -0.919567 2747.7 1294.61)"/><path d="M2308 1310 2258.5 1399 2209 1310Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#1A1AFF" fill-rule="evenodd"/><path d="M2993 2087 3082 2136.5 2993 2186Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#FF1A1A" fill-rule="evenodd"/><path d="M1336 1669C1336 1577.87 1409.87 1504 1501 1504 1592.13 1504 1666 1577.87 1666 1669 1666 1760.13 1592.13 1834 1501 1834 1409.87 1834 1336 1760.13 1336 1669Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#E2F0D9" fill-rule="evenodd"/><text font-family="Arial,Arial_MSFontService,sans-serif" font-weight="400" font-size="110" transform="matrix(1 0 0 1 1394.11 1939)">feed</text><path d="M1452 1503 1501.5 1413 1551 1503Z" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="#FF1A1A" fill-rule="evenodd"/><path d="M0 0 64.9166 0 64.9166 63.499" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="none" fill-rule="evenodd" transform="matrix(-1 1.22465e-16 1.22465e-16 1 1565.92 1350)"/><path d="M3485 2389.76 3485 2465.49 3815.37 2465.49 3815.37 1723 3881.19 1723" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="none" fill-rule="evenodd"/><path d="M3484.37 1056 3484.37 973.499 1502.1 973.499 1502.1 1153.17 1566 1153.17" stroke="#000000" stroke-width="13.75" stroke-miterlimit="8" fill="none" fill-rule="evenodd"/><rect x="1865.5" y="2107.5" width="99" height="99" stroke="#000000" stroke-width="4.58333" stroke-miterlimit="8" fill="#E2F0D9"/><text font-family="Arial,Arial_MSFontService,sans-serif" font-weight="400" font-size="69" transform="matrix(1 0 0 1 1491.03 2174)">component</text><rect x="1865.5" y="2404.5" width="99" height="99" stroke="#000000" stroke-width="4.58333" stroke-miterlimit="8" fill="#1A1AFF"/><rect x="1865.5" y="2251.5" width="99" height="99" stroke="#000000" stroke-width="4.58333" stroke-miterlimit="8" fill="#FF1A1A"/><rect x="0" y="0" width="660" height="462" stroke="#000000" stroke-width="6.875" stroke-miterlimit="8" fill="none" transform="matrix(-1 0 0 1 1996.5 2069.5)"/><text font-family="Arial,Arial_MSFontService,sans-serif" font-weight="400" font-size="69" transform="matrix(1 0 0 1 1380.7 2323)">connection.out</text><text font-family="Arial,Arial_MSFontService,sans-serif" font-weight="400" font-size="69" transform="matrix(1 0 0 1 1425.72 2474)">connection.in</text></g></svg>ng CSTRfig.svg…]()

<p align="center">
  *Figure 1: Illustration of the reactor-separator-recycle system.*
</p>

```julia
using ModelingToolkit, JuMP, EOptInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

# Formulate MTK algebraic System
@connector Stream begin
    @variables begin
        F(t),   [input=true]
        y_A(t), [input=true]
        y_B(t), [input=true]
        y_C(t), [input=true]
    end
    @parameters begin
        V_A = 8.937e-2
        V_B = 1.018e-1
        V_C = 1.13e-1
    end
end
@mtkmodel Influent begin
    @components begin
        out = Stream()
    end
    @parameters begin
        F       # Free design variable
        y_A = 1
        y_B = 0
        y_C = 0
    end
    @equations begin
        out.F ~ F
        out.y_A ~ y_A
        out.y_B ~ y_B
        out.y_C ~ y_C
    end
end
@mtkmodel Mixer begin
    @components begin
        in1 = Stream()
        in2 = Stream()
        out = Stream()
    end
    @equations begin
        out.F ~ in1.F + in2.F
        out.y_A ~ (in1.y_A*in1.F + in2.y_A*in2.F)/(in1.F + in2.F)
        out.y_B ~ (in1.y_B*in1.F + in2.y_B*in2.F)/(in1.F + in2.F)
        out.y_C ~ (in1.y_C*in1.F + in2.y_C*in2.F)/(in1.F + in2.F)
    end
end
@mtkmodel CSTR begin
    @components begin
        in = Stream()
        out = Stream()
    end
    @parameters begin
        V       # Free design variable
        k_1 = 0.4
        k_2 = 0.055
    end
    begin
        r_1 = k_1*out.y_A/(out.y_A*in.V_A + out.y_B*in.V_B + out.y_C*in.V_C)
        r_2 = k_2*out.y_B/(out.y_A*in.V_A + out.y_B*in.V_B + out.y_C*in.V_C)
    end
    @equations begin
        out.F ~ in.F
        out.y_A + out.y_B + out.y_C ~ 1
        out.y_B*out.F ~ in.y_B*in.F + (r_1 - r_2)*V
        out.y_C*out.F ~ in.y_C*in.F + r_2*V
    end
end
@mtkmodel Separator1 begin
    @components begin
        in = Stream()
        outV = Stream()
        outL = Stream()
    end
    @equations begin
        in.F ~ outV.F + outL.F
        in.y_B*in.F ~ outL.y_B*outL.F
        in.y_C*in.F ~ outL.y_C*outL.F
        
        outV.y_A + outV.y_B + outV.y_C ~ 1
        outV.y_C ~ 0
        outV.y_B ~ 0

        outL.y_A + outL.y_B + outL.y_C ~ 1
        outL.y_A ~ 0
    end
end
@mtkmodel Separator2 begin
    @components begin
        in = Stream()
        outV = Stream()
        outL = Stream()
    end
    @equations begin
        in.F ~ outV.F + outL.F
        in.y_B*in.F ~ outV.F

        outV.y_A + outV.y_B + outV.y_C ~ 1
        outV.y_A ~ 0
        outV.y_C ~ 0

        outL.y_A + outL.y_B + outL.y_C ~ 1
        outL.y_A ~ 0
        outL.y_B ~ 0
    end
end
@mtkmodel ReactorSeparatorRecycle begin
    @components begin
        influent = Influent()
        mixer = Mixer()
        cstr = CSTR()
        sep1 = Separator1()
        sep2 = Separator2()
    end
    @equations begin
        connect(influent.out, mixer.in1)
        connect(mixer.out, cstr.in)
        connect(cstr.out, sep1.in)
        connect(sep1.outV, mixer.in2)
        connect(sep1.outL, sep2.in)
    end
end

@mtkcompile s = ReactorSeparatorRecycle()

# Define symbolic expressions of constraints and objective
# Use syntax System.Component.Connector.Variable, System.Component.Component.Parameter, or System.Component.Parameter
exprF5 = s.sep2.outV.F
exprTau = s.cstr.V/(s.cstr.out.F*(s.cstr.out.y_A*s.cstr.in.V_A + s.cstr.out.y_B*s.cstr.in.V_B + s.cstr.out.y_C*s.cstr.in.V_C))
f_CSTR = (25764 + 8178*s.cstr.V)/2.5
s1cap = 132718 + s.cstr.out.F*(369*s.cstr.out.y_A - 1113.9*s.cstr.out.y_B)
s2cap = 25000 + s.sep1.outL.F*(6984.5*s.sep1.outL.y_B - 3869.53*s.sep1.outL.y_C^2)
s1op = s.cstr.out.F*(3+36.11*s.cstr.out.y_A + 7.71*s.cstr.out.y_B)*26.32e-3
s2op = s.sep1.outL.F*(26.21 + 29.45*s.sep1.outL.y_B)*26.32e-3;
f_Sep = (s1cap+s2cap)/2.5 + 0.52*(s1op+s2op)
g1 = 25 - exprF5
g2 = 475/3600 - exprTau
obj = f_CSTR + f_Sep

# Solve using EAGO
using EAGO
model = Model(EAGO.Optimizer)
decision_vars(s) # Displays: sep1₊in₊F(t), sep1₊in₊y_B(t), sep1₊in₊y_C(t), sep1₊outL₊y_C(t), influent₊F, cstr₊V
xL = zeros(6) # lower bound on x
xU = [100, 1, 1, 1, 100, 10] # upper bound on x
@variable(model, xL[i] <= x[i=1:6] <= xU[i]) # ̂x = (̂z,p), ̂z = (z(t),...), p = (influent₊F, cstr₊V)
register_nlsystem(model, s, obj, [g1, g2])
JuMP.optimize!(model)
JuMP.value.(x)

# Obtain observed variable solutions
full_solutions(model, s)
```

## Example Usage: ODE Models

A nonlinear kinetic parameter estimation problem is used to demonstrate the use of `register_odesystem` to formulate a `Model` from an ODE `System`, solved using `Ipopt`.

```julia
using ModelingToolkit, JuMP, EOptInterface
using ModelingToolkit: t_nounits as t, D_nounits as D

# Formulate MTK ODE System
@mtkmodel KineticParameterEstimation begin
    @parameters begin
        T = 273
        K_2 = 46*exp(6500/T-18)
        K_3 = 2*K_2
        k_1 = 53
        k_1s = k_1*1e-6
        k_5 = 1.2e-3
        c_O2 = 2e-3

        k_2f    # Free design variable
        k_3f    # Free design variable
        k_4     # Free design variable
    end
    @variables begin
        x_A(t) = 0.0
        x_B(t) = 0.0
        x_D(t) = 0.0
        x_Y(t) = 0.4
        x_Z(t) = 140.0
        I(t)
    end
    @equations begin
        D(x_A) ~ k_1*x_Z*x_Y - c_O2*(k_2f + k_3f)*x_A + k_2f/K_2*x_D + k_3f/K_3*x_B - k_5*x_A^2
        D(x_B) ~ c_O2*k_3f*x_A - (k_3f/K_3 + k_4)*x_B
        D(x_D) ~ c_O2*k_2f*x_A - k_2f/K_2*x_D
        D(x_Y) ~ -k_1s*x_Z*x_Y
        D(x_Z) ~ -k_1*x_Z*x_Y
        I ~ x_A + 2/21*x_B + 2/21*x_D
    end
end

@mtkcompile o = KineticParameterEstimation()

tspan = (0.0,2.0)
tstep = 0.01
include("kinetic_intensity_data.jl")
intensity(x_A,x_B,x_D) = x_A + 2/21*x_B + 2/21*x_D

# Solve using Ipopt
using Ipopt
model = Model(Ipopt.Optimizer)
decision_vars(o) # Displays: x_Z(t), x_Y(t), x_D(t), x_B(t), x_A(t), k_2f, k_3f, k_4
# FIRST, create discretized differntial state decision variables "z"
N = Int(floor((tspan[2] - tspan[1])/tstep))+1
V = length(unknowns(o))
zL = zeros(V) # lower bound on z
zU = [140.0, 0.4, 140.0, 140.0, 140.0] # upper bound on z
@variable(model, zL[i] <= z[i in 1:V,1:N] <= zU[i]) # ̇z = (x_Z(t), x_Y(t), x_D(t), x_B(t), x_A(t))
# SECOND, create free design decision variables "p"
pL = [10, 10, 0.001] # lower bound on p
pU = [1200, 1200, 40] # upper bound on p
@variable(model, pL[i] <= p[i=1:3] <= pU[i]) # p = (k_2f, k_3f, k_4)
register_odesystem(model, o, tspan, tstep, "EE") # Try changing between "EE" and "IE"!
@objective(model, Min, sum((intensity(z[5,i],z[4,i],z[3,i]) - data[i-1])^2 for i in 2:N))
JuMP.optimize!(model)
```

## References
1. Y. Ma, S. Gowda, R. Anantharaman, C. Laughman, V. Shah, and C. Rackauckas, **ModelingToolkit: A composable graph transformation system for equation-based modeling**, 2021.
2. M. Lubin, O. Dowson, J. Dias Garcia, J. Huchette, B. Legat, and J. P. Vielma, **JuMP 1.0: Recent improvements to a modeling language for mathematical optimization**, *Mathematical Programming Computation*, vol. 15, p. 581-589, 2023.
3. M. Wilhelm and M. Stuber, **EAGO.jl Easy Advanced Global Optimization in Julia**, *Optimization Methods and Software*, vol. 37, no. 2, pp. 425-450, 2022.
4. A. Wächter, & L. T. Biegler, **On the implementation of an interior-point filter line-search algorithm for large-scale nonlinear programming**, *Mathematical programming*, vol. 106, no. 1, pp. 25-57, 2006.

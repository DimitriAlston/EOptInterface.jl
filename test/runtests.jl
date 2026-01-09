using Test, EOptInterface, JuMP, ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@testset "EOptInterface.jl" begin
    @mtkmodel AlgebraicTest begin
        @parameters begin
            k₁ = 0.40
            k₂ = 0.055
            V_A = 8.937e-2
            V_B = 1.018e-1
            V_C = 1.130e-1
            V
            F₁
        end
        @variables begin
            F₂(t); F₃(t); F₄(t); F₅(t); F₆(t); F₇(t)
            y_3A(t); y_3B(t); y_3C(t)
            y_4B(t); y_4C(t)
        end
        begin
            r₁ = (k₁*y_3A)/(y_3A*V_A + y_3B*V_B + y_3C*V_C)
            r₂ = (k₂*y_3B)/(y_3A*V_A + y_3B*V_B + y_3C*V_C)
        end
        @equations begin
            F₅ ~ (y_4B * F₄)
            F₁ + F₇ ~ F₂
            (y_3A * F₃) ~ F₂ - r₁*V
            (y_3B * F₃) ~ (r₁ - r₂)*V
            (y_3C * F₃) ~ r₂*V
            F₃ ~ F₄ + F₇
            (y_3B * F₃) ~ (y_4B * F₄)
            (y_3C * F₃) ~ (y_4C * F₄)
            F₄ ~ F₅ + F₆
            y_3A + y_3B + y_3C ~ 1
            y_4B + y_4C ~ 1
        end
    end
    @mtkcompile n = AlgebraicTest()
    g1 = 25 - n.F₅
    g2 = 475/3600 - n.V/(n.F₃*(n.y_3A*n.V_A + n.y_3B*n.V_B + n.y_3C*n.V_C))
    f_CSTR = (25764 + 8178*n.V)/2.5
    s1cap = 132718 + n.F₃*(369*n.y_3A - 1113.9*n.y_3B)
    s2cap = 25000 + n.F₄*(6984.5*n.y_4B - 3869.53*n.y_4C^2)
    s1op = n.F₃*(3+36.11*n.y_3A + 7.71*n.y_3B)*26.32e-3
    s2op = n.F₄*(26.21 + 29.45*n.y_4B)*26.32e-3;
    f_Sep = (s1cap+s2cap)/2.5 + 0.52*(s1op+s2op)
    obj = f_CSTR + f_Sep
    using Ipopt
    model = JuMP.Model(Ipopt.Optimizer)
    EOptInterface.decision_vars(n)
    xL = zeros(6)
    xU = [1, 100, 1, 1, 10, 100]
    @variable(model, xL[i] <= xvar[i in 1:6] <= xU[i])
    EOptInterface.register_nlsystem(model, n, obj, [g1, g2])
    JuMP.optimize!(model)
    @test JuMP.value.(xvar)[6] == 26.316700011101723
    @test JuMP.termination_status(model) == JuMP.LOCALLY_SOLVED
    @test EOptInterface.full_solutions(model, n)[n.y_3C] == 0.01385685924868818
end

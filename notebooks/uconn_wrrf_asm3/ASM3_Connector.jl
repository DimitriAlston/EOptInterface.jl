using ModelingToolkit

@connector MaterialStream begin
    flow_rate(t), [input = true]
    S_O(t), [input = true]
    S_I(t), [input = true]
    S_S(t), [input = true]
    S_NH(t), [input = true]
    S_N2(t), [input = true]
    S_NO(t), [input = true]
    S_ALK(t), [input = true]
    X_I(t), [input = true]
    X_S(t), [input = true]
    X_H(t), [input = true]
    X_STO(t), [input = true]
    X_A(t), [input = true]
    X_TS(t), [input = true]
end

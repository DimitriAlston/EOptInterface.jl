# reference
# [1] Gujer, W., Henze, M., Mino, T., & Van Loosdrecht, M. (1999). Activated sludge model No. 3. Water science and technology, 39(1), 183-193.

# [2] Henze, M., Gujer, W., Mino, T., & Van Loosedrecht, M. (2006). Activated sludge models ASM1, ASM2, ASM2d and ASM3. IWA publishing.


function theta_T(par::Any)
    return log(par[2]/par[1])/(10.0)
end

function k_fun(par::Any,T::Any)
    return par[2]*exp(theta_T(par)*(T-20))
end
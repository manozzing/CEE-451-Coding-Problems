## Shooting method
function ShootingMethod(x₁::Float64, xₙ₊₁::Float64, f₁::Float64, g₁::Float64, F₁::Float64, G₁::Float64, N::Int64)
    Δx = (xₙ₊₁-x₁)/N
    x = zeros(N+1)
    f = zeros(N+1)
    g = zeros(N+1)
    F = zeros(N+1)
    G = zeros(N+1)

    x[1] = x₁
    f[1] = f₁
    g[1] = g₁
    F[1] = F₁
    G[1] = G₁

    s = g₁

    for i in 1:N
        x[i+1] = x[i] + Δx
        f[i+1] = f[i] + g[i]*Δx
        g[i+1] = g[i] - f[i]*Δx
        F[i+1] = F[i] + G[i]*Δx
        G[i+1] = G[i] - F[i]*Δx
    end
    
    sᴺᴱᵂ = s - (f[N+1]-1)/F[N+1]

    return sᴺᴱᵂ
end

#Error calculation for Shooting method with different discretization
function ϵ_Shooting_N(s::Float64, N::Int64)
    s1 = ShootingMethod(0.0, 1.0, 0.0, s, 0.0, 1.0, N)
    s2 = ShootingMethod(0.0, 1.0, 0.0, s, 0.0, 1.0, N+1)
    ϵₛ = abs(2*(s2-s1)/(s2+s1))
    return ϵₛ, s2
end

ϵ_Shooting(1.0, 10)

## Find N for grid_invariance ϵₛ < δ
function find_N_Shooting(s::Float64, δ::Float64)
    ϵₛ = 1.0
    N = 2
    while (ϵₛ > δ)
        ϵₛ, s_final = ϵ_Shooting_N(s, N)
        N = N + 1
        if N > 10000
            break
        end
    end
    return N, s_final
end

find_N_Shooting(1.0, 0.0001)

#Error calculation for Shooting method with different repeatation
function ϵ_Shooting_M(s::Float64, M::Int64)
    s1 = ShootingMethod(0.0, 1.0, 0.0, s, 0.0, 1.0, N)
    s2 = ShootingMethod(0.0, 1.0, 0.0, s, 0.0, 1.0, N+1)
    ϵₛ = abs(2*(s2-s1)/(s2+s1))
    return ϵₛ
end

## Find M for Repeatation ϵₛ < δ
function find_M_Shooting(s::Float64, δ::Float64, N::Int64)
    ϵₛ = 1.0
    M = 1
    ss = zeros(10000)
    ss[1] = s
    
    for i in 1:10000
        ss[i+1] = ShootingMethod(0.0, 1.0, 0.0, ss[i], 0.0, 1.0, N)

        if abs(2*(ss[i+1]-ss[i])/(ss[i+1]+ss[i])) < δ
            return ss[i+1], ss[i], i
            break
        end
    end
end

find_M_Shooting(1000.0, 0.0000000000000000000000000001, 3)

#println("Execution goes well without syntax error")
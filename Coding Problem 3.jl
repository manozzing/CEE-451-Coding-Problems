## Define f(x,y)
function fxy(x::Float64, y::Float64, k::Float64)
    return -k * sqrt(y) / log(x+1)
end

## Euler step method
function EulerStepMethod(xₐ::Float64, yₐ::Float64, xᵦ::Float64, N::Int64, k::Float64)
    Δx = (xᵦ-xₐ)/N
    x = zeros(N)
    y = zeros(N)
    x[1] = xₐ + Δx
    y[1] = yₐ + fxy(xₐ, yₐ, k)*Δx
    for i in 1:N-1
        x[i+1] = x[i] + Δx
        y[i+1] = y[i] + fxy(x[i], y[i], k)*Δx
    end
    yᵦ = y[N]
    y_c = y[N÷2]
    return yᵦ, y_c
end

println("Execution goes well without syntax error")

EulerStepMethod(1.0, 4.0, 5.0, 8, 0.25)

## Error Calculation
function ϵ_Euler(xₐ::Float64, yₐ::Float64, xᵦ::Float64, N::Int64, k::Float64)
    yᵦ_N1, y_c_N1 = EulerStepMethod(xₐ, yₐ, xᵦ, N, k)
    yᵦ_N2, y_c_N2 = EulerStepMethod(xₐ, yₐ, xᵦ, N+2, k)
    ϵ_c = abs(2*(y_c_N2-y_c_N1)/(y_c_N2+y_c_N1))
    ϵᵦ = abs(2*(yᵦ_N2-yᵦ_N1)/(yᵦ_N2+yᵦ_N1))
    return ϵᵦ, ϵ_c
end

ϵ_Euler(1.0, 4.0, 5.0, 26, 0.25)

## Predictor-Correction Method
function pcm(xₐ::Float64, yₐ::Float64, xᵦ::Float64, N::Int64, k::Float64)
    Δx = (xᵦ-xₐ)/N
    x = zeros(N)
    y = zeros(N)
    yₚ = zeros(N)
    x[1] = xₐ + Δx
    yₚ[1] = yₐ + fxy(xₐ, yₐ, k)*Δx 
    y[1] = yₐ + 0.5(fxy(xₐ, yₐ, k) + fxy(xₐ, yₚ[1], k))*Δx
    for i in 1:N-1
        x[i+1] = x[i] + Δx
        yₚ[i+1] = y[i] + fxy(x[i], y[i], k)*Δx
        y[i+1] = y[i] + 0.5*(fxy(x[i], y[i], k) + fxy(x[i+1], yₚ[i+1], k))*Δx
    end
    yᵦ = y[N]
    y_c = y[N÷2]
    return yᵦ, y_c
end

## Error Calculation for pcm
function ϵ_pcm(xₐ::Float64, yₐ::Float64, xᵦ::Float64, N::Int64, k::Float64)
    yᵦ_N1, y_c_N1 = pcm(xₐ, yₐ, xᵦ, N, k)
    yᵦ_N2, y_c_N2 = pcm(xₐ, yₐ, xᵦ, N+2, k)
    ϵ_c = abs(2*(y_c_N2-y_c_N1)/(y_c_N2+y_c_N1))
    ϵᵦ = abs(2*(yᵦ_N2-yᵦ_N1)/(yᵦ_N2+yᵦ_N1))
    return ϵᵦ, ϵ_c
end

ϵ_Euler(1.0, 4.0, 5.0, 2, 0.25)

## Repeat until ϵ < δ: Euler Step Method
function find_N_Euler(δ::Float64, k::Float64)
    j = 0
    ϵ_c = 1.0
    ϵᵦ = 1.0
    N = 2
    while (ϵ_c > δ || ϵᵦ > δ)
        ϵᵦ, ϵ_c = ϵ_Euler(1.0, 4.0, 5.0, N, k)
        N = N + 2
        if N > 100
            break
        end
    end
    return N
end

find_N_Euler(0.001, 0.25)



## Repeat until ϵ < δ: Predictor-Correction Method
function find_N_pcm(δ::Float64, k::Float64)
    j = 0
    ϵ_c = 1.0
    ϵᵦ = 1.0
    N = 2
    while (ϵ_c > δ || ϵᵦ > δ)
        ϵᵦ, ϵ_c = ϵ_pcm(1.0, 4.0, 5.0, N, k)
        N = N + 2
        if N > 100
            break
        end
    end
    return N
end

find_N_pcm(0.001, 0.25)

## Making matrix for the outputs with repeatation: Euler Step Method
function Euler_Matrix_Filling(k::Float64)
    Euler_matrix = zeros(find_N_Euler(0.001, k)÷2, 5)
    for i in 1:find_N_Euler(0.001, k)÷2
        Euler_matrix[i, 1] = 2*i
        Euler_matrix[i, 2], Euler_matrix[i, 3] = EulerStepMethod(1.0, 4.0, 5.0, 2*i, k)
        Euler_matrix[i, 4], Euler_matrix[i, 5] = ϵ_Euler(1.0, 4.0, 5.0, 2*i, k)
    end
    return Euler_matrix
end

Euler_Matrix_Filling(0.25)
Euler_Matrix_Filling(0.5)

## Making matrix for the outputs with repeatation: Perdictor-Correction Method
function PCM_Matrix_Filling(k::Float64)
    PCM_matrix = zeros(find_N_pcm(0.001, k)÷2, 5)
    for i in 1:find_N_pcm(0.001, k)÷2
        PCM_matrix[i, 1] = 2*i
        PCM_matrix[i, 2], PCM_matrix[i, 3] = pcm(1.0, 4.0, 5.0, 2*i, k)
        PCM_matrix[i, 4], PCM_matrix[i, 5] = ϵ_pcm(1.0, 4.0, 5.0, 2*i, k)
    end
    return PCM_matrix
end

PCM_Matrix_Filling(0.25)
PCM_Matrix_Filling(0.5)
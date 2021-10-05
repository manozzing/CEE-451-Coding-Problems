# Making blank arrayes for I and f
I = zeros(10000)
f = zeros(10000)

# Defining function for integration of exp(-λx)
function integral_1(xₐ::Float64, xᵦ::Float64, λ::Float64, N::Int64)
    
    for n in 2:N
        f_1 = exp(-λ*xₐ)
        f_end = exp(-λ*xᵦ)
        I[1] = 0.5*(f_1 + f_end)*(xᵦ - xₐ)
        dx = (xᵦ - xₐ)/n
        I[n] = 0.5*(f_1 + f_end)*dx
        for i in 1:n-1
            x = xₐ + dx*i
            f = exp(-λ*x)
            I[n] = I[n] + f*dx
        end
    end
    return I[N]
end

# Defining function for finding N when the integration converges
function finding_n(δ::Float64, Nₘ::Int64)
    N = 2
    while (abs((integral_1(0.0,4.0,1.0,N)-integral_1(0.0,4.0,1.0,N-1))/
            (integral_1(0.0,4.0,1.0,N)+integral_1(0.0,4.0,1.0,N-1))) >= δ && N < Nₘ)
            N = N+1
    end
    return N, integral_1(0.0,4.0,1.0,N)
end

# finding n when δ = 0.001, Nₘ = 10000
finding_n(0.001, 10000)

# making array of integration per different N until the integration converges
I_1 = zeros(finding_n(0.001, 10000)[1]-1)
for i in 2:finding_n(0.001, 10000)[1]
    I_1[i-1] = integral_1(0.0, 4.0, 1.0, i)
end
I_1

# Defining function for integration of exp(sin(x)+sqrt(x))
function integral_2(xₐ::Float64, xᵦ::Float64, N::Int64)
    for n in 2:N
        f_1 = exp(sin(xₐ)+sqrt(xₐ))
        f_end = exp(sin(xᵦ)+sqrt(xᵦ))
        I[1] = 0.5*(f_1 + f_end)*(xᵦ - xₐ)
        dx = (xᵦ - xₐ)/n
        I[n] = 0.5*(f_1 + f_end)*dx
        for i in 1:n-1
            x = xₐ + dx*i
            f = exp(sin(x)+sqrt(x))
            I[n] = I[n] + f*dx
        end
    end
    return I[N]
end

# Defining function for finding N when the integration converges
function finding_n_2(δ::Float64, Nₘ::Int64)
    N = 2
    while (abs((integral_2(0.0,4.0,N)-integral_2(0.0,4.0,N-1))/
            (integral_2(0.0,4.0,N)+integral_2(0.0,4.0,N-1))) >= δ && N < Nₘ)
            N = N+1
    end
    return N, integral_2(0.0,4.0,N)
end

# finding n when δ = 0.001, Nₘ = 10000
finding_n_2(0.001, 10000)

# making array of integration per different N until the integration converges
I_2 = zeros(finding_n_2(0.001, 10000)[1]-1)
for i in 2:finding_n_2(0.001, 10000)[1]
    I_2[i-1] = integral_2(0.0, 4.0, i)
end
I_2

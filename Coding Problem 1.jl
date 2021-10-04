### Import required package 
using Plots

### Q1. Defining function for Newton-Rhapson method
function Newton_Raphson(A::Float64, B::Float64, x₀::Float64, Nₘ::Int64, δ::Float64)
    x = x₀
    N = 0
    y = y = A - B*x^2*log(x+1) ::Float64
    while (abs(y) > δ && N < Nₘ)
        y = A - B*x^2*log(x+1)
        dydx = -2*B*x*log(x+1) - B*x^2/(x+1)
        x = x - y/dydx
        N = N+1
    end

    return x, N
end

### Q2. Consider the case A=1, B=0.2, x₀=0.5
Newton_Raphson(1.0, 0.2, 0.5, 10000, 0.000000001)

### Q3. a) Compute x as a function of B over the range 0.2 ≤ B ≤ 2 for A = 4
# Making arrays for B and x
B_Array = zeros(1001)
x_Array = zeros(1001)

# For loop for making 1000-steps of ranges for B and x
for i in 1:1001
    B_Array[i] = 0.2 + 1.8*(i-1)/1000
    x_Array[i] = Newton_Raphson(4.0, 0.2 + 1.8*(i-1)/1000, 0.5, 10000, 0.000000001)[1]
end

# Plotting
plot(B_Array, x_Array, legend = false, title = "The root x per B in Question 3", xlabel = "B", ylabel = "x")

### Q3. b) Compute x as a function of A over the range 1 ≤ A ≤ 6 for A = 0.5
# Making arrays for A and x
A_Array = zeros(1001)
x2_Array = zeros(1001)

# For loop for making 1000-steps of ranges for B and x
for i in 1:1001
    A_Array[i] = 1.0 + 5.0*(i-1)/1000
    x2_Array[i] = Newton_Raphson(1.0 + 5.0*(i-1)/1000, 0.5, 0.5, 10000, 0.000000001)[1]
end

# Plotting
plot(A_Array, x2_Array, legend = false, title = "The root x per A in Question 3", xlabel = "A", ylabel = "x")

### Thank you for your consideration.
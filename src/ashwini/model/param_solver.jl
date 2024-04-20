include("diffusivity.jl")

"""
## param_solver
The key function to get all the parameters required to solve 
for c₂ at a given timestep in the timeseries.

Returns the 16 parameters: n₀(t), n₁(t), n₂(t), n₃(t),
p₀(t), p₁(t), p₂(t), p₃(t), s₀(t), s₁(t), s₂(t), <c₂>ₙ,
<c₂>ₚ , <c₂>ₛ, c₂|ₓ₌₀ and c₂|ₓ₌ₗ in this order

### Params:
    - I: The current at the timestep.
    - Δt: The span of the timestep.
    - c₂ₙₚ: The previous average concentration in negative electrode
    - c₂ₛₚ: The The previous average concentration in separator
    - c₂ₚₚ: The previous average concentration in positive electrode
    - c₂ₙ₀ₚ: The previous concentration in negative electrode @ x=0
    - c₂ₚₗₚ: The previous average concentration in positive electrode @x=L
    - jn_0: The reaction rate @x=0
    - jp_l: The reaction rate @x=L
"""
function param_solver(I, Δt, c₂ₙₚ, c₂ₛₚ, c₂ₚₚ, c₂ₙ₀ₚ, c₂ₚₗₚ, jn_0, jp_l)

    A = zeros(16, 16)
    B = zeros(16, 1)

    D₂ₙ = d₂ₙ(c₂ₙₚ)
    D₂ₛ = d₂ₛ(c₂ₛₚ)
    D₂ₚ = d₂ₚ(c₂ₚₚ)
    
    # Eqn 73
    A[1, 2] = 1
    A[2, 6] = 1
    
    # Eqn 74
    A[3, 1] = -1; A[3, 9] = 1

    # Eqn 75
    A[4, 1] = -1; A[4, 5] = 1; A[4, 10] = -lₛ; A[4, 11] = -lₛ^2

    # Eqn 76
    A[5, 3] = -2*D₂ₙ*lₙ
    A[5, 4] = -3*D₂ₙ*lₙ^2
    A[5, 10] = -D₂ₛ
    A[5, 11] = -2*lₙ*D₂ₛ

    # Eqn 77
    A[6, 7] = 2*D₂ₚ*lₚ
    A[6, 8] = 3*D₂ₚ*lₚ^2
    A[6, 10] = -D₂ₛ
    A[6, 11] = -2*(lₙ+lₛ)*D₂ₛ

    # Eqn 78, 79, 80
    A[7, 1] = 1; A[7, 3] = (2*lₙ^2)/3; A[7, 4] = (3*lₙ^3)/4; A[7, 12] = -1
    A[8, 5] = 1; A[8, 7] = (2*lₚ^2)/3; A[8, 8] = (3*lₚ^3)/4; A[8, 13] = -1
    A[9, 9] = 1; A[9, 10] = lₛ/2; A[9, 11] = (2*lₛ^2)/3; A[9, 14] = -1

    # Eqn 81
    A[10, 12] = lₙ*ϵ₂ₙ; A[10, 13] = lₚ*ϵ₂ₚ; A[10, 14] = lₛ*ϵ₂ₛ

    # Eqn 82 & 83
    A[11, 3] = 2 * D₂ₙ * Δt/ϵ₂ₙ
    A[11, 4] = 3 * lₙ * D₂ₙ * Δt/ϵ₂ₙ
    A[11, 12] = 1

    A[12, 7] = 2 * D₂ₚ * Δt/ϵ₂ₚ
    A[12, 8] = 3 * lₚ * D₂ₚ * Δt/ϵ₂ₚ
    A[12, 13] = 1

    # Eqn 84 & 85
    A[13, 3] = 2*D₂ₙ*Δt/ϵ₂ₙ; A[13, 15] = 1

    A[14, 7] = 2*D₂ₚ*Δt/ϵ₂ₚ; A[14, 16] = 1

    # Eqn 86 & 87
    A[15, 1] = 1; A[15, 3] = lₙ^2; A[15, 4] = lₙ^3; A[15, 15] = -1
    A[16, 5] = 1; A[16, 7] = lₚ^2; A[16, 8] = lₚ^3; A[16, 16] = -1


    # Filling the B matrix to solve the eqn AX = B
    B[10] = c₂₀ * (lₙ*ϵ₂ₙ + lₛ*ϵ₂ₛ + lₚ*ϵ₂ₚ)
    B[11] = c₂ₙₚ + aₙ*(1 - t₊)*j̅ₙ(I)*Δt/ϵ₂ₙ
    B[12] = c₂ₚₚ+ aₚ*(1 - t₊)*j̅ₚ(I)*Δt/ϵ₂ₚ
    B[13] = c₂ₙ₀ₚ + aₙ*(1 - t₊)*jn_0*Δt/ϵ₂ₙ
    B[14] = c₂ₚₗₚ + aₚ*(1 - t₊)*jp_l*Δt/ϵ₂ₚ
    
    # for i = 1:16
    #     println(A[i, :])
    # end

    # for i = 1:16
    #     println(B[i])
    # end

    A\B
end
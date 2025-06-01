using LinearAlgebra

#function that returns the channel for a given ρ, p, and j
function ε_depo(ρ,p,j)
    return (1-p)*ρ+p/(2j+1)*Matrix(1.0I, 2j+1, 2j+1)
end

# Function to rotate a vector v to align the z-axis with n
function rotate_to_axis(v::Vector{Float64}, n::Vector{Float64})
    n = normalize(n)  # Ensure n is a unit vector
    z_axis = [0.0, 0.0, 1.0]  # Reference z-axis

    # Find rotation axis as cross product of z-axis and target n
    axis = cross(z_axis, n)
    axis_norm = norm(axis)

    if axis_norm ≈ 0  # If axis is nearly zero, v is already aligned
        return v
    end

    axis /= axis_norm  # Normalize rotation axis
    θ = acos(dot(z_axis, n))  # Angle between z and n

    # Apply Rodrigues' rotation formula
    v_rot = v * cos(θ) + cross(axis, v) * sin(θ) + axis * dot(axis, v) * (1 - cos(θ))

    return v_rot  # Return as a vector
end

# Function to generate a random unit vector within a solid angle Ω around a given unit vector n
function random_unit_vector(Ω::Float64, n::Vector{Float64})
    # Compute max polar angle from solid angle Ω
    cosθ_max = 1 - Ω / (2π)  # Derived from solid angle fraction
    cosθ = rand() * (1 - cosθ_max) + cosθ_max  # Uniform sampling in [cosθ_max, 1]
    θ = acos(cosθ)

    # Uniform azimuthal angle
    φ = rand() * 2π

    # Convert to Cartesian coordinates in the z-aligned frame
    v = [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]

    # Rotate to align z-axis with n
    return rotate_to_axis(v, n)
end

# Function to generate a mixture of states for the AM channel
function ε_AM(ρ,θ,n,Ω,nvecs,N) 
    J_x=colspinjm("x",N)
    J_y=colspinjm("y",N)
    J_z=colspinjm("z",N)
    n = normalize(n)  # unit vector
    σ = 0.0*ρ
    # Generate random unit vectors
    vectors = [random_unit_vector(Ω, n) for _ in 1:nvecs]
    # Create mixture of states
    for nk in vectors
        Jnk = nk[1]*J_x + nk[2]*J_y + nk[3]*J_z
        Uθ = exp(-im * θ * Jnk)
        σ = σ + Uθ*ρ*Uθ'
    end
    return σ/tr(σ) # trace-1 state
end

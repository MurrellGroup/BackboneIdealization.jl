module BackboneIdealization

using Backboner
using Flux
using Optimisers
using ParameterSchedulers
using Statistics

export idealize

function length_loss(m::Matrix{<:Real}, ideal::Vector{<:Real}; mask::Vector{Bool} = ones(Bool, size(m, 2)-1))
    sum(abs2, (get_bond_lengths(Backbone(m)) .- ideal) .* mask)
end

function angle_loss(m::Matrix{<:Real}, ideal::Vector{<:Real}; mask::Vector{Bool} = ones(Bool, size(m, 2)-2))
    angles = get_bond_angles(Backbone(m))
    sum(abs2, (angles .- ideal) .* mask)
end

function total_loss(
    coords::Matrix{T}, offsets::Matrix{T},
    ideal_lengths::Vector{T}, ideal_angles::Vector{T},
    length_mask::Vector{Bool} = ones(Bool, size(coords, 2)-1),
    angle_mask::Vector{Bool} = ones(Bool, size(coords, 2)-2),
) where T <: Real
    offset_coords = coords + offsets
    length_loss_val = length_loss(offset_coords, ideal_lengths, mask=length_mask)
    angle_loss_val = angle_loss(offset_coords, ideal_angles, mask=angle_mask)
    return length_loss_val + angle_loss_val
end

function idealize(
    backbone::Backbone{T}, ideal_lengths::Vector{T}, ideal_angles::Vector{T};
    mask_tolerance = 0.5, n_iterations = 300,
) where T <: Real
    scale_factor = mean(ideal_lengths)

    coords = backbone.coords / scale_factor
    bonds = ChainedBonds(Backbone(coords))

    ideal_lengths_cycled = [ideal for (_, ideal) in zip(get_bond_lengths(bonds), Iterators.cycle(ideal_lengths / scale_factor))]
    ideal_angles_cycled = [ideal for (_, ideal) in zip(get_bond_angles(bonds), Iterators.cycle(ideal_angles))]

    length_mask = collect(abs.(get_bond_lengths(bonds) .- ideal_lengths_cycled) .< mask_tolerance)
    angle_mask = collect(length_mask[1:end-1] .& length_mask[2:end])
    
    offsets = zeros(T, size(coords))
    state = Flux.setup(Adam(), offsets)
    scheduler = Exp(start=0.03, decay=0.9)

    for (eta, i) in zip(scheduler, 1:n_iterations)
        ∇offsets, = Flux.gradient(offsets) do Δ
            log(total_loss(coords, Δ, ideal_lengths_cycled, ideal_angles_cycled, length_mask, angle_mask))
        end

        Flux.adjust!(state, eta)
        Flux.update!(state, offsets, ∇offsets)
    end

    new_coords = coords .+ offsets
    rescaled_coords = new_coords * scale_factor
    return Backbone(rescaled_coords)
end

end

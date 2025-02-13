module BackboneIdealization

using Backboner
using Flux
using Optimisers
using ParameterSchedulers

export idealize

function idealize(
    backbone::Backbone, ideal_lengths::AbstractVector, ideal_angles::AbstractVector;
    mask_tolerance = 0.5, n_iterations = 300,
)
    scale_factor = sum(ideal_lengths) / length(ideal_lengths)

    coords = backbone.coords / scale_factor
    bonds = ChainedBonds(Backbone(coords))

    ideal_lengths_cycled = [ideal for (_, ideal) in zip(get_bond_lengths(bonds), Iterators.cycle(ideal_lengths / scale_factor))]
    ideal_angles_cycled = [ideal for (_, ideal) in zip(get_bond_angles(bonds), Iterators.cycle(ideal_angles))]

    length_mask = collect(abs.(get_bond_lengths(bonds) - ideal_lengths_cycled) .< mask_tolerance)
    angle_mask = collect(length_mask[1:end-1] .& length_mask[2:end])
    
    Δcoords = zeros(eltype(coords), size(coords))
    state = Flux.setup(Adam(), Δcoords)
    scheduler = Exp(start=0.03, decay=0.9)

    for (η, i) in zip(scheduler, 1:n_iterations)
        ∇Δcoords, = Flux.gradient(Δcoords) do Δcoords
            new_backbone = Backbone(coords + Δcoords)
            length_loss = sum(abs2, (get_bond_lengths(new_backbone) - ideal_lengths_cycled) .* length_mask)
            angle_loss = sum(abs2, (get_bond_angles(new_backbone) - ideal_angles_cycled) .* angle_mask)
            log(length_loss + angle_loss)
        end

        Flux.adjust!(state, η)
        Flux.update!(state, Δcoords, ∇Δcoords)
    end

    return Backbone((coords + Δcoords) * scale_factor)
end

end

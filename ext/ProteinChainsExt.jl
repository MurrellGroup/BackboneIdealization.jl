module ProteinChainsExt

using BackboneIdealization
using Backboner
using ProteinChains

function BackboneIdealization.idealize(
    backbone::Backbone{T},
    geometry::BackboneGeometry = DEFAULT_BACKBONE_GEOMETRY;
    kwargs...
) where T<:Real
    BackboneIdealization.idealize(
        backbone,
        T[geometry.N_Ca_length, geometry.Ca_C_length, geometry.C_N_length],
        T[geometry.N_Ca_C_angle, geometry.Ca_C_N_angle, geometry.C_N_Ca_angle];
        kwargs...
    )
end

# warning: strips non-backbone atoms
function BackboneIdealization.idealize(
    chain::ProteinChain,
    args...; kwargs...
)
    ProteinChain(
        chain.id,
        get_atoms(idealize(Backbone(chain), args...; kwargs...)),
        chain.sequence,
        chain.numbering,
        chain.ins_codes;
        ProteinChains.propertypairs(chain, ProteinChains.NoFields())...
    )
end

function BackboneIdealization.idealize(
    structure::ProteinStructure,
    args...; kwargs...
)
    ProteinStructure(
        structure.name,
        structure.atoms,
        map(chain -> idealize(chain, args...; kwargs...), structure.chains);
        ProteinChains.propertypairs(structure, ProteinChains.NoFields())...
    )
end

end
using BackboneIdealization
using Test

using Backboner
using ProteinChains

@testset "BackboneIdealization.jl" begin
    
    @testset "idealization" begin
        chain = pdb"1qg8"A # has missing residues
        geometry = DEFAULT_BACKBONE_GEOMETRY
        chain_ideal = idealize(chain, geometry)
        @test all(≈(geometry.C_N_length), sort(get_bond_lengths(ChainedBonds(chain_ideal))[3:3:end])[1:end-2]) # two gaps

        scaled_geometry = BackboneGeometry(N_Ca_length=0.1geometry.N_Ca_length, Ca_C_length=0.1geometry.Ca_C_length, C_N_length=0.1geometry.C_N_length)
        backbone_scaled = Backbone(Backbone(chain).coords * 0.1)
        backbone_ideal_scaled = idealize(backbone_scaled, scaled_geometry)
        @test all(≈(scaled_geometry.C_N_length), sort(get_bond_lengths(backbone_ideal_scaled)[3:3:end])[1:end-2]) # two gaps

        @test idealize(pdb"1ASS")["A"] == idealize(pdb"1ASS"A)
    end

end

# BackboneIdealization

[![Build Status](https://github.com/MurrellGroup/BackboneIdealization.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MurrellGroup/BackboneIdealization.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/MurrellGroup/BackboneIdealization.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/BackboneIdealization.jl)

## Installation

```julia
using Pkg
pkg"registry add https://github.com/MurrellGroup/MurrellGroupRegistry"
pkg"add BackboneIdealization"
```

## Usage

```julia
julia> using BackboneIdealization, ProteinChains, Backboner

julia> chain = pdb"1ASS"A;
[ Info: Downloading file from PDB: 1ASS

julia> idealized = idealize(chain, BackboneGeometry());

julia> round.(get_bond_lengths(ChainedBonds(idealized)), digits=12)
455-element Vector{Float64}:
 1.46
 1.52
 1.33
 1.46
 1.52
 ⋮
 1.33
 1.46
 1.52
 1.33
 1.46
 1.52
```

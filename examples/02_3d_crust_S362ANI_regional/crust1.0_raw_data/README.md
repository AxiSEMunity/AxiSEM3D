# CRUST1.0 Raw Data

This directory contains the distributed CRUST1.0 files and the companion Fortran 77 routines used to extract either a single 1D profile, map products, or XYZ tables.

## Included files

- `crust1.bnds`: boundary topography.
- `crust1.vp`: P-wave speed.
- `crust1.vs`: S-wave speed.
- `crust1.rho`: density.

The original CRUST1.0 package also referenced three helper routines:

- `getCN1point.f`: extract a 1D profile at a latitude/longitude point.
- `getCN1maps.f`: generate map products.
- `getCN1xyz.f`: generate XYZ output.

## Model description

CRUST1.0 is an 8-layer model defined on 1x1 degree cells. The published values represent averages over each cell, so queries should use the cell midpoint. For example, the cell spanning 5 to 6 degrees latitude and 150 to 151 degrees longitude should be queried at 5.5 degrees latitude and 150.5 degrees longitude.

This was described by the original authors as a prototype model, with feedback directed to `glaske-at-ucsd.edu` and a planned CRUST1.1 update after the July 15, 2013 initial release.

## Layers

1. Water
2. Ice
3. Upper sediments
4. Middle sediments
5. Lower sediments
6. Upper crystalline crust
7. Middle crystalline crust
8. Lower crystalline crust
9. Sub-Moho properties: `V_Pn`, `V_Sn`, and density

The crustal layers are incomplete in some cells for sedimentary properties. The ninth layer is associated with LLNL model G3Cv3 on continents and a thermal model in the oceans.

## File format

The model spans 89.5 to -89.5 degrees latitude and -179.5 to 179.5 degrees longitude. Longitudes are the inner loop, so each latitude row stores all longitudes before moving to the next latitude. The files start at 89.5 degrees north and 179.5 degrees west.

Each line contains the nine layer values for one parameter in one cell. Every model file therefore contains `360 x 180 = 64800` lines.

## Derived outputs from `getCN1maps.f`

- `map-bd[x]`, where `x = 1..9`:
  1. top of water
  2. bottom of water
  3. bottom of ice
  4. bottom of sediments 1
  5. bottom of sediments 2
  6. bottom of sediments 3
  7. bottom of crystalline crust 1
  8. bottom of crystalline crust 2
  9. bottom of crystalline crust 3, which is the Moho depth
- `map-vp[x]`, where `x = 1..8` gives in-crust `vp` and `x = 9` gives `VPn`
- `map-vs[x]`, following the same convention for `vs`
- `map-ro[x]`, following the same convention for density
- `sedthk`: sediment thickness
- `crsthk`: crustal thickness without water

The `getCN1xyz.f` output uses the same data but emits one line per cell with longitude, latitude, and value.

The generated maps span 89.5 to -89.5 degrees latitude and -179.5 to 179.5 degrees longitude.
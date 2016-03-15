# pyaeromesh
Some incomplete but (semi)usable tools for making OpenFOAM blockmeshes for aerofoils, wings and turbomachinery cascades.
Requires pyFoam.
The wing generation library (aerofoilmesher.py) is the most complete.
It can produce a sort of C-mesh for 2D aerofoils or 3d wings of arbitary aerofoil profile, dihedral, sweep and angle of attack,
but will produce horribly skew cells if any of these angles are too large (probably far beyond stall though).
Run test_afmesh.py to see how it works (wing_template [3D case] or wing2d must be in the same directory). 


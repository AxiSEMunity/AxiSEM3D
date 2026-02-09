# Meshes

- See here details of the mesher used by AxiSEM3D, [Salvus Mesher Lite](https://gitlab.com/swp_ethz/public/SalvusMeshLite)
  .

- The installation page for AxiSEM3D here also contains an excellent
  section on the mesher:
  <https://github.com/AxiSEMunity/AxiSEM3D/wiki/mesher>.

In general, you have a choice between Cartesian (box-shaped) meshes, and
D-shaped (semi-circular) meshes. These obviously solve for the same
physics; and you can actually adapt the Cartesian simulation to account
for the curvature of the Earth’s surface (to move from local to regional
scales) or reduce the D-shaped mesh to a wedge-shaped mesh (to move from
global to regional scales).

It is worth bearing in mind that unless you use a full D-shaped (i.e.,
180°/full semi-circle) mesh, you need to think about the ‘false’
boundaries that you impose at the edges of your box or wedge. There is
more on this in the mesher description in, but in short you can use an
absorbing or reflecting boundary. In the former case you have a few
options depending on what you need, the latter is far simpler and
somewhat cheaper but will also contaminate your seismograms with
reflections from the edge of the box if you let the simulation run on
for too long.
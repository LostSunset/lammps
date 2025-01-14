.. index:: fix surface/local

fix surface/local command
===============

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID surface/local input args input args ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* surface/local = style name of this fix command
* input = zero or more input keywords can be specified

  .. parsed-literal::

       *input* args = source source-args
         *source* = *mol* or *stl*
            *mol* arg = template-ID
               template-ID = ID of molecule template specified in a
               separate :doc:`molecule <molecule>` command, which defines a set of triangles or lines
            *stl* args = stype stlfile
               stype = numeric type assigned to all triangles in STL file
               stlfile = STL filename which defines a set of triangles

* zero or more keyword/value pairs may be appended
* keyword = *flat* or *temperature*

  .. parsed-literal::

       *flat* value = maxangle
         maxangle = maximum angle (degrees) between a pair of connected triangles/lines for a flat connection
       *temperature* value = Tsurf
         Tsurf = surface temperature (degrees Kelvin), required if model with heat is used



Examples
""""""""

.. code-block:: LAMMPS

   read_data data.tris
   fix 1 all surface/local NULL

   molecule tris surf.tri
   fix 1 all surface/local tris

   fix 1 all surface/local surf.tri.stl

Description
"""""""""""

Enable granular surfaces to be used as boundary conditions on
particles in a granular simulation.  Granular surfaces are defined as
a set of triangles (for 3d models) or a set of line segments (for 2d
models).

The :doc:`Howto granular surfaces <Howto_granular_surfaces>` doc page
gives an overview of granular surfaces of two types, *global* and
*local*, and gives guidelines for how they should be defined.

This command is used for models with *local* surfaces.  The :doc:`fix
surface/global <fix_surface_global>` command is used for models with
*global* surfaces.  As explained on the :doc:`Howto granular surfaces
<Howto_granular_surfaces>` doc page, *local* surfaces are most
appropriate when there is a large number of them.  They should
typically use surfaces (triangles/lines) whose size is no more than a
few times larger than the spherical particles used in a granular
model.  There can be as many surfaces as needed to describe one or
more physical surfaces at high resolution.

*Local* triangles or line segments are stored as triangle or line
particles and are distributed across processors in the same manner as
other particles, based on which processor's sub-domain the center of
the triangle or line segment is inside of.

*Local* triangles/lines can be defined in 3 ways.  The latter two
correspond to the 2 options listed above for the *source* argument of
the *input* keyword:

* via a data file, read by the :doc:`read_data <read_data>` command
* via a molecule file(s), read by the :doc:`molecule <molecule>` command
* via an STL file(s), read by this commmand

If triangles/lines were previously read in by the :doc:`read_data
<read_data>` command, then distributed triangles or lines already
exist.  As explained on the :doc:`Howto granular surfaces
<Howto_granular_surfaces>` doc page, these are "particles" as defined
by the :doc:`atom_style tri or line <atom_style>` command, typically
as a sub-style of the :doc:`atom_style hybrid <atom_style>` command.

Use of the *input* keyword adds new triangle/line particles to the
system.

If triangles or lines were previously read in by the :doc:`molecule
<molecule>` command, the *source* keyword of the *input* keyword is
*mol* and tts *template-ID* argument is the molecule template ID used
with the :doc:`molecule <molecule>` command.  Note that a
doc:`molecule <molecule>` command can read and assign serveral
molecule files to the same template-ID.  Each molecule file must
define triangles or lines, not atoms.  For multiple molecule files,
the set of triangle or line particles defined by this input option
will be the union of the triangles and lines from all the molecule
files.  Note that each line/triangle in a molecule file is assigned a
type and molecule ID.

An STL (stereolithography) file defines a set of triangles.  For use
with this command, the *source* argument of the *input* keyword is
*stl*.  The *stype* argument is the numeric type assigned to all the
triangles from the file.  Note that STL files do not contain types or
other flags for each triangle.  The *stlfile* argument is the name of
the STL file.  It can be in text or binary format; this command
auto-detects the format. One triangle particle is created for each
triangle in the STL file(s).  Note that STL files cannot be used for
2d simulations since they only define triangles.  Each triangle
partilce from an STL file is assigned a molecule ID = 1.

This `Wikepedia page
<https://en.wikipedia.org/wiki/STL_(file_format)>`_ describes the
format of both text and binary STL files.  Binary STL files can be
converted to ASCII for editing with the stl_bin2txt tool in the
lammps/tools directory.  Examples of text STL files are included in
the examples/gransurf directory.

Note that this command allows for pre-defined triangle/line particles
read in by the :doc:`read_data <read_data>` command as well as ultiple
uses of the *input* keyword, each with a *source* argument as either
*mol* or *stl*.  The number of triangle/line particles created by this
command will be the union of those already read by the :doc:`read_data
<read_data>` command and those specified by the *input* keywords.

Once all the distributed triangle/line particles are defined, this
command calculates their connectivity.  Two triangles are "connected"
if they have a single corner point in common or an edge in common (2
corner points).  Two line segments are "connected" if the they have an
end point in common.  More technical details on connectivity and its
significance for granular surface simulations is given on :doc:`Howto
granular surfaces <Howto_granular_surfaces>` doc page.  In brief, a
pair of connected surfaces interact with a particle which contacts
both of them simultaneously according to a set of rules which are
designed to generate physically sensible forces on the particle.

Note that there is no requirement that all the surfaces be connected
to one another.  The surfaces can represent the surface of one or more
independent objects.  Particles interact with the surface when they
are close enough to overlap (touch) one or more individual triangles
or lines.  Both sides of a triangle or line interact with particles.
Thus a surface can be infinitely thin, e.g. the blade of a mixer.  See
the :doc:`Howto granular surfaces <Howto_granular_surfaces>` doc page
for restrictions on the geometry of a collection of triangles or
lines.

The nature of individual surface/particle interactions are determined
by the :doc:`pair_coeff <pair_coeff>` command which specifies
interaction parameeters for each pair of particle types.  It is thus
important to specify different types for granular particles and
granular surfaces (triangle/line particles).  Typically a granular
simulation with local surfaces uses the :doc:`pair_style hybrid
<pair_hybrid>` command so that multiple sub-styles can be defined by
the :doc:`pair_coeff <pair_coeff>` command, one or more for
particle/particle interactions, and one or more for particle/surface
interactions.  Similar to the :doc:`fix surface/global
<fix_surface_global>` command this allows multiple particle/surface
granular interaction models to be used, based on the surface particle
types.

----------

These are the optional keywords and values.

The *flat* keyword sets a *maxangle* threshold for the angle (in
degrees) between two connected surfaces (triangles or line segments)
which will be treated as "flat" by the particle/surface interaction
models.  A flat connection means a single force will be applied to the
particle even if it is contact with both surfaces simultaneously.  See
the :doc:`Howto granular surfaces <Howto_granular_surfaces>` doc page
for more details.  The default for *maxangle* is one degree.

The *temperature* keyword is required if any of the granular models
defined by the :doc:`pair_coeff <pair_coeff>` command includes a heat
model which depends on the surface temperature.  Otherwise it is
ignored.  Its *Tsurf* value is the temperature of the surface in
degrees Kelvin.

Note that the *smax* keyword used by the :doc:`fix surface/global
<fix_surface_global>` is not used by this command, because local
triangles and lines are already particles and their type is limited by
the maximum number of particle types.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.  None of the :doc:`fix_modify <fix_modify>` options are
relevant to this fix.  No global or per-atom quantities are stored by
this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`fix surface/global <fix_surface_global>`

Default
"""""""

none

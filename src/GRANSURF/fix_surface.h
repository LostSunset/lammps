/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_FIX_SURFACE_H
#define LMP_FIX_SURFACE_H

#include "fix.h"
#include <map>
#include <tuple>

namespace LAMMPS_NS {

class FixSurface : public Fix {
 public:

  // connectivity info between lines or tris

  struct Connect2d {      // line connectivity

                          // counts, not including self
    int np1,np2;          // # of lines connected to endpts 1/2

                          // pairs of endpoint connections
    int *neigh_p1;        // indices of lines connected to endpt 1
    int *neigh_p2;        // ditto for connections to endpt 2
    int *pwhich_p1;       // which point (0,1) on other line is endpt 1
    int *pwhich_p2;       // ditto for endpt 2
    int *nside_p1;        // consistency of other line normal
    int *nside_p2;        // ditto for endpt 2
                          //   SAME_SIDE = 2 normals are on same side of surf
                          //   OPPOSITE_SIDE = opposite sides of surf
    int *aflag_p1;        // is this line + other line a FLAT,CONCAVE,CONVEX surf
    int *aflag_p2;        // ditto for endpt 2
                          //   surf = on normal side of this line
                          //   aflag = FLAT, CONCAVE, CONVEX
  };

  struct Connect3d {      // tri connectivity

                          // counts, not including self
                          // also not including edge-connected tris for corner pts
    int ne1,ne2,ne3;      // # of tris connected to edges 1,2,3
    int nc1,nc2,nc3;      // # of tris connected to corner pts 1,2,3

                          // pairs of edge connections
    int *neigh_e1;        // indices of tris connected to edge 1
    int *neigh_e2;        // ditto for connections to edge 2
    int *neigh_e3;        // ditto for connections to edge 3
    int *ewhich_e1;       // which edge (0,1,2) on other tri shares edge 1
    int *ewhich_e2;       // ditto for edge 2
    int *ewhich_e3;       // ditto for edge 3
    int *nside_e1;        // consistency of other tri normal
    int *nside_e2;        // ditto for edge 2
    int *nside_e3;        // ditto for edge 3
                          //   SAME_SIDE = 2 normals are on same side of surf
                          //   OPPOSITE_SIDE = opposite sides of surf
    int *aflag_e1;        // is this tri + other tri a FLAT,CONCAVE,CONVEX surf
    int *aflag_e2;        // ditto for edge 2
    int *aflag_e3;        // ditto for edge 3
                          //   surf = on normal side of this tri
                          //   aflag = FLAT, CONCAVE, CONVEX

                          // pairs of corner pt connections
    int *neigh_c1;        // indices of tris connected to corner pt 1
    int *neigh_c2;        // ditto for connections to corner pt 2
    int *neigh_c3;        // ditto for connections to corner pt 3
    int *cwhich_c1;       // which corner pt (0,1,2) on other tri shares corner pt 1
    int *cwhich_c2;       // ditto for corner pt 2
    int *cwhich_c3;       // ditto for corner pt 3
  };

  
  FixSurface(class LAMMPS *, int, char **);
  ~FixSurface() override;

  virtual int setmask() = 0;
  
  virtual void post_constructor() {}

 protected:

  // surfs read from molecule or STL files
  
  struct Point {
    double x[3];
  };

  struct Line {
    int mol,type;           // molID and type of the line
    int p1,p2;              // indices of points in line segment
    double norm[3];         // unit normal to line = Z x (p2-p1)
  };

  struct Tri {
    int mol,type;           // modID and type of the triangle
    int p1,p2,p3;           // indices of points in triangle
    double norm[3];         // unit normal to tri plane = (p2-p1) x (p3-p1)
  };

  // methods common to both global and local surfs

  void extract_from_molecule(char *, std::map<std::tuple<double,double,double>,int> *,
                             int &, int &, Point *&,
                             int &, Line *&, int &, Tri *&);
  void extract_from_stlfile(char *, int, std::map<std::tuple<double,double,double>,int> *,
                            int &, int &, Point *&, int &, Tri *&);

};

}    // namespace LAMMPS_NS

#endif

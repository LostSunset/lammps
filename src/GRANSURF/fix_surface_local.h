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

#ifdef FIX_CLASS

FixStyle(surface/local,FixSurfaceLocal)

#else

#ifndef LMP_FIX_SURFACE_LOCAL_H
#define LMP_FIX_SURFACE_LOCAL_H

#include <stdio.h>
#include "fix_surface.h"
#include <map>
#include "my_pool_chunk.h"

namespace LAMMPS_NS {

class FixSurfaceLocal : public FixSurface {
 public:

  // 2d/3d connectivity

  Connect2d *connect2d;         // 2d connection info
  Connect3d *connect3d;         // 3d connection info
  int nmax_connect;             // allocated size of connect2d/3d

  FixSurfaceLocal(class LAMMPS *, int, char **);
  virtual ~FixSurfaceLocal();
  int setmask();
  void post_constructor();
  
  void setup_pre_neighbor() override;
  void pre_neighbor() override;

  void grow_arrays(int) override;
  void grow_connect();
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  void clear_bonus() override;

  int pack_border(int, int *, double *) override;
  int unpack_border(int, int, double *) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  double memory_usage() override;

 private:
  int dimension,mode;
  double flatthresh;
  double Twall;
  int tvar;
  char *tstr;

  int ninput;
  int *input_modes,*input_stypes;
  char **input_sources;
  
  int flag_complete;    // whether one-time connectivity info has been calculated

  class AtomVecLine *avec_line;
  class AtomVecTri *avec_tri;

  // memory allocation for tagint and int vectors in Connect 2d/3d

  MyPoolChunk<tagint> *tcp;     // allocator for 2d/3d connectivity vecs

  struct Pool2d {
    int neigh_p1,neigh_p2;    // pool indices of neigh_p12 chunks
    int pwhich_p1,pwhich_p2;  // pool indices of pwhich_p12 chunks
    int nside_p1,nside_p2;    // pool indices of nside_p12 chunks
    int aflag_p1,aflag_p2;    // pool indices of aflag_p12 chunks
  };

  struct Pool3d {
    int neigh_e1,neigh_e2,neigh_e3;     // pool indices of neigh_e123 chunks
    int ewhich_e1,ewhich_e2,ewhich_e3;  // pool indices of neigh_e123 chunks
    int nside_e1,nside_e2,nside_e3;     // pool indices of neigh_e123 chunks
    int aflag_e1,aflag_e2,aflag_e3;     // pool indices of neigh_e123 chunks
    int neigh_c1,neigh_c2,neigh_c3;     // pool indices of neigh_c123 chunks
    int cwhich_c1,cwhich_c2,cwhich_c3;  // pool indices of neigh_c123 chunks
  };

  Pool2d *pool2d;               // pool indices of connect2d vectors
  Pool3d *pool3d;               // pool indices of connect3d vectors

  int *atom2connect;       // per-atom index into connect 2d/3d vecs, -1 if none
  int *connect2atom;       // per-connect index into atoms

  // ragged 2d arrays for 2d connectivity for global case

  tagint **neigh_p1,**neigh_p2;

  // ragged 2d arrays for 3d connectivity for global case

  tagint **neigh_e1,**neigh_e2,**neigh_e3;
  tagint **neigh_c1,**neigh_c2,**neigh_c3;

  // size of local/ghost connection info vectors

  int nlocal_connect,nghost_connect;

  // data structs for calculating global connectivity of line/tri particles
  // only used by Rvous comm during setup

  struct InRvous {            // datum for each bin of each endpt in each line
    int proc;                 // owning proc
    int ibin;                 // bin assignment
    int ilocal;               // index of line particle
    int ipoint;               // 0/1 for either endpt
    tagint atomID;            // ID of line particle
    double x[3];              // coords of endpt
  };

  struct OutRvous {           // datum returned to owning proc
    int ilocal;               // index of line particle
    int ipoint;               // 0/1 for either endpt
    tagint atomID;            // ID of other line particle with matching endpt
  };

  double bboxlo[3],bboxhi[3];    // bounding box around all lines/tris

  // data structs for extracting surfs from molecule template or STL file
  // these are global data structs for all surfs, including connectivity
  // only used during setup, then deleted

  Point *points;              // global list of points
  Line *lines;                // global list of lines
  Tri *tris;                  // global list of tris
  int npoints,nlines,ntris;   // count of each
  int maxpoints;
  
  Connect2d *connect2dall;    // global connectivity info
  Connect3d *connect3dall;

  int **plines;               // ragged 2d array for global line end pt lists
  int **etris;                // ragged 2d array for global tri edge lists
  int **ctris;                // ragged 2d array for global tri corner pt lists

  struct OnePt {               // one end/corner point of iline/itri in a bin
    int iatom;                 // local index of the line/tri in atom list
    int iconnect;              // local index of the line/tri in connect list
    int ptwhich;               // 1/2 for two end pts of line
                               // 1/2/3 for three corner pts of tri
  };

  OnePt *pts;                      // list of pts in all bins
  int *bincount;                   // # of pts per bin
  int *binfirst;                   // index into pts of first pt in bin
  int nbins,nbinx,nbiny,nbinz;     // # of total bins and in each dim
  double invbinx,invbiny,invbinz;  // inverse of bin size
  double binlo[3],binhi[3];        // bounds of region that is binned

  double eps,epssq;            // distance tolerance for matching 2 points
                               // for different lines/tris to be connected

  int nmatch;                  // # of line connections
  int nmatch1,nmatch2;         // # of tri connections
  int errormatch;              // # of errors with line connectivity
  int errormatch1,errormatch2; // # of errors with tri connectivity

  int vecflag;            // 0/1 whether tri matching should also
                          // store variable-length vecs of corner connections

  // static variable for Rvous communication callback to access class data
  // callback functions for Rvous communication

  static FixSurfaceLocal *fptr;
  static void linematch(int, char *);
  static void trimatch(int, char *);
  static int point_match(int, char *, int &, int *&, char *&, void *);

  // private methods

  // local connectivity build from line or triangle particles

  void connectivity2d_local();
  void connectivity3d_local();
  void calculate_endpts(double **);
  int overlap_bins_2d(double *, double, int *);
  void calculate_corners(double **);
  int overlap_bins_3d(double *, double, int *);

  // global connectivity build from molecule or STL files

  int check_exist();
  void assign2d();
  void assign3d();

  // finish the connectivity build for both local and global cases

  void connectivity2d_complete();
  void connectivity3d_complete();
  int same_point(double *, double *);
  void epsilon_calculate();
};

}

#endif
#endif

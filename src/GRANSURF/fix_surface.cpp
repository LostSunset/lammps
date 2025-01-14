/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_surface.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "molecule.h"
#include "stl_reader.h"

using namespace LAMMPS_NS;

#define DELTA 128

/* ---------------------------------------------------------------------- */

FixSurface::FixSurface(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

FixSurface::~FixSurface() {}

/* ----------------------------------------------------------------------
   extract lines or tris from a molecule template ID for one or more molecules
   concatenate into single list of points and lines or tris
   identify unique points using hash
------------------------------------------------------------------------- */

void FixSurface::extract_from_molecule(char *molID,
                                       std::map<std::tuple<double,double,double>,int> *hash,
                                       int &npoints, int &maxpoints,
                                       Point *&points, int &nlines, Line *&lines,
                                       int &ntris, Tri *&tris)
{
  int dimension = domain->dimension;
  
  int imol = atom->find_molecule(molID);
  if (imol == -1)
    error->all(FLERR,"Molecule template ID for fix surface does not exist");

  // loop over one or more molecules in molID

  Molecule **onemols = &atom->molecules[imol];
  int nmol = onemols[0]->nset;

  for (int m = 0; m < nmol; m++) {
    if (dimension == 2)
      if (onemols[m]->lineflag == 0)
        error->all(FLERR,"Fix surface molecule must have lines");
    if (dimension == 3)
      if (onemols[m]->triflag == 0)
        error->all(FLERR,"Fix surface molecule must have triangles");

    int nl = onemols[m]->nlines;
    int nt = onemols[m]->ntris;

    nlines += nl;
    ntris += nt;
    lines = (Line *) memory->srealloc(lines,nlines*sizeof(Line),
                                      "surface:lines");
    tris = (Tri *) memory->srealloc(tris,ntris*sizeof(Tri),
                                    "surface:tris");

    // offset line/tri index lists by previous npoints
    // pi,p2,p3 are C-style indices into points vector

    if (dimension == 2) {
      int *molline = onemols[m]->molline;
      int *typeline = onemols[m]->typeline;
      double **epts = onemols[m]->lines;
      int iline = nlines - nl;

      for (int i = 0; i < nl; i++) {
        lines[iline].mol = molline[i];
        lines[iline].type = typeline[i];

        auto key = std::make_tuple(epts[i][0],epts[i][1],0.0);
        if (hash->find(key) == hash->end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface:points");
          }
          (*hash)[key] = npoints;
          points[npoints].x[0] = epts[i][0];
          points[npoints].x[1] = epts[i][1];
          points[npoints].x[2] = 0.0;
          lines[iline].p1 = npoints;
          npoints++;
        } else lines[iline].p1 = (*hash)[key];

        key = std::make_tuple(epts[i][2],epts[i][3],0.0);
        if (hash->find(key) == hash->end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface:points");
          }
          (*hash)[key] = npoints;
          points[npoints].x[0] = epts[i][2];
          points[npoints].x[1] = epts[i][3];
          points[npoints].x[2] = 0.0;
          lines[iline].p2 = npoints;
          npoints++;
        } else lines[iline].p2 = (*hash)[key];

        iline++;
      }
    }

    if (dimension == 3) {
      int *moltri = onemols[m]->moltri;
      int *typetri = onemols[m]->typetri;
      double **cpts = onemols[m]->tris;
      int itri = ntris - nt;

      for (int i = 0; i < nt; i++) {
        tris[itri].mol = moltri[i];
        tris[itri].type = typetri[i];

        auto key = std::make_tuple(cpts[i][0],cpts[i][1],cpts[i][2]);
        if (hash->find(key) == hash->end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface:points");
          }
          (*hash)[key] = npoints;
          points[npoints].x[0] = cpts[i][0];
          points[npoints].x[1] = cpts[i][1];
          points[npoints].x[2] = cpts[i][2];
          tris[itri].p1 = npoints;
          npoints++;
        } else tris[itri].p1 = (*hash)[key];

        key = std::make_tuple(cpts[i][3],cpts[i][4],cpts[i][5]);
        if (hash->find(key) == hash->end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface:points");
          }
          (*hash)[key] = npoints;
          points[npoints].x[0] = cpts[i][3];
          points[npoints].x[1] = cpts[i][4];
          points[npoints].x[2] = cpts[i][5];
          tris[itri].p2 = npoints;
          npoints++;
        } else tris[itri].p2 = (*hash)[key];

        key = std::make_tuple(cpts[i][6],cpts[i][7],cpts[i][8]);
        if (hash->find(key) == hash->end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface:points");
          }
          (*hash)[key] = npoints;
          points[npoints].x[0] = cpts[i][6];
          points[npoints].x[1] = cpts[i][7];
          points[npoints].x[2] = cpts[i][8];
          tris[itri].p3 = npoints;
          npoints++;
        } else tris[itri].p3 = (*hash)[key];

        itri++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   extract triangles from an STL file, can be text or binary
   concatenate into single list of points and tris
   identify unique points using hash
------------------------------------------------------------------------- */

void FixSurface::extract_from_stlfile(char *filename, int stype,
                                      std::map<std::tuple<double,double,double>,int> *hash,
                                      int &npoints, int &maxpoints,
                                      Point *&points, int &ntris, Tri *&tris)
{
  if (domain->dimension == 2)
    error->all(FLERR,"Fix surface cannot use an STL file for 2d simulations");

  // read tris from STL file
  // stltris = tri coords internal to STL reader

  STLReader *stl = new STLReader(lmp);
  double **stltris;
  ntris = stl->read_file(filename,stltris);

  tris = (Tri *) memory->srealloc(tris,ntris*sizeof(Tri),"surface:tris");

  // loop over STL tris
  // populate points and tris data structs
  // for each tri: set molID = 1 and type = stype

  for (int itri = 0; itri < ntris; itri++) {
    tris[itri].mol = 1;
    tris[itri].type = stype;

    auto key = std::make_tuple(stltris[itri][0],stltris[itri][1],stltris[itri][2]);
    if (hash->find(key) == hash->end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface:points");
      }
      (*hash)[key] = npoints;
      points[npoints].x[0] = stltris[itri][0];
      points[npoints].x[1] = stltris[itri][1];
      points[npoints].x[2] = stltris[itri][2];
      tris[itri].p1 = npoints;
      npoints++;
    } else tris[itri].p1 = (*hash)[key];

    key = std::make_tuple(stltris[itri][3],stltris[itri][4],stltris[itri][5]);
    if (hash->find(key) == hash->end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface:points");
      }
      (*hash)[key] = npoints;
      points[npoints].x[0] = stltris[itri][3];
      points[npoints].x[1] = stltris[itri][4];
      points[npoints].x[2] = stltris[itri][5];
      tris[itri].p2 = npoints;
      npoints++;
    } else tris[itri].p2 = (*hash)[key];

    key = std::make_tuple(stltris[itri][6],stltris[itri][7],stltris[itri][8]);
    if (hash->find(key) == hash->end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface:points");
      }
      (*hash)[key] = npoints;
      points[npoints].x[0] = stltris[itri][6];
      points[npoints].x[1] = stltris[itri][7];
      points[npoints].x[2] = stltris[itri][8];
      tris[itri].p3 = npoints;
      npoints++;
    } else tris[itri].p3 = (*hash)[key];
  }

  // delete STL reader

  delete stl;
}


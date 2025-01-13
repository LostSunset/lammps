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

#ifndef RUN_WHAM_H
#define RUN_WHAM_H

#include <QDialog>

class RunWHAM : public QDialog {
    Q_OBJECT

public:
    explicit RunWHAM(QWidget *parent = nullptr);
    ~RunWHAM() = default;

private slots:
    void accept() override;
};

#endif

// Local Variables:
// c-basic-offset: 4
// End:

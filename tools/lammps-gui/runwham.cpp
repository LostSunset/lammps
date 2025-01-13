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

#include "runwham.h"

#include <QDialogButtonBox>
#include <QGridLayout>
#include <QIcon>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSizePolicy>

RunWHAM::RunWHAM(QWidget *parent) : QDialog(parent)
{
    int irow     = 0;
    auto *layout = new QGridLayout;
    auto *top    = new QLabel("Run WHAM Utility:");
    layout->addWidget(top, irow++, 0, 1, 2, Qt::AlignHCenter);

    auto *label = new QLabel("Units: ");
    layout->addWidget(label, irow, 0, 1, 1, Qt::AlignLeft);
    auto *field = new QLineEdit("real");
    layout->addWidget(field, irow++, 1, 1, 1);
    label = new QLabel("Mininum value: ");
    layout->addWidget(label, irow, 0, 1, 1, Qt::AlignLeft);
    field = new QLineEdit("");
    layout->addWidget(field, irow++, 1, 1, 1);
    label = new QLabel("Maximum value: ");
    layout->addWidget(label, irow, 0, 1, 1, Qt::AlignLeft);
    field = new QLineEdit("");
    layout->addWidget(field, irow++, 1, 1, 1);
    label = new QLabel("Number of bins: ");
    layout->addWidget(label, irow, 0, 1, 1, Qt::AlignLeft);
    field = new QLineEdit("50");
    layout->addWidget(field, irow++, 1, 1, 1);
    label = new QLabel("Tolerance: ");
    layout->addWidget(label, irow, 0, 1, 1, Qt::AlignLeft);
    field = new QLineEdit("1.0e-8");
    layout->addWidget(field, irow++, 1, 1, 1);
    label = new QLabel("Temperature: ");
    layout->addWidget(label, irow, 0, 1, 1, Qt::AlignLeft);
    field = new QLineEdit("300.0");
    layout->addWidget(field, irow++, 1, 1, 1);
    label = new QLabel("Metadata file: ");
    layout->addWidget(label, irow, 0, 1, 1, Qt::AlignLeft);
    field = new QLineEdit("");
    layout->addWidget(field, irow++, 1, 1, 1);
    label = new QLabel("Output file: ");
    layout->addWidget(label, irow, 0, 1, 1, Qt::AlignLeft);
    field = new QLineEdit("");
    layout->addWidget(field, irow++, 1, 1, 1);

    auto *buttonBox = new QDialogButtonBox(QDialogButtonBox::Cancel);
    auto *add       = new QPushButton("&Run WHAM");
    add->setObjectName("run_wham");
    buttonBox->addButton(add, QDialogButtonBox::ActionRole);
    connect(add, &QPushButton::released, this, &RunWHAM::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &QDialog::reject);

    layout->addWidget(buttonBox, irow++, 0, 1, 2, Qt::AlignHCenter);
    setLayout(layout);
    setWindowIcon(QIcon(":/icons/lammps-icon-128x128.png"));
    setWindowTitle("LAMMPS-GUI - Run WHAM Utility");
    resize(300, 200);
}

void RunWHAM::accept()
{
    // assemble command line and run WHAM
    fprintf(stderr, "Running WHAM\n");
    QDialog::accept();
}

// Local Variables:
// c-basic-offset: 4
// End:

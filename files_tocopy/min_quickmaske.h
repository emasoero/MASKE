/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef MINIMIZE_CLASS

MinimizeStyle(quickmaske,MinQuickmaske)

#else

#ifndef LMP_MIN_QUICKMASKE_H
#define LMP_MIN_QUICKMASKE_H

#include "min.h"

namespace LAMMPS_NS {

class MinQuickmaske : public Min {
 public:
  MinQuickmaske(class LAMMPS *);
  ~MinQuickmaske() {}
  void init();
  void setup_style();
  void reset_vectors();
  int iterate(int);

 private:
  double dt;
  bigint last_negative;
};

}

#endif
#endif

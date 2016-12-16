/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <math.h>
#include "fix_cwall_lj126.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCWallLJ126::FixCWallLJ126(LAMMPS *lmp, int narg, char **arg) :
  FixWall(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixCWallLJ126::precompute(int m)
{
  coeff1[m] = 48.0 * epsilon[m] * pow(sigma[m],12.0);
  coeff2[m] = 24.0 * epsilon[m] * pow(sigma[m],6.0);
  coeff3[m] = 4.0 * epsilon[m] * pow(sigma[m],12.0);
  coeff4[m] = 4.0 * epsilon[m] * pow(sigma[m],6.0);

  double r2inv = 1.0/(cutoff[m]*cutoff[m]);
  double r6inv = r2inv*r2inv*r2inv;
  offset[m] = r6inv*(coeff3[m]*r6inv - coeff4[m]);
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   fixwall for x and z direction as circle shape.
   coord: circle radius
   which is initialied as z direction.
   One using this program should change the center in x and z direction by hand. such as the line
   67 and line 68. radius cand be changed in input file.
------------------------------------------------------------------------- */


void FixCWallLJ126::wall_particle(int m, int which, double coord)
{
  double delta,rinv,r2inv,r6inv,fwall;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  double centerx = 27.00;
  double centerz = 72.00;
  double deltax,deltaz, deltaxz;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      //if (side < 0) delta = x[i][dim] - coord;
      //else delta = coord - x[i][dim];
	  deltaz = x[i][dim]-centerz;
	  deltax = x[i][dim-2]-centerx;
	  deltaxz = sqrt(deltax*deltax+deltaz*deltaz);
	  delta = coord - deltaxz;
	  
      if (delta >= cutoff[m]) continue;
      if (delta <= 0.0) {
        onflag = 1;
        continue;
      }
      rinv = 1.0/delta;
      r2inv = rinv*rinv;
      r6inv = r2inv*r2inv*r2inv;
      fwall = side * r6inv*(coeff1[m]*r6inv - coeff2[m]) * rinv;
	  
      f[i][dim] -= fwall*deltaz/deltaxz;
	  f[i][dim-2] -= fwall*deltax/deltaxz;
	  
      ewall[0] += r6inv*(coeff3[m]*r6inv - coeff4[m]) - offset[m];
      ewall[m+1] += fwall;
    }

  if (onflag) error->one(FLERR,"Particle on or inside fix cwall surface");
}

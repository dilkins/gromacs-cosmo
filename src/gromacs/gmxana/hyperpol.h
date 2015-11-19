/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef _hyperpol_h
#define _hyperpol_h


#include "typedefs.h"
#include "vec.h"
#include "network.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void mol_unitvectors(const rvec xv1, const rvec xv2, const rvec xv3, rvec u1, rvec u2, rvec u3);

real rotate_beta(real invnormx, real invnormz, const rvec xv2, const rvec xv3, const rvec pout, const rvec pin1, const rvec pin2, real *beta_m);

extern void rotate_beta_theta(real invnormx, real invnormz, const rvec xv2, const rvec xv3, const rvec pin1, const rvec pin2, real *beta_m, real *beta_2, real *beta_1);

void rotate_wave_vec(const rvec wave_vec, const int rot_label, rvec rot_vec); 

void induced_second_order_dipole(real invnormx, real invnormz, const rvec xv2, const rvec xv3, const rvec pout, const rvec pin, real ***betamol, real *mu_ind);


extern void dipole_atom2mol(int *n, int *index, t_block *mols);

extern void double_sum_fade(t_pbc *pbc,real *beta,rvec *x,int n,real rmax2,rvec *qvec,int nbinq,real fade,real inv_width,real *temp_method);

extern void double_sum(t_pbc *pbc,real *beta,rvec *x,int n,real rmax2,rvec *qvec,int nbinq,real *temp_method);

extern void calc_beta(t_pbc *pbc,int natoms,int n,int *ind,rvec *x,real norm_x,real norm_z,rvec pol_out,rvec pol_in1,rvec pol_in2,real *bete_mol_1d,real *beta);

extern void Projected_Scattering_Amplitude(const int nf, const int nt, const int nga, const int nq,
     const rvec xi, const rvec xv2, const rvec xv3,  rvec **arr_scatt_wave_vec,
     rvec ***pout_theta_gamma, rvec ***pin_theta_gamma,
     real *********Onsite_term, real *******Cos_scatt_ampl, real *******Sin_scatt_ampl, 
     real **********Onsite_term_addr, real ********Cos_scatt_ampl_addr, real ********Sin_scatt_ampl_addr); 



extern void Allocate_scattering_amplitude(const int nf, const int nt, const int nga, const int nq, real **********Incoh_term_addr, real ********Cos_scatt_ampl_addr, real ********Sin_scatt_ampl_addr);

extern void Free_scattering_amplitude(const int nf, const int nt, const int nga, real *********Incoh_term, real *******Cos_scatt_ampl, real *******Sin_scatt_ampl);

extern void Allocate_Scattering_Intensity(const int nf, const int nt,  const int nq,
                               real **********tot_tensor_squared_addr, real  **********incoh_tensor_squared_addr);

extern void Scattering_Intensity_t(const int nf, const int nt, const int nga, const int nq,
                               const real inversesize, real *********Incoh_term, real *******Cos_scatt_ampl, real *******Sin_scatt_ampl,
                               real  *********tot_tensor_squared,      real *********incoh_tensor_squared,
                               real **********tot_tensor_squared_addr, real **********incoh_tensor_squared_addr);

extern void Print_tensors(const int nt, const int nframes, const real invgamma, real *theta_vec, real *********tot_tensor_squared,
                           real *********incoh_tensor_squared, const char *fnTENSOR,  const char *fnINCTENSOR ,char **grpname ,const output_env_t oenv);

extern void Print_scattering_pattern(const int nt,  const int nframes, const real invgamma,
                               real bmzxx, real bmzzz, real bmzyy, real *theta_vec,
                               real *********tot_tensor_squared, real *********incoh_tensor_squared,
                               const char *fnTHETA, char **grpname , const output_env_t oenv);

extern void Free_Scattering_Intensity(const int nf, const int nt, 
                               real *********tot_tensor_squared, real *********incoh_tensor_squared); 
           


#ifdef __cplusplus
}
#endif

#endif

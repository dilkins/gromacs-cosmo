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
#include "gromacs/random/random.h"
#include "coulomb.h"
#include "gromacs/math/gmxcomplex.h"
#include "types/inputrec.h"
#include "complex.h"


#ifdef __cplusplus
extern "C" {
#endif

extern void mol_unitvectors(const rvec xv1, const rvec xv2, const rvec xv3, rvec u1, rvec u2, rvec u3);

real rotate_beta(const rvec xv2, const rvec xv3, const rvec pout, const rvec pin1, const rvec pin2, real *beta_m);

extern void rotate_beta_theta( const rvec xv2, const rvec xv3, const rvec pin1, const rvec pin2, real *beta_m, real *beta_2, real *beta_1);

void rotate_wave_vec(const rvec wave_vec, const int rot_label, rvec rot_vec); 

void induced_second_order_dipole( const rvec xv2, const rvec xv3, const rvec pout, const rvec pin, real ***betamol, real *mu_ind);

void induced_second_order_dipole_fluct_beta( const rvec xv2, const rvec xv3, const rvec pout, const rvec pin, real ***betamol, real *mu_ind);


void beta_gaussian_noise(real ***beta_const_mean, real ***beta_const_dev, gmx_rng_t rng, real ****beta_gauss);


extern void dipole_atom2mol(int *n, int *index, t_block *mols);

extern void double_sum_fade(t_pbc *pbc,real *beta,rvec *x,int n,real rmax2,rvec *qvec,int nbinq,real fade,real inv_width,real *temp_method);

extern void double_sum(t_pbc *pbc,real *beta,rvec *x,int n,real rmax2,rvec *qvec,int nbinq,real *temp_method);

extern void Projected_Scattering_Amplitude(const int nf, const int nt, const int nga, const int nq,
     const rvec xi, const rvec xv2, const rvec xv3,  rvec **arr_scatt_wave_vec,
     rvec ***pout_theta_gamma, rvec ***pin_theta_gamma,
     real **********Onsite_term_addr, real ********Cos_scatt_ampl_addr, real ********Sin_scatt_ampl_addr); 

extern void Allocate_scattering_amplitude(const int nf, const int nt, const int nga, const int nq, real **********Incoh_term_addr, real ********Cos_scatt_ampl_addr, real ********Sin_scatt_ampl_addr);

extern void Free_scattering_amplitude(const int nf, const int nt, const int nga, real *********Incoh_term, real *******Cos_scatt_ampl, real *******Sin_scatt_ampl);

extern void Allocate_Scattering_Intensity(const int nf, const int nt,  const int nq, const int nthreads,
                               real **********tot_tensor_squared_addr, real  **********incoh_tensor_squared_addr);

extern void Scattering_Intensity_t(int tid, const int nf, real time, const int nt, const int nga, const int nq, 
                                   const real inversenorm, const real inversesize, const real *theta_vec,
                                   real *********Incoh_term, real *******Cos_scatt_ampl, real *******Sin_scatt_ampl,
                                   real  *********tot_tensor_squared,      real *********incoh_tensor_squared,
                                   real **********tot_tensor_squared_addr, real **********incoh_tensor_squared_addr,
                                   const char *fnTIMEEVOLTENSOR, FILE  *fname);

extern void Print_tensors( const int nfaces, const int nt, const int nframes, const real invgamma, real *theta_vec, int nbinq, rvec **arr_qvec_faces, real *********tot_tensor_squared,
                           real *********incoh_tensor_squared, const char *fnTENSOR,  const char *fnINCTENSOR , const char *fnQSWIPE, char **grpname ,const output_env_t oenv);

extern void Print_scattering_pattern(const int nt,  const int nframes, const real invgamma,
                               real bmzxx, real bmzzz, real bmzyy, real *theta_vec,
                               real *********tot_tensor_squared, real *********incoh_tensor_squared,
                               const char *fnTHETA, char **grpname , const output_env_t oenv);

extern void Free_Scattering_Intensity(const int nf, const int nt, 
                               real *********tot_tensor_squared, real *********incoh_tensor_squared); 
           

typedef struct {
  real beta_gas[27];
  real D[27][212];
} t_Map;

typedef struct {
  real   **krrinput; // input quantity for krr that is read from file and it is krrinput[ndataset][gridpoints] 
  rvec    *grid; // input local grid geometry centred on a given atom
  rvec    *translgrid; //translated and rotated grid for each molecule
  rvec    *rotgrid; //rotated grid for each molecule
  real ****coeff;// coefficients for molecular beta obtained from kernel. pointer of dimension ndataset*3*3*3 in case of krr or gridpoints*3*3*3 in case of scalar fit
  real    *meanquant; //average of feature vector on each grid point
  int      gridpoints; // number of points in the grid
  int      gridcenter; //index of the grid that is taken as reference of the potential or field;
  real     kerndev;
  int      ndataset; //number of points used to train the kernel
  rvec     rspace_grid; // a point in the global grid that fills all unit cell
  real  ***quantity_on_grid; //scalar quantity computed on the global grid, i.e. density
  real  ***quantity_on_grid_x,***quantity_on_grid_y,***quantity_on_grid_z ; // vector quantities computed on global grid, i.e. electric field
  real    *interp_quant_grid; // interpolated quantity from the global onto local grid (i.e. density)
  rvec    *vec_interp_quant_grid; // interpolated vector quantity from the global onto local grid (i.e. electric field)
  real    *weights; // weigh the density based on the distance from central atom
  real    *selfterm;  // self term in the density
} t_Kern;

typedef struct {
  int   atom[1];
  real  q[1];
} t_Ion;


//read input files for the calculation of beta using an electric field map
void readMap(const char *fnMAP, t_Map *Map);
//compute the electrostatic field using a map the form of the field is r^-2
void calc_efield_map(t_pbc *pbc,t_topology *top, t_block *mols, t_Ion *Cation, t_Ion *Anion, int  *molindex[],
                 const int gind , const int isize0, const int ncations, const int nanions,
                 rvec *x, const int imol,  rvec xvec, rvec yvec,  rvec zvec, real electrostatic_cutoff2, real ***field_addr);
//compute beta for a given molecule using the electric field map
void calc_beta_efield_map(t_Map *Map, real **Efield, real ***betamol, real ****betamol_addr);

//compute the component of the second order dipole projected onto the polarization vectors pin and pout
extern void induced_second_order_fluct_dipole( matrix cosdirmat,
                                       const rvec pout, const rvec pin,
                                       real ***betamol,
                                       real *mu_ind);

//compute the component of the second order dipole projected onto the polarization vectors pin and pout
extern void  induced_second_order_fluct_dipole_fluct_beta( matrix cosdirmat,
                                       const rvec pout, const rvec pin,
                                       real ***betamol,
                                       real *mu_ind);


// read input files for the calculation of beta using either a scalar kernel or krr
void readKern(const char *fnCOEFF, const char *fnGRD, const char *fnINPKRR, int kern_order, real std_dev, real kappa, rvec *xref, gmx_bool bAtomCenter, real **betamean, t_Kern *Kern);
//compute molecular beta using kernel ridge regression
void calc_beta_krr(t_Kern *Krr, t_pbc *pbc, t_topology *top, t_block *mols, int  *molindex[], const int gind , const int isize0, rvec *x, rvec xcm_transl, const int imol, real rmax2, real electrostatic_cutoff, real ****betamol);



//compute molecular beta using a scalar kernel
void calc_beta_skern( t_Kern *SKern_rho_O, t_Kern *SKern_rho_H, t_Kern *SKern_E, int kern_order, real *betamean, real ****betamol);
//switching function to bring a function smoothly from a given value to zero
void switch_fn(real r_dist, real electrostatic_cutoff, real rmax, real inv_width, real *sw_coeff);

//function to compute <beta(0)*beta(r)>, where beta is the projected beta component, resulting from pin and pout
extern void  calc_beta_corr( t_pbc *pbc, t_block *mols, int  *molindex[],
                 const int gind , const int isize0, int nbin, const real rmax2,  real invhbinw,
                 rvec *x, real *mu_ind_mols, real **beta_corr);
void calc_ft_beta_corr( t_pbc *pbc, t_block *mols, int  *molindex[],
                 const int gind , const int isize0,  int nbinq, rvec **arr_qvec_faces,  real rmax2,  real invhbinw,
                 rvec *x, real *mu_ind_mols, real **ft_beta_corr);


//compute direction cosine matrix
void calc_cosdirmat(const char *fnREFMOL, t_topology *top, int molsize,  int ind0, rvec *xref, rvec *xmol, matrix *cosdirmat, matrix *invcosdirmat, rvec *xvec, rvec *yvec, rvec *zvec);


void rotate_local_grid_points(t_Kern *SKern_rho_O, t_Kern *SKern_rho_H, t_Kern *SKern_E, int ePBC, matrix box, matrix cosdirmat, rvec xi);

//build global grid in real space used for the calculation of the density and electrostatic potential
//for the specific case of the scalar kernel
void initialize_free_quantities_on_grid(t_Kern *Kern, t_inputrec *ir,  rvec *grid_spacing, rvec *grid_invspacing, gmx_bool bEFIELD, gmx_bool bALLOC, int **gridsize);
//compute the density using the scalar kernel
void calc_dens_on_grid(t_Kern *Kern, t_inputrec *ir, t_pbc *pbc, 
                       t_block *mols, int  *molindex[], int atom_id0, int nspecies ,int isize0, rvec *x,
                       real std_dev_dens, real inv_std_dev_dens, rvec grid_invspacing, int *gridsize, rvec grid_spacing);

//function to interpolate linearly the density from points in the global real space grid
//to the points of the local grid centred on the molecule, whose position is xi
void trilinear_interpolation_kern(t_Kern *Kern, t_inputrec *ir, t_pbc *pbc, rvec xi, rvec grid_invspacing, rvec grid_spacing);

void vec_trilinear_interpolation_kern(t_Kern *Kern, t_inputrec *ir, t_pbc *pbc, matrix invcosdirmat, rvec xi, rvec grid_invspacing, rvec grid_spacing, rvec Emean);


void bspline_efield(t_Kern *Kern, t_inputrec *ir, t_pbc *pbc, matrix invcosdirmat, rvec xi, rvec grid_invspacing, int interp_order, rvec grid_spacing, rvec Emean);

//read in the reference molecule, used to compute the direction cosine matrix
void read_reference_mol(const char *fnREFMOL, rvec **xref);

//check if there are ions in the configuration
int check_ion(t_topology *top, char *name);
//get their indexes, charge and atom names
void identifyIon(t_topology *top,t_Ion *Ion, char *name);

//the following functions are used to compute long range of electrostatic potentials with SPME

//simple factorial function used in splines
unsigned long factorial(unsigned long f);
//simple spline function
real Bspline(real x,int n);

long long combi(int n,int k);

void do_fft(real ***rmatr,t_complex ***kmatr,int *dims, real multiplication_factor, int fwbck);
void setup_ewald_pair_potential(int *grid, int interp_order, int kmax,t_complex ****FT_pair_pot,real invkappa2);
void calculate_spme_efield(t_Kern *Kern, t_inputrec *ir, t_topology *top,  
                          matrix box, real invvol, t_block *mols, int  *molindex[],
                         int *chged_atom_indexes, int n_chged_atoms, int *grid, rvec grid_spacing,
                         int interp_order, rvec *x, int isize0,
                         t_complex ***FT_pair_pot, rvec *Emean, real eps);

int index_wrap(int idx, int wrap);

#ifdef __cplusplus
}
#endif

#endif

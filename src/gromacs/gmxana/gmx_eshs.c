 /*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#include <fftw3.h>

#include "sysstuff.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "xvgr.h"
#include "gromacs/fileio/futil.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/random/random.h"
#include "gromacs/fileio/trxio.h"
#include "physics.h"
#include "index.h"
#include "gromacs/utility/smalloc.h"
#include "calcgrid.h"
#include "nrnb.h"
#include "coulomb.h"
#include "gstat.h"
#include "gromacs/fileio/matio.h"
#include "gmx_ana.h"
#include "hyperpol.h"
#include "names.h"
#include "coulomb.h"
#include "gromacs/math/gmxcomplex.h"
#include "gromacs/fft/fft.h"
#include "gromacs/fft/parallel_3dfft.h"
#include "gromacs/math/do_fit.h"
#include "mtop_util.h"
#include "typedefs.h"
#include "force.h"

#include "gromacs/legacyheaders/gmx_fatal.h"

#define nm2Bohr 10.0/0.529177
#define ANG2NM 0.1
#define AU2VNM 5.14220652e2

static void do_eshs(t_topology *top,  const char *fnTRX,
                   const char *fnSFACT, const char *fnTHETA, const real angle_corr,
                   const char *fnVCOEFF, const char *fnVGRD, const char *fnVINP,
                   const char *fnRGRDO, const char *fnCOEFFO,
                   const char *fnRGRDH, const char *fnCOEFFH, const char *fnMAP,
                   const char *fnBETACORR, const char *fnFTBETACORR, const char *fnREFMOL,
                   const char *method, const char *kern,
                   gmx_bool bIONS, char *catname, char *anname, gmx_bool bPBC, 
                   int qbin, int nbinq, int kern_order, real std_dev_dens, real fspacing,
                   real binwidth, int nbintheta, int nbingamma, real pin_angle, real pout_angle,
                   real cutoff_field, real maxelcut, real kappa, int interp_order, int kmax, real kernstd,
                   int *isize, int  *molindex[], char **grpname, int ng,
                   const output_env_t oenv, real eps, real *sigma_vals, real ecorrcut)
{
    FILE          *fp, *fpn;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, nanions, ncations, i, j, k, qq, n, c, tt, rr, nframes, nfaces, gr_ind, nbin, aa, bb, cc;
    real         **s_method, **s_method_coh, **s_method_incoh, *temp_method, ****s_method_t, ****s_method_coh_t, ****s_method_incoh_t, ***mu_sq_t, ***coh_temp;
    real           qnorm, maxq, incoh_temp = 0.0, tot_temp = 0.0, gamma = 0.0 ,theta0 = 5.0, check_pol;
    real          *cos_t, *sin_t, ****cos_tq, ****sin_tq,   mu_ind =0.0, mu_sq =0.0, mod_f ;
    real         **field_ad, electrostatic_cutoff2, electrostatic_cutoff, max_spacing, maxelcut2,  invkappa2,  ***beta_mol, *betamean ,*mu_ind_mols, ****mu_ind_t, *beta_corr, *ft_beta_corr;
    int            max_i, isize0, ind0, indj, ind1;
    real           t, rmax2, rmax,  r, r_dist, r2, q_xi, dq;
    real          *inv_segvol, normfac, segvol, spherevol, prev_spherevol, invsize0, invgamma, invhbinw, inv_width,  theta=0, *theta_vec;
    rvec          *x, xcm, xcm_transl, dx,  *x_i1, xi, x01, x02, *arr_qvec, **arr_qvec_faces ,vec_polin, vec_polout, ***vec_pout_theta_gamma, ***vec_pin_theta_gamma;
    rvec           pol_perp, pol_par,vec_kout, vec_2kin, pol_in1, pol_in2, vec_kout_2kin ;
    rvec           xvec, yvec, zvec, *xmol, *xref, Emean;
    real          *qref;
    matrix         cosdirmat,invcosdirmat; 
    real           invvol, invvol_sum;
    t_Map         *Map=NULL;
    t_Kern        *Krr = NULL;
    t_Kern        *SKern_rho_O = NULL;
    t_Kern        *SKern_rho_H = NULL;
    t_Kern        *SKern_E = NULL;
    t_Ion         *Cation=NULL, *Anion=NULL;
    t_inputrec    *ir=NULL;
    t_complex   ***FT_pair_pot;
    matrix         box, box_pbc;
    rvec           grid_spacing, grid_invspacing;
    real           inv_std_dev_dens, dens_deb, inv_tot_npoints_local_grid;
    int            *gridsize;
    int            ePBC = -1, ePBCrdf = -1;
    int            nplots = 1;
    t_block       *mols = NULL;
    t_atom        *atom = NULL;
    t_pbc          pbc;
    gmx_rmpbc_t    gpbc = NULL;
    gmx_rng_t      rng = NULL;
    int            mol, a, molsize;
    int            atom_id_0, nspecies_0, atom_id_1, nspecies_1;
    int           *chged_atom_indexes, n_chged_atoms;
    real	   rr2, rnorm,fac,beta;
    rvec	   xj, deltar;
    real ecorrcut2;

    fprintf(stderr,"Initialize number of atoms, get charge indexes, the number of atoms in each molecule and get the reference molecule\n");
    atom = top->atoms.atom;
    mols = &(top->mols);
    isize0 = isize[0];
//    isize0 = 2;
    molsize = mols->index[molindex[0][1]] - mols->index[molindex[0][0]];
    snew(xmol,molsize);
    snew(xref,molsize);
    snew(qref,molsize);
    nfaces = 6;
    invsize0 = 1.0/isize0;
    invgamma = 1.0/nbingamma;
    ecorrcut2 = ecorrcut*ecorrcut;
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    // need to know the index of oxygen atoms in a molecule and of the hydrogens
    // also need to know the number of oxygen and hydrogens in each molecule
    // if all works we need to make this better and use topology and related gromacs functions
    atom_id_0 = 0;
    nspecies_0 = 1;
    atom_id_1 = 1;
    nspecies_1 = 2;
    n_chged_atoms = 0;
    for (i = 0; i < molsize; i++)
    {
        qref[i] = top->atoms.atom[i].q;
        if (top->atoms.atom[i].q != 0.0 )
        {
           n_chged_atoms ++;
        }
    }
    snew(chged_atom_indexes, n_chged_atoms);
    ind0 = 0;
    for (i = 0; i < molsize; i++)
    {
        if (top->atoms.atom[i].q != 0.0 )
        {
          chged_atom_indexes[ind0] = i;
          ind0++;
        }
    }

    if (ePBC == -1)
    {
        ePBC = guess_ePBC(box);
    }
    copy_mat(box, box_pbc);
    ePBCrdf = ePBC;
    set_pbc(&pbc, ePBCrdf, box_pbc);

    if (fnREFMOL)
    {
      read_reference_mol(fnREFMOL,&xref);
    }
    else
    {
       xref[0][XX] = 0.0; xref[0][YY] = 0.0; xref[0][ZZ] = 0.0;
       xref[1][XX] = 0.075695; xref[1][YY] = 0.0; xref[1][ZZ] = 0.0585882;
       xref[2][XX] = -0.075695; xref[2][YY] = 0.0; xref[2][ZZ] = 0.0585882;
       xref[3][XX] = 0.0; xref[3][YY] = 0.0; xref[3][ZZ] = 0.01546;
    }

    fprintf(stderr,"\nNumber of atoms %d\n",natoms);
    fprintf(stderr,"\nNumber of molecules %d\n",isize0);
    fprintf(stderr,"\nNumber of atoms in molecule %d\n",molsize);
    fprintf(stderr,"\nName of group %s\n",grpname[0]);
    fprintf(stderr,"\nNumber of charged species in molecule %d\n",n_chged_atoms);

    // read electrostatic fit map input file
    if (kern[0] == 'm' || kern[0] == 'n' )
    {
       Map=(t_Map *)calloc(1,sizeof(t_Map));
       readMap(fnMAP, Map);
       if (kern[0] == 'm')
       {
          fprintf(stderr,"initialized electric field map to compute beta\n");
       }
       else
       {
          fprintf(stderr,"read the constant beta, beta will not fluctuate with the environment\n");
       }
    }
    else if (kern[0] == 'k')
    {
       Krr = (t_Kern *)calloc(1,sizeof(t_Kern));
       readKern(fnVCOEFF, fnVGRD, fnVINP, 0, 0, 0, NULL, FALSE, NULL ,Krr);
       Krr->kerndev = 0.5/((kernstd*kernstd));
       fprintf(stderr,"initialized kernel ridge regression to compute beta with standard dev = %f\n", kernstd);
    }
    else if (kern[0] == 's')
    {
       fprintf(stderr,"about to initialize scalar kernel \n");
       SKern_E = (t_Kern *)calloc(1,sizeof(t_Kern));
       readKern(fnVCOEFF, fnVGRD, NULL, kern_order, 0, kappa, xref, FALSE, &betamean, SKern_E);
       fprintf(stderr,"scalar kernel coefficients for the electric field read\n");
       SKern_rho_O = (t_Kern *)calloc(1,sizeof(t_Kern));
       readKern(fnCOEFFO, fnRGRDO, NULL, kern_order, std_dev_dens, 0, xref, TRUE, NULL,SKern_rho_O);
       fprintf(stderr,"scalar kernel coefficients for the oxygen density read\n");
       SKern_rho_H = (t_Kern *)calloc(1,sizeof(t_Kern));
       readKern(fnCOEFFH, fnRGRDH, NULL, kern_order, std_dev_dens, 0, xref, FALSE, NULL,SKern_rho_H);
       fprintf(stderr,"scalar kernel coefficients for the hydrogen density read\n");
       fprintf(stderr,"initialized scalar kernel \n");
       fprintf(stderr,"the density for the scalar kernel for each grid point i and for all atomic species j\n");
       fprintf(stderr,"is computed using gaussians rho_i = sum_j exp(-(x_i-x_j)^2/(2*std_dev^2))*weight_i\n");
       fprintf(stderr,"with a standard deviation of %f\n",std_dev_dens);
       fprintf(stderr,"the total number of training points is %d\n",SKern_rho_H->gridpoints + SKern_rho_O->gridpoints + SKern_E->gridpoints);
       inv_tot_npoints_local_grid = 1.0/(SKern_rho_H->gridpoints + SKern_rho_O->gridpoints + SKern_E->gridpoints);
    }

    if (bIONS )
    {
       // search for cations
       if (catname)
       {
          ncations = check_ion(top, catname);
          Cation = (t_Ion *)calloc(ncations,sizeof(t_Ion));
          identifyIon(top,Cation,catname);
       }
       else
       {
          gmx_fatal(FARGS,"Wrong cation name or cation not specified\n. Specify correct name of cation\n");
       }
   
       // search for anions
       if (anname)
       {
          nanions=check_ion(top, anname);
          Anion=(t_Ion *)calloc(nanions,sizeof(t_Ion));
          identifyIon(top, Anion, anname);
       }
       else
       {
          gmx_fatal(FARGS,"Wrong anion name or anion not specified\n. Specify correct name of anion\n");
       }

       if (nanions > 0)
       {
          if (Anion[0].q[0] > 0)
          {
             gmx_fatal(FARGS,"wrong definition of anion, check topology\n");
          }
          fprintf(stderr,"Name of anion %s\n Number of anions %d\n",anname, nanions );       
       }
       if (ncations > 0)
       {
          if (Cation[0].q[0] < 0)
          {
             gmx_fatal(FARGS,"wrong definition of cation, check topology\n");
          }
          fprintf(stderr,"Name of cation %s\n Number of cations %d\n",catname, ncations );
       }
    }

    if (bPBC)
    {
        rmax2   =  max_cutoff2(FALSE ? epbcXY : epbcXYZ, box_pbc);
        nbin     = (int)(sqrt(rmax2)/binwidth);
        invhbinw = 1.0 / binwidth;

        electrostatic_cutoff2 =  min(cutoff_field*cutoff_field, rmax2) ;
        fprintf(stderr, "rmax2 = %f\n", rmax2);
        maxelcut2 = maxelcut*maxelcut; 
        if (fnBETACORR)
        {
            fprintf(stderr, "number of bins for <beta(0)*beta(r)> = %d\n", nbin);
            nfaces = 1;
            if (nbingamma >1 || nbintheta >1  )
            {
                gmx_fatal(FARGS, "when computing <beta(0)*beta(r)> choose nplanes = 1 and nbintheta = 1");
            }
        } 
        inv_width = ( method[0] != 'd') ? 1.0 : M_PI*0.5/(sqrt(maxelcut2)-sqrt(electrostatic_cutoff2)) ;
        electrostatic_cutoff = sqrt(electrostatic_cutoff2);
        if ((method[0] == 'd') && (nbingamma >1 || fnBETACORR))
        {
            gmx_fatal(FARGS, "when using the double summation method choose nplanes = 1\n Also, choose the single method to compute <beta(0)*beta(r)>");
        }
        if ((method[0] == 'd') && (electrostatic_cutoff > maxelcut || maxelcut > sqrt(rmax2)  || electrostatic_cutoff > sqrt(rmax2))  )
        {
           fprintf(stderr,"electrostatic_cutoff =%f maxcutoff=%f rmax=%f\n",electrostatic_cutoff,maxelcut,sqrt(rmax2));
           gmx_fatal(FARGS, "wrong choice of cutoffs to truncate the potential or to compute the double sum, choose cutoffs appropriately\n");
        }
        if (kern[0] == 's')
        {
           ir=(t_inputrec *)calloc(1,sizeof(t_inputrec));
           ir->nkx = roundf(box[XX][XX]/fspacing); ir->nky = roundf(box[YY][YY]/fspacing); ir->nkz = roundf(box[ZZ][ZZ]/fspacing); ir->fourier_spacing = fspacing;
//           fprintf(flog,"Compute the intensity using ewald sums\n"); 
//           fprintf(flog, "screening kappa parameter used in ewald = %f\n", kappa);
//           fprintf(flog, "hard core smoothening parameter = %f\n", core_term);
//           if (ir->nkx > 0 && ir->nky > 0 && ir->nkz > 0)
//           {
              /* Mark fourier_spacing as not used */
//              ir->fourier_spacing = 0;
//           }
//           else if (ir->nkx != 0 && ir->nky != 0 && ir->nkz != 0)
//           {
//               gmx_fatal(FARGS, "Some of the Fourier grid sizes are set, but all of them need to be set.");
//           }
//           printf("about to build grid\n");
//           max_spacing = calc_grid(flog, box, ir->fourier_spacing,
//                                &(ir->nkx), &(ir->nky), &(ir->nkz));

           fprintf(stderr,"f spacing %f\n",ir->fourier_spacing);
           fprintf(stderr,"make the grid\n");
           snew(gridsize, DIM);
           initialize_free_quantities_on_grid(SKern_rho_O, ir, &grid_spacing,&grid_invspacing, FALSE, TRUE, &gridsize); 
           initialize_free_quantities_on_grid(SKern_rho_H, ir, &grid_spacing,&grid_invspacing, FALSE, TRUE, &gridsize);
           initialize_free_quantities_on_grid(SKern_E, ir, &grid_spacing, &grid_invspacing, TRUE, TRUE, &gridsize);
           inv_std_dev_dens = 0.5/(std_dev_dens*std_dev_dens);
           fprintf(stderr,"grid made and quantities on global grid allocated\n");
           fprintf(stderr,"initialize ewald pair potential\n");
           invkappa2 = 1.0/(sqr((box[XX][XX]+box[YY][YY]+box[ZZ][ZZ])/3.0)*kappa*kappa);
           //invkappa2=1.0/(kappa*kappa);
           //kmax is the upper limit in the summation over wave-vectors in the reciprocal term of spme
           //kmax is defined as kmax=int(walpha*max(box_lengths))
           //walpha=pi/rcut and rcut =min(box_lengths)*min(0.5,1.2*natoms^(-1/6))
           //kmax = (int)(max(max(box[XX][XX],box[YY][YY]),box[ZZ][ZZ])*M_PI/(min(min(box[XX][XX],box[YY][YY]),box[ZZ][ZZ])*min(0.5,1.2*pow(natoms,-1.0/6.0))));
           fprintf(stderr,"kappa is %f kappa^-2 (in fractional coords) is %f\n",kappa,invkappa2);
           fprintf(stderr,"number of fourier components for spme is 4pi/3 * the cube of %d\n",kmax);
           setup_ewald_pair_potential(gridsize,interp_order,kmax,&FT_pair_pot,invkappa2);
           fprintf(stderr,"ewald pair potential has been set-up\n");
        
        }
        if (method[0] =='d')
        {
           fprintf(stderr, "cutoff for electric field calculation = %f\n", sqrt(electrostatic_cutoff2));
           fprintf(stderr, "switching parameter (PI/2)*1/(max_cutoff - electrostatic_cutoff) = %f\n", inv_width);
        }
        
    }
    else
    {
        rmax2   = sqr(3*max(box[XX][XX], max(box[YY][YY], box[ZZ][ZZ])));
    }
    if (debug)
    {
        fprintf(debug, "rmax2 = %g\n", rmax2);
    }

    set_pbc(&pbc, ePBCrdf, box_pbc);
    rmax     = sqrt(rmax2);
    //initialize beta tensor
    snew(beta_mol, DIM);
    snew(mu_ind_mols, isize0);
    
    snew(beta_corr, nbin+1);
    snew(ft_beta_corr,nbinq);
    for (i = 0; i < DIM; i++)
    {
        snew(beta_mol[i], DIM);
    }
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            snew(beta_mol[i][j], DIM);
        }
    }

    snew(s_method, ng);
    snew(s_method_coh, ng);
    snew(s_method_incoh, ng);
    snew(s_method_t, ng);
    snew(s_method_coh_t, ng);
    snew(s_method_incoh_t, ng);
    max_i = isize[0];

    /*allocate memory for beta in lab frame and initialize beta in mol frame*/

    for (g = 0; g < ng; g++)
    {
        /*allocate memory for s_method array */
           snew(s_method[g], nbinq);
           snew(s_method_coh[g], nbinq);
           snew(s_method_incoh[g],nbinq);
           snew(temp_method, nbinq);
           snew(arr_qvec,nbinq);
           /*initialize incoming and outcoming wave-vectors*/
           vec_kout[XX] = 1.0; 
           vec_kout[YY] = 0.0;
           vec_kout[ZZ] = 0.0;
           vec_2kin[XX] = 0.0;
           vec_2kin[YY] = 0.0;
           vec_2kin[ZZ] = 1.0;
           /* compute the polarization vectors for outcoming beam*/     
           cprod(vec_kout, vec_2kin, pol_perp);
           cprod(vec_kout, pol_perp, pol_par);
           //fprintf(stderr,"polarization vector perpendicular to outcoming beam %f %f %f\n",pol_perp[XX], pol_perp[YY], pol_perp[ZZ]);
           //fprintf(stderr,"polarization vector parallel to Outcoming beam %f %f %f\n",pol_par[XX], pol_par[YY], pol_par[ZZ]);
           svmul(sin(M_PI/180.0*pout_angle), pol_perp, pol_perp);
           svmul(cos(M_PI/180.0*pout_angle), pol_par,  pol_par);
           rvec_add(pol_perp, pol_par, vec_polout);
           /* compute the polarization vectors for incoming beam*/                  
           cprod(vec_kout, vec_2kin, pol_perp);
           cprod(vec_2kin, pol_perp, pol_par);
           //fprintf(stderr,"polarization vector perpendicular to incoming beam %f %f %f\n",pol_perp[XX], pol_perp[YY], pol_perp[ZZ]);
           //fprintf(stderr,"polarization vector parallel to incoming %f %f %f\n",pol_par[XX], pol_par[YY], pol_par[ZZ]);
           svmul(sin(M_PI/180.0*pin_angle), pol_perp, pol_perp);
           svmul(cos(M_PI/180.0*pin_angle), pol_par , pol_par);
           rvec_add(pol_perp, pol_par, vec_polin);
           /*normalize kout, kin, vec_polout, vec_polin */
           unitv(vec_kout, vec_kout);
           unitv(vec_2kin, vec_2kin);
           rvec_sub(vec_kout, vec_2kin, vec_kout_2kin);
           unitv(vec_kout_2kin, vec_kout_2kin);
           unitv(vec_polin, vec_polin);
           unitv(vec_polout, vec_polout);
           /*----------------------------------------------------------------------*/         
 
           snew(cos_t, nbinq);
           snew(sin_t, nbinq);
           theta0  = theta0*M_PI/180.0 ;
           if (method[0] == 's')
           {
              qnorm = M_PI*2.0/(rmax*2.0)*qbin;
           }
           else
           {
              qnorm = M_PI*2.0/(rmax);
           }
           printf("----INITIALIZE DIRECTION OF INCOMING AND OUTCOMING WAVE-VECTORS-------------\n");
           printf("direction of incoming wave-vector is %f %f %f\n", vec_2kin[XX], vec_2kin[YY], vec_2kin[ZZ]);
           printf("direction of outcoming wave-vector is %f %f %f\n", vec_kout[XX], vec_kout[YY], vec_kout[ZZ]);
           printf("----INITIALIZE DIRECTION OF INCOMING AND OUTCOMING POLARIZATION VECTORS-----\n");
           printf("polarization of incoming wave-vector is %f %f %f\n", vec_polin[XX], vec_polin[YY], vec_polin[ZZ]);
           printf("polarization of outcoming wave-vector is %f %f %f\n", vec_polout[XX], vec_polout[YY], vec_polout[ZZ]);
           printf("----INITIALIZE DIRECTION OF SCATTERED WAVE VECTOR: q=kout-2kin -------------\n");
           printf("direction of scattered wave-vector is %f %f %f\n", vec_kout[XX]-vec_2kin[XX], vec_kout[YY]-vec_2kin[YY], vec_kout[ZZ]-vec_2kin[ZZ]);
           if (method[0] == 's')
           { 
              printf("minimum wave-vector is (2pi/L)*qbin = %f\n", qnorm);
              printf("maximum wave-vector is (2pi/L)*qbin*nbinq = %f\n", qnorm*nbinq);
           }
           else
           {
              printf("minimum wave-vector is (2pi/(L/2)) = %f\n", qnorm);
              printf("maximum wave-vector is (2pi/(L/2)) +(2pi/(L/2))/qbin * nbinq = %f\n", qnorm + qnorm/qbin*nbinq);
           }
           for (qq = 0; qq< nbinq; qq++)
           {
              // the magnitude of the scattered wave-vector has to be sqrt(2)*2pi/L*n, so we multiply by sqrt(2) since vec_kout_2kin is normalized
              arr_qvec[qq][XX] = sqrt(2.0)*(qnorm + qnorm*qq)*(vec_kout_2kin[XX]) ;
              arr_qvec[qq][YY] = sqrt(2.0)*(qnorm + qnorm*qq)*(vec_kout_2kin[YY]) ;
              arr_qvec[qq][ZZ] = sqrt(2.0)*(qnorm + qnorm*qq)*(vec_kout_2kin[ZZ]) ;
           }
           if (arr_qvec[0][YY] != 0.0 && arr_qvec[0][XX] != sqrt(2.0)*(qnorm + qnorm*qq)*1.0 && arr_qvec[0][ZZ] != sqrt(2.0)*(qnorm + qnorm*qq)*(-1.0))
           {
              gmx_fatal(FARGS,"the direction of the scattered wave-vector is not x - z.\n It is sufficient to compute the intensity using this direction of the scattered wave-vector, choose directions of incoming and outcoming wave-vectors such that q = kout -2kin = (x - z)*2Pi/L*n, where n is an integer  \n");
           }
  
           snew(s_method_t[g], nfaces);
           snew(s_method_coh_t[g], nfaces);
           snew(s_method_incoh_t[g], nfaces);
           snew(arr_qvec_faces,nfaces);
           snew(vec_pout_theta_gamma,nfaces);
           snew(vec_pin_theta_gamma, nfaces);
           snew(theta_vec, nbintheta);
           // this loop is to get the scattering wave-vector and polarization vectors for different faces of the cube
           printf("\n----COMPUTE THE POLARIZATION VECTORS AT DIFFERENT SCATTERING ANGLES THETA.---------------------\n");
           printf("----THETA IS DEFINED AS THE ANGLE BETWEEN THE INCOMING AND OUTCOMING WAVE-VECTORS.---------------\n");  
           printf("----THE POLARIZATION VECTORS ARE ALSO COMPUTED AT DIFFERENT ANGLES GAMMA.------------------------\n");
           printf("----GAMMA IS THE ANGLE DEFINED BY A PLANE THAT GOES THROUGH THE SCATTERING WAVE-VECTOR AND-------\n");
           printf("----THE PLANE PARALLEL TO THE CHOSEN FACE OF THE SIMULATION BOX----------------------------------\n \n");
           for (rr = 0; rr< nfaces; rr++)
           {
              snew(s_method_t[g][rr], nbinq);
              snew(s_method_coh_t[g][rr], nbinq);
              snew(s_method_incoh_t[g][rr], nbinq);
              snew(arr_qvec_faces[rr],nbinq);
              for (qq = 0; qq< nbinq; qq++)
              {
                   //now the magnitude of the wave-vector is indeed sqrt(2)*2pi/L*n so we don't need to multiply by sqrt(2).
                 if (method[0] == 's')
                 {
                    arr_qvec_faces[rr][qq][XX] = (qnorm + qnorm*qq)*(1.0) ;
                    arr_qvec_faces[rr][qq][YY] = 0.0 ;
                    arr_qvec_faces[rr][qq][ZZ] = (qnorm + qnorm*qq)*(-1.0) ;
                 }
                 else
                 {
                    arr_qvec_faces[rr][qq][XX] = (qnorm + qnorm*qq/qbin)*(1.0) ;
                    arr_qvec_faces[rr][qq][YY] = 0.0 ;
                    arr_qvec_faces[rr][qq][ZZ] = (qnorm + qnorm*qq/qbin)*(-1.0) ;
                 }
                 rotate_wave_vec(arr_qvec_faces[rr][qq], rr, arr_qvec_faces[rr][qq]);
                 snew(s_method_t[g][rr], nbintheta);
                 snew(s_method_coh_t[g][rr], nbintheta);
                 snew(s_method_incoh_t[g][rr], nbintheta);          
                 snew(vec_pout_theta_gamma[rr], nbintheta);
                 snew(vec_pin_theta_gamma[rr], nbintheta);
 
                 for(tt = 0 ; tt < nbintheta; tt++)
                 {
                     snew(s_method_t[g][rr][tt],nbinq);
                     snew(s_method_coh_t[g][rr][tt],nbinq);
                     snew(s_method_incoh_t[g][rr][tt],nbinq);
                     snew(vec_pout_theta_gamma[rr][tt], nbingamma);
                     snew(vec_pin_theta_gamma[rr][tt], nbingamma);
                 }
                 for (tt = 0; tt < nbintheta; tt++)
                 {
                         //theta is the angle between the outcoming wave-vector and the scattered wave-vector
                         //the experimental_theta is the angle between the incoming and outcoming wave-vectors and it is twice this value
                         theta_vec[tt] = ( tt< nbintheta*0.5 ) ? /*5.0/(6.0*2.0)*/ 0.5*M_PI*(-1.0 + tt*2.0/nbintheta) : /*5.0/(6.0*2.0)**/ 0.5*M_PI*(-1.0 + 2.0/nbintheta + tt*2.0/nbintheta) ;
                         // we get very close to ±90 degrees because kin and kout -> infinity at theta=±90 degrees
                         theta_vec[0] = 0.5*M_PI*(-1.0)*0.999;
                         theta_vec[nbintheta-1] = 0.5*M_PI*(-1.0 + 2.0/nbintheta + (nbintheta-1)*2.0/nbintheta)*0.999 ;
                         if (fnBETACORR || fnFTBETACORR)
                         {
                            theta_vec[tt] = angle_corr*M_PI/(2.0*180.0);
                         }
                         // this loop is to  rotate the scattering plane wrt the scattering wave-vector
                         for (c = 0; c < nbingamma; c++)
                         {
                             gamma = c*M_PI*invgamma;
                             // compute the outcoming wave-vector as a function of the angles gamma and theta, and the corresponding polarization vector                       
                             vec_kout[XX] = 0.5*(1.0 + tan(theta_vec[tt])*cos(gamma));
                             vec_kout[YY] = 0.5*(sqrt(2)*tan(theta_vec[tt])*sin(gamma));
                             vec_kout[ZZ] = 0.5*(-1.0 + tan(theta_vec[tt])*cos(gamma));
                             rotate_wave_vec(vec_kout, rr, vec_kout);
                             cprod(vec_kout, arr_qvec_faces[rr][0], pol_perp);
                             cprod(pol_perp, vec_kout, pol_par);
                             svmul(sin(M_PI/180.0*pout_angle), pol_perp, pol_perp);
                             //fprintf(stderr,"polarization vector perpendicolar to outcoming beam  %f %f %f\n",pol_perp[XX], pol_perp[YY], pol_perp[ZZ]);
                             svmul(cos(M_PI/180.0*pout_angle), pol_par,  pol_par);
                             //fprintf(stderr,"polarization vector parallel to outcoming beam %f %f %f\n",pol_par[XX], pol_par[YY], pol_par[ZZ]);
                             rvec_add(pol_perp, pol_par, vec_pout_theta_gamma[rr][tt][c]);
                             unitv( vec_pout_theta_gamma[rr][tt][c], vec_pout_theta_gamma[rr][tt][c]);
         
                             // compute the incoming wave-vector as a function of the angles gamma and theta, and the corresponding polarization vector
                             vec_2kin[XX] = -0.5*(1.0 - tan(theta_vec[tt])*cos(gamma));
                             vec_2kin[YY] = -0.5*(-sqrt(2)*tan(theta_vec[tt])*sin(gamma));
                             vec_2kin[ZZ] = -0.5*(-1.0 - tan(theta_vec[tt])*cos(gamma));
                             rotate_wave_vec(vec_2kin, rr, vec_2kin);         
                             cprod(vec_2kin, arr_qvec_faces[rr][0], pol_perp);
                             cprod(pol_perp, vec_2kin, pol_par);
                             svmul(sin(M_PI/180.0*pin_angle), pol_perp, pol_perp);
                             //fprintf(stderr,"polarization vector perpdicular to incoming beam %f %f %f\n",pol_perp[XX], pol_perp[YY], pol_perp[ZZ]);
                             svmul(cos(M_PI/180.0*pin_angle), pol_par , pol_par);
                             //fprintf(stderr,"polarization vector parallel to incoming %f %f %f\n",pol_par[XX], pol_par[YY], pol_par[ZZ]);
                             rvec_add(pol_perp, pol_par, vec_pin_theta_gamma[rr][tt][c]);
                             unitv(vec_pin_theta_gamma[rr][tt][c], vec_pin_theta_gamma[rr][tt][c]);
                             printf("polarization vectors at angles theta_expt = %f gamma = %f and face index %d \n",M_PI+2.0*theta_vec[tt]*180.0/M_PI, gamma*180.0/M_PI, rr);
                             printf("incoming polarization vector = %f %f %f \n",vec_pin_theta_gamma[rr][tt][c][XX], vec_pin_theta_gamma[rr][tt][c][YY], vec_pin_theta_gamma[rr][tt][c][ZZ]);
                             printf("outcoming polarization vector = %f %f %f \n",vec_pout_theta_gamma[rr][tt][c][XX], vec_pout_theta_gamma[rr][tt][c][YY], vec_pout_theta_gamma[rr][tt][c][ZZ]);
                             printf("direction of scattered wave-vector = %f %f %f \n", vec_kout[XX] -vec_2kin[XX], vec_kout[YY] - vec_2kin[YY], vec_kout[ZZ] -vec_2kin[ZZ]);
                             check_pol = iprod(vec_pin_theta_gamma[rr][tt][c],vec_pout_theta_gamma[rr][tt][c]);
                             printf("(incoming polarization vec) dot (outcoming polarization vec) = %.17g\n",check_pol);
                         }
                 }
              }
              if ((rr == 0) || (rr == 2) || (rr == 4))
              {
                 printf("smallest scattering wavevector along one of the diagonals of the faces of the simulation box at |q| = %f nm^-1\n",norm(arr_qvec_faces[rr][0]));
              }
              else if ((rr == 1) || (rr == 3) || (rr == 5))
              {
                 printf("smallest scattering wavevector along one of the diagonals of sides of the simulation box at |q| = %f nm^-1\n",norm(arr_qvec_faces[rr][0]));
              }
              printf("q = %f %f %f\n", arr_qvec_faces[rr][0][XX], arr_qvec_faces[rr][0][YY], arr_qvec_faces[rr][0][ZZ]);
           }
    }
    
    snew(x_i1, max_i);
    nframes    = 0;
    invvol_sum = 0;
    if (bPBC && (NULL != top))
    {
        gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms);
    }

    rng=gmx_rng_init(gmx_rng_make_seed());

    if (method[0] == 's' )
    {
        fprintf(stderr,"using a single sum to compute intensity\n");
        do
        {
            copy_mat(box, box_pbc);
            if (top != NULL)
            {
                gmx_rmpbc(gpbc, natoms, box, x);
            }
            set_pbc(&pbc, ePBCrdf, box_pbc);
            invvol      = 1/det(box_pbc);
            invvol_sum += invvol;
            
            for (g = 0; g < ng; g++)
            {
                snew(cos_tq, nfaces);
                snew(sin_tq, nfaces);
                snew(mu_sq_t, nfaces); 

                for (rr = 0; rr < nfaces; rr++)
                {
                   snew(cos_tq[rr], nbintheta);
                   snew(sin_tq[rr], nbintheta);
                   snew(mu_sq_t[rr], nbintheta);
                   for (tt = 0; tt < nbintheta; tt++)
                   {
                       snew(cos_tq[rr][tt],nbingamma);
                       snew(sin_tq[rr][tt],nbingamma);
                       snew(mu_sq_t[rr][tt], nbingamma);
                       for (c = 0; c <nbingamma; c++)
                       {
                           snew(cos_tq[rr][tt][c],nbinq);
                           snew(sin_tq[rr][tt][c],nbinq);
                       }
                   }
                }
                 
                if (kern[0] == 's')
                {
                   fprintf(stderr,"put atoms in box\n");
                   put_atoms_in_box(ePBC, box, natoms, x);
                   fprintf(stderr,"about to compute density on a grid\n");
                   calc_dens_on_grid(SKern_rho_O, ir, &pbc,  mols, molindex, 
                                     atom_id_0, nspecies_0, isize0, x, std_dev_dens,
                                     inv_std_dev_dens, grid_invspacing, gridsize,grid_spacing);
                   fprintf(stderr,"computed O dens\n");
                   calc_dens_on_grid(SKern_rho_H, ir, &pbc,
                                     mols, molindex, atom_id_1, nspecies_1, isize0, x, std_dev_dens,
                                     inv_std_dev_dens, grid_invspacing, gridsize, grid_spacing);
                   fprintf(stderr,"computed H dens\n");
                   
                   calculate_spme_efield(SKern_E, ir, top, box, invvol, mols, molindex,
                                       chged_atom_indexes,n_chged_atoms,
                                       gridsize, grid_spacing, interp_order, x, isize0, FT_pair_pot, &Emean,eps);
                   fprintf(stderr,"computed electric field with spme\n");
                   calc_efield_correction(SKern_E, ir, top, &pbc, box, invvol, mols, molindex,
                                       chged_atom_indexes, n_chged_atoms,
                                       gridsize, grid_spacing, grid_invspacing, x, isize0, sigma_vals, ecorrcut,
                                       ecorrcut2, gridsize);
                   fprintf(stderr,"computed real-space correction to electric field\n");

                   //fprintf(stderr,"average field %f %f %f\n", Emean[XX], Emean[YY], Emean[ZZ]);
                }
                for (i = 0; i < isize0; i++)
                {
                    ind0  = mols->index[molindex[g][i]];
                    copy_rvec(x[ind0], xi);
                    for (aa = 0; aa < molsize; aa++)
                    {
                       pbc_dx(&pbc,x[ind0+aa],x[ind0],xmol[aa]);
                       //printf("O %f %f %f\n",x[ind0+aa][XX]*10.0,x[ind0+aa][YY]*10.0,x[ind0+aa][ZZ]*10.0);
                    }
                    calc_cosdirmat( fnREFMOL, top, molsize, ind0,  xref, xmol, &cosdirmat, &invcosdirmat, &xvec, &yvec, &zvec );
                    if (kern[0] == 'm')
                    {
                        calc_efield_map(&pbc, top, mols, Cation, Anion, molindex, g, isize0, ncations, nanions,
                                x, ind0, xvec, yvec, zvec, electrostatic_cutoff2, &field_ad);
                        calc_beta_efield_map(Map, field_ad, beta_mol, &beta_mol );
                    }
                    else if (kern[0] == 's' )
                    {

                        rotate_local_grid_points(SKern_rho_O, SKern_rho_H, SKern_E, ePBC, box, cosdirmat, xi );
                        //fprintf(stderr,"finished rotating and translating grid\n");
                        trilinear_interpolation_kern(SKern_rho_O, ir, &pbc, xi, grid_invspacing, grid_spacing);
                        //fprintf(stderr,"finished interpolation O kern\n");
                        trilinear_interpolation_kern(SKern_rho_H, ir, &pbc, xi, grid_invspacing, grid_spacing);
                        //fprintf(stderr,"finished interpolation H kern\n");
                        if (debug)
                        {
                              //vec_trilinear_interpolation_kern(SKern_E, ir, &pbc, invcosdirmat, xi,
//                                                         grid_invspacing, grid_spacing, Emean);
                            bspline_efield(SKern_E, ir, &pbc, invcosdirmat, xi,
                                           grid_invspacing, interp_order, grid_spacing, Emean);


			    // Now loop over all other molecules and find the real-space electric field at this position.
/*			    for (j=0;j<isize0;j++)
			    {
				if (i!=j)
				{
					ind1 = mols->index[molindex[g][j]];
					copy_rvec(x[ind1],xj);
					pbc_dx(&pbc,xi,xj,deltar);
					rnorm = sqrt(deltar[0]*deltar[0] + deltar[1]*deltar[1] + deltar[2]*deltar[2]);
					rr2 = rnorm*rnorm;
					beta = kappa;
					fac = (2.0 * beta / sqrt(M_PI))*exp(-beta*beta*rr2);
					fac += erfc(beta*rnorm)/rnorm;
					SKern_E->vec_interp_quant_grid[0][0] += fac*deltar[0]/rr2;
					SKern_E->vec_interp_quant_grid[0][1] += fac*deltar[1]/rr2;
					SKern_E->vec_interp_quant_grid[0][2] += fac*deltar[2]/rr2;
				}
			    }
                            printf("force electrostatic %d %f %f %f\n",i,SKern_E->vec_interp_quant_grid[0][0],
                                             SKern_E->vec_interp_quant_grid[0][1],SKern_E->vec_interp_quant_grid[0][2]);
*/
                        }
                        else
                        {

                        vec_trilinear_interpolation_kern(SKern_E, ir, &pbc, invcosdirmat, xi,
                                                         grid_invspacing, grid_spacing, Emean);
                        }
                        //fprintf(stderr,"finished interpolation E kern\n");

                        calc_beta_skern(SKern_rho_O, SKern_rho_H, SKern_E, kern_order, betamean, &beta_mol);
                        
/*                        if (debug)
                        {
                           for (aa = 0; aa < DIM; aa++)
                           {
                               for (bb = 0; bb < DIM; bb++)
                               {
                                  for (cc = 0; cc < DIM ; cc++)
                                  {
                                     printf("beta %d %d %d %f\n",aa, bb, cc, beta_mol[aa][bb][cc]);
                                  }
                               }
                           }
                        }
*/

//                        if (debug )
//                        {
//                           gmx_fatal(FARGS,"end\n");
//                        }

                    }
                    else if (kern[0] == 'k')
                    {
                        for (gr_ind =0; gr_ind < Krr->gridpoints; gr_ind++)
                        {
                           mvmul(cosdirmat,Krr->grid[gr_ind],SKern_E->rotgrid[gr_ind]);
                        }
                       //copy_rvec(x[ind0],xcm_transl);
                       //we use this only temporarily, because the centre water molecule has been centred
                       //in the centre of core charge
                       svmul(0.0117176,zvec,xcm);
                       rvec_add(xcm,x[ind0],xcm_transl);
                       calc_beta_krr(Krr, &pbc, top,  mols, molindex, g, isize0, x, xcm_transl, ind0, rmax2, electrostatic_cutoff,&beta_mol );
                    }
                    else if (kern[0] == 'n')
                    {
                      for (aa = 0; aa < DIM; aa++)
                      {
                          for (bb = 0; bb < DIM; bb++)
                          {
                              for (cc = 0; cc < DIM; cc++)
                              {
                                 beta_mol[aa][bb][cc] = Map->beta_gas[aa*9+ bb*3+ cc];
                              }
                          }
                      }
                    }
                    for (rr = 0; rr < nfaces; rr++)
                    {
                        for (tt = 0; tt < nbintheta; tt++ )
                        {
                            for (c  = 0; c < nbingamma; c++)
                            {
                                if (kern[0] == 'n')
                                {
                                   induced_second_order_fluct_dipole(cosdirmat, 
                                                                     vec_pout_theta_gamma[rr][tt][c], vec_pin_theta_gamma[rr][tt][c], 
                                                                     beta_mol, &mu_ind);
                                }
                                else
                                {
                                   induced_second_order_fluct_dipole_fluct_beta(cosdirmat,
                                                                     vec_pout_theta_gamma[rr][tt][c], vec_pin_theta_gamma[rr][tt][c],
                                                                     beta_mol, &mu_ind);

                                }
                                mu_sq_t[rr][tt][c] += mu_ind*mu_ind;
                                mu_ind_mols[i] = mu_ind;
                                for (qq = 0; qq < nbinq; qq++)
                                {
                                    q_xi = iprod(arr_qvec_faces[rr][qq],xi);
                                    cos_tq[rr][tt][c][qq] += mu_ind*cos(q_xi);
                                    sin_tq[rr][tt][c][qq] += mu_ind*sin(q_xi);
                                }
                            }
                        }
                     }
                }
                for (rr = 0; rr < nfaces; rr++)
                {
                   for (tt = 0; tt < nbintheta; tt++)
                   {
                       for (c = 0; c < nbingamma; c++)
                       {
                           incoh_temp = mu_sq_t[rr][tt][c]*invsize0;
                           for (qq = 0; qq < nbinq; qq++)
                           {
                              tot_temp = (cos_tq[rr][tt][c][qq]*cos_tq[rr][tt][c][qq] + sin_tq[rr][tt][c][qq]*sin_tq[rr][tt][c][qq])*invsize0;
                              s_method_t[g][rr][tt][qq] +=  tot_temp  ;
                              s_method_coh_t[g][rr][tt][qq] += tot_temp - incoh_temp;
                              s_method_incoh_t[g][rr][tt][qq] += incoh_temp ;
                           }
                       }
                   }
                }
                if (fnBETACORR)
                {
                 calc_beta_corr( &pbc,  mols, molindex, g, isize0, nbin, rmax2, invhbinw, x, mu_ind_mols, &beta_corr);
                }
                else if (fnFTBETACORR)
                {
                 calc_ft_beta_corr( &pbc,  mols, molindex, g, isize0,  nbinq, arr_qvec_faces, rmax2, invhbinw, x, mu_ind_mols, &ft_beta_corr);
                }
            }
            for (rr = 0; rr < nfaces; rr++)
            {
               for (tt  = 0; tt < nbintheta; tt++)
               {
                   for (c = 0; c < nbingamma; c++)
                   {
                      sfree(cos_tq[rr][tt][c]);
                      sfree(sin_tq[rr][tt][c]);
                   }
                   sfree(cos_tq[rr][tt]);
                   sfree(sin_tq[rr][tt]);
                   sfree(mu_sq_t[rr][tt]);
               }
               sfree(cos_tq[rr]);
               sfree(sin_tq[rr]);
               sfree(mu_sq_t[rr]);
            }
            sfree(cos_tq);
            sfree(sin_tq);
            sfree(mu_sq_t);
            nframes++;
        }
        while (read_next_x(oenv, status, &t, x, box));
    }

    else if ( method[0] =='d')
    {
        fprintf(stderr,"do a double loop over atoms, more expensive but you can smoothen the intensity at low q values\n");
        fprintf(stderr,"the intensity is smoothened using a switching function between %f nm and %f nm\n",electrostatic_cutoff, maxelcut);
        fprintf(stderr, "switching parameter (PI/2)*1/(max_cutoff - electrostatic_cutoff) = %f\n", inv_width);

        snew(mu_ind_t, isize0);
        for (i = 0; i < isize0; i++)
        {
            snew(mu_ind_t[i],nfaces);
            for (rr = 0; rr < nfaces; rr ++) 
            {
                snew(mu_ind_t[i][rr],nbintheta);
                for (tt = 0; tt < nbintheta; tt++)
                {
                   snew(mu_ind_t[i][rr][tt],nbingamma);
                }
            }
        }
        do
        {
            // Must init pbc every step because of pressure coupling 
            copy_mat(box, box_pbc);
            if (top != NULL)
            {
                gmx_rmpbc(gpbc, natoms, box, x);
            }
            set_pbc(&pbc, ePBCrdf, box_pbc);
            invvol      = 1/det(box_pbc);
            invvol_sum += invvol;

            for (g = 0; g < ng; g++)
            {
                snew(mu_sq_t, nfaces);
                snew(coh_temp,nfaces);
                for (rr = 0; rr < nfaces; rr++)
                {
                   snew(mu_sq_t[rr], nbintheta);
                   snew(coh_temp[rr],nbintheta);
                   for ( tt = 0; tt < nbintheta; tt++)
                   {
                      snew(coh_temp[rr][tt],nbinq);
                      snew(mu_sq_t[rr][tt],nbingamma);
                   }
                }
                for (i = 0; i < isize0; i++)
                {
                    ind0  = mols->index[molindex[g][i]];
                    copy_rvec(x[ind0], xi);
                    for (aa = 0; aa < molsize; aa++)
                    {
                       pbc_dx(&pbc,x[ind0+aa],x[ind0],xmol[aa]);
                    }
                    calc_cosdirmat( fnREFMOL,  top, molsize, ind0,  xref, xmol, &cosdirmat, &invcosdirmat, &xvec, &yvec, &zvec );
                    if (kern[0] == 'm')
                    {
                        calc_efield_map(&pbc, top, mols, Cation, Anion, molindex, g, isize0, ncations, nanions,
                                x, ind0, xvec, yvec, zvec, electrostatic_cutoff2, &field_ad);
                        calc_beta_efield_map(Map, field_ad, beta_mol, &beta_mol );
                    }
                    else if (kern[0] == 'k')
                    {
                        for (gr_ind =0; gr_ind < Krr->gridpoints; gr_ind++)
                        {
                           mvmul(cosdirmat,Krr->grid[gr_ind],Krr->rotgrid[gr_ind]);
                        }
                        //copy_rvec(x[ind0],xcm_transl);
                        //we use this only temporarily, because the centre water molecule has been centred in the centre of core charge
                        svmul(0.0117176,zvec,xcm);
                        rvec_add(xcm,x[ind0],xcm_transl);
                        calc_beta_krr(Krr, &pbc, top,  mols, molindex, g, isize0, x, xcm_transl, ind0, rmax2, electrostatic_cutoff,&beta_mol );

                    }
                    else if (kern[0] == 'n')
                    {
                      for (aa = 0; aa < DIM; aa++)
                      {
                          for (bb = 0; bb < DIM; bb++)
                          {
                              for (cc = 0; cc < DIM; cc++)
                              {
                                 beta_mol[aa][bb][cc] = Map->beta_gas[aa*9+ bb*3+ cc];
                              }
                          }
                      }
                    }
                    for (rr = 0; rr < nfaces; rr++)
                    {
                        for (tt = 0; tt < nbintheta; tt++ )
                        {
                            for (c = 0 ; c < nbingamma; c++)
                            {

                                if (kern[0] == 'n')
                                {
                                   induced_second_order_fluct_dipole(cosdirmat,
                                                                     vec_pout_theta_gamma[rr][tt][c], vec_pin_theta_gamma[rr][tt][c],
                                                                     beta_mol, &mu_ind);
                                }
                                else
                                {
                                   induced_second_order_fluct_dipole_fluct_beta(cosdirmat,
                                                                     vec_pout_theta_gamma[rr][tt][c], vec_pin_theta_gamma[rr][tt][c],
                                                                     beta_mol, &mu_ind);

                                }
                                mu_ind_t[i][rr][tt][c] = mu_ind;
                                mu_sq_t[rr][tt][c] += mu_ind*mu_ind;
                               //printf("mu_ind_t %f \n", mu_ind_t[i][rr][tt][c]);
                            }
                        }
                    }
                }
                for (i = 0; i < isize0 -1; i++)
                {
                    ind0  = mols->index[molindex[g][i]];
                    copy_rvec(x[ind0], xi);
                    //printf("xi %f %f %f\n",xi[XX],xi[YY],xi[ZZ]);
                    for (j = i + 1; j < isize0 ; j++)
                    {
                        ind0 = mols->index[molindex[g][j]];
                        pbc_dx(&pbc, xi, x[ind0], dx);
                        r2 = iprod(dx, dx);
                        if (r2 <= maxelcut2 )
                        {
                           r_dist = sqrt(r2);
                           if ( r_dist <= electrostatic_cutoff)
                           {
                              for (rr = 0; rr < nfaces; rr++)
                              {
                                  for (tt = 0; tt < nbintheta; tt++)
                                  {
                                      for (c = 0; c < nbingamma; c++)
                                      {  mod_f = (mu_ind_t[i][rr][tt][c] * mu_ind_t[j][rr][tt][c]);
                                         //mod_f = (mu_ind_t[i][rr][tt][c] + mu_ind_t[j][rr][tt][c]);
                                         //mod_f *= mod_f;
                                         //fprintf(stderr,"mod_f %f i %d rr %d tt %d c %d \n",mod_f,i,rr,tt,c);
                                         for (qq = 0; qq < nbinq; qq++)
                                         {
                                             coh_temp[rr][tt][qq] += mod_f*cos(iprod(arr_qvec_faces[rr][qq],dx));
                                         }
                                      }
                                  }
                              }
                           }
                           else
                           {
                              for (rr = 0; rr < nfaces; rr++)
                              {
                                  for (tt = 0; tt < nbintheta; tt++)
                                  {
                                      for (c = 0; c < nbingamma; c++)
                                      {
                                          mod_f = (mu_ind_t[i][rr][tt][c] * mu_ind_t[j][rr][tt][c])*sqr(cos((r_dist-electrostatic_cutoff)*inv_width)) ;
                                          //mod_f = (mu_ind_t[i][rr][tt][c] + mu_ind_t[j][rr][tt][c])*cos((r_dist-electrostatic_cutoff)*inv_width) ;
                                          //mod_f *= mod_f;
                                          //fprintf(stderr,"mod_f2 %f i %d rr %d tt %d c %d \n",mod_f,i,rr,tt,c);
                                          for (qq = 0; qq < nbinq; qq++)
                                          {
                                             coh_temp[rr][tt][qq] += mod_f*cos(iprod(arr_qvec_faces[rr][qq],dx));
                                          //   fprintf(stderr,"coh_temp %f rr %d tt %d qq %d\n",coh_temp[rr][tt][qq],rr,tt,qq);
                                          }
                                      }
                                  }
                              }
                           }
                        }
                    }
                }
                for (rr = 0; rr < nfaces; rr++)
                {
                   for (tt = 0; tt < nbintheta; tt++)
                   {
                      for (c = 0; c < nbingamma; c++)
                      {
                         incoh_temp = mu_sq_t[rr][tt][c]*invsize0;
                         for (qq = 0; qq < nbinq; qq++)
                         {
                            s_method_t[g][rr][tt][qq] +=  2.0*coh_temp[rr][tt][qq]*invsize0 + incoh_temp  ;
                            s_method_coh_t[g][rr][tt][qq] += 2.0*coh_temp[rr][tt][qq]*invsize0 ;
                            s_method_incoh_t[g][rr][tt][qq] += incoh_temp ;
                            //printf("incoh_temp %f\n",s_method_incoh_t[g][rr][tt][qq]);
                         }
                      }
                   }
                }
            }
            for (rr = 0; rr < nfaces; rr++)
            {
               for (tt  = 0; tt < nbintheta; tt++)
               {
                   sfree(mu_sq_t[rr][tt]);
                   sfree(coh_temp[rr][tt]);
               }
               sfree(mu_sq_t[rr]);
               sfree(coh_temp[rr]);
            }
            sfree(mu_sq_t);
            sfree(coh_temp);
            nframes++;
        }
        while (read_next_x(oenv, status, &t, x, box));
    }
    
    fprintf(stderr, "\n");
    if (bPBC && (NULL != top))
    {
        gmx_rmpbc_done(gpbc);
    }
    
    close_trj(status);

    sfree(x);

    /* Average volume */
    invvol = invvol_sum/nframes;
    for (g = 0; g < ng; g++)
    {
        for (qq = 0; qq < nbinq ; qq++)
        {
           s_method[g][qq] = s_method[g][qq]/(nframes)  ;
           s_method_coh[g][qq] = s_method_coh[g][qq]/(nframes)  ;
           s_method_incoh[g][qq] = s_method_incoh[g][qq]/(nframes) ;
        }
    }

    sprintf(gtitle, "Non-linear optical scattering ");
    fp = xvgropen(fnTHETA, "S(theta)", "theta", "S(theta)", oenv);
    sprintf(refgt, "%s", "");
    fprintf(fp, "@    s%d legend \"incoherent\"\n",nplots);
    fprintf(fp, "@target G0.S%d\n",nplots);
    fprintf(fp, "@type xy\n");
    nplots++ ;
    for (tt = 0; tt < nbintheta ; tt++)
    {
       theta = (M_PI+2.0*theta_vec[tt])*180.0/M_PI;
       for (rr = 0; rr< nfaces; rr++)
       {
           fprintf(fp, "%10g", theta);
           fprintf(fp, " %10g", s_method_incoh_t[0][rr][tt][0]/nframes*invgamma);
           fprintf(fp, "\n");
       }
    }
    fprintf(fp,"&\n");
    for (qq = 0; qq < nbinq ; qq++)
    {
       for  (rr = 0; rr < nfaces; rr++)
       {
          fprintf(fp, "@    s%d legend \" coherent q=%g\" \n",nplots,norm(arr_qvec_faces[rr][qq]));
          fprintf(fp, "@target G0.S%d\n",nplots);
          fprintf(fp, "@type xy\n");
          nplots++ ;
          for (tt = 0; tt < nbintheta ; tt++)
          {
             theta = (M_PI+2.0*theta_vec[tt])*180.0/M_PI;
             fprintf(fp, "%10g", theta);
             fprintf(fp, " %10g", s_method_coh_t[0][rr][tt][qq]/nframes*invgamma);
             fprintf(fp, "\n");
          }
          fprintf(fp,"&\n");
          fprintf(fp, "@    s%d legend \"total q=%g\"\n",nplots,norm(arr_qvec_faces[rr][qq]));
          fprintf(fp, "@target G0.S%d\n",nplots);
          nplots++ ;
          fprintf(fp, "@type xy\n");
          for (tt = 0; tt < nbintheta  ; tt++) 
          { 
             theta = (M_PI+2.0*theta_vec[tt])*180.0/M_PI;
             fprintf(fp, "%10g", theta);
             fprintf(fp, " %10g", s_method_t[0][rr][tt][qq]/nframes*invgamma);
             fprintf(fp, "\n");
          }
          fprintf(fp,"&\n");
       }
    }
    gmx_ffclose(fp);
    
    if (!fnBETACORR )
    {
    // print the nonlinear scattering intensity as a function of wave-vector only if you don't compute the <beta(0)*beta(r)>
    nplots = 1;
    sprintf(gtitle, "Non-linear optical scattering ");
    fpn = xvgropen(fnSFACT, "S(q)", "q (nm^-1)", "S(q)", oenv);
    sprintf(refgt, "%s", "");
    for (tt = 0; tt < nbintheta ; tt++)
    {
       theta = (M_PI+ 2.0*theta_vec[tt])*180.0/M_PI;
       if ((round(abs(theta)) == 45.0) || (round(abs(theta)) == 30.0 ) || (round(abs(theta)) == 60.0 ) || (round(abs(theta)) == 90.0)
          || (round(abs(theta)) == 150.0)  || (round(abs(theta)) == 120.0) || (round(abs(theta)) == 10.0) || (tt ==  0) || (tt == nbintheta/2))
       {
           fprintf(fpn, "@    s%d legend \"incoherent theta=%g all faces\"\n",nplots,theta);
           fprintf(fpn, "@target G0.S%d\n",nplots);
           fprintf(fpn, "@type xy\n");
           for (qq = 0; qq< nbinq; qq++)
           {
               fprintf(fpn, "%10g", norm(arr_qvec_faces[1][qq]) );
               fprintf(fpn, " %10g", (s_method_incoh_t[0][1][tt][qq]+ s_method_incoh_t[0][3][tt][qq] + s_method_incoh_t[0][5][tt][qq])/nframes*invgamma/3.0);
               fprintf(fpn, "\n");
               fprintf(fpn, "%10g", norm(arr_qvec_faces[0][qq]) );
               fprintf(fpn, " %10g", (s_method_incoh_t[0][0][tt][qq]+ s_method_incoh_t[0][2][tt][qq] + s_method_incoh_t[0][4][tt][qq])/nframes*invgamma/3.0);
               fprintf(fpn, "\n");
           }
           fprintf(fpn,"&\n");
           nplots++;
           fprintf(fpn, "@    s%d legend \"coherent theta=%g all faces\"\n",nplots,theta);
           fprintf(fpn, "@target G0.S%d\n",nplots);
           fprintf(fpn, "@type xy\n");
           for (qq = 0; qq< nbinq; qq++)
           {
               fprintf(fpn, "%10g", norm(arr_qvec_faces[1][qq]) );
               fprintf(fpn, " %10g", (s_method_coh_t[0][1][tt][qq]+ s_method_coh_t[0][3][tt][qq] + s_method_coh_t[0][5][tt][qq])/nframes*invgamma/3.0);
               fprintf(fpn, "\n");
               fprintf(fpn, "%10g", norm(arr_qvec_faces[0][qq]) );
               fprintf(fpn, " %10g", (s_method_coh_t[0][0][tt][qq]+ s_method_coh_t[0][2][tt][qq] + s_method_coh_t[0][4][tt][qq])/nframes*invgamma/3.0);
               fprintf(fpn, "\n");
           }
           fprintf(fpn,"&\n");
           nplots++;
           fprintf(fpn, "@    s%d legend \"total theta=%g all faces\"\n",nplots,theta);
           fprintf(fpn, "@target G0.S%d\n",nplots);
           fprintf(fpn, "@type xy\n");
           for (qq = 0; qq< nbinq; qq++)
           {
               fprintf(fpn, "%10g", norm(arr_qvec_faces[1][qq]) );
               fprintf(fpn, " %10g", (s_method_t[0][1][tt][qq]+ s_method_t[0][3][tt][qq] + s_method_t[0][5][tt][qq])/nframes*invgamma/3.0);
               fprintf(fpn, "\n");
               fprintf(fpn, "%10g", norm(arr_qvec_faces[0][qq]) );
               fprintf(fpn, " %10g", (s_method_t[0][0][tt][qq]+ s_method_t[0][2][tt][qq] + s_method_t[0][4][tt][qq])/nframes*invgamma/3.0);
               fprintf(fpn, "\n");
           }
           fprintf(fpn,"&\n");
           nplots++;
           for (rr = 0; rr< nfaces; rr++)
           {
               fprintf(fpn, "@    s%d legend \"incoherent theta=%4g face index=%d \"\n",nplots,theta,rr);
               fprintf(fpn, "@target G0.S%d\n",nplots);
               fprintf(fpn, "@type xy\n");
               for (qq = 0; qq< nbinq; qq++)
               {
                   fprintf(fpn, "%10g", norm(arr_qvec_faces[rr][qq]) );
                   fprintf(fpn, " %10g", s_method_incoh_t[0][rr][tt][qq]/nframes*invgamma);
                   fprintf(fpn, "\n");
               }
               fprintf(fpn,"&\n");
               nplots++;
               fprintf(fpn, "@    s%d legend \"coherent theta=%4g face index=%d \"\n",nplots,theta,rr);
               fprintf(fpn, "@target G0.S%d\n",nplots);
               fprintf(fpn, "@type xy\n");
               for (qq = 0; qq< nbinq; qq++)
               {
                   fprintf(fpn, "%10g", norm(arr_qvec_faces[rr][qq]) );
                   fprintf(fpn, " %10g", s_method_coh_t[0][rr][tt][qq]/nframes*invgamma);
                   fprintf(fpn, "\n");
               }
               fprintf(fpn,"&\n");
               nplots++;
               fprintf(fpn, "@    s%d legend \"total theta=%4g face index=%d \"\n",nplots,theta,rr);
               fprintf(fpn, "@target G0.S%d\n",nplots);
               fprintf(fpn, "@type xy\n");
               for (qq = 0; qq< nbinq; qq++)
               {
                   fprintf(fpn, "%10g", norm(arr_qvec_faces[rr][qq]) );
                   fprintf(fpn, " %10g", s_method_t[0][rr][tt][qq]/nframes*invgamma);
                   fprintf(fpn, "\n");
               }
               fprintf(fpn,"&\n");
               nplots++;
           }
       }
    }
    do_view(oenv, fnSFACT, NULL);
    }

    if (fnBETACORR)
    {
       /* Calculate volume of sphere segments or length of circle segments */
       snew(inv_segvol, (nbin+1));
       prev_spherevol = 0;
       inv_segvol[0] = 1.0;
       normfac = 1.0/(nframes*invvol*isize0*(isize[0]-1));
       for (i = 1; (i < (nbin+1)); i++)
       {
           r = i*binwidth;
           spherevol = (4.0/3.0)*M_PI*r*r*r;
           segvol         = spherevol-prev_spherevol;
           inv_segvol[i]  = 1.0/segvol;
           prev_spherevol = spherevol;
       }
      
       sprintf(gtitle, "Non-linear optical scattering ");
       fpn = xvgropen(fnBETACORR, "hyperpolarizability spatial correlation", "r [nm]", "beta(0) beta(r)", oenv);
       sprintf(refgt, "%s", "");
       fprintf(fpn, "@type xy\n");

       for (i = 0; i < nbin+1; i++)
       {
           fprintf(fpn, "%10g %10g\n", i*binwidth, beta_corr[i]*normfac*inv_segvol[i] );
       }
       gmx_ffclose(fpn);
    }
    else if (fnFTBETACORR)
    {
       sprintf(gtitle, "FT of beta-beta correlation ");
       fpn = xvgropen(fnFTBETACORR, "hyperpolarizability spatial correlation", "r [nm]", "beta(0) beta(r)", oenv);
       sprintf(refgt, "%s", "");
       fprintf(fpn, "@type xy\n");
       for (qq = 0; qq < nbinq; qq++)
       {   
           fprintf(fpn, "%10g %10g\n", norm(arr_qvec_faces[0][qq]), ft_beta_corr[qq]/nframes );
       }
       gmx_ffclose(fpn);
    }

    for (g = 0; g < ng; g++)
    {
       sfree(s_method[g]);
       sfree(s_method_coh[g]);
       sfree(s_method_incoh[g]);
    }
    if(method[0] == 'm')
    {
        for (i = 0; rr < isize0; rr++)
        {
           for (rr = 0; rr < nfaces; rr++)
           {
               for (tt = 0; tt < nbintheta; tt++)
               {
                   sfree(mu_ind_t[i][rr][tt]);
               }
               sfree(mu_ind_t[i][rr]);
           }
           sfree(mu_ind_t[i]);
        }
        sfree(mu_ind_t);
    }

    sfree(s_method);
    sfree(s_method_coh);
    sfree(arr_qvec);
    sfree(cos_t);
    sfree(sin_t);
    sfree(temp_method);

    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            sfree(beta_mol[i][j]);
        }
        sfree(beta_mol[i]);
    }
    sfree(beta_mol);
    sfree(beta_corr);
    sfree(mu_ind_mols);
    if (kern[0] == 's')
    {
        initialize_free_quantities_on_grid(SKern_rho_O, ir, &grid_spacing,&grid_invspacing, FALSE, FALSE, &gridsize);
        initialize_free_quantities_on_grid(SKern_rho_H, ir, &grid_spacing,&grid_invspacing, FALSE, FALSE, &gridsize);
        initialize_free_quantities_on_grid(SKern_E, ir, &grid_spacing, &grid_invspacing, TRUE, FALSE, &gridsize);
        sfree(gridsize);
        fprintf(stderr,"quantities computed on global grid freed\n");
        for ( i = 0; i< SKern_E->ndataset; i++)
        {
           for ( aa = 0; aa < DIM; aa++)
           {
               for( bb = 0; bb < DIM; bb++)
               {
                   sfree(SKern_E->coeff[i][aa][bb]);
               }
               sfree(SKern_E->coeff[i][aa]);
           }
           sfree(SKern_E->coeff[i]);
        }
        sfree(SKern_E->coeff);
        for ( i = 0; i< SKern_rho_O->ndataset; i++)
        {
           for ( aa = 0; aa < DIM; aa++)
           {
               for( bb = 0; bb < DIM; bb++)
               {
                   sfree(SKern_rho_O->coeff[i][aa][bb]);
                   sfree(SKern_rho_H->coeff[i][aa][bb]);
               }
               sfree(SKern_rho_O->coeff[i][aa]);
               sfree(SKern_rho_H->coeff[i][aa]);
           }
           sfree(SKern_rho_O->coeff[i]);
           sfree(SKern_rho_H->coeff[i]);
        }
        sfree(SKern_rho_O->coeff);
        sfree(SKern_rho_H->coeff);

        sfree(SKern_E->grid);                 sfree(SKern_rho_O->grid);          sfree(SKern_rho_H->grid);       
        sfree(SKern_E->rotgrid);              sfree(SKern_rho_O->rotgrid);       sfree(SKern_rho_H->rotgrid);
        sfree(SKern_E->translgrid);           sfree(SKern_rho_O->translgrid);    sfree(SKern_rho_H->translgrid);
        sfree(SKern_E->meanquant);            sfree(SKern_rho_O->meanquant);     sfree(SKern_rho_H->meanquant);
        sfree(SKern_E->weights);              sfree(SKern_rho_O->weights);       sfree(SKern_rho_H->weights);     
        sfree(SKern_E->selfterm);              sfree(SKern_rho_O->selfterm);     sfree(SKern_rho_H->selfterm);            
        
   } 
}

void read_reference_mol(const char *fnREFMOL, rvec **xref)
{
      char refatomname[256], commentline[256];
      FILE *fp;
      int n_outputs, nats, i;
      real temp;

      fp = gmx_ffopen(fnREFMOL, "r");
      printf("reference molecule file\n");
      n_outputs = fscanf(fp,"%d ",&nats);
      printf("%d\n",nats);
      n_outputs = fscanf(fp,"%s ",commentline);
      printf("%s\n",commentline);
      for (i = 0; i< nats ; i++)
      {
          n_outputs = fscanf(fp,"%s ",refatomname);
          n_outputs = fscanf(fp,"%f ", &temp);
          (*xref)[i][XX] = temp*ANG2NM;
          n_outputs = fscanf(fp,"%f ", &temp);
          (*xref)[i][YY] = temp*ANG2NM;
          n_outputs = fscanf(fp,"%f ", &temp);
          (*xref)[i][ZZ] = temp*ANG2NM;
          printf("%s %f %f %f\n",refatomname, (*xref)[i][XX],(*xref)[i][YY],(*xref)[i][ZZ]);
      }
      fclose(fp);
}


void readMap(const char *fnMAP, t_Map *Map)
{
    FILE *fm;
    int i,j;
    int n_outputs;
     
    fm = gmx_ffopen(fnMAP, "r");
    for (i = 0; i<27; i++)
    {
        n_outputs = fscanf(fm,"%f ",&Map->beta_gas[i]);
    }

    // read coefficients
    for (i=0;i<27;i++) {  // 27 elements
      for (j=0;j<36;j++) {  // 212 components
        n_outputs = fscanf(fm,"%f ",&Map->D[i][j]);
      }
    }
  
    fclose(fm);
}

void readKern(const char *fnCOEFF, const char *fnGRD, const char *fnINPKRR, int kern_order, real std_dev, real kappa, rvec *xref, gmx_bool bAtomCenter, real **betamean, t_Kern *Kern)
{
    FILE *fk, *fg, *fp;
    int i, j, ch = 0, n_outputs ;
    int a, b, c;
    rvec dx;
    real temp;

    fp = (fnINPKRR == NULL) ? NULL : gmx_ffopen(fnINPKRR, "r");
    fk = gmx_ffopen(fnCOEFF, "r");
    fg = gmx_ffopen(fnGRD, "r");

    while(!feof(fk))
    {
      ch = fgetc(fk);
      if(ch == '\n')
      {
         Kern->ndataset ++;
      }
    }
    Kern->ndataset ++;
    rewind(fk);

    if (kappa != 0.0)
    {
       Kern->ndataset-= DIM*DIM*DIM;
    }

    while(!feof(fg))
    {
      ch = fgetc(fg);
      if(ch == '\n')
      {
         Kern->gridpoints ++;
      }
    }
    rewind(fg);
    Kern->gridpoints ++;
  
    //allocate all pointers needed for kernel operations
    if (fnINPKRR)
    {
       snew(Kern->krrinput, Kern->ndataset);
       snew(Kern->grid, Kern->gridpoints);
       snew(Kern->rotgrid, Kern->gridpoints);
       snew(Kern->translgrid, Kern->gridpoints);
       snew(Kern->coeff, Kern->ndataset);
       //allocate potential and krr coefficients
       for ( i = 0; i< Kern->ndataset; i++)
       {
          snew(Kern->krrinput[i], Kern->gridpoints);
          snew(Kern->coeff[i],DIM);
          for ( a = 0; a < DIM; a++)
          {
              snew(Kern->coeff[i][a],DIM);
              for( b = 0; b < DIM; b++)
              { 
                  snew(Kern->coeff[i][a][b],DIM);
              }
          }
       }
       //fprintf(stderr,"allocated coefficients\n");
       for (i = 0; i < Kern->ndataset; i++)
       {
           for ( j = 0; j < Kern->gridpoints; j++)
           {
               n_outputs = fscanf(fp, "%f", &Kern->krrinput[i][j]) ;
               //printf("vij %d %d %f\n",i,j,Kern->krrinput[i][j]);
           }
           for (a = 0; a < DIM; a++)
           {
               for (b = 0; b < DIM; b++)
               {
                   for (c = 0; c < DIM; c++)
                   {
                        n_outputs = fscanf(fk, "%f", &Kern->coeff[i][a][b][c]);
                        //fprintf(stderr,"coefficient %d %d %d %d %f\n",i,a,b,c, Kern->coeff[i][a][b][c]);
                   }
               }
           }
       }
       for ( j = 0; j < Kern->gridpoints; j++)
       {
           n_outputs = fscanf(fg, "%f", &Kern->grid[j][XX]);
           n_outputs = fscanf(fg, "%f", &Kern->grid[j][YY]);
           n_outputs = fscanf(fg, "%f", &Kern->grid[j][ZZ]);
           if (Kern->krrinput[0][j] == 0.0)
           {
              Kern->gridcenter = j;
           }
       }
    }
    else
    {
        snew(Kern->grid, Kern->gridpoints);
        snew(Kern->rotgrid, Kern->gridpoints);
        snew(Kern->translgrid, Kern->gridpoints);
        snew(Kern->meanquant,Kern->ndataset);
        snew(Kern->coeff, Kern->ndataset);
        for ( i = 0; i< Kern->ndataset; i++)
        {
           snew(Kern->coeff[i],DIM);
           for ( a = 0; a < DIM; a++)
           {
               snew(Kern->coeff[i][a],DIM);
               for( b = 0; b < DIM; b++)
               {
                   snew(Kern->coeff[i][a][b],DIM);
               }
           }
        }

        if (std_dev !=0.0)
        {
           snew(Kern->interp_quant_grid,Kern->gridpoints);
        }
        else
        {
           snew(Kern->vec_interp_quant_grid,Kern->gridpoints);
        }
    
        for (i = 0; i < Kern->ndataset; i++)
        {
           for (a = 0; a < DIM; a++)
           {
               for (b = 0; b < DIM; b++)
               {
                   for (c = 0; c < DIM; c++)
                   {
                           n_outputs = fscanf(fk, "%f", &Kern->coeff[i][a][b][c]);
                   }
               }
           }
           n_outputs = fscanf(fk, "%f", &Kern->meanquant[i]);
           printf(" %f\n", Kern->meanquant[i]);
        }
        if (kappa != 0.0)
        {
           snew(*betamean,DIM*DIM*DIM);
           for (i = 0; i < DIM*DIM*DIM; i++)
           {
               n_outputs = fscanf(fk, "%f", &temp);
             (*betamean)[i] = temp;
           }
        }
        for ( j = 0; j < Kern->gridpoints; j++)
        {
            n_outputs = fscanf(fg, "%f", &Kern->grid[j][XX]);
            n_outputs = fscanf(fg, "%f", &Kern->grid[j][YY]);
            n_outputs = fscanf(fg, "%f", &Kern->grid[j][ZZ]);
        }

        snew(Kern->weights,Kern->gridpoints);
        snew(Kern->selfterm,Kern->gridpoints);
        for (j = 0; j < Kern->gridpoints; j ++)
        {
           if (std_dev != 0.0)
           {
              Kern->weights[j] =  exp(-0.5*norm2(Kern->grid[j])/(4.0*std_dev*4.0*std_dev))  ;
              if (bAtomCenter)
              {
                 Kern->selfterm[j] = Kern->weights[j]*exp(-0.5*norm2(Kern->grid[j])/(std_dev*std_dev));
              }
              else
              {
                 rvec_sub(Kern->grid[j],xref[1],dx);
                 Kern->selfterm[j] = Kern->weights[j]*exp(-0.5*norm2(dx)/(std_dev*std_dev));
                 rvec_sub(Kern->grid[j],xref[2],dx); 
                 Kern->selfterm[j] += Kern->weights[j]*exp(-0.5*norm2(dx)/(std_dev*std_dev));
              }
           }
           else
           {
              Kern->weights[j] = 1.0;
/*              rvec_sub(Kern->grid[j],xref[0],dx);
              Kern->selfterm[j] = qref[0]*gmx_erf(norm(dx)*kappa)/norm(dx);
              rvec_sub(Kern->grid[j],xref[1],dx);
              Kern->selfterm[j] += qref[1]*gmx_erf(norm(dx)*kappa)/norm(dx);
              rvec_sub(Kern->grid[j],xref[2],dx);
              Kern->selfterm[j] += qref[2]*gmx_erf(norm(dx)*kappa)/norm(dx);
              rvec_sub(Kern->grid[j],xref[3],dx);
              Kern->selfterm[j] += qref[3]*gmx_erf(norm(dx)*kappa)/norm(dx);
*/
              Kern->selfterm[j] = 0.0;
           }
        }

    if (Kern->gridpoints*kern_order != Kern->ndataset && std_dev != 0.0)
    {
       gmx_fatal(FARGS,"check input files for the scalar kernel of the density.\n the number of lines in the weights file has to be Ngridpoints*kern_order");
    }
    else if (Kern->gridpoints*kern_order*3 != Kern->ndataset && std_dev == 0.0)
    {
       gmx_fatal(FARGS,"check input files for the scalar kernel of the electric field.\n the number of lines in the weights file has to be 3*Ngridpoints*kern_order");
    }

    }
    fclose(fk);
    if (fnINPKRR)
    {
      fclose(fp);
    }
    fclose(fg);
}

void identifyIon(t_topology *top,t_Ion *Ion,char *name)
{
  int i,n=0;

  for(i=0;i<top->atoms.nr;i++){
    if (!strcmp(*top->atoms.atomname[i],name)){
      Ion[n].atom[0]=i;
      Ion[n].q[0]=top->atoms.atom[i].q;
      n++;
    }
  }
  return;
}

int check_ion(t_topology *top,char *name) 
{
  int i,n=0;

  for (i=0;i<top->atoms.nr;i++) {
    if (!strcmp(*top->atoms.atomname[i],name)) {
      n+=1;
    }
  }
  return n;
}

void induced_second_order_fluct_dipole_fluct_beta(matrix cosdirmat,
                                                  const rvec pout, const rvec pin,
                                                  real ***betamol, real *mu_ind)
{
    int i, j, k;
    int l, m, n;
    //betal is the hyperpolarizability in lab frame
    real ***betal,  mu_temp = 0.0;
    //initialize beta
    snew(betal, DIM);
    *mu_ind = 0.0;
    for (i = 0; i < DIM; i++)
    {
        snew(betal[i], DIM);
        for (j = 0; j < DIM; j++)
        {
            snew(betal[i][j], DIM);
        }
    }
    //compute beta_lab as beta_lab(i,j,k)=sum(p,q,r) c_i,p c_j,q c_k,r * beta_mol_p,q,r
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
           for (k = 0; k < DIM; k++)
           {
               for (l = 0; l < DIM; l++)
               {
                    for (m = 0; m < DIM; m++)
                    {
                         for ( n = 0; n < DIM; n++)
                         {
                              betal[i][j][k] +=cosdirmat[i][l]*cosdirmat[j][m]*cosdirmat[k][n]*betamol[l][m][n];
                         }
                    }
               }
           }
        }
    }

    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
           for (k = 0; k < DIM; k++)
           {
               *mu_ind += betal[i][j][k]*pout[i]*pin[j]*pin[k];
       
           }
        }
    }

    //FREE beta lab
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            sfree(betal[i][j]);
        }
        sfree(betal[i]);
    }
    sfree(betal);


}



void induced_second_order_fluct_dipole(matrix cosdirmat, 
                                       const rvec pout, const rvec pin,
                                       real ***betamol, real *mu_ind)
{
    int i, j, k;
    int p, l, m;
    //betal is the hyperpolarizability in lab frame
    real ***betal,  mu_temp = 0.0;


    //initialize beta
    snew(betal, DIM);
    for (i = 0; i < DIM; i++)
    {
        snew(betal[i], DIM);
        for (j = 0; j < DIM; j++)
        {
            snew(betal[i][j], DIM);
        }
    }
    //compute beta_lab as beta_lab(i,j,k)=sum(p,q,r) c_i,p c_j,q c_k,r * beta_mol_p,q,r
    // the loops over p,q,r have been contracted using mathematica FullSimplify

    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
           for (k = 0; k < DIM; k++)
           {
               betal[i][j][k] = cosdirmat[i][ 0]*(cosdirmat[j][ 0]*(betamol[0][ 0][ 0]*cosdirmat[k][ 0] +
                                   betamol[0][ 0][ 1]*cosdirmat[k][ 1] + betamol[0][ 0][ 2]*cosdirmat[k][ 2]) +
                                   cosdirmat[j][ 1]*(betamol[0][ 1][ 0]*cosdirmat[k][ 0] + betamol[0][ 1][ 1]*cosdirmat[k][ 1] +
                                   betamol[0][ 1][ 2]*cosdirmat[k][ 2]) + cosdirmat[j][ 2]*(betamol[0][2][ 0]*cosdirmat[k][ 0] +
                                   betamol[0][ 2][ 1]*cosdirmat[k][ 1] + betamol[0][ 2][ 2]*cosdirmat[k][ 2])) +
                                   cosdirmat[i][ 1]*(cosdirmat[j][ 0]*(betamol[1][ 0][ 0]*cosdirmat[k][ 0] +
                                   betamol[1][ 0][ 1]*cosdirmat[k][ 1] + betamol[1][ 0][ 2]*cosdirmat[k][ 2]) +
                                   cosdirmat[j][ 1]*(betamol[1][ 1][ 0]*cosdirmat[k][ 0] + betamol[1][ 1][ 1]*cosdirmat[k][ 1] +
                                   betamol[1][ 1][ 2]*cosdirmat[k][ 2]) + cosdirmat[j][ 2]*(betamol[1][ 2][ 0]*cosdirmat[k][ 0] +
                                   betamol[1][ 2][ 1]*cosdirmat[k][ 1] + betamol[1][ 2][ 2]*cosdirmat[k][ 2])) +
                                   cosdirmat[i][ 2]*(cosdirmat[j][ 0]*(betamol[2][ 0][ 0]*cosdirmat[k][ 0] +
                                   betamol[2][ 0][ 1]*cosdirmat[k][ 1] + betamol[2][ 0][ 2]*cosdirmat[k][ 2]) +
                                   cosdirmat[j][ 1]*(betamol[2][  1][ 0]*cosdirmat[k][ 0] + betamol[2][ 1][ 1]*cosdirmat[k][ 1] +
                                   betamol[2][ 1][ 2]*cosdirmat[k][ 2]) + cosdirmat[j][ 2]*(betamol[2][ 2][ 0]*cosdirmat[k][ 0] +
                                   betamol[2][ 2][ 1]*cosdirmat[k][ 1] + betamol[2][ 2][ 2]*cosdirmat[k][ 2])) ;
           }
        }
    }

    //compute induced dipole component as mu = sum (i,j,k) (eout dot i) (ein dot j) (ein dot k) beta_lab_i,j,k
    // the loops over i,j,k have been contracted using mathematica FullSimplify
    *mu_ind =  betal[0][ 0][ 0]*pin[0]*pin[0]*pout[0] + betal[0][ 0][ 1]*pin[0]*pin[1]*pout[0] +
               betal[0][ 1][ 0]*pin[0]*pin[1]*pout[0] + betal[0][ 1][ 1]*pin[1]*pin[1]*pout[0] +
               betal[0][ 0][ 2]*pin[0]*pin[2]*pout[0] + betal[0][ 2][ 0]*pin[0]*pin[2]*pout[0] +
               betal[0][ 1][ 2]*pin[1]*pin[2]*pout[0] + betal[0][ 2][ 1]*pin[1]*pin[2]*pout[0] +
               betal[0][ 2][ 2]*pin[2]*pin[2]*pout[0] + betal[1][ 0][ 0]*pin[0]*pin[0]*pout[1] +
               betal[1][ 0][ 1]*pin[0]*pin[1]*pout[1] + betal[1][ 1][ 0]*pin[0]*pin[1]*pout[1] +
               betal[1][ 1][ 1]*pin[1]*pin[1]*pout[1] + betal[1][ 0][ 2]*pin[0]*pin[2]*pout[1] +
               betal[1][ 2][ 0]*pin[0]*pin[2]*pout[1] + betal[1][ 1][ 2]*pin[1]*pin[2]*pout[1] +
               betal[1][ 2][ 1]*pin[1]*pin[2]*pout[1] + betal[1][ 2][ 2]*pin[2]*pin[2]*pout[1] +
               (betal[2][ 0][ 0]*pin[0]*pin[0] + pin[1]*((betal[2][ 0][ 1] + betal[2][ 1][ 0])*pin[0] +
               betal[2][ 1][ 1]*pin[1]) + ((betal[2][ 0][ 2] + betal[2][ 2][ 0])*pin[0] +
               (betal[2][ 1][ 2] + betal[2][ 2][ 1])*pin[1])*pin[2] + betal[2][ 2][ 2]*pin[2]*pin[2])*pout[2] ;


/*    for (i = 0; i< DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
           for (k = 0; k < DIM; k++)
           {
               mu_temp += betal[i][j][k]*pout_temp[i]*pin_temp[j]*pin_temp[k];
           }
        }
    }
    *mu_ind = mu_temp;
*/

    //FREE beta lab
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            sfree(betal[i][j]);
        }
        sfree(betal[i]);
    }
    sfree(betal);
}


void calc_beta_skern( t_Kern *SKern_rho_O, t_Kern *SKern_rho_H, t_Kern *SKern_E, int kern_order, real *betamean, real ****betamol)
{
     int gr_ind, a, b, c, kern_ind;
     real feature_vec, feature_vec_x, feature_vec_y, feature_vec_z;
     int ind_ex, ind_ey, ind_ez, ind_rho;

     for (a = 0; a < DIM ; a++)
     {
         for (b = 0; b < DIM ; b++)
         {
            for (c = 0; c < DIM ; c++)
            {
               (*betamol)[a][b][c] = 0.0;
            }
         }
     }

     for (kern_ind = 0; kern_ind < kern_order; kern_ind ++)
     {
         for (gr_ind = 0; gr_ind < SKern_E->gridpoints; gr_ind++)
         {
            ind_ex = DIM*(gr_ind + SKern_E->gridpoints*kern_ind) + XX;
            ind_ey = ind_ex + YY;
            ind_ez = ind_ex + ZZ;
            feature_vec_x = pow(SKern_E->vec_interp_quant_grid[gr_ind][XX] - SKern_E->meanquant[ind_ex], kern_ind+1);
            feature_vec_y = pow(SKern_E->vec_interp_quant_grid[gr_ind][YY] - SKern_E->meanquant[ind_ey], kern_ind+1);
            feature_vec_z = pow(SKern_E->vec_interp_quant_grid[gr_ind][ZZ] - SKern_E->meanquant[ind_ez], kern_ind+1);
            
/*
            printf("index of coeff %d\n", ind_ex);
            printf("index of coeff %d\n", ind_ey);
            printf("index of coeff %d\n", ind_ez);
*/
//            printf("kern_ind %d gr_ind %d feature_vec_x %f\n",kern_ind, gr_ind, feature_vec_x);
//            printf("kern_ind %d gr_ind %d feature_vec_y %f\n",kern_ind, gr_ind, feature_vec_y);
//            printf("kern_ind %d gr_ind %d feature_vec_z %f\n",kern_ind, gr_ind, feature_vec_z);


            if (debug)
            {
            printf("electric_field= %f %f %f\n",feature_vec_x,feature_vec_y,feature_vec_z);

            printf("feature_vec_x %f \n",feature_vec_x);
            printf("feature_vec_y %f \n",feature_vec_y);
            printf("feature_vec_z %f \n",feature_vec_z);

            printf("predicted_vec_1 %f\n",SKern_E->coeff[ind_ex][0][0][0]*feature_vec_x);
            printf("predicted_vec_1 %f\n",SKern_E->coeff[ind_ey][0][0][0]*feature_vec_y);
            printf("predicted_vec_1 %f\n",SKern_E->coeff[ind_ez][0][0][0]*feature_vec_z);
            printf("predicted_vec_2 %f\n",SKern_E->coeff[ind_ex][2][2][2]*feature_vec_x);
            printf("predicted_vec_2 %f\n",SKern_E->coeff[ind_ey][2][2][2]*feature_vec_y);
            printf("predicted_vec_2 %f\n",SKern_E->coeff[ind_ez][2][2][2]*feature_vec_z);
            }

/*
            printf("coeff_vec_2 %f\n",SKern_E->coeff[ind_ex][2][2][2]);
            printf("coeff_vec_2 %f\n",SKern_E->coeff[ind_ey][2][2][2]);
            printf("coeff_vec_2 %f\n",SKern_E->coeff[ind_ez][2][2][2]);
*/

            for (a = 0; a < DIM ; a++)
            {
                for (b = 0; b < DIM ; b++)
                {
                   for (c = 0; c < DIM ; c++)
                   {
                      (*betamol)[a][b][c] += SKern_E->coeff[ind_ex][a][b][c]*feature_vec_x;
//                      printf("beta_%d_%d_%d %f coeff %f\n",a,b,c, (*betamol)[a][b][c],SKern_E->coeff[ind_ex][a][b][c] );
                      (*betamol)[a][b][c] += SKern_E->coeff[ind_ey][a][b][c]*feature_vec_y;
//                      printf("beta_%d_%d_%d %f coeff %f\n",a,b,c, (*betamol)[a][b][c],SKern_E->coeff[ind_ey][a][b][c] );
                      (*betamol)[a][b][c] += SKern_E->coeff[ind_ez][a][b][c]*feature_vec_z;
//                      printf("beta_%d_%d_%d %f coeff %f\n",a,b,c, (*betamol)[a][b][c],SKern_E->coeff[ind_ez][a][b][c] );

                   }
                }
            }
         }
         for (gr_ind = 0; gr_ind < SKern_rho_O->gridpoints; gr_ind++)
         {
            ind_rho = gr_ind + SKern_rho_O->gridpoints*kern_ind;
            feature_vec = pow(SKern_rho_O->interp_quant_grid[gr_ind] - SKern_rho_O->meanquant[ind_rho] , kern_ind+1);

//            printf("kern_ind %d gr_ind %d feature_vec_rho_O %f\n",kern_ind, gr_ind, feature_vec);

            if (debug)
            {
            printf("feature_vec %f\n",feature_vec);
            printf("predicted_vec_1 %f\n",SKern_rho_O->coeff[ind_rho][0][0][0]*feature_vec);
            printf("predicted_vec_2 %f\n",SKern_rho_O->coeff[ind_rho][2][2][2]*feature_vec);
            printf("coeff_vec_2 %f\n",SKern_rho_O->coeff[ind_rho][2][2][2]);
            }

            for (a = 0; a < DIM ; a++)
            {
                for (b = 0; b < DIM ; b++)
                {
                   for (c = 0; c < DIM ; c++)
                   {
                      (*betamol)[a][b][c] += SKern_rho_O->coeff[ind_rho][a][b][c]*feature_vec;
//                      printf("beta_%d_%d_%d %f coeff %f\n",a,b,c, (*betamol)[a][b][c],SKern_rho_O->coeff[ind_rho][a][b][c] );


                   }
                }
            }
         }
         for (gr_ind = 0; gr_ind < SKern_rho_H->gridpoints; gr_ind++)
         {
            ind_rho = gr_ind + SKern_rho_H->gridpoints*kern_ind;
            feature_vec = pow(SKern_rho_H->interp_quant_grid[gr_ind] - SKern_rho_H->meanquant[ind_rho] , kern_ind+1);

//            printf("kern_ind %d gr_ind %d feature_vec_rho_H %f\n",kern_ind, gr_ind, feature_vec);


            if (debug)
            {
            printf("feature_vec %f\n",feature_vec);
            printf("predicted_vec_1 %f\n",SKern_rho_H->coeff[ind_rho][0][0][0]*feature_vec);
            printf("predicted_vec_2 %f\n",SKern_rho_H->coeff[ind_rho][2][2][2]*feature_vec);
            printf("coeff_vec_2 %f\n",SKern_rho_H->coeff[ind_rho][2][2][2]);
            }

            for (a = 0; a < DIM ; a++)
            {
                for (b = 0; b < DIM ; b++)
                {
                   for (c = 0; c < DIM ; c++)
                   {
                      (*betamol)[a][b][c] += SKern_rho_H->coeff[ind_rho][a][b][c]*feature_vec;
//                       printf("beta_%d_%d_%d %f coeff %f\n",a,b,c, (*betamol)[a][b][c],SKern_rho_H->coeff[ind_rho][a][b][c] );
                   }
                }
            }
         }
     }
     for (a = 0; a < DIM ; a++)
     {
         for (b = 0; b < DIM ; b++)
         {
            for (c = 0; c < DIM ; c++)
            {
               (*betamol)[a][b][c] += betamean[c+DIM*b+DIM*DIM*a];
        
//                printf("beta_final_%d_%d_%d %f betamean %f\n",a,b,c, (*betamol)[a][b][c], betamean[c+DIM*b+DIM*DIM*a] );
            }
         }
     }
//     if (debug)
//     {
//        gmx_fatal(FARGS,"consider only one loop\n");
//     }   
     
}

void calc_beta_krr(t_Kern *Krr, t_pbc *pbc, t_topology *top, t_block *mols, int  *molindex[],
                 const int gind , const int isize0, rvec *x, rvec xcm_transl, const int imol, real rmax2, real electrostatic_cutoff ,real ****betamol)
{
   int    grid_ind, set_ind, j, m, indj, a, b, c;
   real   d2, delV, kernel_fn, *Vij, v_t, vreference, r_dist, el_cut2, sw_coeff;
   rvec   vecij, vecr, *vecs, *vecs_deb=NULL;
   matrix inverted_mat;

   //initialize molecular beta
   for (a = 0; a < DIM; a++)
   {
       for (b = 0; b < DIM; b++)
       {
          for (c = 0; c < DIM; c++)
          {
             (*betamol)[a][b][c] = 0.0 ;
          }
       }
   }
   //initialize potential at a point in the grid in the vicinity of molecule imol due to all molecules in a cutoff
   snew(Vij,Krr->gridpoints);
   snew(vecs,4);
   vreference = 0.0;
   v_t = 0.0;
/*   if (bEWALD)
   {
       for (j = 0; j < isize0; j++)
       {
           indj = mols->index[molindex[gind][j]];
           pbc_dx(pbc, xcm_transl, x[indj] ,vecij); 
           d2 = norm2(vecij);
           if (d2 < rmax2 && imol != indj)
           {
               for (m = 1; m < 4; m++)
               {
                  pbc_dx(pbc,x[indj+m],x[imol],vecs[m]);
               }

               for (grid_ind = 0; grid_ind < Krr->gridpoints; grid_ind++)
               {
                   for (m = 1; m < 4; m++)
                   {
                       rvec_sub(vecs[m], Krr->rotgrid[grid_ind],vecr);
                       r_dist = norm(vecr);
                       v_t = (top->atoms.atom[m].q)*exp(-invkappa2*r_dist*r_dist)/(core_term + r_dist) ;
                       Vij[grid_ind] += v_t;
                       vreference += v_t;
//                       printf("within cutoff potential %d %d %f %f %f\n", m, grid_ind, sqrt(d2) ,r_dist, Vij[grid_ind]);
                   }
               }
           }
       }
//       gmx_fatal(FARGS,"end\n");
   }

   if (bFADE)
   {
       el_cut2 = electrostatic_cutoff*electrostatic_cutoff;
       for (j = 0; j < isize0; j++)
       {
           indj = mols->index[molindex[gind][j]];
           pbc_dx(pbc, xcm_transl, x[indj] ,vecij);
           d2 = norm2(vecij);
           //fprintf(stderr,"distance squared %d %f\n",imol, d2);
           if (d2 < el_cut2 && indj != imol)
           {
               for (m = 1; m < 4; m++)
               {
                  pbc_dx(pbc,x[indj+m],x[imol],vecs[m]);
               }
               for (grid_ind = 0; grid_ind < Krr->gridpoints; grid_ind++)
               {
                   for (m = 1; m < 4; m++)
                   {
                       rvec_sub(vecs[m], Krr->rotgrid[grid_ind],vecr);
                       v_t = (top->atoms.atom[m].q)*invnorm(vecr) ;
                       Vij[grid_ind] += v_t;
                       vreference += v_t;
                       //printf("within cutoff potential %d %d %f %f %f\n",indj, m, sqrt(d2) ,r_dist, Vij[grid_ind]);
                   }
               }
           }
           else  if (d2 > el_cut2 && d2 < rmax2 )
           {
               r_dist = sqrt(d2);
               sw_coeff = cos((r_dist-electrostatic_cutoff)*inv_width);
               sw_coeff *= sw_coeff;
               for (m = 1; m < 4; m++)
               {
                  pbc_dx(pbc,x[indj+m],x[imol],vecs[m]);
               }
               for (grid_ind = 0; grid_ind < Krr->gridpoints; grid_ind++)
               {
                  for (m = 1; m < 4; m++)
                  {
                     rvec_sub(vecs[m], Krr->rotgrid[grid_ind],vecr);
                     v_t =  (top->atoms.atom[m].q)*sw_coeff*invnorm(vecr);
                     Vij[grid_ind] += v_t;
                     vreference += v_t;
                    //printf("above cutoff %d %d %f %f %f\n",indj, m, sqrt(d2), norm(vecr), Vij[grid_ind]);
                  }
               }
           }
       }
   }

   else
   {
*/
   el_cut2 = electrostatic_cutoff*electrostatic_cutoff;
   for (j = 0; j < isize0; j++)
   {
       indj = mols->index[molindex[gind][j]];
       pbc_dx(pbc, xcm_transl, x[indj] ,vecij);
       d2 = norm2(vecij);
       if (d2 < rmax2 && indj != imol)
       {
           for (m = 1; m < 4; m++)
           {
               pbc_dx(pbc,x[indj+m],x[imol],vecs[m]);
           }
           for (grid_ind = 0; grid_ind < Krr->gridpoints; grid_ind++)
           {
              for (m = 1; m < 4; m++)
              {
                  rvec_sub(vecs[m], Krr->rotgrid[grid_ind],vecr);
                  v_t = (top->atoms.atom[m].q)*invnorm(vecr);
                  Vij[grid_ind] += v_t;
                  vreference += v_t;
                  //printf("distance potential %d %d %f %f\n",indj, m, norm(vecr), Vij[grid_ind]);
                  //gmx_fatal(FARGS, "end\n");
              }
           }
       }
   }
/*
   }
*/
   vreference /= Krr->gridpoints;
   for (grid_ind = 0; grid_ind < Krr->gridpoints; grid_ind++)
   {
       Vij[grid_ind] -= vreference;
//       printf("v final imol grid_ind %d %d %f\n",imol,grid_ind,Vij[grid_ind]);
   }
//   gmx_fatal(FARGS,"end\n");
   for (set_ind = 0; set_ind < Krr->ndataset; set_ind++)
   {
        delV = 0.0;
        for (grid_ind = 0; grid_ind < Krr->gridpoints; grid_ind++)
        {
            delV -= sqr( (Vij[grid_ind]) - (Krr->krrinput[set_ind][grid_ind]) );
//            printf("del V %f %d %d %d\n",sqr( (Vij[grid_ind]) - (Krr->krrinput[set_ind][grid_ind]) ), grid_ind, set_ind, imol );
        }
        kernel_fn = exp(delV*(Krr->kerndev));
        //printf("kernel_fn %d %d %f %f\n", set_ind, imol, delV, kernel_fn);
        for (a = 0; a < DIM; a++)
        {
            for (b = 0; b < DIM; b++)
            {
                for (c = 0; c < DIM; c++)
                {
                   (*betamol)[a][b][c] += (Krr->coeff[set_ind][a][b][c])*kernel_fn;
//                   printf("betamol krr %d %d %d %d %f %f %f\n",a,b,c,set_ind,(Krr->coeff[set_ind][a][b][c])*kernel_fn, kernel_fn, (Krr->coeff[set_ind][a][b][c]));
                }
            }
        }
   }
//   for (a = 0; a < DIM; a++)
//   {
//       for (b = 0; b < DIM; b++)
//       {
//           for (c = 0; c < DIM; c++)
//           {
//               printf("betamol final %d %d %d %f\n",a,b,c,betamol[a][b][c]);
//           }
//       }
//   }
   sfree(Vij);
   sfree(vecs);
//  *betamol_addr = betamol;
}

void switch_fn(real r_dist, real electrostatic_cutoff, real rmax, real inv_width, real *sw_coeff)
{
   real temp;

   if (r_dist <= electrostatic_cutoff)
   {
      *sw_coeff = 1.0;
   }
   else if (r_dist >= rmax)
   {
     *sw_coeff = 0.0;
   }
   else
   {
      temp = cos((r_dist-electrostatic_cutoff)*inv_width);
     *sw_coeff = temp*temp ;
   }
}

void calc_beta_corr( t_pbc *pbc, t_block *mols, int  *molindex[],
                 const int gind , const int isize0,  int nbin, const real rmax2,  real invhbinw,
                 rvec *x, real *mu_ind_mols, real **beta_corr)
{

   rvec   dx;
   real   d2;
   real  *beta_t_corr;
   int    i, j, indj, imol, bin_ind;
   int    *count;
   
   snew(beta_t_corr,nbin+1);
   snew(count,nbin+1);

   for (i = 0; i < isize0 -1; i ++)
   {
      imol = mols->index[molindex[gind][i]];
      for (j = i+1 ; j < isize0 ; j++)
      {
          indj = mols->index[molindex[gind][j]];
          pbc_dx(pbc, x[imol], x[indj] ,dx);
          d2 = iprod(dx, dx);
          if ( d2 < rmax2 )
          {
             bin_ind = sqrt(d2)*invhbinw;
             beta_t_corr[bin_ind] += mu_ind_mols[i]*mu_ind_mols[j];
             count[bin_ind] += 2;
          }
      }
   }

   for ( i = 0; i < nbin +1; i++)
   {
       if (count[i] != 0)
       {
          (*beta_corr)[i] += beta_t_corr[i];
       }
   }
   sfree(beta_t_corr);
   sfree(count);
}

void calc_ft_beta_corr( t_pbc *pbc, t_block *mols, int  *molindex[],
                 const int gind , const int isize0,  int nbinq, rvec **arr_qvec_faces, real rmax2,  real invhbinw,
                 rvec *x, real *mu_ind_mols, real **ft_beta_corr)
{
   rvec   dx;
   real   d2, qr, dist, invdist;
   int    i, j, qq,indj, imol;
   real  *arr_qvec_norm, *invq;
   
   snew(arr_qvec_norm,nbinq);
   snew(invq,nbinq);
   for (qq = 0; qq < nbinq; qq ++)
   {
       arr_qvec_norm[qq] = norm(arr_qvec_faces[0][qq]);
       invq[qq] = 1.0/arr_qvec_norm[qq] ;
   }   

   for (i = 0; i < isize0 -1; i ++)
   {
      imol = mols->index[molindex[gind][i]];
      for (j = i+1 ; j < isize0 ; j++)
      {
          indj = mols->index[molindex[gind][j]];
          pbc_dx(pbc, x[imol], x[indj] ,dx);
          d2 = iprod(dx, dx);
          if ( d2 < rmax2 )
          {
             dist = sqrt(d2);
             invdist = 1.0/dist;
             for (qq = 0; qq < nbinq; qq++)
             {
                qr = arr_qvec_norm[qq]*dist;
                (*ft_beta_corr)[qq] += mu_ind_mols[i]*mu_ind_mols[j]*sin(qr)*invdist*invq[qq] ;
             }
          }
      }
   }
   sfree(invq);
   sfree(arr_qvec_norm);
}



void calc_efield_map(t_pbc *pbc,t_topology *top, t_block *mols, t_Ion *Cation, t_Ion *Anion, int  *molindex[], 
                 const int gind , const int isize0, const int ncations, const int nanions,
                 rvec *x, const int imol,  rvec xvec, rvec yvec,  rvec zvec, real electrostatic_cutoff2, real ***field_addr)
{
   int    j, m, n, indj;
   real  chg, d2, r, r2, ir, ir3, ir5;
   real   vx, vy, vz, chr3, chr5;
   real **Efield;
   rvec   vecij, vecr;

   snew(Efield,3);
   for (n = 0; n < 3; n++)
   {
       snew(Efield[n],9);
   }
   for (j = 0; j < isize0; j++)
   {
       indj = mols->index[molindex[gind][j]];
       pbc_dx(pbc, x[imol],x[indj],vecij);
       d2 = norm2(vecij);
       if (d2 < electrostatic_cutoff2 && indj != imol)
       {
          for (n = 0; n < 3; n++)
          {
              for (m = 1; m < 4; m++)
              {
                  pbc_dx(pbc,x[imol+n],x[indj+m],vecr);
                  svmul(nm2Bohr,vecr,vecr);
                  r=norm(vecr);
                  r2=r*r;
                  ir=1.0/r;
                  ir3=ir*ir*ir;
                  ir5=ir3*ir*ir;
                  chg = top->atoms.atom[m].q;
                  vx = iprod(vecr,xvec);
                  vy = iprod(vecr,yvec);
                  vz = iprod(vecr,zvec);
                  chr3 = chg*ir3;
                  chr5 = chg*ir5;
                  Efield[n][0] += vx*chr3 ;
                  Efield[n][1] += vy*chr3 ;
                  Efield[n][2] += vz*chr3 ;
                  // electric field gradient
                  Efield[n][3] += chr5*(r2 - 3*vx*vx );
                  Efield[n][4] += chr5*(r2 - 3*vy*vy );
                  Efield[n][5] += chr5*(r2 - 3*vz*vz );
                  Efield[n][6] -= chr5*( 3*vx*vy );
                  Efield[n][7] -= chr5*( 3*vy*vz );
                  Efield[n][8] -= chr5*( 3*vz*vx );
              }
          }
       }
   }
   for (j = 0; j < ncations; j++)
   {
       pbc_dx(pbc, x[imol],x[Cation[j].atom[0]],vecij);
       d2 = norm2(vecij);
       if (d2 < electrostatic_cutoff2 )
       {
          for (n = 0; n < 3; n++)
          {
             for (m = 0; m < 1; m++)
             {
                 pbc_dx(pbc,x[imol+n],x[Cation[j].atom[m]],vecr);
                 svmul(nm2Bohr,vecr,vecr);
                 r=norm(vecr);
                 r2=r*r;
                 ir=1.0/r;
                 ir3=ir*ir*ir;
                 ir5=ir3*ir*ir;
                 chg =Cation[j].q[m];
                 vx = iprod(vecr,xvec);
                 vy = iprod(vecr,yvec);
                 vz = iprod(vecr,zvec);
                 chr3 = chg*ir3;
                 chr5 = chg*ir5;
                 Efield[n][0] += vx*chr3 ;
                 Efield[n][1] += vy*chr3 ;
                 Efield[n][2] += vz*chr3 ;
                 // electric field gradient
                 Efield[n][3] += chr5*(r2 - 3*vx*vx );
                 Efield[n][4] += chr5*(r2 - 3*vy*vy );
                 Efield[n][5] += chr5*(r2 - 3*vz*vz );
                 Efield[n][6] -= chr5*( 3*vx*vy );
                 Efield[n][7] -= chr5*( 3*vy*vz );
                 Efield[n][8] -= chr5*( 3*vz*vx );
             }
          }
       }
   }
   for (j = 0; j < nanions; j++)
   {
       pbc_dx(pbc, x[imol],x[Anion[j].atom[0]],vecij);
       d2 = norm2(vecij);
       if (d2 < electrostatic_cutoff2)
       {
          for (n = 0; n < 3; n++)
          {
             for (m = 0; m < 1; m++)
             {
                 pbc_dx(pbc,x[imol+n],x[Anion[j].atom[m]],vecr);
                 //printf("anion position  %f %f %f\n" , x[Anion[j].atom[m]][0], x[Anion[j].atom[m]][1], x[Anion[j].atom[m]][2]);
                 //printf("vector distance %f %f %f\n", vecr[0], vecr[1], vecr[2]);
                 svmul(nm2Bohr,vecr,vecr);
                 r=norm(vecr);
                 r2=r*r;
                 ir=1.0/r;
                 ir3=ir*ir*ir;
                 ir5=ir3*ir*ir;
                 chg =Anion[j].q[m];
                 vx = iprod(vecr,xvec);
                 vy = iprod(vecr,yvec);
                 vz = iprod(vecr,zvec);
                 chr3 = chg*ir3;
                 chr5 = chg*ir5;
                 Efield[n][0] += vx*chr3 ;
                 Efield[n][1] += vy*chr3 ;
                 Efield[n][2] += vz*chr3 ;
                 // electric field gradient
                 Efield[n][3] += chr5*(r2 - 3*vx*vx );
                 Efield[n][4] += chr5*(r2 - 3*vy*vy );
                 Efield[n][5] += chr5*(r2 - 3*vz*vz );
                 Efield[n][6] -= chr5*( 3*vx*vy );
                 Efield[n][7] -= chr5*( 3*vy*vz );
                 Efield[n][8] -= chr5*( 3*vz*vx );
             }
          }
       }
   }
   *field_addr = Efield;
    //fprintf(stderr,"efield %f %f %f\n",Efield[0], Efield[1], Efield[2]);
}

void calc_beta_efield_map(t_Map *Map, real **Efield, real ***betamol, real ****betamol_addr)
{

  int i,index, n;
  int a,b,c;
  real delta;
  real betatemp[27] = {0.0};

  for (i=0;i<27;i++)
  {

    delta=0;
    index=0;

   // Electric field components
   for (n=0;n<3;n++) {
     delta+=Efield[0][n]*Map->D[i][index*3]+
            Efield[1][n]*Map->D[i][index*3+1]+
            Efield[2][n]*Map->D[i][index*3+2];
     index++;
   }

    // E^2 components   
    for (n=0;n<3;n++) {
      delta+=Efield[0][n]*Efield[0][n]*Map->D[i][index*3]+
             Efield[1][n]*Efield[1][n]*Map->D[i][index*3+1]+
             Efield[2][n]*Efield[2][n]*Map->D[i][index*3+2];
      index++;
    }

    // Gradient components
    for (n=3;n<9;n++) {
      delta+=Efield[0][n]*Map->D[i][index*3]+
             Efield[1][n]*Map->D[i][index*3+1]+
             Efield[2][n]*Map->D[i][index*3+2];
      index++;
    }

    // Calculate beta_liquid
    betatemp[i]=Map->beta_gas[i]+delta;
  }

  for (a = 0; a < DIM; a++)
  {
      for (b = 0; b < DIM; b++)
      {
          for (c = 0; c < DIM; c++)
          {
                  betamol[a][b][c] = betatemp[a*9+ b*3+ c];
          }
      }
  }
  //printf("bzxx %f bzyy %f bzzz %f\n",betamol[2][0][0], betamol[2][1][1], betamol[2][2][2]);
  //printf("betatemp_zxx %f betatemp_zyy %f betatemp_zzz %f\n",betatemp[2*9], betatemp[2*9 + 1*3 +1], betatemp[26]);
  *betamol_addr = betamol;

   for (n = 0; n < 3; n++)
   {
     sfree(Efield[n]);
   }
   sfree(Efield);
}


void calc_cosdirmat(const char *fnREFMOL, t_topology *top, int molsize,  int ind0, rvec *xref, rvec *xmol,
                    matrix *cosdirmat, matrix *invcosdirmat, rvec *xvec, rvec *yvec, rvec *zvec)
{
    rvec *xref_t=NULL ;
    real *w_rls =NULL;
    matrix cosdirmat_t;
    int i;

    rvec_add( xmol[1], xmol[2], *zvec);
//    printf("zvec %f %f %f\n",(*zvec)[0],(*zvec)[1],(*zvec)[2]);
    cprod(*zvec,xmol[1], *yvec);
//    printf("yvec %f %f %f\n",(*yvec)[0],(*yvec)[1],(*yvec)[2]);
    unitv(*yvec,*yvec);
    unitv(*zvec,*zvec);
//    printf("zvec %f %f %f\n",(*zvec)[0],(*zvec)[1],(*zvec)[2]);
//    printf("yvec %f %f %f\n",(*yvec)[0],(*yvec)[1],(*yvec)[2]);
    cprod(*yvec,*zvec,*xvec);

    if (!fnREFMOL)
    {
       (*cosdirmat)[0][0] = (*xvec)[0]; (*cosdirmat)[0][1] = (*yvec)[0]; (*cosdirmat)[0][2] = (*zvec)[0];
       (*cosdirmat)[1][0] = (*xvec)[1]; (*cosdirmat)[1][1] = (*yvec)[1]; (*cosdirmat)[1][2] = (*zvec)[1];
       (*cosdirmat)[2][0] = (*xvec)[2]; (*cosdirmat)[2][1] = (*yvec)[2]; (*cosdirmat)[2][2] = (*zvec)[2];
/*       printf("xmol[1] %f %f %f\n",xmol[1][XX],xmol[1][YY],xmol[1][ZZ]);
       printf("xmol[2] %f %f %f\n",xmol[2][XX],xmol[2][YY],xmol[2][ZZ]);
       printf("cosdirmat %13.7f %13.7f %13.7f\n",(*cosdirmat)[0][0],(*cosdirmat)[0][1],(*cosdirmat)[0][2]);
       printf("cosdirmat %13.7f %13.7f %13.7f\n",(*cosdirmat)[1][0],(*cosdirmat)[1][1],(*cosdirmat)[1][2]);
       printf("cosdirmat %13.7f %13.7f %13.7f\n",(*cosdirmat)[2][0],(*cosdirmat)[2][1],(*cosdirmat)[2][2]);
*/
       m_inv(*cosdirmat,*invcosdirmat);
    }
    else
    {
      snew(w_rls,molsize);
      snew(xref_t,molsize);
      for (i = 0; i < molsize; i++)
      {
          w_rls[i] = top->atoms.atom[ind0+i].m;
          copy_rvec(xref[i],xref_t[i]);
      }

      reset_x_ndim(DIM, molsize, NULL , molsize, NULL, xref_t, w_rls);
      reset_x_ndim(DIM, molsize, NULL, molsize, NULL, xmol, w_rls);
      calc_fit_R(DIM, molsize, w_rls, xref_t, xmol, cosdirmat_t);
      copy_mat(cosdirmat_t, *invcosdirmat);
      m_inv(cosdirmat_t, *cosdirmat);

/*
      printf("cosdirmat %13.7f %13.7f %13.7f\n",(*cosdirmat)[0][0],(*cosdirmat)[0][1],(*cosdirmat)[0][2]);
      printf("cosdirmat %13.7f %13.7f %13.7f\n",(*cosdirmat)[1][0],(*cosdirmat)[1][1],(*cosdirmat)[1][2]);
      printf("cosdirmat %13.7f %13.7f %13.7f\n",(*cosdirmat)[2][0],(*cosdirmat)[2][1],(*cosdirmat)[2][2]);
*/

      sfree(w_rls);
      sfree(xref_t);
    }
/*
    printf("invcosdirmat %13.7f %13.7f %13.7f\n",(*invcosdirmat)[0][0],(*invcosdirmat)[0][1],(*invcosdirmat)[0][2]);
    printf("invcosdirmat %13.7f %13.7f %13.7f\n",(*invcosdirmat)[1][0],(*invcosdirmat)[1][1],(*invcosdirmat)[1][2]);
    printf("invcosdirmat %13.7f %13.7f %13.7f\n",(*invcosdirmat)[2][0],(*invcosdirmat)[2][1],(*invcosdirmat)[2][2]);
*/
}

void rotate_local_grid_points(t_Kern *SKern_rho_O, t_Kern *SKern_rho_H, t_Kern *SKern_E, int ePBC, matrix box, matrix cosdirmat,rvec xi)
{
  int gr_ind;
   
  for (gr_ind =0; gr_ind < SKern_rho_O->gridpoints; gr_ind++)
  {
     mvmul(cosdirmat,SKern_rho_O->grid[gr_ind],SKern_rho_O->rotgrid[gr_ind]);
     rvec_add(xi,SKern_rho_O->rotgrid[gr_ind],SKern_rho_O->translgrid[gr_ind]);
     put_atoms_in_box(ePBC, box, 1, &SKern_rho_O->translgrid[gr_ind]);

     copy_rvec(SKern_rho_O->translgrid[gr_ind],SKern_rho_H->translgrid[gr_ind]);
  }

  for (gr_ind = 0; gr_ind < SKern_E->gridpoints; gr_ind++)
  {
     mvmul(cosdirmat,SKern_E->grid[gr_ind],SKern_E->rotgrid[gr_ind]);    
//     rvec_add(xi,SKern_E->grid[gr_ind],SKern_E->translgrid[gr_ind]);
     rvec_add(xi,SKern_E->rotgrid[gr_ind],SKern_E->translgrid[gr_ind]);
     put_atoms_in_box(ePBC, box, 1, &SKern_E->translgrid[gr_ind]);
//     printf("X %f %f %f\n",10.0*SKern_E->translgrid[gr_ind][XX], 10.0*SKern_E->translgrid[gr_ind][YY], 10.0*SKern_E->translgrid[gr_ind][ZZ]);
  }
}

  
void initialize_free_quantities_on_grid(t_Kern *Kern, t_inputrec *ir, rvec *grid_spacing, rvec *grid_invspacing, gmx_bool bEFIELD,  gmx_bool bALLOC, int **gridsize)
{
     int ix, iy, iz;

     if (bALLOC)
     {
        (*grid_invspacing)[XX] = 1.0/ir->fourier_spacing;
        (*grid_invspacing)[YY] = 1.0/ir->fourier_spacing;
        (*grid_invspacing)[ZZ] = 1.0/ir->fourier_spacing;
   
        (*grid_spacing)[XX] = ir->fourier_spacing;
        (*grid_spacing)[YY] = ir->fourier_spacing;
        (*grid_spacing)[ZZ] = ir->fourier_spacing;
   
        srenew(*gridsize,DIM);
        (*gridsize)[XX] = ir->nkx; (*gridsize)[YY] = ir->nky; (*gridsize)[ZZ] = ir->nkz;     
   
        if (!bEFIELD)
        {
           snew(Kern->quantity_on_grid,ir->nkx);
           for (ix = 0; ix < ir->nkx; ix++)
           {
               snew(Kern->quantity_on_grid[ix], ir->nky);
               for (iy = 0; iy < ir->nky ; iy ++)
               {
                  snew(Kern->quantity_on_grid[ix][iy], ir->nkz);
               }
           }
        }
        else
        {
           snew(Kern->quantity_on_grid_x,ir->nkx);
           snew(Kern->quantity_on_grid_y,ir->nkx);
           snew(Kern->quantity_on_grid_z,ir->nkx);
           for (ix = 0; ix < ir->nkx; ix++)
           {
               snew(Kern->quantity_on_grid_x[ix], ir->nky);
               snew(Kern->quantity_on_grid_y[ix], ir->nky);
               snew(Kern->quantity_on_grid_z[ix], ir->nky);
               for (iy = 0; iy < ir->nky ; iy ++)
               {
                  snew(Kern->quantity_on_grid_x[ix][iy], ir->nkz);
                  snew(Kern->quantity_on_grid_y[ix][iy], ir->nkz);
                  snew(Kern->quantity_on_grid_z[ix][iy], ir->nkz);
               }
           }     
        }
     }
     else
     {
        if (!bEFIELD)
        {
           for (ix = 0; ix < ir->nkx; ix++)
           {
               for (iy = 0; iy < ir->nky ; iy ++)
               {
                  sfree(Kern->quantity_on_grid[ix][iy]);
               }
               sfree(Kern->quantity_on_grid[ix]);
           }
           sfree(Kern->quantity_on_grid);
        }
        else
        {
           for (ix = 0; ix < ir->nkx; ix++)
           {
               for (iy = 0; iy < ir->nky ; iy ++)
               {
                  sfree(Kern->quantity_on_grid_x[ix][iy]);
                  sfree(Kern->quantity_on_grid_y[ix][iy]);
                  sfree(Kern->quantity_on_grid_z[ix][iy]);
               }
               sfree(Kern->quantity_on_grid_x[ix]);
               sfree(Kern->quantity_on_grid_y[ix]);
               sfree(Kern->quantity_on_grid_z[ix]);
           }
           sfree(Kern->quantity_on_grid_x);
           sfree(Kern->quantity_on_grid_y);
           sfree(Kern->quantity_on_grid_z);
        }
     }
}

void calc_efield_correction(t_Kern *Kern, t_inputrec *ir, t_topology *top, t_pbc *pbc, 
                          matrix box, real invvol, t_block *mols, int  *molindex[],
                         int *chged_atom_indexes, int n_chged_atoms, int *grid, rvec grid_spacing, rvec grid_invspacing,
                         rvec *x, int isize0, real *sigma_vals,real dxcut, real dxcut2, int *gridsize)
{
	// Calculate (in real space) the correction to the reciprocal-space part of the electric field,
	// using a different cutoff radius.

	int ix,iy,iz,m,i,j,ind0,n,size_near2;
	int ind_x, ind_y, ind_z, mx;
	rvec dx,xi;
	real charge,dx2,ef0,invdx2,dx2s,dx2b,dxs,dxb,invdx,scfc;
	int *bin_ind0;
  int **relevant_grid_points,*half_size_grid_points,*size_nearest_grid_points;

	scfc = 2.0 / sqrt(M_PI);

	// Note: sigma_vals should be a 4-d array:
	// 0: kappa2
	// 1: kappa
	// 2: kappa2^2
	// 3: kappa^2
	// kappa is the value that we use for PME calculations, and kappa2 is the value that we would like
	// to use for the output electric field.

	// Firstly, check whether or not a user-defined value has been given to the "small" (i.e., target)
	// sigma. If not (i.e., it's still at its default value of -1.0), then we don't do this part of the
	// program.
	if (sigma_vals[0] > 0.0)
	{
	
		// We don't free and re-initialize Kern->quantity_on_grid_XYZ,
		// because we're going to append to our previous values.
	
		snew(bin_ind0,DIM);
		snew(half_size_grid_points,DIM);
		snew(size_nearest_grid_points,DIM);

		mx = 0;
		for (ix=0;ix<DIM;ix++)
		{
			half_size_grid_points[ix] = floor(dxcut*grid_invspacing[ix]) + 1;
			size_nearest_grid_points[ix] = 2 * half_size_grid_points[ix];
			if (size_nearest_grid_points[ix]>mx){mx = size_nearest_grid_points[ix];}
		}
	
		snew(relevant_grid_points,mx);
		for (i=0;i<mx;i++)
		{
			snew(relevant_grid_points[i],DIM);
		}

		// Loop over molecules.
		for (n=0;n<isize0;n++)
		{
			for (m = 0;m<n_chged_atoms;m++)
			{
				ind0 = mols->index[molindex[0][n]] + chged_atom_indexes[m] ;
				copy_rvec(x[ind0],xi);
				charge = top->atoms.atom[ind0].q;

				// Work out what the closest grid point is to this molecule.
				for (ix=0;ix<DIM;ix++)
				{
					bin_ind0[ix] = roundf(xi[ix]*grid_invspacing[ix]);
					if (bin_ind0[ix] == gridsize[ix]){bin_ind0[ix]=0;}
					if (bin_ind0[ix] == -1){bin_ind0[ix]=gridsize[ix]-1;}
				}
	
				// Now find out which grid points should be checked.
				for (ix=0;ix<DIM;ix++)
				{
					for (j=0;j<size_nearest_grid_points[ix];j++)
					{
						relevant_grid_points[j][ix] = j - half_size_grid_points[ix] + bin_ind0[ix];
						if (relevant_grid_points[j][ix] >= gridsize[ix]){relevant_grid_points[j][ix] -= gridsize[ix];}
						if (relevant_grid_points[j][ix] <  0){relevant_grid_points[j][ix] += gridsize[ix];}
					}
				}
	
				for (ix=0;ix < size_nearest_grid_points[XX];ix++)
				{
					ind_x = relevant_grid_points[ix][XX];
					Kern->rspace_grid[XX] = ind_x * grid_spacing[XX];
					for (iy=0;iy < size_nearest_grid_points[YY];iy++)
					{
						ind_y = relevant_grid_points[iy][YY];
						Kern->rspace_grid[YY] = ind_y * grid_spacing[YY];
						for (iz=0;iz < size_nearest_grid_points[ZZ];iz++)
						{
							ind_z = relevant_grid_points[iz][ZZ];
							Kern->rspace_grid[ZZ] = ind_z * grid_spacing[ZZ];
							pbc_dx(pbc,xi,Kern->rspace_grid,dx);
							dx2 = norm2(dx);
							if (dx2<=dxcut2)
							{
								// Calculate electric field correction terms.
								invdx2 = 1.0/dx2;
								dx2s = dx2 * sigma_vals[2];
								dx2b = dx2 * sigma_vals[3];
								ef0 = sigma_vals[1]*exp(-1.0*dx2b) - sigma_vals[0]*exp(-1.0*dx2s);
								ef0 *= scfc;
								invdx = sqrt(invdx2);
								dxs = sqrt(dx2s);
								dxb = sqrt(dx2b);
								ef0 += invdx*( gmx_erf(dxs) - gmx_erf(dxb));
								ef0 *= charge*invdx2;
								Kern->quantity_on_grid_x[ind_x][ind_y][ind_z] += ef0 * dx[XX];
								Kern->quantity_on_grid_y[ind_x][ind_y][ind_z] += ef0 * dx[YY];
								Kern->quantity_on_grid_z[ind_x][ind_y][ind_z] += ef0 * dx[ZZ];
							}
						}
					}
				}
			}
		}

	}

}

void calc_dens_on_grid(t_Kern *Kern, t_inputrec *ir, t_pbc *pbc,
                       t_block *mols, int  *molindex[], int atom_id0, int nspecies ,int isize0, rvec *x,
                       real std_dev_dens, real inv_std_dev_dens, rvec grid_invspacing, int *gridsize, rvec grid_spacing)
{
  int   m, i, j, ix, iy, iz, ind0;
  int   ind_x, ind_y, ind_z;
  rvec  dx;
  real  inv_std;
  int  *bin_ind0;
  int  *bin_ind_std;
  int  *bin_indmax;
  int  *bin_indmin;
//  int  *rspace_gridsize;
  int   size_nearest_grid_points, half_size_grid_points;
//  real  ***dens_deb;
  int **relevant_grid_points;

  size_nearest_grid_points = roundf(std_dev_dens*12.0*grid_invspacing[XX]);
  half_size_grid_points = roundf(size_nearest_grid_points*0.5);
  snew(bin_ind0,DIM);
  snew(bin_ind_std,DIM);
  snew(bin_indmax,DIM);
  snew(bin_indmin,DIM);

//  rspace_gridsize[XX] = ir->nkx  ; rspace_gridsize[YY] = ir->nky  ; rspace_gridsize[ZZ] = ir->nkz  ;

  snew(relevant_grid_points,size_nearest_grid_points);
  for (i = 0; i < size_nearest_grid_points; i++)
  {
      snew(relevant_grid_points[i],DIM);
  }


  for (ix = 0; ix < ir->nkx; ix++)
  {
      for (iy = 0; iy < ir->nky ; iy ++)
      {
        sfree(Kern->quantity_on_grid[ix][iy]);
      }
      sfree(Kern->quantity_on_grid[ix]);
  }
  sfree(Kern->quantity_on_grid);

  snew(Kern->quantity_on_grid,ir->nkx);
  //snew(dens_deb,ir->nkx);
  for (ix = 0; ix < ir->nkx; ix++)
  {
      snew(Kern->quantity_on_grid[ix], ir->nky);
  //    snew(dens_deb[ix],ir->nky);
      for (iy = 0; iy < ir->nky ; iy ++)
      {
        snew(Kern->quantity_on_grid[ix][iy], ir->nkz);
   //     snew(dens_deb[ix][iy], ir->nkz);
      }
  }

  //fprintf(stderr,"start looping over atoms\n");
  for (i = 0; i < isize0; i++)
  {
      ind0 = mols->index[molindex[0][i]];
      for (m = atom_id0; m < nspecies + atom_id0 ; m++)
      {
         for (ix = 0; ix < DIM; ix ++)
         {
            bin_ind0[ix] = roundf((x[ind0+m][ix] /*- 0.5*grid_spacing[ix]*/)*grid_invspacing[ix] );
            if (bin_ind0[ix] == gridsize[ix])
            {
               bin_ind0[ix] = 0;
            }
            if (bin_ind0[ix] == -1)
            {
               bin_ind0[ix] = gridsize[ix] -1 ;
            }
         }
         for ( j = 0; j < size_nearest_grid_points; j ++)
         {
            for (ix = 0; ix < DIM; ix ++)
            {
               relevant_grid_points[j][ix] = j - half_size_grid_points + bin_ind0[ix];
       	       if (relevant_grid_points[j][ix] >= gridsize[ix])
               {
                  relevant_grid_points[j][ix] -= gridsize[ix];
               }
               if (relevant_grid_points[j][ix] < 0)
               {
                  relevant_grid_points[j][ix] += gridsize[ix];
               }
            }
         } 
   
         for (ix = 0; ix < size_nearest_grid_points; ix++)
         {
             ind_x = relevant_grid_points[ix][XX];
             Kern->rspace_grid[XX]=ind_x*grid_spacing[XX];
             for (iy = 0; iy < size_nearest_grid_points ; iy ++)
             {
                ind_y = relevant_grid_points[iy][YY];
                Kern->rspace_grid[YY]=ind_y*grid_spacing[YY];
                for (iz = 0; iz < size_nearest_grid_points ; iz ++)
                {
                   ind_z = relevant_grid_points[iz][ZZ];
                   Kern->rspace_grid[ZZ]=ind_z*grid_spacing[ZZ];  
                   pbc_dx(pbc,x[ind0+m],Kern->rspace_grid,dx);
                   Kern->quantity_on_grid[ind_x][ind_y][ind_z] += exp(-inv_std_dev_dens*norm2(dx)) ;
                }   
             }
         }
      }   
  }

  sfree(bin_ind0);
  sfree(bin_ind_std);
  sfree(bin_indmax);
  sfree(bin_indmin);

  for (i = 0; i < size_nearest_grid_points; i++)
  {
      sfree(relevant_grid_points[i]);
  }
  sfree(relevant_grid_points);
}

void trilinear_interpolation_kern(t_Kern *Kern, t_inputrec *ir, t_pbc *pbc, rvec xi, rvec grid_invspacing, rvec grid_spacing)
{
     int i, ix, iy, iz, ind0;
     int bin_indx0, bin_indy0, bin_indz0, bin_indx1, bin_indy1, bin_indz1;
     real xd, yd, zd, sel_dist;
     rvec delx, delxdeb;

     sfree(Kern->interp_quant_grid);
     snew(Kern->interp_quant_grid,Kern->gridpoints);

     for (i = 0; i < Kern->gridpoints; i++)
     {
         bin_indx1 = ceil((Kern->translgrid[i][XX] )*grid_invspacing[XX] );
         bin_indy1 = ceil((Kern->translgrid[i][YY] )*grid_invspacing[YY] );
         bin_indz1 = ceil((Kern->translgrid[i][ZZ] )*grid_invspacing[ZZ] );
         if (bin_indx1 > ir->nkx-1)
         {
            bin_indx1 = 0;
         }
         if (bin_indz1 > ir->nkz-1)
         {
            bin_indz1 = 0;
         }
         if (bin_indy1 > ir->nky-1)
         {
            bin_indy1 = 0;
         }

         bin_indx0 = (bin_indx1 > 0 ) ? bin_indx1  -1 : ir->nkx-1 ; 
         bin_indy0 = (bin_indy1 > 0 ) ? bin_indy1  -1 : ir->nky-1 ;
         bin_indz0 = (bin_indz1 > 0 ) ? bin_indz1  -1 : ir->nkz-1 ;

         Kern->rspace_grid[XX]=grid_spacing[XX]*bin_indx0;
         Kern->rspace_grid[YY]=grid_spacing[YY]*bin_indy0;
         Kern->rspace_grid[ZZ]=grid_spacing[ZZ]*bin_indz0;

         pbc_dx(pbc,Kern->translgrid[i],Kern->rspace_grid,delx);

         xd = fabs(delx[XX])*grid_invspacing[XX];
         yd = fabs(delx[YY])*grid_invspacing[YY];
         zd = fabs(delx[ZZ])*grid_invspacing[ZZ];

         pbc_dx(pbc,xi,Kern->translgrid[i],delx);

         Kern->interp_quant_grid[i] = Kern->quantity_on_grid[bin_indx0][bin_indy0][bin_indz0]*(1.0 - xd)*(1.0 - zd)*(1 - yd) +
                                     Kern->quantity_on_grid[bin_indx1][bin_indy0][bin_indz0]*xd*(1 - yd)*(1 - zd) +
                                     Kern->quantity_on_grid[bin_indx0][bin_indy1][bin_indz0]*(1.0 - xd)*yd*(1 - zd) +
                                     Kern->quantity_on_grid[bin_indx0][bin_indy0][bin_indz1]*(1.0 - xd)*(1 - yd)*zd + 
                                     Kern->quantity_on_grid[bin_indx1][bin_indy0][bin_indz1]*xd*(1 - yd)*zd +
                                     Kern->quantity_on_grid[bin_indx0][bin_indy1][bin_indz1]*(1 - xd)*yd*zd +
                                     Kern->quantity_on_grid[bin_indx1][bin_indy1][bin_indz0]*xd*yd*(1 - zd) +
                                     Kern->quantity_on_grid[bin_indx1][bin_indy1][bin_indz1]*xd*yd*zd;
         Kern->interp_quant_grid[i] *= Kern->weights[i];
//         Kern->interp_quant_grid[i] -= Kern->selfterm[i];
     }
}


void bspline_efield(t_Kern *Kern, t_inputrec *ir, t_pbc *pbc, matrix invcosdirmat, rvec xi, rvec grid_invspacing, int interp_order, rvec grid_spacing, rvec Emean)
{
     int i, a, b, c, kk1, kk2, kk3;
     int bin_indx0, bin_indy0, bin_indz0, bin_indx1, bin_indy1, bin_indz1;
     real ***Ematr_x, ***Ematr_y, ***Ematr_z;
     rvec bsplcoeff, vec_t, du, Efield_temp;

     sfree(Kern->vec_interp_quant_grid);
     snew(Kern->vec_interp_quant_grid,Kern->gridpoints);

     snew(Ematr_x,ir->nkx);
     snew(Ematr_y,ir->nkx);
     snew(Ematr_z,ir->nkx);

     for (a = 0; a < ir->nkx; a++)
     {
         snew(Ematr_x[a],ir->nky);
         snew(Ematr_y[a],ir->nky);
         snew(Ematr_z[a],ir->nky);
         for (b = 0; b < ir->nky; b++)
         {
             snew(Ematr_x[a][b],ir->nkz);
             snew(Ematr_y[a][b],ir->nkz);
             snew(Ematr_z[a][b],ir->nkz);
         }
     }

     for (i = 0; i < Kern->gridpoints; i++)
     {

         bin_indx0 = floor((Kern->translgrid[i][XX] )*grid_invspacing[XX] );
         bin_indy0 = floor((Kern->translgrid[i][YY] )*grid_invspacing[YY] );
         bin_indz0 = floor((Kern->translgrid[i][ZZ] )*grid_invspacing[ZZ] );

         if (bin_indx0 == ir->nkx  )
         {
            bin_indx0 = 0;
         }
         if (bin_indy0 == ir->nky )
         {
            bin_indy0 = 0;
         }
         if (bin_indz0 == ir->nkz )
         {
            bin_indz0 = 0;
         }

/*         fprintf(stderr,"binx %d biny %d binz %d NEW i %d\n",bin_indx0, bin_indy0, bin_indz0,i);


         bin_indx1 = ceil((Kern->translgrid[i][XX] )*grid_invspacing[XX] );
         bin_indy1 = ceil((Kern->translgrid[i][YY] )*grid_invspacing[YY] );
         bin_indz1 = ceil((Kern->translgrid[i][ZZ] )*grid_invspacing[ZZ] );

         if (bin_indx1 > ir->nkx-1)
         {
            bin_indx1 = 0;
         }
         if (bin_indz1 > ir->nkz-1)
         {
            bin_indz1 = 0;
         }
         if (bin_indy1 > ir->nky-1)
         {
            bin_indy1 = 0;
         }

         bin_indx0 = (bin_indx1 > 0 ) ? bin_indx1  -1 : ir->nkx-1 ;
         bin_indy0 = (bin_indy1 > 0 ) ? bin_indy1  -1 : ir->nky-1 ;
         bin_indz0 = (bin_indz1 > 0 ) ? bin_indz1  -1 : ir->nkz-1 ;
         
         fprintf(stderr,"binx %d biny %d binz %d\n",bin_indx0, bin_indy0, bin_indz0);

*/
         vec_t[XX] = 0.0;
         vec_t[YY] = 0.0;
         vec_t[ZZ] = 0.0;
       
         for (a = bin_indx0 - interp_order +1; a <= bin_indx0 + interp_order; a++)
         {
             du[XX] = Kern->translgrid[i][XX]*grid_invspacing[XX] -a;
             kk1 = index_wrap(a,ir->nkx);
             bsplcoeff[XX] = Bspline(du[XX],interp_order);
             for (b = bin_indy0 - interp_order + 1; b <= bin_indy0 + interp_order; b++)
             {
                 du[YY] = Kern->translgrid[i][YY]*grid_invspacing[YY]  -b;
                 kk2 = index_wrap(b,ir->nky);
                 bsplcoeff[YY] = Bspline(du[YY],interp_order);
                 for (c = bin_indz0 - interp_order + 1; c <= bin_indz0 + interp_order; c++)
                 {
                     du[ZZ] = Kern->translgrid[i][ZZ]*grid_invspacing[ZZ] -c;
                     kk3 = index_wrap(c,ir->nkz);
                     bsplcoeff[ZZ] = Bspline(du[ZZ],interp_order);
                     
                     vec_t[XX] += bsplcoeff[ZZ]*bsplcoeff[YY]*bsplcoeff[XX]*(Kern->quantity_on_grid_x[kk1][kk2][kk3]-Emean[XX]);
                     vec_t[YY] += bsplcoeff[ZZ]*bsplcoeff[YY]*bsplcoeff[XX]*(Kern->quantity_on_grid_y[kk1][kk2][kk3]-Emean[YY]);
                     vec_t[ZZ] += bsplcoeff[ZZ]*bsplcoeff[YY]*bsplcoeff[XX]*(Kern->quantity_on_grid_z[kk1][kk2][kk3]-Emean[ZZ]);
//                     printf("vec_t %d %d %d %f %f %f\n",a,b,c,vec_t[XX],vec_t[YY],vec_t[ZZ]);
                     
                 }
             }
         }
	 
         /*for (a = bin_indx0 - interp_order +1; a <= bin_indx0 + interp_order; a++)
         {
                kk1 = index_wrap(a,ir->nkx);
		for (b = bin_indy0 - interp_order + 1; b <= bin_indy0 + interp_order; b++) 
                {
                        kk2 = index_wrap(b,ir->nky);
			for (c = bin_indz0 - interp_order + 1; c <= bin_indz0 + interp_order; c++)
                        {
                                kk3 = index_wrap(c,ir->nkz);
				vec_t[XX] += Ematr_x[kk1][kk2][kk3];
				vec_t[YY] += Ematr_y[kk1][kk2][kk3];
				vec_t[ZZ] += Ematr_z[kk1][kk2][kk3];
//                                rvec_inc(Kern->vec_interp_quant_grid[i],Efield_temp);
			}
		}	
	 }
         */
//         mvmul(invcosdirmat,vec_t,Kern->vec_interp_quant_grid[i]);
	// Instead of rotating the electric field into the "molecular frame", here we leave it in the lab frame;
	// it's assumed that we're going to have only one grid point at the position of the oxygen atom, and thus
	// will calculate the true (SPME) electric field.
	Kern->vec_interp_quant_grid[i][XX] = vec_t[XX];
	Kern->vec_interp_quant_grid[i][YY] = vec_t[YY];
	Kern->vec_interp_quant_grid[i][ZZ] = vec_t[ZZ];
        printf("mol_grid_espl %f %f\n",Kern->translgrid[i][ZZ],Kern->vec_interp_quant_grid[i][ZZ]);
        printf("four_grid_espl %f %f\n",grid_spacing[ZZ]*bin_indz0,Kern->quantity_on_grid_z[bin_indx0][bin_indy0][bin_indz0]-Emean[ZZ]);
     }

     for ( i = 0; i < ir->nkz-1; i++)
     {
         printf("grid_efield_z %f %f\n",grid_spacing[ZZ]*i,Kern->quantity_on_grid_z[0][0][i]);
     }


     for (a = 0; a < ir->nkx; a++)
     {
         for (b = 0; b < ir->nky; b++)
         {
             sfree(Ematr_x[a][b]);
             sfree(Ematr_y[a][b]);
             sfree(Ematr_z[a][b]);
         }
         sfree(Ematr_x[a]);
         sfree(Ematr_y[a]);
         sfree(Ematr_z[a]);
     }
     sfree(Ematr_x);
     sfree(Ematr_y);
     sfree(Ematr_z);

//     fprintf(stderr,"freed everything\n");  
}

void vec_lagrange_interpolation_kern(t_Kern *Kern, t_inputrec *ir, t_pbc *pbc, matrix invcosdirmat, rvec xi, rvec grid_invspacing, rvec grid_spacing, rvec Emean)
{
     int i, j, k, l, ix, iy, iz, d, ind0;
     int bin_indx0, bin_indy0, bin_indz0, bin_indx1, bin_indy1, bin_indz1;
     real xd, yd, zd, sel_dist;
     rvec delx, delxdeb, vec_t;

     sfree(Kern->vec_interp_quant_grid);
     snew(Kern->vec_interp_quant_grid,Kern->gridpoints);

     for (i = 0; i < Kern->gridpoints; i++)
     {


         bin_indx1 = ceil((Kern->translgrid[i][XX] )*grid_invspacing[XX] );
         bin_indy1 = ceil((Kern->translgrid[i][YY] )*grid_invspacing[YY] );
         bin_indz1 = ceil((Kern->translgrid[i][ZZ] )*grid_invspacing[ZZ] );

         if (bin_indx1 > ir->nkx-1)
         {
            bin_indx1 = 0;
         }
         if (bin_indz1 > ir->nkz-1)
         {
            bin_indz1 = 0;
         }
         if (bin_indy1 > ir->nky-1)
         {
            bin_indy1 = 0;
         }

         bin_indx0 = (bin_indx1 > 0 ) ? bin_indx1  -1 : ir->nkx-1 ;
         bin_indy0 = (bin_indy1 > 0 ) ? bin_indy1  -1 : ir->nky-1 ;
         bin_indz0 = (bin_indz1 > 0 ) ? bin_indz1  -1 : ir->nkz-1 ;

         Kern->rspace_grid[XX]=grid_spacing[XX]*bin_indx0;
         Kern->rspace_grid[YY]=grid_spacing[YY]*bin_indy0;
         Kern->rspace_grid[ZZ]=grid_spacing[ZZ]*bin_indz0;

         pbc_dx(pbc,Kern->translgrid[i],Kern->rspace_grid,delx);

         int npoints = 1;

					// Get a list of the points that we will be using for interpolation.
				  int *xlist,*ylist,*zlist;
					snew (xlist,2*npoints+1);
					snew (ylist,2*npoints+1);
					snew (zlist,2*npoints+1);
					for (j=1;j<=2*npoints;j++)
					{
						xlist[j] = index_wrap(bin_indx0 - npoints + j - 1,ir->nkx);
						ylist[j] = index_wrap(bin_indy0 - npoints + j - 1,ir->nky);
						zlist[j] = index_wrap(bin_indz0 - npoints + j - 1,ir->nkz);
					}

					// Get arrays of the x1, x2, x3 values, as floats. We will take the spatial grid to start from 1.0, and go to the decimal version of
					// 2*npoints. This will avoid any potential unpleasantness with periodic boundary conditions (wherein the actual coordinates of the
					// grid points might have a discontinuity).
					float *x1a,*x2a,*x3a;
					snew(x1a,2*npoints+1);
					snew(x2a,2*npoints+1);
					snew(x3a,2*npoints+1);
					for (j=1;j<=2*npoints;j++)
					{
						x1a[j] = (float)j;
						x2a[j] = (float)j;
						x3a[j] = (float)j;
					}

					// Now get an array of the electric field components on the grid points.
					float ***efield_x, ***efield_y, ***efield_z;
					snew(efield_x,2*npoints+1);snew(efield_y,2*npoints+1);snew(efield_z,2*npoints+1);
					for (j=0;j<2*npoints+1;j++){snew(efield_x[j],2*npoints+1);snew(efield_y[j],2*npoints+1);snew(efield_z[j],2*npoints+1);}
					for (j=0;j<2*npoints+1;j++){for (k=0;k<2*npoints+1;k++){snew(efield_x[j][k],2*npoints+1);snew(efield_y[j][k],2*npoints+1);snew(efield_z[j][k],2*npoints+1);}}
					for (j=1;j<=2*npoints;j++)
					{
						for (k=1;k<=2*npoints;k++)
						{
							for (l=1;l<=2*npoints;l++)
							{
								efield_x[j][k][l] = Kern->quantity_on_grid_x[xlist[j]][ylist[k]][zlist[l]];
								efield_y[j][k][l] = Kern->quantity_on_grid_y[xlist[j]][ylist[k]][zlist[l]];
								efield_z[j][k][l] = Kern->quantity_on_grid_z[xlist[j]][ylist[k]][zlist[l]];
							}
						}
					}

					// Finally, in terms of the coordinates chosen, what is the point on which we wish to interpolate the electric field?
					

//Kern->quantity_on_grid_x[bin_indx0][bin_indy0][bin_indz0]

/*
	 float *x1a,*x2a,*x3a;
	 snew(x1a,2*npoints +1);
	 snew(x2a,2*npoints +1);
	 snew(x3a,2*npoints +1);

	 j = 1;
	 for (i=bin_indx0-npoints-1;i<bin_indx0+npoints;i++)
	 {
	 	x1a[j] = Kern->grid_spacing[XX]*index_wrap(i,ir->nkx);
		j++;
	 }
	 j = 1;
	 for (i=bin_indy0-npoints-1;i<bin_indy0+npoints;i++)
	 {
		x2a[j] = Kern->grid_spacing[YY]*index_wrap(i,ir->nky);
		j++;
	 }
         j = 1;
         for (i=bin_indz0-npoints-1;i<bin_indz0+npoints;i++)
         {
                x3a[j] = Kern->grid_spacing[ZZ]*index_wrap(i,ir->nkz);
                j++;
         }
*/

         xd = fabs(delx[XX])*grid_invspacing[XX];
         yd = fabs(delx[YY])*grid_invspacing[YY];
         zd = fabs(delx[ZZ])*grid_invspacing[ZZ];



         copy_rvec(vec_t,Kern->vec_interp_quant_grid[i]);
	}
}

void vec_trilinear_interpolation_kern(t_Kern *Kern, t_inputrec *ir, t_pbc *pbc, matrix invcosdirmat, rvec xi, rvec grid_invspacing, rvec grid_spacing, rvec Emean)
{
     int i, ix, iy, iz, d, ind0;
     int bin_indx0, bin_indy0, bin_indz0, bin_indx1, bin_indy1, bin_indz1;
     real xd, yd, zd, sel_dist;
     rvec delx, delxdeb, vec_t;

     sfree(Kern->vec_interp_quant_grid);
     snew(Kern->vec_interp_quant_grid,Kern->gridpoints);

     for (i = 0; i < Kern->gridpoints; i++)
     {

/*         if (debug)
         {
           Kern->translgrid[i][XX]= 0.0;
           Kern->translgrid[i][YY]= 0.0;
           Kern->translgrid[i][ZZ]= i*0.11;
           Emean[ZZ]=0.0;
         }
*/

         bin_indx1 = ceil((Kern->translgrid[i][XX] )*grid_invspacing[XX] );
         bin_indy1 = ceil((Kern->translgrid[i][YY] )*grid_invspacing[YY] );
         bin_indz1 = ceil((Kern->translgrid[i][ZZ] )*grid_invspacing[ZZ] );

         if (bin_indx1 > ir->nkx-1)
         {
            bin_indx1 = 0;
         }
         if (bin_indz1 > ir->nkz-1)
         {
            bin_indz1 = 0;
         }
         if (bin_indy1 > ir->nky-1)
         {
            bin_indy1 = 0;
         }

         bin_indx0 = (bin_indx1 > 0 ) ? bin_indx1  -1 : ir->nkx-1 ;
         bin_indy0 = (bin_indy1 > 0 ) ? bin_indy1  -1 : ir->nky-1 ;
         bin_indz0 = (bin_indz1 > 0 ) ? bin_indz1  -1 : ir->nkz-1 ;

         Kern->rspace_grid[XX]=grid_spacing[XX]*bin_indx0;
         Kern->rspace_grid[YY]=grid_spacing[YY]*bin_indy0;
         Kern->rspace_grid[ZZ]=grid_spacing[ZZ]*bin_indz0;

         pbc_dx(pbc,Kern->translgrid[i],Kern->rspace_grid,delx);

         xd = fabs(delx[XX])*grid_invspacing[XX];
         yd = fabs(delx[YY])*grid_invspacing[YY];
         zd = fabs(delx[ZZ])*grid_invspacing[ZZ];

//         pbc_dx(pbc,xi,Kern->translgrid[i],delx);

                                vec_t[XX] = (Kern->quantity_on_grid_x[bin_indx0][bin_indy0][bin_indz0]-Emean[XX])*(1.0 - xd)*(1.0 - zd)*(1 - yd) +
                                         (Kern->quantity_on_grid_x[bin_indx1][bin_indy0][bin_indz0]-Emean[XX])*xd*(1 - yd)*(1 - zd) +
                                         (Kern->quantity_on_grid_x[bin_indx0][bin_indy1][bin_indz0]-Emean[XX])*(1.0 - xd)*yd*(1 - zd) +
                                         (Kern->quantity_on_grid_x[bin_indx0][bin_indy0][bin_indz1]-Emean[XX])*(1.0 - xd)*(1 - yd)*zd +
                                         (Kern->quantity_on_grid_x[bin_indx1][bin_indy0][bin_indz1]-Emean[XX])*xd*(1 - yd)*zd +
                                         (Kern->quantity_on_grid_x[bin_indx0][bin_indy1][bin_indz1]-Emean[XX])*(1 - xd)*yd*zd +
                                         (Kern->quantity_on_grid_x[bin_indx1][bin_indy1][bin_indz0]-Emean[XX])*xd*yd*(1 - zd) +
                                         (Kern->quantity_on_grid_x[bin_indx1][bin_indy1][bin_indz1]-Emean[XX])*xd*yd*zd;

                                vec_t[YY] = (Kern->quantity_on_grid_y[bin_indx0][bin_indy0][bin_indz0]-Emean[YY])*(1.0 - xd)*(1.0 - zd)*(1 - yd) +
                                         (Kern->quantity_on_grid_y[bin_indx1][bin_indy0][bin_indz0]-Emean[YY])*xd*(1 - yd)*(1 - zd) +
                                         (Kern->quantity_on_grid_y[bin_indx0][bin_indy1][bin_indz0]-Emean[YY])*(1.0 - xd)*yd*(1 - zd) +
                                         (Kern->quantity_on_grid_y[bin_indx0][bin_indy0][bin_indz1]-Emean[YY])*(1.0 - xd)*(1 - yd)*zd +
                                         (Kern->quantity_on_grid_y[bin_indx1][bin_indy0][bin_indz1]-Emean[YY])*xd*(1 - yd)*zd +
                                         (Kern->quantity_on_grid_y[bin_indx0][bin_indy1][bin_indz1]-Emean[YY])*(1 - xd)*yd*zd +
                                         (Kern->quantity_on_grid_y[bin_indx1][bin_indy1][bin_indz0]-Emean[YY])*xd*yd*(1 - zd) +
                                         (Kern->quantity_on_grid_y[bin_indx1][bin_indy1][bin_indz1]-Emean[YY])*xd*yd*zd;

                                vec_t[ZZ] = (Kern->quantity_on_grid_z[bin_indx0][bin_indy0][bin_indz0]-Emean[ZZ])*(1.0 - xd)*(1.0 - zd)*(1 - yd) +
                                         (Kern->quantity_on_grid_z[bin_indx1][bin_indy0][bin_indz0]-Emean[ZZ])*xd*(1 - yd)*(1 - zd) +
                                         (Kern->quantity_on_grid_z[bin_indx0][bin_indy1][bin_indz0]-Emean[ZZ])*(1.0 - xd)*yd*(1 - zd) +
                                         (Kern->quantity_on_grid_z[bin_indx0][bin_indy0][bin_indz1]-Emean[ZZ])*(1.0 - xd)*(1 - yd)*zd +
                                         (Kern->quantity_on_grid_z[bin_indx1][bin_indy0][bin_indz1]-Emean[ZZ])*xd*(1 - yd)*zd +
                                         (Kern->quantity_on_grid_z[bin_indx0][bin_indy1][bin_indz1]-Emean[ZZ])*(1 - xd)*yd*zd +
                                         (Kern->quantity_on_grid_z[bin_indx1][bin_indy1][bin_indz0]-Emean[ZZ])*xd*yd*(1 - zd) +
                                         (Kern->quantity_on_grid_z[bin_indx1][bin_indy1][bin_indz1]-Emean[ZZ])*xd*yd*zd;

//                                 mvmul(invcosdirmat,vec_t,Kern->vec_interp_quant_grid[i]);
              copy_rvec(vec_t,Kern->vec_interp_quant_grid[i]);
//            Kern->vec_interp_quant_grid[i][d] *= Kern->weights[i];
//            Kern->vec_interp_quant_grid[i][d] -= Kern->selfterm[i];

         printf("xi %f %f %f\n",xi[XX],xi[YY],xi[ZZ]);
         printf("grid %f %f %f\n",Kern->translgrid[i][XX],Kern->translgrid[i][YY],Kern->translgrid[i][ZZ]);
         printf("grid index %d %d %d\n",bin_indx0,bin_indy0,bin_indz0);
         printf("nearest grid 0 0 0 %f %f %f %f\n",Kern->rspace_grid[XX],Kern->rspace_grid[YY],Kern->rspace_grid[ZZ],Kern->quantity_on_grid_z[bin_indx0][bin_indy0][bin_indz0]-Emean[ZZ]);
         printf("nearest grid 1 1 1 %f %f %f %f\n",grid_spacing[XX]*bin_indx1,grid_spacing[XX]*bin_indy1,grid_spacing[XX]*bin_indz1,Kern->quantity_on_grid_z[bin_indx1][bin_indy1][bin_indz1] -Emean[ZZ]);
         printf("molecular_efield_z %f %f\n",Kern->translgrid[i][ZZ],Kern->vec_interp_quant_grid[i][ZZ]);
         printf("distance grid to x %f \n",norm(delx));
         printf("d= %f electric_field_bare= %f %f %f\n",norm(delx),vec_t[XX], vec_t[YY], vec_t[ZZ]);
         printf("d= %f electric field rotated= %f %f %f\n",norm(delx),Kern->vec_interp_quant_grid[i][XX], Kern->vec_interp_quant_grid[i][YY], Kern->vec_interp_quant_grid[i][ZZ]);
         pbc_dx(pbc,Kern->rspace_grid,xi,delxdeb);
         printf("xd %f yd %f zd %f\n",xd,yd,zd);

     }
     for ( i = 0; i < ir->nkz-1; i++)
     {
         printf("grid_efield_z %f %f\n",grid_spacing[ZZ]*i,Kern->quantity_on_grid_z[0][0][i]);
     }
     gmx_fatal(FARGS,"exit from elfield\n");
}

//
// functions to compute the long range electrostatic potential from David Wilkins
//

void setup_ewald_pair_potential(int *grid, int interp_order, int kmax,t_complex ****FT_pair_pot,real invkappa2)
{
	real ***pair_potential,***Bmatr,***Mmatr;
	real **bx,**by,**bz, pi, expfac, scalfac, tpi;
	t_complex ***pK, ***mK;
	real fx, fy, fz, ***cx, ***cy, ***cz, *cx_temp, *cy_temp, *cz_temp,c, s, f, arg, arg_mult;
	int i,j,k,m, kx, ky, kz, k2;
        

	pi = M_PI ;

	// Initialize pair potential matrix, M and B matrices.
	snew(pair_potential,grid[0]);
	snew(Mmatr,grid[0]);
        snew(pK,grid[0]);
        snew(mK,grid[0]);
        snew(*FT_pair_pot,grid[0]);
	for (i = 0; i < grid[0]; i++)
	{
		snew(pair_potential[i],grid[1]);
		snew(Mmatr[i],grid[1]);
                snew((*FT_pair_pot)[i],grid[1]);
                snew(pK[i],grid[1]);
                snew(mK[i],grid[1]);
		for (j = 0; j < grid[1]; j++)
		{
			snew(pair_potential[i][j],grid[2]);
			snew(Mmatr[i][j],grid[2]);
                        snew((*FT_pair_pot)[i][j],grid[2]);
                        snew(pK[i][j],grid[2]/2 + 1);
                        snew(mK[i][j],grid[2]/2 + 1);
		}
	}
        fprintf(stderr,"allocated potential and matrices\n");
	snew(Bmatr,kmax+1);
	for (i = 0;i <= kmax; i++)
	{
		snew(Bmatr[i],kmax+1);
		for (j = 0;j <= kmax; j++)
		{
			snew(Bmatr[i][j],kmax+1);
		}
	}
        fprintf(stderr,"allocated Bmatrix\n");
	// Calculate the M matrix. This does not show up in the original SPME paper, but is convenient if we want to calculate
	// the electrostatic potential.
	for (i = 0; i < grid[0]; i++)
	{
		for (j = 0; j < grid[1]; j++)
		{
			for (k = 0;k < grid[2]; k++)
			{
				Mmatr[i][j][k] = Bspline(i,interp_order)*Bspline(j,interp_order)*Bspline(k,interp_order);
			}
		}
	}

	// Calculate b coefficients for the B matrix.
	snew(bx,2);
	snew(by,2);
	snew(bz,2);
	snew(bx[0],kmax+1);
	snew(bx[1],kmax+1);
	snew(by[0],kmax+1);
	snew(by[1],kmax+1);
	snew(bz[0],kmax+1);
	snew(bz[1],kmax+1);
	for (i = 0; i <= kmax; i++)
	{
		for (m = 0; m<interp_order-1;m++)
		{
			bx[0][i] += Bspline(m+1,interp_order)*cos(2.0*pi*m*i/grid[0]);
			bx[1][i] += Bspline(m+1,interp_order)*sin(2.0*pi*m*i/grid[0]);
		}
		bx[0][i] = 1.0/(bx[0][i]*bx[0][i] + bx[1][i]*bx[1][i]);
	}

	for (j = 0; j <= kmax; j++)
	{
		for (m = 0; m<interp_order-1;m++)
		{
			by[0][j] += Bspline(m+1,interp_order)*cos(2.0*pi*m*j/grid[1]);
			by[1][j] += Bspline(m+1,interp_order)*sin(2.0*pi*m*j/grid[1]);
		}
		by[0][j] = 1.0/(by[0][j]*by[0][j] + by[1][j]*by[1][j]);
	}
	for (k = 0; k <= kmax; k++)
	{
		for (m = 0; m<interp_order-1;m++)
		{
			bz[0][k] += Bspline(m+1,interp_order)*cos(2.0*pi*m*k/grid[2]);
			bz[1][k] += Bspline(m+1,interp_order)*sin(2.0*pi*m*k/grid[2]);
		}
		bz[0][k] = 1.0/(bz[0][k]*bz[0][k] + bz[1][k]*bz[1][k]);
	}
	// Calculate the B matrix (this is required for calculating the pair potential).
	for (i = 0; i <= kmax; i++)
	{
		for (j = 0; j <= kmax; j++)
		{
			for (k = 0; k <= kmax; k++)
			{
				Bmatr[i][j][k] = bx[0][i]*by[0][j]*bz[0][k];
			}
		}
	}
	// Calculate the pair potential.
	expfac = -1.0 * M_PI * M_PI * invkappa2;
	tpi = 2.0 * M_PI;
	scalfac = 1.0 / ( grid[0]*grid[1]*grid[2]*M_PI);  

	snew(bx,2);
	snew(by,2);
	snew(bz,2);
	snew(bx[0],grid[0]);
	snew(bx[1],grid[0]);
	snew(by[0],grid[1]);
	snew(by[1],grid[1]);
	snew(bz[0],grid[2]);
	snew(bz[1],grid[2]);

	snew(cx,2);
	snew(cy,2);
	snew(cz,2);
	snew(cx[0],grid[0]);
	snew(cx[1],grid[0]);
	snew(cy[0],grid[1]);
	snew(cy[1],grid[1]);
	snew(cz[0],grid[2]);
	snew(cz[1],grid[2]);
        snew(cx_temp,grid[0]);
        snew(cy_temp,grid[1]);
        snew(cz_temp,grid[2]);
	for (i=0;i<grid[0];i++)
	{
		snew(cx[0][i],kmax+1);
		snew(cx[1][i],kmax+1);
	}
	for (j=0;j<grid[1];j++)
	{
		snew(cy[0][j],kmax+1);
		snew(cy[1][j],kmax+1);
	}
	for (k=0;k<grid[2];k++)
	{
		snew(cz[0][k],kmax+1);
		snew(cz[1][k],kmax+1);
	}

	for (i = 0;i < grid[0]; i++)
	{
		cx[0][i][0] = 1.0;
		cx[1][i][0] = 0.0;
		fx = (real)(i) / grid[0];
		c = cos(tpi * fx);
		s = sin(tpi * fx);
		for (m = 1;m <= kmax;m++)
		{
			cx[0][i][m] = c*cx[0][i][m-1] - s*cx[1][i][m-1];
			cx[1][i][m] = s*cx[0][i][m-1] + c*cx[1][i][m-1];
		}
	}

	for (j = 0;j < grid[1]; j++)
	{
		cy[0][j][0] = 1.0;
		cy[1][j][0] = 0.0;
		fy = (real)(j) / grid[1];
		c = cos(tpi * fy);
		s = sin(tpi * fy);
		for (m = 1;m <= kmax;m++)
		{
			cy[0][j][m] = c*cy[0][j][m-1] - s*cy[1][j][m-1];
			cy[1][j][m] = s*cy[0][j][m-1] + c*cy[1][j][m-1];
		}
	}
	for (k = 0;k < grid[2]; k++)
	{
		cz[0][k][0] = 1.0;
		cz[1][k][0] = 0.0;
		fz = (real)(k) / grid[2];
		c = cos(tpi * fz);
		s = sin(tpi * fz);
		for (m = 1;m <= kmax;m++)
		{
			cz[0][k][m] = c*cz[0][k][m-1] - s*cz[1][k][m-1];
			cz[1][k][m] = s*cz[0][k][m-1] + c*cz[1][k][m-1];
		}
	}

	f = 1.0;
	for (kx = 0;kx <= kmax; kx++)
	{
		if (kx==1){f = 2.0;}
		ky = 0;
		kz = 0;
		k2 = kx*kx;

                for (i = 0; i < grid[0]; i ++)
                {
                    cx_temp[i]=cx[0][i][kx];
                }

		// Contributions to the pair potential from ky = 0 and kz = 0
		if (k2>0)
		{
			arg = exp(k2 * expfac)*f*Bmatr[kx][ky][kz]/k2;
			for (i = 0;i<grid[0];i++)
			{
                                arg_mult=arg*cx_temp[i];
				for (j = 0;j<grid[1];j++)
				{
					for (k = 0;k<grid[2];k++)
					{
						pair_potential[i][j][k] += arg_mult;
					}
				}
			}
		}

		// Contributions from ky = 0 and |kz| > 0
		for (kz = 1;kz <= kmax;kz++)
		{
			ky = 0;
			k2 = kx*kx + kz*kz;
			arg = exp(k2 * expfac)*2.0*f*Bmatr[kx][ky][kz]/k2;
                        for (k = 0; k < grid[2]; k++)
                        {
                            cz_temp[k]=cz[0][k][kz];
                        }
			for (i = 0;i<grid[0];i++)
			{
                                arg_mult=arg*cx_temp[i];
				for (j = 0;j<grid[1];j++)
				{
					for (k = 0;k<grid[2];k++)
					{
						pair_potential[i][j][k] += arg_mult*cz_temp[k];
					}
				}
			}
		}

		for (ky = 1;ky <= kmax;ky++)
		{
			// Contributions from |ky| > 0 and kz = 0
			kz = 0;
			k2 = kx*kx + ky*ky;
			arg = exp(k2*expfac) *2.0*f*Bmatr[kx][ky][kz]/ k2;
                        for (j = 0; j < grid[1]; j++)
                        {
                            cy_temp[j]=cy[0][j][ky];
                        }
			for (i = 0;i<grid[0];i++)
			{
				for (j = 0;j<grid[1];j++)
				{
                                        arg_mult=arg*cx_temp[i]*cy_temp[j];
					for (k = 0;k<grid[2];k++)
					{
						pair_potential[i][j][k] += arg_mult;
					}
				}
			}

			// Contributions from |ky| > 0 and |kz| > 0
			for (kz = 1;kz <= kmax;kz++)
			{
				k2 = kx*kx + ky*ky + kz*kz;
				arg = exp(k2*expfac)*4.0*f*Bmatr[kx][ky][kz] / k2;
                                for (k = 0; k <grid[2]; k++)
                                {
                                    cz_temp[k]=cz[0][k][kz];
                                }
				for (i = 0;i<grid[0];i++)
				{
					for (j = 0;j<grid[1];j++)
					{       arg_mult=arg*cx_temp[i]*cy_temp[j];
						for (k = 0;k<grid[2];k++)
						{
							pair_potential[i][j][k] += arg_mult*cz_temp[k];
						}
					}
				}
			}
		}
	}

	// Scale the pair potential here to avoid having to do so later on.
	for (i=0;i<grid[0];i++)
	{
		for (j=0;j<grid[1];j++)
		{
			for (k=0;k<grid[2];k++)
			{
				pair_potential[i][j][k] *= scalfac;
			}
		}
	}


	// Now Fourier transform the pair potential and the M matrix.

	do_fft(pair_potential,pK,grid,1.0,GMX_FFT_REAL_TO_COMPLEX);
	do_fft(Mmatr,mK,grid,1.0,GMX_FFT_REAL_TO_COMPLEX);

//	// Multiply the FTed pair potential and element matrix element-wise to get the theta matrix.
	for (i = 0;i < grid[0];i++)
	{
		for (j = 0;j < grid[1];j++)
		{
			for (k=0;k < grid[2]/2 + 1;k++)
			{
				(*FT_pair_pot)[i][j][k].re = (pK[i][j][k].re*mK[i][j][k].re - pK[i][j][k].im*mK[i][j][k].im);
				(*FT_pair_pot)[i][j][k].im = (pK[i][j][k].re*mK[i][j][k].im + pK[i][j][k].im*mK[i][j][k].re);
			}
		}
	}

        for (i = 0; i < grid[0]; i++)
        {
                for (j = 0; j < grid[1]; j++)
                {
                        sfree(pair_potential[i][j]);
                        sfree(Mmatr[i][j]);
                        sfree(pK[i][j]);
                        sfree(mK[i][j]);
                }
                sfree(pair_potential[i]);
                sfree(Mmatr[i]);
                sfree(pK[i]);
                sfree(mK[i]);
        }
        sfree(pair_potential);
        sfree(Mmatr);        
        sfree(pK);
        sfree(mK);

        sfree(bx[0]);
        sfree(bx[1]);
        sfree(by[0]);
        sfree(by[1]);
        sfree(bz[0]);
        sfree(bz[1]);
        sfree(bx); sfree(by); sfree(bz);

        for (i=0;i<grid[0];i++)
        {
                sfree(cx[0][i]);
                sfree(cx[1][i]);
        }
        for (j=0;j<grid[1];j++)
        {
                sfree(cy[0][j]);
                sfree(cy[1][j]);
        }
        for (k=0;k<grid[2];k++)
        {
                sfree(cz[0][k]);
                sfree(cz[1][k]);
        }
        sfree(cx[0]);
        sfree(cx[1]);
        sfree(cy[0]); 
        sfree(cy[1]); 
        sfree(cz[0]);
        sfree(cz[1]);
        sfree(cx); sfree(cy); sfree(cz);
        sfree(cx_temp); sfree(cy_temp); sfree(cz_temp);

}

void calculate_spme_efield(t_Kern *Kern, t_inputrec *ir, t_topology *top,
                         matrix box, real invvol, t_block *mols, int  *molindex[],
                         int *chged_atom_indexes, int n_chged_atoms, int *grid, rvec grid_spacing,
                         int interp_order, rvec *x, int isize0,
                         t_complex ***FT_pair_pot, rvec *Emean, real eps)
{

        // Here we use the calculated FT_pair potential to find the electrostatic potential on the grid points.
        // Firstly, we must create the charge matrix [as in JCP 103 (19), 1995]. We Fourier transform this and
        // then multiply it by FT_pair_potential. Inverse FTing the result will give the potential on these points.
        // We require that the input coordinates be fractional coordinates (though this should be easy to alter).
        real charge, ***Qmatr_z, ***Qmatr_y, ***Qmatr_x;
        real mult_fac, scalfac;
        t_complex ****qF, ****convF, ***convF_z, ***qF_z, ***convF_y, ***qF_y, ***convF_x, ***qF_x;
        int n, i, j, k, m, d,points[3], kk1, kk2, kk3, idx;
        int ix, iy;
        int ind0;
        rvec xi, invbox, du;
        rvec ***Qmatr,bsplcoeff, d_bsplcoeff;
	rvec dielectric;

        invbox[XX] =1.0/ box[XX][XX]; invbox[YY] = 1.0/box[YY][YY]; invbox[ZZ] = 1.0/box[ZZ][ZZ];
        mult_fac = pow(invvol/(grid_spacing[XX]*grid_spacing[YY]*grid_spacing[ZZ]),1.0/3.0);
        for (ix = 0; ix < ir->nkx; ix++)
        {
            for (iy = 0; iy < ir->nky ; iy ++)
            {
              sfree(Kern->quantity_on_grid_x[ix][iy]);
              sfree(Kern->quantity_on_grid_y[ix][iy]);
              sfree(Kern->quantity_on_grid_z[ix][iy]);
            }
            sfree(Kern->quantity_on_grid_x[ix]);
            sfree(Kern->quantity_on_grid_y[ix]);
            sfree(Kern->quantity_on_grid_z[ix]);
        }
        sfree(Kern->quantity_on_grid_x);
        sfree(Kern->quantity_on_grid_y);
        sfree(Kern->quantity_on_grid_z);

        snew(Kern->quantity_on_grid_x,ir->nkx);
        snew(Kern->quantity_on_grid_y,ir->nkx);
        snew(Kern->quantity_on_grid_z,ir->nkx);
        for (ix = 0; ix < ir->nkx; ix++)
        {
            snew(Kern->quantity_on_grid_x[ix], ir->nky);
            snew(Kern->quantity_on_grid_y[ix], ir->nky);
            snew(Kern->quantity_on_grid_z[ix], ir->nky);
            for (iy = 0; iy < ir->nky ; iy ++)
            {
              snew(Kern->quantity_on_grid_x[ix][iy], ir->nkz);
              snew(Kern->quantity_on_grid_y[ix][iy], ir->nkz);
              snew(Kern->quantity_on_grid_z[ix][iy], ir->nkz);
            }
        }

        //fprintf(stderr,"allocated potential on grid\n ");

        snew(Qmatr_x,grid[0]);
        snew(qF_x,grid[0]);
        snew(convF_x,grid[0]);
        snew(Qmatr_y,grid[0]);
        snew(qF_y,grid[0]);
        snew(convF_y,grid[0]);
        snew(Qmatr_z,grid[0]);
        snew(qF_z,grid[0]);
        snew(convF_z,grid[0]);
        for (i = 0;i<grid[0];i++)
        {
                snew(Qmatr_x[i],grid[1]);
                snew(qF_x[i],grid[1]);
                snew(convF_x[i],grid[1]);
                snew(Qmatr_y[i],grid[1]);
                snew(qF_y[i],grid[1]);
                snew(convF_y[i],grid[1]);
                snew(Qmatr_z[i],grid[1]);
                snew(qF_z[i],grid[1]);
                snew(convF_z[i],grid[1]);
                for (j = 0;j<grid[1];j++)
                {
                        snew(Qmatr_x[i][j],grid[2]);
                        snew(qF_x[i][j],grid[2]/2 + 1);
                        snew(convF_x[i][j],grid[2]/2 + 1);
                        snew(Qmatr_y[i][j],grid[2]);
                        snew(qF_y[i][j],grid[2]/2 + 1);
                        snew(convF_y[i][j],grid[2]/2 + 1);
                        snew(Qmatr_z[i][j],grid[2]);
                        snew(qF_z[i][j],grid[2]/2 + 1);
                        snew(convF_z[i][j],grid[2]/2 + 1);
                }
        }

        //fprintf(stderr,"allocated charge matrix on grid\n ");

        // NOTE: as currently written, we also require that periodic boundary conditions have been applied to
        // the coordinates.
	for (i=0;i<3;i++)
	{
		// Zero the electric-field correction due to dielectric boundary conditions.
		dielectric[i] = 0.0;
	}
        for (n = 0;n < isize0 ;n++)
        {
                for (m = 0; m <  n_chged_atoms  ; m++)
                {
                       ind0 = mols->index[molindex[0][n]]   + chged_atom_indexes[m] ;
                       copy_rvec(x[ind0],xi);
                        charge = top->atoms.atom[ind0].q; 
//                       charge= 1.0;
                       //fprintf(stderr,"charge =%f\n",charge);
                       //fprintf(stderr,"chged atoms=%d chged_atom_indexes %d atom ind %d\n",n_chged_atoms,chged_atom_indexes[m],ind0);
                       //Find the nearest grid points to charge n.
                       for (i = 0;i < 3;i++)
                       {
			       if (eps > 0.0)
			       {
				       // If we have as input a negative dielectric constant, this means that we don't want to include these
				       // boundary conditions.
				       dielectric[i] -= charge*xi[i];
			       }
                               xi[i] *= invbox[i];                            
                               points[i] = floor(xi[i]*grid[i]);
                       }
                       // Now, loop over all the grid points onto which this charge is allocated. In each dimension,
                       // there are 2*interp_order.
                       for (i = points[0]-interp_order+1;i<=points[0]+interp_order;i++)
                       {
                               du[XX] = xi[XX]*grid[XX] - i;
                               kk1 = index_wrap(i,grid[XX]);
                               bsplcoeff[XX] = Bspline(du[XX],interp_order);
                               d_bsplcoeff[XX] = (Bspline(du[XX],interp_order-1) - Bspline(du[XX]-1.0,interp_order-1));
                               for (j = points[1]-interp_order+1;j<=points[1]+interp_order;j++)
                               {
                                       du[YY] = xi[YY]*grid[YY] - j;
                                       kk2 = index_wrap(j,grid[YY]);
                                       bsplcoeff[YY] =Bspline(du[YY],interp_order);
                                       d_bsplcoeff[YY] = (Bspline(du[YY],interp_order-1) - Bspline(du[YY]-1.0,interp_order-1));
                                       for (k = points[2]-interp_order+1;k<=points[2]+interp_order;k++)
                                       {
                                               du[ZZ] = xi[ZZ]*grid[ZZ] - k;
                                               kk3 = index_wrap(k,grid[ZZ]);
                                               bsplcoeff[ZZ] = Bspline(du[ZZ],interp_order);
                                               //Qmatr_x[kk1][kk2][kk3] += charge*bsplcoeff[XX]*bsplcoeff[YY]* bsplcoeff[ZZ];
                                               d_bsplcoeff[ZZ] = (Bspline(du[ZZ],interp_order-1)-Bspline(du[ZZ]-1.0,interp_order-1));
                                               //Qmatr_x[kk1][kk2][kk3] += charge*bsplcoeff[XX]*bsplcoeff[YY]* bsplcoeff[ZZ];
                                               Qmatr_x[kk1][kk2][kk3] -= charge*d_bsplcoeff[XX]*bsplcoeff[YY]*bsplcoeff[ZZ];
                                               Qmatr_y[kk1][kk2][kk3] -= charge*bsplcoeff[XX]*d_bsplcoeff[YY]*bsplcoeff[ZZ];
                                               Qmatr_z[kk1][kk2][kk3] -= charge*bsplcoeff[XX]*bsplcoeff[YY]*d_bsplcoeff[ZZ];
                                       }
                               }
                       }
                }
        }

        //fprintf(stderr,"computed Qmatr\n");
        // Now Fourier transform the charge matrix.
        do_fft(Qmatr_z,qF_z,grid,1.0,GMX_FFT_REAL_TO_COMPLEX);
        do_fft(Qmatr_y,qF_y,grid,1.0,GMX_FFT_REAL_TO_COMPLEX);
        do_fft(Qmatr_x,qF_x,grid,1.0,GMX_FFT_REAL_TO_COMPLEX);

        // Multiply the FTed charge matrix by the FTed Ewald pair potential.
        // NOTE: the "pair potential" here is not exacly the same as the JCP paper cited
        // at the beginning of this function. It is equal to the potential of this paper, convoluted
        // with the matrix of Bspline functions. Using this modified pair potential, we can obtain very
        // straightforwardly the electrostatic potential on the grid.
        for (i = 0;i < grid[0];i++)
        {
                for (j = 0;j < grid[1];j++)
                {
                        for (k = 0;k < grid[2]/2 + 1;k++)
                        {
                                convF_z[i][j][k].re = qF_z[i][j][k].re*FT_pair_pot[i][j][k].re - qF_z[i][j][k].im*FT_pair_pot[i][j][k].im;
                                convF_z[i][j][k].im = qF_z[i][j][k].re*FT_pair_pot[i][j][k].im + qF_z[i][j][k].im*FT_pair_pot[i][j][k].re;
                                convF_x[i][j][k].re = qF_x[i][j][k].re*FT_pair_pot[i][j][k].re - qF_x[i][j][k].im*FT_pair_pot[i][j][k].im;
                                convF_x[i][j][k].im = qF_x[i][j][k].re*FT_pair_pot[i][j][k].im + qF_x[i][j][k].im*FT_pair_pot[i][j][k].re;
                                convF_y[i][j][k].re = qF_y[i][j][k].re*FT_pair_pot[i][j][k].re - qF_y[i][j][k].im*FT_pair_pot[i][j][k].im;
                                convF_y[i][j][k].im = qF_y[i][j][k].re*FT_pair_pot[i][j][k].im + qF_y[i][j][k].im*FT_pair_pot[i][j][k].re;
                        }
                }
        }
        // FT it back to get the convolution in real space. This convolution is the electrostatic potential on the grid points.
      
        do_fft(Kern->quantity_on_grid_x,convF_x,grid,mult_fac,  GMX_FFT_COMPLEX_TO_REAL);
        do_fft(Kern->quantity_on_grid_y,convF_y,grid,mult_fac,  GMX_FFT_COMPLEX_TO_REAL);
        do_fft(Kern->quantity_on_grid_z,convF_z,grid,mult_fac,  GMX_FFT_COMPLEX_TO_REAL);

        (*Emean)[XX] = 0.0;
        (*Emean)[YY] = 0.0;
        (*Emean)[ZZ] = 0.0;
        for (i = 0;i<grid[0];i++)
        {
                for (j = 0;j<grid[1];j++)
                {
                        sfree(Qmatr_x[i][j]);
                        sfree(qF_x[i][j]);
                        sfree(convF_x[i][j]);
                        sfree(Qmatr_y[i][j]);
                        sfree(qF_y[i][j]);
                        sfree(convF_y[i][j]);
                        sfree(Qmatr_z[i][j]);
                        sfree(qF_z[i][j]);
                        sfree(convF_z[i][j]);
//                        for (k = 0; k < grid[2];k++)
//                        { 
//                            (*Emean)[XX] += Kern->quantity_on_grid_x[i][j][k];
//                            (*Emean)[YY] += Kern->quantity_on_grid_y[i][j][k];
//                            (*Emean)[ZZ] += Kern->quantity_on_grid_z[i][j][k];
//                            printf("emean %f %f %f\n",(*Emean)[XX],(*Emean)[YY],(*Emean)[ZZ]);
//                        }

                }
                sfree(Qmatr_x[i]);
                sfree(qF_x[i]);
                sfree(convF_x[i]);
                sfree(Qmatr_y[i]);
                sfree(qF_y[i]);
                sfree(convF_y[i]);
                sfree(Qmatr_z[i]);
                sfree(qF_z[i]);
                sfree(convF_z[i]);
        }
        sfree(Qmatr_x);
        sfree(qF_x);
        sfree(convF_x);
        sfree(Qmatr_y);
        sfree(qF_y);
        sfree(convF_y);
        sfree(Qmatr_z);
        sfree(qF_z);
        sfree(convF_z);

	// The correction due to dielectric boundary conditions must now be multiplied by some constants.
	if (eps > 0.0)
	{
		// A negative epsilon (which is unphysical) just means that we don't want to apply these
		// boundary conditions.
		scalfac = 4.0*M_PI * invvol / (1.0 + 2.0*eps);
		for (i=0;i<3;i++)
		{
			dielectric[i] = scalfac*dielectric[i];
		}
		// We've calculated this correction; add it to the calculated electric field.
		for (i=0;i<grid[0];i++)
		{
			for (j=0;j<grid[1];j++)
			{
				for (k=0;k<grid[2];k++)
				{
					Kern->quantity_on_grid_x[i][j][k] += dielectric[0];
					Kern->quantity_on_grid_y[i][j][k] += dielectric[1];
					Kern->quantity_on_grid_z[i][j][k] += dielectric[2];
				}
			}
		}
	}


//        (*Emean)[XX] /=(grid[0]*grid[1]*grid[2]);
//        (*Emean)[YY] /=(grid[0]*grid[1]*grid[2]);
//        (*Emean)[ZZ] /=(grid[0]*grid[1]*grid[2]);


//compute field and subtract its mean value. now just as a test
/*
        for (i = 0;i<grid[0];i++)
        {
                for (j = 0;j<grid[1];j++)
                {
                        for (k = 0; k < grid[2];k++)
                        {
                            Kern->quantity_on_grid_x[i][j][k] -= (*Emean)[XX];
                            Kern->quantity_on_grid_y[i][j][k] -= (*Emean)[YY];
                            Kern->quantity_on_grid_z[i][j][k] -= (*Emean)[ZZ];
//                            printf("emean %f %f %f\n",(*Emean)[XX],(*Emean)[YY],(*Emean)[ZZ]);
                        }

                }
        }
*/
}

int index_wrap(int idx, int wrap)
{
        int ans;
        ans = idx;
        while (ans <0){ans +=wrap;}
        while (ans >=wrap){ans-=wrap;}
        return ans;
}

real Bspline(real x,int n)
{
	// Calculate B-spline function of order n.

	// We save the values of the coefficients between calls of this function to avoid recalculating them.
	// NBP: This may not be best practice, but currently is just as ported from Fortran.
	static real **coeffs;
	int i;
	static int init;
	static int order;
	real ans = 0, xpr;

	if (x > (real)n || x < 0.0) {return ans;}

	/* For SPME, we may need two types of spline function: if order n is used for calculating the energy (or potential), then
	we would need order n-1 for calculating the forces (or electric field). We initialize both sets of coefficients to be on
	the safe side. */
	if (init==0)
	{
		order = n;
		snew(coeffs,2);
		snew(coeffs[1],n+1);
		snew(coeffs[0],n);
		for (i = 0; i <= n; i++)
		{
			// Coefficients for order n
			coeffs[1][i] = pow(-1.0,i);
			coeffs[1][i] *= combi(n,i);
			coeffs[1][i] /= factorial(n-1);
		}
		for (i = 0; i < n; i++)
		{
			// Coefficients for order n-1
			coeffs[0][i] = pow(-1.0,i);
			coeffs[0][i] *= combi(n-1,i);
			coeffs[0][i] /= factorial(n-2);
		}
		init = 1;
	}

	if (n == order)
	{
		for (i = 0; i<=n;i++)
		{
			xpr = x - i;
			if (xpr<0.0){xpr = 0.0;}
			ans += coeffs[1][i]*pow(xpr,n-1);
		}
	}
	else if (n == order-1)
	{
		for (i = 0; i<=n;i++)
		{
			xpr = x - i;
			if (xpr<0.0){xpr = 0.0;}
			ans += coeffs[0][i]*pow(xpr,n-1);
		}
	}
	else
	{
		fprintf(stderr,"Incorrect order passed: %i. Should be %i or %i!\n",n,order,order-1);
		exit(-1);
	}
	return ans;
}

unsigned long factorial(unsigned long f)
{
    if ( f == 0 ) 
        return 1;
    return(f * factorial(f - 1));
}

long long combi(int n,int k)
{
    long long ans=1;
    k=k>n-k?n-k:k;
    int j=1;
    for(;j<=k;j++,n--)
    {
        if(n%j==0)
        {
            ans*=n/j;
        }else
        if(ans%j==0)
        {
            ans=ans/j*n;
        }else
        {
            ans=(ans*n)/j;
        }
    }
    return ans;
}

void do_fft(real ***rmatr,t_complex ***kmatr,int *dims, real multiplication_factor, int fwbck)
{
	// Do fast Fourier transform. We pass in a real array and a complex array, and the integer fwbck tells us
	// whether we are transforming backwards or forwards (i.e., it tells us which is the target array and which
	// is input).
	// OK, quite clearly gromacs doesn't want to cooperate. No problem. Let's do it my way.
	static gmx_fft_t fftF_2d, fftF_1d;
	static int plan_made = 0;
	real *rdat;
	t_complex *kdat,*k_in, *k_out, *kmatr_temp;
	int i,j,k,m, kdims[3] = {dims[0], dims[1], dims[2]/2 + 1};

	rdat = (real*) malloc(sizeof(real)*dims[1]*dims[2]);
	kdat = (t_complex*) malloc(sizeof(t_complex)*kdims[1]*kdims[2]);
	k_in = (t_complex*) malloc(sizeof(t_complex) * kdims[0]);
	k_out= (t_complex*) malloc(sizeof(t_complex) * kdims[0]);
        kmatr_temp = (t_complex*) malloc(sizeof(t_complex) * kdims[0]*kdims[1]*kdims[2]);
	if (plan_made == 0)
	{
		gmx_fft_init_2d_real(&fftF_2d,dims[1],dims[2],GMX_FFT_FLAG_NONE);
		gmx_fft_init_1d(&fftF_1d,dims[0],GMX_FFT_FLAG_NONE);
		plan_made = 1;
	}

	if (fwbck == GMX_FFT_REAL_TO_COMPLEX)
	{
		// Do real-to-complex FFT.
		// We start off by doing multiple 2d real-to-complex FFTs.

		for (i=0;i<dims[0];i++)
		{
			for (j=0;j<dims[1];j++) {
				for (k=0;k<dims[2];k++) {rdat[j*dims[2] + k] = rmatr[i][j][k];}
			}
			// Now we have this slice of the 3d array, we would like to carry out a real-to-complex transform on it.
			gmx_fft_2d_real(fftF_2d,GMX_FFT_REAL_TO_COMPLEX,rdat,kdat);
	
			// Put the result into the larger matrix (which will eventually contain the final answer).
			m=0;
			for (j=0;j<kdims[1];j++) {
				for (k=0;k<kdims[2];k++) {
                                        kmatr_temp[i*kdims[1]*kdims[2] + j*kdims[2] + k].re = kdat[m].re;
                                        kmatr_temp[i*kdims[1]*kdims[2] + j*kdims[2] + k].im = kdat[m].im;
					m++;
				}
			}
		}

		// OK, now we do a complex-to-complex FT along the remaining direction to get the final result.

		for (j=0;j<kdims[1];j++)
		{
			for (k=0;k<kdims[2];k++)
			{
				for (i=0;i<kdims[0];i++)
                                {
                                   k_in[i].re = kmatr_temp[i*kdims[1]*kdims[2] + j*kdims[2] + k].re;
                                   k_in[i].im = kmatr_temp[i*kdims[1]*kdims[2] + j*kdims[2] + k].im;
                                }
				// OK, now we have a one-dimensional array to FFT.
				gmx_fft_1d(fftF_1d,GMX_FFT_FORWARD,k_in,k_out);

				// Great - now, put the result back into the larger matrix.
				for (i=0;i<kdims[0];i++) 
                                {
                                    kmatr[i][j][k].re = multiplication_factor*k_out[i].re;
                                    kmatr[i][j][k].im = multiplication_factor*k_out[i].im;
                                }
			}
		}
	} else if (fwbck == GMX_FFT_COMPLEX_TO_REAL)
	{
		// Do complex-to-real FFT.
		// Firstly, go along one direction and do complex-to-complex FTs.

		for (j=0;j<kdims[1];j++)
		{
			for (k=0;k<kdims[2];k++)
			{
				for (i=0;i<kdims[0];i++) {k_in[i].re = kmatr[i][j][k].re;k_in[i].im = kmatr[i][j][k].im;}
				// OK, now we have a one-dimensional array to FFT.
				int status = gmx_fft_1d(fftF_1d,GMX_FFT_BACKWARD,k_in,k_out);
	
				// Great - now, put the result back into the larger matrix.
				for (i=0;i<kdims[0];i++) 
                                {
                                    kmatr_temp[i*kdims[1]*kdims[2] + j*kdims[2] + k].re = k_out[i].re;
                                    kmatr_temp[i*kdims[1]*kdims[2] + j*kdims[2] + k].im = k_out[i].im;
                                }
			}
		}

		// We finish off by doing multiple 2d complex-to-real FFTs.
		for (i=0;i<kdims[0];i++)
		{
			for (j=0;j<kdims[1];j++) {
				for (k=0;k<kdims[2];k++) {
					kdat[j*kdims[2] + k].re = kmatr_temp[i*kdims[1]*kdims[2] + j*kdims[2] + k].re;
					kdat[j*kdims[2] + k].im = kmatr_temp[i*kdims[1]*kdims[2] + j*kdims[2] + k].im;
				}
			}
			// Now we have this slice of the 3d array, we would like to carry out a complex-to-real transform on it.
			gmx_fft_2d_real(fftF_2d,GMX_FFT_COMPLEX_TO_REAL,kdat,rdat);

			// Put the result into the larger matrix (which will now contain the final answer).
			m=0;
			for (j=0;j<dims[1];j++) {
				for (k=0;k<dims[2];k++) {
					rmatr[i][j][k] = multiplication_factor*rdat[m];
					m++;
				}
			}
		}
	} else
	{
		fprintf(stderr,"Incorrect flag passed to do_fft: %i\n",fwbck);
		exit(-1);
	}
        sfree(rdat); sfree(kdat), sfree(kmatr_temp);
        sfree(k_in);  sfree(k_out);
}
 
int gmx_eshs(int argc, char *argv[])
{
    const char        *desc[] = {
        "The structure of liquids can be studied by elastic second harmonic scattering.",
        "[THISMODULE] calculates the non-linear optical scattering intensity in 2 different ways.",
        "The simplest method (single, default) is to use 1/N<|sum_i beta_IJK(i) exp[iq dot r_i]|^2>.",
        "This however converges slowly with the simulation time and is more prone to noise.[PAR]",
        "A method that converges more quickly wrt simulation time (double), but nmol times more expensive, is",
        "I(q)=1/N<sum_i beta_IJK(i)^2> + 1/N< sum_i sum_{i!=j} beta_IJK(i)beta_IJK(j) cos(q dot r_ij) >.",
        "I(q)=incoherent term + coherent term. For more details see Bersohn, et al. JCP 45, 3184 (1966)",
        "pout, pin, are the angles formed between the polarization vectors and the scattering plane.",
        "Common polarization combinations are PSS, PPP, SPP, SSS . [PAR]",
    };
    static gmx_bool          bPBC = TRUE, bIONS = FALSE;
    static real              electrostatic_cutoff = 1.2, maxelcut = 2.0, kappa = 5.0,  kernstd = 10.0 ;
    static real              fspacing = 0.01, pout_angle = 0.0 , pin_angle = 0.0, std_dev_dens = 0.05;
    static real              binwidth = 0.002, angle_corr = 90.0, eps = -1.0 , kmax_spme = 4.0;
    static int               ngroups = 1, nbintheta = 10, nbingamma = 2 ,qbin = 1, nbinq = 10 ;
    static int               nkx = 0, nky = 0, nkz = 0, kern_order = 2, interp_order = 4, kmax = 0;
		static real							 kappa2 = -1.0, ecorrcut = -1.0;

    static const char *methodt[] = {NULL, "single", "double" ,NULL };
    static const char *kernt[] = {NULL, "krr", "scalar", "none", "map", NULL};
    static char *catname = NULL;
    static char *anname =  NULL;

    t_pargs            pa[] = {
        { "-nbintheta",     FALSE, etINT, {&nbintheta},
        "number of bins over scattering angle theta chosen between -pi/2 and + pi/2 (available only with thetaswipe)" },
        { "-nplanes",       FALSE, etINT, {&nbingamma},
        "number of scattering planes that lie on the scattered wave-vector to average over, -PI/2< gamma< PI/2" },
        { "-angle_corr",       FALSE, etREAL, {&angle_corr},
        "angle at which to compute <beta(0)beta(r)>" },
        { "-qbin",          FALSE, etINT, {&qbin},
        "choose wave-vector to sample given by 2pi/box-length*qbin" },
        { "-nbinq",         FALSE, etINT, {&nbinq},
        "how many bins in the reciprocal space" },
        { "-stddens",       FALSE, etREAL, {&std_dev_dens}, "standard deviation to compute density on a grid. Use only with scalar kernel [nm]."},
        { "-fourierspacing",          FALSE, etREAL, {&fspacing}, "fourier spacing [nm] gives lower bound for number of wave vectors to use in each direction with Ewald, overridden by nkx,nky,nkz " },
        { "-binw",          FALSE, etREAL, {&binwidth}, "width of bin to compute <beta_lab(0) beta_lab(r)> " },
        { "-pout",          FALSE, etREAL, {&pout_angle}, "polarization angle of outcoming beam in degrees. For P choose 0, for S choose 90" },
        { "-pin",           FALSE, etREAL, {&pin_angle}, "polarization angle of incoming beam in degrees. For P choose 0, for S choose 90" },
        { "-cutoff",        FALSE, etREAL, {&electrostatic_cutoff}, "cutoff for the calculation of electrostatics around a molecule and/or for method=double" },
        { "-maxcutoff",        FALSE, etREAL, {&maxelcut}, "cutoff to smoothly truncate the calculation of the double sum" },
        { "-kappa",        FALSE, etREAL, {&kappa}, "screening parameter for the ewald term, i.e. erf(r*kappa)/r, in nm^-1" },
        { "-kernorder",        FALSE, etINT, {&kern_order}, "kernel order, where beta = sum_i sum_(kern_ind) c[i*kern_ind] * (feature_vec[i]-mean(feature_vec))^kern_ind " },     
        { "-splorder",        FALSE, etINT, {&interp_order}, "interpolation order for b-splines" },
        { "-kmax_spme",        FALSE, etREAL, {&kmax_spme}, "maximum wave-vector to use when performing th SPME in fourier space. This gives the number of images used to compute spme " },
        { "-kernstd",       FALSE, etREAL, {&kernstd}, "standard deviation of kernel function. only makes sense if kernel ridge regression is used"},

        { "-method",     FALSE, etENUM, {methodt}, "I(q) using a single summation O(N) or a double summation, O(N*N)" },
        { "-kern",   FALSE, etENUM, {kernt}, "what method to use to compute beta"},
        { "-ions",   FALSE, etBOOL, {&bIONS}, "compute molecular hyperpolarizability when ions are present"},
        { "-cn",     FALSE, etSTR, {&catname}, "name of cation"},
        { "-an",     FALSE, etSTR, {&anname}, "name of anion"},
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances. Always use, results without PBC not tested." },
        { "-ng",       FALSE, etINT, {&ngroups}, 
          "Number of secondary groups, not available for now. Only tip4p water implemented." },
        { "-eps",	FALSE, etREAL, {&eps}, "dielectric constant"},
        { "-kappa2", FALSE, etREAL, {&kappa2}, "large kappa for real-space correction."},
				{ "-ecorrcut", FALSE, etREAL, {&ecorrcut}, "cutoff length for electric field correction."},
    };
#define NPA asize(pa)
    const char        *fnTPS, *fnNDX , *fnBETACORR = NULL, *fnFTBETACORR= NULL, *fnREFMOL = NULL;
    const char        *fnVCOEFF=NULL, *fnVGRD=NULL, *fnVINP=NULL;
    const char        *fnRGRDO=NULL, *fnRGRDH=NULL, *fnCOEFFO=NULL, *fnCOEFFH=NULL;
    const char        *fnMAP=NULL;
    output_env_t       oenv;
    int           *gnx;
    int            nFF[2];
    atom_id      **grpindex;
    char         **grpname = NULL;
    /*gmx_bool       bGkr, bMU, bSlab;*/

    t_filenm           fnm[] = {
        { efTRX, "-f",  NULL,     ffREAD },
        { efMAP, "-emap",    "static.map",   ffOPTRD },
        { efDAT, "-vinp",    "vinput.dat", ffOPTRD},
        { efDAT, "-vgrid",   "vgrid.dat", ffOPTRD},
        { efDAT, "-vcoeff",  "vcoeff.dat", ffOPTRD},
        { efDAT, "-rhogridO",   "rhogridO.dat", ffOPTRD},
        { efDAT, "-rhogridH",   "rhogridH.dat", ffOPTRD},
        { efDAT, "-rhocoeffH",   "rhocoeffH.dat", ffOPTRD},
        { efDAT, "-rhocoeffO",   "rhocoeffO.dat", ffOPTRD},
        { efDAT, "-refmol",  "refmol.dat", ffOPTRD},
        { efDAT, "-ewlog",   "ewaldlog.dat",ffOPTWR},
        { efTPS, NULL,  NULL,     ffREAD },
        { efNDX, NULL,  NULL,     ffOPTRD },
        { efXVG, "-o",  "non_linear_sfact",    ffWRITE },
        { efXVG, "-otheta", "non_linear_sfact_vs_theta", ffOPTWR },
        { efXVG, "-betacorr", "beta_correlation", ffOPTWR },
        { efXVG, "-ftbetacorr", "FT_beta_correlation", ffOPTWR },

    };
#define NFILE asize(fnm)
    int            npargs;
    t_pargs       *ppa;
    t_topology    *top;
    int            ePBC;
    int            k, natoms;
    matrix         box;
    
    npargs = asize(pa);
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                           NFILE, fnm, NPA, pa, asize(desc), desc, 0, NULL, &oenv))
    {
        return 0;
    }

    fnTPS = ftp2fn_null(efTPS, NFILE, fnm);
    fnNDX = ftp2fn_null(efNDX, NFILE, fnm);
    fnMAP = opt2fn_null("-emap", NFILE,fnm);
    fnVCOEFF = opt2fn_null("-vcoeff", NFILE,fnm);
    fnVGRD = opt2fn_null("-vgrid", NFILE,fnm);
    fnVINP = opt2fn_null("-vinp", NFILE,fnm);
    fnRGRDO = opt2fn_null("-rhogridO", NFILE,fnm);
    fnRGRDH = opt2fn_null("-rhogridH", NFILE,fnm);
    fnCOEFFO = opt2fn_null("-rhocoeffO", NFILE,fnm);
    fnCOEFFH = opt2fn_null("-rhocoeffH", NFILE,fnm);
    fnBETACORR = opt2fn_null("-betacorr", NFILE,fnm);
    fnFTBETACORR = opt2fn_null("-ftbetacorr", NFILE,fnm);
    fnREFMOL = opt2fn_null("-refmol", NFILE, fnm);


    if (!fnTPS && !fnNDX)
    {
        gmx_fatal(FARGS, "Neither index file nor topology file specified\n"
                  "Nothing to do!");
    }

    if ((*kernt)[0] == 's')
    {
       if (!fnVCOEFF || !fnVGRD || !fnRGRDO || !fnRGRDH || !fnCOEFFO || !fnCOEFFH)
       {
          gmx_fatal(FARGS, "specify all files for scalar kernel using -vcoeff, -vgrid, -vinp, -rhogridO, -rhogridH, -rhocoeffO, rhocoeffH\n");
       }
    }
    else if ((*kernt)[0] == 'm')
    {
       if (!fnMAP )
       {
          gmx_fatal(FARGS, "specify map file with -emap\n");
       }
    }
    else if ((*kernt)[0] == 'k' )
    {
       if (!fnVCOEFF || !fnVGRD || !fnVINP)
       {
          gmx_fatal(FARGS, "specify all files for krrr using -vcoeff, -vgrid, -vinp\n");
       }
    }

    snew(top, ngroups+1);
    ePBC = read_tpx_top(ftp2fn(efTPS, NFILE, fnm), NULL, box,
                        &natoms, NULL, NULL, NULL, top);

    snew(gnx, ngroups+1);
    snew(grpname, ngroups+1);
    snew(grpindex, ngroups+1);
    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm),
             ngroups +1 , gnx, grpindex, grpname);

      // If no kmax is specified, then one is chosen according to the "optimizing" formula of Frenkel and Smit.
      //kmax = (int)(M_PI/(min(0.5,1.2*pow(natoms,-1.0/6.0))));
      //kmax = (int)(2.5*box[XX][XX]*kappa/M_PI);
    kmax = roundf(kmax_spme*box[XX][XX]);
    fprintf(stderr,"\n setting kmax_spme = %i\n",kmax);

    fprintf(stderr," Start indexing the atoms to each molecule\n");
    dipole_atom2mol(&gnx[0], grpindex[0], &(top->mols));

		if (ecorrcut < 0.0)
		{
			// No cutoff length has been chosen for the electric field correction term. We will choose one here.
			ecorrcut = 1.1 * (sqrt(1.0/kappa) + sqrt(1.0/kappa2));
			fprintf(stderr,"Electric field cutoff chosen: %f nm\n",ecorrcut);
//		ecorrcut = 1.120944;
//		The above value (in nm) was suitable. This is given fairly well by the formula above.
		}

		real *sigma_vals;
		snew(sigma_vals,4);
		sigma_vals[0] = kappa2;
		sigma_vals[1] = kappa;
		sigma_vals[2] = kappa2*kappa2;
		sigma_vals[3] = kappa*kappa;

    do_eshs(top, ftp2fn(efTRX, NFILE, fnm),
            opt2fn("-o", NFILE, fnm), opt2fn("-otheta", NFILE, fnm), angle_corr,
           fnVCOEFF, fnVGRD, fnVINP, fnRGRDO, fnCOEFFO,
           fnRGRDH, fnCOEFFH, fnMAP, fnBETACORR, fnFTBETACORR,
           fnREFMOL, methodt[0], kernt[0], bIONS, catname, anname, bPBC,  qbin, nbinq,
           kern_order, std_dev_dens, fspacing, binwidth,
           nbintheta, nbingamma, pin_angle, pout_angle, 
           electrostatic_cutoff, maxelcut, kappa, interp_order, kmax, kernstd, gnx, grpindex, grpname, ngroups, oenv, eps, sigma_vals, ecorrcut);
    return 0;
}

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
#include "gromacs/math/do_fit.h"
#include "mtop_util.h"
#include "typedefs.h"
#include "force.h"

#include "gromacs/legacyheaders/gmx_fatal.h"

#define nm2Bohr 10.0/0.529177
#define ANG2NM 0.1
#define AU2VNM 5.14220652e2



static void do_shsint(t_topology *top,  const char *fnTRX,
                   const char *fnSFACT, const char *fnTHETA, const real angle_corr,
                   const char *fnMAP, const char *fnBETAINP,
                   const char *fnBETACORR, const char *fnFTBETACORR, const char *fnREFMOL,
                   const char *method, const char *kern,
                   gmx_bool bIONS, char *catname, char *anname, gmx_bool bPBC, 
                   int qbin, int nbinq, 
                   real binwidth, int nbintheta, int nbingamma, real pin_angle, real pout_angle,
                   int *isize, int  *molindex[], char **grpname, int ng,
                   const output_env_t oenv)
{
    FILE          *fp, *fpn, *fnINP;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, nanions, ncations, i, j, k, qq, n, c, tt, rr, bytes_beta, nframes_betainp, nframes, nfaces, gr_ind, nbin, aa, bb, cc;
    real         **s_method, **s_method_coh, **s_method_incoh, *temp_method, ****s_method_t, ****s_method_coh_t, ****s_method_incoh_t, ***mu_sq_t, ***coh_temp;
    real           qnorm, maxq, incoh_temp = 0.0, tot_temp = 0.0, gamma = 0.0 ,theta0 = 5.0, check_pol;
    real          *cos_t, *sin_t, ****cos_tq, ****sin_tq,   mu_ind =0.0, mu_sq =0.0, mod_f ;
    real         **field_ad, electrostatic_cutoff2, electrostatic_cutoff, max_spacing, maxelcut2,  invkappa2, betainp_value, ***beta_mol, *betamean ,*mu_ind_mols, ****mu_ind_t, *beta_corr, *ft_beta_corr;
    int            max_i, isize0, ind0, indj;
    real           t, rmax2, rmax,  r, r_dist, r2, q_xi, dq;
    real          *inv_segvol, normfac, segvol, spherevol, prev_spherevol, invsize0, invgamma, invhbinw, inv_width,  theta=0, *theta_vec;
    rvec          *x, xcm, xcm_transl, dx,  *x_i1, xi, x01, x02, *arr_qvec, **arr_qvec_faces ,vec_polin, vec_polout, ***vec_pout_theta_gamma, ***vec_pin_theta_gamma;
    rvec           pol_perp, pol_par,vec_kout, vec_2kin, pol_in1, pol_in2, vec_kout_2kin ;
    rvec           xvec, yvec, zvec, *xmol, *xref, Emean;
    real          *qref;
    matrix         cosdirmat,invcosdirmat; 
    real           invvol, invvol_sum;
    t_Map         *Map=NULL;
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
    if ( kern[0] == 'n' )
    {
       Map=(t_Map *)calloc(1,sizeof(t_Map));
       readMap(fnMAP, Map);
       fprintf(stderr,"initialized electric field map to compute beta\n");
       
    }
    else if (kern[0] == 's')
    {
       fprintf(stderr,"you have chosen the scalar kernel\nnow I will read the file containing the beta values for each molecule\n");
       fnINP=fopen(fnBETAINP,"rb");
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

        fprintf(stderr, "rmax2 = %f\n", rmax2);
        if (fnBETACORR)
        {
            fprintf(stderr, "number of bins for <beta(0)*beta(r)> = %d\n", nbin);
            nfaces = 1;
            if (nbingamma >1 || nbintheta >1  )
            {
                gmx_fatal(FARGS, "when computing <beta(0)*beta(r)> choose nplanes = 1 and nbintheta = 1");
            }
        } 
        inv_width = 1.0  ;
        
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
    nframes_betainp = 0; 
    betainp_value = 0;

//    nframes_betainp = (int *)calloc(1,sizeof(int));
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
            if (kern[0] == 's')
            {
               bytes_beta = fread(&nframes_betainp,sizeof(int),1,fnINP);                        fprintf(stderr,"frame n. of betafile %d\n",nframes_betainp);
            }
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
                    if (kern[0] == 's' )
                    {
                      for (aa = 0; aa < DIM; aa++)
                      {
                          for (bb = 0; bb < DIM; bb++)
                          {
                              for (cc = 0; cc < DIM; cc++)
                              {
                                 bytes_beta = fread(&betainp_value,sizeof(float),1,fnINP);
                                 beta_mol[aa][bb][cc] = betainp_value ;
                                 printf("beta %d %d %d %f\n",aa, bb, cc, beta_mol[aa][bb][cc]);
                              }
                          }
                      }

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
}

int gmx_shsint(int argc, char *argv[])
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
    static real              pout_angle = 0.0 , pin_angle = 0.0;
    static real              binwidth = 0.002, angle_corr = 90.0;
    static int               ngroups = 1, nbintheta = 10, nbingamma = 2 ,qbin = 1, nbinq = 10 ;

    static const char *methodt[] = {NULL, "single", NULL };
    static const char *kernt[] = {NULL, "scalar", "none", NULL};
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
        { "-binw",          FALSE, etREAL, {&binwidth}, "width of bin to compute <beta_lab(0) beta_lab(r)> " },
        { "-pout",          FALSE, etREAL, {&pout_angle}, "polarization angle of outcoming beam in degrees. For P choose 0, for S choose 90" },
        { "-pin",           FALSE, etREAL, {&pin_angle}, "polarization angle of incoming beam in degrees. For P choose 0, for S choose 90" },

        { "-method",     FALSE, etENUM, {methodt}, "I(q) using a single summation O(N) or a double summation, O(N*N)" },
        { "-kern",   FALSE, etENUM, {kernt}, "what method to use to compute beta"},
        { "-ions",   FALSE, etBOOL, {&bIONS}, "compute molecular hyperpolarizability when ions are present"},
        { "-cn",     FALSE, etSTR, {&catname}, "name of cation"},
        { "-an",     FALSE, etSTR, {&anname}, "name of anion"},
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances. Always use, results without PBC not tested." },
        { "-ng",       FALSE, etINT, {&ngroups}, 
          "Number of secondary groups, not available for now. Only tip4p water implemented." },
    };
#define NPA asize(pa)
    const char        *fnTPS, *fnNDX , *fnBETACORR = NULL, *fnFTBETACORR= NULL, *fnREFMOL = NULL;
    const char        *fnBETAINP=NULL;
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
        { efDAT, "-refmol",  "refmol.dat", ffOPTRD},
        { efDAT, "-betainp",  "betafile.dat", ffOPTRD},
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
    fnBETAINP = opt2fn_null("-betainp", NFILE,fnm);
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
       if (!fnBETAINP)
       {
          gmx_fatal(FARGS, "specify input file containing the trajectory of molecular betas with -betainp\n");
       }
    }
    else if ((*kernt)[0] == 'n')
    {
       if (!fnMAP )
       {
          gmx_fatal(FARGS, "specify map containing the constant beta values\n");
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

    fprintf(stderr,"Start indexing the atoms to each molecule\n");
    dipole_atom2mol(&gnx[0], grpindex[0], &(top->mols));

    do_shsint(top, ftp2fn(efTRX, NFILE, fnm),
           opt2fn("-o", NFILE, fnm), opt2fn("-otheta", NFILE, fnm), angle_corr,
           fnMAP, fnBETAINP, fnBETACORR, fnFTBETACORR,
           fnREFMOL, methodt[0], kernt[0], bIONS, catname, anname, bPBC,  qbin, nbinq,
           binwidth, nbintheta, nbingamma, pin_angle, pout_angle, 
           gnx, grpindex, grpname, ngroups, oenv);
    return 0;
}

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
#include "typedefs.h"
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

#include "gromacs/legacyheaders/gmx_fatal.h"

static void do_nonlinearopticalscatteringtheta(t_topology *top, /*const char *fnNDX, const char *fnTPS,*/ const char *fnTRX,
                   const char *fnSFACT, const char *fnTHETA, const char *method,
                   gmx_bool bPBC, gmx_bool bKleinmannsymm, gmx_bool bSpectrum , gmx_bool bRandbeta, gmx_bool bThetaswipe,
                   int qbin, int nbinq, real koutx, real kouty, real koutz,
                   real kinx, real kiny, real kinz, real  bmzxx, real  bmzyy, real bmzzz,
                   int nbintheta, int nbingamma ,real pin_angle, real pout_angle ,
                   real fade, int *isize, int  *molindex[], char **grpname, int ng,
                   const output_env_t oenv)
{
    FILE          *fp;
    FILE          *fpn;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, i, j, k, qq, n, c, tt, rr, nframes, nfaces;
    real         **s_method, **s_method_coh, **s_method_incoh, *temp_method, ****s_method_t, ****s_method_coh_t, ****s_method_incoh_t, ***mu_sq_t ;
    real           qnorm, maxq, coh_temp = 0.0,  incoh_temp = 0.0, tot_temp = 0.0, gamma = 0.0 ,theta0 = 5.0, check_pol;
    real          *cos_t, *sin_t, ****cos_tq, ****sin_tq, *cq, *sq , *c_square, *s_square, *cs_sc ,***beta_mol, ***beta_mol_const, ***beta_const_dev, *beta_mol_1d, *beta_lab_2_t, *beta_lab_1_t, beta_fact, mu_ind =0.0, mu_sq =0.0 ;
    real           beta_lab_sq_2 = 0.0, beta_lab_sq_1 = 0.0, beta_lab_1_2 =0.0, beta_lab_2 = 0.0, beta_lab_1 = 0.0, b22 = 0.0, b21 = 0.0, b12 = 0.0, b11 = 0.0;
    int            max_i, isize0, ind0;
    real           t, rmax2, rmax,  r, r_dist, r2, q_xi, dq, invhbinw, normfac, norm_x, norm_z, mod_f, inv_width;
    real           segvol, spherevol, prev_spherevol, invsize0, invgamma, theta=0, *theta_vec;
    rvec          *x, dx,  *x_i1, xi, x01, x02, *arr_qvec, **arr_qvec_faces ,vec_polin, vec_polout, ***vec_pout_theta_gamma, ***vec_pin_theta_gamma;
    rvec           pol_perp, pol_par,vec_kout, vec_2kin, pol_in1, pol_in2, vec_kout_2kin ; 
    real           invvol, invvol_sum, rho;
    matrix         box, box_pbc;
    int            ePBC = -1, ePBCrdf = -1;
    t_block       *mols = NULL;
    t_atom        *atom = NULL;
    t_pbc          pbc;
    gmx_rmpbc_t    gpbc = NULL;
    gmx_rng_t      rng = NULL;
    int            mol, a;

    atom = top->atoms.atom;
    mols = &(top->mols);
    isize0 = isize[0];
    nfaces = 6;
    invsize0 = 1.0/isize0;
    invgamma = 1.0/nbingamma;
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    
    fprintf(stderr,"\nnumber of atoms %d\n",natoms);

    /* initialize some handy things */
    if (ePBC == -1)
    {
        ePBC = guess_ePBC(box);
    }
    copy_mat(box, box_pbc);
    ePBCrdf = ePBC;
    if (bPBC)
    {
        rmax2   =  max_cutoff2(FALSE ? epbcXY : epbcXYZ, box_pbc);
        fprintf(stderr, "rmax2 = %f\n", rmax2);
        
    }
    else
    {
        rmax2   = sqr(3*max(box[XX][XX], max(box[YY][YY], box[ZZ][ZZ])));
    }
    if (debug)
    {
        fprintf(debug, "rmax2 = %g\n", rmax2);
    }
    
    /*here find the normalization of the OH bond-length and of the dipole vector*/
    set_pbc(&pbc, ePBCrdf, box_pbc);
    ind0  = mols->index[molindex[0][0]];
    pbc_dx(&pbc, x[ind0+1], x[ind0], x01);
    pbc_dx(&pbc, x[ind0+2], x[ind0], x02);
    rvec_sub( x01, x02, xi);
    norm_x = gmx_invsqrt(norm2(xi));
    rvec_add( x01, x02, xi);
    norm_z = gmx_invsqrt(norm2(xi));
    rmax     = sqrt(rmax2);

    snew(s_method, ng);
    snew(s_method_coh, ng);
    snew(s_method_incoh, ng);
    snew(s_method_t, ng);
    snew(s_method_coh_t, ng);
    snew(s_method_incoh_t, ng);
    max_i = isize[0];

    /*allocate memory for beta in lab frame and initialize beta in mol frame*/
    snew(beta_lab_2_t, isize[0]);
    snew(beta_lab_1_t, isize[0]);
    snew(beta_mol, DIM);
    snew(beta_mol_const, DIM);
    snew(beta_const_dev, DIM);
    snew(beta_mol_1d, 7);
    for (i = 0; i < DIM; i++)
    {
        snew(beta_mol[i], DIM);
        snew(beta_mol_const[i], DIM);
        snew(beta_const_dev[i], DIM);
    }
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            snew(beta_mol[i][j], DIM);
            snew(beta_mol_const[i][j], DIM);
            snew(beta_const_dev[i][j], DIM);
        }
    }

    /*beta_mol parameters in a.u. taken from Kusalik Mol. Phys. 99 1107-1120, (2001)*/
    /*convert from a.u. to [e^3]*[nm^3]/[Hartree^2]*/
    beta_mol[2][0][0] = bmzxx ; //0.22 ; //5.7 /**0.000148184711*/;
    beta_mol[2][2][2] = bmzzz ; //12.8 ; //31.6 /**0.000148184711*/;
    beta_mol[2][1][1] = bmzyy ; //7.5  ; //10.9 /**0.000148184711*/;
    beta_mol_const[2][0][0] = bmzxx ; //0.22 ; //5.7 /**0.000148184711*/;
    beta_mol_const[2][2][2] = bmzzz ; //12.8 ; //31.6 /**0.000148184711*/;
    beta_mol_const[2][1][1] = bmzyy ; //7.5  ; //10.9 /**0.000148184711*/;

    fprintf(stderr,"beta_mol[z][x][x] %f\n",beta_mol[2][0][0]);
    fprintf(stderr,"beta_mol[z][z][z] %f\n",beta_mol[2][2][2]);
    fprintf(stderr,"beta_mol[z][y][y] %f\n",beta_mol[2][1][1]);

    beta_mol_1d[0] = beta_mol[2][0][0];
    beta_mol_1d[1] = beta_mol[2][1][1];
    beta_mol_1d[2] = beta_mol[2][2][2];

    if (bKleinmannsymm)
    {
       beta_mol[0][0][2] = bmzxx ;
       beta_mol[0][2][0] = bmzxx ;
       beta_mol[1][1][2] = bmzyy;
       beta_mol[1][2][1] = bmzyy;
       beta_mol_const[0][0][2] = bmzxx ;
       beta_mol_const[0][2][0] = bmzxx ;
       beta_mol_const[1][1][2] = bmzyy;
       beta_mol_const[1][2][1] = bmzyy;
       beta_mol_1d[3] = beta_mol[0][0][2];
       beta_mol_1d[4] = beta_mol[0][2][0];
       beta_mol_1d[5] = beta_mol[1][1][2];
       beta_mol_1d[6] = beta_mol[1][2][1];
    }
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            for (k = 0; k < DIM; k++)
            {
                if ((i == 0 &&  j == 2 && k == 0) || 
                  (i == 0 &&  j == 0 && k == 2) ||
                  (i == 1 &&  j == 1 && k == 2) ||
                  (i == 1 &&  j == 2 && k == 1) ||
                  (i == 2 &&  j == 1 && k == 1) ||
                  (i == 2 &&  j == 0 && k == 0) ||
                  (i == 2 &&  j == 2 && k == 2))
                  {
                     beta_const_dev[i][j][k] = beta_mol_const[i][j][k];
                  }
                else
                  {
                     beta_const_dev[i][j][k] = 1.0 ;
                  }
            }
        }
    }

    
    for (g = 0; g < ng; g++)
    {
        /* this is THE array */
        /*allocate memory for s_method array */
        snew(cq,nbintheta+1);
        snew(sq,nbintheta+1);
        snew(c_square,nbintheta+1);
        snew(s_square,nbintheta+1);
        snew(cs_sc, nbintheta+1);
        inv_width = (fade == 0.0 ) ? 1.0 : M_PI*0.5/(rmax-fade) ;
        if (bSpectrum == FALSE)
        {
            snew(s_method[g], nbintheta+1);
            snew(s_method_coh[g], nbintheta+1);
            snew(s_method_incoh[g], nbintheta+1);
            snew(temp_method, nbintheta+1);
            snew(arr_qvec, nbintheta+1);
            snew(cos_t, nbintheta+1);
            snew(sin_t, nbintheta+1);
            qnorm = M_PI*2.0/(rmax*2.0)*qbin;
            fprintf(stderr,"|q| = %f\n", qnorm);
            for (qq = 0; qq <= nbintheta; qq++)
            {
                theta  = M_PI/(2*nbintheta)*qq -M_PI*0.5 ;
                arr_qvec[qq][XX] = qnorm*sin(theta);  // /( sqrt(-2*cos(theta) + 2.0)) ;
                arr_qvec[qq][YY] = 0.0 ;
                arr_qvec[qq][ZZ] = qnorm*(cos(theta) -1.0);  // /( sqrt(-2*cos(theta) + 2.0))  ; 
                cq[qq] = cos(theta) ;
                sq[qq] = sin(theta) ;
                c_square[qq] = cq[qq]*cq[qq] ;
                cs_sc[qq] = cq[qq]*sq[qq] ;
                s_square[qq] = sq[qq]*sq[qq] ;
            }
        }
        else
        {
           snew(s_method[g], nbinq);
           snew(s_method_coh[g], nbinq);
           snew(s_method_incoh[g],nbinq);
           snew(temp_method, nbinq);
           snew(arr_qvec,nbinq);
           /*initialize incoming and outcoming wave-vectors*/
           vec_kout[XX] = koutx; 
           vec_kout[YY] = kouty;
           vec_kout[ZZ] = koutz;
           vec_2kin[XX] = kinx;
           vec_2kin[YY] = kiny;
           vec_2kin[ZZ] = kinz;
           /* compute the polarization vectors for outcoming beam*/     
           cprod(vec_kout, vec_2kin, pol_perp);
           cprod(vec_kout, pol_perp, pol_par);
           fprintf(stderr,"polarization vector perpendicular to outcoming beam %f %f %f\n",pol_perp[XX], pol_perp[YY], pol_perp[ZZ]);
           fprintf(stderr,"polarization vector parallel to Outcoming beam %f %f %f\n",pol_par[XX], pol_par[YY], pol_par[ZZ]);
           svmul(sin(M_PI/180.0*pout_angle), pol_perp, pol_perp);
           svmul(cos(M_PI/180.0*pout_angle), pol_par,  pol_par);
           rvec_add(pol_perp, pol_par, vec_polout);
           /* compute the polarization vectors for incoming beam*/                  
           cprod(vec_kout, vec_2kin, pol_perp);
           cprod(vec_2kin, pol_perp, pol_par);
           fprintf(stderr,"polarization vector perpendicular to incoming beam %f %f %f\n",pol_perp[XX], pol_perp[YY], pol_perp[ZZ]);
           fprintf(stderr,"polarization vector parallel to incoming %f %f %f\n",pol_par[XX], pol_par[YY], pol_par[ZZ]);
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
           qnorm = M_PI*2.0/(rmax*2.0)*qbin;
           printf("----INITIALIZE DIRECTION OF INCOMING AND OUTCOMING WAVE-VECTORS----\n");
           printf("direction of incoming wave-vector is %f %f %f\n", vec_2kin[XX], vec_2kin[YY], vec_2kin[ZZ]);
           printf("direction of outcoming wave-vector is %f %f %f\n", vec_kout[XX], vec_kout[YY], vec_kout[ZZ]);
           printf("----INITIALIZE DIRECTION OF INCOMING AND OUTCOMING POLARIZATION VECTORS----\n");
           printf("polarization of incoming wave-vector is %f %f %f\n", vec_polin[XX], vec_polin[YY], vec_polin[ZZ]);
           printf("polarization of outcoming wave-vector is %f %f %f\n", vec_polout[XX], vec_polout[YY], vec_polout[ZZ]);
           printf("----INITIALIZE DIRECTION OF SCATTERED WAVE VECTOR: q=kout-2kin ----\n");
           printf("direction of scattered wave-vector is %f %f %f\n", vec_kout[XX]-vec_2kin[XX], vec_kout[YY]-vec_2kin[YY], vec_kout[ZZ]-vec_2kin[ZZ]);
           printf("minimum wave-vector is (2pi/L)*qbin = %f\n", qnorm);
           printf("maximum wave-vector is (2pi/L)*qbin*nbinq = %f\n", qnorm*nbinq);
           cq[0] = cos(theta0) ;
           sq[0] = sin(theta0) ;
           c_square[0] = cq[0]*cq[0] ;
           cs_sc[0] = cq[0]*sq[0] ;
           s_square[0] = sq[0]*sq[0] ;

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
  
           if (bThetaswipe == TRUE)
           {
              snew(s_method_t[g], nfaces);
              snew(s_method_coh_t[g], nfaces);
              snew(s_method_incoh_t[g], nfaces);
              snew(arr_qvec_faces,nfaces);
              snew(vec_pout_theta_gamma,nfaces);
              snew(vec_pin_theta_gamma, nfaces);
              snew(theta_vec, nbintheta);
              // this loop is to get the scattering wave-vector and polarization vectors for different faces of the cube
              printf("\n----COMPUTE THE POLARIZATION VECTORS AT DIFFERENT SCATTERING ANGLES THETA.---------------------\n");
              printf("----THETA IS DEFINED AS THE ANGLE BETWEEN THE INCOMING AND OUTCOMING WAVE-VECTORS.-------------\n");              
              printf("----THE POLARIZATION VECTORS ARE ALSO COMPUTED AT DIFFERENT ANGLES GAMMA.----------------------\n");
              printf("----GAMMA IS THE ANGLE DEFINED BY A PLANE THAT GOES THROUGH THE SCATTERING WAVE-VECTOR AND-----\n");
              printf("----THE PLANE PARALLEL TO THE CHOSEN FACE OF THE SIMULATION BOX--------------------------------\n \n");
              for (rr = 0; rr< nfaces; rr++)
              {
                 snew(s_method_t[g][rr], nbinq);
                 snew(s_method_coh_t[g][rr], nbinq);
                 snew(s_method_incoh_t[g][rr], nbinq);
                 snew(arr_qvec_faces[rr],nbinq);
                 for (qq = 0; qq< nbinq; qq++)
                 {
                      //now the magnitude of the wave-vector is indeed sqrt(2)*2pi/L*n so we don't need to multiply by sqrt(2).
                    arr_qvec_faces[rr][qq][XX] = (qnorm + qnorm*qq)*(1.0) ;
                    arr_qvec_faces[rr][qq][YY] = 0.0 ;
                    arr_qvec_faces[rr][qq][ZZ] = (qnorm + qnorm*qq)*(-1.0) ;
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
                                printf("polarization vectors at angles theta_expt = %f gamma = %f and face index %d \n",2.0*theta_vec[tt]*180.0/M_PI, gamma*180.0/M_PI, rr);
                                printf("incoming polarization vector = %f %f %f \n",vec_pin_theta_gamma[rr][tt][c][XX], vec_pin_theta_gamma[rr][tt][c][YY], vec_pin_theta_gamma[rr][tt][c][ZZ]);
                                printf("outcoming polarization vector = %f %f %f \n",vec_pout_theta_gamma[rr][tt][c][XX], vec_pout_theta_gamma[rr][tt][c][YY], vec_pout_theta_gamma[rr][tt][c][ZZ]);
                                printf("direction of scattered wave-vector = %f %f %f \n", vec_kout[XX] -vec_2kin[XX], vec_kout[YY] - vec_2kin[YY], vec_kout[ZZ] -vec_2kin[ZZ]);
                                check_pol = iprod(vec_pin_theta_gamma[rr][tt][c],vec_pout_theta_gamma[rr][tt][c]);
                                printf("(incoming polarization vec) dot (outcoming polarization vec) = %.17g\n",check_pol);
                            }
                        }
                 }
                 printf("wave_vec %f %f %f\n", arr_qvec_faces[rr][0][XX], arr_qvec_faces[rr][0][YY], arr_qvec_faces[rr][0][ZZ]);
              }
           }
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

    if (method[0] == 'm' && fade != 0.0 && bSpectrum == TRUE)
    {
       fprintf(stderr,"modified direct method (modsumexp) with spectrum\n");
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
                beta_lab_sq_1 = 0.0;
                beta_lab_sq_2 = 0.0;
                beta_lab_1_2 = 0.0 ;
                snew(temp_method, nbinq);
                for (i = 0; i < isize0; i++)
                {
                    ind0  = mols->index[molindex[g][i]];
                    copy_rvec(x[ind0], x_i1[i]);
                    pbc_dx(&pbc, x[ind0+1], x_i1[i], x01);
                    pbc_dx(&pbc, x[ind0+2], x_i1[i], x02);
                    rotate_beta_theta(x01, x02, pol_in1, pol_in2, beta_mol_1d, &beta_lab_2, &beta_lab_1 );
                    beta_lab_sq_2 += beta_lab_2*beta_lab_2 ;
                    beta_lab_sq_1 += beta_lab_1*beta_lab_1 ;
                    beta_lab_1_2  += beta_lab_1*beta_lab_2 ;
                    beta_lab_2_t[i] = beta_lab_2;
                    beta_lab_1_t[i] = beta_lab_1;
                }
                for (i = 0; i < isize0 -1.0; i++)
                {
                    /* Real rdf between points in space */
                    ind0 = mols->index[molindex[g][i]];
                    copy_rvec(x[ind0], xi);
                    for (j = i + 1; j < isize0 ; j++)
                    {
                        pbc_dx(&pbc, xi, x_i1[j], dx);
                        r2 = iprod(dx, dx);
                        if (r2 <= rmax2 )
                        {
                            r_dist = sqrt(r2);
                            b22 = beta_lab_2_t[i]*beta_lab_2_t[j]  ;
                            b11 = beta_lab_1_t[i]*beta_lab_1_t[j]  ;
                            b21 = beta_lab_2_t[i]*beta_lab_1_t[j]  ;
                            b12 = beta_lab_1_t[i]*beta_lab_2_t[j]  ;

                            if ( r_dist <= fade)
                            {
                                mod_f = b22*s_square[0] + b11*c_square[0] - (b12 + b21)*cs_sc[0] ;
                                for (qq = 0; qq < nbinq; qq++)
                                {
                                   temp_method[qq] += mod_f*cos(iprod(arr_qvec[qq],dx));
                                }
                            }
                            else
                            {
                                mod_f = (b22*s_square[0] + b11*c_square[0] - (b12 + b21)*cs_sc[0])*sqr(cos((r_dist-fade)*inv_width)) ;
                                for (qq = 0; qq < nbinq; qq++)
                                {
                                   /*fprintf(stderr,"s_method[g][qq] %f\n", mod_f*cos((minq+dq*qq)*iprod(arr_qvec,dx)) );*/
                                   temp_method[qq] += mod_f*cos(iprod(arr_qvec[qq],dx));
                                }
                            }
                        }
                    }
                }
                incoh_temp = (beta_lab_sq_1*c_square[0] +  beta_lab_sq_2*s_square[0] - 2.0*beta_lab_1_2*cs_sc[0])*invsize0;
                for (qq = 0; qq < nbinq; qq++)
                {
                   coh_temp = 2.0*temp_method[qq]*invsize0;
                   s_method_coh[g][qq] += coh_temp ;
                   s_method[g][qq] += coh_temp + incoh_temp ;
                   s_method_incoh[g][qq] += incoh_temp;
                }
            }
            nframes++;
        }
        while (read_next_x(oenv, status, &t, x, box));
    }

    else if (method[0] == 's' && bSpectrum == TRUE && fade == 0.0 && bThetaswipe == FALSE )
    {
        fprintf(stderr,"loop with sumexp method and spectrum \n");
        do
        {
            /* Must init pbc every step because of pressure coupling */
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
                snew(cos_t,nbinq);
                snew(sin_t,nbinq);
                mu_sq = 0.0;
                for (i = 0; i < isize0; i++)
                {
                    ind0  = mols->index[molindex[g][i]];
                    copy_rvec(x[ind0], xi);
                    pbc_dx(&pbc, x[ind0+1], xi, x01);
                    pbc_dx(&pbc, x[ind0+2], xi, x02);
                    induced_second_order_dipole(x01, x02, vec_polout, vec_polin, beta_mol, &mu_ind);
                    mu_sq += mu_ind*mu_ind; 
                    for (qq = 0; qq < nbinq; qq++)
                    {
                       q_xi = iprod(arr_qvec[qq],xi);
                       cos_t[qq] += mu_ind*cos(q_xi);
                       sin_t[qq] += mu_ind*sin(q_xi);
                    }
                }
                incoh_temp = (mu_sq)*invsize0;
                for (qq = 0; qq < nbinq; qq++)
                {
                   tot_temp = (cos_t[qq]*cos_t[qq] + sin_t[qq]*sin_t[qq])*invsize0;
                   s_method[g][qq] +=  tot_temp  ;
                   s_method_coh[g][qq] += tot_temp - incoh_temp;
                   s_method_incoh[g][qq] += incoh_temp ;
                }
            }
            nframes++;
        }
        while (read_next_x(oenv, status, &t, x, box));
    }
    else if (method[0] == 's' && bSpectrum == TRUE && fade == 0.0 && bThetaswipe == TRUE )
    {
        fprintf(stderr,"loop with sumexp method, spectrum and swiping over the scattering angle\n");
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
                       snew(mu_sq_t[rr][tt], nbintheta);
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
                    pbc_dx(&pbc, x[ind0+1], xi, x01);
                    pbc_dx(&pbc, x[ind0+2], xi,  x02);

                    if (bRandbeta)
                    {
                         beta_gaussian_noise(beta_mol_const, beta_const_dev, rng, &beta_mol);
//                         beta_mol[2][0][0] = bmzxx*gmx_rng_gaussian_real(rng) + bmzxx;
//                         beta_mol[2][1][1] = bmzyy*gmx_rng_gaussian_real(rng) + bmzyy;
//                         beta_mol[2][2][2] = bmzzz*gmx_rng_gaussian_real(rng) + bmzzz;
//                         beta_mol[0][2][0] = bmzxx*gmx_rng_gaussian_real(rng) + bmzxx;
//                         beta_mol[0][0][2] = beta_mol[0][2][0];
//                         beta_mol[1][2][1] = bmzyy*gmx_rng_gaussian_real(rng) + bmzyy;
//                         beta_mol[1][1][2] = beta_mol[1][2][1];
                         //printf("%f %f %f %f %f\n",beta_mol[2][0][0], beta_mol[2][1][1], beta_mol[2][2][2], beta_mol[0][2][0], beta_mol[1][2][1]);
                    }
                    for (rr = 0; rr < nfaces; rr++)
                    {
                        for (tt = 0; tt < nbintheta; tt++ )
                        {
                            for (c  = 0; c < nbingamma; c++)
                            {
                                induced_second_order_dipole(x01, x02, vec_pout_theta_gamma[rr][tt][c], vec_pin_theta_gamma[rr][tt][c], beta_mol, &mu_ind);
                                mu_sq_t[rr][tt][c] += mu_ind*mu_ind;
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

    else if ((fade == 0.0  && method[0]=='m') || (method[0] == 's' && fade != 0.0) || bSpectrum == FALSE  )
    {
        gmx_fatal(FARGS,"Combination of method %c, spectrum %d and fading %d  unavailable", method[0], bSpectrum, fade);
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
    if (bSpectrum == FALSE)
    {
       for (g = 0; g < ng; g++)
       {
           for (qq = 0; qq <= nbintheta ; qq++)
           {
              s_method[g][qq] = s_method[g][qq]/(nframes)  ;
              s_method_coh[g][qq] = s_method_coh[g][qq]/(nframes)  ;
              s_method_incoh[g][qq] = s_method_incoh[g][qq]/(nframes) ;
           }  
       }
    }
    else 
    {
       for (g = 0; g < ng; g++)
       {
           for (qq = 0; qq < nbinq ; qq++)
           {
              s_method[g][qq] = s_method[g][qq]/(nframes)  ;
              s_method_coh[g][qq] = s_method_coh[g][qq]/(nframes)  ;
              s_method_incoh[g][qq] = s_method_incoh[g][qq]/(nframes) ;
           }
       }
    }

    if (bThetaswipe == TRUE )
    {
       int      nplots = 1;
       sprintf(gtitle, "Non-linear optical scattering ");
       fp = xvgropen(fnTHETA, "S(theta)", "theta", "S(theta)", oenv);
       sprintf(refgt, "%s", "");
       fprintf(fp, "@    s%d legend \"incoherent\"\n",nplots);
       fprintf(fp, "@target G0.S%d\n",nplots);
       fprintf(fp, "@type xy\n");
       nplots++ ;
       for (tt = 0; tt < nbintheta ; tt++)
       {
          theta = 2.0*theta_vec[tt]*180.0/M_PI;
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
                theta = 2.0*theta_vec[tt]*180.0/M_PI;
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
                theta = 2.0*theta_vec[tt]*180.0/M_PI;
                fprintf(fp, "%10g", theta);
                fprintf(fp, " %10g", s_method_t[0][rr][tt][qq]/nframes*invgamma);
                fprintf(fp, "\n");
             }
             fprintf(fp,"&\n");
          }
       }
       gmx_ffclose(fp);

       nplots = 1;
       sprintf(gtitle, "Non-linear optical scattering ");
       fpn = xvgropen(fnSFACT, "S(q)", "q (nm^-1)", "S(q)", oenv);
       sprintf(refgt, "%s", "");
       for (tt = 0; tt < nbintheta ; tt++)
       {
          theta = 2.0*theta_vec[tt]*180.0/M_PI;
          if ((round(theta) == -45.0) || (round(theta) == -30.0 ) || (round(theta) == -60.0 ) || (round(theta) == -90.0)
             || (round(theta) == -150.0)  || (round(theta) == -120.0) || (tt ==  0))
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
       gmx_ffclose(fpn);
    }
    else
    {
        sprintf(gtitle, "Non-linear optical scattering ");
        fp = xvgropen(fnSFACT, gtitle, "q(nm^-1)", "S(q)", oenv);
        sprintf(refgt, "%s", "");
        fprintf(fp, "@    s0 legend \"coherent\" \n");
        fprintf(fp, "@target G0.S0\n");
        if (ng == 1)
        {
            if (output_env_get_print_xvgr_codes(oenv))
            {
                fprintf(fp, "@ subtitle \"%s%s - %s\"\n", grpname[0], refgt, grpname[1]);
            }
        }
        else
        {
            if (output_env_get_print_xvgr_codes(oenv))
            {
                fprintf(fp, "@ subtitle \"reference %s%s\"\n", grpname[0], refgt);
            }
            xvgr_legend(fp, ng, (const char**)(grpname+1), oenv);
        }
        for (qq = 0; qq < nbinq  ; qq++)
        {
            fprintf(fp, "%10g", norm(arr_qvec[qq]));
            for (g = 0; g < ng; g++)
            {
                fprintf(fp, " %10g", s_method_coh[g][qq]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp,"&\n");
        fprintf(fp, "@    s1 legend \"incoherent\"\n");
        fprintf(fp, "@target G0.S1\n");
        fprintf(fp, "@type xy\n");
        for (qq = 0; qq < nbinq  ; qq++)
        {
            fprintf(fp, "%10g", norm(arr_qvec[qq]));
            for (g = 0; g < ng; g++)
            {
                fprintf(fp, " %10g", s_method_incoh[g][qq]);
            }
            fprintf(fp, "\n");
        }
        fprintf(fp,"&\n");
        fprintf(fp, "@    s2 legend \"total\"\n");
        fprintf(fp, "@target G0.S2\n");
        fprintf(fp, "@type xy\n");
        for (qq = 0; qq < nbinq  ; qq++)
        {
            fprintf(fp, "%10g", norm(arr_qvec[qq]));
            for (g = 0; g < ng; g++)
            {
                fprintf(fp, " %10g", s_method[g][qq]);
            }
            fprintf(fp, "\n");
        }
        gmx_ffclose(fp);
    }
    do_view(oenv, fnSFACT, NULL);

    for (g = 0; g < ng; g++)
    {
       sfree(s_method[g]);
       sfree(s_method_coh[g]);
       sfree(s_method_incoh[g]);
    }
    sfree(beta_mol_1d);
    sfree(s_method);
    sfree(s_method_coh);
    sfree(arr_qvec);
    sfree(cos_t);
    sfree(sin_t);
    sfree(temp_method);
    sfree(cq);
    sfree(sq);
    sfree(c_square);
    sfree(s_square);
    sfree(cs_sc);
    sfree(beta_lab_1_t);
    sfree(beta_lab_2_t);
    
    gmx_rng_destroy(rng);

    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            sfree(beta_mol[i][j]);
            sfree(beta_mol_const[i][j]);
            sfree(beta_const_dev[i][j]);
        }
    }
    for (j = 0; j < DIM; j++)
    {
        sfree(beta_mol[j]);
        sfree(beta_mol_const[j]);
        sfree(beta_const_dev[j]);
    }
    sfree(beta_mol);
    sfree(beta_mol_const);
    sfree(beta_const_dev);
}


void rotate_beta_theta( const rvec xv2, const rvec xv3, const rvec pin1, const rvec pin2, real *beta_m, real *beta_2, real *beta_1)
{
    rvec xvec, yvec, zvec;
    real x_in1, y_in1, z_in1;
    real x_in2, y_in2, z_in2;

    rvec_add( xv2, xv3, zvec);
    cprod(zvec,xv2, yvec);
    unitv(zvec,zvec);
    unitv(yvec,yvec);
    cprod(yvec,zvec,xvec);

    x_in1 = iprod(xvec, pin1);
    y_in1 = iprod(yvec, pin1);
    z_in1 = iprod(zvec, pin1);
    x_in2 = iprod(xvec, pin2);
    y_in2 = iprod(yvec, pin2);
    z_in2 = iprod(zvec, pin2);

    /*should be faster than triple loop over i,j,k and uses only the 7 non zero component of water (C2v symmetry) */
    *beta_2 = beta_m[0]*zvec[2]*x_in1*x_in2 +
              beta_m[1]*zvec[2]*y_in1*y_in2 +
              beta_m[2]*zvec[2]*z_in1*z_in2 +
              beta_m[3]*xvec[2]*x_in1*z_in2 +
              beta_m[4]*xvec[2]*z_in1*x_in2 +
              beta_m[5]*yvec[2]*y_in1*z_in2 +
              beta_m[6]*yvec[2]*z_in1*y_in2 ;

    *beta_1 =  beta_m[0]*zvec[0]*x_in1*x_in2 +
              beta_m[1]*zvec[0]*y_in1*y_in2 +
              beta_m[2]*zvec[0]*z_in1*z_in2 +
              beta_m[3]*xvec[0]*x_in1*z_in2 +
              beta_m[4]*xvec[0]*z_in1*x_in2 +
              beta_m[5]*yvec[0]*y_in1*z_in2 +
              beta_m[6]*yvec[0]*z_in1*y_in2 ;
}

void rotate_wave_vec(const rvec wave_vec, const int rot_label, rvec rot_vec)
{
     real rot_angle = -M_PI/4.0;
     rvec rotated_vec;
     matrix rot_matrix;
     rvec rot_vec_tmp;

     if (rot_label == 0)
     {
        svmul(1.0, wave_vec, rot_vec);
     }
     if (rot_label >= 1)
     {
       //first rotate the diagonal of -45 degrees wrt to the y axis to get the vector parallel to the x axis
       rot_matrix[0][0] = cos(rot_angle);
       rot_matrix[0][2] = sin(rot_angle);
       rot_matrix[1][1] = 1.0;
       rot_matrix[2][0] = -sin(rot_angle);
       rot_matrix[2][2] = cos(rot_angle);
       mvmul(rot_matrix, wave_vec, rotated_vec);
       svmul(1/sqrt(2.0), rotated_vec, rot_vec);
       copy_rvec(rot_vec, rot_vec_tmp);
       clear_mat(rot_matrix);
       
     }
     if (rot_label >= 2)
     {
       //now rotate the vector of -45 degrees wrt to the z axis to get the vector parallel to the diagonal of the x and y axis
       rot_matrix[0][0] = cos(-rot_angle);
       rot_matrix[0][1] = -sin(-rot_angle);
       rot_matrix[1][0] = sin(-rot_angle);
       rot_matrix[1][1] = cos(-rot_angle);
       rot_matrix[2][2] = 1.0;
       mvmul(rot_matrix, rot_vec_tmp, rotated_vec);
       svmul(sqrt(2.0), rotated_vec, rot_vec);
       copy_rvec(rot_vec, rot_vec_tmp);
       clear_mat(rot_matrix);      
     }
     if (rot_label >= 3)
     {
       // now rotate the vector of -45 degrees wrt the z axis to get the vector parallel to the y axis
       rot_matrix[0][0] = cos(-rot_angle);
       rot_matrix[0][1] = -sin(-rot_angle);
       rot_matrix[1][0] = sin(-rot_angle);
       rot_matrix[1][1] = cos(-rot_angle);
       rot_matrix[2][2] = 1.0;
       mvmul(rot_matrix, rot_vec_tmp, rotated_vec);
       svmul(1.0/(sqrt(2.0)), rotated_vec, rot_vec);
       copy_rvec(rot_vec, rot_vec_tmp);
       clear_mat(rot_matrix);
     }
     if (rot_label >= 4)
     {
       // now rotate the vector of 45 degrees wrt the x axis to get the diagonal of the y and z axes
       rot_matrix[0][0] = 1.0;
       rot_matrix[1][1] = cos(rot_angle);
       rot_matrix[1][2] = -sin(rot_angle);
       rot_matrix[2][1] = sin(rot_angle);
       rot_matrix[2][2] = cos(rot_angle);
       mvmul(rot_matrix, rot_vec_tmp, rotated_vec);
       svmul(sqrt(2.0), rotated_vec, rot_vec);
       copy_rvec(rot_vec, rot_vec_tmp);
       clear_mat(rot_matrix);
     }
     if (rot_label >= 5)
     {
       //now rotate the vector of 45 degrees wrt the x axis to get the vector parallel to the z axis
       rot_matrix[0][0] = 1.0;
       rot_matrix[1][1] = cos(rot_angle);
       rot_matrix[1][2] = -sin(rot_angle);
       rot_matrix[2][1] = sin(rot_angle);
       rot_matrix[2][2] = cos(rot_angle);
       mvmul(rot_matrix, rot_vec_tmp, rotated_vec);
       svmul(1.0/(sqrt(2.0)), rotated_vec, rot_vec);
       copy_rvec(rot_vec, rot_vec_tmp);
       clear_mat(rot_matrix);
     }

}


void beta_gaussian_noise(real ***beta_const_mean, real ***beta_const_dev, gmx_rng_t rng, real ****beta_gauss)
{
    int a, b, c;

    for (a = 0; a < DIM; a++)
    {
        for (b = 0; b < DIM; b++)
        {
            for (c = 0; c < DIM; c++)
            {
                (*beta_gauss)[a][b][c] = beta_const_dev[a][b][c]*gmx_rng_gaussian_real(rng) + beta_const_mean[a][b][c];
                //printf("beta_mol %d %d %d %f\n", a, b, c, (*beta_gauss)[a][b][c]);
            }
        }
    }  


}

void induced_second_order_dipole(const rvec xv2, const rvec xv3, const rvec pout, const rvec pin, real ***betamol, real *mu_ind)
{
    int i, j, k;
    int p, l, m;
    rvec xvec, yvec, zvec;
    //betal is the hyperpolarizability in lab frame
    real ***betal, **cosdirmat, mu_temp = 0.0;

    //initialize beta
    snew(betal, DIM);
    snew(cosdirmat, DIM);
    for (i = 0; i < DIM; i++)
    {
        snew(betal[i], DIM);
        snew(cosdirmat[i], DIM);
    }
    
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            snew(betal[i][j], DIM);
        }
    }

    rvec_add( xv2, xv3, zvec);
    cprod(zvec,xv2, yvec);
    unitv(zvec,zvec);
    unitv(yvec,yvec);
    cprod(yvec,zvec,xvec);


    cosdirmat[0][0] = xvec[0]; cosdirmat[0][1] = yvec[0]; cosdirmat[0][2] = zvec[0];
    cosdirmat[1][0] = xvec[1]; cosdirmat[1][1] = yvec[1]; cosdirmat[1][2] = zvec[1];
    cosdirmat[2][0] = xvec[2]; cosdirmat[2][1] = yvec[2]; cosdirmat[2][2] = zvec[2];

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

    //FREE beta lab
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            sfree(betal[i][j]);
        }
    }
    
    for (j = 0; j < DIM; j++)
    {
        sfree(betal[j]);
        sfree(cosdirmat[j]);
    }
    
    sfree(betal);
    sfree(cosdirmat);
}

void dipole_atom2mol(int *n, int *index, t_block *mols)
{
    int nmol, i, j, m;

    nmol = 0;
    i    = 0;
    while (i < *n)
    {
        m = 0;
        while (m < mols->nr && index[i] != mols->index[m])
        {
            m++;
        }
        if (m == mols->nr)
        {
            gmx_fatal(FARGS, "index[%d]=%d does not correspond to the first atom of a molecule", i+1, index[i]+1);
        }
        for (j = mols->index[m]; j < mols->index[m+1]; j++)
        {
            if (i >= *n || index[i] != j)
            {
                gmx_fatal(FARGS, "The index group is not a set of whole molecules");
            }
            i++;
        }
        index[nmol++] = m;
        //fprintf(stderr,"index[nmol]: index[%d]=%d\n",nmol, index[nmol]);
    }
    printf("There are %d molecules in the selection\n", nmol);
    *n = nmol;
}


int gmx_nonlinearopticalscatteringtheta(int argc, char *argv[])
{
    const char        *desc[] = {
        "The structure of liquids can be studied by elastic second harmonic scattering.",
        "[THISMODULE] calculates the non-linear optical scattering intensity per molecule in 2 different ways.",
        "The simplest method (sumexp) is to use 1/N<|sum_i beta_IJK(i) exp[iq dot r_i]|^2>.",
        "This however converges slowly with the simulation time and is more prone to noise.[PAR]",
        "A method that converges more quickly wrt simulation time (modsumexp, default), but nmol times more expensive is to do",
        "I(q)=1/N<sum_i beta_IJK(i)^2> + 1/N< sum_i sum_{i!=j} beta_IJK(i)beta_IJK(j) cos(q dot r_ij) >.",
        "I(q)=incoherent term + coherent term. For more details see Bersohn, et al. JCP 45, 3184 (1966)",
        "The values for the hyperpolarizability for a water molecule beta_IJK are taken by",
        "A. V. Gubskaya et al., Mol. Phys. 99, 13 1107 (2001) (computed at MP4 level with mean field liquid water).",
        "pout, pin1, pin2 are the polarization directions of the three beams.",
        "Common polarization combinations are PSS, PPP, SPP, SSS .",
        "Under Kleinmann symmetry (default) beta_ijj = beta_jij = beta_jji otherwise beta_ijj = beta_jij. [PAR]",
    };
    static gmx_bool    bPBC = TRUE, bKleinmannsymm = TRUE, bSpectrum = TRUE, bThetaswipe = FALSE, bRandbeta = FALSE;
    static real        fade = 0.0;
    static real        koutx = 1.0, kouty = 0.0 , koutz = 0.0, kinx = 0.0, kiny = 0.0, kinz = 1.0, pout_angle = 0.0 , pin_angle = 0.0;
    static real        bmzxx = 5.7 , bmzyy = 10.9 , bmzzz = 31.6 ;
    static int         ngroups = 1, nbintheta = 10, nbingamma = 2 ,qbin = 1, nbinq = 10 ;

    static const char *methodt[] = { NULL, "modsumexp", "sumexp" ,NULL }; 

    t_pargs            pa[] = {
        { "-nbintheta",      FALSE, etINT, {&nbintheta},
        "number of bins over scattering angle theta chosen between -pi/2 and + pi/2 (available only with thetaswipe)" },
        { "-nplanes",      FALSE, etINT, {&nbingamma},
        "number of scattering planes that lie on the scattered wave-vector to average over, -PI/2< gamma< PI/2" },
        { "-qbin",      FALSE, etINT, {&qbin},
        "which wave-vector to sample which is 2pi/box-length*qbin" },
        { "-nbinq",      FALSE, etINT, {&nbinq},
        "how many bins in the reciprocal space" },
        { "-bmzxx",         FALSE, etREAL, {&bmzxx}, "component of beta in zxx " },
        { "-bmzyy",         FALSE, etREAL, {&bmzyy}, "component of beta in zyy " },
        { "-bmzzz",         FALSE, etREAL, {&bmzzz}, "component of beta in zzz " },
        { "-koutx",         FALSE, etREAL, {&koutx}, "direction of kout in x direction " },
        { "-kouty",         FALSE, etREAL, {&kouty}, "direction of kout in y direction " },
        { "-koutz",         FALSE, etREAL, {&koutz}, "direction of kout in z direction " },
        { "-kinx",         FALSE, etREAL, {&kinx}, "direction of kin in x direction " },
        { "-kiny",         FALSE, etREAL, {&kiny}, "direction of kin in y direction " },
        { "-kinz",         FALSE, etREAL, {&kinz}, "direction of kin in z direction " },
        { "-pout",         FALSE, etREAL, {&pout_angle}, "polarization angle of outcoming beam in degrees. For P choose 0, for S choose 90" },
        { "-pin",         FALSE, etREAL, {&pin_angle}, "polarization angle of incoming beam in degrees. For P choose 0, for S choose 90" },
        { "-method",     FALSE, etENUM, {methodt},
          "I(q) using the different methods" },
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances. Without PBC the maximum range will be three times the largest box edge." },
        { "-spectrum",    FALSE, etBOOL, {&bSpectrum}, "compute spectrum at fixed angle of 5 degrees (only sumexp)" },
        { "-randbeta",    FALSE, etBOOL, {&bRandbeta}, "for each molecule and at each step add gaussian noise to the value of beta determined by bmzxx, bmzyy, bmzzz" },        
        { "-thetaswipe",    FALSE, etBOOL, {&bThetaswipe}, "compute spectrum swiping over the scattering angle theta" },
        { "-klein",    FALSE, etBOOL, {&bKleinmannsymm}, "Assume Kleinmann symmetry" },
        { "-ng",       FALSE, etINT, {&ngroups}, 
          "Number of secondary groups to compute RDFs around a central group" },
        { "-fade",     FALSE, etREAL, {&fade},
          "In the method the modification function cos((rij-fade)*pi/(2*(L/2-fade))) is used in the fourier transform."
          " If fade is 0.0 nothing is done." },

    };
#define NPA asize(pa)
    const char        *fnTPS, *fnNDX;
    output_env_t       oenv;
    int           *gnx;
    int            nFF[2];
    atom_id      **grpindex;
    char         **grpname = NULL;
    /*gmx_bool       bGkr, bMU, bSlab;*/

    t_filenm           fnm[] = {
        { efTRX, "-f",  NULL,     ffREAD },
        { efTPS, NULL,  NULL,     ffREAD },
        { efNDX, NULL,  NULL,     ffOPTRD },
        { efXVG, "-o",  "non_linear_sfact",    ffWRITE },
        { efXVG, "-otheta", "non_linear_sfact_vs_theta", ffOPTWR },

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

    if (!fnTPS && !fnNDX)
    {
        gmx_fatal(FARGS, "Neither index file nor topology file specified\n"
                  "Nothing to do!");
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
   
    do_nonlinearopticalscatteringtheta(top, ftp2fn(efTRX, NFILE, fnm),
           opt2fn("-o", NFILE, fnm), opt2fn("-otheta", NFILE, fnm), methodt[0],  bPBC, bKleinmannsymm, bSpectrum, bRandbeta, bThetaswipe, qbin, nbinq,
           koutx, kouty, koutz, kinx, kiny, kinz  , bmzxx, bmzyy, bmzzz ,nbintheta, nbingamma, pin_angle, pout_angle, 
           fade, gnx, grpindex, grpname, ngroups, oenv);

    return 0;
}

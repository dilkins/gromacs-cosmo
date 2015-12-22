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

#define nm2Bohr 10.0/0.529177
#define AU2VNM 5.14220652e2


static void do_eshs(t_topology *top,  const char *fnTRX, const char *fnMAP, const char *fnKRR,
                   const char *fnGRD, const char *fnPOT,
                   const char *fnSFACT, const char *fnTHETA, const char *method,
                   gmx_bool bIONS, char *catname, char *anname, gmx_bool bPBC, 
                   int qbin, int nbinq, real koutx, real kouty, real koutz,
                   real kinx, real kiny, real kinz,
                   int nbintheta, int nbingamma, real pin_angle, real pout_angle,
                   real cutoff_field, real kernstd,
                   int *isize, int  *molindex[], char **grpname, int ng,
                   const output_env_t oenv)
{
    FILE          *fp;
    FILE          *fpn;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, nanions, ncations, i, j, k, qq, n, c, tt, rr, nframes, nfaces;
    real         **s_method, **s_method_coh, **s_method_incoh, *temp_method, ****s_method_t, ****s_method_coh_t, ****s_method_incoh_t, ***mu_sq_t ;
    real           qnorm, maxq, coh_temp = 0.0,  incoh_temp = 0.0, tot_temp = 0.0, gamma = 0.0 ,theta0 = 5.0, check_pol;
    real          *cos_t, *sin_t, ****cos_tq, ****sin_tq,   mu_ind =0.0, mu_sq =0.0 ;
    real          *field_ad, electrostatic_cutoff, ***beta_mol;
    int            max_i, isize0, ind0;
    real           t, rmax2, rmax,  r, r_dist, r2, q_xi, dq;
    real           segvol, spherevol, prev_spherevol, invsize0, invgamma, theta=0, *theta_vec;
    rvec          *x, dx,  *x_i1, xi, x01, x02, *arr_qvec, **arr_qvec_faces ,vec_polin, vec_polout, ***vec_pout_theta_gamma, ***vec_pin_theta_gamma;
    rvec           pol_perp, pol_par,vec_kout, vec_2kin, pol_in1, pol_in2, vec_kout_2kin ;
    rvec           xvec, yvec, zvec;
    matrix         cosdirmat; 
    real           invvol, invvol_sum;
    t_Map         *Map=NULL;
    t_Kern        *Krr = NULL;
    t_Ion         *Cation=NULL,*Anion=NULL;
    matrix         box, box_pbc;
    int            ePBC = -1, ePBCrdf = -1;
    int            nplots = 1; 
    t_block       *mols = NULL;
    t_atom        *atom = NULL;
    t_pbc          pbc;
    gmx_rmpbc_t    gpbc = NULL;
    gmx_rng_t      rng = NULL;
    int            mol, a;

    // read electrostatic fit map input file
    if (fnMAP)
    {
       Map=(t_Map *)calloc(1,sizeof(t_Map));
       readMap(fnMAP, Map);
    }
    else
    {
       Krr = (t_Kern *)calloc(1,sizeof(t_Kern));
       readKern(fnKRR, fnGRD, fnPOT, Krr);
       Krr->kerndev = -0.5/((kernstd*kernstd));
       fprintf(stderr,"kernel of standard dev %f\n", Krr->kerndev);
       fprintf(stderr,"first 2 and last coefficients %f %f %f\n",Krr->coeff[0][XX][XX][XX], Krr->coeff[0][ZZ][ZZ][ZZ], Krr->coeff[Krr->ndataset-1][ZZ][ZZ][ZZ]);
       fprintf(stderr,"last potential %f\n",Krr->V[Krr->gridpoints-1][Krr->ndataset-1]);
    }

    atom = top->atoms.atom;
    mols = &(top->mols);
    isize0 = isize[0];
    nfaces = 6;
    invsize0 = 1.0/isize0;
    invgamma = 1.0/nbingamma;
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);

    fprintf(stderr,"\nNumber of atoms %d\n",natoms);
    fprintf(stderr,"\nName of molecules %s\n",grpname[0]);

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
        electrostatic_cutoff =  min(cutoff_field*cutoff_field, rmax2);
        fprintf(stderr, "rmax2 = %f\n", rmax2);
        fprintf(stderr, "cutoff for electric field calculation = %f\n", sqrt(electrostatic_cutoff));
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
        /* this is THE array */
        /*allocate memory for s_method array */
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
                    pbc_dx(&pbc, x[ind0+2], xi, x02);
 
                    rvec_add( x01, x02, zvec);
                    cprod(zvec,x01, yvec);
                    unitv(yvec,yvec);
                    unitv(zvec,zvec);
                    cprod(yvec,zvec,xvec);
                    cosdirmat[0][0] = xvec[0]; cosdirmat[0][1] = yvec[0]; cosdirmat[0][2] = zvec[0];
                    cosdirmat[1][0] = xvec[1]; cosdirmat[1][1] = yvec[1]; cosdirmat[1][2] = zvec[1];
                    cosdirmat[2][0] = xvec[2]; cosdirmat[2][1] = yvec[2]; cosdirmat[2][2] = zvec[2];                    
                    if (fnMAP)
                    {
                        calc_efield_map(&pbc, top, mols, Cation, Anion, molindex, g, isize0, ncations, nanions,
                                x, ind0, xvec, yvec, zvec, electrostatic_cutoff, &field_ad);
                        calc_beta_efield_map(Map, field_ad, beta_mol, &beta_mol );
                    }
                    else
                    {
                        calc_beta_krr(Krr, &pbc, top, mols, molindex, g, isize0, x, ind0, cosdirmat  ,electrostatic_cutoff, beta_mol, &beta_mol );
                        printf("beta %d %d %d %f\n",0, 0, 0, beta_mol[0][0][0]);
                        printf("beta %d %d %d %f\n",2, 2, 2, beta_mol[2][2][2]);
                        printf("beta %d %d %d %f\n",2, 0, 0, beta_mol[2][0][0]);
                        printf("beta %d %d %d %f\n",0, 2, 0, beta_mol[0][2][0]);
                        printf("beta %d %d %d %f\n",0, 0, 2, beta_mol[0][0][2]);
                    }
                    for (rr = 0; rr < nfaces; rr++)
                    {
                        for (tt = 0; tt < nbintheta; tt++ )
                        {
                            for (c  = 0; c < nbingamma; c++)
                            {

                                induced_second_order_fluct_dipole(cosdirmat, 
                                                                  vec_pout_theta_gamma[rr][tt][c], vec_pin_theta_gamma[rr][tt][c], 
                                                                  beta_mol, &mu_ind);
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

    else if ( method[0]!='s')
    {
        gmx_fatal(FARGS,"Method %c\n", method[0]);
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
    do_view(oenv, fnSFACT, NULL);

    for (g = 0; g < ng; g++)
    {
       sfree(s_method[g]);
       sfree(s_method_coh[g]);
       sfree(s_method_incoh[g]);
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
    }
    for (j = 0; j < DIM; j++)
    {
        sfree(beta_mol[j]);
    }
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
    for (i = 0; i<81; i++)
    {
       n_outputs = fscanf(fm,"%f ",&Map->gamma_gas[i]);

    }
    // read coefficients
    for (i=0;i<27;i++) {  // 27 elements
      for (j=0;j<212;j++) {  // 212 components
        n_outputs = fscanf(fm,"%f ",&Map->D[i][j]);
      }
    }
  
    fclose(fm);
}

void readKern(const char *fnKRR, const char *fnGRD, const char *fnPOT, t_Kern *Krr)
{
    FILE *fk, *fg, *fp;
    int i, j, ch = 0, n_outputs ;
    int a, b, c;

    fp = gmx_ffopen(fnPOT, "r");
    fk = gmx_ffopen(fnKRR, "r");
    fg = gmx_ffopen(fnGRD, "r");

    while(!feof(fp))
    {
      ch = fgetc(fp);
      if(ch == '\n')
      {
         Krr->ndataset ++;
      }
    }
    Krr->ndataset ++;
    rewind(fp);

    while(!feof(fg))
    {
      ch = fgetc(fg);
      if(ch == '\n')
      {
         Krr->gridpoints ++;
      }
    }
    rewind(fg);
    Krr->gridpoints ++;

    //allocate all pointers needed for krr operations
    snew(Krr->V, Krr->ndataset);
    snew(Krr->Vgrid, Krr->gridpoints);
    snew(Krr->translgrid, Krr->gridpoints);
    snew(Krr->coeff, Krr->ndataset);
    //allocate potential and krr coefficients
    for ( i = 0; i< Krr->ndataset; i++)
    {
        snew(Krr->V[i], Krr->gridpoints);
        snew(Krr->coeff[i],DIM);
    }
    
    //allocate coefficients for dimensions of beta;
    //fprintf(stderr,"allocated potential\n");
    for ( i = 0; i< Krr->ndataset; i++)
    {
       for ( a = 0; a < DIM; a++)
       {
           snew(Krr->coeff[i][a],DIM);
       }
    }

    for ( i = 0; i< Krr->ndataset; i++)
    {
       for ( a = 0; a < DIM; a++)
       {
           for( b = 0; b < DIM; b++)
           { 
               snew(Krr->coeff[i][a][b],DIM);
           }
       }
    }
    //fprintf(stderr,"allocated coefficients\n");

    for (i = 0; i < Krr->ndataset; i++)
    {
        for ( j = 0; j < Krr->gridpoints; j++)
        {
            n_outputs = fscanf(fp, "%f", &Krr->V[i][j]) ;
        }
        for (a = 0; a < DIM; a++)
        {
            for (b = 0; b < DIM; b++)
            {
                for (c = 0; c < DIM; c++)
                {
                     n_outputs = fscanf(fk, "%f", &Krr->coeff[i][a][b][c]);
                     //fprintf(stderr,"coefficient %d %d %d %d %f\n",i,a,b,c, Krr->coeff[i][a][b][c]);
                }
            }
        }
    }
    for ( j = 0; j < Krr->gridpoints; j++)
    {
        n_outputs = fscanf(fg, "%f", &Krr->Vgrid[j][XX]);
        n_outputs = fscanf(fg, "%f", &Krr->Vgrid[j][YY]);
        n_outputs = fscanf(fg, "%f", &Krr->Vgrid[j][ZZ]);
        if (Krr->Vgrid[j][XX] == 0.0 && Krr->Vgrid[j][YY] == 0.0 &&  Krr->Vgrid[j][ZZ] == 0.0)
        {
           if (Krr->V[0][j] == 0.0)
           {
              Krr->gridcenter = j;
           }
           else
           {
              gmx_fatal(FARGS,"could not find the correct grid center (0.0 ,0.0 ,0.0) for the reference of the electrostatic potential\n check the input files");
           }
        }
    }
    fclose(fk);
    fclose(fp);
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
    }

    for (i = 0; i < DIM; i++)
    {
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
    }

    sfree(betal);
}

void calc_beta_krr(t_Kern *Krr, t_pbc *pbc,t_topology *top, t_block *mols, int  *molindex[],
                 const int gind , const int isize0,
                 rvec *x, const int imol, matrix cosdirmat, 
                 real electrostatic_cutoff , real ***betamol, real ****betamol_addr)
{
   int    grid_ind, set_ind, j, m, indj, a, b, c;
   real   d2, r, delV, kernel_fn, *Vij;
   rvec   vecij, vecrotated, vecr;
   matrix cosdirmatj;

   //initialize molecular beta
   for (a = 0; a < DIM; a++)
   {
       for (b = 0; b < DIM; b++)
       {
          for (c = 0; c < DIM; c++)
          {
             betamol[a][b][c] = 0.0 ;
          }
       }
   }
   //initialize potential at a point in the grid in the vicinity of molecule imol due to all molecules in a cutoff
   snew(Vij,Krr->gridpoints);
   for (j = 0; j < isize0; j++)
   {
       indj = mols->index[molindex[gind][j]];
       pbc_dx(pbc, x[imol], x[indj] ,vecij);
       d2 = norm2(vecij);
       if (d2 < electrostatic_cutoff && indj != imol) 
       {
           for (grid_ind = 0; grid_ind < Krr->gridpoints; grid_ind++)
           {
              //translate each point of the grid to the point of molecule i (not used)
              rvec_add(x[imol], Krr->Vgrid[grid_ind], Krr->translgrid[grid_ind]);
              for (m = 1; m < 4; m++)
              {
                  //compute displacement vector between molecule i and each atom of molecule j
                  //printf("vec i %f %f %f\n", x[imol][XX], x[imol][YY], x[imol][ZZ]);
                  pbc_dx(pbc,x[indj+m],x[imol],vecij);
                  //printf("vec j+m %f %f %f j+m %d\n", x[indj+m][XX], x[indj+m][YY], x[indj+m][ZZ],indj+m);
                  //printf("displecement vecij %f %f %f\n",vecij[XX], vecij[YY], vecij[ZZ]);
                  //rotate molecule j according to the direction cosines of molecule i
                  mvmul(cosdirmat,vecij,vecrotated);
                  //printf("rotated displacement vecij %f %f %f\n",vecrotated[XX], vecrotated[YY], vecrotated[ZZ]);
                  //now compute displacement vector between a point in the grid and the atom of molecule j rotated according to molecule i
                  rvec_sub(vecrotated,Krr->Vgrid[grid_ind],vecr);
                  //printf("translated grid %f %f %f grid_ind %d\n",Krr->translgrid[grid_ind][XX], Krr->translgrid[grid_ind][YY], Krr->translgrid[grid_ind][ZZ], grid_ind);
                  //printf("final displacement vector %f %f %f\n",vecr[XX],vecr[YY],vecr[ZZ]);
                  r=10.0*norm(vecr);
                  Vij[grid_ind] += (top->atoms.atom[m].q)/r;
              }
           }
       }
   }
   for (grid_ind = 0; grid_ind < Krr->gridpoints; grid_ind++)
   {
       Vij[grid_ind] -=  Vij[Krr->gridcenter];
       //printf("Vij %f grid label %d\n",Vij[grid_ind], grid_ind);
   }
   for (set_ind = 0; set_ind < Krr->ndataset; set_ind++)
   {
        delV = 0.0;
        for (grid_ind = 0; grid_ind < Krr->gridpoints; grid_ind++)
        {
            delV += sqr( (Vij[grid_ind]) - (Krr->V[set_ind][grid_ind]) );
        }
        kernel_fn = exp(delV*(Krr->kerndev));
        for (a = 0; a < DIM; a++)
        {
            for (b = 0; b < DIM; b++)
            {
                for (c = 0; c < DIM; c++)
                {
                   betamol[a][b][c] += (Krr->coeff[set_ind][a][b][c])*kernel_fn;
                   //fprintf(stderr,"beta %d %d %d  %f\n", a, b, c, betamol[a][b][c]);
                   //*(betamol_addr)[a][b][c] += (Krr->coeff[set_ind][a][b][c])*kernel_fn; 
                   //fprintf(stderr,"beta %d %d %d  %f\n", a, b, c, *(betamol_addr)[a][b][c]);
                }
            }
        }
   }
   sfree(Vij);
  *betamol_addr = betamol;
}

void calc_efield_map(t_pbc *pbc,t_topology *top, t_block *mols, t_Ion *Cation, t_Ion *Anion, int  *molindex[], 
                 const int gind , const int isize0, const int ncations, const int nanions,
                 rvec *x, const int imol,  rvec xvec, rvec yvec,  rvec zvec, real electrostatic_cutoff, real **field_addr)
{
   int    j, m, indj;
   real  chg, d2, r, r2, ir, ir3, ir5;
   real   vx, vy, vz, chr3, chr5;
   real  *Efield;
   rvec   vecij, vecr;

   snew(Efield,9);
   for (j = 0; j < isize0; j++)
   {
       indj = mols->index[molindex[gind][j]];
       pbc_dx(pbc, x[imol],x[indj],vecij);
       d2 = norm2(vecij);
       if (d2 < electrostatic_cutoff && indj != imol)
       {
          for (m = 1; m < 4; m++)
          {
              pbc_dx(pbc,x[imol],x[indj+m],vecr);
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
              Efield[0] += vx*chr3 ;
              Efield[1] += vy*chr3 ;
              Efield[2] += vz*chr3 ;
              // electric field gradient
              Efield[3] += chr5*(r2 - 3*vx*vx );
              Efield[4] += chr5*(r2 - 3*vy*vy );
              Efield[5] += chr5*(r2 - 3*vz*vz );
              Efield[6] -= chr5*( 3*vx*vy );
              Efield[7] -= chr5*( 3*vy*vz );
              Efield[8] -= chr5*( 3*vz*vx );
          }
       }
   }
   for (j = 0; j < ncations; j++)
   {
       pbc_dx(pbc, x[imol],x[Cation[j].atom[0]],vecij);
       d2 = norm2(vecij);
       if (d2 < electrostatic_cutoff )
       {
          for (m = 0; m < 1; m++)
          {
              pbc_dx(pbc,x[imol],x[Cation[j].atom[m]],vecr);
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
              Efield[0] += vx*chr3 ;
              Efield[1] += vy*chr3 ;
              Efield[2] += vz*chr3 ;
              // electric field gradient
              Efield[3] += chr5*(r2 - 3*vx*vx );
              Efield[4] += chr5*(r2 - 3*vy*vy );
              Efield[5] += chr5*(r2 - 3*vz*vz );
              Efield[6] -= chr5*( 3*vx*vy );
              Efield[7] -= chr5*( 3*vy*vz );
              Efield[8] -= chr5*( 3*vz*vx );
          }
       }
   }
   for (j = 0; j < nanions; j++)
   {
       pbc_dx(pbc, x[imol],x[Anion[j].atom[0]],vecij);
       d2 = norm2(vecij);
       if (d2 < electrostatic_cutoff)
       {
          for (m = 0; m < 1; m++)
          {
              pbc_dx(pbc,x[imol],x[Anion[j].atom[m]],vecr);
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
              Efield[0] += vx*chr3 ;
              Efield[1] += vy*chr3 ;
              Efield[2] += vz*chr3 ;
              // electric field gradient
              Efield[3] += chr5*(r2 - 3*vx*vx );
              Efield[4] += chr5*(r2 - 3*vy*vy );
              Efield[5] += chr5*(r2 - 3*vz*vz );
              Efield[6] -= chr5*( 3*vx*vy );
              Efield[7] -= chr5*( 3*vy*vz );
              Efield[8] -= chr5*( 3*vz*vx );
          }
       }
   }
   *field_addr = Efield;
    //fprintf(stderr,"efield %f %f %f\n",Efield[0], Efield[1], Efield[2]);
}

void calc_beta_efield_map(t_Map *Map, real *Efield, real ***betamol, real ****betamol_addr)
{

  int i,index,N=1,M=4;
  int a,b,c,d,e,f,g;
  real delta,term;
  real betatemp[27] = {0.0};

  for (i=0;i<27;i++)
  {

    delta=0;
    index=0;

   // Electric field components
    for (a=0;a<=N;a++) {
      for (b=0;b<=N;b++) {
        for (c=0;c<=N;c++) {
          if (a+b+c>0 && a+b+c<=N) {
            term=pow(Efield[0],a)*pow(Efield[1],b)*pow(Efield[2],c);
            delta+=term*Map->D[i][index];
            index++;
          }
        }
      }
    }
    // Electric field gradient components
    for (a=0;a<=M;a++) {
      for (b=0;b<=M;b++) {
        for (c=0;c<=M;c++) {
          for (d=0;d<=M;d++) {
            for (e=0;e<=M;e++) {
              for (f=0;f<=M;f++) {
                if (a+b+c+d+e+f>0 && a+b+c+d+e+f<=M) {
                  term=pow(Efield[3],a)*pow(Efield[4],b)*pow(Efield[5],c)*pow(Efield[6],d)*pow(Efield[7],e)*pow(Efield[8],f);
                  delta+=term*Map->D[i][index];
                  index++;
                }
              }
            }
          }
        }
      }
    }
    // Calculate molecular beta sorrounded by electric field of molecules
    betatemp[i]=Map->beta_gas[i]+
                  (Map->gamma_gas[i*3]*Efield[0]+Map->gamma_gas[i*3+1]*Efield[1]+Map->gamma_gas[i*3+2]*Efield[2])+delta;
//    printf("beta temp %f map beta gas %f gamma term %f delta %f Ex %f\n",betatemp[i], Map->beta_gas[i], (Map->gamma_gas[i*3]*Efield[0]+Map->gamma_gas[i*3+1]*Efield[1]+Map->gamma_gas[i*3+2]*Efield[2]), delta, Efield[0] );
  }
  for (a = 0; a < DIM; a++)
  {
      for (b = 0; b < DIM; b++)
      {
          for (c = 0; c < DIM; c++)
          {
              //remove this if when you want to compute the full tensor 
              //if ((a == 0 &&  b == 2 && c == 0) || 
              //    (a == 0 &&  b == 0 && c == 2) ||
              //    (a == 1 &&  b == 1 && c == 2) ||
              //    (a == 1 &&  b == 2 && c == 1) ||
              //    (a == 2 &&  b == 1 && c == 1) ||
              //   (a == 2 &&  b == 0 && c == 0) ||
              //    (a == 2 &&  b == 2 && c == 2)) 
              //{
                  betamol[a][b][c] = betatemp[a*9+ b*3+ c];
              //}
          }
      }
  }
  //printf("bzxx %f bzyy %f bzzz %f\n",betamol[2][0][0], betamol[2][1][1], betamol[2][2][2]);
  //printf("betatemp_zxx %f betatemp_zyy %f betatemp_zzz %f\n",betatemp[2*9], betatemp[2*9 + 1*3 +1], betatemp[26]);
  *betamol_addr = betamol;
   sfree(Efield);
}

    
 
int gmx_eshs(int argc, char *argv[])
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
    static gmx_bool          bPBC = TRUE, bIONS = FALSE;
    static real              electrostatic_cutoff = 2.0, kernstd = 0.9 ;
    static real              koutx = 1.0, kouty = 0.0 , koutz = 0.0, kinx = 0.0, kiny = 0.0, kinz = 1.0, pout_angle = 0.0 , pin_angle = 0.0;
    static int               ngroups = 1, nbintheta = 10, nbingamma = 2 ,qbin = 1, nbinq = 10 ;

    static const char *methodt[] = {NULL, "sumexp", NULL }; 
    static char *catname = NULL;
    static char *anname =  NULL;

    t_pargs            pa[] = {
        { "-nbintheta",     FALSE, etINT, {&nbintheta},
        "number of bins over scattering angle theta chosen between -pi/2 and + pi/2 (available only with thetaswipe)" },
        { "-nplanes",       FALSE, etINT, {&nbingamma},
        "number of scattering planes that lie on the scattered wave-vector to average over, -PI/2< gamma< PI/2" },
        { "-qbin",          FALSE, etINT, {&qbin},
        "which wave-vector to sample which is 2pi/box-length*qbin" },
        { "-nbinq",         FALSE, etINT, {&nbinq},
        "how many bins in the reciprocal space" },
        { "-koutx",         FALSE, etREAL, {&koutx}, "direction of kout in x direction " },
        { "-kouty",         FALSE, etREAL, {&kouty}, "direction of kout in y direction " },
        { "-koutz",         FALSE, etREAL, {&koutz}, "direction of kout in z direction " },
        { "-kinx",          FALSE, etREAL, {&kinx}, "direction of kin in x direction " },
        { "-kiny",          FALSE, etREAL, {&kiny}, "direction of kin in y direction " },
        { "-kinz",          FALSE, etREAL, {&kinz}, "direction of kin in z direction " },
        { "-pout",          FALSE, etREAL, {&pout_angle}, "polarization angle of outcoming beam in degrees. For P choose 0, for S choose 90" },
        { "-pin",           FALSE, etREAL, {&pin_angle}, "polarization angle of incoming beam in degrees. For P choose 0, for S choose 90" },
        { "-cutoff",        FALSE, etREAL, {&electrostatic_cutoff}, "cutoff for the calculation of electric field or electrostatic potential around a molecule" },
        { "-kernstd",       FALSE, etREAL, {&kernstd}, "standard deviation of kernel function. only makes sense if kernel ridge reduction is used"},
        { "-method",     FALSE, etENUM, {methodt}, "I(q) using the different methods" },
        { "-ions",   FALSE, etBOOL, {&bIONS}, "compute molecular hyperpolarizability when ions are present"},
        { "-cn",     FALSE, etSTR, {&catname}, "name of cation"},
        { "-an",     FALSE, etSTR, {&anname}, "name of anion"},
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances. Without PBC the maximum range will be three times the largest box edge." },
        { "-ng",       FALSE, etINT, {&ngroups}, 
          "Number of secondary groups to compute RDFs around a central group" },
    };
#define NPA asize(pa)
    const char        *fnTPS, *fnNDX, *fnMAP , *fnKRR , *fnGRD , *fnPOT ;
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
        { efDAT, "-krrpar",  "krr_param.dat", ffOPTRD},
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
    fnMAP = opt2fn_null("-emap", NFILE,fnm);
    fnKRR = opt2fn_null("-krrpar", NFILE,fnm);
    fnGRD = opt2fn_null("-vgrid", NFILE,fnm);
    fnPOT = opt2fn_null("-vinp", NFILE,fnm);

    if (!fnTPS && !fnNDX)
    {
        gmx_fatal(FARGS, "Neither index file nor topology file specified\n"
                  "Nothing to do!");
    }
 
    if (!fnMAP && !fnKRR)
    {
       gmx_fatal(FARGS, "neither electrostatic map file or potential kernel file specified\n");
    }

    if (fnMAP && fnKRR)
    {
       gmx_fatal(FARGS, "specify either kernel ridge reduction files or electrostatic map files\n");
    }

    if (fnKRR || fnGRD || fnPOT)
    {
       if (!fnKRR || !fnGRD || !fnPOT)
       {
           gmx_fatal(FARGS, "specify all the kernel ridge reduction files\n");
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
   
    do_eshs(top, ftp2fn(efTRX, NFILE, fnm),
           fnMAP, fnKRR, fnGRD, fnPOT, opt2fn("-o", NFILE, fnm), opt2fn("-otheta", NFILE, fnm), methodt[0], 
           bIONS, catname, anname, bPBC, qbin, nbinq,
           koutx, kouty, koutz, kinx, kiny, kinz, nbintheta, nbingamma, pin_angle, pout_angle, 
           electrostatic_cutoff, kernstd, gnx, grpindex, grpname, ngroups, oenv);

    return 0;
}

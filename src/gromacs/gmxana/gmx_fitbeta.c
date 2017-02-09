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


static void do_fitbeta(t_topology *top, /*const char *fnNDX, const char *fnTPS,*/ const char *fnTRX,
                   const char *fnTENSOR, const char *fnINCTENSOR, const char *fnTIMEEVOLTENSOR, 
                   const char *fnTHETA, const char *fnQSWIPE,  const char *method,
                   gmx_bool bPBC, int nbinq, int nbintheta, int nbingamma ,real pin_angle, real pout_angle ,
                   real bmzxx, real bmzyy, real bmzzz,
                   int *isize, int  *molindex[], char **grpname, int ng,
                   const output_env_t oenv, const int nthr)
{
    t_trxstatus   *status, **trxin;
    t_trxframe     *fr;
    gmx_bool       bHaveFrame;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, i, j, k, qq, n, c, tt, rr, nframes, nfaces;
    real           qnorm, maxq, gamma = 0.0 ;
    int            max_i, isize0, ind0;
    real           *********onsite_term, *******cos_scattering_ampl, *******sin_scattering_ampl;
    real           **dummy_ptr_bast, *********tot_tensor_squared, *********incoh_tensor_squared;
    real           t, rmax2, rmax,  r_dist, r2,  dq, invhbinw, invnormfac, norm_x, norm_z, mod_f, inv_width;
    real           invsize0, invgamma, theta=0, *theta_vec;
    rvec          *x, dx,  *x_i1, xi, x01, x02,  **arr_qvec_faces ,vec_polin, vec_polout, ***vec_pout_theta_gamma, ***vec_pin_theta_gamma;
    rvec           pol_perp, pol_par,vec_kout, vec_2kin, pol_in1, pol_in2, vec_kout_2kin ; 
    real           invvol, invvol_sum, rho;
    matrix         box, box_pbc;
    int            ePBC = -1, ePBCrdf = -1;
    t_block       *mols = NULL;
    t_atom        *atom = NULL;
    t_pbc          pbc;
    gmx_rmpbc_t    gpbc = NULL;
    int            mol, a;
    FILE          *fptime;


    atom = top->atoms.atom;
    mols = &(top->mols);
    isize0 = isize[0];
    nfaces = 6;
    invsize0 = 1.0/isize0;
    invgamma = 1.0/nbingamma;
    invnormfac = invsize0*invgamma/3.0;
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);

    if (nthr > 1)
    {
       snew(fr, nthr);
       snew(trxin,nthr); 
    }

    for (j = 0; j < nthr; j++)
    {
       read_first_frame(oenv, &trxin[j],fnTRX, &fr[j],TRX_READ_X) ;
    }

    for (i = 0; i <nthr; i++)
    {
      for (j = 0; j < i; j++)
      {
          bHaveFrame=read_next_frame(oenv, trxin[i], &fr[i]);
      }
    }

    fprintf(stderr,"next x %f %f %f\n time %f\n",fr[0].x[0][0],fr[0].x[0][1],fr[0].x[0][2],fr[0].time);
    fprintf(stderr,"next x %f %f %f\n time %f\n",fr[1].x[0][0],fr[1].x[0][1],fr[1].x[0][2],fr[1].time);



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
    
    /*here find the normalization of the OH bond-length and of the dipole vector*/
    set_pbc(&pbc, ePBCrdf, box_pbc);
    ind0  = mols->index[molindex[0][0]];
    pbc_dx(&pbc, x[ind0], x[ind0+1], x01);
    pbc_dx(&pbc, x[ind0], x[ind0+2], x02);
    rvec_sub( x01, x02, xi);
    norm_x = gmx_invsqrt(norm2(xi));
    rvec_add( x01, x02, xi);
    norm_z = gmx_invsqrt(norm2(xi));
    rmax     = sqrt(rmax2);

    max_i = isize[0];

    for (g = 0; g < ng; g++)
    {
        qnorm = M_PI*2.0/(rmax*2.0);
        snew(arr_qvec_faces,nfaces);
        snew(vec_pout_theta_gamma,nfaces);
        snew(vec_pin_theta_gamma, nfaces);
        snew(theta_vec, nbintheta);

        // this loop is to get the scattering wave-vector and polarization vectors for different faces of the cube
        printf("\n----INITIALIZING DIRECTION OF SCATTERED WAVE VECTOR: q=kout-2kin------------------------------\n");
        printf("----q=(2pi/L)*(Nx,Ny,Nz), where Nx, Ny and Nz are integers-------------------------------------\n");
        printf("----COMPUTE THE POLARIZATION VECTORS AT DIFFERENT SCATTERING ANGLES THETA.---------------------\n");
        printf("----THETA IS DEFINED AS THE ANGLE BETWEEN THE INCOMING AND OUTCOMING WAVE-VECTORS.-------------\n");              
        printf("----THE POLARIZATION VECTORS ARE ALSO COMPUTED AT DIFFERENT ANGLES GAMMA.----------------------\n");
        printf("----GAMMA IS THE ANGLE DEFINED BY A PLANE THAT GOES THROUGH THE SCATTERING WAVE-VECTOR AND-----\n");
        printf("----THE PLANE PARALLEL TO THE CHOSEN FACE OF THE SIMULATION BOX--------------------------------\n \n");
        for (rr = 0; rr< nfaces; rr++)
        {
           snew(arr_qvec_faces[rr],nbinq);
           for (qq = 0; qq< nbinq; qq++)
           {
              //now the magnitude of the wave-vector is indeed sqrt(2)*2pi/L*n so we don't need to multiply by sqrt(2).
              arr_qvec_faces[rr][qq][XX] = (qnorm + qnorm*qq)*(1.0) ;
              arr_qvec_faces[rr][qq][YY] = 0.0 ;
              arr_qvec_faces[rr][qq][ZZ] = (qnorm + qnorm*qq)*(-1.0) ;
              rotate_wave_vec(arr_qvec_faces[rr][qq], rr, arr_qvec_faces[rr][qq]);
              snew(vec_pout_theta_gamma[rr], nbintheta);
              snew(vec_pin_theta_gamma[rr], nbintheta);
 
              for(tt = 0 ; tt < nbintheta; tt++)
              {
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
                          svmul(cos(M_PI/180.0*pout_angle), pol_par,  pol_par);
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
                          svmul(cos(M_PI/180.0*pin_angle), pol_par , pol_par);
                          rvec_add(pol_perp, pol_par, vec_pin_theta_gamma[rr][tt][c]);
                          unitv(vec_pin_theta_gamma[rr][tt][c], vec_pin_theta_gamma[rr][tt][c]);
                          printf("polarization vectors at angles theta_expt = %f gamma = %f and face index %d \n",(M_PI + 2.0*theta_vec[tt])*180.0/M_PI, gamma*180.0/M_PI, rr);
                          printf("incoming polarization vector = %f %f %f \n",vec_pin_theta_gamma[rr][tt][c][XX], vec_pin_theta_gamma[rr][tt][c][YY], vec_pin_theta_gamma[rr][tt][c][ZZ]);
                          printf("outcoming polarization vector = %f %f %f \n",vec_pout_theta_gamma[rr][tt][c][XX], vec_pout_theta_gamma[rr][tt][c][YY], vec_pout_theta_gamma[rr][tt][c][ZZ]);
                          printf("direction of scattered wave-vector = %f %f %f \n", vec_kout[XX] -vec_2kin[XX], vec_kout[YY] - vec_2kin[YY], vec_kout[ZZ] -vec_2kin[ZZ]);
                      }
                  }
           }
           printf("wave_vec %f %f %f\n", arr_qvec_faces[rr][0][XX], arr_qvec_faces[rr][0][YY], arr_qvec_faces[rr][0][ZZ]);
        }
    
    Allocate_Scattering_Intensity(nfaces,  nbintheta,  nbinq,
                                  &tot_tensor_squared, &incoh_tensor_squared);
    }
    
    nframes    = 0;
    invvol_sum = 0;
    if (bPBC && (NULL != top))
    {
        gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms);
    }

    if (fnTIMEEVOLTENSOR)
    {
      fptime = xvgropen(fnTIMEEVOLTENSOR, " time (ps) theta (degrees) <Scat_ampl_ijk*Scatt_ampl_lmn(time)>", " ", " ", oenv);
      sprintf(refgt,"%s", "");
      fprintf(fptime, "@ subtitle \"%s%s - %s\"\n", grpname[0], refgt, grpname[1]);
    }

    if (method[0] == 's')
    {
        fprintf(stderr,"use sumexp method, compute the 729 elements of the tensor containing orientational correlations of the molecules\n");
        do
        {
            for (nt = 0; nt < nthr; nt++)
            {
               // Must init pbc every step because of pressure coupling 
               copy_mat(box, box_pbc);
               if (top != NULL)
               {
                   gmx_rmpbc(gpbc, natoms, fr[0].box, fr[0].x);
               }
               set_pbc(&pbc, ePBCrdf, box_pbc);
               invvol      = 1/det(box_pbc);
               invvol_sum += invvol;
               for (g = 0; g < ng; g++)
               {
                   Allocate_scattering_amplitude( nfaces, nbintheta, nbingamma, nbinq,  &onsite_term, &cos_scattering_ampl, &sin_scattering_ampl);
                   for (i = 0; i < isize0; i++)
                   {
                       ind0  = mols->index[molindex[g][i]];
                       copy_rvec(fr[0].x[ind0], xi);
                       pbc_dx(&pbc, xi, fr[0].x[ind0+1], x01);
                       pbc_dx(&pbc, xi, fr[0].x[ind0+2], x02);
                       //this function computes cos_scattering_ampl and sin_scattering_ampl and the incoherent term
                       //summed up over the molecules
                       Projected_Scattering_Amplitude(nfaces, nbintheta, nbingamma, nbinq,
                                                      xi, x01, x02, arr_qvec_faces,
                                                      vec_pout_theta_gamma, vec_pin_theta_gamma, 
                                                      onsite_term, cos_scattering_ampl, sin_scattering_ampl,
                                                      &onsite_term, &cos_scattering_ampl, &sin_scattering_ampl);
                   }
                   // this function computes tot_tensor_squared += cos_scattering_ampl(i,j,k)*cos_scattering_ampl(iprime,jprime,kprime) + 
                   // sin_scattering_ampl(i,j,k)*sin_scattering_ampl(iprime,jprime,kprime) +
                   // and also incoh_tensor_squared += incoherent_term
                   // the sum runs over frames
                   Scattering_Intensity_t(nfaces, t ,nbintheta,  nbingamma ,nbinq, invnormfac, invsize0, theta_vec,
                                          onsite_term, cos_scattering_ampl, sin_scattering_ampl,
                                          tot_tensor_squared, incoh_tensor_squared,
                                          &tot_tensor_squared, &incoh_tensor_squared, fnTIMEEVOLTENSOR, fptime);          
                   Free_scattering_amplitude(nfaces, nbintheta, nbingamma, onsite_term, cos_scattering_ampl, sin_scattering_ampl);
               }
            }
            nframes++;
            do
            {
               bHaveFrame = read_next_frame(oenv, trxin[nt], &fr[nt]);
               k++;
            }
            while(bHaveFrame || k < nthr);
            //read_next_x(oenv, status, &t, x, box)
        }
        while (bHaveFrame);
    }

    if (fnTIMEEVOLTENSOR)
    {
        gmx_ffclose(fptime); 
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
 
    Print_tensors( nbintheta,    nframes, invgamma, theta_vec, nbinq, arr_qvec_faces,
                  tot_tensor_squared, incoh_tensor_squared,  fnTENSOR, fnINCTENSOR, fnQSWIPE, grpname ,oenv);
 
    if (fnTHETA)
    {
        Print_scattering_pattern( nbintheta,  nframes, invgamma,
                                 bmzxx, bmzzz,  bmzyy,  theta_vec,
                                tot_tensor_squared, incoh_tensor_squared,
                                fnTHETA, grpname, oenv);
    }

    Free_Scattering_Intensity( nfaces, nbintheta,
                              tot_tensor_squared, incoh_tensor_squared); 

}

void Print_scattering_pattern(const int nt,  const int nframes, const real invgamma,
                               real bmzxx, real bmzzz, real bmzyy, real *theta_vec,
                               real *********tot_tensor_squared, real *********incoh_tensor_squared, 
                               const char *fnTHETA,  char **grpname, const output_env_t oenv)
{
     int pm, qm, sm, ppm, qpm, spm, qq;
     int i, j, k, p, nplots = 1;
     int rr, tt;
     real ***beta_mol, intensity_theta, theta, intensity_inc_theta;
     FILE *fp;
     char title[STRLEN], gtitle[STRLEN], refgt[30];

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

     beta_mol[2][0][0] = bmzxx ; 
     beta_mol[2][2][2] = bmzzz ;
     beta_mol[2][1][1] = bmzyy ; 
     beta_mol[0][0][2] = bmzxx ;
     beta_mol[0][2][0] = bmzxx ;
     beta_mol[1][1][2] = bmzyy ;
     beta_mol[1][2][1] = bmzyy ;
 
     sprintf(gtitle,"Second harmonic scattering pattern");
     fp = xvgropen(fnTHETA, "I(theta)", "theta (degree)", "I(theta)", oenv);
     sprintf(refgt, "%s", "");
     fprintf(fp, "@ subtitle \"%s%s - %s\"\n", grpname[0], refgt, grpname[1]);
     fprintf(fp, "@    s%d legend \"total\"\n",nplots);
     fprintf(fp, "@target G0.S%d\n",nplots);
     fprintf(fp, "@type xy\n");

     for (p = 0; p < 2; p++)
     {     
         for (tt = 0; tt < nt; tt++ )
         {
             intensity_theta = 0.0;
             intensity_inc_theta = 0.0 ;
             theta = (M_PI + 2.0*theta_vec[tt])*180.0/M_PI;
             for (pm = 0; pm < DIM; pm++)
             {
                 for (qm = 0; qm < DIM; qm++)
                 {
                     for (sm = 0; sm < DIM; sm++)
                     {
                          for (ppm = 0; ppm < DIM; ppm++)
                          {
                              for (qpm = 0; qpm < DIM; qpm++)
                              {
                                  for (spm = 0; spm < DIM; spm++)
                                  {
                                      if (p == 1 )
                                      {
                                          //face 1,3,5 correspond to the setup where the scattering wave-vector lies along the three orthogonal sides of the simulation box
                                          //face 0,2,4 correspond to the setup where the scattering wave-vector lies along the three orthogonal diagonals of the faces the simulation box
                                          intensity_theta += beta_mol[pm][qm][sm]*beta_mol[ppm][qpm][spm]*(tot_tensor_squared[1][tt][pm][qm][sm][ppm][qpm][spm][0]
                                                             + tot_tensor_squared[3][tt][pm][qm][sm][ppm][qpm][spm][0] +
                                                             tot_tensor_squared[5][tt][pm][qm][sm][ppm][qpm][spm][0]);
                                      }
                                      else
                                      {
                                         intensity_inc_theta += beta_mol[pm][qm][sm]*beta_mol[ppm][qpm][spm]*(incoh_tensor_squared[1][tt][pm][qm][sm][ppm][qpm][spm][0]
                                                           + incoh_tensor_squared[3][tt][pm][qm][sm][ppm][qpm][spm][0] +
                                                           incoh_tensor_squared[5][tt][pm][qm][sm][ppm][qpm][spm][0]);

                                      }
                                  }
                              }
                          }
                     }
                 }
             }
             fprintf(fp, "%10g", theta);
             if (p == 1)
             {
                 fprintf(fp, " %10g", intensity_theta/(nframes*3.0)*invgamma);
             }
             else
             {
                 fprintf(fp, " %10g", intensity_inc_theta/(nframes*3.0)*invgamma);
             }

             fprintf(fp, "\n");
         }
         if (p == 0)
         {
            nplots++;
            fprintf(fp, "&\n");
            fprintf(fp, "@    s%d legend \"incoherent\"\n",nplots);
            fprintf(fp, "@target G0.S%d\n",nplots);
            fprintf(fp, "@type xy\n");
         }
     }
     fprintf(fp, "&\n");
     gmx_ffclose(fp);

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

void Print_tensors(const int nt,  const int nframes, const real invgamma, real *theta_vec ,int nbinq, 
                   rvec **arr_qvec_faces, real *********tot_tensor_squared, 
                   real *********incoh_tensor_squared, const char *fnTENSOR, 
                   const char *fnINCTENSOR , const char *fnQSWIPE, char **grpname ,const output_env_t oenv)
{
    int pm, qm, sm, ppm, qpm, spm, qq;
    int rr, tt, thetain;
    char refgt[30],str[80],thetastr[80],strin[80];
    real theta;
    FILE *fpn, *fpp, *fq, *fqinc;

    fpn = xvgropen(fnTENSOR, "<Scat_ampl_ijk*Scatt_ampl_lmn>", " ", " ", oenv);
    sprintf(refgt,"%s", "");
    fprintf(fpn, "@ subtitle \"%s%s - %s\"\n", grpname[0], refgt, grpname[1]);
    if (fnINCTENSOR)
    {
       fpp = xvgropen(fnINCTENSOR, "<Inc_Scatt_ampl_ijk*Inc_Scatt_ampl_lmn>", " " , " ", oenv);
       sprintf(refgt,"%s", "");
       fprintf(fpp, "@ subtitle \"%s%s - %s\"\n", grpname[0], refgt, grpname[1]);
    }
    // average over the intensities computed at a wave-vector that lies along the diagonal of the three faces of the simulation box
    for (tt = 0; tt < nt; tt++ )
    {
        theta = (M_PI + 2.0*theta_vec[tt])*180.0/M_PI;
        fprintf(fpn, "%10g ",theta);
        if (fnINCTENSOR)
        {
           fprintf(fpp, "%10g ",theta);
        }
        for (pm = 0; pm < DIM; pm++)
        {
            for (qm = 0; qm < DIM; qm++)
            {
                for (sm = 0; sm < DIM; sm++)
                {
                     for (ppm = 0; ppm < DIM; ppm++)
                     {
                         for (qpm = 0; qpm < DIM; qpm++)
                         {
                             for (spm = 0; spm < DIM; spm++)
                             {
                                 for (qq = 0; qq < 1; qq++)
                                 {
                                    fprintf(fpn, "%10g ",(tot_tensor_squared[1][tt][pm][qm][sm][ppm][qpm][spm][qq]
                                                        + tot_tensor_squared[3][tt][pm][qm][sm][ppm][qpm][spm][qq]+
                                                        tot_tensor_squared[5][tt][pm][qm][sm][ppm][qpm][spm][qq])/(nframes*3.0)*invgamma);
                                    if (fnINCTENSOR)
                                    {
                                       fprintf(fpp, "%10g ",(incoh_tensor_squared[1][tt][pm][qm][sm][ppm][qpm][spm][qq]
                                                        + incoh_tensor_squared[3][tt][pm][qm][sm][ppm][qpm][spm][qq]+
                                                        incoh_tensor_squared[5][tt][pm][qm][sm][ppm][qpm][spm][qq])/(nframes*3.0)*invgamma);
                                    }
                                 }
                             }
                         }
                     }
                }
            }
        }
        fprintf(fpn, "\n");
        if (fnINCTENSOR)
           {
              fprintf(fpp, "\n");
           }
    }
    gmx_ffclose(fpn);
    if (fnINCTENSOR)
    {
       gmx_ffclose(fpp);
    }
    if (fnQSWIPE)
    {
        for (tt = 0; tt < nt; tt++ )
        {
           theta = (M_PI + 2.0*theta_vec[tt])*180.0/M_PI;
           thetain=roundf(theta);
           sprintf(thetastr,"%d",thetain);
           strcpy(str,fnQSWIPE);
           str[strlen(str)-1]=0;
           str[strlen(str)-1]=0;
           str[strlen(str)-1]=0;
           str[strlen(str)-1]=0;
           strcpy(strin,str);
           strcat(strin,"_Inctheta");
           strcat(strin,thetastr);
           strcat(strin,".xvg");
           strcat(str,"_theta");
           strcat(str,thetastr);
           strcat(str,".xvg");
           fq = xvgropen(str, "<Scat_ampl_ijk*Scatt_ampl_lmn> as a function of q"," ", " ", oenv);
           fqinc = xvgropen(strin, "<Inc_Scat_ampl_ijk*Inc_Scatt_ampl_lmn> as a function of q"," ", " ", oenv);
           sprintf(refgt,"%s", "");
           fprintf(fq, "@ subtitle \"%s%s - %s\"\n", grpname[0], refgt, grpname[1]);
           fprintf(fqinc, "@ subtitle \"%s%s - %s\"\n", grpname[0], refgt, grpname[1]);

           for (qq = 0; qq < nbinq; qq++)
           {
               fprintf(fq, "%10g ",norm(arr_qvec_faces[1][qq]));
               fprintf(fqinc, "%10g ",norm(arr_qvec_faces[1][qq]));

               for (pm = 0; pm < DIM; pm++)
               {
                   for (qm = 0; qm < DIM; qm++)
                   {
                       for (sm = 0; sm < DIM; sm++)
                       {
                            for (ppm = 0; ppm < DIM; ppm++)
                            {
                                for (qpm = 0; qpm < DIM; qpm++)
                                {
                                    for (spm = 0; spm < DIM; spm++)
                                    {
                                        fprintf(fq, "   %10g", (tot_tensor_squared[1][tt][pm][qm][sm][ppm][qpm][spm][qq]
                                                             + tot_tensor_squared[3][tt][pm][qm][sm][ppm][qpm][spm][qq]+
                                                              tot_tensor_squared[5][tt][pm][qm][sm][ppm][qpm][spm][qq])/(nframes*3.0)*invgamma); 

                                        fprintf(fqinc, "   %10g", (incoh_tensor_squared[1][tt][pm][qm][sm][ppm][qpm][spm][qq]
                                                             + incoh_tensor_squared[3][tt][pm][qm][sm][ppm][qpm][spm][qq]+
                                                              incoh_tensor_squared[5][tt][pm][qm][sm][ppm][qpm][spm][qq])/(nframes*3.0)*invgamma);

                                    }
                                }
                            }
                       }
                   }
               }
               fprintf(fq, "\n");
               fprintf(fqinc, "\n");
           }
           gmx_ffclose(fq); 
           gmx_ffclose(fqinc);          
        }
    }



}

void Projected_Scattering_Amplitude(const int nf, const int nt, const int nga, const int nq,
     const rvec xi, const rvec xv2, const rvec xv3,  rvec **arr_scatt_wave_vec,
     rvec ***pout_theta_gamma, rvec ***pin_theta_gamma, 
     real *********Onsite_term, real *******Cos_scatt_ampl, real *******Sin_scatt_ampl, 
     real **********Onsite_term_addr, real ********Cos_scatt_ampl_addr, real ********Sin_scatt_ampl_addr)
{
    // this function computes the 27 elements: Scos[pm][qm][sm] = sum_ijk pout[i]*pin[j]*pin[k] c[i][pm] c[j][qm] c[k][sm]*cos(q*r)
    // Ssin[pm][qm][sm] = sum_ijk pout[i]*pin[j]*pin[k] c[i][pm] c[j][qm] c[k][sm]*sin(q*r)
    // for each configuration defined by the angles theta, gamma, and cube index (face or diagonal)
    // and for each combination of incoming and outcoming polarization vectors
    // i,j,k are indexes of rotation matrix in lab frame
    // it takes as input the inverse of the norm of the length of the H1-H2 vector and of the O-H1 vector (it has to be changed for polarizable FFs)
    // also it takes as input the direction of OH1 and of OH2
    // and the polarization combinations
    // and the number of bins over gamma, theta and cube indexes (always 6)
    // rr, tt, c are indexes relating to the given side or diagonal of the cube that is sampled (rr),
    // or to the scattering angle (tt), and to the scattering plane (c)

    int i, j, k, pm, qm, sm, rr, tt, c, qq;
    rvec xvec, yvec, zvec;
    matrix cosdirmat, cosdirmat_t;
    real ***rot_element_matrix_prime;
    real co_q_xi, si_q_xi, q_xi, rot_element_matrix_pqs ;

    rvec_add( xv2, xv3, zvec);
    cprod(zvec,xv2, yvec);
    unitv(zvec,zvec);
    unitv(yvec,yvec);
    cprod(yvec,zvec,xvec);

    cosdirmat[0][0] = xvec[0]; cosdirmat[0][1] = yvec[0]; cosdirmat[0][2] = zvec[0];
    cosdirmat[1][0] = xvec[1]; cosdirmat[1][1] = yvec[1]; cosdirmat[1][2] = zvec[1];
    cosdirmat[2][0] = xvec[2]; cosdirmat[2][1] = yvec[2]; cosdirmat[2][2] = zvec[2];

    //to compute the scattering amplitude we use the transpose of the rotation matrix
    //in order to have loops with the indexes in the right order efficiencywise
    transpose(cosdirmat, cosdirmat_t);

    snew(rot_element_matrix_prime,DIM);
    for (i = 0; i < DIM; i++)
    {
        snew(rot_element_matrix_prime[i],DIM);
        for ( j = 0; j < DIM; j++)
        {
           snew(rot_element_matrix_prime[i][j],DIM);
        }
    }
 
    for (rr = 0; rr < nf; rr++)
    {
        for (tt = 0; tt < nt; tt++ )
        {
            for (c  = 0; c < nga; c++)
            {
                for (pm = 0; pm < DIM; pm++)
                {
                    for (qm = 0; qm < DIM; qm++)
                    {
                        for (sm = 0; sm < DIM; sm++)
                        {    
                             rot_element_matrix_pqs = 0.0;                    
                             for (i = 0; i < DIM; i++)
                             {
                                 for (j = 0; j < DIM; j++)
                                 {
                                     for (k = 0; k < DIM; k++)
                                     {              
                                          rot_element_matrix_pqs += pout_theta_gamma[rr][tt][c][i]*pin_theta_gamma[rr][tt][c][j]*pin_theta_gamma[rr][tt][c][k]
                                                                                  *cosdirmat_t[pm][i]*cosdirmat_t[qm][j]*cosdirmat_t[sm][k] ;
                                     }
                                 }
                             }
                             rot_element_matrix_prime[pm][qm][sm] = rot_element_matrix_pqs;
                             for (qq = 0; qq < nq; qq++)
                             {
                                 q_xi = iprod(xi,arr_scatt_wave_vec[rr][qq]);
                                 co_q_xi = cos(q_xi);
                                 si_q_xi = sin(q_xi);
                                 Cos_scatt_ampl[rr][tt][c][pm][qm][sm][qq] += rot_element_matrix_pqs*co_q_xi ;
                                 Sin_scatt_ampl[rr][tt][c][pm][qm][sm][qq] += rot_element_matrix_pqs*si_q_xi ;
                             }
                        }
                    }
                }
                for (pm = 0; pm < DIM; pm++)
                {
                    for (qm = 0; qm < DIM; qm++)
                    {
                        for (sm = 0; sm < DIM; sm++)
                        {
                             for (i = 0; i < DIM; i++)
                             {
                                 for (j = 0; j < DIM; j++)
                                 {
                                     for (k = 0; k < DIM; k++)
                                     {
                                          Onsite_term[rr][tt][c][pm][qm][sm][i][j][k] += rot_element_matrix_prime[pm][qm][sm]*rot_element_matrix_prime[i][j][k] ;
                                     }
                                 }
                             }

                        }
                    }
                }

            }
        }
    }
    *Onsite_term_addr = Onsite_term;
    *Cos_scatt_ampl_addr = Cos_scatt_ampl;
    *Sin_scatt_ampl_addr = Sin_scatt_ampl; 
    for (i = 0; i < DIM; i++)
    {
        for ( j = 0; j < DIM; j++)
        {
           sfree(rot_element_matrix_prime[i][j]);
        }
        sfree(rot_element_matrix_prime[i]);
    }
    sfree(rot_element_matrix_prime);

    clear_mat(cosdirmat);
    clear_mat(cosdirmat_t);

}


void Allocate_scattering_amplitude(const int nf, const int nt, const int nga, const int nq, real **********Incoh_term_addr, real ********Cos_scatt_ampl_addr, real ********Sin_scatt_ampl_addr)
{

    int rr, tt, c, qq;
    int pm, qm, sm, ppm, qpm;
    real *********Incoh_term, *******Cos_scatt_ampl, *******Sin_scatt_ampl;
 
    snew(Incoh_term, nf);
    snew(Cos_scatt_ampl, nf);
    snew(Sin_scatt_ampl, nf);
    for (rr = 0; rr < nf; rr++)
    {
        snew(Incoh_term[rr],nt);
        snew(Cos_scatt_ampl[rr], nt);
        snew(Sin_scatt_ampl[rr], nt);
        for (tt = 0; tt < nt; tt++ )
        {
            snew(Incoh_term[rr][tt],nga);
            snew(Cos_scatt_ampl[rr][tt], nga);
            snew(Sin_scatt_ampl[rr][tt], nga);
            for (c  = 0; c < nga; c++)
            {
                snew(Incoh_term[rr][tt][c],DIM);
                snew(Cos_scatt_ampl[rr][tt][c], DIM);
                snew(Sin_scatt_ampl[rr][tt][c], DIM);
                for (pm = 0; pm < DIM; pm++)
                {
                    snew(Incoh_term[rr][tt][c][pm],DIM);
                    snew(Cos_scatt_ampl[rr][tt][c][pm], DIM);
                    snew(Sin_scatt_ampl[rr][tt][c][pm], DIM);
                    for (qm = 0; qm < DIM; qm++)
                    {
                        snew(Incoh_term[rr][tt][c][pm][qm],DIM);
                        snew(Cos_scatt_ampl[rr][tt][c][pm][qm], DIM);
                        snew(Sin_scatt_ampl[rr][tt][c][pm][qm], DIM);
                        for (sm = 0; sm < DIM; sm++)
                        {
                           snew(Cos_scatt_ampl[rr][tt][c][pm][qm][sm], nq);
                           snew(Sin_scatt_ampl[rr][tt][c][pm][qm][sm], nq);
                           snew(Incoh_term[rr][tt][c][pm][qm][sm],DIM);
                           for (ppm = 0; ppm < DIM; ppm++)
                           {
                                snew(Incoh_term[rr][tt][c][pm][qm][sm][ppm],DIM);
                                for (qpm = 0; qpm < DIM; qpm++)
                                {
                                    snew(Incoh_term[rr][tt][c][pm][qm][sm][ppm][qpm],DIM);
                                }
                           }
                        }
                    }
                }
            }
        }
    }
    *Cos_scatt_ampl_addr = Cos_scatt_ampl;
    *Sin_scatt_ampl_addr = Sin_scatt_ampl;
    *Incoh_term_addr = Incoh_term;
}


void Free_scattering_amplitude(const int nf, const int nt, const int nga, real *********Incoh_term, real *******Cos_scatt_ampl, real *******Sin_scatt_ampl)
{

    int rr, tt, c;
    int pm, qm, sm;
    int ppm, qpm;

    for (rr = 0; rr < nf; rr++)
    {
        for (tt = 0; tt < nt; tt++ )
        {
            for (c  = 0; c < nga; c++)
            {
                for (pm = 0; pm < DIM; pm++)
                {
                    for (qm = 0; qm < DIM; qm++)
                    {
                        for (sm = 0; sm < DIM; sm++)
                        {
                           for (ppm = 0; ppm < DIM; ppm++)
                           {
                               for (qpm = 0; qpm < DIM; qpm++)
                               {
                                  sfree(Incoh_term[rr][tt][c][pm][qm][sm][ppm][qpm]);
                               }
                               sfree(Incoh_term[rr][tt][c][pm][qm][sm][ppm]);
                           }
                           sfree(Cos_scatt_ampl[rr][tt][c][pm][qm][sm]);
                           sfree(Sin_scatt_ampl[rr][tt][c][pm][qm][sm]);
                        }
                        sfree(Incoh_term[rr][tt][c][pm][qm]);
                        sfree(Cos_scatt_ampl[rr][tt][c][pm][qm]);
                        sfree(Sin_scatt_ampl[rr][tt][c][pm][qm]);
                    }
                    sfree(Incoh_term[rr][tt][c][pm]);
                    sfree(Cos_scatt_ampl[rr][tt][c][pm]);
                   sfree(Sin_scatt_ampl[rr][tt][c][pm]);
                }
                sfree(Incoh_term[rr][tt][c]);
                sfree(Cos_scatt_ampl[rr][tt][c]);
                sfree(Sin_scatt_ampl[rr][tt][c]);
           }
            sfree(Incoh_term[rr][tt]);
            sfree(Cos_scatt_ampl[rr][tt]);
            sfree(Sin_scatt_ampl[rr][tt]);
        }
        sfree(Incoh_term[rr]);
        sfree(Cos_scatt_ampl[rr]);
        sfree(Sin_scatt_ampl[rr]);
    }
    sfree(Incoh_term);
    sfree(Cos_scatt_ampl);
    sfree(Sin_scatt_ampl);

}


void Scattering_Intensity_t(const int nf, real time, const int nt, const int nga, const int nq,
                               const real inversenorm, const real inversesize,  const  real *theta_vec,
                               real *********Incoh_term, real *******Cos_scatt_ampl, real *******Sin_scatt_ampl,
                               real  *********tot_tensor_squared,      real *********incoh_tensor_squared,
                               real **********tot_tensor_squared_addr, real **********incoh_tensor_squared_addr, 
                               const char *fnTIMEEVOLTENSOR, FILE  *fname)
{ 
    // this function computes the 27*27 elements cos_scatt_ampl[pm][qm][sm]*cos_scatt_ampl[ppm][qpm][spm] + sin_scatt_ampl[pm][qm][sm]*sin_scatt_ampl[ppm][qpm][spm]
    int i, j, k, pm, qm, sm, ppm, qpm, spm, rr, tt, c, qq, pl;
    real co_q_xi, si_q_xi, q_xi;
    real incoh_temp, tot_temp ;
    real print_temp = 0.0;
    real theta;


    if (fnTIMEEVOLTENSOR)
    {
       fprintf(fname, "%10g ", time);
    }
    for (rr = 0; rr < nf; rr++)
    {
        for (tt = 0; tt < nt; tt++ )
        {
            if (fnTIMEEVOLTENSOR && rr==0)
            {
               theta = (M_PI + 2.0*theta_vec[tt])*180.0/M_PI;
               fprintf(fname, "%10g ",theta);
            }
            for (c  = 0; c < nga; c++)
            {
                for (pm = 0; pm < DIM; pm++)
                {
                    for (qm = 0; qm < DIM; qm++)
                    {
                        for (sm = 0; sm < DIM; sm++)
                        {
                             for (ppm = 0; ppm < DIM; ppm++)
                             {
                                 for (qpm = 0; qpm < DIM; qpm++)
                                 {
                                     for (spm = 0; spm < DIM; spm++)
                                     {
                                         if ((fnTIMEEVOLTENSOR) && (c==0) && (rr==0) )
                                         {
                                            print_temp = 0.0;
                                            for (pl = 0; pl < nga; pl++)
                                            {
                                              print_temp += (Cos_scatt_ampl[1][tt][pl][pm][qm][sm][0]*Cos_scatt_ampl[1][tt][pl][ppm][qpm][spm][0] +
                                                       Sin_scatt_ampl[1][tt][pl][pm][qm][sm][0]*Sin_scatt_ampl[1][tt][pl][ppm][qpm][spm][qq]);
                                              print_temp +=   (Cos_scatt_ampl[3][tt][pl][pm][qm][sm][0]*Cos_scatt_ampl[3][tt][pl][ppm][qpm][spm][0] +
                                                       Sin_scatt_ampl[3][tt][pl][pm][qm][sm][0]*Sin_scatt_ampl[3][tt][pl][ppm][qpm][spm][0]);
                                              print_temp += (Cos_scatt_ampl[5][tt][pl][pm][qm][sm][0]*Cos_scatt_ampl[5][tt][pl][ppm][qpm][spm][0] +
                                                       Sin_scatt_ampl[5][tt][pl][pm][qm][sm][0]*Sin_scatt_ampl[5][tt][pl][ppm][qpm][spm][0]);
                                            }
                                            fprintf(fname, "%10g ", print_temp*inversenorm);
                                             
                                         }
                                         for (qq = 0; qq < nq; qq++)
                                         {
                                           tot_temp = (Cos_scatt_ampl[rr][tt][c][pm][qm][sm][qq]*Cos_scatt_ampl[rr][tt][c][ppm][qpm][spm][qq] + 
                                                       Sin_scatt_ampl[rr][tt][c][pm][qm][sm][qq]*Sin_scatt_ampl[rr][tt][c][ppm][qpm][spm][qq])*inversesize;
                                           tot_tensor_squared[rr][tt][pm][qm][sm][ppm][qpm][spm][qq]  +=  tot_temp;
                                           incoh_tensor_squared[rr][tt][pm][qm][sm][ppm][qpm][spm][qq]  +=  Incoh_term[rr][tt][c][pm][qm][sm][ppm][qpm][spm]*inversesize;
                                         }
                                    
                                     }
                                 }
                             }
                        }
                    }
                }
            }
        }
    }

    if (fnTIMEEVOLTENSOR)
    {
       fprintf(fname, "\n");
    }

    *tot_tensor_squared_addr = tot_tensor_squared;
    *incoh_tensor_squared_addr = incoh_tensor_squared;

}

void Allocate_Scattering_Intensity(const int nf, const int nt,  const int nq,
                               real **********tot_tensor_squared_addr, real **********incoh_tensor_squared_addr)
{
    int rr, tt, qq;
    int pm, qm, sm;
    int ppm, qpm, spm;
    real *********tot_tensor_squared, *********incoh_tensor_squared;

    snew(tot_tensor_squared, nf);
    snew(incoh_tensor_squared, nf);
    for (rr = 0; rr < nf; rr++)
    {
        snew(tot_tensor_squared[rr], nt);
        snew(incoh_tensor_squared[rr], nt);
        for (tt = 0; tt < nt; tt++ )
        {
            snew(tot_tensor_squared[rr][tt], DIM);
            snew(incoh_tensor_squared[rr][tt], DIM);
            for (pm = 0; pm < DIM; pm++)
            {
                snew(tot_tensor_squared[rr][tt][pm], DIM);
                snew(incoh_tensor_squared[rr][tt][pm], DIM);
                for (qm = 0; qm < DIM; qm++)
                {
                    snew(tot_tensor_squared[rr][tt][pm][qm], DIM);
                    snew(incoh_tensor_squared[rr][tt][pm][qm], DIM);
                    for (sm = 0; sm < DIM; sm++)
                    {
                         snew(tot_tensor_squared[rr][tt][pm][qm][sm], DIM);
                         snew(incoh_tensor_squared[rr][tt][pm][qm][sm], DIM);
                         for (ppm = 0; ppm < DIM; ppm++)
                         {
                             snew(tot_tensor_squared[rr][tt][pm][qm][sm][ppm], DIM);
                             snew(incoh_tensor_squared[rr][tt][pm][qm][sm][ppm], DIM);
                             for (qpm = 0; qpm < DIM; qpm++)
                             {
                                 snew(tot_tensor_squared[rr][tt][pm][qm][sm][ppm][qpm], DIM);
                                 snew(incoh_tensor_squared[rr][tt][pm][qm][sm][ppm][qpm], DIM);
                                 for (spm = 0; spm < DIM; spm++)
                                 {
                                     snew(tot_tensor_squared[rr][tt][pm][qm][sm][ppm][qpm][spm], nq);
                                     snew(incoh_tensor_squared[rr][tt][pm][qm][sm][ppm][qpm][spm], nq);
                                 }
                             }
                         }
                    }
                }
            }
        }
    }
    *tot_tensor_squared_addr = tot_tensor_squared;
    *incoh_tensor_squared_addr = incoh_tensor_squared;
}

void Free_Scattering_Intensity(const int nf, const int nt,
                               real *********tot_tensor_squared, real *********incoh_tensor_squared)
{
    int rr, tt;
    int pm, qm, sm;
    int ppm, qpm, spm;
    for (rr = 0; rr < nf; rr++)
    {
        for (tt = 0; tt < nt; tt++ )
        {
            for (pm = 0; pm < DIM; pm++)
            {
                for (qm = 0; qm < DIM; qm++)
                {
                    for (sm = 0; sm < DIM; sm++)
                    {
                         for (ppm = 0; ppm < DIM; ppm++)
                         {
                             for (qpm = 0; qpm < DIM; qpm++)
                             {
                                 for (spm = 0; spm < DIM; spm++)
                                 {
                                     sfree(tot_tensor_squared[rr][tt][pm][qm][sm][ppm][qpm][spm]);
                                     sfree(incoh_tensor_squared[rr][tt][pm][qm][sm][ppm][qpm][spm]);
                                 }
                                 sfree(tot_tensor_squared[rr][tt][pm][qm][sm][ppm][qpm]);
                                 sfree(incoh_tensor_squared[rr][tt][pm][qm][sm][ppm][qpm]);
                             }
                             sfree(tot_tensor_squared[rr][tt][pm][qm][sm][ppm]);
                             sfree(incoh_tensor_squared[rr][tt][pm][qm][sm][ppm]);
                         }
                         sfree(tot_tensor_squared[rr][tt][pm][qm][sm]);
                         sfree(incoh_tensor_squared[rr][tt][pm][qm][sm]);
                    }
                    sfree(tot_tensor_squared[rr][tt][pm][qm]);
                    sfree(incoh_tensor_squared[rr][tt][pm][qm]);
                }
                sfree(tot_tensor_squared[rr][tt][pm]);
                sfree(incoh_tensor_squared[rr][tt][pm]);
            }
            sfree(tot_tensor_squared[rr][tt]);
            sfree(incoh_tensor_squared[rr][tt]);
        }
        sfree(tot_tensor_squared[rr]);
        sfree(incoh_tensor_squared[rr]);
    }
    sfree(tot_tensor_squared);
    sfree(incoh_tensor_squared);
}


int gmx_fitbeta(int argc, char *argv[])
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
    static gmx_bool    bPBC = TRUE;
    static real        pout_angle = 0.0 , pin_angle = 0.0;
    static real        bmzxx = 5.7 , bmzyy = 10.9 , bmzzz = 31.6 ;
    static int         ngroups = 1, nbintheta = 10, nbingamma = 2 , nbinq = 10, nthr = 1;

    static const char *methodt[] = { NULL, "sumexp" ,NULL }; 

    t_pargs            pa[] = {
        { "-nbintheta",      FALSE, etINT, {&nbintheta},
        "number of bins over scattering angle theta chosen between -pi/2 and + pi/2 (available only with thetaswipe)" },
        { "-nplanes",      FALSE, etINT, {&nbingamma},
        "number of scattering planes that lie on the scattered wave-vector to average over, -PI/2< gamma< PI/2" },
        { "-nbinq",      FALSE, etINT, {&nbinq},
        "how many bins in the reciprocal space" },
        { "-pout",         FALSE, etREAL, {&pout_angle}, "polarization angle of outcoming beam in degrees. For P choose 0, for S choose 90" },
        { "-pin",         FALSE, etREAL, {&pin_angle}, "polarization angle of incoming beam in degrees. For P choose 0, for S choose 90" },
        { "-method",     FALSE, etENUM, {methodt},
          "I(q) using the different methods" },
        { "-bmzxx",         FALSE, etREAL, {&bmzxx}, "component of beta in zxx " },
        { "-bmzyy",         FALSE, etREAL, {&bmzyy}, "component of beta in zyy " },
        { "-bmzzz",         FALSE, etREAL, {&bmzzz}, "component of beta in zzz " },
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances. Without PBC the maximum range will be three times the largest box edge." },
        { "-ng",       FALSE, etINT, {&ngroups}, 
          "Number of secondary groups to compute RDFs around a central group" },
        { "-nthr",       FALSE, etINT, {&nthr},
          "number of openMP threads" },
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
        { efNDX, NULL,  NULL,     ffREAD },
        { efXVG, "-o",  "tot_intensity",    ffWRITE },
        { efXVG, "-oin",  "inc_intensity",    ffOPTWR },
        { efXVG, "-oev",  "time_evol_intensity",    ffOPTWR },
        { efXVG, "-otheta", "non_linear_sfact_vs_theta", ffOPTWR },
        { efXVG, "-oqplot", "non_linear_sfact_vs_q", ffOPTWR },


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

        if (ngroups  != 1)
    {
        gmx_fatal(FARGS, "Intensity for ng>1 not yet implemented rerun with ng = 1 (default) \n");
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
   
    do_fitbeta(top, ftp2fn(efTRX, NFILE, fnm),
           opt2fn("-o", NFILE, fnm), opt2fn_null("-oin", NFILE, fnm), 
           opt2fn_null("-oev", NFILE, fnm), opt2fn_null("-otheta", NFILE, fnm),
           opt2fn_null("-oqplot", NFILE, fnm),
           methodt[0],  bPBC, nbinq,
           nbintheta, nbingamma, pin_angle, pout_angle,
           bmzxx, bmzyy, bmzzz, 
           gnx, grpindex, grpname, ngroups, oenv,nthr);

    return 0;
}

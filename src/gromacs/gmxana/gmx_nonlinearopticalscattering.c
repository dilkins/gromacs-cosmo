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

static void do_nonlinearopticalscattering(t_topology *top, /*const char *fnNDX, const char *fnTPS,*/ const char *fnTRX,
                   const char *fnSFACT, const char *fnOSRDF, const char *fnORDF, /*const char *fnHQ, */
                   const char *method,
                   gmx_bool bPBC, gmx_bool bNormalize,
                   real maxq, real minq, int nbinq, real kx, real ky, real kz, real binwidth,
                   real fade, real faderdf, int *isize, int  *molindex[], char **grpname, int ng,
                   const output_env_t oenv)
{
    FILE          *fp;
    FILE          *fpn;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, i, ii, j, k, nbin, qq, j0, j1, n, n_j ,nframes;
    int          **count;
    real         **s_method, **s_method_g_r, *analytical_integral, *arr_q, qel, *temp_method, *cos_q, *sin_q, beta_lab;
    /*char         **grpname;*/
    int            isize_cm = 0, nrdf = 0, max_i, isize0, isize_g;
    atom_id      /* **index, */ *index_cm = NULL;
    gmx_int64_t   *sum;
    real           t, rmax2, rmax,  r, r_dist, r2, r2ii, q_xi, dq, invhbinw, normfac, mod_f, inv_width;
    real           segvol, spherevol, prev_spherevol, **rdf;
    rvec          *x, dx, *x0 = NULL, *x_i1, xi, arr_qvec, uv1, uv2, uv3;
    real          *inv_segvol, invvol, invvol_sum, rho;
    gmx_bool       bClose, *bExcl, bTop, bNonSelfExcl;
    matrix         box, box_pbc;
    int          **npairs;
    atom_id        ix, jx, ***pairs;
    /*t_topology    *top  = NULL;*/
    int            ePBC = -1, ePBCrdf = -1;
    t_block       *mols = NULL;
    t_blocka      *excl;
    t_atom        *atom = NULL;
    t_pbc          pbc;
    gmx_rmpbc_t    gpbc = NULL;
    int           *is   = NULL, **coi = NULL, cur, mol, i1, res, a;

    excl = NULL;

    bClose = FALSE ; 

    atom = top->atoms.atom;
    mols = &(top->mols);
    isize0 = isize[0];
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    fprintf(stderr,"natoms %d\n",natoms);
    fprintf(stderr,"molindex %d\n",molindex[0][0]);
    fprintf(stderr,"isize[0] %d\n",isize[0]);
    fprintf(stderr,"ng %d\n",ng);

    /* initialize some handy things */
    if (ePBC == -1)
    {
        ePBC = guess_ePBC(box);
    }
    copy_mat(box, box_pbc);
    ePBCrdf = ePBC;
    if (bPBC)
    {
        rmax2   = /*0.99*0.99* */ max_cutoff2(FALSE ? epbcXY : epbcXYZ, box_pbc);
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

    /* We use the double amount of bins, so we can correctly
     * write the rdf and rdf_cn output at i*binwidth values.
     */
    nbin     = (int)(sqrt(rmax2) * 2 / binwidth);
    invhbinw = 2.0 / binwidth;
    rmax     = sqrt(rmax2);

    snew(count, ng);
    snew(pairs, ng);
    snew(npairs, ng);
    snew(s_method, ng);
    snew(s_method_g_r, ng);
    snew(bExcl, natoms);
    max_i = 0;
    for (g = 0; g < ng; g++)
    {
        if (isize[g+1] > max_i)
        {
            max_i = isize[g+1];
        }

        /* this is THE array */
        snew(count[g], nbin+1);
        /* make pairlist array for groups and exclusions */
        snew(pairs[g], isize[0]);
        snew(npairs[g], isize[0]);
        /*allocate memory for s_method array */
        snew(s_method[g], nbinq);
        snew(s_method_g_r[g], nbinq);
        snew(arr_q,nbinq);
        snew(temp_method,nbinq);
        snew(cos_q,nbinq);
        snew(sin_q,nbinq);
        snew(analytical_integral,nbinq);
        normfac = 1.0/sqrt(kx*kx + ky*ky + kz*kz) ;
        arr_qvec[XX] = kx*normfac;
        arr_qvec[YY] = ky*normfac;
        arr_qvec[ZZ] = kz*normfac;
        inv_width = (fade == 0.0 ) ? 1.0 : M_PI*0.5/(rmax-fade) ; 
        dq=(maxq-minq)/nbinq ;       
        for (qq = 0; qq< nbinq; qq++)
        {
            arr_q[qq]=minq+dq*qq;
            if (fade == 0.0)
            {
               analytical_integral[qq]=((sin(arr_q[qq]*rmax) - arr_q[qq]*rmax*cos(arr_q[qq]*rmax))/(arr_q[qq]*arr_q[qq]*arr_q[qq]))*4.0*M_PI*isize0;
            }
            else
            {
               qel = arr_q[qq] ;
               analytical_integral[qq] = 4.0*M_PI*isize0*(-fade*qel*cos(fade*qel) + sin(fade*qel))/pow(qel,3.0); 
               analytical_integral[qq] +=  4.0*M_PI*isize0*((fade*qel*(M_PI + qel*(fade - rmax))*(pow(M_PI,2) - 2*pow(qel,2)*pow(fade - rmax,2))*\
               (M_PI + qel*(-fade + rmax))*cos(fade*qel) -\
               pow(M_PI,2)*qel*(M_PI + qel*(fade - rmax))*rmax*(M_PI + qel*(-fade + rmax))*cos(qel*rmax) + \
               (-pow(M_PI,4) + pow(M_PI,2)*pow(qel,2)*pow(fade - rmax,2) - 2*pow(qel,4)*pow(fade - rmax,4))*sin(fade*qel) + \
               pow(M_PI,2)*(pow(M_PI,2) - 3*pow(qel,2)*pow(fade - rmax,2))*sin(qel*rmax))/\
               (2.*pow(qel,3)*pow(M_PI + qel*(fade - rmax),2)*pow(M_PI + qel*(-fade + rmax),2)));
            }                                                                 
        }
        for (i = 0; i < isize[0]; i++)
        {   
            /* We can only have exclusions with atomic rdfs */
            ix = molindex[0][i];
            /*fprintf(stderr,"ix isize[g+1] isize[0] g %d %d %d %d\n",ix, isize[g+1], isize[0], g);*/
            for (j = 0; j < natoms; j++)
            {
                bExcl[j] = FALSE;
            }
            k = 0;
            snew(pairs[g][i], isize[g+1]);
            bNonSelfExcl = FALSE;
            for (j = 0; j < isize[g+1]; j++)
            {
                jx = molindex[g+1][j];
                /*fprintf(stderr,"ix jx %d %d\n",ix,jx);*/
                if (!bExcl[jx])
                {
                    pairs[g][i][k++] = jx;
                }
                else if (ix != jx)
                {
                    bNonSelfExcl = TRUE;
                }
            }
            if (bNonSelfExcl)
            {
                npairs[g][i] = k;
                srenew(pairs[g][i], npairs[g][i]);
            }
            else
            {
                npairs[g][i] = -1;
                sfree(pairs[g][i]);
            }
        }
    }
    sfree(bExcl);

    snew(x_i1, max_i);
    nframes    = 0;
    invvol_sum = 0;
    if (bPBC && (NULL != top))
    {
        gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms);
    }
    if (method[0] == 'c')
    {
        do
        {
            /* Must init pbc every step because of pressure coupling */
            copy_mat(box, box_pbc);
            if (bPBC)
            {
                if (top != NULL)
                {
                    gmx_rmpbc(gpbc, natoms, box, x);
                }
                set_pbc(&pbc, ePBCrdf, box_pbc);
    
            }
            invvol      = 1/det(box_pbc);
            invvol_sum += invvol;
    
            for (g = 0; g < ng; g++)
            {
                for (i = 0; i < isize[g+1]; i++)
                {
                    int ind0 ;
                    copy_rvec(x[molindex[g+1][i]], x_i1[i]);
                    ind0  = mols->index[molindex[g+1][i]];
                    mol_unitvectors(x[ind0], x[ind0+1], x[ind0+2], uv1, uv2, uv3);
                    rotate_beta(x[ind0], x[ind0+1], x[ind0+2], arr_qvec ,/*const real beta_m[DIM][DIM][DIM],*/  beta_lab);
                    gmx_fatal(FARGS,"beta_lab %f\n",beta_lab); 
                }
                for (i = 0; i < isize0; i++)
                {
                    if (bClose)
                    {
                        /* Special loop, since we need to determine the minimum distance
                         * over all selected atoms in the reference molecule/residue. */
                        isize_g = isize[g+1];
                        for (j = 0; j < isize_g; j++)
                        {
                            r2 = 1e30;
                            /* Loop over the selected atoms in the reference molecule */
                            for (ii = coi[0][i]; ii < coi[0][i+1]; ii++)
                            {
                                if (bPBC)
                                {
                                    pbc_dx(&pbc, x[molindex[0][ii]], x_i1[j], dx);
                                }
                                else
                                {
                                    rvec_sub(x[molindex[0][ii]], x_i1[j], dx);
                                }
                                r2ii = iprod(dx, dx);
                                if (r2ii < r2)
                                {
                                    r2 = r2ii;
                                }
                            }
                            if (r2 > 0.0 && r2 <= rmax2)
                            {
                                count[g][(int)(sqrt(r2)*invhbinw)]++;
                            }
                        }
                    }
                    else
                    {
                        /* Real rdf between points in space */
                        copy_rvec(x[molindex[0][i]], xi);
                        if ( npairs[g][i] >= 0)
                        {
                            /* Expensive loop, because of indexing */
                            n_j=0.0;
                            snew(temp_method,nbinq);
                            for (j = 0; j < npairs[g][i]; j++)
                            {
                                jx = pairs[g][i][j];
                                if (bPBC)
                                {
                                    pbc_dx(&pbc, xi, x[jx], dx);
                                }
                                else
                                {
                                    rvec_sub(xi, x[jx], dx);
                                }
    
                                r2 = iprod(dx, dx);
                                if (r2 > 0.0 &&  r2 <= rmax2)
                                {
                                    r_dist = sqrt(r2);
                                    count[g][(int)(r_dist*invhbinw)]++;
                                    n_j ++ ;
                                    if ((fade == 0.0) || (r_dist <= fade))
                                    {
                                        mod_f = 1.0 /r_dist ;
                                        for (qq = 0; qq < nbinq; qq++)
                                        {
                                            temp_method[qq] += mod_f*sin(arr_q[qq]*r_dist)  ;
                                        }
                                    }
                                    else
                                    {
                                        mod_f = sqr(cos((r_dist-fade)*inv_width))/r_dist ;
                                        for (qq = 0; qq < nbinq; qq++)
                                        {
                                            temp_method[qq] += mod_f*sin(arr_q[qq]*r_dist)  ;
                                        }                                    
                                    }
                                }
                            }
                            for (qq = 0; qq < nbinq; qq++)
                            {
                                 s_method[g][qq] += temp_method[qq]/arr_q[qq] - analytical_integral[qq]*invvol ;
                            }
                            sfree(temp_method);
                        }
                        else
                        {
                            /* Cheaper loop, no exclusions */
                            isize_g = isize[g+1];
                            n_j = 0.0;
                            snew(temp_method,nbinq);
                            for (j = 0; j < isize_g; j++)
                            {
                                if (bPBC)
                                {
                                    pbc_dx(&pbc, xi, x_i1[j], dx);
                                }
                                else
                                {
                                    rvec_sub(xi, x_i1[j], dx);
                                }
                                r2 = iprod(dx, dx);
                                if (r2 > 0.0 && r2 <= rmax2)
                                {   
                                    r_dist = sqrt(r2);
                                    count[g][(int)(r_dist*invhbinw)]++;
                                    n_j ++ ;
                                    if ((fade == 0.0) || (r_dist <= fade))
                                    {
                                        mod_f = 1.0 /r_dist ;
                                        for (qq = 0; qq < nbinq; qq++)
                                        {
                                            temp_method[qq] += mod_f*sin(arr_q[qq]*r_dist)  ;
                                        }
                                    }
                                    else
                                    {
                                        mod_f = sqr(cos((r_dist-fade)*inv_width))/r_dist ;
                                        for (qq = 0; qq < nbinq; qq++)
                                        {
                                            temp_method[qq] += mod_f*sin(arr_q[qq]*r_dist)  ;
                                        }
                                    }
                                }
                            }
                            for (qq = 0; qq < nbinq; qq++)
                            {
                                 s_method[g][qq] += temp_method[qq]/arr_q[qq] - analytical_integral[qq]*invvol ;
                            }
                            sfree(temp_method) ;
                        }
                    }
                }
            }
            nframes++;
        }
        while (read_next_x(oenv, status, &t, x, box));
    }
    else if (method[0] == 's')
    {   
        fprintf(stderr,"loop with sumexp method \n");
        do
        {
            /* Must init pbc every step because of pressure coupling */
            copy_mat(box, box_pbc);
            if (bPBC)
            {
                if (top != NULL)
                {
                    gmx_rmpbc(gpbc, natoms, box, x);
                }
                set_pbc(&pbc, ePBCrdf, box_pbc);
    
            }
            invvol      = 1/det(box_pbc);
            invvol_sum += invvol;
            for (g = 0; g < ng; g++)
            {
                snew(cos_q,nbinq);
                snew(sin_q,nbinq);
                for (i = 0; i < isize0; i++)
                {
                    copy_rvec(x[molindex[0][i]], xi);
                    /*isize_g = isize[g+1];*/
                    for (qq = 0; qq < nbinq; qq++)
                    {
                        q_xi=(minq+dq*qq)*iprod(arr_qvec,xi);
                        cos_q[qq] += cos(q_xi);
                        sin_q[qq] += sin(q_xi);
                    }
                }
                for (qq = 0; qq < nbinq; qq++)
                {
                    s_method[g][qq] += sqr(cos_q[qq]) + sqr(sin_q[qq]);
                }
                sfree(cos_q);
                sfree(sin_q);
            }
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
    if (method[0]=='c')
    {
       /* Calculate volume of sphere segments or length of circle segments */
       snew(inv_segvol, (nbin+1)/2);
       prev_spherevol = 0;
       for (i = 0; (i < (nbin+1)/2); i++)
       {
           r = (i + 0.5)*binwidth;
           spherevol = (4.0/3.0)*M_PI*r*r*r;
           segvol         = spherevol-prev_spherevol;
           inv_segvol[i]  = 1.0/segvol;
           prev_spherevol = spherevol;
       }
   
       snew(rdf, ng);
       for (g = 0; g < ng; g++)
       {
           /* We have to normalize by dividing by the number of frames */
           normfac = 1.0/(nframes*invvol*isize0*isize[g+1]);
           /* Do the normalization */
           nrdf = max((nbin+1)/2, 1+2*faderdf/binwidth);
           snew(rdf[g], nrdf);
           for (i = 0; i < (nbin+1)/2; i++)
           {
               r = i*binwidth;
               if (i == 0)
               {
                   j = count[g][0];
               }
               else
               {
                   j = count[g][i*2-1] + count[g][i*2];
               }
               if ( (faderdf > 0) && (r >= faderdf) )
               {
                   rdf[g][i] = 1 + (j*inv_segvol[i]*normfac-1)*exp(-16*sqr(r/faderdf-1));
               }
               else
               {
                   if (bNormalize)
                   {
                       rdf[g][i] = j*inv_segvol[i]*normfac;
                   }
                   else
                   {
                       rdf[g][i] = j/(binwidth*isize0*nframes);
                   }
               }
           }
           for (; (i < nrdf); i++)
           {
               rdf[g][i] = 1.0;
           }
       }
       for (g = 0; g < ng; g++)
       {
           /* compute analytical integral for the cosmo method and the S(q), and also the S(q) for the conventional g(r) method*/
           for (qq = 0; qq < nbinq ; qq++)
           {
               s_method[g][qq] = 1.0 + s_method[g][qq]/(nframes*isize0) ;
               for (i = 0; i< (nbin+1)/2 ; i++)
               {
                   r = i*binwidth;
                   if ((fade == 0) || (r <= fade ))
                   {
                       s_method_g_r[g][qq] += binwidth*r*sin(arr_q[qq]*r)*(rdf[g][i]-1.0)/arr_q[qq] ;
                   }
                   else
                   {
                       mod_f = sqr(cos((r - fade)*inv_width)) ; 
                       s_method_g_r[g][qq] += mod_f*binwidth*r*sin(arr_q[qq]*r)*(rdf[g][i]-1.0)/arr_q[qq] ;
                   }
               }               
               s_method_g_r[g][qq] = s_method_g_r[g][qq]*4.0*M_PI*isize0*invvol + 1.0;
           }
       }
    }
    else if (method[0]=='s')
    {
       for (g = 0; g < ng; g++)
       {
           for (qq = 0; qq < nbinq ; qq++)
           {
               s_method[g][qq] = s_method[g][qq]/(nframes*isize0)  ;
           }
       }
    }

    sprintf(gtitle, "Structure factor");
    fp = xvgropen(fnSFACT, gtitle, "q", "S(q)", oenv);
    sprintf(refgt, "%s", "");
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
        fprintf(fp, "%10g", arr_q[qq]);
        for (g = 0; g < ng; g++)
        {
            fprintf(fp, " %10g", s_method[g][qq]);
        }
        fprintf(fp, "\n");
    }
    gmx_ffclose(fp);

    do_view(oenv, fnSFACT, NULL);

    if ((fnOSRDF || fnORDF) && method[0]!='s')
    {
        if (fnOSRDF) {fp = xvgropen(fnOSRDF, "S(q) evaluated from g(r)", "q", "S(q)", oenv);}
        if (fnORDF)  {fpn = xvgropen(fnORDF, "Radial distribution function", "r", "g(r)", oenv);}
        if (ng == 1)
        {
            if (output_env_get_print_xvgr_codes(oenv))
            {
                if (fnOSRDF) {fprintf(fp, "@ subtitle \"%s-%s\"\n", grpname[0], grpname[1]);}
                if (fnORDF)  {fprintf(fpn, "@ subtitle \"%s-%s\"\n", grpname[0], grpname[1]);}
            }
        }
        else
        {
            if (output_env_get_print_xvgr_codes(oenv))
            {
                if (fnOSRDF) {fprintf(fp, "@ subtitle \"reference %s\"\n", grpname[0]);}
                if (fnORDF)  {fprintf(fpn, "@ subtitle \"reference %s\"\n", grpname[0]);}
            }
            if (fnOSRDF) {xvgr_legend(fp, ng, (const char**)(grpname+1), oenv);}
            if (fnORDF)  {xvgr_legend(fpn, ng, (const char**)(grpname+1), oenv);}
        }
        if (fnOSRDF)
        {
           for (qq = 0; qq < nbinq  ; qq++)
           {
               fprintf(fp, "%10g", arr_q[qq]);
               for (g = 0; g < ng; g++)
               {
                   fprintf(fp, " %10g", s_method_g_r[g][qq]);
               }
               fprintf(fp, "\n");
           }
           gmx_ffclose(fp);
           do_view(oenv, fnOSRDF, NULL);
        }
        if (fnORDF)
        {
           for (i = 0; (i < nrdf)  ; i++)
           {
               fprintf(fpn, "%10g", i*binwidth);
               for (g = 0; g < ng; g++)
               {
                   fprintf(fpn, " %10g", rdf[g][i]);
               }
               fprintf(fpn, "\n");
           }
           gmx_ffclose(fpn);
           do_view(oenv, fnORDF, NULL);
        }
    }
    if (method[0] == 'c')
    {
       for (g = 0; g < ng; g++)
       {
          sfree(rdf[g]);
          sfree(s_method[g]);
          sfree(s_method_g_r[g]);
       }
       sfree(rdf);
       sfree(s_method);
       sfree(analytical_integral);
       sfree(s_method_g_r);
       sfree(arr_q);
       sfree(cos_q);
       sfree(sin_q);
    }
    else if (method[0] == 's')  
    {
       for (g = 0; g < ng; g++)
       {
          sfree(s_method[g]);
       }
       sfree(s_method);
       sfree(analytical_integral);
       sfree(arr_q);
    }
}

void mol_unitvectors(const rvec xv1, const rvec xv2, const rvec xv3, rvec u1, rvec u2, rvec u3)
{
    rvec_sub( xv2, xv3, u2);
    unitv(u2, u2);
    svmul(2.0, xv1, u1);
    rvec_sub( u1, xv2, u1);
    rvec_sub( u1, xv2, u1);
    unitv(u1, u1);
    cprod(u1 , u2, u3);
    
}

void rotate_beta(const rvec xv1, const rvec xv2, const rvec xv3, const rvec qv, /*const real beta_m[DIM][DIM][DIM],*/ real beta_l)
{
    int i, k, q;
    matrix um;
    real c1;
    rvec c2;
    clear_rvec(c2);
    rvec_sub( xv2, xv3, um[1]);
    unitv(um[1], um[1]);
    svmul(2.0, xv1, um[0]);
    rvec_sub( um[0], xv2, um[0]);
    rvec_sub( um[0], xv2, um[0]);
    unitv(um[0], um[0]);
    cprod(um[0] , um[1], um[2]);
    
   for ( i = 0; i < DIM; i++)
   {
       c1 = iprod(um[i], qv);
       for ( k = 0; k < DIM; k++)
       {
           c2[k] = iprod(um[k], qv);
           for ( q = 0; q < DIM; q++)
           {
               fprintf(stderr,"%f %f %f\n", um[i][XX], um[k][XX], um[q][XX]);
               fprintf(stderr,"%f %f %f\n", c1, c2[k], iprod(um[q],qv));
               beta_l += c1*c2[k]*iprod(um[q],qv)*2.0;
           }
       }
   }

}

/*void rotate_beta(const rvec uv[], const rvec qv, beta_m, beta_l)
{
   int  i, k, q;
   real c1;
   rvec c2;
   clear_rvec(beta_m);
   for ( i = 0; i < DIM; i++)
   {
       c1 = iprod(uv[i], qv);
       for ( j = 0; j < DIM; j++)
       {
           c2[j] = iprod(xv2, qv);
           for ( k = 0; k < DIM; k++)
           {
               beta_l += c1*c2[j]*iprod(xv3,qv)
           }
       }
   }

    iprod 
}*/

void dipole_atom2molindex(int *n, int *index, t_block *mols)
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
        /* Modify the index in place */
        index[nmol++] = m;
        /*fprintf(stderr, "index nmol m %d %d %d\n", index[nmol], nmol, m);*/
    }
    printf("There are %d molecules in the selection\n", nmol);
    *n = nmol;
}

int gmx_nonlinearopticalscattering(int argc, char *argv[])
{
    const char        *desc[] = {
        "The structure of liquids can be studied by either neutron or X-ray scattering",
        "[THISMODULE] calculates the non-linear optical scattering intensity per molecule in different ways.",
        "The simplest method (sumexp) is to use 1/N|sum_i beta_IJK(i) exp[iqr_i]|^2.",
        "This however converges slowly with the box size for small values of q.[PAR]",
        "The method cosmo (default) uses the following expression to compute I(q):",
        "I(q)=1+ 1/N<sum_{ij} sin(q r_{ij})/(q r_{ij})> -4pi rho int r sin q r dr (to edit).",
    };
    static gmx_bool    bPBC = TRUE, bNormalize = TRUE;
    static real        binwidth = 0.002, maxq=100.0, minq=2.0*M_PI/1000.0, fade = 0.0, faderdf = 0.0;
    static real        kx = 1, ky = 0, kz = 0;
    static int         ngroups = 1, nbinq = 100;

    static const char *methodt[] = { NULL, "cosmo",  "sumexp",  NULL }; 

    t_pargs            pa[] = {
        { "-maxq",      FALSE, etREAL, {&maxq},
        "max wave-vector (1/nm)" },
        { "-minq",      FALSE, etREAL, {&minq},
        "min wave-vector (1/nm)" },
        { "-nbinq",      FALSE, etINT, {&nbinq},
        "number of bins over wave-vector" },
        { "-qx",         FALSE, etREAL, {&kx}, "direction of q-vector in x (1 or 0), use with sumexp method" },
        { "-qy",         FALSE, etREAL, {&ky}, "direction of q-vector in y (1 or 0), use with sumexp method"},
        { "-qz",         FALSE, etREAL, {&kz}, "direction of q-vector in z (1 or 0), use with sumexp method" },
        { "-bin",      FALSE, etREAL, {&binwidth},
          "Binwidth for g(r) (nm)" },
        { "-method",     FALSE, etENUM, {methodt},
          "I(q) using the different methods" },
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances. Without PBC the maximum range will be three times the largest box edge." },
        { "-norm",     FALSE, etBOOL, {&bNormalize},
          "Normalize for volume and density" },
        { "-ng",       FALSE, etINT, {&ngroups},
          "Number of secondary groups to compute RDFs around a central group" },
        { "-fade",     FALSE, etREAL, {&fade},
          "In the cosmo method the modification function cos^2((rij-fade)*pi/(2*(L/2-fade))) is used in the fourier transform."
          " If fade is 0.0 nothing is done." },
        { "-faderdf",     FALSE, etREAL, {&faderdf},
          "From this distance onwards the RDF is tranformed by g'(r) = 1 + [g(r)-1] exp(-(r/faderdf-1)^2 to make it go to 1 smoothly. "
          " If faderdf is 0.0 nothing is done." },

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
        { efXVG, "-o",  "sfact",    ffWRITE },
        { efXVG, "-osrdf", "sfact_rdf", ffOPTWR },
        { efXVG, "-ordf", "rdf", ffOPTWR },
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

    dipole_atom2molindex(&gnx[0], grpindex[0], &(top->mols));
   
    do_nonlinearopticalscattering(top, /*fnNDX,*/ /*fnTPS,*/ ftp2fn(efTRX, NFILE, fnm),
           opt2fn("-o", NFILE, fnm), opt2fn_null("-osrdf", NFILE, fnm),
           opt2fn_null("-ordf", NFILE, fnm),
           methodt[0],  bPBC, bNormalize,  maxq, minq, nbinq, kx, ky, kz, binwidth,
           fade, faderdf, gnx, grpindex, grpname, ngroups, oenv);

    return 0;
}

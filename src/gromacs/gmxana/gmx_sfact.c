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
/*#include "sfactor_func.h"*/
#include "names.h"

#include "gromacs/legacyheaders/gmx_fatal.h"

static void do_sfact(const char *fnNDX, const char *fnTPS, const char *fnTRX,
                   const char *fnSFACT, const char *fnOSRDF, const char *fnORDF, /*const char *fnHQ, */
                   const char *method,
                   gmx_bool bPBC, gmx_bool bNormalize,
                   real cutoff, real maxq, real minq, int nbinq, real kx, real ky, real kz, real binwidth, real fade, int ng,
                   const output_env_t oenv)
{
    FILE          *fp;
    FILE          *fpn;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, i, ii, j, k, nbin, qq, j0, j1, n, nframes;
    int          **count;
    real         **s_method, **s_method_g_r, analytical_integral, *arr_q, *cos_q, *sin_q;
    char         **grpname;
    int           *isize, isize_cm = 0, nrdf = 0, max_i, isize0, isize_g;
    atom_id      **index, *index_cm = NULL;
    gmx_int64_t   *sum;
    real           t, rmax2, rmax, cut2, r, r_dist, r2, r2ii, q_xi, dq, invhbinw, normfac;
    real           segvol, spherevol, prev_spherevol, **rdf;
    rvec          *x, dx, *x0 = NULL, *x_i1, xi, arr_qvec ;
    real          *inv_segvol, invvol, invvol_sum, rho;
    gmx_bool       bClose, *bExcl, bTop, bNonSelfExcl;
    matrix         box, box_pbc;
    int          **npairs;
    atom_id        ix, jx, ***pairs;
    t_topology    *top  = NULL;
    int            ePBC = -1, ePBCrdf = -1;
    t_block       *mols = NULL;
    t_blocka      *excl;
    t_atom        *atom = NULL;
    t_pbc          pbc;
    gmx_rmpbc_t    gpbc = NULL;
    int           *is   = NULL, **coi = NULL, cur, mol, i1, res, a;

    excl = NULL;

    bClose = FALSE ; /*(method[0] != 'c'); */
    if (fnTPS)
    {
        snew(top, 1);
        bTop = read_tps_conf(fnTPS, title, top, &ePBC, &x, NULL, box, TRUE);
        if (bTop )
        {
            /* get exclusions from topology */
            excl = &(top->excls);
        }
    }
    snew(grpname, ng+1);
    snew(isize, ng+1);
    snew(index, ng+1);
    fprintf(stderr, "\nSelect a reference group and %d group%s\n",
            ng, ng == 1 ? "" : "s");
    if (fnTPS)
    {
        get_index(&(top->atoms), fnNDX, ng+1, isize, index, grpname);
        atom = top->atoms.atom;
    }
    else
    {
        rd_index(fnNDX, ng+1, isize, index, grpname);

    }

    if (bClose )
    {
        isize0 = is[0];
        snew(x0, isize0);
    }
    else
    {
        isize0 = isize[0];
    }
    natoms = read_first_x(oenv, &status, fnTRX, &t, &x, box);
    if (!natoms)
    {
        gmx_fatal(FARGS, "Could not read coordinates from statusfile\n");
    }
    if (fnTPS)
    {
        /* check with topology */
        if (natoms > top->atoms.nr)
        {
            gmx_fatal(FARGS, "Trajectory (%d atoms) does not match topology (%d atoms)",
                      natoms, top->atoms.nr);
        }
    }
    /* check with index groups */
    for (i = 0; i < ng+1; i++)
    {   
        for (j = 0; j < isize[i]; j++)
        {
            if (index[i][j] >= natoms)
            {
                gmx_fatal(FARGS, "Atom index (%d) in index group %s (%d atoms) larger "
                          "than number of atoms in trajectory (%d atoms)",
                          index[i][j], grpname[i], isize[i], natoms);
            }
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
        rmax2   = 0.99*0.99*max_cutoff2(FALSE ? epbcXY : epbcXYZ, box_pbc);
        fprintf(stderr,"rmax2 %f\n",rmax2);
        fprintf(stderr,"box_pbc %f\n",box_pbc[XX][XX]);
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
    cut2     = sqr(cutoff);
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
        snew(cos_q,nbinq);
        snew(sin_q,nbinq);
        normfac = 1.0/sqrt(kx*kx + ky*ky + kz*kz) ;
        arr_qvec[XX] = kx*normfac;
        arr_qvec[YY] = ky*normfac;
        arr_qvec[ZZ] = kz*normfac;
        dq=(maxq-minq)/nbinq ;       
        for (qq = 0; qq< nbinq; qq++)
        {
            arr_q[qq]=minq+dq*qq;
        }
        for (i = 0; i < isize[0]; i++)
        {   
            /* We can only have exclusions with atomic rdfs */
            ix = index[0][i];
            for (j = 0; j < natoms; j++)
            {
                bExcl[j] = FALSE;
            }
            if (excl)
            {
               for (j = excl->index[ix]; j < excl->index[ix+1]; j++)
               {
                   bExcl[excl->a[j]] = TRUE;
               }
            }
            k = 0;
            snew(pairs[g][i], isize[g+1]);
            bNonSelfExcl = FALSE;
            for (j = 0; j < isize[g+1]; j++)
            {
                jx = index[g+1][j];
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
                    copy_rvec(x[index[g+1][i]], x_i1[i]);
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
                                    pbc_dx(&pbc, x[index[0][ii]], x_i1[j], dx);
                                }
                                else
                                {
                                    rvec_sub(x[index[0][ii]], x_i1[j], dx);
                                }
                                r2ii = iprod(dx, dx);
                                if (r2ii < r2)
                                {
                                    r2 = r2ii;
                                }
                            }
                            if (r2 > cut2 && r2 <= rmax2)
                            {
                                count[g][(int)(sqrt(r2)*invhbinw)]++;
                            }
                        }
                    }
                    else
                    {
                        /* Real rdf between points in space */
                        copy_rvec(x[index[0][i]], xi);
                        if ( npairs[g][i] >= 0)
                        {
                            /* Expensive loop, because of indexing */
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
                                if (r2 > cut2 && r2 <= rmax2)
                                {
                                    r_dist = sqrt(r2);
                                    count[g][(int)(r_dist*invhbinw)]++;
                                    for (qq = 0; qq < nbinq; qq++)
                                    {
                                        s_method[g][qq]+=sin(arr_q[qq]*r_dist)/(arr_q[qq]*r_dist);
                                    }
    
                                }
                            }
                        }
                        else
                        {
                            /* Cheaper loop, no exclusions */
                            isize_g = isize[g+1];
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
                                if (r2 > cut2 && r2 <= rmax2)
                                {   
                                    r_dist = sqrt(r2);
                                    count[g][(int)(r_dist*invhbinw)]++;
                                    for (qq = 0; qq < nbinq; qq++)
                                    {
                                        s_method[g][qq]+=sin(arr_q[qq]*r_dist)/(arr_q[qq]*r_dist);
                                    }
    
                                }
                            }
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
                for (i = 0; i < isize[g+1]; i++)
                {
                    copy_rvec(x[index[g+1][i]], x_i1[i]);
                }
                snew(cos_q,nbinq);
                snew(sin_q,nbinq);
                for (i = 0; i < isize0; i++)
                {
                    copy_rvec(x[index[0][i]], xi);
                    isize_g = isize[g+1];
                    for (qq = 0; qq < nbinq; qq++)
                    {
                        q_xi=(minq+dq*qq)*iprod(arr_qvec,xi);
                        cos_q[qq] += cos(q_xi);
                        sin_q[qq] += sin(q_xi);
                    }
                }
                for (qq = 0; qq < nbinq; qq++)
                {
                    s_method[g][qq]+= cos_q[qq]*cos_q[qq] + sin_q[qq]*sin_q[qq];
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
           nrdf = max((nbin+1)/2, 1+2*fade/binwidth);
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
               if ((fade > 0) && (r >= fade))
               {
                   rdf[g][i] = 1 + (j*inv_segvol[i]*normfac-1)*exp(-16*sqr(r/fade-1));
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
               analytical_integral=(sin(arr_q[qq]*rmax) - arr_q[qq]*rmax*cos(arr_q[qq]*rmax))/(arr_q[qq]*arr_q[qq]*arr_q[qq]);
               s_method[g][qq] = 1.0 + s_method[g][qq]/(nframes*isize0) - 4.0*M_PI*isize0*invvol*analytical_integral;
               for (i = 0; i< (nbin+1)/2 ; i++)
               {
                   r = i*binwidth;
                   s_method_g_r[g][qq] += binwidth*r*sin(arr_q[qq]*r)*(rdf[g][i]-1.0)/arr_q[qq] ;
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

    if (fnOSRDF && method[0]!='s')
    {
        fp = xvgropen(fnOSRDF, "S(q) evaluated from g(r)", "q", "S(q)", oenv);
        fpn = xvgropen(fnORDF, "Radial distribution function", "r", "g(r)", oenv);
        if (ng == 1)
        {
            if (output_env_get_print_xvgr_codes(oenv))
            {
                fprintf(fp, "@ subtitle \"%s-%s\"\n", grpname[0], grpname[1]);
                fprintf(fpn, "@ subtitle \"%s-%s\"\n", grpname[0], grpname[1]);
            }
        }
        else
        {
            if (output_env_get_print_xvgr_codes(oenv))
            {
                fprintf(fp, "@ subtitle \"reference %s\"\n", grpname[0]);
                fprintf(fpn, "@ subtitle \"reference %s\"\n", grpname[0]);
            }
            xvgr_legend(fp, ng, (const char**)(grpname+1), oenv);
            xvgr_legend(fpn, ng, (const char**)(grpname+1), oenv);
        }
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
    /*
    if (fnHQ)
    {
        int   nhq = 401;
        real *hq, *integrand, Q;

        rho = isize[1]*invvol;
        snew(hq, nhq);
        snew(integrand, nrdf);
        for (i = 0; (i < nhq); i++)
        {
            Q            = i*0.5;
            integrand[0] = 0;
            for (j = 1; (j < nrdf); j++)
            {
                r             = j*binwidth;
                integrand[j]  = (Q == 0) ? 1.0 : sin(Q*r)/(Q*r);
                integrand[j] *= 4.0*M_PI*rho*r*r*(rdf[0][j]-1.0);
            }
            hq[i] = print_and_integrate(debug, nrdf, binwidth, integrand, NULL, 0);
        }
        fp = xvgropen(fnHQ, "h(Q)", "Q(/nm)", "h(Q)", oenv);
        for (i = 0; (i < nhq); i++)
        {
            fprintf(fp, "%10g %10g\n", i*0.5, hq[i]);
        }
        gmx_ffclose(fp);
        do_view(oenv, fnHQ, NULL);
        sfree(hq);
        sfree(integrand);
    }
    */
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
       sfree(s_method_g_r);
       sfree(arr_q);
    }
    else if (method[0] == 's')  
    {
       for (g = 0; g < ng; g++)
       {
          sfree(s_method[g]);
       }
       sfree(s_method);
       sfree(arr_q);
    }
}


int gmx_sfact(int argc, char *argv[])
{
    const char        *desc[] = {
        "The structure of liquids can be studied by either neutron or X-ray scattering",
        "[THISMODULE] calculates the structure factor S(q) in different ways.",
        "The simplest method (sumexp) is to use 1/N|sum_j exp(iqrj)|^2.",
        "This however converges slowly with the box size for small values of q.[PAR]",
        "The option rdf is based on the definition of S(q) for an isotropic system c.f",
        "M.P. Allen and D.J. Tildesley pp. 58.[PAR]",
        "The method cosmo (default) uses the following expression to compute S(q):",
        "S(q)=1+ 1/N<sum_{ij} sin(qr_{ij})/(qr_{ij})> -4pi rho int r sin qr dr.[PAR]",
    };
    static gmx_bool    /*bCM     = FALSE,*/ bPBC = TRUE, bNormalize = TRUE;
    static real        cutoff  = 0, binwidth = 0.002, maxq=100.0, minq=2.0*M_PI/1000.0, fade = 0.0;
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
        { "-kx",         FALSE, etREAL, {&kx}, "direction of k-vector in x (1 or 0)" },
        { "-ky",         FALSE, etREAL, {&ky}, "direction of k-vector in y (1 or 0)"},
        { "-kz",         FALSE, etREAL, {&kz}, "direction of k-vector in z (1 or 0)" },
        { "-bin",      FALSE, etREAL, {&binwidth},
          "Binwidth for g(r) (nm)" },
/*        { "-com",      FALSE, etBOOL, {&bCM},
          "RDF with respect to the center of mass of first group" }, */
        { "-method",     FALSE, etENUM, {methodt},
          "S(q) using the different methods" },
/*        { "-rdf",   FALSE, etENUM, {rdft},
          "RDF type" },
*/
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances. Without PBC the maximum range will be three times the largest box edge." },
        { "-norm",     FALSE, etBOOL, {&bNormalize},
          "Normalize for volume and density" },
        { "-cut",      FALSE, etREAL, {&cutoff},
          "Shortest distance (nm) to be considered"},
        { "-ng",       FALSE, etINT, {&ngroups},
          "Number of secondary groups to compute RDFs around a central group" },
        { "-fade",     FALSE, etREAL, {&fade},
          "From this distance onwards the RDF is tranformed by g'(r) = 1 + [g(r)-1] exp(-(r/fade-1)^2 to make it go to 1 smoothly. If fade is 0.0 nothing is done." }
    };
#define NPA asize(pa)
    const char        *fnTPS, *fnNDX;
    output_env_t       oenv;

    t_filenm           fnm[] = {
        { efTRX, "-f",  NULL,     ffREAD },
        { efTPS, NULL,  NULL,     ffOPTRD },
        { efNDX, NULL,  NULL,     ffOPTRD },
        { efXVG, "-o",  "sfact",    ffWRITE },
        { efXVG, "-osrdf", "sfact_rdf", ffOPTWR },
        { efXVG, "-ordf", "rdf", ffOPTWR },
       /* { efXVG, "-hq", "hq",     ffOPTWR },*/
    };
#define NFILE asize(fnm)
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
   
    do_sfact(fnNDX, fnTPS, ftp2fn(efTRX, NFILE, fnm),
           opt2fn("-o", NFILE, fnm), opt2fn_null("-osrdf", NFILE, fnm),
           opt2fn_null("-ordf", NFILE, fnm),
           /*opt2fn_null("-hq", NFILE, fnm),*/
           /*bCM,*/ methodt[0],  bPBC, bNormalize, cutoff, maxq, minq, nbinq, kx, ky, kz, binwidth, fade, ngroups,
           oenv);

    return 0;
}

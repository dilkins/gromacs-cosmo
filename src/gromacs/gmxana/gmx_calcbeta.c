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




static void do_calcbeta(t_topology *top,  const char *fnTRX,
                   const char *fnBETA,
                   const char *fnVCOEFF, const char *fnVGRD,
                   const char *fnRGRDO, const char *fnCOEFFO,
                   const char *fnRGRDH, const char *fnCOEFFH, const char *fnREFMOL,
                   const char *method, const char *kern,
                   gmx_bool bPBC, int kern_order, real std_dev_dens, real fspacing,
                   real kappa, int interp_order, int kmax,
                   int *isize, int  *molindex[], char **grpname, int ng,
                   const output_env_t oenv, real eps)
{
    FILE          *fBETA;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, nanions, ncations, i, j, k, qq, n, c, tt, rr, nframes_bla[1], nframes, gr_ind, nbin, aa, bb, cc;
    real         **field_ad, electrostatic_cutoff2, electrostatic_cutoff, max_spacing, maxelcut2,  invkappa2,  ***beta_mol, *betamean ,*mu_ind_mols, ****mu_ind_t, *beta_corr, *ft_beta_corr;
    int            max_i, isize0, ind0, indj;
    real           t, rmax2, rmax,  r, r_dist, r2, q_xi, dq;
    real          *inv_segvol, normfac, segvol, spherevol, prev_spherevol, invsize0, invgamma, invhbinw, inv_width,  theta=0, *theta_vec;
    rvec          *x, xcm, xcm_transl, dx,  *x_i1, xi, x01, x02, *arr_qvec, **arr_qvec_faces ,vec_polin, vec_polout, ***vec_pout_theta_gamma, ***vec_pin_theta_gamma;
    rvec           xvec, yvec, zvec, *xmol, *xref, Emean;
    real          *qref;
    matrix         cosdirmat,invcosdirmat; 
    real           invvol, invvol_sum;
    t_Map         *Map=NULL;
    t_Kern        *Krr = NULL;
    t_Kern        *SKern_rho_O = NULL;
    t_Kern        *SKern_rho_H = NULL;
    t_Kern        *SKern_E = NULL;
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
    fprintf(stderr,"open binary file to write beta_mol for the trajectory\n");
    fBETA=fopen(fnBETA,"wb");
    atom = top->atoms.atom;
    mols = &(top->mols);
    isize0 = isize[0];
//    isize0 = 2;
    molsize = mols->index[molindex[0][1]] - mols->index[molindex[0][0]];
    snew(xmol,molsize);
    snew(xref,molsize);
    snew(qref,molsize);
    invsize0 = 1.0/isize0;
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
    if (kern[0] == 's')
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

    if (bPBC)
    {
        rmax2   =  max_cutoff2(FALSE ? epbcXY : epbcXYZ, box_pbc);

        fprintf(stderr, "rmax2 = %f\n", rmax2);
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
           fprintf(stderr,"max wave-vector for ewald sum is %d\n",kmax);
           setup_ewald_pair_potential(gridsize,interp_order,kmax,&FT_pair_pot,invkappa2);
           fprintf(stderr,"ewald pair potential has been set-up\n");
        
        }
        
    }

    set_pbc(&pbc, ePBCrdf, box_pbc);
    rmax     = sqrt(rmax2);
    //initialize beta tensor
    snew(beta_mol, DIM);
    snew(mu_ind_mols, isize0);
    
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

    max_i = isize[0];

    /*allocate memory for beta in lab frame and initialize beta in mol frame*/

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

                    rotate_local_grid_points(SKern_rho_O, SKern_rho_H, SKern_E, ePBC, box, cosdirmat, xi );
                    //fprintf(stderr,"finished rotating and translating grid\n");
                    trilinear_interpolation_kern(SKern_rho_O, ir, &pbc, xi, grid_invspacing, grid_spacing);
                    //fprintf(stderr,"finished interpolation O kern\n");
                    trilinear_interpolation_kern(SKern_rho_H, ir, &pbc, xi, grid_invspacing, grid_spacing);
                    //fprintf(stderr,"finished interpolation H kern\n");
                    vec_trilinear_interpolation_kern(SKern_E, ir, &pbc, invcosdirmat, xi,
                                                     grid_invspacing, grid_spacing, Emean);
                    //fprintf(stderr,"finished interpolation E kern\n");

                    calc_beta_skern(SKern_rho_O, SKern_rho_H, SKern_E, kern_order, betamean, &beta_mol);

                    //char abc[10]="abcdefghij";
                    //fwrite(abc,sizeof(abc[10]),sizeof(abc)/sizeof(abc[0]),fBETA);
                    fwrite(&nframes,sizeof(int),1,fBETA);
                    
                    for (aa = 0; aa < DIM; aa++)
                    {
                        for (bb = 0; bb < DIM; bb++)
                        {
                           for (cc = 0; cc < DIM ; cc++)
                           {
                              printf("beta %d %d %d %f\n",aa, bb, cc, beta_mol[aa][bb][cc]);
                              fwrite(&beta_mol[aa][bb][cc],sizeof(float),1,fBETA);
                           }
                        }
                    }

//                    if (debug )
//                    {
//                       gmx_fatal(FARGS,"end\n");
//                    }

                }
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
    fclose(fBETA);
    sfree(x);

    /* Average volume */
    invvol = invvol_sum/nframes;

    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            sfree(beta_mol[i][j]);
        }
        sfree(beta_mol[i]);
    }
    sfree(beta_mol);
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

 
int gmx_calcbeta(int argc, char *argv[])
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
    static gmx_bool          bPBC = TRUE;
    static real              electrostatic_cutoff = 1.2, kappa = 5.0 ;
    static real              fspacing = 0.01, std_dev_dens = 0.05;
    static real              eps = -1.0 ;
    static int               ngroups = 1 ;
    static int               nkx = 0, nky = 0, nkz = 0, kern_order = 2, interp_order = 4, kmax =0;

    static const char *methodt[] = {NULL, "single", "double" ,NULL };
    static const char *kernt[] = {NULL, "krr", "scalar", "none", "map", NULL};
    static char *catname = NULL;
    static char *anname =  NULL;

    t_pargs            pa[] = {
        { "-stddens",       FALSE, etREAL, {&std_dev_dens}, "standard deviation to compute density on a grid. Use only with scalar kernel [nm]."},
        { "-fourierspacing",          FALSE, etREAL, {&fspacing}, "fourier spacing [nm] gives lower bound for number of wave vectors to use in each direction with Ewald, overridden by nkx,nky,nkz " },
        { "-kappa",        FALSE, etREAL, {&kappa}, "screening parameter for the ewald term, i.e. erf(r*kappa)/r, in nm^-1" },
        { "-kernorder",        FALSE, etINT, {&kern_order}, "kernel order, where beta = sum_i sum_(kern_ind) c[i*kern_ind] * (feature_vec[i]-mean(feature_vec))^kern_ind " },     
        { "-splorder",        FALSE, etINT, {&interp_order}, "interpolation order for b-splines" },
        { "-kmax_spme",        FALSE, etINT, {&kmax}, "number of fourier components within cutoff in fourier space for SPME (actual number is 4pi/3*kmax_spme^3)" },
        { "-method",     FALSE, etENUM, {methodt}, "I(q) using a single summation O(N) or a double summation, O(N*N)" },
        { "-kern",   FALSE, etENUM, {kernt}, "what method to use to compute beta"},
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances. Always use, results without PBC not tested." },
        { "-ng",       FALSE, etINT, {&ngroups}, 
          "Number of secondary groups, not available for now. Only tip4p water implemented." },
	{ "-eps",	FALSE, etREAL, {&eps}, "dielectric constant"},
    };
#define NPA asize(pa)
    const char        *fnTPS, *fnNDX , *fnBETA = NULL;
    const char        *fnVCOEFF=NULL, *fnVGRD=NULL, *fnREFMOL=NULL;
    const char        *fnRGRDO=NULL, *fnRGRDH=NULL, *fnCOEFFO=NULL, *fnCOEFFH=NULL;
    output_env_t       oenv;
    int           *gnx;
    int            nFF[2];
    atom_id      **grpindex;
    char         **grpname = NULL;

    t_filenm           fnm[] = {
        { efTRX, "-f",  NULL,     ffREAD },
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
        { efDAT, "-o",  "betafile.dat",    ffWRITE },

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
    fnVCOEFF = opt2fn_null("-vcoeff", NFILE,fnm);
    fnVGRD = opt2fn_null("-vgrid", NFILE,fnm);
    fnRGRDO = opt2fn_null("-rhogridO", NFILE,fnm);
    fnRGRDH = opt2fn_null("-rhogridH", NFILE,fnm);
    fnCOEFFO = opt2fn_null("-rhocoeffO", NFILE,fnm);
    fnCOEFFH = opt2fn_null("-rhocoeffH", NFILE,fnm);
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
    else if ((*kernt)[0] != 's')
    {
       gmx_fatal(FARGS, "use only scalar kernel to compute molecular beta\n");
    }
 
    snew(top, ngroups+1);
    ePBC = read_tpx_top(ftp2fn(efTPS, NFILE, fnm), NULL, box,
                        &natoms, NULL, NULL, NULL, top);

    snew(gnx, ngroups+1);
    snew(grpname, ngroups+1);
    snew(grpindex, ngroups+1);
    get_index(&top->atoms, ftp2fn_null(efNDX, NFILE, fnm),
             ngroups +1 , gnx, grpindex, grpname);

    if (kmax == 0)
    {
      // If no kmax is specified, then one is chosen according to the "optimizing" formula of Frenkel and Smit.
      //kmax = (int)(M_PI/(min(0.5,1.2*pow(natoms,-1.0/6.0))));
      kmax = (int)(2.5*box[XX][XX]*kappa/M_PI);
      fprintf(stderr,"\nNo kmax_spme specified; setting kmax_spme = %i\n",kmax);
    }


    fprintf(stderr,"Start indexing the atoms to each molecule\n");
    dipole_atom2mol(&gnx[0], grpindex[0], &(top->mols));

    do_calcbeta(top, ftp2fn(efTRX, NFILE, fnm),
           opt2fn("-o", NFILE, fnm),
           fnVCOEFF, fnVGRD, fnRGRDO, fnCOEFFO,
           fnRGRDH, fnCOEFFH, fnREFMOL, methodt[0], kernt[0], bPBC, 
           kern_order, std_dev_dens, fspacing, 
           kappa, interp_order, kmax,  gnx, grpindex, grpname, ngroups, oenv, eps);
    return 0;
}

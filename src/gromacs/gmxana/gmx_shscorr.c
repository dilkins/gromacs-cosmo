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

static void do_shscorr(t_topology *top,  const char *fnTRX,
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
                   const output_env_t oenv, int nmax, real intheta)

{
    FILE          *fp, *fpn;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, nanions, ncations, i, j, k, qq, n, c, tt, rr, nframes, nfaces, gr_ind, nbin, aa, bb, cc;
    real         ***s_method, ***s_method_coh, ***s_method_incoh, **temp_method, ****s_method_t, ****s_method_coh_t, ****s_method_incoh_t, ***mu_sq_t, ***coh_temp;
    real           qnorm, maxq, incoh_temp = 0.0, tot_temp = 0.0, gamma = 0.0 ,theta0 = 5.0, check_pol;
    real          *cos_t, *sin_t, ****cos_tq, ****sin_tq,   mu_ind =0.0, mu_sq =0.0, mod_f ;
    real         **field_ad, electrostatic_cutoff2, electrostatic_cutoff, max_spacing, maxelcut2,  invkappa2,  ***beta_mol, *betamean ,*mu_ind_mols, ****mu_ind_t, *beta_corr, *ft_beta_corr;
    int            max_i, isize0, ind0, indj;
    real           t, rmax2, rmax,  r, r_dist, r2, q_xi, dq;
    real          *inv_segvol, normfac, segvol, spherevol, prev_spherevol, invsize0, invgamma, invhbinw, inv_width,  theta=0, *theta_vec;
    rvec          *x, xcm, xcm_transl, dx,  *x_i1, xi, x01, x02, qvec, **arr_qvec, **arr_qvec_faces ,vec_polin, vec_polout, ***vec_pout_theta_gamma, ***vec_pin_theta_gamma;
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
    int		       nx,ny,nz,maxnpoint,**narray,nmax2,n_used,n2,ii,jj,kk,l,ff,*num_count,nn,nxa,nya,nza,nxb,nyb,nzb;
    real	      **kvec,***u_vec,***v_vec,**basis,*coeff,saout,sain,caout,cain,*****beta_lab,st,ct,zt;
	rvec			***all_r;
	real			**intens_total,**intens_cohrt,**intens_incoh,****s_array,****c_array,rval,**snt,**cst;

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
       fprintf(stderr,"the density on the grid points is weighted with gaussians that are centred on the molecule\n");
       fprintf(stderr,"with a standard deviation of %f\n",4.0*std_dev_dens);
       fprintf(stderr,"if you want to change the standard deviation for these weights you have to modify the source\n");
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

	// Read in all frames.
	nframes = 0;
	do{nframes++;} while (read_next_x(oenv,status,&t,x,box));
	read_first_x(oenv, &status, fnTRX, &t, &x, box);

	snew(all_r,nframes);
	for (i=0;i<nframes;i++)
	{
		snew(all_r[i],isize0);
	}
	for (i=0;i<nframes;i++)
	{
		for (j=0;j<isize0;j++)
		{
			snew(all_r[i][j],molsize);
		}
	}

	// We already have the first frame.
	for (ii=0;ii<isize0;ii++)
	{
		ind0 = mols->index[molindex[0][ii]];
		for (j=0;j<molsize;j++)
		{
			copy_rvec(x[ind0+j],all_r[0][ii][j]);
		}
	}

	for (i=1;i<nframes;i++)
	{
		read_next_x(oenv,status,&t,x,box);
		for (ii=0;ii<isize0;ii++)
		{
			ind0 = mols->index[molindex[0][ii]];
			for (j=0;j<molsize;j++)
			{
				copy_rvec(x[ind0+j],all_r[i][ii][j]);
			}
		}
		
	}
//	// DO SOME QUICK TESTING.
//	fprintf(stderr,"%i\n",nframes);
//	fprintf(stderr,"%f %f %f\n",all_r[0][0][0][0],all_r[0][0][0][1],all_r[0][0][0][2]);
//	fprintf(stderr,"%f %f %f\n",all_r[0][0][1][0],all_r[0][0][1][1],all_r[0][0][1][2]);
//	fprintf(stderr,"%f %f %f\n",all_r[0][0][2][0],all_r[0][0][2][1],all_r[0][0][2][2]);
//	fprintf(stderr,"%f %f %f\n",all_r[0][0][3][0],all_r[0][0][3][1],all_r[0][0][3][2]);
//	fprintf(stderr,"%f %f %f\n",all_r[99][0][0][0],all_r[99][0][0][1],all_r[99][0][0][2]);
//	fprintf(stderr,"%f %f %f\n",all_r[99][0][1][0],all_r[99][0][1][1],all_r[99][0][1][2]);
//	fprintf(stderr,"%f %f %f\n",all_r[99][0][2][0],all_r[99][0][2][1],all_r[99][0][2][2]);
//	fprintf(stderr,"%f %f %f\n",all_r[99][0][3][0],all_r[99][0][3][1],all_r[99][0][3][2]);
//	exit(0);

	// Initialize sine and cosine arrays.
	snew(s_array,nframes);
	snew(c_array,nframes);
	for (ff=0;ff<nframes;ff++)
	{
		snew(s_array[ff],isize0);
		snew(c_array[ff],isize0);
	}
	for (ff=0;ff<nframes;ff++)
	{
		for (i=0;i<isize0;i++)
		{
			snew(s_array[ff][i],3);
			snew(c_array[ff][i],3);
		}
	}
	for (ff=0;ff<nframes;ff++)
	{
		for (i=0;i<isize0;i++)
		{
			for (j=0;j<3;j++)
			{
				snew(s_array[ff][i][j],nmax+1);
				snew(c_array[ff][i][j],nmax+1);
			}
		}
	}

	// Now, fill these arrays.
	for (ff=0;ff<nframes;ff++)
	{
		for (i=0;i<isize0;i++)
		{
			for (j=0;j<3;j++)
			{
				rval = all_r[ff][i][0][j];
				s_array[ff][i][j][0] = 0.0;
				c_array[ff][i][j][0] = 1.0;
				q_xi = 2.0 * M_PI * rval / box[j][j];
				// CHECK THIS. AND TAKE OUTSIDE OF LOOPS.
				s_array[ff][i][j][1] = sin(q_xi);
				c_array[ff][i][j][1] = cos(q_xi);
//				fprintf(stderr,"HERE %i %i %i %f %f %f %f\n",ff,i,j,rval,q_xi,sin(q_xi),cos(q_xi));
				for (nn=2;nn<nmax+1;nn++)
				{
					s_array[ff][i][j][nn] = s_array[ff][i][j][nn-1]*c_array[ff][i][j][1] + c_array[ff][i][j][nn-1]*s_array[ff][i][j][1];
					c_array[ff][i][j][nn] = c_array[ff][i][j][nn-1]*c_array[ff][i][j][1] - s_array[ff][i][j][nn-1]*s_array[ff][i][j][1];
				}
			}
		}
	}

//	exit(0);

    //initialize beta tensor

    snew(beta_mol, DIM);
	snew(beta_lab, DIM);
    snew(mu_ind_mols, isize0);
    
    snew(beta_corr, nbin+1);
    snew(ft_beta_corr,nbinq);
    for (i = 0; i < DIM; i++)
    {
        snew(beta_mol[i], DIM);
		snew(beta_lab[i], DIM);
    }
    for (i = 0; i < DIM; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            snew(beta_mol[i][j], DIM);
            snew(beta_lab[i][j], DIM);
        }
    }
	for (i = 0;i < DIM; i++)
	{
		for (j = 0;j < DIM; j++)
		{
			for (k = 0;k < DIM; k++)
			{
				snew(beta_lab[i][j][k],isize0);
			}
		}
	}
	for (i = 0;i < DIM; i++)
	{
		for (j = 0;j < DIM; j++)
		{
			for (k = 0;k < DIM; k++)
			{
				for (l = 0;l < isize0;l++)
				{
					snew(beta_lab[i][j][k][l],nframes);
				}
			}
		}
	}

    if (method[0] == 's')
    {
       qnorm = M_PI*2.0/(rmax*2.0);
    }
    else
    {
      qnorm = M_PI*2.0/(rmax);
    }


	// We want to find all (nx,ny,nz) values such that nx^2 + ny^2 + nz^2 <= nmax^2.
	// Firstly, allocate an array to hold all of these values.
	// maxnpoint is the maximum number of n points we might want. We will need fewer
	// than this.
	maxnpoint = (2*nmax + 1)*(2*nmax + 1)*(2*nmax + 1);
	snew(narray,maxnpoint);
	for (i = 0;i < maxnpoint; i++)
	{
		snew(narray[i],4);
	}
	
	nmax2 = nmax*nmax;
	n_used = 0;
	for (nx = -nmax;nx <= nmax;nx++)
	{
		for (ny = -nmax;ny <= nmax;ny++)
		{
			for (nz = -nmax;nz <= nmax;nz++)
			{
				n2 = nx*nx + ny*ny + nz*nz;
				if ((n2 != 0) && (n2 <= nmax2))
				{
					narray[n_used][0] = nx;
					narray[n_used][1] = ny;
					narray[n_used][2] = nz;
					narray[n_used][3] = n2;
					n_used += 1;
				}
			}
		}
	}
	
	snew(kvec,n_used);
	for (qq=0;qq<n_used;qq++)
	{
		snew(kvec[qq],4);
//		for (j=0;j<4;j++)
//		{
		kvec[qq][0] = qnorm * narray[qq][0];
		kvec[qq][1] = qnorm * narray[qq][1];
		kvec[qq][2] = qnorm * narray[qq][2];
		kvec[qq][3] = qnorm * sqrt(1.0*narray[qq][3]);
//		}
	}
	
	// We now have a number of (nx,ny,nz) vectors; these will determine our scattering
	// wavevectors. For each one, let's go through and work out the vectors u and v.
	snew(u_vec,n_used);
	snew(v_vec,n_used);
	snew(basis,3);
	for (i=0;i<3;i++)
	{
		snew(basis[i],3);
	}
	for (i=0;i<n_used;i++)
	{
		snew(u_vec[i],nbingamma);
		snew(v_vec[i],nbingamma);
		for (j=0;j<nbingamma;j++)
		{
			snew(u_vec[i][j],3);
			snew(v_vec[i][j],3);
		}
	}
	
	saout=sin(M_PI/180.0*pout_angle);
	caout=cos(M_PI/180.0*pout_angle);
	sain=sin(M_PI/180.0*pin_angle);
	cain=cos(M_PI/180.0*pin_angle);
	
	snew(coeff,3);
	intheta = intheta * M_PI/180.0;
	coeff[0] = sin(intheta * 0.5);
	
	for (c=0;c<nbingamma;c++)
	{
		gamma = c*M_PI*invgamma;
		coeff[1] = cos(intheta * 0.5) * sin(gamma);
		coeff[2] = cos(intheta * 0.5) * cos(gamma);
		for (i=0;i<n_used;i++)
		{
			// Given this gamma and this q vector, let's figure out the scattering plane.
			// intheta gives us the scattering angle; that should be all we need!
			// 1. Create an orthonormal basis. The first vector is this kvec element.
			for (j=0;j<3;j++)
			{
				basis[0][j] = kvec[i][j];
			}
			// To find the second vector, we take another vector that's definitely different
			// from the first one, and take the cross product.
			basis[2][0] = basis[0][0] + 0.4;
			basis[2][1] = basis[0][1] + 0.3;
			basis[2][2] = basis[0][2] + 0.1;
			cprod(basis[0],basis[2],basis[1]);
			// The third vector is the cross product of the first with the second.
			cprod(basis[0],basis[1],basis[2]);

			unitv(basis[0],basis[0]);
			unitv(basis[1],basis[1]);
			unitv(basis[2],basis[2]);

			// Now, given this 3-dimensional basis, find k_in and k_out. We take gamma to be the angle made with the third
			// basis vector.
			for (j=0;j<3;j++)
			{
				vec_2kin[j] = coeff[0]*basis[0][j] - coeff[1]*basis[1][j] - coeff[2]*basis[2][j];
				vec_kout[j] = coeff[0]*basis[0][j] + coeff[1]*basis[1][j] + coeff[2]*basis[2][j];
			}
			unitv(vec_2kin,vec_2kin);
			unitv(vec_kout,vec_kout);

			// Polarization vectors for outgoing beam.
			cprod(vec_kout,vec_2kin,pol_perp);
			cprod(vec_kout,pol_perp,pol_par);
			svmul(saout,pol_perp,pol_perp);
			svmul(caout,pol_par,pol_par);
			rvec_add(pol_perp,pol_par,vec_polout);
	
			// Polarization vectors for incoming beam.
			cprod(vec_kout,vec_2kin,pol_perp);
			cprod(vec_2kin,pol_perp,pol_par);
			svmul(sain,pol_perp,pol_perp);
			svmul(cain,pol_par,pol_par);
			rvec_add(pol_perp,pol_par,vec_polin);

			// Normalize everything.
			unitv(vec_polout,vec_polout);
			unitv(vec_polin,vec_polin);

			// Store the u and v vectors.
			for (j=0;j<3;j++)
			{
				u_vec[i][c][j] = vec_polout[j];
				v_vec[i][c][j] = vec_polin[j];
			}

		}
	}
	
	fprintf(stderr,"Using %i points in reciprocal space\n",n_used);

	// We now have a certain number of q vectors, along with the U and V vectors corresponding to
	// them, for a certain number of scattering planes each.

/* USE THESE LINES TO CHECK:

	for (i =0;i<n_used;i++)
	{
		fprintf(stderr,"%i %i %i %i\n",narray[i][0],narray[i][1],narray[i][2],narray[i][3]);
		fprintf(stderr,"%f %f %f %f\n",kvec[i][0],kvec[i][1],kvec[i][2],kvec[i][3]);
		fprintf(stderr,"U1 %f %f %f\n",u_vec[i][0][0],u_vec[i][0][1],u_vec[i][0][2]);
		fprintf(stderr,"U2 %f %f %f\n",u_vec[i][1][0],u_vec[i][1][1],u_vec[i][1][2]);
		fprintf(stderr,"V1 %f %f %f\n",v_vec[i][0][0],v_vec[i][0][1],v_vec[i][0][2]);
		fprintf(stderr,"V2 %f %f %f\n",v_vec[i][1][0],v_vec[i][1][1],v_vec[i][1][2]);
	}

*/

	// Next: for each molecule, rotate the molecular hyperpolarizability tensor into the lab frame (later on,
	// we will multiply by the elements of the polarization vectors).
	for (ii=0;ii<DIM;ii++)
	{
		for (jj=0;jj<DIM;jj++)
		{
			for (kk=0;kk<DIM;kk++)
			{
//				beta_mol[ii][jj][kk] = Map->beta_gas[ii*9 + jj*3 + kk];
				beta_mol[ii][jj][kk] = 0.0;
			}
		}
	}

	beta_mol[2][2][2] = 2.28499883794;
	// DMW: This should be changed later on!

	for (ff=0;ff<nframes;ff++)
	{
		for (i = 0;i<isize0;i++)
		{
			// For this frame and this molecule, calculate the direction cosines. This is the only time
			// we'll need them.

			for (ii=0;ii<molsize;ii++)
			{
				pbc_dx(&pbc,all_r[ff][i][ii],all_r[ff][i][0],xmol[ii]);
			}

			calc_cosdirmat( fnREFMOL, top, molsize, ind0,  xref, xmol, &cosdirmat, &invcosdirmat, &xvec, &yvec, &zvec );

			for (aa=0;aa<DIM;aa++)
			{
				for (bb=0;bb<DIM;bb++)
				{
					for (cc=0;cc<DIM;cc++)
					{
						beta_lab[aa][bb][cc][i][ff] = 0.0;
						for (ii=0;ii<DIM;ii++)
						{
							for (jj=0;jj<DIM;jj++)
							{
								for (kk=0;kk<DIM;kk++)
								{
									beta_lab[aa][bb][cc][i][ff] += beta_mol[ii][jj][kk]*cosdirmat[aa][ii]*cosdirmat[bb][jj]*cosdirmat[cc][kk];
								}
							}
						}
					}
				}
			}
		}
	}

	if (kern[0] != 'n')
	{
		fprintf(stderr,"WARNING: This code has been written for kern[0] = n, so may not work well otherwise!\n");
	}

	snew(snt,nframes);
	snew(cst,nframes);
	for (ff=0;ff<nframes;ff++)
	{
		snew(snt[ff],isize0);
		snew(cst[ff],isize0);
	}

    snew(s_method, ng);
    snew(s_method_coh, ng);
    snew(s_method_incoh, ng);
    snew(s_method_t, ng);
    snew(s_method_coh_t, ng);
    snew(s_method_incoh_t, ng);
    max_i = isize[0];

    for (g = 0; g < ng; g++)
    {
   	    /*allocate memory for s_method array */
   	       snew(s_method[g], n_used);
   	       snew(s_method_coh[g], n_used);
   	       snew(s_method_incoh[g],n_used);
   	       snew(temp_method, n_used);
   	       snew(arr_qvec,n_used);

			for (i=0;i<n_used;i++)
			{
				snew(s_method[g][i],nbingamma);
				snew(s_method_coh[g][i],nbingamma);
				snew(s_method_incoh[g][i],nbingamma);
				snew(temp_method[i],nbingamma);
				snew(arr_qvec[i],nbingamma);
			}

		// For each q vector, let's calculate the intensity.
		for (qq=0;qq<n_used;qq++)
		{
			fprintf(stderr,"Doing q-point number %i of %i\n",qq+1,n_used);
			// To recap (for my own sake: so that I don't forget!):
			// For this group (of which there is only one), we are going
			// to take every different q-vector that we have generated,
			// of which there are n_used, and we will take each value of gamma
			// for each q-vector. This will give us n_used*nbingamma different
			// scattering planes, each with its own u and v vectors (i.e., n_used*nbingamma
			// different experimental setups).
			// For each of these experimental setups, we will calculate the intensity of
			// ESHS light scattered. Each of these setups has a corresponding value of q.

			nx = narray[qq][0];
			ny = narray[qq][1];
			nz = narray[qq][2];
			nxa = abs(nx);
			nya = abs(ny);
			nza = abs(nz);
			// The sine matrices are odd in their arguments; we must take this into account.
			nxb = 2*(nx>0) - 1;
			nyb = 2*(ny>0) - 1;
			nzb = 2*(nz>0) - 1;

			for (ff=0;ff<nframes;ff++)
			{
				for (i=0;i<isize0;i++)
				{
					cst[ff][i] = ( c_array[ff][i][0][nxa]*c_array[ff][i][1][nya]*c_array[ff][i][2][nza] - nxb*nyb*s_array[ff][i][0][nxa]*s_array[ff][i][1][nya]*c_array[ff][i][2][nza] - nxb*nzb*s_array[ff][i][0][nxa]*c_array[ff][i][1][nya]*s_array[ff][i][2][nza] - nyb*nzb*c_array[ff][i][0][nxa]*s_array[ff][i][1][nya]*s_array[ff][i][2][nza]);
					snt[ff][i] = ( nxb*s_array[ff][i][0][nxa]*c_array[ff][i][1][nya]*c_array[ff][i][2][nza] + nyb*c_array[ff][i][0][nxa]*s_array[ff][i][1][nya]*c_array[ff][i][2][nza] + nzb*c_array[ff][i][0][nxa]*c_array[ff][i][1][nya]*s_array[ff][i][2][nza] - nxb*nyb*nzb*s_array[ff][i][0][nxa]*s_array[ff][i][1][nya]*s_array[ff][i][2][nza]);
				}
			}

			for (c=0;c<nbingamma;c++)
			{
				s_method[g][qq][c] = 0.0;
				s_method_coh[g][qq][c] = 0.0;
				s_method_incoh[g][qq][c] = 0.0;

				// For this given group,
				// q-vector,
				// and value of gamma,
				// we are going to calculate the intensity.
				for (ff=0;ff<nframes;ff++)
				{
					st = 0.0;
					ct = 0.0;
					zt = 0.0;
					for (i=0;i<isize0;i++)
					{
						mu_ind = 0.0;
						for (aa=0;aa<DIM;aa++)
						{
							for (bb=0;bb<DIM;bb++)
							{
								mu_ind += beta_lab[aa][bb][bb][i][ff]*u_vec[qq][c][aa]*v_vec[qq][c][bb]*v_vec[qq][c][bb];
								for (cc=bb+1;cc<DIM;cc++)
								{
									mu_ind += 2.0*beta_lab[aa][bb][cc][i][ff]*u_vec[qq][c][aa]*v_vec[qq][c][bb]*v_vec[qq][c][cc];
								}
							}
						}

/*						mu_ind =beta_lab[0][0][0][i][ff]*u_vec[qq][c][0]*v_vec[qq][c][0]*v_vec[qq][c][0] + beta_lab[0][0][1][i][ff]*u_vec[qq][c][0]*v_vec[qq][c][0]*v_vec[qq][c][1] +
								beta_lab[0][0][2][i][ff]*u_vec[qq][c][0]*v_vec[qq][c][0]*v_vec[qq][c][2] + beta_lab[0][1][0][i][ff]*u_vec[qq][c][0]*v_vec[qq][c][1]*v_vec[qq][c][0] +
								beta_lab[0][1][1][i][ff]*u_vec[qq][c][0]*v_vec[qq][c][1]*v_vec[qq][c][1] + beta_lab[0][1][2][i][ff]*u_vec[qq][c][0]*v_vec[qq][c][1]*v_vec[qq][c][2] +
								beta_lab[0][2][0][i][ff]*u_vec[qq][c][0]*v_vec[qq][c][2]*v_vec[qq][c][0] + beta_lab[0][2][1][i][ff]*u_vec[qq][c][0]*v_vec[qq][c][2]*v_vec[qq][c][1] + 
								beta_lab[0][2][2][i][ff]*u_vec[qq][c][0]*v_vec[qq][c][2]*v_vec[qq][c][2] +
								beta_lab[1][0][0][i][ff]*u_vec[qq][c][1]*v_vec[qq][c][0]*v_vec[qq][c][0] + beta_lab[1][0][1][i][ff]*u_vec[qq][c][1]*v_vec[qq][c][0]*v_vec[qq][c][1] +
								beta_lab[1][0][2][i][ff]*u_vec[qq][c][1]*v_vec[qq][c][0]*v_vec[qq][c][2] + beta_lab[1][1][0][i][ff]*u_vec[qq][c][1]*v_vec[qq][c][1]*v_vec[qq][c][0] +
								beta_lab[1][1][1][i][ff]*u_vec[qq][c][1]*v_vec[qq][c][1]*v_vec[qq][c][1] + beta_lab[1][1][2][i][ff]*u_vec[qq][c][1]*v_vec[qq][c][1]*v_vec[qq][c][2] +
								beta_lab[1][2][0][i][ff]*u_vec[qq][c][1]*v_vec[qq][c][2]*v_vec[qq][c][0] + beta_lab[1][2][1][i][ff]*u_vec[qq][c][1]*v_vec[qq][c][2]*v_vec[qq][c][1] + 
								beta_lab[1][2][2][i][ff]*u_vec[qq][c][1]*v_vec[qq][c][2]*v_vec[qq][c][2] +
								beta_lab[2][0][0][i][ff]*u_vec[qq][c][2]*v_vec[qq][c][0]*v_vec[qq][c][0] + beta_lab[2][0][1][i][ff]*u_vec[qq][c][2]*v_vec[qq][c][0]*v_vec[qq][c][1] +
								beta_lab[2][0][2][i][ff]*u_vec[qq][c][2]*v_vec[qq][c][0]*v_vec[qq][c][2] + beta_lab[2][1][0][i][ff]*u_vec[qq][c][2]*v_vec[qq][c][1]*v_vec[qq][c][0] +
								beta_lab[2][1][1][i][ff]*u_vec[qq][c][2]*v_vec[qq][c][1]*v_vec[qq][c][1] + beta_lab[2][1][2][i][ff]*u_vec[qq][c][2]*v_vec[qq][c][1]*v_vec[qq][c][2] +
								beta_lab[2][2][0][i][ff]*u_vec[qq][c][2]*v_vec[qq][c][2]*v_vec[qq][c][0] + beta_lab[2][2][1][i][ff]*u_vec[qq][c][2]*v_vec[qq][c][2]*v_vec[qq][c][1] + 
								beta_lab[2][2][2][i][ff]*u_vec[qq][c][2]*v_vec[qq][c][2]*v_vec[qq][c][2];*/

						// For this molecule, mu_ind is the component of beta that we're interested in, after all transformations are done.
						ct += mu_ind*cst[ff][i];
						st += mu_ind*snt[ff][i];
						zt += mu_ind*mu_ind;
					}
					// Add to the average the intensity for a single frame.
					s_method[g][qq][c] += st*st + ct*ct;
					s_method_incoh[g][qq][c] += zt;
				}
				s_method_coh[g][qq][c] = s_method[g][qq][c] - s_method_incoh[g][qq][c];
			}

		}

		// For each q-value, average over all values of gamma.
		for (g=0;g<ng;g++)
		{
			for (qq=0;qq<n_used;qq++)
			{
				for (c=1;c<nbingamma;c++)
				{
					s_method[g][qq][0] += s_method[g][qq][c];
					s_method_coh[g][qq][0] += s_method_coh[g][qq][c];
					s_method_incoh[g][qq][0] += s_method_incoh[g][qq][c];
				}
			}
		}

    }

	// DMW: EDITING FROM HERE
	// Still to do:
	// 1. Output intensity.
	// 2. Check that results are correct.

	// Now, we have the intensity for many different wavevectors. What we really want is
	// the intensity as a function of the magnitude |q|. The largest possible square
	// magnitude is nmax2.
	snew(intens_total,ng);
	snew(intens_cohrt,ng);
	snew(intens_incoh,ng);
	for (g=0;g<ng;g++)
	{
		snew(intens_total[g],nmax2+1);
		snew(intens_cohrt[g],nmax2+1);
		snew(intens_incoh[g],nmax2+1);
	}
	snew(num_count,nmax2);

	// num_count[qq] will be the number of n vectors with n^2 = qq. This will
	// help us to calculate the average.
	for (g=0;g<ng;g++)
	{
		for (qq=0;qq<n_used;qq++)
		{
			n2 = narray[qq][3];
			intens_total[g][n2] += s_method[g][qq][0];
			intens_cohrt[g][n2] += s_method_coh[g][qq][0];
			intens_incoh[g][n2] += s_method_incoh[g][qq][0];
			num_count[n2] += 1;
		}
		for (qq=0;qq<nmax2+1;qq++)
		{
			if (num_count[qq]>0)
			{
				intens_total[g][qq] /= num_count[qq];
				intens_cohrt[g][qq] /= num_count[qq];
				intens_incoh[g][qq] /= num_count[qq];
			}
		}
	}

	// Now we have the intensity for each type of scattering, averaged over frames and
	// values of gamma; we have the average intensity for each possible modulus |q|.
	// (Or rather, as it stands, nx^2 + ny^2 + nz^2).
	// All that remains, I think, is to print it out!
    fpn = xvgropen("total_intensity.out", "S(q)", "q [1/nm]", "S(q)", oenv);
    fprintf(fpn, "@type xy\n");
	for (qq=0;qq<nmax2+1;qq++)
	{
		if (num_count[qq]>0)
		{
			fprintf(fpn, "%10g %10g\n",sqrt(qq)*qnorm,intens_total[0][qq]*invsize0*invgamma/(nframes));
		}
	}
	gmx_ffclose(fpn);

    fpn = xvgropen("cohrt_intensity.out", "S(q)", "q [1/nm]", "S(q)", oenv);
    fprintf(fpn, "@type xy\n");
	for (qq=0;qq<nmax2+1;qq++)
	{
		if (num_count[qq]>0)
		{
			fprintf(fpn, "%10g %10g\n",sqrt(qq)*qnorm,intens_cohrt[0][qq]*invsize0*invgamma/(nframes));
		}
	}
	gmx_ffclose(fpn);

    fpn = xvgropen("incoh_intensity.out", "S(q)", "q [1/nm]", "S(q)", oenv);
    fprintf(fpn, "@type xy\n");
	for (qq=0;qq<nmax2+1;qq++)
	{
		if (num_count[qq]>0)
		{
			fprintf(fpn, "%10g %10g\n",sqrt(qq)*qnorm,intens_incoh[0][qq]*invsize0*invgamma/(nframes));
		}
	}
	gmx_ffclose(fpn);

	exit(0);
    
    snew(x_i1, max_i);
    nframes    = 0;
    invvol_sum = 0;
    if (bPBC && (NULL != top))
    {
        gpbc = gmx_rmpbc_init(&top->idef, ePBC, natoms);
    }

    rng=gmx_rng_init(gmx_rng_make_seed());
    
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
/*           s_method[g][qq] = s_method[g][qq]/(nframes)  ;
           s_method_coh[g][qq] = s_method_coh[g][qq]/(nframes)  ;
           s_method_incoh[g][qq] = s_method_incoh[g][qq]/(nframes) ;*/
	// DMW: HAD TO COMMENT OUT.
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
    
    if (!fnBETACORR )
    {
    // print the nonlinear scattering intensity as a function of wave-vector only if you don't compute the <beta(0)*beta(r)>
    nplots = 1;
    sprintf(gtitle, "Non-linear optical scattering ");
    fpn = xvgropen(fnSFACT, "S(q)", "q (nm^-1)", "S(q)", oenv);
    sprintf(refgt, "%s", "");
    for (tt = 0; tt < nbintheta ; tt++)
    {
       theta = 2.0*theta_vec[tt]*180.0/M_PI;
       if ((round(theta) == -45.0) || (round(theta) == -30.0 ) || (round(theta) == -60.0 ) || (round(theta) == -90.0)
          || (round(theta) == -150.0)  || (round(theta) == -120.0) || (tt ==  0) || (tt == nbintheta/2))
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
 
int gmx_shscorr(int argc, char *argv[])
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
    static real              binwidth = 0.002, angle_corr = 90.0 ;
    static int               ngroups = 1, nbintheta = 10, nbingamma = 2 ,qbin = 1, nbinq = 10 ;
    static int               nkx = 0, nky = 0, nkz = 0, kern_order = 2, interp_order = 4, kmax =20;
    static int		     	 nmax = 10;
	static real				 intheta = 90;

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
		{ "-nmax",	    	FALSE, etINT, {&nmax}, "maximum modulus of n vector"},
		{ "-intheta",		FALSE, etREAL, {&intheta}, "theta value"},
        { "-cutoff",        FALSE, etREAL, {&electrostatic_cutoff}, "cutoff for the calculation of electrostatics around a molecule and/or for method=double" },
        { "-maxcutoff",        FALSE, etREAL, {&maxelcut}, "cutoff to smoothly truncate the calculation of the double sum" },
        { "-kappa",        FALSE, etREAL, {&kappa}, "screening parameter for the ewald term, i.e. erf(r*kappa)/r, in nm^-1" },
        { "-kernorder",        FALSE, etINT, {&kern_order}, "kernel order, where beta = sum_i sum_(kern_ind) c[i*kern_ind] * (feature_vec[i]-mean(feature_vec))^kern_ind " },     
        { "-splorder",        FALSE, etINT, {&interp_order}, "interpolation order for b-splines" },
        { "-kmax_spme",        FALSE, etINT, {&kmax}, "max wave vector defining the images to use in SPME" },
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

    fprintf(stderr,"Start indexing the atoms to each molecule\n");
    dipole_atom2mol(&gnx[0], grpindex[0], &(top->mols));

    do_shscorr(top, ftp2fn(efTRX, NFILE, fnm),
            opt2fn("-o", NFILE, fnm), opt2fn("-otheta", NFILE, fnm), angle_corr,
           fnVCOEFF, fnVGRD, fnVINP, fnRGRDO, fnCOEFFO,
           fnRGRDH, fnCOEFFH, fnMAP, fnBETACORR, fnFTBETACORR,
           fnREFMOL, methodt[0], kernt[0], bIONS, catname, anname, bPBC,  qbin, nbinq,
           kern_order, std_dev_dens, fspacing, binwidth,
           nbintheta, nbingamma, pin_angle, pout_angle, 
           electrostatic_cutoff, maxelcut, kappa, interp_order, kmax, kernstd, gnx, grpindex, grpname, ngroups, oenv, nmax, intheta);
    return 0;
}

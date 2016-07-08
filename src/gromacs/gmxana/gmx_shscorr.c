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
                   const output_env_t oenv, int nmax, int n2max, real intheta)

{
    FILE          *fp, *fpn;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, nanions, ncations, i, j, k, qq, n, c, tt, rr, nframes, nfaces, gr_ind, nbin, aa, bb, cc;
    real         ***s_method, ***s_method_coh, ***s_method_incoh, **temp_method, ****s_method_t, ****s_method_coh_t, ****s_method_incoh_t, ***mu_sq_t, ***coh_temp;
    real           qnorm, maxq, incoh_temp = 0.0, tot_temp = 0.0, gamma = 0.0 ,theta0 = 5.0, check_pol;
    real          *cos_t, *sin_t, ****cos_tq, ****sin_tq, mu_sq =0.0, mod_f ;
    real         **field_ad, electrostatic_cutoff2, electrostatic_cutoff, max_spacing, maxelcut2,  invkappa2,  ***beta_mol, *betamean, ****mu_ind_t;
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
    int           *chged_atom_indexes, n_chged_atoms,nrp;
    int		       nx,ny,nz,maxnpoint,**narray,**narray2,nmax2,n2max2,n_used,n2,n2_2,ii,jj,kk,l,ff,*num_count,nn,nxa,nya,nza,nxb,nyb,nzb,*repeat_list,*num_repeats,**to_repeat,nmx,ss,n2a;
    real	      **kvec,***u_vec,***v_vec,**basis,*coeff,saout,sain,caout,cain,*****beta_lab,*st,*ct,*zt;
    rvec			**all_r;
    real			***intens_total,***intens_cohrt,***intens_incoh,****s_array,****c_array,rval,dotp,mu_ind = 0.0,*induced_mu,xx,yy,zz;

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
    qnorm = M_PI*2.0/(rmax*2.0);

	if (kern[0] != 'n'){fprintf(stderr,"WARNING: This code has been written for kern[0] = n, so may not work well otherwise!\n");}

/*******************************MOLECULAR-FRAME BETA************************************************************************************************************/

	// Initialize molecular-frame beta tensor
	snew(beta_mol,DIM);
	for (i=0;i<DIM;i++){snew(beta_mol[i],DIM);}
	for (i=0;i<DIM;i++){for (j=0;j<DIM;j++){snew(beta_mol[i][j],DIM);}}

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

	beta_mol[2][2][2] = 1.0;
	// DMW: This should be changed later on!

/*******************************K-VECTORS***********************************************************************************************************************/

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
	n2max2 = n2max*n2max;
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

	snew(narray2,n_used);
	for (i=0;i<n_used;i++)
	{
		snew(narray2[i],4);
	}
	for (i=0;i<n_used;i++)
	{
		for (j=0;j<4;j++)
		{
			narray2[i][j] = narray[i][j];
		}
	}
	sfree(narray);

	// OK, let's go through and see how many q values we don't have to do because their q-vectors are multiples
	// of existing ones.
	snew(repeat_list,n_used);
	for (qq=0;qq<n_used;qq++)
	{
		repeat_list[qq] = 0;
	}
	for (qq=0;qq<n_used;qq++)
	{
		for (rr=qq+1;rr<n_used;rr++)
		{
			if (repeat_list[rr]==0)
			{
				if ( (narray2[qq][1]*narray2[rr][2]) == (narray2[qq][2]*narray2[rr][1]))
				{
					if ( (narray2[qq][0]*narray2[rr][2]) == (narray2[qq][2]*narray2[rr][0]))
					{
						if ( (narray2[qq][0]*narray2[rr][1]) == (narray2[qq][1]*narray2[rr][0]))
						{
							repeat_list[rr]=qq;
						}
					}
				}
			}
		}
	}

	// Now for each of the groups that are related by scale factors, find the shortest vector.
	snew(num_repeats,n_used);
	for (qq=0;qq<n_used;qq++){num_repeats[qq]=0;}
	for (qq=0;qq<n_used;qq++)
	{
		// We only want to look at the vectors that are "originals".
		if (repeat_list[qq]==0)
		{
			rr = qq;
			for (ss=qq+1;ss<n_used;ss++)
			{
				if (repeat_list[ss] == qq)
				{
					if (narray2[ss][3] < narray2[rr][3]){rr = ss;}
				}
			}
			// Now, rr is the vector in the same direction as qq with the smallest magnitude. This is the one we want to keep.
			num_repeats[rr] = 1;
		}
	}

	nrp = 0;
	for (rr=0;rr<n_used;rr++)
	{
		if (num_repeats[rr]==1){nrp++;}
	}
	snew(narray,nrp);
	for (i=0;i<nrp;i++)
	{
		snew(narray[i],4);
	}
	nmx = 0;
	for (rr=0;rr<n_used;rr++)
	{
		if (num_repeats[rr]==1)
		{
			for (j=0;j<4;j++)
			{
				narray[nmx][j] = narray2[rr][j];
			}
			nmx++;
		}
	}
	sfree(narray2);
	n_used = nrp;

	snew(kvec,n_used);
	for (qq=0;qq<n_used;qq++)
	{
		snew(kvec[qq],4);
		kvec[qq][0] = qnorm * narray[qq][0];
		kvec[qq][1] = qnorm * narray[qq][1];
		kvec[qq][2] = qnorm * narray[qq][2];
		kvec[qq][3] = qnorm * sqrt(1.0*narray[qq][3]);
	}

	sfree(num_repeats);
	sfree(repeat_list);

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
	
	fprintf(stderr,"\nUsing %i points in reciprocal space\n",n_used);

	// We now have a certain number of q vectors, along with the U and V vectors corresponding to
	// them, for a certain number of scattering planes each.

/**************************INITIALIZE ARRAYS DEPENDING ON NUMBER OF FRAMES**************************************************************************************/

	// Find number of frames.
	nframes = 0;
	do{nframes++;} while (read_next_x(oenv,status,&t,x,box));

	// Initialize coordinate array.
	snew(all_r,isize0);
	for (i=0;i<isize0;i++){snew(all_r[i],molsize);}

	// Initialize sine and cosine arrays.
	snew(s_array,nframes);
	snew(c_array,nframes);
	for (ff=0;ff<nframes;ff++){snew(s_array[ff],isize0);snew(c_array[ff],isize0);}
	for (ff=0;ff<nframes;ff++){for (i=0;i<isize0;i++){snew(s_array[ff][i],3);snew(c_array[ff][i],3);}}
	for (ff=0;ff<nframes;ff++){for (i=0;i<isize0;i++){for (j=0;j<3;j++){snew(s_array[ff][i][j],nmax+1);snew(c_array[ff][i][j],nmax+1);}}}

    //initialize beta tensor

	snew(beta_lab, DIM);
    for (i = 0; i < DIM; i++){snew(beta_lab[i], DIM);}
    for (i = 0; i < DIM; i++){for (j = 0; j < DIM; j++){snew(beta_lab[i][j], DIM);}}
	for (i = 0;i < DIM; i++){for (j = 0;j < DIM; j++){for (k = 0;k < DIM; k++){snew(beta_lab[i][j][k],isize0);}}}
	for (i = 0;i < DIM; i++){for (j = 0;j < DIM; j++){for (k = 0;k < DIM; k++){for (ii = 0;ii < isize0;ii++){snew(beta_lab[i][j][k][ii],nframes);}}}}

/**************************INITIALIZE ARRAYS USED TO CALCULATE INTENSITIES**************************************************************************************/

	snew(ct,nmax2+1);
	snew(st,nmax2+1);
	snew(zt,nmax2+1);

	// We will be calculating the intensity as a function of the magnitude |q|. The largest
	// possible square magnitude is n2max2.
	snew(intens_total,ng);
	snew(intens_cohrt,ng);
	snew(intens_incoh,ng);
	for (g=0;g<ng;g++)
	{
		snew(intens_total[g],nframes);
		snew(intens_cohrt[g],nframes);
		snew(intens_incoh[g],nframes);
	}
	for (g=0;g<ng;g++)
	{
		for (ff=0;ff<nframes;ff++)
		{
			snew(intens_total[g][ff],n2max2+1);
			snew(intens_cohrt[g][ff],n2max2+1);
			snew(intens_incoh[g][ff],n2max2+1);
		}
	}
	snew(num_count,n2max2+1);

/**************************READ FRAMES IN AND FILL ARRAYS*******************************************************************************************************/

	// Coefficients for sine and cosine arrays.
	for (j=0;j<3;j++)
	{
		coeff[j] = 2.0 * M_PI / box[j][j];
	}

	// Read in frames.

	for (ff=0;ff<nframes;ff++)
	{
		if (ff==0){read_first_x(oenv, &status, fnTRX, &t, &x, box);}
		else{read_next_x(oenv,status,&t,x,box);}
		for (ii=0;ii<isize0;ii++)
		{
			ind0 = mols->index[molindex[0][ii]];
			for (j=0;j<molsize;j++)
			{
				copy_rvec(x[ind0+j],all_r[ii][j]);
			}
		}

		// Fill sine and cosine arrays.
		for (i=0;i<isize0;i++)
		{
			for (j=0;j<3;j++)
			{
				rval = all_r[i][0][j];
				s_array[ff][i][j][0] = 0.0;
				c_array[ff][i][j][0] = 1.0;
				q_xi = rval * coeff[j];
				s_array[ff][i][j][1] = sin(q_xi);
				c_array[ff][i][j][1] = cos(q_xi);
				for (nn=2;nn<nmax+1;nn++)
				{
					s_array[ff][i][j][nn] = s_array[ff][i][j][nn-1]*c_array[ff][i][j][1] + c_array[ff][i][j][nn-1]*s_array[ff][i][j][1];
					c_array[ff][i][j][nn] = c_array[ff][i][j][nn-1]*c_array[ff][i][j][1] - s_array[ff][i][j][nn-1]*s_array[ff][i][j][1];
				}
			}
		}

		// For each molecule, rotate the molecular hyperpolarizability tensor into the lab frame (later on,
		// we will multiply by the elements of the polarization vectors).
		for (i = 0;i<isize0;i++)
		{
			// For this frame and this molecule, calculate the direction cosines.
			for (ii=0;ii<molsize;ii++)
			{
				pbc_dx(&pbc,all_r[i][ii],all_r[i][0],xmol[ii]);
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

	// We no longer need our r array, so let's free it up now.
	sfree(all_r);

/*************************CALCULATE INTENSITY*******************************************************************************************************************/

	fprintf(stderr,"\n");

	// For each group we take every frame, and for each frame we take each different q-vector generated,
	// and each value of gamma for each q-vector. This will give us several different scattering planes,
	// each with their own u and v vectors, i.e. n_used*nbingamma diffrent experimental setups.
	// For each setup, we calculate the intensity of ESHS light scattered.

//    for (g = 0; g < ng; g++)
	for (ff=0;ff<nframes;ff++)
    {
		// For each q vector, let's calculate the intensity.
//		for (qq=0;qq<n_used;qq++)
//		for (ff=0;ff<nframes;ff++)
		for (g = 0;g<ng;g++)
		{
//			for (ff=0;ff<nframes;ff++)
			for (qq=0;qq<n_used;qq++)
			{
//				fprintf(stderr,"Doing q-point number %i of %i\n",qq+1,n_used);
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

				n2 = narray[qq][3];
				// If n2a is greater than 1, then this means that we can use this k-vector to get the intensities for higher k-vectors.
				// (and the number of these vectors we can get is equal to n2a).
				n2a = (int)sqrt((float)nmax2/n2);

				for (c=0;c<nbingamma;c++)
				{
					// For this given group,
					// q-vector,
					// and value of gamma,
					// we are going to calculate the intensity.
					for (j=0;j<n2a;j++)
					{
						st[j] = 0.0;
						ct[j] = 0.0;
						zt[j] = 0.0;
					}
					for (i=0;i<isize0;i++)
					{
						mu_ind = 0.0;

						for (aa=0;aa<DIM;aa++)
						{
							mu_ind += beta_lab[aa][0][0][i][ff]*u_vec[qq][c][aa]*v_vec[qq][c][0]*v_vec[qq][c][0] +
									  beta_lab[aa][1][1][i][ff]*u_vec[qq][c][aa]*v_vec[qq][c][1]*v_vec[qq][c][1] +
									  beta_lab[aa][2][2][i][ff]*u_vec[qq][c][aa]*v_vec[qq][c][2]*v_vec[qq][c][2];
							for (bb=0;bb<DIM;bb++)
							{
								for (cc=bb+1;cc<DIM;cc++)
								{
									mu_ind += 2.0*beta_lab[aa][bb][cc][i][ff]*u_vec[qq][c][aa]*v_vec[qq][c][bb]*v_vec[qq][c][cc];
								}
							}
						}

						// For this molecule, mu_ind is the component of beta that we're interested in, after all transformations are done.
//						q_xi = kvec[qq][0]*all_r[ff][i][0][0] + kvec[qq][1]*all_r[ff][i][0][1] + kvec[qq][2]*all_r[ff][i][0][2];

						for (j=0;j<n2a;j++)
						{
							nx = nxa*(j+1);
							ny = nya*(j+1);
							nz = nza*(j+1);
							ct[j] += mu_ind*( c_array[ff][i][0][nx]*c_array[ff][i][1][ny]*c_array[ff][i][2][nz] - 
								nxb*nyb*s_array[ff][i][0][nx]*s_array[ff][i][1][ny]*c_array[ff][i][2][nz] - 
								nxb*nzb*s_array[ff][i][0][nx]*c_array[ff][i][1][ny]*s_array[ff][i][2][nz] - 
								nyb*nzb*c_array[ff][i][0][nx]*s_array[ff][i][1][ny]*s_array[ff][i][2][nz]);
							st[j] += mu_ind*( nxb*s_array[ff][i][0][nx]*c_array[ff][i][1][ny]*c_array[ff][i][2][nz] + 
								nyb*c_array[ff][i][0][nx]*s_array[ff][i][1][ny]*c_array[ff][i][2][nz] + 
								nzb*c_array[ff][i][0][nx]*c_array[ff][i][1][ny]*s_array[ff][i][2][nz] - 
								nxb*nyb*nzb*s_array[ff][i][0][nx]*s_array[ff][i][1][ny]*s_array[ff][i][2][nz]);
							zt[j] += mu_ind*mu_ind;
						}
						
					}
					// Add to the average the intensity for a single frame and gamma.
					for (j=0;j<n2a;j++)
					{
						n2_2 = n2*(j+1)*(j+1);
						intens_total[g][0][n2_2] += st[j]*st[j] + ct[j]*ct[j];
						intens_incoh[g][0][n2_2] += zt[j];
						intens_cohrt[g][0][n2_2] += intens_total[g][0][n2_2] - intens_incoh[g][0][n2_2];
						num_count[n2_2] += 1;
					}
				}
			}

		}

    }

/****************************************************************** PRINT OUT ********************************************************************************************/

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
			fprintf(fpn, "%10g %10g\n",sqrt(qq)*qnorm,intens_total[0][0][qq]*invsize0/num_count[qq]);
		}
	}
	gmx_ffclose(fpn);

    fpn = xvgropen("cohrt_intensity.out", "S(q)", "q [1/nm]", "S(q)", oenv);
    fprintf(fpn, "@type xy\n");
	for (qq=0;qq<nmax2+1;qq++)
	{
		if (num_count[qq]>0)
		{
			fprintf(fpn, "%10g %10g\n",sqrt(qq)*qnorm,intens_cohrt[0][0][qq]*invsize0/num_count[qq]);
		}
	}
	gmx_ffclose(fpn);

    fpn = xvgropen("incoh_intensity.out", "S(q)", "q [1/nm]", "S(q)", oenv);
    fprintf(fpn, "@type xy\n");
	for (qq=0;qq<nmax2+1;qq++)
	{
		if (num_count[qq]>0)
		{
			fprintf(fpn, "%10g %10g\n",sqrt(qq)*qnorm,intens_incoh[0][0][qq]*invsize0/num_count[qq]);
		}
	}
	gmx_ffclose(fpn);
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
    static int		     	 nmax = 10,n2max = 20;
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
		{ "-nmax",	    	FALSE, etINT, {&nmax}, "maximum modulus of n vector for sphere"},
		{ "-n2max",			FALSE, etINT, {&n2max}, "maximum modulus of n vector for cutoff"},
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
           electrostatic_cutoff, maxelcut, kappa, interp_order, kmax, kernstd, gnx, grpindex, grpname, ngroups, oenv, nmax, n2max, intheta);
    return 0;
}

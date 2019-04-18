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
                   const char *fnSFACT,
                   const char *fnBETACORR, const char *fnREFMOL,
                   const char *kern,
                   gmx_bool bIONS, char *catname, char *anname, gmx_bool bPBC, 
                   real binwidth, int nbingamma, real pin_angle, real pout_angle,
                   int *isize, int  *molindex[], char **grpname,
                   const output_env_t oenv, int nmax, int n2max, real intheta, int skip, char *betafile)

{
    FILE          *fp, *fpn;
    t_trxstatus   *status;
    char           outf1[STRLEN], outf2[STRLEN];
    char           title[STRLEN], gtitle[STRLEN], refgt[30];
    int            g, natoms, nanions, ncations, i, j, k, qq, n, c, tt, rr, nframes, nfaces, gr_ind, nbin, aa, bb, cc;
    real         ***s_method, ***s_method_coh, ***s_method_incoh, **temp_method, ****s_method_t, ****s_method_coh_t, ****s_method_incoh_t, ***mu_sq_t, ***coh_temp;
    real           qnorm, maxq, incoh_temp = 0.0, tot_temp = 0.0, gamma = 0.0 ,theta0 = 5.0, check_pol;
    real          *cos_t, *sin_t, ****cos_tq, ****sin_tq, mu_sq =0.0, mod_f ;
    real         **field_ad, max_spacing, ***beta_mol, *betamean, ****mu_ind_t;
    int            max_i, isize0, ind0, indj;
    real           t, rmax2, rmax,  r, r_dist, r2, q_xi, dq;
    real          *inv_segvol, normfac, segvol, spherevol, prev_spherevol, invsize0, invgamma, invhbinw, inv_width,  theta=0, *theta_vec;
    rvec          *x, xcm, xcm_transl, dx,  *x_i1, xi, x01, x02, qvec, **arr_qvec, **arr_qvec_faces ,vec_polin, vec_polout, ***vec_pout_theta_gamma, ***vec_pin_theta_gamma;
    rvec           pol_perp, pol_par,vec_kout, vec_2kin, pol_in1, pol_in2, vec_kout_2kin ;
    rvec           xvec, yvec, zvec, *xmol, *xref, Emean;
    real          *qref;
    matrix         cosdirmat,invcosdirmat; 
    real           invvol, invvol_sum;
//    t_Map         *Map=NULL;
    t_Kern        *Krr = NULL;
    t_Kern        *SKern_rho_O = NULL;
    t_Kern        *SKern_rho_H = NULL;
    t_Kern        *SKern_E = NULL;
    t_Ion         *Cation=NULL, *Anion=NULL;
    t_inputrec    *ir=NULL;
    t_complex   ***FT_pair_pot;
    matrix         box, box_pbc;
    rvec           grid_spacing, grid_invspacing;
    real           dens_deb, inv_tot_npoints_local_grid;
    int            *gridsize;
    int            ePBC = -1, ePBCrdf = -1;
    int            nplots = 1;
    t_block       *mols = NULL;
    t_atom        *atom = NULL;
    t_pbc          pbc;
    gmx_rmpbc_t    gpbc = NULL;
    gmx_rng_t      rng = NULL;
    int            mol, a, molsize,nfr;
    int            atom_id_0, nspecies_0, atom_id_1, nspecies_1;
    int           *chged_atom_indexes, n_chged_atoms,nrp;
    int		       nx,ny,nz,maxnpoint,**narray,**narray2,nmax2,n2max2,n_used,n2,n2_2,ii,jj,kk,l,ff,*num_count,nn,nxa,nya,nza,nxb,nyb,nzb,*repeat_list,*num_repeats,**to_repeat,nmx,ss,n2a;
    real	      ***u_vec,***v_vec,**basis,*coeff,saout,sain,caout,cain,**beta_lab,*st,*ct,*zt;
	real			***uvv;
    real			*intens_total,*intens_cohrt,*intens_incoh,***s_array,***c_array,rval,dotp,mu_ind = 0.0,*induced_mu,xx,yy,zz;
    char commentline[256];
    int n_outputs;

    fprintf(stderr,"Initialize number of atoms, get charge indexes, the number of atoms in each molecule and get the reference molecule\n");
    atom = top->atoms.atom;
    mols = &(top->mols);
    isize0 = isize[0];
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

//        fprintf(stderr, "rmax2 = %f\n", rmax2);
//        maxelcut2 = maxelcut*maxelcut; 
        if (fnBETACORR)
        {
            fprintf(stderr, "number of bins for <beta(0)*beta(r)> = %d\n", nbin);
            nfaces = 1;
            if (nbingamma >1 )
            {
                gmx_fatal(FARGS, "when computing <beta(0)*beta(r)> choose nplanes = 1 and nbintheta = 1");
            }
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

/*******************************OPEN BETA FILE******************************************************************************************************************/

	FILE *all_betas;
	if (strcmp(betafile,""))
	{
		all_betas = gmx_ffopen(betafile, "r");
	} else {
		fprintf(stderr,"Using vector beta\n");
	}

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
				beta_mol[ii][jj][kk] = 0.0;
			}
		}
	}

	beta_mol[2][2][2] = 1.0;

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
	for (nx = 0;nx <= nmax;nx+=(1 + (nx/skip)))
	{
		for (ny = 0;ny <= nmax;ny+=(1 + (ny/skip)))
		{
			for (nz = 0;nz <= nmax;nz+=(1 + (nz/skip)))
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
		snew(narray[i],5);
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
			narray[nmx][4] = (int)sqrt((float)n2max2/narray[nmx][3]);
			nmx++;
		}
	}
	sfree(narray2);
	n_used = nrp;

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

	snew(uvv,DIM*DIM*DIM);
	for (i=0;i<DIM*DIM*DIM;i++)
	{
		snew(uvv[i],n_used);
	}
	for (i=0;i<DIM*DIM*DIM;i++)
	{
		for (qq=0;qq<n_used;qq++)
		{
			snew(uvv[i][qq],nbingamma);
		}
	}
//	snew(uvv,DIM);
//	for (ii=0;ii<DIM;ii++){snew(uvv[ii],DIM);}
//	for (ii=0;ii<DIM;ii++){for (jj=0;jj<DIM;jj++){snew(uvv[ii][jj],DIM);}}
//	for (ii=0;ii<DIM;ii++){for (jj=0;jj<DIM;jj++){for (kk=0;kk<DIM;kk++){snew(uvv[ii][jj][kk],n_used);}}}
//	for (ii=0;ii<DIM;ii++){for (jj=0;jj<DIM;jj++){for (kk=0;kk<DIM;kk++){for (qq=0;qq<n_used;qq++){snew(uvv[ii][jj][kk][qq],nbingamma);}}}}

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
				basis[0][j] = (float)narray[i][j];
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

	// uvv matrix.
	rr = 0;
	for (ii=0;ii<3;ii++)
	{
		for (jj=0;jj<3;jj++)
		{
			for (kk=jj;kk<3;kk++)
			{
				for (qq=0;qq<n_used;qq++)
				{
					for (c=0;c<nbingamma;c++)
					{
						uvv[rr][qq][c] = u_vec[qq][c][ii]*v_vec[qq][c][jj]*v_vec[qq][c][kk];
						if (kk!=jj){uvv[rr][qq][c]*=2.0;}
					}
				}
				rr++;
			}
		}
	}
	sfree(u_vec);
	sfree(v_vec);
	
	fprintf(stderr,"\nUsing %i points in reciprocal space\n",n_used);

	// We now have a certain number of q vectors, along with the U and V vectors corresponding to
	// them, for a certain number of scattering planes each.

/**************************INITIALIZE ARRAYS DEPENDING ON NUMBER OF FRAMES**************************************************************************************/

	// Find number of frames.
	nframes = 0;
	do{nframes++;} while (read_next_x(oenv,status,&t,x,box));
	fprintf(stderr,"\n\nTotal number of frames: %i\n\n",nframes);

	// Initialize sine and cosine arrays.
	snew(s_array,isize0);
	snew(c_array,isize0);
	for (i=0;i<isize0;i++){snew(s_array[i],3);snew(c_array[i],3);}
	for (i=0;i<isize0;i++){for (j=0;j<3;j++){snew(s_array[i][j],n2max+1);snew(c_array[i][j],n2max+1);}}

    //initialize beta tensor

	snew(beta_lab, DIM*DIM*DIM);
	for (i=0;i<DIM*DIM*DIM;i++)
	{
		snew(beta_lab[i],isize0);
	}

/**************************INITIALIZE ARRAYS USED TO CALCULATE INTENSITIES**************************************************************************************/

	snew(ct,n2max2+1);
	snew(st,n2max2+1);
	snew(zt,n2max2+1);

	// We will be calculating the intensity as a function of the magnitude |q|. The largest
	// possible square magnitude is n2max2.
	snew(intens_total,n2max2+1);
	snew(intens_cohrt,n2max2+1);
	snew(intens_incoh,n2max2+1);
	snew(num_count,n2max2+1);

/**************************READ FRAMES IN AND FILL ARRAYS*******************************************************************************************************/
/**************************CALCULATE INTENSITY******************************************************************************************************************/

	// Coefficients for sine and cosine arrays.
	for (j=0;j<3;j++)
	{
		coeff[j] = 2.0 * M_PI / box[j][j];
	}

	// Read in frames.
	read_first_x(oenv, &status, fnTRX, &t, &x, box);
	do
	{

		// Read in the comment line of the beta file
		if (strcmp(betafile,"")) {
			n_outputs = fscanf(all_betas, "%s",commentline);
		}

		// Fill sine and cosine arrays.
		for (i=0;i<isize0;i++)
		{
			ind0 = mols->index[molindex[0][i]];
			for (j=0;j<3;j++)
			{
				rval = x[ind0][j];
				s_array[i][j][0] = 0.0;
				c_array[i][j][0] = 1.0;
				q_xi = rval * coeff[j];
				s_array[i][j][1] = sin(q_xi);
				c_array[i][j][1] = cos(q_xi);
				for (nn=2;nn<n2max+1;nn++)
				{
					s_array[i][j][nn] = s_array[i][j][nn-1]*c_array[i][j][1] + c_array[i][j][nn-1]*s_array[i][j][1];
					c_array[i][j][nn] = c_array[i][j][nn-1]*c_array[i][j][1] - s_array[i][j][nn-1]*s_array[i][j][1];
				}
			}

			// For this molecule, read in the laboratory-frame hyperpolarizability tensor from an external file
			if (strcmp(betafile,"")) {
				for (rr=0;rr<DIM*DIM*DIM;rr++)
				{
					n_outputs = fscanf(all_betas,"%f ",&beta_lab[rr][i]);
				}
//				fprintf(stderr,"BETA %f %f %f\n",beta_lab[0][i],beta_lab[14][i],beta_lab[26][i]);
			} else {

				// For each molecule, rotate the molecular hyperpolarizability tensor into the lab frame (later on,
				// we will multiply by the elements of the polarization vectors).
				for (ii=0;ii<molsize;ii++)
				{
					pbc_dx(&pbc,x[ind0+ii],x[ind0],xmol[ii]);
				}

				calc_cosdirmat( fnREFMOL, top, molsize, ind0,  xref, xmol, &cosdirmat, &invcosdirmat, &xvec, &yvec, &zvec );

				rr = 0;
				for (aa=0;aa<DIM;aa++)
				{
					for (bb=0;bb<DIM;bb++)
					{
						for (cc=bb;cc<DIM;cc++)
						{
							beta_lab[rr][i] = 0.0;
							for (ii=0;ii<DIM;ii++)
							{
								for (jj=0;jj<DIM;jj++)
								{
									for (kk=0;kk<DIM;kk++)
									{
										beta_lab[rr][i] += beta_mol[ii][jj][kk]*cosdirmat[aa][ii]*cosdirmat[bb][jj]*cosdirmat[cc][kk];
									}
								}
							}
							rr++;
						}
					}
				}
			}
		}

		// For this frame we now take each different q-vector generated, and each value of gamma
		// for each q-vector. This will give us several different scattering planes, each with
		// their own u and v vectors, i.e. n_used*nbingamma different experimental setups.
		// For each setup, we calculate the intensity of ESHS light scattered.

		for (qq=0;qq<n_used;qq++)
		{
			nxa = narray[qq][0];
			nya = narray[qq][1];
			nza = narray[qq][2];
//			nxa = abs(nx);
//			nya = abs(ny);
//			nza = abs(nz);
			// The sine matrices are odd in their arguments; we must take this into account.
//			nxb = 2*(nx>0) - 1;
//			nyb = 2*(ny>0) - 1;
//			nzb = 2*(nz>0) - 1;

			n2 = narray[qq][3];
			// If n2a is greater than 1, then this means that we can use this k-vector to get the intensities for higher k-vectors.
			// (and the number of these vectors we can get is equal to n2a).
			n2a = narray[qq][4];

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
					mu_ind = beta_lab[0][i]*uvv[0][qq][c]   + beta_lab[1][i]*uvv[1][qq][c]   +
						 beta_lab[2][i]*uvv[2][qq][c]   + beta_lab[3][i]*uvv[3][qq][c]   +
						 beta_lab[4][i]*uvv[4][qq][c]   + beta_lab[5][i]*uvv[5][qq][c]   +
						 beta_lab[6][i]*uvv[6][qq][c]   + beta_lab[7][i]*uvv[7][qq][c]   +
						 beta_lab[8][i]*uvv[8][qq][c]   + beta_lab[9][i]*uvv[9][qq][c]   +
						 beta_lab[10][i]*uvv[10][qq][c] + beta_lab[11][i]*uvv[11][qq][c] +
						 beta_lab[12][i]*uvv[12][qq][c] + beta_lab[13][i]*uvv[13][qq][c] +
						 beta_lab[14][i]*uvv[14][qq][c] + beta_lab[15][i]*uvv[15][qq][c] +
						 beta_lab[16][i]*uvv[16][qq][c] + beta_lab[17][i]*uvv[17][qq][c];

					// For this molecule, mu_ind is the component of beta that we're interested in, after all transformations are done.

					for (j=0;j<n2a;j++)
					{
						nx = nxa*(j+1);
						ny = nya*(j+1);
						nz = nza*(j+1);
						ct[j] += mu_ind*( c_array[i][0][nx]*(c_array[i][1][ny]*c_array[i][2][nz] - 
								  		     s_array[i][1][ny]*s_array[i][2][nz]) -
								  s_array[i][0][nx]*(s_array[i][1][ny]*c_array[i][2][nz] + 
								  		     c_array[i][1][ny]*s_array[i][2][nz])); 
						st[j] += mu_ind*( s_array[i][0][nx]*(c_array[i][1][ny]*c_array[i][2][nz] - 
								  		     s_array[i][1][ny]*s_array[i][2][nz]) +
								  c_array[i][0][nx]*(s_array[i][1][ny]*c_array[i][2][nz] + 
								  		     c_array[i][1][ny]*s_array[i][2][nz])); 

						zt[j] += mu_ind*mu_ind;
					}
					
				}
				// Add to the average the intensity for a single frame and gamma.
				for (j=0;j<n2a;j++)
				{
					n2_2 = n2*(j+1)*(j+1);
					intens_total[n2_2] += st[j]*st[j] + ct[j]*ct[j];
					intens_incoh[n2_2] += zt[j];
					num_count[n2_2] += 1;
				}
			}
		}
    } while (read_next_x(oenv,status,&t,x,box));
/****************************************************************** PRINT OUT ********************************************************************************************/

	// Now we have the intensity for each type of scattering, averaged over frames and
	// values of gamma; we have the average intensity for each possible modulus |q|.
	// (Or rather, as it stands, nx^2 + ny^2 + nz^2).
	// All that remains, I think, is to print it out!
    fpn = xvgropen(fnSFACT, "S(q)", "q [1/nm]", "S(q)", oenv);
    fprintf(fpn, "@type xy\n");
	for (qq=0;qq<n2max2+1;qq++)
	{
		if (num_count[qq]>0)
		{
			fprintf(fpn, "%10g %10g %10g\n",sqrt(qq)*qnorm,intens_total[qq]*invsize0/num_count[qq],intens_incoh[qq]*invsize0/num_count[qq]);
		}
	}
	gmx_ffclose(fpn);

/*    fpn = xvgropen("cohrt_intensity.out", "S(q)", "q [1/nm]", "S(q)", oenv);
    fprintf(fpn, "@type xy\n");
	for (qq=0;qq<n2max2+1;qq++)
	{
		if (num_count[qq]>0)
		{
			fprintf(fpn, "%10g %10g\n",sqrt(qq)*qnorm,(intens_total[qq]-intens_incoh[qq])*invsize0/num_count[qq]);
		}
	}
	gmx_ffclose(fpn);

    fpn = xvgropen("incoh_intensity.out", "S(q)", "q [1/nm]", "S(q)", oenv);
    fprintf(fpn, "@type xy\n");
	for (qq=0;qq<n2max2+1;qq++)
	{
		if (num_count[qq]>0)
		{
			fprintf(fpn, "%10g %10g\n",sqrt(qq)*qnorm,intens_incoh[qq]*invsize0/num_count[qq]);
		}
	}
	gmx_ffclose(fpn);*/
}
/****************************************************************** MAIN SUBROUTINE **************************************************************************************/
 
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
    static real              pout_angle = 0.0 , pin_angle = 0.0;
    //, std_dev_dens = 0.05;
    static real              binwidth = 0.002;
    static int               ngroups = 1, nbingamma = 2 ;
    static int               nkx = 0, nky = 0, nkz = 0;
    static int		     	 nmax = 10,n2max = 20,skip = -1;
    static real				 intheta = 90;

    static const char *kernt[] = {NULL, "krr", "scalar", "none", "map", NULL};
    static char *catname = NULL;
    static char *anname =  NULL;
    static char *betafile = "";

    t_pargs            pa[] = {
        { "-nplanes",       FALSE, etINT, {&nbingamma},
        "number of scattering planes that lie on the scattered wave-vector to average over, -PI/2< gamma< PI/2" },
        { "-binw",          FALSE, etREAL, {&binwidth}, "width of bin to compute <beta_lab(0) beta_lab(r)> " },
        { "-pout",          FALSE, etREAL, {&pout_angle}, "polarization angle of outcoming beam in degrees. For P choose 0, for S choose 90" },
        { "-pin",           FALSE, etREAL, {&pin_angle}, "polarization angle of incoming beam in degrees. For P choose 0, for S choose 90" },
	{ "-nmax",	    	FALSE, etINT, {&nmax}, "maximum modulus of n vector for sphere"},
	{ "-n2max",			FALSE, etINT, {&n2max}, "maximum modulus of n vector for cutoff"},
	{ "-intheta",		FALSE, etREAL, {&intheta}, "theta value"},
	{ "-skip",			FALSE, etINT, {&skip}, "q-vector skip"},
        { "-kern",   FALSE, etENUM, {kernt}, "what method to use to compute beta"},
        { "-ions",   FALSE, etBOOL, {&bIONS}, "compute molecular hyperpolarizability when ions are present"},
        { "-cn",     FALSE, etSTR, {&catname}, "name of cation"},
        { "-an",     FALSE, etSTR, {&anname}, "name of anion"},
        { "-pbc",      FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances. Always use, results without PBC not tested." },
        { "-ng",       FALSE, etINT, {&ngroups}, 
          "Number of secondary groups, not available for now. Only tip4p water implemented." },
	{ "-beta", FALSE, etSTR, {&betafile}, "file containing beta tensors"},
    };
#define NPA asize(pa)
    const char        *fnTPS, *fnNDX , *fnBETACORR = NULL, *fnREFMOL = NULL;
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
    fnBETACORR = opt2fn_null("-betacorr", NFILE,fnm);
    fnREFMOL = opt2fn_null("-refmol", NFILE, fnm);


    if (!fnTPS && !fnNDX)
    {
        gmx_fatal(FARGS, "Neither index file nor topology file specified\n"
                  "Nothing to do!");
    }

//    if ((*kernt)[0] == 's')
//    {
//       if (!fnVCOEFF || !fnVGRD || !fnRGRDO || !fnRGRDH || !fnCOEFFO || !fnCOEFFH)
//       {
//          gmx_fatal(FARGS, "specify all files for scalar kernel using -vcoeff, -vgrid, -vinp, -rhogridO, -rhogridH, -rhocoeffO, rhocoeffH\n");
//       }
//    }
//    else if ((*kernt)[0] == 'm')
//    {
//       if (!fnMAP )
//       {
//          gmx_fatal(FARGS, "specify map file with -emap\n");
//       }
//    }
//    else if ((*kernt)[0] == 'k' )
//    {
//       if (!fnVCOEFF || !fnVGRD || !fnVINP)
//       {
//          gmx_fatal(FARGS, "specify all files for krrr using -vcoeff, -vgrid, -vinp\n");
//       }
//    }
 
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

	if (n2max < nmax){n2max=nmax;}
	if (skip < 0){skip = 2*nmax;}

    do_shscorr(top, ftp2fn(efTRX, NFILE, fnm),
            opt2fn("-o", NFILE, fnm), 
           fnBETACORR,
           fnREFMOL, kernt[0], bIONS, catname, anname, bPBC, 
           binwidth,
           nbingamma, pin_angle, pout_angle, 
           gnx, grpindex, grpname, oenv, nmax, n2max, intheta, skip, betafile);
    return 0;
}

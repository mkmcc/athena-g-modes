#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:  */

/* Custom boundary conditions for the x1 direction*/

static void bc_ix1(GridS *pGrid);

void problem(DomainS *pDomain)
{
    GridS *pGrid=(pDomain->Grid);
    int i, is = pGrid->is, ie = pGrid->ie;
    int j, js = pGrid->js, je = pGrid->je;
    int k, ks = pGrid->ks, ke = pGrid->ke;
    Real x1,x2,x3;
    Real density, pressure, vel, lx, ly;
    Real boost = par_getd("problem", "boost");

    
    
    /* size of the domain (in physical coordinates) */
    lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
    ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
    
    /* if prat=0.8, vx = -1.8965 (t&e find other vals)*/
    
    
    /* setup reigon of high pressure and low pressure */
    
    
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
                cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
                
            
                if (x1 > 0.0) {
                    pressure = 0.01;
                } else {
                    pressure = 1.0;
                    
                }
                
                
                vel = boost;
                    
                }
                
                if (x1 > 0.0) {
                    density= 1.0;
                } else {
                    density= 1.0;
                    
                }
                
                
                pGrid->U[0][j][i].d = density;
                pGrid->U[0][j][i].M1 = density*vel;
                pGrid->U[0][j][i].M2 = 0.0;
#ifndef ISOTHERMAL
                pGrid->U[0][j][i].E = pressure/Gamma_1 + (SQR(pGrid->U[0][j][i].M1) + SQR(pGrid->U[0][j][i].M2))/(2.0*(pGrid->U[0][j][i].d));
                
#endif
            }
        }
    }
    
    /* enroll special functions */
    bvals_mhd_fun(pDomain,left_x1, bc_ix1);
    
    return;
    
}
void problem_write_restart(MeshS *pM, FILE *fp)
{
    return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
    return;
}

ConsFun_t get_usr_expr(const char *expr)
{
    
    return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
    return NULL;
}



void Userwork_in_loop(MeshS *pM)
{
    return ;
}


void Userwork_after_loop(MeshS *pM)
{
    return ;
}

/*---------------------------------------------------------------------------------- Constant pressure boundary condition for the inner x1 boundary*/



static void bc_ix1(GridS *pGrid)
{
    int is = pGrid->is;
    int js = pGrid->js, je = pGrid->je;
    int ks = pGrid->ks, ke = pGrid->ke;
    int i,j,k;
#ifdef MHD
    int ju, ku; /* j-upper, k-upper */
#endif
    Real P=1.0;                /* set a constant pressure */
    
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=1; i<=nghost; i++) {
                pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
                
                pGrid->U[k][j][is-i].E = P/Gamma_1 + (SQR(pGrid->U[0][j][i].M1) + SQR(pGrid->U[0][j][i].M2))/(2*(pGrid->U[0][j][i].d));
            }
        }
    }
    
    
    return;
}





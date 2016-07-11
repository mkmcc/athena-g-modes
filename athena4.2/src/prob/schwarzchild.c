/* simulates an atmosphere stable to schwarzchild convection.  ie, g-modes. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

static void reflect_ix2(GridS *pGrid);
static void reflect_ox2(GridS *pGrid);
static Real grav_pot2(const Real x1, const Real x2, const Real x3);


void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j;
  int is,ie,js,je;
  Real amp,x1,x2,x3,lx,ly;
  Real rho, P, drho, dP;
  Real alpha, beta, n, omega, kx, dv, f;

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;

  lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

  n = sqrt(9.0/10.0);           /* brunt        */
  kx = 4.0 * PI / lx;           /* wave vector  */
  omega = 0.5 * n;              /* wave freq    */

  alpha = 1.0 / (2.0 * 5.0/3.0 * kx);
  beta  = sqrt(SQR(n/omega)-1.0);

  amp = 1.0e-3;

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,0,&x1,&x2,&x3);

      /* background */
      rho = pow(1.0 + x2/2, -3.0);
      P   = pow(1.0 + x2/2, -2.0);

      /* eigenmodes */
      dv =  amp * sin(beta*kx*x2) / beta;
      f  = -amp * (cos(beta*kx*x2) - alpha/beta * sin(beta*kx*x2));

      dv *= exp(alpha * kx * x2);
      f  *= exp(alpha * kx * x2);

      /* perturbations */
      drho = -SQR(n)/omega * dv * sin(kx*x1);
      dP   =   omega/kx    * f  * sin(kx*x1);

      /* write to grid */
      pGrid->U[0][j][i].d = rho * (1.0 + drho);
      pGrid->U[0][j][i].E =   P * (1.0 + dP)/Gamma_1;

      pGrid->U[0][j][i].M1 =  f * sin(kx*x1) * rho;
      pGrid->U[0][j][i].M2 = dv * cos(kx*x1) * rho;
    }
  }

  StaticGravPot = grav_pot2;
  bvals_mhd_fun(pDomain, left_x2,  reflect_ix2);
  bvals_mhd_fun(pDomain, right_x2, reflect_ox2);

  return;
}

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd;

  StaticGravPot = grav_pot2;
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x2,  reflect_ix2);
      bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x2, reflect_ox2);
    }
  }

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
}

void Userwork_after_loop(MeshS *pM)
{
}



/*------------------------------------------------------------------------------
 * reflect_ix2: special reflecting boundary functions in x2 for 2D sims
 */
static void reflect_ix2(GridS *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][js-j][i]    =  pGrid->U[k][js+(j-1)][i];
        pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js-j][i].M2; /* reflect 2-mom. */
        pGrid->U[k][js-j][i].E +=
          pGrid->U[k][js+(j-1)][i].d*1.0*(2*j-1)*pGrid->dx2/Gamma_1;
      }
    }
  }

  return;
}


/*------------------------------------------------------------------------------
 * reflect_ox2: special reflecting boundary functions in x2 for 2D sims
 */
static void reflect_ox2(GridS *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke, ku;
  int i,j,k,il,iu,jl,ju; /* i/j-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][je+j][i]    =  pGrid->U[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je+j][i].M2; /* reflect 2-mom. */
        pGrid->U[k][je+j][i].E -=
          pGrid->U[k][je-(j-1)][i].d*1.0*(2*j-1)*pGrid->dx2/Gamma_1;
      }
    }
  }

  return;
}

static Real grav_pot2(const Real x1, const Real x2, const Real x3)
{
  return x2;
}

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"



/* problem() function defines the initial condition */
void problem(DomainS *pDomain)
{
  /* define variables */
  GridS *pGrid = pDomain->Grid;
  int i,j,is,ie,js,je;
  Real x1,x2,x3,lx,ly;
  Real rho,P,vx,vy;
  Real vb,face;
 
  face = par_getd("problem", "interface");
  vb = par_getd("problem", "boost"); 
  
  /* limits for loops (integer coordinates) */
  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;

  /* size of the domain (in physical coordinates) */
  lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];


  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,0,&x1,&x2,&x3);

      vy = 0.0;
      vx = vb;
      rho = 1.0;

      if (x1 > face){
         P = 0.01;
      }
      else {
         P = 1;
      }
     

      pGrid->U[0][j][i].d = rho;
      pGrid->U[0][j][i].E = P/Gamma_1 + 0.5*rho*(SQR(vx) + SQR(vy));

      pGrid->U[0][j][i].M1 = vx * rho;
      pGrid->U[0][j][i].M2 = vy * rho;

    }
  }


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
   int nl, nd, ntot;
   GridS *pGrid;
   int is, ie, js, je, ks, ke;
   int i,j,k;
   Real KE,TE;


   for (nl=0; nl<=(pM->NLevels)-1; nl++) {
      for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++) {
         if (pM->Domain[nl][nd].Grid != NULL) {
            pGrid = pM->Domain[nl][nd].Grid;

            is = pGrid->is; ie = pGrid->ie;
            js = pGrid->js; je = pGrid->je;
            ks = pGrid->ks; ke = pGrid->ke;

            for (j=js; j<=je; j++){
               for (i=is; i<=ie; i++){
                  KE = (SQR(pGrid->U[k][j][i].M1)+SQR(pGrid->U[k][j][i].M2))/(2*pGrid->U[k][j][i].d);
                  TE = pGrid->U[k][j][i].E - KE;
                  TE -= SQR(pGrid->U[k][j][i].d)*(1/SQR(TE))*pGrid->dt;
                  if (TE < 3*pGrid->U[k][j][i].d*0.01/2){
                     TE = 3*pGrid->U[k][j][i].d*0.01/2;
                  }
                  pGrid->U[k][j][i].E = KE + TE;
               }
            }

         }
      }
   }

   return;
}

void Userwork_after_loop(MeshS *pM)
{
}




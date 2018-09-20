#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <tools.hpp>

using namespace bemtool;
using namespace std;

int main(int argc, char* argv[]){

  // Loading of the mesh
  Geometry node("torus.msh");
  cout << "nb node: " << NbNode(node) << endl;
  
  Mesh2D mesh;
  mesh.Load(node,1); Orienting(mesh);
  int nb_elt = NbElt(mesh);
  cout << "nb_elt:\t" << nb_elt << endl;
  
  Real meshwidth = 0.;
  for(int j=0; j<nb_elt; j++){
    if(sqrt(Vol(mesh[j]))>meshwidth){
      meshwidth=sqrt(Vol(mesh[j]));
    }
  }
  cout << "meshwidth = " << meshwidth << endl;

  
  // Degrees of freedom
  Dof<RT0_2D> dof0(mesh);  
  int nb_dof0 = NbDof(dof0);
  cout << "nb_dof0 = " << nb_dof0 << endl;
  Dof<P0_2D> dof1(mesh);  
  int nb_dof1 = NbDof(dof1);  
  cout << "nb_dof1 = " << nb_dof1 << endl;  
  int nb_dof_tot = nb_dof0+nb_dof1;
  int offset[2] = {0,nb_dof0};
  
  cout << "nb_dof_tot = " << nb_dof_tot << endl;
  
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //    Assembly of matrices      //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  BIOp<LA_SL_3D_RT0xRT0>       S(mesh,mesh,1.);
  BIOp<LA_SL_3D_DivRT0xDivRT0> R(mesh,mesh,1.);
  BIOp<LA_SL_3D_P0xDivRT0>     T(mesh,mesh,1.);
  BIOp<LA_SL_3D_P0xP0>         S2(mesh,mesh,1.);  
  
  DenseMatrix<Cplx>  A(nb_dof_tot,nb_dof_tot);
  DenseMatrix<Cplx>  M(nb_dof_tot,nb_dof_tot);
  DenseMatrix<Cplx>  P(nb_dof_tot,nb_dof_tot);
  
  progress bar("Assemblage",nb_elt);
  for(int j=0; j<nb_elt; j++){bar++;
    for(int k=0; k<nb_elt; k++){      
      
      {//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	const N3& jj = offset[0]+dof0[j];
	const N3& kk = offset[0]+dof0[k];
	A(jj,kk) += S(j,k);
	P(jj,kk) += S(j,k);
	P(jj,kk) += R(j,k);	
      }//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
      
      {//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	const N1& jj = offset[1]+dof1[j];
	const N3& kk = offset[0]+dof0[k];
	C1x3 Tloc = T(j,k);
	A(jj,kk) += Tloc;
	A(kk,jj) += transpose(Tloc);
	M(kk,jj) += transpose(Tloc);
      }//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
      
      {//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
	const N1& jj = offset[1]+dof1[j];
	const N1& kk = offset[1]+dof1[k];
	P(jj,kk) += S2(j,k);	
      }//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//      
    }
  }
  bar.end();  
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //    Export des matrices      //
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  ofstream file;
  
  file.open("A.txt");
  for(int j=0; j<nb_dof_tot; j++){
    for(int k=0; k<nb_dof_tot; k++){
      file << A(j,k).real() << " ";
    }
  }
  file.close();
  
  file.open("P.txt");
  for(int j=0; j<nb_dof_tot; j++){
    for(int k=0; k<nb_dof_tot; k++){
      file << P(j,k).real() << " ";
    }
  }
  file.close();
  
}

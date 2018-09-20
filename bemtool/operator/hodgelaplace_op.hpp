//===================================================================
//
//  Copyright 2017 Xavier Claeys
//
//  bemtool is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  bemtool is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with bemtool.  If not, see <http://www.gnu.org/licenses/>
//
//====================================================================
#ifndef BEMTOOL_OPERATOR_HODGELAPLACE_HPP
#define BEMTOOL_OPERATOR_HODGELAPLACE_HPP

#include "operator.hpp"

namespace bemtool {
  
  
  /*=====================================
    Simple couche Laplace P1 x P1
    =====================================*/

  template <typename PhiX, typename PhiY>
  class BIOpKernel<HL,SL_OP,3,PhiX,PhiY>{
    
  public:
    typedef BIOpKernelTraits<HL,SL_OP,3,PhiX,PhiY> Trait;
    
  private:
  const typename Trait::MeshX&      meshx;
  const typename Trait::MeshY&      meshy;
        typename Trait::MatType     inter;
        typename Trait::JacX        dx;
        typename Trait::JacY        dy;
                        PhiX        phix;
                        PhiY        phiy;
  const                 std::vector<R3>&  normalx;
  const                 std::vector<R3>&  normaly;
  const                 Real        kappa;
                        R3          x0_y0,x_y,nx,ny;
                        Real        h,r;
                        Cplx        ker,val,val2;

  public:
    BIOpKernel<HL,SL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				     const typename Trait::MeshY& my,
				     const Real& k):
    meshx(mx), phix(mx), meshy(my), phiy(my), 
      normalx(NormalTo(mx)), normaly(NormalTo(my)),
      kappa(k) {};
    
    
    inline void Assign(const int& ix, const int& iy){
      const typename Trait::EltX& ex=meshx[ix];
      const typename Trait::EltY& ey=meshy[iy];
      phix.Assign(ix);
      phiy.Assign(iy);
      h     = DetJac(ex)*DetJac(ey);
      x0_y0 = ex[0]-ey[0];
      dx    = MatJac(ex);
      dy    = MatJac(ey);
      nx    = normalx[ix];
      ny    = normaly[iy];
    }
    

    inline const typename Trait::MatType&
    operator()(const typename Trait::Rdx& tx,
	       const typename Trait::Rdy& ty){
      x_y = x0_y0 + dx*tx-dy*ty;
      r   = norm2(x_y);
      ker = h/(4*pi*r);
      for(int j=0; j<Trait::nb_dof_x; j++){
	for(int k=0; k<Trait::nb_dof_y; k++){
	  val = ker*phix(j,tx)*phiy(k,ty)*(nx,ny);
	  inter(j,k) = val;
	}
      }
      return inter;
    }
    
    
    inline const Cplx&
    operator()(const typename Trait::Rdx& tx,
	       const typename Trait::Rdy& ty,
	       const int& kx, const int& ky){
      x_y = x0_y0 + dx*tx-dy*ty;
      r   = norm2(x_y);
      ker = h/(4*pi*r);
      return val2 = ker*phix(kx,tx)*phiy(ky,ty)*(nx,ny);
    }
    
    
  };
    
  typedef BIOpKernel<HL,SL_OP,3,P1_2D,P1_2D> HL_SL_3D_P1xP1;
  

  /*=====================================
    Simple couche Laplace RT0 x Curl(P1)
    =====================================*/
  
  
  template <typename PhiX, typename PhiY>
  class BIOpKernel<HL,DL_OP,3,PhiX,PhiY>{
    
  public:
    typedef BIOpKernelTraits<HL,DL_OP,3,PhiX,PhiY> Trait;
    
  private:
    const typename Trait::MeshX&     meshx;
    const typename Trait::MeshY&     meshy;
          typename Trait::MatType    inter;
          typename Trait::JacX       dx;
          typename Trait::JacY       dy;
          typename Trait::GradPhiY   grad_phiy;    
                          PhiX       phix;
                          PhiY       phiy;
    const        std::vector<R3>&    normaly;    
    const                 Real       kappa;
                          R3         x0_y0,x_y,nx,ny;
                          Real       h,r;
                          Cplx       ker,val,val2;

  public:
    BIOpKernel<HL,DL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				     const typename Trait::MeshY& my,
				     const Real& k):
    meshx(mx), phix(mx), meshy(my), phiy(my), grad_phiy(my),
      normaly(NormalTo(my)), kappa(k) {};
    
    inline void Assign(const int& ix, const int& iy){
      const typename Trait::EltX& ex=meshx[ix];
      const typename Trait::EltY& ey=meshy[iy];
      phix.Assign(ix);
      phiy.Assign(iy);
      grad_phiy.Assign(iy);
      h     = DetJac(ex)*DetJac(ey);
      x0_y0 = ex[0]-ey[0];
      dx    = MatJac(ex);
      dy    = MatJac(ey);
      ny    = normaly[iy];
    }

    inline const typename Trait::MatType&
    operator()(const typename Trait::Rdx& tx,
	       const typename Trait::Rdy& ty){
      x_y = x0_y0 + dx*tx-dy*ty;
      r   = norm2(x_y);
      ker = h/(4*pi*r);
      for(int j=0; j<Trait::nb_dof_x; j++){
	for(int k=0; k<Trait::nb_dof_y; k++){
	  inter(j,k) = ( phix(j,tx),vprod(ny,grad_phiy(k,ty)) )*ker;
	}
      }
      return inter;
    }
    
    inline const Cplx&
    operator()(const typename Trait::Rdx& tx,
	       const typename Trait::Rdy& ty,
	       const int& kx, const int& ky){
      x_y = x0_y0 + dx*tx-dy*ty;
      r   = norm2(x_y);
      ker = h/(4*pi*r);
      return val2 = ( phix(kx,tx),vprod(ny,grad_phiy(ky,ty)) )*ker;
    }
        
  };
  
  typedef BIOpKernel<HL,DL_OP,3,RT0_2D,P1_2D> HL_SL_3D_RT0xCurlP1;

  /*=====================================
    Simple couche Laplace DivRT0 x DivRT0
    =====================================*/
  
  
  template <typename PhiX, typename PhiY>
  class BIOpKernel<HL,HS_OP,3,PhiX,PhiY>{
    
  public:
    typedef BIOpKernelTraits<HL,HS_OP,3,PhiX,PhiY> Trait;
    
  private:
    const typename Trait::MeshX&     meshx;
    const typename Trait::MeshY&     meshy;
          typename Trait::MatType    inter;
          typename Trait::JacX       dx;
          typename Trait::JacY       dy;
          typename Trait::DivPhiX    div_phix;        
          typename Trait::DivPhiY    div_phiy;    
                          PhiX       phix;
                          PhiY       phiy;
    const                 Real       kappa;
                          R3         x0_y0,x_y,nx,ny;
                          Real       h,r;
                          Cplx       ker,val,val2;

  public:
    BIOpKernel<HL,HS_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				     const typename Trait::MeshY& my,
				     const Real& k):
    meshx(mx), phix(mx), div_phix(mx),
      meshy(my), phiy(my), div_phiy(my), kappa(k) {};
    
    
    inline void Assign(const int& ix, const int& iy){
      const typename Trait::EltX& ex=meshx[ix];
      const typename Trait::EltY& ey=meshy[iy];
      phix.Assign(ix);
      phiy.Assign(iy);
      div_phix.Assign(ix);
      div_phiy.Assign(iy);
      h     = DetJac(ex)*DetJac(ey);
      x0_y0 = ex[0]-ey[0];
      dx    = MatJac(ex);
      dy    = MatJac(ey);
    }
    

    inline const typename Trait::MatType&
    operator()(const typename Trait::Rdx& tx,
	       const typename Trait::Rdy& ty){
      x_y = x0_y0 + dx*tx-dy*ty;
      r   = norm2(x_y);
      ker = h/(4*pi*r);
      for(int j=0; j<Trait::nb_dof_x; j++){
	for(int k=0; k<Trait::nb_dof_y; k++){
	  inter(j,k) = div_phix(j,tx)*div_phix(k,ty)*ker;
	}
      }
      return inter;
    }
    
    
    inline const Cplx&
    operator()(const typename Trait::Rdx& tx,
	       const typename Trait::Rdy& ty,
	       const int& kx, const int& ky){
      x_y = x0_y0 + dx*tx-dy*ty;
      r   = norm2(x_y);
      ker = h/(4*pi*r);
      return val2 = div_phix(kx,tx)*div_phix(ky,ty)*ker;
    }
        
  };
  
  typedef BIOpKernel<HL,HS_OP,3,RT0_2D,RT0_2D> HL_SL_3D_DivRT0xDivRT0;


    
  /*=====================================
    Simple couche Laplace RT0 x RT0
    =====================================*/

  template <typename PhiX, typename PhiY>
  class BIOpKernel<HL,TDL_OP,3,PhiX,PhiY>{
    
  public:
    typedef BIOpKernelTraits<HL,TDL_OP,3,PhiX,PhiY> Trait;
    
  private:
  const typename Trait::MeshX&      meshx;
  const typename Trait::MeshY&      meshy;
        typename Trait::MatType     inter;
        typename Trait::JacX        dx;
        typename Trait::JacY        dy;
        typename Trait::DivPhiY     div_phiy;
                        PhiX        phix;
                        PhiY        phiy;
  const                 Real        kappa, inv_kappa2;
                        R3          x0_y0,x_y,nx,ny;
                        Real        h,r;
                        Cplx        ker,val,val2;

  public:
    BIOpKernel<HL,TDL_OP,3,PhiX,PhiY>(const typename Trait::MeshX& mx,
				     const typename Trait::MeshY& my,
				     const Real& k):
    meshx(mx), phix(mx), meshy(my), phiy(my), div_phiy(my),
      kappa(k), inv_kappa2(1.) {};
    
    
    inline void Assign(const int& ix, const int& iy){
      const typename Trait::EltX& ex=meshx[ix];
      const typename Trait::EltY& ey=meshy[iy];
      phix.Assign(ix);
      phiy.Assign(iy);
      div_phiy.Assign(iy);
      h     = DetJac(ex)*DetJac(ey);
      x0_y0 = ex[0]-ey[0];
      dx    = MatJac(ex);
      dy    = MatJac(ey);
    }
    

    inline const typename Trait::MatType&
    operator()(const typename Trait::Rdx& tx,
	       const typename Trait::Rdy& ty){
      x_y = x0_y0 + dx*tx-dy*ty;
      r   = norm2(x_y);
      ker = h/(4*pi*r);
      for(int j=0; j<Trait::nb_dof_x; j++){
	for(int k=0; k<Trait::nb_dof_y; k++){
	  val =  ker*( phix(j,tx),phiy(k,ty) );
	  inter(j,k) = val;
	}
      }
      return inter;
    }
    
    
    inline const Cplx&
    operator()(const typename Trait::Rdx& tx,
	       const typename Trait::Rdy& ty,
	       const int& kx, const int& ky){
      x_y = x0_y0 + dx*tx-dy*ty;
      r   = norm2(x_y);
      ker = h/(4*pi*r);
      return val2 = ( phix(kx,tx),phiy(ky,ty) )*ker;
    }
    
    
  };
    
  typedef BIOpKernel<HL,TDL_OP,3,RT0_2D,RT0_2D> LA_SL_3D_RT0xRT0;


  
}

















#endif

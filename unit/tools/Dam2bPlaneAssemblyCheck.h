#pragma once
#include "PlaneStressCoupling.h"
#include "RankFourTensor.h"
#include "MooseTypes.h"
#include "ADUtils.h"
#include <map>

// Checks only constitutive cross-blocks on original one-QP QUAD4 geometry.
// AD differentiates the frozen algorithmic stress linearization, not the nonlinear CDP map.
class Dam2bPlaneAssemblyCheck
{
public:
  using Tensor = AbaqusCDPStateIntegrator::SymmetricTensor;
  struct Geometry { double area; std::array<RealGradient,4> grad; };
  std::map<unsigned int,Geometry> geometry;
  std::ofstream output;

  Dam2bPlaneAssemblyCheck(const char * input_path, const char * output_path) : output(output_path)
  {
    std::ifstream input(input_path);
    if (!input || !output) throw std::runtime_error("Cannot open geometry/assembly output");
    std::string line;
    std::getline(input,line);
    while (std::getline(input,line))
    {
      if (!line.empty() && line.back()=='\r') line.pop_back();
      std::istringstream stream(line); std::string field; std::vector<double> x;
      while (std::getline(stream,field,','))
      {
        std::size_t end; double v=std::stod(field,&end);
        if (end!=field.size() || !std::isfinite(v)) throw std::runtime_error("Invalid geometry field");
        x.push_back(v);
      }
      if (x.size()!=11 || !(x[2]>0) || x[0]!=std::floor(x[0])) throw std::runtime_error("Invalid geometry row");
      Geometry g; g.area=x[2];
      for (unsigned int j=0;j<4;++j) g.grad[j]=RealGradient(x[3+2*j],x[4+2*j],0.);
      if (!geometry.emplace(static_cast<unsigned int>(x[0]),g).second) throw std::runtime_error("Duplicate geometry");
    }
    output << std::setprecision(17) << "abaqus_element,mode,row_dof,column_dof,h,old_J,new_J,frozen_tensor_actual_AD,FD_center,FD_plus,FD_minus,base_S33,plus_substeps,minus_substeps\n";
  }

  void observe(unsigned int label, unsigned int mode, const AbaqusCDPSubstepIntegrator & steps,
               const Tensor & eold, const Tensor & enew, double dt,
               const AbaqusCDPStateIntegrator::State & old,
               const AbaqusCDPSubstepIntegrator::LinearizedResult & result)
  {
    const auto & g=geometry.at(label);
    const auto & stress0=result.result.final_result.cauchy_stress;
    const std::array<std::array<unsigned int,2>,6> pairs={{{0,0},{1,1},{2,2},{0,1},{1,2},{0,2}}};
    RankFourTensor C;
    for (unsigned int row=0;row<6;++row)
      for (unsigned int col=0;col<6;++col)
      {
        auto i=pairs[row][0],j=pairs[row][1],k=pairs[col][0],l=pairs[col][1];
        double v=result.algorithmic_tangent[col][row]*(k==l ? 1. : .5);
        C(i,j,k,l)=C(j,i,k,l)=C(i,j,l,k)=C(j,i,l,k)=v;
      }
    ADReal e[3][3]={};
    for (unsigned int j=0;j<4;++j)
    {
      ADReal ux=0.,uy=0.,z=0.;
      Moose::derivInsert(ux.derivatives(),2*j,1.);
      Moose::derivInsert(uy.derivatives(),2*j+1,1.);
      Moose::derivInsert(z.derivatives(),8+j,1.);
      e[0][0]+=g.grad[j](0)*ux; e[1][1]+=g.grad[j](1)*uy;
      e[0][1]+=.5*(g.grad[j](1)*ux+g.grad[j](0)*uy);
      e[2][2]+=.25*z;
    }
    e[1][0]=e[0][1];
    ADReal stress[3][3]={};
    for (unsigned int i=0;i<3;++i)
      for (unsigned int j=0;j<3;++j)
        for (unsigned int k=0;k<3;++k)
          for (unsigned int l=0;l<3;++l) stress[i][j]+=C(i,j,k,l)*e[k][l];
    std::array<ADReal,12> residual={};
    for (unsigned int i=0;i<4;++i)
    {
      for (unsigned int a=0;a<2;++a)
        residual[2*i+a]=g.area*.95*(stress[a][0]*g.grad[i](0)+stress[a][1]*g.grad[i](1));
      residual[8+i]=g.area*.25*stress[2][2];
    }
    auto force=[&](const Tensor & s,unsigned int row) {
      if (row>=8) return g.area*.25*s[2];
      const unsigned int a=row%2,i=row/2;
      return g.area*.95*(s[a]*g.grad[i](a)+s[3]*g.grad[i](1-a));
    };
    for (unsigned int col=0;col<12;++col)
    {
      Tensor direction={};
      if (col>=8) direction[2]=.25;
      else { unsigned int a=col%2,j=col/2; direction[a]=g.grad[j](a); direction[3]=.5*g.grad[j](1-a); }
      for (double h : {1e-8,1e-9,1e-10})
      {
        auto ep=enew,em=enew;
        for (unsigned int k=0;k<6;++k) { ep[k]+=h*direction[k]; em[k]-=h*direction[k]; }
        const auto plus=steps.integrate(eold,ep,dt,old),minus=steps.integrate(eold,em,dt,old);
        for (unsigned int row=0;row<12;++row)
        {
          if ((col>=8)==(row>=8)) continue;
          double original,correction;
          if (row<8)
          {
            unsigned int a=row%2,i=row/2;
            original=g.area*.95*C(a,a,2,2)*g.grad[i](a)*.25;
            correction=g.area*PlaneStressCoupling::displacementShearColumn(C,a,g.grad[i],.25,-.05);
          }
          else
          {
            unsigned int a=col%2,j=col/2;
            original=g.area*.25*C(2,2,a,a)*g.grad[j](a);
            correction=g.area*PlaneStressCoupling::thicknessShearRow(C,a,.25,g.grad[j]);
          }
          const double fp=force(plus.final_result.cauchy_stress,row),fm=force(minus.final_result.cauchy_stress,row),f0=force(stress0,row);
          output << label << ',' << mode << ',' << row << ',' << col << ',' << h << ','
                 << original << ',' << original+correction << ',' << residual[row].derivatives()[col] << ','
                 << (fp-fm)/(2*h) << ',' << (fp-f0)/h << ',' << (f0-fm)/h << ',' << stress0[2] << ','
                 << plus.accepted_substeps << ',' << minus.accepted_substeps << '\n';
        }
      }
    }
  }
};

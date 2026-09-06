#pragma once
#include <array>
// Full Bbar(test) : C : Bbar(trial); tensor shear, no major-symmetry assumption.
namespace CDPBbar
{
template <typename Tensor>
double contract(const Tensor & C, unsigned int a, unsigned int b,
                const std::array<double,3> & g, const std::array<double,3> & h,
                double avg_g, double avg_h)
{
  double out=0;
  for (unsigned int i=0;i<3;++i)
    for (unsigned int j=0;j<3;++j)
    {
      const double test=0.5*((i==a ? g[j]:0)+(j==a ? g[i]:0))
                       +(i==j ? (avg_g-g[a])/3.0:0);
      for (unsigned int k=0;k<3;++k)
        for (unsigned int l=0;l<3;++l)
        {
          const double trial=0.5*((k==b ? h[l]:0)+(l==b ? h[k]:0))
                            +(k==l ? (avg_h-h[b])/3.0:0);
          out+=test*C(i,j,k,l)*trial;
        }
    }
  return out;
}
}

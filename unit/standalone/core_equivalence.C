// Before/after software-kernel comparison; no mesh, equilibrium solve or LIMS Job.
#include "AbaqusCDPSubstepIntegrator.h"
#include "CDPDiagnostics.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char ** argv)
{
  const unsigned int repeats = argc > 1 ? std::stoul(argv[1]) : 1;
  const std::string directory = "test/tests/cdp_material_table/data/";
  const double E = 3.04e10;
  const CDPMaterialTable table(directory + "compression_hardening.csv",
                               directory + "compression_damage.csv",
                               directory + "tension_stiffening.csv",
                               directory + "tension_damage.csv", E);
  const double ft = table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION, 0).stress.value;
  const double fc = table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::COMPRESSION, 0).stress.value;
  const std::array<AbaqusCDPLocalIntegrator::SymmetricTensor, 4> strains = {{
      {1.05*ft/E,-0.21*ft/E,-0.21*ft/E,0,0,0},
      {-1.05*fc/E,0.21*fc/E,0.21*fc/E,0,0,0},
      {1.08*ft/E,-0.216*ft/E,-0.216*ft/E,0.08*ft/E,0,0},
      {-1.05*fc/E,0.21*fc/E,0.21*fc/E,0.07*fc/E,-0.03*fc/E,0}}};
  std::cout << std::setprecision(17);
  unsigned int id=0;
  volatile double checksum=0;
  for (bool ad : {true,false})
  {
    AbaqusCDPLocalIntegrator::Parameters parameters{E,0.2,36.31,0.1,1.16,0.667,40,1e-9,1e-7,1e-6};
    parameters.use_automatic_differentiation_jacobian=ad;
    const AbaqusCDPLocalIntegrator local(table,parameters);
    const AbaqusCDPStateIntegrator state(local,{0.0,1.0,5e-4,1e-12});
    for (const auto & strain : strains)
      for (bool aggregate : {false,true})
      {
        ++id;
        const AbaqusCDPSubstepIntegrator integrator(state,{16,std::abs(strain[0])/3,1e-8,aggregate,aggregate});
        CDPDiagnostics::Counters costs={}; CDPDiagnostics::Context context;context.counters=&costs;
        const auto evaluate=[&](){return integrator.integrateLinearized({},strain,1e-3,{});};
        const auto result=[&](){CDPDiagnostics::Binding binding(&context);return evaluate();}();
        const auto emit=[&](const char *label,const auto & values){std::size_t i=0;for(const auto v:values)std::cout<<"value,"<<id<<','<<label<<','<<i++<<','<<v<<'\n';};
        emit("stress",result.result.final_result.cauchy_stress);
        emit("plastic",result.result.final_result.state.backbone.plastic_strain);
        emit("viscous_plastic",result.result.final_result.state.viscous_plastic_strain);
        const auto & s=result.result.final_result.state;
        emit("history",std::array<double,4>{s.backbone.tensile_equivalent_plastic_strain,s.backbone.compressive_equivalent_plastic_strain,s.viscous_tension_damage,s.viscous_compression_damage});
        for(std::size_t c=0;c<6;++c)emit(("tangent"+std::to_string(c)).c_str(),result.algorithmic_tangent[c]);
        const auto local_result=local.integrateLinearized(strain,{});
        for(std::size_t c=0;c<14;++c)emit(("transition"+std::to_string(c)).c_str(),local_result.derivative[c]);
        std::cout<<"value,"<<id<<",substeps,0,"<<result.result.accepted_substeps<<'\n';
        std::cout<<"cost,"<<id<<",spectrum,"<<costs[CDPDiagnostics::SPECTRUM].calls<<'\n';
        std::cout<<"cost,"<<id<<",backsolve,"<<costs[CDPDiagnostics::BACKSOLVE].calls<<'\n';
        const auto start=std::chrono::steady_clock::now();
        for(unsigned int n=0;n<repeats;++n)checksum+=evaluate().result.final_result.cauchy_stress[0];
        const double us=std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"timing,"<<id<<",microseconds,"<<us<<'\n';
      }
  }
  std::cerr<<"checksum="<<checksum<<'\n';
}

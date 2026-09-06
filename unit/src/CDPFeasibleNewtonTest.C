#include "AbaqusCDPLocalIntegrator.h"
#include "gtest/gtest.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
namespace {
CDPMaterialTable projectedTable() {
 const std::string d="test/tests/cdp_material_table/uniaxial_tension_20260830/";
 return CDPMaterialTable(d+"compression_hardening.csv",d+"compression_damage.csv",d+"tension_stiffening.csv",d+"tension_damage.csv",29791500000.0);
}
}
TEST(CDPFeasibleNewton, D01Capture0) {
 auto table=projectedTable();AbaqusCDPLocalIntegrator::Parameters p{29791500000.0,.2,36,.1,1.16,.667};
 AbaqusCDPLocalIntegrator::State old;
old.plastic_strain={-0.0001700362842915049,-0.00017003628428836188,0.00084379512253861754,1.216360837019724e-05,-1.7356925963159499e-05,-1.7356925977610017e-05};old.tensile_equivalent_plastic_strain=0.0008444160985219005;old.compressive_equivalent_plastic_strain=0;
AbaqusCDPLocalIntegrator::SymmetricTensor target={-9.0365343635225316e-05,-9.036534363203684e-05,0.00093614078695955288,1.2324873920195368e-05,-1.7592418891230148e-05,-1.7592418905879295e-05};
 const auto saved=old;
 EXPECT_THROW(AbaqusCDPLocalIntegrator(table,p).integrate(target,old),std::runtime_error);
 p.project_failed_newton_step=true;
 const auto result=AbaqusCDPLocalIntegrator(table,p).integrate(target,old);
 EXPECT_LE(result.residual_norm,p.residual_tolerance);EXPECT_LE(result.iterations,40u);
 EXPECT_GE(result.plastic_multiplier,0);EXPECT_GE(result.state.tensile_equivalent_plastic_strain,old.tensile_equivalent_plastic_strain);EXPECT_GE(result.state.compressive_equivalent_plastic_strain,old.compressive_equivalent_plastic_strain);
 const auto lin=AbaqusCDPLocalIntegrator(table,p).integrateLinearized(target,old);
 EXPECT_EQ(result.effective_stress,lin.result.effective_stress);
 for(const auto & column:lin.derivative)for(double x:column)EXPECT_TRUE(std::isfinite(x));
 // Tight tolerance isolates derivative accuracy from stopping error in this software test.
 auto derivative_parameters=p;derivative_parameters.residual_tolerance=1e-12;
 const AbaqusCDPLocalIntegrator derivative_core(table,derivative_parameters);
 const auto tangent=derivative_core.integrateLinearized(target,old);
 auto plus=target,minus=target;const double h=1e-10;plus[2]+=h;minus[2]-=h;
 const auto rp=derivative_core.integrate(plus,old),rm=derivative_core.integrate(minus,old);
 double derivative_error=0,derivative_scale=1;
 for(unsigned i=0;i<6;++i){const double fd=(rp.effective_stress[i]-rm.effective_stress[i])/(2*h);derivative_error=std::max(derivative_error,std::abs(fd-tangent.derivative[2][i]));derivative_scale=std::max(derivative_scale,std::abs(fd));}
 EXPECT_LT(derivative_error/derivative_scale,1e-3);
 EXPECT_EQ(old.plastic_strain,saved.plastic_strain);EXPECT_DOUBLE_EQ(old.tensile_equivalent_plastic_strain,saved.tensile_equivalent_plastic_strain);EXPECT_DOUBLE_EQ(old.compressive_equivalent_plastic_strain,saved.compressive_equivalent_plastic_strain);
}
TEST(CDPFeasibleNewton, D02Capture0) {
 auto table=projectedTable();AbaqusCDPLocalIntegrator::Parameters p{29791500000.0,.2,36,.1,1.16,.667};
 AbaqusCDPLocalIntegrator::State old;
old.plastic_strain={0.0014872687374944899,0.0014872687371291821,-0.0014620339967781621,-0.00020665604500810279,0.00016646232338713882,0.00016646232341545116};old.tensile_equivalent_plastic_strain=0.00012298361418062488;old.compressive_equivalent_plastic_strain=0.0014613947559259319;
AbaqusCDPLocalIntegrator::SymmetricTensor target={0.0015240141114345868,0.0015240141110805254,-0.0013975892305307224,-0.00021690460335515682,0.00017610941571070507,0.0001761094157224077};
 const auto saved=old;
 EXPECT_THROW(AbaqusCDPLocalIntegrator(table,p).integrate(target,old),std::runtime_error);
 p.project_failed_newton_step=true;
 const auto result=AbaqusCDPLocalIntegrator(table,p).integrate(target,old);
 EXPECT_LE(result.residual_norm,p.residual_tolerance);EXPECT_LE(result.iterations,40u);
 EXPECT_GE(result.plastic_multiplier,0);EXPECT_GE(result.state.tensile_equivalent_plastic_strain,old.tensile_equivalent_plastic_strain);EXPECT_GE(result.state.compressive_equivalent_plastic_strain,old.compressive_equivalent_plastic_strain);
 const auto lin=AbaqusCDPLocalIntegrator(table,p).integrateLinearized(target,old);
 EXPECT_EQ(result.effective_stress,lin.result.effective_stress);
 for(const auto & column:lin.derivative)for(double x:column)EXPECT_TRUE(std::isfinite(x));
 // Tight tolerance isolates derivative accuracy from stopping error in this software test.
 auto derivative_parameters=p;derivative_parameters.residual_tolerance=1e-12;
 const AbaqusCDPLocalIntegrator derivative_core(table,derivative_parameters);
 const auto tangent=derivative_core.integrateLinearized(target,old);
 auto plus=target,minus=target;const double h=1e-10;plus[2]+=h;minus[2]-=h;
 const auto rp=derivative_core.integrate(plus,old),rm=derivative_core.integrate(minus,old);
 double derivative_error=0,derivative_scale=1;
 for(unsigned i=0;i<6;++i){const double fd=(rp.effective_stress[i]-rm.effective_stress[i])/(2*h);derivative_error=std::max(derivative_error,std::abs(fd-tangent.derivative[2][i]));derivative_scale=std::max(derivative_scale,std::abs(fd));}
 EXPECT_LT(derivative_error/derivative_scale,1e-3);
 EXPECT_EQ(old.plastic_strain,saved.plastic_strain);EXPECT_DOUBLE_EQ(old.tensile_equivalent_plastic_strain,saved.tensile_equivalent_plastic_strain);EXPECT_DOUBLE_EQ(old.compressive_equivalent_plastic_strain,saved.compressive_equivalent_plastic_strain);
}
TEST(CDPFeasibleNewton, D02Capture1) {
 auto table=projectedTable();AbaqusCDPLocalIntegrator::Parameters p{29791500000.0,.2,36,.1,1.16,.667};
 AbaqusCDPLocalIntegrator::State old;
old.plastic_strain={0.00089439568582129524,0.00089439568360668502,8.5302688470749991e-05,0.00061818102327210037,-0.0013700165132083122,-0.0013700165147215513};old.tensile_equivalent_plastic_strain=0.002195768637365299;old.compressive_equivalent_plastic_strain=0.0004327354920335101;
AbaqusCDPLocalIntegrator::SymmetricTensor target={0.0009531850571030951,0.0009531850548886397,0.00015815542001929389,0.00061822768025922,-0.0013701250945158156,-0.0013701250960286499};
 const auto saved=old;
 EXPECT_THROW(AbaqusCDPLocalIntegrator(table,p).integrate(target,old),std::runtime_error);
 p.project_failed_newton_step=true;
 const auto result=AbaqusCDPLocalIntegrator(table,p).integrate(target,old);
 EXPECT_LE(result.residual_norm,p.residual_tolerance);EXPECT_LE(result.iterations,40u);
 EXPECT_GE(result.plastic_multiplier,0);EXPECT_GE(result.state.tensile_equivalent_plastic_strain,old.tensile_equivalent_plastic_strain);EXPECT_GE(result.state.compressive_equivalent_plastic_strain,old.compressive_equivalent_plastic_strain);
 const auto lin=AbaqusCDPLocalIntegrator(table,p).integrateLinearized(target,old);
 EXPECT_EQ(result.effective_stress,lin.result.effective_stress);
 for(const auto & column:lin.derivative)for(double x:column)EXPECT_TRUE(std::isfinite(x));
 // Tight tolerance isolates derivative accuracy from stopping error in this software test.
 auto derivative_parameters=p;derivative_parameters.residual_tolerance=1e-12;
 const AbaqusCDPLocalIntegrator derivative_core(table,derivative_parameters);
 const auto tangent=derivative_core.integrateLinearized(target,old);
 auto plus=target,minus=target;const double h=1e-10;plus[2]+=h;minus[2]-=h;
 const auto rp=derivative_core.integrate(plus,old),rm=derivative_core.integrate(minus,old);
 double derivative_error=0,derivative_scale=1;
 for(unsigned i=0;i<6;++i){const double fd=(rp.effective_stress[i]-rm.effective_stress[i])/(2*h);derivative_error=std::max(derivative_error,std::abs(fd-tangent.derivative[2][i]));derivative_scale=std::max(derivative_scale,std::abs(fd));}
 EXPECT_LT(derivative_error/derivative_scale,1e-3);
 EXPECT_EQ(old.plastic_strain,saved.plastic_strain);EXPECT_DOUBLE_EQ(old.tensile_equivalent_plastic_strain,saved.tensile_equivalent_plastic_strain);EXPECT_DOUBLE_EQ(old.compressive_equivalent_plastic_strain,saved.compressive_equivalent_plastic_strain);
}
TEST(CDPFeasibleNewton, D02Capture2) {
 auto table=projectedTable();AbaqusCDPLocalIntegrator::Parameters p{29791500000.0,.2,36,.1,1.16,.667};
 AbaqusCDPLocalIntegrator::State old;
old.plastic_strain={0.00089792970657001636,0.00089792969212419726,0.00011291180885880218,0.00062554735417117379,-0.0013831757902785732,-0.0013831757919230423};old.tensile_equivalent_plastic_strain=0.002249274742578067;old.compressive_equivalent_plastic_strain=0.0004327354920335101;
AbaqusCDPLocalIntegrator::SymmetricTensor target={0.00096133133028604694,0.00096133131583884689,0.00018357006345837883,0.00062554901815871085,-0.0013831792766431663,-0.0013831792782876442};
 const auto saved=old;
 EXPECT_THROW(AbaqusCDPLocalIntegrator(table,p).integrate(target,old),std::runtime_error);
 p.project_failed_newton_step=true;
 const auto result=AbaqusCDPLocalIntegrator(table,p).integrate(target,old);
 EXPECT_LE(result.residual_norm,p.residual_tolerance);EXPECT_LE(result.iterations,40u);
 EXPECT_GE(result.plastic_multiplier,0);EXPECT_GE(result.state.tensile_equivalent_plastic_strain,old.tensile_equivalent_plastic_strain);EXPECT_GE(result.state.compressive_equivalent_plastic_strain,old.compressive_equivalent_plastic_strain);
 const auto lin=AbaqusCDPLocalIntegrator(table,p).integrateLinearized(target,old);
 EXPECT_EQ(result.effective_stress,lin.result.effective_stress);
 for(const auto & column:lin.derivative)for(double x:column)EXPECT_TRUE(std::isfinite(x));
 // Tight tolerance isolates derivative accuracy from stopping error in this software test.
 auto derivative_parameters=p;derivative_parameters.residual_tolerance=1e-12;
 const AbaqusCDPLocalIntegrator derivative_core(table,derivative_parameters);
 const auto tangent=derivative_core.integrateLinearized(target,old);
 auto plus=target,minus=target;const double h=1e-10;plus[2]+=h;minus[2]-=h;
 const auto rp=derivative_core.integrate(plus,old),rm=derivative_core.integrate(minus,old);
 double derivative_error=0,derivative_scale=1;
 for(unsigned i=0;i<6;++i){const double fd=(rp.effective_stress[i]-rm.effective_stress[i])/(2*h);derivative_error=std::max(derivative_error,std::abs(fd-tangent.derivative[2][i]));derivative_scale=std::max(derivative_scale,std::abs(fd));}
 EXPECT_LT(derivative_error/derivative_scale,1e-3);
 EXPECT_EQ(old.plastic_strain,saved.plastic_strain);EXPECT_DOUBLE_EQ(old.tensile_equivalent_plastic_strain,saved.tensile_equivalent_plastic_strain);EXPECT_DOUBLE_EQ(old.compressive_equivalent_plastic_strain,saved.compressive_equivalent_plastic_strain);
}
TEST(CDPFeasibleNewton, D02Capture3) {
 auto table=projectedTable();AbaqusCDPLocalIntegrator::Parameters p{29791500000.0,.2,36,.1,1.16,.667};
 AbaqusCDPLocalIntegrator::State old;
old.plastic_strain={0.00090746960819255583,0.00090746987188716798,0.00012445839532105251,0.00063690805969194729,-0.0014009128986647683,-0.0014009126611421736};old.tensile_equivalent_plastic_strain=0.0022945280968675206;old.compressive_equivalent_plastic_strain=0.0004327354920335101;
AbaqusCDPLocalIntegrator::SymmetricTensor target={0.000972388249845218,0.00097238851354146155,0.00019442585526969455,0.0006369081687145761,-0.0014009130930480956,-0.0014009128555236452};
 const auto saved=old;
 EXPECT_THROW(AbaqusCDPLocalIntegrator(table,p).integrate(target,old),std::runtime_error);
 p.project_failed_newton_step=true;
 const auto result=AbaqusCDPLocalIntegrator(table,p).integrate(target,old);
 EXPECT_LE(result.residual_norm,p.residual_tolerance);EXPECT_LE(result.iterations,40u);
 EXPECT_GE(result.plastic_multiplier,0);EXPECT_GE(result.state.tensile_equivalent_plastic_strain,old.tensile_equivalent_plastic_strain);EXPECT_GE(result.state.compressive_equivalent_plastic_strain,old.compressive_equivalent_plastic_strain);
 const auto lin=AbaqusCDPLocalIntegrator(table,p).integrateLinearized(target,old);
 EXPECT_EQ(result.effective_stress,lin.result.effective_stress);
 for(const auto & column:lin.derivative)for(double x:column)EXPECT_TRUE(std::isfinite(x));
 // Tight tolerance isolates derivative accuracy from stopping error in this software test.
 auto derivative_parameters=p;derivative_parameters.residual_tolerance=1e-12;
 const AbaqusCDPLocalIntegrator derivative_core(table,derivative_parameters);
 const auto tangent=derivative_core.integrateLinearized(target,old);
 auto plus=target,minus=target;const double h=1e-10;plus[2]+=h;minus[2]-=h;
 const auto rp=derivative_core.integrate(plus,old),rm=derivative_core.integrate(minus,old);
 double derivative_error=0,derivative_scale=1;
 for(unsigned i=0;i<6;++i){const double fd=(rp.effective_stress[i]-rm.effective_stress[i])/(2*h);derivative_error=std::max(derivative_error,std::abs(fd-tangent.derivative[2][i]));derivative_scale=std::max(derivative_scale,std::abs(fd));}
 EXPECT_LT(derivative_error/derivative_scale,1e-3);
 EXPECT_EQ(old.plastic_strain,saved.plastic_strain);EXPECT_DOUBLE_EQ(old.tensile_equivalent_plastic_strain,saved.tensile_equivalent_plastic_strain);EXPECT_DOUBLE_EQ(old.compressive_equivalent_plastic_strain,saved.compressive_equivalent_plastic_strain);
}
TEST(CDPFeasibleNewton, OriginalSuccessfulStatesAndTangentsUnchanged) {
 auto table=projectedTable();AbaqusCDPLocalIntegrator::Parameters p{29791500000.0,.2,36,.1,1.16,.667};auto q=p;q.project_failed_newton_step=true;
 const double ft=table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::TENSION,0).stress.value;
 const double fc=table.responseByEquivalentPlasticStrain(CDPMaterialTable::Branch::COMPRESSION,0).stress.value;
 for(const double stress : {1000.0,1.01*ft,-1.01*fc}){
 const double e=stress/p.youngs_modulus;const AbaqusCDPLocalIntegrator::SymmetricTensor strain={-.2*e,-.2*e,e,0,0,0};
 const auto base=AbaqusCDPLocalIntegrator(table,p).integrateLinearized(strain,{});const auto result=AbaqusCDPLocalIntegrator(table,q).integrateLinearized(strain,{});
 EXPECT_EQ(base.result.effective_stress,result.result.effective_stress);EXPECT_EQ(base.result.state.plastic_strain,result.result.state.plastic_strain);EXPECT_DOUBLE_EQ(base.result.state.tensile_equivalent_plastic_strain,result.result.state.tensile_equivalent_plastic_strain);EXPECT_DOUBLE_EQ(base.result.state.compressive_equivalent_plastic_strain,result.result.state.compressive_equivalent_plastic_strain);EXPECT_EQ(base.derivative,result.derivative);
 }
}

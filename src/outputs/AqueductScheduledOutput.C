#include "AqueductScheduledOutput.h"
#include "AqueductOutputSchedule.h"
#include "MooseApp.h"

registerMooseObject("DamSafetyApp",AqueductScheduledExodus);
registerMooseObject("DamSafetyApp",AqueductScheduledCSV);

template <typename Base>
InputParameters AqueductScheduledOutput<Base>::validParams()
{
  auto p=Base::validParams();
  p.template addRequiredParam<Real>("phase_origin","Total time of the ELC origin");
  p.template addRequiredParam<Real>("phase_end","Total time of the ELC endpoint");
  p.template addRequiredRangeCheckedParam<Real>("mark_interval","mark_interval>0","Source output interval");
  p.template addRangeCheckedParam<Real>("mark_tolerance",1e-10,"mark_tolerance>=0",
                                       "Rounding tolerance for output marks only; never a solve tolerance");
  return p;
}
template <typename Base>
AqueductScheduledOutput<Base>::AqueductScheduledOutput(const InputParameters & p)
 :Base(p),_origin(this->template getParam<Real>("phase_origin")),
 _end(this->template getParam<Real>("phase_end")),
 _period(this->template getParam<Real>("mark_interval")),
 _mark_tolerance(this->template getParam<Real>("mark_tolerance")),
 _last_mark(this->template declareRestartableData<long long>("aq_last_output_mark",-1))
{
  AqueductOutputSchedule::slot(_origin,_origin,_end,_period,_mark_tolerance);
  if(this->template getParam<Real>("min_simulation_time_interval")!=0. ||
     this->template getParam<unsigned int>("time_step_interval")!=1 ||
     !this->template getParam<std::vector<Real>>("sync_times").empty() ||
     this->isParamValid("sync_times_object"))
    this->paramError("mark_interval","Do not combine mark scheduling with other interval or sync controls");
}
template <typename Base>
void AqueductScheduledOutput<Base>::outputStep(const ExecFlagType & type)
{
  if((!this->_allow_output && type!=EXEC_FORCED) ||
     (type==EXEC_INITIAL && this->_app.isRecovering()) ||
     (!this->_execute_on.isValueSet(type) && type!=EXEC_FORCED))return;
  const auto mark=AqueductOutputSchedule::slot(this->_time,_origin,_end,_period,_mark_tolerance);
  if(mark<0 || mark<=_last_mark)return;
  Base::outputStep(type);
  _last_mark=mark;
}
template class AqueductScheduledOutput<Exodus>;
template class AqueductScheduledOutput<CSV>;

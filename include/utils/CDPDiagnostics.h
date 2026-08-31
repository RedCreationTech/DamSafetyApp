#pragma once

#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>

/** Optional thread-local diagnostics. Never owns or changes constitutive state. */
namespace CDPDiagnostics
{
enum Kind { MATERIAL, PARTITION, LOCAL, NEWTON_STEP, LOCAL_LINEARIZED, STATE, STATE_LINEARIZED,
            RESIDUAL, AD_JACOBIAN, FD_JACOBIAN, SPECTRUM, AD_SPECTRUM, FACTOR,
            BACKSOLVE, STATE_CHAIN, TRACE_IO, RECOMPUTE, COUNT };
inline constexpr std::array<const char *, COUNT> names = {
    "material", "partition", "local", "newton_step", "local_linearized", "state", "state_linearized",
    "residual", "ad_jacobian", "fd_jacobian", "spectrum", "ad_spectrum", "factor",
    "backsolve", "state_chain", "trace_io", "diagnostic_recompute"};
struct Counter { std::uint64_t calls=0, failed=0; double microseconds=0, failed_microseconds=0; };
using Counters = std::array<Counter, COUNT>;
struct Context
{
  Counters * counters=nullptr;
  std::ostream * stream=nullptr;
  unsigned int * failure_samples=nullptr;
  unsigned int maximum_failure_samples=4;
  bool trace=false;
  double time=0, dt=0;
  std::uint64_t step=0, element=0, qp=0, call=0;
  unsigned int partition=0, substep=0;
};
extern thread_local Context * current;
extern thread_local double excluded_microseconds;
class Binding
{
public:
  explicit Binding(Context * context) : _previous(current), _excluded(current && !context)
  { if(_excluded) _start=std::chrono::steady_clock::now(); current=context; }
  ~Binding()
  {
    if(_excluded)
    {
      const double us=std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-_start).count();
      excluded_microseconds+=us;
      if(_previous->counters)
      { auto & c=(*_previous->counters)[RECOMPUTE];++c.calls;c.microseconds+=us; }
    }
    current=_previous;
  }
  Binding(const Binding &)=delete;
private:
  Context * _previous;
  bool _excluded;
  std::chrono::steady_clock::time_point _start;
};
class Scope
{
public:
  explicit Scope(Kind kind)
    : _counter(current && current->counters ? &(*current->counters)[kind] : nullptr),
      _exceptions(std::uncaught_exceptions())
  { if (_counter) { ++_counter->calls; _start=std::chrono::steady_clock::now(); _excluded_start=excluded_microseconds; } }
  ~Scope()
  {
    if (!_counter) return;
    const double elapsed=std::chrono::duration<double,std::micro>(std::chrono::steady_clock::now()-_start).count();
    const double us=std::max(0.0,elapsed-(excluded_microseconds-_excluded_start));
    _counter->microseconds+=us;
    if (_failed || std::uncaught_exceptions()>_exceptions) { ++_counter->failed; _counter->failed_microseconds+=us; }
  }
  void failed() { _failed=true; }
  Scope(const Scope &)=delete;
private:
  Counter * _counter;
  int _exceptions;
  bool _failed=false;
  double _excluded_start=0;
  std::chrono::steady_clock::time_point _start;
};
inline void json(std::ostream & out, double value)
{ if (std::isfinite(value)) out<<std::setprecision(17)<<value; else out<<"null"; }
template<class T, std::size_t N> void json(std::ostream & out,const std::array<T,N> & values)
{
  out<<'[';
  for (std::size_t i=0;i<N;++i) { if(i) out<<','; json(out,values[i]); }
  out<<']';
}
template<class T> void stateJson(std::ostream & out,const T & state)
{
  out<<"{\"plastic\":";json(out,state.backbone.plastic_strain);
  out<<",\"kappa_t\":";json(out,state.backbone.tensile_equivalent_plastic_strain);
  out<<",\"kappa_c\":";json(out,state.backbone.compressive_equivalent_plastic_strain);
  out<<",\"viscous_plastic\":";json(out,state.viscous_plastic_strain);
  out<<",\"damage_t\":";json(out,state.viscous_tension_damage);
  out<<",\"damage_c\":";json(out,state.viscous_compression_damage);out<<'}';
}
inline void event(const std::string & kind,const std::string & payload)
{
  if (!current || !current->stream) return;
  Scope timer(TRACE_IO);
  auto & out=*current->stream;
  out<<std::setprecision(17)<<"{\"event\":\""<<kind<<"\",\"time\":"<<current->time
     <<",\"dt\":"<<current->dt<<",\"step\":"<<current->step<<",\"element\":"<<current->element
     <<",\"qp\":"<<current->qp<<",\"call\":"<<current->call<<",\"partition\":"<<current->partition
     <<",\"substep\":"<<current->substep<<','<<payload<<"}\n";
  out.flush();
}
inline bool sampleFailure()
{
  if (!current || !current->stream || !current->failure_samples ||
      *current->failure_samples>=current->maximum_failure_samples) return false;
  ++*current->failure_samples;
  return true;
}
}

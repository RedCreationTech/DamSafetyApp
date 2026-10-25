#include "AbaqusCDPPointPacket.h"

#include <cctype>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace
{
const char * names[4] = {"compression_hardening", "compression_damage",
                         "tension_stiffening", "tension_damage"};
const char * headers[4] = {"stress_pa,inelastic_strain", "damage_c,inelastic_strain",
                           "stress_pa,cracking_strain", "damage_t,cracking_strain"};

void bits(std::ostream & out, double value)
{
  static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559,
                "CDP point packets require IEEE-754 binary64");
  std::uint64_t raw;
  std::memcpy(&raw, &value, sizeof(raw));
  std::ostringstream token;
  token << std::hex << std::setfill('0') << std::setw(16) << raw;
  out << ' ' << token.str();
}
double number(std::istream & input)
{
  std::string token;
  if (!(input >> token) || token.size()!=16)
    throw std::runtime_error("Invalid CDP packet binary64 token");
  for (unsigned char c : token)
    if (!std::isxdigit(c))
      throw std::runtime_error("Invalid CDP packet binary64 token");
  const auto raw = static_cast<std::uint64_t>(std::stoull(token, nullptr, 16));
  double value;
  std::memcpy(&value, &raw, sizeof(value));
  return value;
}
template <typename UInt> UInt integer(std::istream & input)
{
  std::string token;
  if (!(input >> token) || token.empty())
    throw std::runtime_error("Missing CDP packet unsigned integer");
  for (unsigned char c : token)
    if (c < '0' || c > '9')
      throw std::runtime_error("Invalid CDP packet unsigned integer");
  const auto value = std::stoull(token);
  if (value > std::numeric_limits<UInt>::max())
    throw std::runtime_error("CDP packet unsigned integer out of range");
  return static_cast<UInt>(value);
}
void expect(std::istream & input, const char * key)
{
  std::string actual;
  if (!(input >> actual) || actual!=key)
    throw std::runtime_error(std::string("Missing CDP packet field: ")+key);
}
template <typename Array> void array(std::ostream & out, const char * key, const Array & values)
{
  out << key;
  for (double value : values) bits(out, value);
  out << '\n';
}
template <typename Array> void array(std::istream & input, const char * key, Array & values)
{
  expect(input,key);
  for (double & value : values) value=number(input);
}
}

void AbaqusCDPPointPacket::setTables(const CDPMaterialTable & table)
{
  for (unsigned int i=0;i<4;++i)
  {
    const auto branch=i<2 ? CDPMaterialTable::Branch::COMPRESSION : CDPMaterialTable::Branch::TENSION;
    const auto & x=i%2 ? table.damageAbscissa(branch) : table.stressAbscissa(branch);
    const auto & y=i%2 ? table.damageValues(branch) : table.stressValues(branch);
    tables[i].clear();
    for (std::size_t j=0;j<x.size();++j) tables[i].emplace_back(y[j],x[j]);
  }
}
void AbaqusCDPPointPacket::write(std::ostream & out) const
{
  out << "AQ_CDP_POINT_V1\ncontext " << element << ' ' << qp << ' ' << side << ' ' << rank << ' '
      << thread << ' ' << step;
  bits(out,time);bits(out,dt);for(double v:position)bits(out,v);out << '\n';
  out << "material " << std::quoted(material) << '\n';
  out << "local";
  for(double v:{local.youngs_modulus,local.poissons_ratio,local.dilation_angle_degrees,
                local.eccentricity,local.biaxial_to_uniaxial_compression_ratio,
                local.tensile_meridian_ratio})bits(out,v);
  out << ' ' << local.maximum_iterations;
  for(double v:{local.residual_tolerance,local.finite_difference_step,local.minimum_line_search})bits(out,v);
  out << ' ' << static_cast<unsigned int>(local.use_automatic_differentiation_jacobian) << '\n';
  out << "state";
  for(double v:{state.tension_recovery,state.compression_recovery,state.relaxation_time,
                state.state_tolerance})bits(out,v);
  out << "\nsubsteps " << substeps.maximum_substeps;
  bits(out,substeps.maximum_strain_increment);bits(out,substeps.tangent_perturbation);out << '\n';
  array(out,"old_strain",old_strain);array(out,"new_strain",new_strain);
  array(out,"old_backbone",old.backbone.plastic_strain);
  out << "old_kappa";bits(out,old.backbone.tensile_equivalent_plastic_strain);
  bits(out,old.backbone.compressive_equivalent_plastic_strain);out << '\n';
  array(out,"old_viscous",old.viscous_plastic_strain);
  out << "old_damage";bits(out,old.viscous_tension_damage);bits(out,old.viscous_compression_damage);
  out << '\n';
  for(unsigned int i=0;i<4;++i)
  {
    out << names[i] << ' ' << tables[i].size() << '\n';
    for(const auto & row:tables[i]){bits(out,row.first);bits(out,row.second);out << '\n';}
  }
  out << "expected_error " << std::quoted(expected_error) << "\nEND\n";
  if(!out)throw std::runtime_error("CDP point packet write failed");
}
AbaqusCDPPointPacket AbaqusCDPPointPacket::read(std::istream & input)
{
  AbaqusCDPPointPacket p;
  expect(input,"AQ_CDP_POINT_V1");expect(input,"context");
  p.element=integer<std::uint64_t>(input);p.qp=integer<unsigned int>(input);
  p.side=integer<unsigned int>(input);p.rank=integer<unsigned int>(input);
  p.thread=integer<unsigned int>(input);p.step=integer<std::uint64_t>(input);
  p.time=number(input);p.dt=number(input);for(double & v:p.position)v=number(input);
  expect(input,"material");input>>std::quoted(p.material);
  expect(input,"local");
  p.local.youngs_modulus=number(input);p.local.poissons_ratio=number(input);
  p.local.dilation_angle_degrees=number(input);p.local.eccentricity=number(input);
  p.local.biaxial_to_uniaxial_compression_ratio=number(input);p.local.tensile_meridian_ratio=number(input);
  p.local.maximum_iterations=integer<unsigned int>(input);p.local.residual_tolerance=number(input);
  p.local.finite_difference_step=number(input);p.local.minimum_line_search=number(input);
  const auto ad=integer<unsigned int>(input);
  if(ad>1)throw std::runtime_error("Invalid CDP packet AD flag");
  p.local.use_automatic_differentiation_jacobian=ad;
  expect(input,"state");p.state.tension_recovery=number(input);p.state.compression_recovery=number(input);
  p.state.relaxation_time=number(input);p.state.state_tolerance=number(input);
  expect(input,"substeps");p.substeps.maximum_substeps=integer<unsigned int>(input);
  p.substeps.maximum_strain_increment=number(input);p.substeps.tangent_perturbation=number(input);
  array(input,"old_strain",p.old_strain);array(input,"new_strain",p.new_strain);
  array(input,"old_backbone",p.old.backbone.plastic_strain);
  expect(input,"old_kappa");p.old.backbone.tensile_equivalent_plastic_strain=number(input);
  p.old.backbone.compressive_equivalent_plastic_strain=number(input);
  array(input,"old_viscous",p.old.viscous_plastic_strain);
  expect(input,"old_damage");p.old.viscous_tension_damage=number(input);
  p.old.viscous_compression_damage=number(input);
  for(unsigned int i=0;i<4;++i)
  {
    expect(input,names[i]);const auto size=integer<std::size_t>(input);
    if(size<2 || size>100000)throw std::runtime_error("Invalid CDP packet table size");
    for(std::size_t j=0;j<size;++j)
    {
      const double value=number(input),abscissa=number(input);
      p.tables[i].emplace_back(value,abscissa);
    }
  }
  expect(input,"expected_error");input>>std::quoted(p.expected_error);expect(input,"END");
  std::string trailing;
  if(input>>trailing)throw std::runtime_error("Trailing CDP packet data");
  return p;
}
std::array<std::string,4> AbaqusCDPPointPacket::writeTables(const std::string & directory) const
{
  std::filesystem::create_directories(directory);
  std::array<std::string,4> paths;
  for(unsigned int i=0;i<4;++i)
  {
    paths[i]=(std::filesystem::path(directory)/(std::string(names[i])+".csv")).string();
    std::ostringstream content;content.imbue(std::locale::classic());
    content << headers[i] << '\n' << std::setprecision(17);
    for(const auto & row:tables[i])content << row.first << ',' << row.second << '\n';
    if(std::filesystem::exists(paths[i]))
    {
      std::ifstream existing(paths[i]);std::ostringstream bytes;bytes<<existing.rdbuf();
      if(bytes.str()!=content.str())throw std::runtime_error("Refusing to overwrite different replay table");
    }
    else
    {
      std::ofstream out(paths[i]);out<<content.str();
      if(!out)throw std::runtime_error("Replay table write failed");
    }
  }
  return paths;
}

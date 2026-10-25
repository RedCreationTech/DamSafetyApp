#include "AbaqusCDPPointPacket.h"
#include "gtest/gtest.h"

#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

TEST(AbaqusCDPPointPacket, LosslessRoundTripIncludingSignedZeroAndAllTables)
{
  AbaqusCDPPointPacket p;p.element=19;p.step=101;p.qp=2;p.time=1.2113125;p.dt=.0025;
  p.old_strain[0]=-0.;p.new_strain[5]=std::nextafter(1.,2.);
  p.local={30e9,.2,36.31,.1,1.16,.6667,40,1e-9,1e-7,1e-6};
  p.state={0.,1.,.05,1e-12};p.substeps={256,2.5e-5,1e-8};
  p.material="test_material";p.expected_error="example: \"quote\" and \\ slash";
  for(auto & table:p.tables)table={{0.,0.},{1.,.1}};
  std::ostringstream bytes;p.write(bytes);std::istringstream input(bytes.str());
  const auto q=AbaqusCDPPointPacket::read(input);std::ostringstream again;q.write(again);
  EXPECT_EQ(bytes.str(),again.str());EXPECT_TRUE(std::signbit(q.old_strain[0]));
  EXPECT_EQ(q.expected_error,p.expected_error);
}
TEST(AbaqusCDPPointPacket, RejectsSignedAndOverflowingContextIdentifiers)
{
  AbaqusCDPPointPacket p;p.element=19;
  for(auto & table:p.tables)table={{0.,0.},{1.,.1}};
  std::ostringstream bytes;p.write(bytes);
  for(const std::string context:{"context -1", "context +1", "context 18446744073709551616",
                                 "context 19 4294967296"})
  {
    auto bad=bytes.str();
    bad.replace(bad.find("context 19"),std::string("context 19").size(),context);
    std::istringstream input(bad);
    EXPECT_THROW(AbaqusCDPPointPacket::read(input),std::exception)<<context;
  }
}
TEST(AbaqusCDPPointPacket, RejectsTruncatedWrongVersionAndTrailingData)
{
  for(const auto data:{"AQ_CDP_POINT_V2", "AQ_CDP_POINT_V1\ncontext 1", "bad"})
  {std::istringstream input(data);EXPECT_THROW(AbaqusCDPPointPacket::read(input),std::runtime_error);}
  AbaqusCDPPointPacket p;for(auto & t:p.tables)t={{0.,0.},{1.,.1}};
  std::ostringstream bytes;p.write(bytes);std::istringstream input(bytes.str()+"unexpected");
  EXPECT_THROW(AbaqusCDPPointPacket::read(input),std::runtime_error);
}

// Explicit opt-in: no made-up failure state and no silent mechanics surrogate.
TEST(AbaqusCDPPointReplay, ExternalCapturedFailure)
{
  const char * capture=std::getenv("CDP_REPLAY_PACKET");
  if(!capture)GTEST_SKIP()<<"No real captured packet provided";
  const char * directory=std::getenv("CDP_REPLAY_OUTPUT");
  ASSERT_NE(directory,nullptr)<<"An isolated replay output directory is required";
  ASSERT_FALSE(std::string(directory).empty());
  ASSERT_FALSE(std::filesystem::exists(directory))<<"Refusing to overwrite replay evidence";
  std::ifstream input(capture);ASSERT_TRUE(input.good());
  const auto p=AbaqusCDPPointPacket::read(input);
  ASSERT_FALSE(p.expected_error.empty());
  const auto tables=p.writeTables((std::filesystem::path(directory)/"tables").string());
  CDPMaterialTable table(tables[0],tables[1],tables[2],tables[3],p.local.youngs_modulus);
  std::ofstream trace_file((std::filesystem::path(directory)/"local-trace.csv").string());
  ASSERT_TRUE(trace_file.good());
  AbaqusCDPLocalIntegrator::DiagnosticTrace trace;trace.stream=&trace_file;
  std::string plain_error,traced_error;
  // Both runs start from the packet's immutable old state. The trace cannot seed the second run.
  for(bool traced:{false,true})
  {
    AbaqusCDPLocalIntegrator local(table,p.local,traced?&trace:nullptr);
    AbaqusCDPStateIntegrator state(local,p.state);
    AbaqusCDPSubstepIntegrator substep(state,p.substeps);
    try{substep.integrateLinearized(p.old_strain,p.new_strain,p.dt,p.old);}
    catch(const std::exception & e){(traced?traced_error:plain_error)=e.what();}
  }
  trace_file.flush();
  std::ofstream receipt((std::filesystem::path(directory)/"replay-receipt.txt").string());
  receipt<<"expected "<<std::quoted(p.expected_error)<<"\nplain "<<std::quoted(plain_error)
      <<"\ntraced "<<std::quoted(traced_error)<<"\ntrace_written "<<trace.written
      <<"\ntrace_dropped "<<trace.dropped<<"\ntrace_io_failed "<<trace.io_failed<<'\n';
  EXPECT_EQ(plain_error,p.expected_error)<<"Exact failure reproduction is a gate, not an assumption";
  EXPECT_EQ(traced_error,plain_error);EXPECT_FALSE(trace.io_failed);EXPECT_TRUE(trace_file.good());
  EXPECT_EQ(trace.dropped,0u)<<"Trace bound reached: failure reproduction may match, but diagnosis is incomplete";
  EXPECT_GT(trace.written,0u);EXPECT_TRUE(receipt.good());
}

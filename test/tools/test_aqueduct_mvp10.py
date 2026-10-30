import math
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'tools/abaqus'))
from build_aqueduct_mvp10 import parse_blocks,render,compose,build

ROOT=Path(__file__).resolve().parents[2]

def fixture():
    return parse_blocks('''[Executioner]
type = Transient
nl_abs_tol = 1e-5
nl_rel_tol = 1e-8
nl_max_its = 35
[TimeIntegrator]
type = AqueductHHT
hht_start_time = 1
inactive_tsteps = 1
[]
[TimeStepper]
type = FunctionDT
function = dt_clock
[]
[]
[Functions/dynamic_on]
type = ParsedFunction
expression = 't>1'
[]
[Functions/dt_clock]
type = ParsedFunction
expression = 'if(t<1,1,.0025)'
[]
[Functions/motion_ax]
type = IntegratedAcceleration
start_time = 1
[]
[Physics/SolidMechanics/QuasiStatic/concrete]
block = '2 3'
[]
[Kernels/mass_ux]
type = AqueductHHTInertia
hht_start_time = 1
[]
[Kernels/shell_mass_rx]
type = AqueductFlatShellInertia
hht_start_time = 1
[]
[AuxKernels/accel_x]
type = NewmarkAccelAux
[]
[AuxKernels/velocity_x]
type = NewmarkVelAux
[]
[Materials/cdp_update]
type = AbaqusCDPStressUpdate
viscosity = .05
maximum_substeps = 256
[]
[Constraints/contact_law_1]
type = AqueductElasticSlipConstraint
mu = .55
critical_slip = .00025331
[]
[Outputs]
csv = true
[]
[Outputs/fields]
type = Exodus
min_simulation_time_interval = .05
[]
''')

class MVP10Tests(unittest.TestCase):
    def compose(self):
        b=fixture();before={k:dict(v) for k,v in b.items()}
        stats=compose(b,[dict(name='P',libmesh_node='42')],[dict(exodus_element='100')],
                      '[Kernels/gravity]\ntype = Gravity\nvalue = -9.8\n[]\n')
        return before,b,stats

    def test_strict_parser_roundtrip_and_rejects_silent_loss(self):
        b=fixture();self.assertEqual(parse_blocks(render(b)),b)
        for text in ['[]','[A]\nx=1\nx=2\n[]','[A]\nx=1',
                     '[A]\nx=1\n[]\n[A]\nx=2\n[]','x=1']:
            with self.assertRaises(ValueError):parse_blocks(text)

    def test_continuous_state_handoff_is_explicit_in_every_rate_path(self):
        _,b,_=self.compose()
        self.assertEqual(b['Executioner']['end_time'],'12')
        self.assertEqual(b['Executioner/TimeIntegrator']['inactive_tsteps'],'0')
        self.assertEqual(b['Executioner/TimeIntegrator']['quasistatic_initialization'],'true')
        for path in ['Kernels/mass_ux','Kernels/shell_mass_rx','Executioner/TimeIntegrator']:
            self.assertEqual(b[path]['hht_start_time'],'2')
        self.assertEqual(b['Functions/motion_ax']['start_time'],'2')
        for path in ['AuxKernels/accel_x','AuxKernels/velocity_x']:
            self.assertEqual(b[path]['static_end_time'],'2')
        self.assertEqual(b['Functions/dynamic_on']['expression'],"'t>2'")
        self.assertNotIn('reset',render(b).lower())

    def test_material_and_contact_inputs_unchanged_and_run_gate_active(self):
        before,b,_=self.compose()
        for path in ['Materials/cdp_update','Constraints/contact_law_1']:
            self.assertEqual(b[path],before[path])
        self.assertEqual(b['UserObjects/aq_expert_pending']['expression'],"'1'")
        self.assertEqual(b['UserObjects/aq_expert_pending']['execute_on'],'INITIAL')
        self.assertEqual(b['Executioner']['nl_abs_tol'],'1e-5')
        self.assertEqual(b['Executioner']['nl_rel_tol'],'1e-8')

    def test_output_intervals_regions_and_no_artificial_solve_sync(self):
        _,b,stats=self.compose()
        expected={'overview':'.1','concrete_damage':'.5','critical_cloud':'.05',
                  'critical_native':'.05','controls':'.01','support_reaction':'.1'}
        for name,dt in expected.items():
            p=b['Outputs/'+name];self.assertEqual(p['mark_interval'],dt)
            self.assertEqual(p['phase_origin'],'2');self.assertNotIn('sync_times',p)
            self.assertNotIn('min_simulation_time_interval',p)
        self.assertEqual(stats['control_channels'],9)
        self.assertEqual(b['VectorPostprocessors/aq_critical_native']['element_ids'],"'99'")
        self.assertEqual(b['Outputs/critical_cloud']['file'],'critical-output.msh')
        self.assertNotIn('Outputs/fields',b)
        self.assertEqual(b['Outputs']['csv'],'false')

    def test_phase_knots_initial_steps_and_minimums(self):
        _,b,_=self.compose()
        clock=b['Functions/dt_clock']
        self.assertEqual(clock['direction'],'LEFT_INCLUSIVE')
        x=list(map(float,clock['x'].strip("'").split()))
        y=list(map(float,clock['y'].strip("'").split()))
        for t in [0.,1.,2.]:self.assertEqual(y[x.index(t)],.001)
        p=b['Executioner/TimeStepper']
        self.assertEqual(p['static_dtmin'],'1e-10');self.assertEqual(p['dynamic_dtmin'],'1e-8')

    def test_default_legacy_integration_does_not_enable_quasistatic_option(self):
        text=(ROOT/'src/timeintegrators/AqueductHHT.C').read_text()
        self.assertIn('"quasistatic_initialization",false',text)
        self.assertIn('_sys.solutionUDot()->zero()',text)
        self.assertIn('_sys.solutionUDotDot()->zero()',text)
        self.assertIn('_du_dotdot_du=0.;computeDuDotDu()',text)
        self.assertIn('quasiStatic() || _t_step<=_inactive_tsteps',text)

    def test_wrong_source_rejected_before_creating_output(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'source').write_text('wrong');(p/'base').write_text('wrong')
            with self.assertRaisesRegex(ValueError,'identity'):
                build(p,p/'base',p/'source',p/'out')
            self.assertFalse((p/'out').exists())

if __name__=='__main__':unittest.main()

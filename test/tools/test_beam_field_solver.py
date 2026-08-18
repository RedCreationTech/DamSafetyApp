import importlib.util
import unittest
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix


PROJECT_DIR = Path(__file__).resolve().parents[2]
TOOL_PATH = PROJECT_DIR / 'tools' / 'abaqus' / 'beam_field_solver.py'
SPEC = importlib.util.spec_from_file_location('beam_field_solver', TOOL_PATH)
SOLVER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SOLVER)


def integrate_sdof(mass, damping, stiffness, dt, nsteps, gamma, beta):
    """单自由度 Newmark 自由振动 (u0=1, v0=0), 返回各步位移。"""
    a0, a1, a2, a3, a4, a5 = SOLVER.newmark_coefficients(gamma, beta, dt)
    effective = a0 * mass + a1 * damping + stiffness
    u, v, acc = 1.0, 0.0, -stiffness / mass
    history = [u]
    for _ in range(nsteps):
        b = (mass * (a0 * u + a2 * v + a3 * acc) +
             damping * (a1 * u + a4 * v + a5 * acc))
        u1 = b / effective
        acc1 = a0 * (u1 - u) - a2 * v - a3 * acc
        v1 = v + dt * ((1 - gamma) * acc + gamma * acc1)
        u, v, acc = u1, v1, acc1
        history.append(u)
    return np.array(history)


class BeamFieldSolverTest(unittest.TestCase):
    def test_circular_section_limits(self):
        self.assertEqual(
            SOLVER.section_fiber_limits('CIRC', [25.0]),
            (25.0, 25.0))

    def test_newmark_defaults_unchanged(self):
        """CLI 默认值保持常加速度法, 既有复现能力不回归。"""
        self.assertEqual(SOLVER.NEWMARK_GAMMA, 0.5)
        self.assertEqual(SOLVER.NEWMARK_BETA, 0.25)
        args = SOLVER.build_parser().parse_args([
            '--report', 'r.json', '--mesh', 'm.e',
            '--displacement-out', 'd.e', '--acceleration-out', 'a.e',
            '--stress-out', 's.e'])
        self.assertEqual(args.newmark_gamma, 0.5)
        self.assertEqual(args.newmark_beta, 0.25)
        # 默认系数与旧内联公式逐项一致
        dt = 0.01
        a0, a1, a2, a3, a4, a5 = SOLVER.newmark_coefficients(0.5, 0.25, dt)
        self.assertAlmostEqual(a0, 1.0 / (0.25 * dt * dt))
        self.assertAlmostEqual(a1, 0.5 / (0.25 * dt))
        self.assertAlmostEqual(a2, 1.0 / (0.25 * dt))
        self.assertAlmostEqual(a3, 1.0 / (2 * 0.25) - 1.0)
        self.assertAlmostEqual(a4, 0.5 / 0.25 - 1.0)
        self.assertAlmostEqual(a5, dt * (0.5 / (2 * 0.25) - 1.0))

    def test_undamped_default_conserves_amplitude(self):
        """γ=0.5/β=0.25 无算法阻尼: 高频自由振动幅值守恒。"""
        dt = 0.01
        omega = 50.0                    # ω·dt = 0.5, 高频
        history = integrate_sdof(1.0, 0.0, omega ** 2, dt, 2000, 0.5, 0.25)
        amplitude = np.abs(history[-400:]).max()
        self.assertAlmostEqual(amplitude, 1.0, delta=0.02)

    def test_damped_newmark_attenuates_high_frequency(self):
        """γ=0.55/β=0.2756 (等价 HHT-α α=-0.05) 衰减高频人为振荡。"""
        dt = 0.01
        omega = 50.0                    # ω·dt = 0.5, 高频
        history = integrate_sdof(
            1.0, 0.0, omega ** 2, dt, 2000, 0.55, 0.2756)
        amplitude = np.abs(history[-400:]).max()
        self.assertLess(amplitude, 0.5)

    def test_axial_stress_recovery(self):
        elastic_modulus = 200.0
        area = 2.0
        length = 1.0
        iy = iz = j = 1.0
        stiffness = SOLVER.beam_ke_local(
            elastic_modulus, 80.0, area, iy, iz, j, length, area, area)
        displacement = np.zeros((1, 2, 6))
        displacement[0, 1, 0] = 0.01
        records = {'beam': [{
            'nodes': (0, 1),
            'transform': np.eye(12),
            'ke': stiffness,
            'keep': list(range(12)),
            'area': area,
            'iy': iy,
            'iz': iz,
            'j': j,
            'cy': 1.0,
            'cz': 1.0,
            'e': elastic_modulus,
        }]}

        fields = SOLVER.recover_beam_stress(displacement, records)

        self.assertAlmostEqual(fields['axial_stress']['beam'][0, 0], 2.0)
        self.assertAlmostEqual(fields['bending_stress']['beam'][0, 0], 0.0)
        self.assertAlmostEqual(fields['torsional_shear']['beam'][0, 0], 0.0)
        self.assertAlmostEqual(fields['vonmises_stress']['beam'][0, 0], 2.0)
        # 单轴应力: 主应力 (2, 0, 0), 轴向应变 σ/E
        self.assertAlmostEqual(fields['s_max_principal']['beam'][0, 0], 2.0)
        self.assertAlmostEqual(fields['s_mid_principal']['beam'][0, 0], 0.0)
        self.assertAlmostEqual(fields['s_min_principal']['beam'][0, 0], 0.0)
        self.assertAlmostEqual(
            fields['axial_strain']['beam'][0, 0], 2.0 / elastic_modulus)
        self.assertAlmostEqual(fields['bending_strain']['beam'][0, 0], 0.0)

    def test_principal_stress_with_torsion(self):
        """σ 与 τ 组合: 主值 σ/2±sqrt(σ²/4+τ²) 与 vonmises 口径自洽。"""
        elastic_modulus = 200.0
        area = 2.0
        length = 1.0
        iy = iz = 1.0
        j = 2.0
        shear_modulus = 80.0
        stiffness = SOLVER.beam_ke_local(
            elastic_modulus, shear_modulus, area, iy, iz, j, length,
            area, area)
        displacement = np.zeros((1, 2, 6))
        displacement[0, 1, 0] = 0.01    # 轴向
        displacement[0, 1, 3] = 0.02    # 扭转
        records = {'beam': [{
            'nodes': (0, 1),
            'transform': np.eye(12),
            'ke': stiffness,
            'keep': list(range(12)),
            'area': area,
            'iy': iy,
            'iz': iz,
            'j': j,
            'cy': 1.0,
            'cz': 1.0,
            'e': elastic_modulus,
        }]}

        fields = SOLVER.recover_beam_stress(displacement, records)

        sigma = fields['axial_stress']['beam'][0, 0]
        tau = shear_modulus * j / length * 0.02 * 1.0 / j
        radius = np.sqrt(sigma ** 2 / 4.0 + tau ** 2)
        self.assertAlmostEqual(fields['torsional_shear']['beam'][0, 0], tau)
        self.assertAlmostEqual(
            fields['s_max_principal']['beam'][0, 0], sigma / 2.0 + radius)
        self.assertAlmostEqual(
            fields['s_min_principal']['beam'][0, 0], sigma / 2.0 - radius)
        self.assertAlmostEqual(
            fields['vonmises_stress']['beam'][0, 0],
            np.sqrt(sigma ** 2 + 3.0 * tau ** 2))

    def test_base_reaction_recovery(self):
        """约束自由度反力 = M·a+C·v+K·u, 非约束自由度为零。"""
        k3 = np.array([
            [4.0, -1.0, 0.0],
            [-1.0, 4.0, -1.0],
            [0.0, -1.0, 4.0]])
        stiffness = csr_matrix(np.block([
            [k3, np.zeros((3, 3))], [np.zeros((3, 3)), 2.0 * k3]]))
        mass = csr_matrix(np.diag([2.0, 3.0, 5.0, 1.0, 1.0, 1.0]))
        nfr = 2
        u6 = np.zeros((nfr, 1, 6))
        v6 = np.zeros((nfr, 1, 6))
        a6 = np.zeros((nfr, 1, 6))
        u6[:, 0, :] = [[1.0, 2.0, 3.0, 0.0, 0.1, 0.2],
                       [0.5, 1.5, 2.5, 0.1, 0.0, 0.3]]
        v6[:, 0, :] = [[0.1, 0.2, 0.3, 0.0, 0.0, 0.0],
                       [0.4, 0.5, 0.6, 0.0, 0.0, 0.0]]
        a6[:, 0, :] = [[10.0, 20.0, 30.0, 1.0, 2.0, 3.0],
                       [15.0, 25.0, 35.0, 4.0, 5.0, 6.0]]
        alpha, beta = 0.5, 0.25

        reactions = SOLVER.recover_base_reactions(
            stiffness, mass, alpha, beta, u6, v6, a6, [0, 5])

        damping = alpha * mass + beta * stiffness
        for frame in range(nfr):
            expected = (mass @ a6[frame, 0, :] +
                        damping @ v6[frame, 0, :] +
                        stiffness @ u6[frame, 0, :])
            # dof 0 为约束平动自由度, 恢复反力
            self.assertAlmostEqual(reactions[frame, 0, 0], expected[0])
        # 约束转动自由度 dof 5 不输出 (rf 仅 x/y/z)
        # 非约束平动自由度 (dof 1, 2) 全部为零
        self.assertTrue(np.all(reactions[:, 0, 1] == 0.0))
        self.assertTrue(np.all(reactions[:, 0, 2] == 0.0))


if __name__ == '__main__':
    unittest.main()

#!/usr/bin/env python3
"""Render the 2D dam added-mass dynamic prototype with ParaView.

Usage:
  pvpython tools/abaqus/render_dam_dynamic.py [result.e] [output_dir]

The video contains only solver time states present in Exodus. Displacement is
visually amplified and the scalar color is element von Mises stress.
"""

import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np


PROJECT_DIR = Path(__file__).resolve().parents[2]
CASE_DIR = PROJECT_DIR / '.build' / 'cases' / 'abaqus-2d-dam-p0'
RESULT = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else (
    CASE_DIR / 'results' / 'dam_dynamic_added_mass_smoke.e')
OUTPUT_DIR = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else (
    CASE_DIR / 'renders' / 'dam_dynamic_added_mass_prototype')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RESOLUTION = [1600, 900]
FPS = 25
SECONDS_PER_STATE = 1


def exodus_stats(path):
    with netCDF4.Dataset(path) as dataset:
        names = [str(value).strip() for value in netCDF4.chartostring(
            dataset.variables['name_nod_var'][:])]
        indices = {name: index + 1 for index, name in enumerate(names)}
        dx = np.asarray(dataset.variables[
            f"vals_nod_var{indices['disp_x']}"][:])
        dy = np.asarray(dataset.variables[
            f"vals_nod_var{indices['disp_y']}"][:])
        max_disp = float(np.sqrt(dx * dx + dy * dy).max())

        elem_names = [str(value).strip() for value in netCDF4.chartostring(
            dataset.variables['name_elem_var'][:])]
        vm_index = elem_names.index('vonmises_stress') + 1
        vm_values = []
        for name, variable in dataset.variables.items():
            if name.startswith(f'vals_elem_var{vm_index}eb'):
                vm_values.append(np.asarray(variable[:]))
        max_vm = max(float(np.nanmax(values)) for values in vm_values)
        times = np.asarray(dataset.variables['time_whole'][:]).tolist()
        x = np.asarray(dataset.variables['coordx'][:])
        y = np.asarray(dataset.variables['coordy'][:])
        bounds = (float(x.min()), float(x.max()), float(y.min()), float(y.max()))
    return max_disp, max_vm, times, bounds


MAX_DISP, MAX_VM, TIMES, BOUNDS = exodus_stats(RESULT)
X_MIN, X_MAX, Y_MIN, Y_MAX = BOUNDS
HEIGHT = Y_MAX - Y_MIN
WARP_SCALE = 0.04 * HEIGHT / max(MAX_DISP, 1e-15)

from paraview.simple import (  # noqa: E402
    ColorBy,
    ExodusIIReader,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetColorTransferFunction,
    GetScalarBar,
    Render,
    SaveScreenshot,
    Show,
    Text,
    WarpByVector,
)


def render():
    print(f'[input] {RESULT}')
    print(f'[stats] states={len(TIMES)}, max|disp|={MAX_DISP:.6e} m, '
          f'max_vonmises={MAX_VM:.6e} Pa, warp={WARP_SCALE:.6e}')

    reader = ExodusIIReader(FileName=[str(RESULT)])
    reader.PointVariables = ['disp_']
    reader.ElementVariables = ['vonmises_stress']
    reader.ElementBlocks = list(reader.ElementBlocks.Available)

    warp = WarpByVector(Input=reader)
    warp.Vectors = ['POINTS', 'disp_']
    warp.ScaleFactor = WARP_SCALE

    view = GetActiveViewOrCreate('RenderView')
    display = Show(warp, view)
    display.Representation = 'Surface With Edges'
    display.EdgeColor = [0.18, 0.18, 0.18]
    ColorBy(display, ('CELLS', 'vonmises_stress'))
    lut = GetColorTransferFunction('vonmises_stress')
    lut.ApplyPreset('Cool to Warm (Extended)', True)
    lut.RescaleTransferFunction(0.0, max(MAX_VM, 1.0))

    view.Background = [0.08, 0.09, 0.12]
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1
    view.CameraFocalPoint = [(X_MIN + X_MAX) / 2,
                             (Y_MIN + Y_MAX) / 2, 0.0]
    view.CameraPosition = [(X_MIN + X_MAX) / 2,
                           (Y_MIN + Y_MAX) / 2, 3.0 * HEIGHT]
    view.CameraViewUp = [0.0, 1.0, 0.0]
    view.ResetCamera()

    scalar_bar = GetScalarBar(lut, view)
    scalar_bar.Title = 'von Mises stress (Pa)'
    scalar_bar.ComponentTitle = ''
    scalar_bar.TitleFontSize = 10
    scalar_bar.LabelFontSize = 9
    scalar_bar.ScalarBarLength = 0.42
    scalar_bar.Visibility = 1

    title = Text()
    title.Text = 'DamSafetyApp | 2D dam | added-mass dynamic prototype'
    title_display = Show(title, view)
    title_display.FontSize = 8
    title_display.Color = [0.92, 0.92, 0.92]
    title_display.WindowLocation = 'Upper Center'

    time_label = Text()
    time_display = Show(time_label, view)
    time_display.FontSize = 7
    time_display.Color = [0.92, 0.92, 0.92]
    time_display.WindowLocation = 'Lower Left Corner'

    scene = GetAnimationScene()
    scene.PlayMode = 'Sequence'
    for index, time_value in enumerate(TIMES):
        scene.AnimationTime = float(time_value)
        time_label.Text = f't = {time_value:.3f} s | state {index + 1}/{len(TIMES)}'
        Render(view)
        frame = OUTPUT_DIR / f'frame-{index:04d}.png'
        SaveScreenshot(str(frame), view, ImageResolution=RESOLUTION)
        print(f'[frame] {frame}')

    video = OUTPUT_DIR / 'dam_dynamic_added_mass_prototype.mp4'
    subprocess.run([
        'ffmpeg', '-y', '-framerate', str(1 / SECONDS_PER_STATE),
        '-i', str(OUTPUT_DIR / 'frame-%04d.png'),
        '-c:v', 'libx264', '-r', str(FPS), '-pix_fmt', 'yuv420p',
        '-movflags', '+faststart', str(video),
    ], check=True)
    print(f'[video] {video} ({video.stat().st_size} bytes)')


if __name__ == '__main__':
    render()

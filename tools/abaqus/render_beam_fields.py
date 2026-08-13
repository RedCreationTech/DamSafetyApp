#!/usr/bin/env python3
"""Render displacement, beam stress and acceleration videos from Exodus.

Usage:
  pvpython tools/abaqus/render_beam_fields.py \
    displacement.e stress.e acceleration.e output_directory
"""

import shutil
import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np


DISPLACEMENT = Path(sys.argv[1]).resolve()
STRESS = Path(sys.argv[2]).resolve()
ACCELERATION = Path(sys.argv[3]).resolve()
OUTPUT_DIR = Path(sys.argv[4]).resolve()
MODE = sys.argv[5] if len(sys.argv) > 5 else 'all'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

FRAME_STRIDE = 33
FPS = 20
RESOLUTION = [1600, 900]
HEIGHT = 17360.0
CAMERA_POSITION = [14000.0, -19000.0, 11000.0]
CAMERA_FOCAL_POINT = [0.0, 0.0, HEIGHT / 2.0]


def names_of(dataset, key):
    return [str(value).strip() for value in netCDF4.chartostring(
        dataset.variables[key][:])]


def vector_max(path, components):
    with netCDF4.Dataset(path) as dataset:
        names = names_of(dataset, 'name_nod_var')
        indices = {name: index + 1 for index, name in enumerate(names)}
        values = [np.asarray(dataset.variables[
            f'vals_nod_var{indices[name]}'][:]) for name in components]
    return float(np.sqrt(sum(value * value for value in values)).max())


def element_max(path, field):
    maximum = 0.0
    with netCDF4.Dataset(path) as dataset:
        names = names_of(dataset, 'name_elem_var')
        field_index = names.index(field) + 1
        prefix = f'vals_elem_var{field_index}eb'
        for name, variable in dataset.variables.items():
            if name.startswith(prefix):
                maximum = max(maximum, float(np.nanmax(variable[:])))
    return maximum


def make_relative_displacement(source, destination):
    """Remove mean base translation for deformation visualization only."""
    shutil.copy(source, destination)
    with netCDF4.Dataset(destination, 'r+') as dataset:
        names = names_of(dataset, 'name_nod_var')
        indices = {name: index + 1 for index, name in enumerate(names)}
        base = np.asarray(dataset.variables['coordz'][:]) < 1.0
        for component in ('disp_x', 'disp_y', 'disp_z'):
            variable = dataset.variables[f'vals_nod_var{indices[component]}']
            values = np.asarray(variable[:])
            variable[:] = values - values[:, base].mean(axis=1)[:, None]


RELATIVE_DISPLACEMENT = OUTPUT_DIR / 'displacement-relative-for-render.e'
make_relative_displacement(DISPLACEMENT, RELATIVE_DISPLACEMENT)
MAX_DISPLACEMENT = vector_max(
    RELATIVE_DISPLACEMENT, ('disp_x', 'disp_y', 'disp_z'))
MAX_ACCELERATION = vector_max(
    ACCELERATION, ('accel_x', 'accel_y', 'accel_z'))
MAX_STRESS = element_max(STRESS, 'vonmises_stress')
WARP_SCALE = 0.08 * HEIGHT / max(MAX_DISPLACEMENT, 1e-12)

from paraview.simple import (  # noqa: E402
    ColorBy,
    Delete,
    ExodusIIReader,
    GetActiveViewOrCreate,
    GetAnimationScene,
    GetColorTransferFunction,
    GetSources,
    GetScalarBar,
    Render,
    ResetSession,
    SaveScreenshot,
    Show,
    Text,
    WarpByVector,
)


def encode(frame_directory, destination):
    subprocess.run([
        'ffmpeg', '-y', '-framerate', str(FPS),
        '-i', str(frame_directory / 'frame-%04d.png'),
        '-c:v', 'libx264', '-pix_fmt', 'yuv420p',
        '-movflags', '+faststart', str(destination),
    ], check=True)


def configure_view(display, lut, title, scalar_title):
    view = GetActiveViewOrCreate('RenderView')
    view.Background = [0.08, 0.09, 0.12]
    view.CameraPosition = CAMERA_POSITION
    view.CameraFocalPoint = CAMERA_FOCAL_POINT
    view.CameraViewUp = [0.0, 0.0, 1.0]
    view.OrientationAxesVisibility = 0

    scalar_bar = GetScalarBar(lut, view)
    scalar_bar.Title = scalar_title
    scalar_bar.ComponentTitle = ''
    scalar_bar.TitleFontSize = 10
    scalar_bar.LabelFontSize = 9
    scalar_bar.ScalarBarLength = 0.42
    scalar_bar.Visibility = 1

    heading = Text()
    heading.Text = title
    heading_display = Show(heading, view)
    heading_display.FontSize = 8
    heading_display.Color = [0.92, 0.92, 0.92]
    heading_display.WindowLocation = 'Upper Center'

    time_label = Text()
    time_display = Show(time_label, view)
    time_display.FontSize = 7
    time_display.Color = [0.92, 0.92, 0.92]
    time_display.WindowLocation = 'Lower Left Corner'
    return view, time_label


def render_scene(name, source_path, association, field, maximum,
                 title, scalar_title, warp=False):
    reader = ExodusIIReader(FileName=[str(source_path)])
    reader.ElementBlocks = list(reader.ElementBlocks.Available)
    if association == 'POINTS':
        reader.PointVariables = [field]
    else:
        reader.ElementVariables = [field]

    source = reader
    if warp:
        source = WarpByVector(Input=reader)
        source.Vectors = ['POINTS', field]
        source.ScaleFactor = WARP_SCALE

    view = GetActiveViewOrCreate('RenderView')
    display = Show(source, view)
    display.Representation = 'Wireframe'
    display.LineWidth = 4.0
    ColorBy(display, (association, field))
    lut = GetColorTransferFunction(field)
    lut.ApplyPreset('Cool to Warm (Extended)', True)
    lut.RescaleTransferFunction(0.0, max(maximum, 1e-12))
    view, time_label = configure_view(
        display, lut, title, scalar_title)

    times = list(reader.TimestepValues)[::FRAME_STRIDE]
    frame_directory = OUTPUT_DIR / f'{name}-frames'
    frame_directory.mkdir(exist_ok=True)
    scene = GetAnimationScene()
    scene.PlayMode = 'Sequence'
    print(f'[render] {name}: {len(times)} frames')
    for index, time_value in enumerate(times):
        scene.AnimationTime = float(time_value)
        time_label.Text = f't = {time_value:.2f} s'
        Render(view)
        frame = frame_directory / f'frame-{index:04d}.png'
        SaveScreenshot(str(frame), view, ImageResolution=RESOLUTION)
        if index % max(1, len(times) // 5) == 0:
            print(f'  {index + 1}/{len(times)}')
    video = OUTPUT_DIR / f'PR-RG-400gal-X_{name}.mp4'
    encode(frame_directory, video)
    print(f'[video] {video} ({video.stat().st_size} bytes)')
    for source_proxy in list(GetSources().values()):
        Delete(source_proxy)
    ResetSession()


if __name__ == '__main__':
    print(f'[range] displacement={MAX_DISPLACEMENT:.6g} mm '
          f'acceleration={MAX_ACCELERATION:.6g} mm/s^2 '
          f'stress={MAX_STRESS:.6g} MPa warp={WARP_SCALE:.6g}')
    if MODE in ('all', 'displacement'):
        render_scene(
            'displacement', RELATIVE_DISPLACEMENT, 'POINTS', 'disp_',
            MAX_DISPLACEMENT,
            'PR-RG-400gal-X | relative displacement field',
            '|displacement| (mm)', warp=True)
    if MODE in ('all', 'stress'):
        render_scene(
            'stress', STRESS, 'CELLS', 'vonmises_stress', MAX_STRESS,
            'PR-RG-400gal-X | beam equivalent stress envelope',
            'von Mises envelope (MPa)')
    if MODE in ('all', 'acceleration'):
        render_scene(
            'acceleration', ACCELERATION, 'POINTS', 'accel_',
            MAX_ACCELERATION,
            'PR-RG-400gal-X | absolute acceleration field',
            '|acceleration| (mm/s^2)')

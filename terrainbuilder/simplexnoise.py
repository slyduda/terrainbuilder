#############################################################################
#  Copyright (C) 2020 - 2021 Sylvester Duda
#
#  This is code is derivative of the Simple Pseudo-random Number Generators
#  Available on GitHub https://github.com/cmcqueen/simplerandom
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, see http://www.gnu.org/licenses/
############################################################################

from numpy import pi, sin, cos
from opensimplex import OpenSimplex

simplex = OpenSimplex()


def octave_noise_2d(octaves: float, persistence: float, scale: float, x: float, y: float):
    total = 0
    frequency = scale
    amplitude = 1.0

    # We have to keep track of the largest possible amplitude,
    # because each octave adds more, and we need a value in [-1, 1].
    max_amplitude = 0.0

    for i in range(octaves):
        total += simplex.noise2d(x * frequency, y * frequency) * amplitude

        frequency *= 2
        max_amplitude += amplitude
        amplitude *= persistence

    return total / max_amplitude


def octave_noise_3d(octaves: float, persistence: float, scale: float, x: float, y: float, z: float):
    total = 0
    frequency = scale
    amplitude = 1.0

    # We have to keep track of the largest possible amplitude,
    # because each octave adds more, and we need a value in [-1, 1].
    max_amplitude = 0.0

    for i in range(octaves):
        total += simplex.noise3d(x * frequency,
                                 y * frequency, z * frequency) * amplitude

        frequency *= 2
        max_amplitude += amplitude
        amplitude *= persistence

    return total / max_amplitude


def octave_noise_4d(octaves: float, persistence: float, scale: float, x: float, y: float, z: float, w: float):
    total = 0
    frequency = scale
    amplitude = 1.0

    # We have to keep track of the largest possible amplitude,
    # because each octave adds more, and we need a value in [-1, 1].
    max_amplitude = 0.0

    for i in range(octaves):
        total += simplex.noise4d(x * frequency,
                                 y * frequency, z * frequency, w * frequency) * amplitude

        frequency *= 2
        max_amplitude += amplitude
        amplitude *= persistence

    return total / max_amplitude


def scaled_octave_noise_2d(octaves: float,  persistence: float, scale: float, loBound: float, hiBound: float, x: float, y: float):
    """2D Scaled Multi-octave Simplex noise.

    Returned value will be between loBound and hiBound.

    Args:
        octaves (float): [description]
        persistence (float): [description]
        scale (float): [description]
        loBound (float): [description]
        hiBound (float): [description]
        x (float): [description]
        y (float): [description]

    Returns:
        [type]: [description]
    """
    return octave_noise_2d(octaves, persistence, scale, x, y) * (hiBound - loBound) / 2 + (hiBound + loBound) / 2


def scaled_octave_noise_3d(octaves: float,  persistence: float, scale: float, loBound: float, hiBound: float, x: float, y: float, z: float):
    """3D Scaled Multi-octave Simplex noise.

    Returned value will be between loBound and hiBound.

    Args:
        octaves (float): [description]
        persistence (float): [description]
        scale (float): [description]
        loBound (float): [description]
        hiBound (float): [description]
        x (float): [description]
        y (float): [description]
        z (float): [description]

    Returns:
        [type]: [description]
    """
    return octave_noise_3d(octaves, persistence, scale, x, y, z) * (hiBound - loBound) / 2 + (hiBound + loBound) / 2


def scaled_octave_noise_4d(octaves: float,  persistence: float, scale: float, loBound: float, hiBound: float, x: float, y: float, z: float, w: float):
    """4D Scaled Multi-octave Simplex noise.

    Returned value will be between loBound and hiBound.

    Args:
        octaves (float): [description]
        persistence (float): [description]
        scale (float): [description]
        loBound (float): [description]
        hiBound (float): [description]
        x (float): [description]
        y (float): [description]
        z (float): [description]
        w (float): [description]

    Returns:
        [type]: [description]
    """
    return octave_noise_4d(octaves, persistence, scale, x, y, z, w) * (hiBound - loBound) / 2 + (hiBound + loBound) / 2


def fast_floor(x: float):
    return int(x) if x > 0 else int(x) - 1


def dot(g: list, x: float, y: float, z: float = 0.0, w: float = 0.0):
    """[summary]

    Args:
        g (list): [description]
        x (float): [description]
        y (float): [description]
        z (float, optional): [description]. Defaults to 0.0.
        w (float, optional): [description]. Defaults to 0.0.

    Returns:
        [type]: [description]
    """
    if w:
        return g[0]*x + g[1]*y + g[2]*z + g[3]*w
    if z:
        return g[0]*x + g[1]*y + g[2]*z
    return g[0]*x + g[1]*y


def simplexnoise(seed, map, width: int, height: int, roughness: float):
    ka = 256/seed
    kb = seed*567 % 256
    kc = (seed*seed) % 256
    kd = (567-seed) % 256
    for y in range(height):
        for x in range(width):
            fNX = x/float(width)  # we let the x-offset define the circle
            fNY = y/float(height)  # we let the x-offset define the circle
            fRdx = fNX*2*pi  # a full circle is two pi radians
            fRdy = fNY*4*pi  # a full circle is two pi radians
            # fYSin = sin(fRdy) Removed because redundant
            fRdsSin = 1.0
            noise_scale = 0.593
            a = fRdsSin*sin(fRdx)
            b = fRdsSin*cos(fRdx)
            c = fRdsSin*sin(fRdy)
            d = fRdsSin*cos(fRdy)
            v = scaled_octave_noise_4d(64.0, roughness, 2.0, 0.0, 1.0,
                                       ka+a*noise_scale,
                                       kb+b*noise_scale,
                                       kc+c*noise_scale,
                                       kd+d*noise_scale)
            if (map[y * width + x] == 0.0):
                map[y * width + x] = v
    return map

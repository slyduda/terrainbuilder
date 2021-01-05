#############################################################################
#  Copyright (C) 2020 - 2021 Sylvester Duda
#
#  This is code is derivative of the C++ Plate Tectonics implementation
#  Available on GitHub https://github.com/Mindwerks/plate-tectonics
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

from numpy.random import RandomState, SeedSequence
from numpy import array, int32


def opposite_of(a: int or float):
    a = bool(a)
    a = not a
    a = int(a)
    return a


def normalize(arr, size: int):
    min = arr[0]
    max = arr[0]

    for i in range(size):
        min = min if min < arr[i] else arr[i]
        max = max if max > arr[i] else arr[i]

    diff = max - min

    if diff > 0:
        for i in range(size):
            arr[i] = (arr[i] - min) / diff

    return arr


def sqrdmd(seed: SeedSequence, map: list, size: int, rgh: float):
    _randstate = RandomState(seed=seed)

    full_size = size * size

    i = 0
    x, y, dx, dy = (0, 0, 0, 0)
    x0, x1, y0, y1 = (0, 0, 0, 0)
    p0, p1, p2, p3 = (0, 0, 0, 0)
    step, line_jump, masked = (0, 0, 0)
    slope, sum, center_sum = (0.0, 0.0, 0.0)

    temp = size - 1
    if temp & (temp - 1) or temp & 3:
        assert('ERROR Side should be 2**n +1')

    temp = size
    slope = rgh
    step = size & ~1

    def CALC_SUM(a, b, c, d, rnd):
        sum = ((a) + (b) + (c) + (d)) * 0.25
        return sum + slope * rnd

    def SAVE_SUM(a):
        is_zero = int(map[a]) == 0
        if is_zero:
            map[a] = sum

    dy = step * size
    CALC_SUM(map[0], map[step], map[dy], map[dy + step],
             _randstate.random_integers(1, int32))
    SAVE_SUM(i)
    center_sum = sum

    p0 = step >> 1
    CALC_SUM(map[0], map[step], center_sum, center_sum,
             _randstate.random_integers(1, int32))
    SAVE_SUM(p0)

    p1 = p0 * size
    CALC_SUM(map[0], map[dy], center_sum, center_sum,
             _randstate.random_integers(1, int32))
    SAVE_SUM(p1)
    map[full_size + p0 - size] = map[p0]  # Copy top val into btm row.
    map[p1 + size - 1] = map[p1]  # Copy left value into right column.
    slope *= rgh
    step >>= 1

    while step > 1:  # Enter the main loop.
        # *************************************************************
        # Calc midpoint of sub squares on the map ("diamond step"). *
        # ************************************************************/
        dx = step
        dy = step * size
        i = (step >> 1) * (size + 1)
        line_jump = step * size + 1 + step - size
        y0 = 0
        y1 = dy
        while y1 < size * size:
            x0 = 0
            x1 = dx
            while x1 < size:
                sum = (map[y0+x0] + map[y0+x1] +
                       map[y1+x0] + map[y1+x1]) * 0.25
                sum = sum + slope * _randstate.random_integers(1, int32)
                masked = int(not bool(int(map[i])))
                map[i] = map[i] * opposite_of(masked) + sum * masked
                # There's additional step taken at the end of last
                # valid loop. That step actually isn't valid because
                # the row ends right then. Thus we are forced to
                # manually remove it after the loop so that 'i'
                # points again to the index accessed last.
                # /
                x0 += dx
                x1 += dx
                i += step

            y0 += dy
            y1 += dy
            i += line_jump - step

        # **************************************************************
            # Calculate each sub diamonds' center point ("square step").
            # Diamond gets its left and right vertices from the square
            # corners of last iteration and its top and bottom vertices
            # from the "diamond step" we just performed.
            # ************************************************************/
        i = step >> 1
        p0 = step  # right */
        p1 = i * size + i  # bottom */
        p2 = 0  # left */
        p3 = full_size + i - (i + 1) * size  # top (wrapping edges) */

        # Calculate "diamond" values for top row in map. */
        while p0 < size:
            sum = (map[p0] + map[p1] + map[p2] + map[p3]) * 0.25
            sum = sum + slope * _randstate.random_sample()
            masked = opposite_of((int(map[i])))
            map[i] = map[i] * opposite_of(masked) + sum * masked
            # Copy it into bottom row. */
            map[full_size + i - size] = map[i]
            p0 += step
            p1 += step
            p2 += step
            p3 += step
            i += step

            # Now that top row's values are calculated starting from
            # 'y = step >> 1' both saves us from recalculating same things
            # twice and guarantees that data will not be read beyond top
            # row of map. 'size - (step >> 1)' guarantees that data will
            # not be read beyond bottom row of map.

        y = step >> 1
        temp = 0
        while y < size - (step >> 1):
            p0 = step >> 1  # right
            p1 = p0 * size  # bottom
            p2 = -p0  # left
            p3 = -p1  # top
            # For even rows add step/2. Otherwise add nothing. */
            x = i = p0 * temp  # Init 'x' while it's easy. */
            i += y * size  # Move 'i' into correct row. */
            p0 += i
            p1 += i
            # For odd rows p2 (left) wraps around map edges. */
            p2 += i + (size - 1) * opposite_of(temp)
            p3 += i
            # size - (step >> 1) guarantees that data will not be
            # read beyond rightmost column of map. */
            while x < size - (step >> 1):
                sum = (map[p0] + map[p1] + map[p2] + map[p3]) * 0.25
                sum = sum + slope * _randstate.random_sample()
                masked = opposite_of(int(map[i]))
                map[i] = map[i] * opposite_of(masked) + sum * masked
                p0 += step
                p1 += step
                p2 += step
                p3 += step
                i += step
                # if we start from leftmost column -> left
                # point (p2) is going over the right border ->
                # wrap it around into the beginning of
                # previous rows left line. */
                p2 -= (size - 1) * opposite_of(x)

                # End of while loop
                x += step

                # copy rows first element into its last */
                i = y * size
                map[i + size - 1] = map[i]

            # End of while loop
            y += step >> 1
            temp = opposite_of(temp)

            slope *= rgh  # reduce amount of randomness for next round */
            step = step >> 1  # split squares and diamonds in half */
        return 0

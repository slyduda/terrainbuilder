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

def sqrdmd(seed, map, size: int, rgh: float):
    pass

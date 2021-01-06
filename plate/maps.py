#############################################################################
#  Copyright (C) 2020 - 2021 Sylvester Duda
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
#############################################################################


from numpy import ones

# TODO. I can change these to zeros too to make plates set crust function easier


class HeightMap(object):
    '''
        A class for constructing a matrix for height maps.
    '''

    def __init__(self, x: int, y: int):
        if x == 0 or y == 0:
            assert(x == 0 or y == 0), 'width and height of {} should be greater than 0'.format(
                self.__class__.__name__)
        self.width = x
        self.height = y
        self._data = ones((x, y), dtype=int)
        self._area = x * y

    def set_all(self, v: int):
        for x in self._data:
            for y in x:
                y = v

    def area(self):
        return self.width * self.height

    def x_mod(self, x: int):
        return x % self.width

    def y_mod(self, y: int):
        return y % self.height


class MassMap(object):
    '''
        A class for constructing a matrix for mass maps.
    '''

    def __init__(self, x: int, y: int):
        if x == 0 or y == 0:
            assert(x == 0 or y == 0), 'width and height of {} should be greater than 0'.format(
                self.__class__.__name__)
        self.width = x
        self.height = y
        self._data = ones((x, y), dtype=float)
        self._area = x * y

    def set_all(self, v: float):
        for x in self._data:
            for y in x:
                y = v

    def area(self):
        return self.width * self.height

    def x_mod(self, x: int):
        return x % self.width

    def y_mod(self, y: int):
        return y % self.height


class AgeMap(object):
    '''
        A class for constructing a matrix for age maps.
    '''

    def __init__(self, x: int, y: int):
        if x == 0 or y == 0:
            assert(x == 0 or y == 0), 'width and height of {} should be greater than 0'.format(
                self.__class__.__name__)
        self.width = x
        self.height = y
        self._data = ones((x, y), dtype=int)
        self._area = x * y

    def set_all(self, v: int):
        for x in self._data:
            for y in x:
                y = v

    def area(self):
        return self.width * self.height

    def x_mod(self, x: int):
        return x % self.width

    def y_mod(self, y: int):
        return y % self.height

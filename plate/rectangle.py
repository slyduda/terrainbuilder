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
#############################################################################
import math


class WorldDimension(object):
    def __init__(self, width, height):
        self.width = width
        self.height = height

    def get_width(self):
        return self.width

    def get_height(self):
        return self.height

    def get_max(self):
        if self.width > self.height:
            return self.width
        else:
            return self.height

    def get_area(self):
        return self.width * self.height

    def x_mod(self, x: int):
        return x % self.get_width()

    def y_mod(self, y: int):
        return y % self.get_height()

    def contains(self, x: int, y: int):
        return (x >= 0 and x < self.width and y >= 0 and y < self.height)

    def normalize(self, x: int, y: int):
        print('UNSURE OF THE MATH ON THIS ONE BECAUSE IT TURNS INTO A FLOAT VALUE')
        self.height = int(self.height / y)
        self.width = int(self.width / x)

    def index_of(self, x: int, y: int):
        return y * self.get_width() + x

    def line_index(self, y: int):
        if y < 0 or y >= self.height:
            assert(y < 0 or y >=
                   self.height), "WorldDimension::line_index: y is not valid"
        return self.index_of(0, y)

    def y_from_index(self, index: int):
        return math.floor(index / self.width)

    def x_from_index(self, index: int):
        return index % self.width

    def normalized_index_of(self, x: int, y: int):
        return self.index_of(self.x_mod(x), self.y_mod(y))

    def x_cap(self, x: int):
        if x < self.width:
            return x
        else:
            return self.width - 1

    def y_cap(self, y: int):
        if y < self.height:
            return y
        else:
            return self.height - 1

    # def larger_size():
        # Use the getMax function


class Rectangle(object):
    def __init__(self, wd: WorldDimension, left: int, right: int, top: int, bottom: int):
        self.wd = wd
        self.left = left
        self.right = right
        self.top = top
        self.bottom = bottom

    def get_map_index(self, px: int, py: int):
        x = px % self.wd.get_width()
        y = py % self.wd.get_height()

        ilft = int(self.left)
        itop = int(self.top)
        irgt = int(self.right) + (self.wd.get_width()
                                  if self.right < ilft else 0)
        ibtm = int(self.bottom) + (self.wd.get_height()
                                   if self.bottom < itop else 0)
        width = irgt - ilft

        if width < 0:
            raise Exception("Failed because width is less than zero")

        xOkA = (x >= ilft) and (x < irgt)
        xOkB = (x + self.wd.get_width() >=
                ilft) and (x + self.wd.get_width() < irgt)
        xOk = xOkA or xOkB

        yOkA = (y >= itop) and (y < ibtm)
        yOkB = (y + self.wd.get_height() >=
                itop) and (y + self.wd.get_height() < ibtm)
        yOk = yOkA or yOkB

        x += self.wd.get_width() if (x < ilft) else 0
        y += self.wd.get_height() if (y < itop) else 0

        x -= ilft  # Calculate offset within local map.
        y -= itop

        if x < 0:
            raise Exception("Failed because x is less than 0")

        if y < 0:
            raise Exception("Failed because y is less than 0")

        if xOk and yOk:
            # px = x These were pointers
            # py = y
            return (y * width + x)
        else:
            return -1

    def get_map_indeces(self, px: int, py: int):
        x = px % self.wd.get_width()
        y = py % self.wd.get_height()

        ilft = int(self.left)
        itop = int(self.top)
        irgt = int(self.right) + (self.wd.get_width()
                                  if self.right < ilft else 0)
        ibtm = int(self.bottom) + (self.wd.get_height()
                                   if self.bottom < itop else 0)
        width = irgt - ilft

        if width < 0:
            raise Exception("Failed because width is less than zero")

        xOkA = (x >= ilft) and (x < irgt)
        xOkB = (x + self.wd.get_width() >=
                ilft) and (x + self.wd.get_width() < irgt)
        xOk = xOkA or xOkB

        yOkA = (y >= itop) and (y < ibtm)
        yOkB = (y + self.wd.get_height() >=
                itop) and (y + self.wd.get_height() < ibtm)
        yOk = yOkA or yOkB

        x += self.wd.get_width() if (x < ilft) else 0
        y += self.wd.get_height() if (y < itop) else 0

        x -= ilft  # Calculate offset within local map.
        y -= itop

        if x < 0:
            raise Exception("Failed because x is less than 0")

        if y < 0:
            raise Exception("Failed because y is less than 0")

        if xOk and yOk:
            return (x, y)
        else:
            return (-1, -1)

    def enlarge_to_contain(self, x: int, y: int):
        if y < self.top:
            self.top = y

        if y > self.bottom:
            self.bottom = y

        if x < self.left:
            self.left = x

        if x > self.right:
            self.right = x

    def get_left(self):
        return self.left

    def get_right(self):
        return self.right

    def get_top(self):
        return self.top

    def get_bottom(self):
        return self.bottom

    def set_left(self, v: int):
        self.left = v

    def set_right(self, v: int):
        self.right = v

    def set_top(self, v: int):
        self.top = v

    def set_bottom(self, v: int):
        self.bottom = v

    def shift(self, dx: int, dy: int):
        self.left += dx
        self.right += dx
        self.top += dy
        self.bottom += dy

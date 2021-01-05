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

from lithosphere import Lithosphere


class platec_api_list_elem(object):
    def __init__(self, _id: int, _data):
        self.data = _data
        self.id = _id


lithospheres = []
last_id = 1


def platec_api_create(seed: int, width: int, height: int, sea_level: float,
                      erosion_period: int, folding_ratio: float,
                      aggr_overlap_abs: int, aggr_overlap_rel: float,
                      cycle_count: int, num_plates: int):

    litho = Lithosphere(seed, width, height, sea_level,
                        erosion_period, folding_ratio, aggr_overlap_abs,
                        aggr_overlap_rel, cycle_count)
    litho.create_plates(num_plates)

    elem = platec_api_list_elem(last_id+1, litho)
    lithospheres.append(elem)

    return litho


def platec_api_destroy():
    lithospheres = []


def platec_api_step(litho: Lithosphere):
    litho.update()

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
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, see http://www.gnu.org/licenses/
#############################################################################


'''
A lithosphere is the rigid, outermost shell of a terrestrial-type planet or natural satellite.

On Earth, it is composed of the crust and the portion of the upper mantle that behaves elastically
on time scales of thousands of years or greater. The crust and upper mantle are distinguished on
the basis of chemistry and mineralogy.
'''

import sys
from numpy import array, zeros, int32
from numpy.random import RandomState, SeedSequence, MT19937

from .rectangle import WorldDimension
from .noise import generate_noise
from .plate import Plate
from .maps import HeightMap, AgeMap


SUBDUCT_RATIO = 0.5

BUOYANCY_BONUS_X = 3
MAX_BUOYANCY_AGE = 20
MULINV_MAX_BUOYANCY_AGE = 1.0 / MAX_BUOYANCY_AGE

RESTART_ENERGY_RATIO = 0.15
RESTART_SPEED_LIMIT = 2.0
RESTART_ITERATIONS = 600
NO_COLLISION_TIME_LIMIT = 10

# NEEDS WORK
CONTINENTAL_BASE = 2.0
OCEANIC_BASE = .1
FLT_EPSILON = sys.float_info.epsilon
BOOL_REGENERATE_CRUST = 1


def findBound(*map, length, x0, y0, dx, dy):
    pass


def findPlate(*plates, x, y, num_plates):
    pass


class PlateCollision(object):
    def __init__(self, _index: int, x: int, y: int, z: float):
        """Container for collision details between two plates.

        In simulation there's usually 2-5% collisions of the entire map
        area. In a 512*512 map that means 5000-13000 collisions.

        When plate collisions are recorded and processed pair-by-pair, some
        of the information is lost if more than two plates collide at the
        same point (there will be no record of the two lower plates
        colliding together, just that they both collided with the tallest
        plate) ONLY IF ALL the collisions between ANY TWO plates of that
        group always include a third, taller/higher  plate. This happens
        most often when plates have long, sharp spikes i.e. in the
        beginning.


        Args:
            _index (int): [description]
            x (int): [description]
            y (int): [description]
            z (float): [description]
        """
        self._index = _index
        self.wx = x
        self.wy = y
        self.crust = z


class PlateArea(object):
    def __init__(self):
        self.btm = None
        self.lft = None
        self.rgt = None
        self.top = None
        self.wdt = None
        self.hgt = None
        self.border = []


class Lithosphere(object):

    def __init__(self, seed: SeedSequence, width: int, height: int, sea_level: float, _erosion_period: int, _folding_ratio: float, aggr_ratio_abs: int, aggr_ration_rel: float, num_cycles: int):
        """Initialize system's height map i.e. topography.

        Args:
            seed (int):
            width (int): Square height map's width in pixels.
            height (int): Square height map's height in pixels.
            sea_level (float): Amount of surface area that becomes oceanic crust.
            _erorsion_period (int): # of iterations between global erosion.
            _folding_ratio (float): Percent of overlapping crust that's folded.
            aggr_ratio_abs (int): # of overlapping points causing aggregation.
            aggr_ration_rel (float): % of overlapping area causing aggregation.
            num_cycles (int):  Number of times system will be restarted.
        """

        self.wd = WorldDimension(width, height)
        self._randstate = RandomState(MT19937(seed))
        self._steps = 0

        self.erosion_period = _erosion_period
        self.folding_ratio = _folding_ratio
        self.aggr_overlap_abs = aggr_ratio_abs
        self.aggr_overlap_rel = aggr_ration_rel

        self.cycle_count = 0
        self.max_cycles = num_cycles
        self.max_plates = 0
        self.num_plates = 0
        self.iter_count = 0

        self.hmap = HeightMap(width, height)
        self.amap = AgeMap(width, height)
        self.plates = []
        self.plate_areas = []

        tmp_dim = WorldDimension(width, height)
        A = tmp_dim.get_area()
        tmp = [0.0 for _ in range(A)]

        tmp = self.create_noise(tmp, tmp_dim, seed, True)

        lowest = tmp[0]
        highest = tmp[0]

        for i in range(1, A):
            lowest = lowest if lowest < tmp[i] else tmp[i]
            highest = highest if highest > tmp[i] else tmp[i]

        for i in range(A):
            tmp[i] = (tmp[i] - lowest) / (highest - lowest)

        sea_threshold = 0.5
        th_step = 0.5

        while (th_step > 0.01):
            count = 0
            for i in range(A):
                count += (tmp[i] < sea_threshold)

            th_step *= 0.5
            if (count / float(A) < sea_level):
                sea_threshold += th_step
            else:
                sea_threshold -= th_step

        sea_level = sea_threshold
        for i in range(A):
            cont = int(tmp[i] > sea_level) * (tmp[i] + CONTINENTAL_BASE)
            ocea = int(tmp[i] <= sea_level) * OCEANIC_BASE
            tmp[i] = cont + ocea

        for i in range(len(tmp)):
            y, x = (i // width, i % width)
            self.hmap._data[y, x] = tmp[i]

        self.peak_ek = 0.0
        self.last_coll_count = 0

        self.collisions = None
        self.subductions = None

        self.imap = [None for _ in range(self.wd.get_area())]

    def create_noise(self, tmp: list, tmp_dim: WorldDimension, seed: SeedSequence, use_simplex: bool = False):
        map = generate_noise(tmp, tmp_dim, seed, use_simplex)
        return map

    def grow_plates(self):
        #  "Grow" plates from their origins until surface is fully populated.
        owner = self.imap  # Create an alias.
        max_border = 1
        i = 0
        while max_border:
            max_border = 0
            for i in range(self.num_plates):
                area = self.plate_areas[i]
                N = len(area.border)
                max_border = max_border if max_border > N else N

                if N == 0:
                    continue

                j = self._randstate.randint(1, high=2147483648) % N
                p = area.border[j]
                cy = self.wd.y_from_index(p)
                cx = self.wd.x_from_index(p)

                width = self.wd.get_width()
                height = self.wd.get_height()
                lft = cx - 1 if cx > 0 else width - 1
                rgt = cx + 1 if cx < width - 1 else 0
                top = cy - 1 if cy > 0 else height - 1
                btm = cy + 1 if cy < height - 1 else 0

                n = top * width + cx  # North.
                s = btm * width + cx  # South.
                w = cy * width + lft  # West.
                e = cy * width + rgt  # East.

                if owner[n] >= self.num_plates:
                    owner[n] = i
                    area.border.append(n)

                    if area.top == self.wd.y_mod(top + 1):
                        area.top = top
                        area.hgt += 1

                if owner[s] >= self.num_plates:
                    owner[s] = i
                    area.border.append(s)

                    if btm == self.wd.y_mod(area.btm + 1):
                        area.btm = btm
                        area.hgt += 1

                if owner[w] >= self.num_plates:
                    owner[w] = i
                    area.border.append(w)

                    if area.lft == self.wd.x_mod(lft + 1):
                        area.lft = lft
                        area.wdt += 1

                if owner[e] >= self.num_plates:
                    owner[e] = i
                    area.border.append(e)

                    if rgt == self.wd.x_mod(area.rgt + 1):
                        area.rgt = rgt
                        area.wdt += 1

                # Overwrite processed point with unprocessed one.
                area.border[j] = area.border[-1]
                area.border.pop()

    def create_plates(self, num_plates):
        map_area = self.wd.get_area()
        self.max_plates = num_plates
        self.num_plates = num_plates

        self.imap = [255 for i in range(map_area)]
        self.collisions = [[] for _ in range(num_plates)]
        self.subductions = [[] for _ in range(num_plates)]

        used = []
        self.plate_areas = [PlateArea() for _ in range(num_plates)]
        for i in range(num_plates):
            area = self.plate_areas[i]
            # Randomly select an unused plate origin.
            p = self._randstate.randint(1, high=2147483648) % (map_area)
            while p in used:
                p = self._randstate.randint(1, high=2147483648) % (map_area)

            y = self.wd.y_from_index(p)
            x = self.wd.x_from_index(p)

            area.lft = area.rgt = x  # Save origin...
            area.top = area.btm = y
            area.wdt = area.hgt = 1

            area.border.append(p)  # ...and mark it as border.
            used.append(p)
            self.imap[p] = i

        self.grow_plates()

        self.plates = [None for _ in range(num_plates)]

        # Extract and create plates from initial terrain.
        for i in range(num_plates):
            area = self.plate_areas[i]
            area.wdt = self.wd.x_cap(area.wdt)
            area.hgt = self.wd.y_cap(area.hgt)

            x0 = area.lft
            x1 = 1 + x0 + area.wdt
            y0 = area.top
            y1 = 1 + y0 + area.hgt
            width = x1 - x0
            height = y1 - y0
            plt = [0.0 for _ in range(width * height)]

            # Copy plate's height data from global map into local map.
            j = 0
            for y in range(y0, y1):
                for x in range(x0, x1):
                    k = self.wd.normalized_index_of(x, y)
                    xf = self.wd.x_from_index(k)
                    yf = self.wd.y_from_index(k)
                    plt[j] = self.hmap._data[yf, xf] * (self.imap[k] == i)
                    j += 1

            self.plates[i] = Plate(self._randstate.randint(1, high=2147483648), plt, width,
                                   height, x0, y0, i, self.wd)

        self.iter_count = num_plates + MAX_BUOYANCY_AGE
        self.peak_ek = 0
        self.last_coll_count = 0

    def get_cycle_count(self):
        return self.cycle_count

    def get_iteration_count(self):
        return self.iter_count

    def get_world_dimension(self):
        return self.wd

    def get_plate_count(self):
        return self.num_plates

    def get_age_map(self):
        return self.amap

    def get_topography(self):
        return self.hmap

    def get_plates_map(self):
        return self.imap

    def update(self):
        self._steps += 1
        total_velocity = 0
        system_kinetic_energy = 0

        for i in range(self.num_plates):
            total_velocity += self.plates[i].get_velocity()
            system_kinetic_energy += self.plates[i].get_momentum()

        if system_kinetic_energy > self.peak_ek:
            self.peak_ek = system_kinetic_energy

        # If there's no continental collisions during past iterations,
        # then interesting activity has ceased and we should restart.
        # Also if the simulation has been going on for too long already,
        # restart, because interesting stuff has most likely ended.
        if total_velocity < RESTART_SPEED_LIMIT or system_kinetic_energy / self.peak_ek < RESTART_ENERGY_RATIO or self.last_coll_count > NO_COLLISION_TIME_LIMIT or self.iter_count > 600:
            print('RESTARTING')
            self.restart()
            return

        map_area = self.wd.get_area()
        prev_imap = list(self.imap)
        self.imap = [255 for _ in range(map_area)]

        # Realize accumulated external forces to each plate.
        for i in range(self.num_plates):
            self.plates[i].reset_segments()

            if self.erosion_period > 0 and self.iter_count % self.erosion_period == 0:
                self.plates[i].erode(CONTINENTAL_BASE)

            self.plates[i].move()

        oceanic_collisions = 0
        continental_collisions = 0

        # Update height and plate index maps.
        # Doing it plate by plate is much faster than doing it index wise:
        # Each plate's map's memory area is accessed sequentially and only
        # once as opposed to calculating "num_plates" indices within plate
        # maps in order to find out which plate(s) own current location.
        self.hmap.set_all(0)

        for i in range(self.num_plates):
            x0 = int(self.plates[i].get_left())
            y0 = int(self.plates[i].get_top())
            x1 = x0 + self.plates[i].get_width()
            y1 = y0 + self.plates[i].get_height()

            plate_map = self.plates[i].get_map(True, False)
            plate_age = self.plates[i].get_map(False, True)

            # Copy first part of plate onto world map.
            j = 0
            for y in range(y0, y1):
                for x in range(x0, x1):
                    y_index, x_index = plate_map.get_indeces(j)
                    x_mod = self.wd.x_mod(x)
                    y_mod = self.wd.y_mod(y)

                    k = self.wd.index_of(x_mod, y_mod)

                    # No crust here...
                    if plate_map._data[y_index, x_index] < 2 * FLT_EPSILON:
                        j += 1
                        continue

                    if self.imap[k] >= self.num_plates:  # No one here yet?
                        # This plate becomes the "owner" of current location
                        # if it is the first plate to have crust on it.
                        self.hmap._data[y_mod,
                                        x_mod] = plate_map._data[y_index, x_index]
                        self.imap[k] = i
                        self.amap._data[y_mod,
                                        x_mod] = plate_age._data[y_index, x_index]
                        j += 1
                        continue

                    # DO NOT ACCEPT HEIGHT EQUALITY! Equality leads to subduction
                    # of shore that 's barely above sea level. It's a lot less
                    # serious problem to treat very shallow waters as continent...
                    prev_is_oceanic = self.hmap._data[y_mod,
                                                      x_mod] < CONTINENTAL_BASE
                    this_is_oceanic = plate_map._data[y_index,
                                                      x_index] < CONTINENTAL_BASE

                    prev_timestamp = self.plates[self.imap[k]
                                                 ].get_crust_timestamp(x_mod, y_mod)
                    this_timestamp = plate_age._data[y_index, x_index]
                    prev_is_bouyant = (self.hmap._data[y_mod, x_mod] > plate_map._data[y_index, x_index]) | ((self.hmap._data[y_mod, x_mod] + 2 * FLT_EPSILON > plate_map._data[y_index, x_index]) & (
                        self.hmap._data[y_mod, x_mod] < 2 * FLT_EPSILON + plate_map._data[y_index, x_index]) & (prev_timestamp >= this_timestamp))

                    # Handle subduction of oceanic crust as special case.
                    if this_is_oceanic and prev_is_bouyant:
                        # This plate will be the subducting one.
                        # The level of effect that subduction has
                        # is directly related to the amount of water
                        # on top of the subducting plate.
                        sediment = SUBDUCT_RATIO * OCEANIC_BASE * \
                            (CONTINENTAL_BASE -
                             plate_map._data[y_index, x_index]) / CONTINENTAL_BASE

                        # Save collision to the receiving plate's list.
                        coll = PlateCollision(i, x_mod, y_mod, sediment)
                        self.subductions[self.imap[k]].append(coll)
                        oceanic_collisions += 1

                        # Remove subducted oceanic lithosphere from plate.
                        # This is crucial for
                        # a) having correct amount of colliding crust(below)
                        # b) protecting subducted locations from receiving
                        # crust from other subductions/collisions.
                        self.plates[i].set_crust(
                            x_mod, y_mod, plate_map._data[y_index, x_index] - OCEANIC_BASE, this_timestamp)

                        if plate_map._data[y_index, x_index] <= 0:
                            j += 1
                            continue  # Nothing more to collide.

                    elif prev_is_oceanic:
                        sediment = SUBDUCT_RATIO * OCEANIC_BASE * \
                            (CONTINENTAL_BASE -
                             self.hmap._data[y_mod, x_mod]) / CONTINENTAL_BASE

                        coll = PlateCollision(
                            self.imap[k], x_mod, y_mod, sediment)
                        self.subductions[i].append(coll)
                        oceanic_collisions += 1

                        self.plates[self.imap[k]].set_crust(
                            x_mod, y_mod, self.hmap._data[y_mod, x_mod] - OCEANIC_BASE, prev_timestamp)
                        self.hmap._data[y_mod, x_mod] -= OCEANIC_BASE

                        if self.hmap._data[y_mod, x_mod] <= 0:
                            self.imap[k] = i
                            self.hmap._data[y_mod,
                                            x_mod] = plate_map._data[y_index, x_index]
                            self.amap._data[y_mod,
                                            x_mod] = plate_age._data[y_index, x_index]
                            j += 1
                            continue

                    # Record collisions to both plates. This also creates
                    # continent segment at the collided location to plates.
                    this_area = self.plates[i].add_collision(x_mod, y_mod)
                    prev_area = self.plates[self.imap[k]
                                            ].add_collision(x_mod, y_mod)

                    # At least two plates are at same location.
                    # Move some crust from the SMALLER plate onto LARGER one.
                    if this_area < prev_area:
                        coll = PlateCollision(
                            self.imap[k], x_mod, y_mod, plate_map._data[y_index, x_index] * self.folding_ratio)

                        # Give some...
                        self.hmap._data[y_mod, x_mod] += coll.crust
                        self.plates[self.imap[k]].set_crust(
                            x_mod, y_mod, self.hmap._data[y_mod, x_mod], plate_age._data[y_index, x_index])

                        # And take some.
                        self.plates[i].set_crust(
                            x_mod, y_mod, plate_map._data[y_index, x_index] * (1.0 - self.folding_ratio), plate_age._data[y_index, x_index])

                        # Add collision to the earlier plate's list.
                        self.collisions[i].append(coll)
                        continental_collisions += 1
                    else:
                        coll = PlateCollision(
                            i, x_mod, y_mod, self.hmap._data[y_mod, x_mod] * self.folding_ratio)

                        self.plates[i].set_crust(
                            x_mod, y_mod, plate_map._data[y_index, x_index] + coll.crust, self.amap._data[y_mod, x_mod])

                        self.plates[self.imap[k]].set_crust(
                            x_mod, y_mod, self.hmap._data[y_mod, x_mod] * (1.0 - self.folding_ratio), self.amap._data[y_mod, x_mod])

                        self.collisions[self.imap[k]].append(coll)
                        continental_collisions += 1

                        # Give the location to the larger plate.
                        self.hmap._data[y_mod,
                                        x_mod] = plate_map._data[y_index, x_index]
                        self.imap[k] = i
                        self.amap._data[y_mod,
                                        x_mod] = plate_age._data[y_index, x_index]

                    j += 1

        self.last_coll_count = (self.last_coll_count + 1) & - \
            (continental_collisions == 0)

        for i in range(self.num_plates):
            for j in range(len(self.subductions[i])):
                coll = self.subductions[i][j]

                # ifdef DEBUG
                if i == coll._index:
                    print("when subducting: SRC == DEST!")
                    exit(1)
                # endif

                # Do not apply friction to oceanic plates.
                # This is a very cheap way to emulate slab pull.
                # Just perform subduction and on our way we go!
                self.plates[i].add_crust_by_subduction(
                    coll.wx, coll.wy, coll.crust, self.iter_count, self.plates[coll._index].get_velocity_x(), self.plates[coll._index].get_velocity_y())

            self.subductions[i] = []

        for i in range(self.num_plates):
            for j in range(len(self.collisions[i])):
                coll = self.collisions[i][j]
                coll_count = 0
                coll_count_i = 0
                coll_count_j = 0
                coll_ratio = 0.0
                coll_ratio_i = 0.0
                coll_ratio_j = 0.0

                # ifdef DEBUG
                if i == coll._index:
                    print("when colliding: SRC == DEST!")
                    exit(1)
                # endif

                # Collision causes friction. Apply it to both plates.
                self.plates[i].apply_friction(coll.crust)
                self.plates[coll._index].apply_friction(coll.crust)

                coll_count_i, coll_ratio_i = self.plates[i].get_collision_info(
                    coll.wx, coll.wy)
                coll_count_j, coll_ratio_j = self.plates[coll._index].get_collision_info(
                    coll.wx, coll.wy)

                # Find the minimum count of collisions between two
                # continents on different plates.
                # It's minimum because large plate will get collisions
                # from all over where as smaller plate will get just
                # a few. It's those few that matter between these two
                # plates, not what the big plate has with all the
                # other plates around it.
                coll_count = int(coll_count_i)
                coll_count -= (coll_count - coll_count_j) & - \
                    (coll_count > coll_count_j)

                # Find maximum amount of collided surface area between
                # two continents on different plates.
                # Like earlier, it's the "experience" of the smaller
                # plate that matters here.
                coll_ratio = float(coll_ratio_i)
                coll_ratio += (coll_ratio_j - coll_ratio) * \
                    (coll_ratio_j > coll_ratio)

                if (coll_count > self.aggr_overlap_abs) or (coll_ratio > self.aggr_overlap_rel):
                    amount = self.plates[i].aggregate_crust(
                        self.plates[coll._index], coll.wx, coll.wy)

                    # Calculate new direction and speed for the
                    # merged plate system, that is , for the
                    # receiving plate!
                    self.plates[coll._index].collide(
                        self.plates[i], coll.wx, coll.wy, amount)

            self.collisions[i] = []

        index_found = [0 for _ in range(self.num_plates)]

        # Fill divergent boundaries with new crustal material, molten magma.
        i = 0
        for y in range(BOOL_REGENERATE_CRUST * self.wd.get_height()):
            for x in range(self.wd.get_width()):
                if self.imap[i] >= self.num_plates:
                    # The owner of this new crust is that neighbour plate
                    # who was located at this point before plates moved.
                    self.imap[i] = prev_imap[i]

                    # ifdef DEBUG
                    if self.imap[i] >= self.num_plates:
                        print("Previous index map has no owner!")
                        exit(1)
                    # endif

                    # If this is oceanic crust then add buoyancy to it.
                    # Magma that has just crystallized into oceanic crust
                    # is more buoyant than that which has had a lot of
                    # time to cool down and become more dense.
                    self.hmap._data[y, x] = self.iter_count
                    self.hmap._data[y, x] = OCEANIC_BASE * BUOYANCY_BONUS_X

                    self.plates[self.imap[i]].set_crust(
                        x, y, OCEANIC_BASE, self.iter_count)

                else:
                    index_found[self.imap[i]] += 1
                    if self.hmap._data[y, x] <= 0.0:
                        print("Occupied point has no land mass!")
                        exit(1)
                i += 1

        # Remove empty plates from the system.
        i = 0
        while i < (self.num_plates):
            if (self.num_plates == 1):
                print("Only one plate left!")

            elif (index_found[i] == 0):
                self.plates[i] = self.plates.pop()
                index_found[i] = index_found.pop()

                # Life is seldom as simple as seems at first.
                # Replace the moved plate's index in the index map
                # to match its current position in the array!
                for j in range(self.wd.get_area()):
                    if self.imap[j] == self.num_plates - 1:
                        self.imap[j] = i

                self.num_plates -= 1
                i -= 1
            i += 1

        for i in range((BUOYANCY_BONUS_X > 0) * map_area):
            # Calculate the inverted age of this piece of crust.
            # Force result to be minimum between inv. age and
            # max buoyancy bonus age.
            y_index, x_index = self.amap.get_indeces(i)
            crust_age = int(self.iter_count -
                            self.amap._data[y_index, x_index])
            crust_age = MAX_BUOYANCY_AGE - crust_age
            crust_age = crust_age & -(crust_age <= MAX_BUOYANCY_AGE)

            self.hmap._data[y_index, x_index] += (self.hmap._data[y_index, x_index] < CONTINENTAL_BASE) * BUOYANCY_BONUS_X * \
                OCEANIC_BASE * crust_age * MULINV_MAX_BUOYANCY_AGE

        self.iter_count += 1

    def get_width(self):
        return self.wd.get_width()

    def get_height(self):
        return self.wd.get_height()

    def is_finished(self):
        return

    def restart(self):
        """Replace plates with a new population.
        """
        map_area = self.wd.get_area()

        # No increment if running for ever.
        self.cycle_count += self.max_cycles > 0
        if self.cycle_count > self.max_cycles:
            return

        # Update height map to include all recent changes.
        self.hmap.set_all(0)
        for i in range(self.num_plates):
            x0 = int(self.plates[i].get_left())
            y0 = int(self.plates[i].get_top())
            x1 = x0 + self.plates[i].get_width()
            y1 = y0 + self.plates[i].get_height()

            plate_map = self.plates[i].get_map(True, False)
            plate_age = self.plates[i].get_map(False, True)

            # Copy first part of plate onto world map.
            j = 0
            for y in range(y0, y1):
                for x in range(x0, x1):
                    x_mod = self.wd.x_mod(x)
                    y_mod = self.wd.y_mod(y)
                    h0 = self.hmap._data[y_mod, x_mod]
                    h1 = plate_map._data[y, x]
                    a0 = self.amap._data[x_mod, y_mod]
                    a1 = plate_age._data[y, x]

                    self.amap._data[x_mod, y_mod] = (
                        h0 * a0 + h1 * a1) / (h0 + h1)
                    self.hmap._data[y_mod, x_mod] += plate_map._data[y, x]
                    j += 1

        # Delete plates.
        self.plates = []
        self.num_plates = 0

        # create new plates IFF there are cycles left to run!
        # However, if max cycle count is "ETERNITY", then 0 < 0 + 1 always.
        eternity = int(not bool(self.max_cycles))
        if (self.cycle_count < self.max_cycles + eternity):
            self.create_plates(num_plates=self.max_plates)

            # Restore the ages of plates' points of crust!
            for i in range(self.num_plates):
                x0 = int(self.plates[i].get_left())
                y0 = int(self.plates[i].get_top())
                x1 = x0 + self.plates[i].get_width()
                y1 = y0 + self.plates[i].get_height()

                plate_map = self.plates[i].get_map(True, False)
                plate_age = self.plates[i].get_map(False, True)
                plate_age_const = self.plates[i].get_map(False, True)

                j = 0
                for y in range(y0, y1):
                    for x in range(x0, x1):
                        x_mod = self.wd.x_mod(x)
                        y_mod = self.wd.y_mod(y)

                        plate_age[j] = self.amap._data[y_mod, x_mod]

            return

        # Add some "virginity buoyancy" to all pixels for a visual boost.
        for i in range((BUOYANCY_BONUS_X > 0) * map_area):
            y_index, x_index = self.amap.get_indeces(i)
            crust_age = self.iter_count - self.amap._data[y_index, x_index]
            crust_age = MAX_BUOYANCY_AGE - crust_age
            crust_age &= -(crust_age <= MAX_BUOYANCY_AGE)

            self.hmap._data[y_index, x_index] += (self.hmap._data[y_index, x_index] < CONTINENTAL_BASE) * BUOYANCY_BONUS_X * \
                OCEANIC_BASE * crust_age * MULINV_MAX_BUOYANCY_AGE

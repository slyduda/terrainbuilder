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

from numpy import array, zeros, int32
from rectangle import WorldDimension
from numpy.random import RandomState, SeedSequence
from noise import generate_noise
from plate import Plate

from maps import HeightMap, AgeMap

SUBDUCT_RATIO = 0.5

BUOYANCY_BONUS_X = 3
MAX_BUOYANCY_AGE = 20
MULINV_MAX_BUOYANCY_AGE = 1.0 / MAX_BUOYANCY_AGE

RESTART_ENERGY_RATIO = 0.15
RESTART_SPEED_LIMIT = 2.0
NO_COLLISION_TIME_LIMIT = 10

# NEEDS WORK
CONTINENTAL_BASE = 0
OCEANIC_BASE = 0
FLT_EPSILON = 0
BOOL_REGENERATE_CRUST = 0


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

    def __init__(self, seed: int, width: int, height: int, sea_level: float, _erosion_period: int, _folding_ratio: float, aggr_ratio_abs: int, aggr_ration_rel: float, num_cycles: int):
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
        self._randstate = RandomState(seed=seed)
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

        tmp_dim = WorldDimension(width, height)
        A = tmp_dim.get_area()
        tmp = [0.0 for _ in range(A)]

        tmp = self.create_noise(tmp, tmp_dim, True)

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
            tmp[i] = (tmp[i] > sea_level) * (tmp[i] + CONTINENTAL_BASE) + \
                (tmp[i] <= sea_level) * OCEANIC_BASE

        self.peak_ek = 0.0
        self.last_coll_count = 0

        self.collisions = None
        self.subductions = None

        self.imap = [0 for _ in range(self.wd.get_area())]

    def create_noise(self, tmp: list, tmp_dim: WorldDimension, use_simplex: bool = False):
        return generate_noise(tmp, tmp_dim, self._randstate, use_simplex)

    def create_plates(self, num_plates):
        map_area = self.wd.get_area()
        self.max_plates = num_plates
        self.num_plates = num_plates

        vec = []  # for plate collision objects
        vec_size = self.wd.get_max()*4  # Not going to use this yet

        self.collisions = []
        self.subductions = []

        for i in range(num_plates):
            self.collisions.append(vec)
            self.subductions.append(vec)

        for i in range(map_area):
            self.imap[i] = i

        area = [PlateArea() for _ in num_plates]
        for i in range(num_plates):
            # Randomly select an unused plate origin.
            p = self.imap[self._randstate.randint(1, int32) % (map_area - i)]
            y = self.wd.y_from_index(p)
            x = self.wd.x_from_index(p)

            area[i].lft = area[i].rgt = x  # Save origin...
            area[i].top = area[i].btm = y
            area[i].wdt = area[i].hgt = 1

            area[i].border.append(p)  # ...and mark it as border.

            # Overwrite used entry with last unused entry in array.
            self.imap[p] = self.imap[map_area - i - 1]

        owner = self.imap  # Create an alias.
        # TODO Might have to copy the imap

        # "Grow" plates from their origins until surface is fully populated.
        max_border = 1
        i = 0
        while max_border:
            max_border = i
            for i in range(num_plates):
                N = len(area[i].border)
                max_border = max_border if max_border > N else N

                if N == 0:
                    continue

                j = self._randstate.randint(1, int32) % N
                p = area[i].border[j]
                cy = self.wd.y_from_index(p)
                cx = self.wd.x_from_index(p)

                lft = cx - 1 if cx > 0 else self.wd.get_width() - 1
                rgt = cx + 1 if cx < self.wd.get_width() - 1 else 0
                top = cy - 1 if cy > 0 else self.wd.get_height() - 1
                btm = cy + 1 if cy < self.wd.get_height() - 1 else 0

                n = top * self.wd.get_width() + cx  # North.
                s = btm * self.wd.get_width() + cx  # South.
                w = cy * self.wd.get_width() + lft  # West.
                e = cy * self.wd.get_width() + rgt  # East.

                if owner[n] >= num_plates:
                    owner[n] = i
                    area[i].border.push_back(n)

                    if area[i].top == self.wd.y_mod(top + 1):
                        area[i].top = top
                        area[i].hgt += 1

                if owner[s] >= num_plates:
                    owner[s] = i
                    area[i].border.push_back(s)

                    if btm == self.wd.y_mod(area[i].btm + 1):
                        area[i].btm = btm
                        area[i].hgt += 1

                if owner[w] >= num_plates:
                    owner[w] = i
                    area[i].border.push_back(w)

                    if area[i].lft == self.wd.x_mod(lft + 1):
                        area[i].lft = lft
                        area[i].wdt += 1

                if owner[e] >= num_plates:
                    owner[e] = i
                    area[i].border.push_back(e)

                    if rgt == self.wd.x_mod(area[i].rgt + 1):
                        area[i].rgt = rgt
                        area[i].wdt += 1

                # Overwrite processed point with unprocessed one.
                area[i].border[j] = area[i].border.back()
                area[i].border.pop_back()

        self.plates = [None for _ in num_plates]

        # Extract and create plates from initial terrain.
        for i in range(num_plates):
            area[i].wdt = self.wd.x_cap(area[i].wdt)
            area[i].hgt = self.wd.y_cap(area[i].hgt)

            x0 = area[i].lft
            x1 = 1 + x0 + area[i].wdt
            y0 = area[i].top
            y1 = 1 + y0 + area[i].hgt
            width = x1 - x0
            height = y1 - y0
            plt = [0.0 for _ in width * height]

            # Copy plate's height data from global map into local map.
            j = 0
            for y in range(y0, y1):
                for x in range(x0, x1):
                    k = self.wd.normalized_index_of(x, y)
                    plt[j] = self.hmap[k] * (owner[k] == i)
                    j += 1

            # Create plate. TODO Fix the SeedSeq
            self.plates[i] = Plate(self._randstate.randint(1, int32), plt, width,
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
        self.imap = [0 for _ in map_area]

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
                    x_mod = self.wd.x_mod(x)
                    y_mod = self.wd.y_mod(y)

                    k = self.wd.index_of(x_mod, y_mod)

                    if plate_map[j] < 2 * FLT_EPSILON:  # No crust here...
                        j += 1
                        continue

                    if self.imap[k] >= self.num_plates:  # No one here yet?
                        # This plate becomes the "owner" of current location
                        # if it is the first plate to have crust on it.
                        self.hmap[k] = plate_map[j]
                        self.imap[k] = i
                        self.amap[k] = plate_age[j]
                        j += 1
                        continue

                    # DO NOT ACCEPT HEIGHT EQUALITY! Equality leads to subduction
                    # of shore that 's barely above sea level. It's a lot less
                    # serious problem to treat very shallow waters as continent...
                    prev_is_oceanic = self.hmap[k] < CONTINENTAL_BASE
                    this_is_oceanic = plate_map[j] < CONTINENTAL_BASE

                    prev_timestamp = self.plates[self.imap[k]
                                                 ].get_crust_timestamp(x_mod, y_mod)
                    this_timestamp = plate_age[j]
                    prev_is_bouyant = (self.hmap[k] > plate_map[j]) | ((self.hmap[k] + 2 * FLT_EPSILON > plate_map[j]) & (
                        self.hmap[k] < 2 * FLT_EPSILON + plate_map[j]) & (prev_timestamp >= this_timestamp))

                    # Handle subduction of oceanic crust as special case.
                    if this_is_oceanic and prev_is_bouyant:
                        # This plate will be the subducting one.
                        # The level of effect that subduction has
                        # is directly related to the amount of water
                        # on top of the subducting plate.
                        sediment = SUBDUCT_RATIO * OCEANIC_BASE * \
                            (CONTINENTAL_BASE -
                             plate_map[j]) / CONTINENTAL_BASE

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
                            x_mod, y_mod, plate_map[j] - OCEANIC_BASE, this_timestamp)

                        if plate_map[j] <= 0:
                            j += 1
                            continue  # Nothing more to collide.

                    elif prev_is_oceanic:
                        sediment = SUBDUCT_RATIO * OCEANIC_BASE * \
                            (CONTINENTAL_BASE -
                             self.hmap[k]) / CONTINENTAL_BASE

                        coll = PlateCollision(
                            self.imap[k], x_mod, y_mod, sediment)
                        self.subductions[i].push_back(coll)
                        oceanic_collisions += 1

                        self.plates[self.imap[k]].set_crust(
                            x_mod, y_mod, self.hmap[k] - OCEANIC_BASE, prev_timestamp)
                        self.hmap[k] -= OCEANIC_BASE

                        if self.hmap[k] <= 0:
                            self.imap[k] = i
                            self.hmap[k] = plate_map[j]
                            self.amap[k] = plate_age[j]
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
                            self.imap[k], x_mod, y_mod, plate_map[j] * self.folding_ratio)

                        # Give some...
                        self.hmap[k] += coll.crust
                        self.plates[self.imap[k]].set_crust(
                            x_mod, y_mod, self.hmap[k], plate_age[j])

                        # And take some.
                        self.plates[i].set_crust(
                            x_mod, y_mod, plate_map[j] * (1.0 - self.folding_ratio), plate_age[j])

                        # Add collision to the earlier plate's list.
                        self.collisions[i].append(coll)
                        continental_collisions += 1
                    else:
                        coll = PlateCollision(
                            i, x_mod, y_mod, self.hmap[k] * self.folding_ratio)

                        self.plates[i].set_crust(
                            x_mod, y_mod, plate_map[j] + coll.crust, self.amap[k])

                        self.plates[self.imap[k]].set_crust(
                            x_mod, y_mod, self.hmap[k] * (1.0 - self.folding_ratio), self.amap[k])

                        self.collisions[self.imap[k]].append(coll)
                        continental_collisions += 1

                        # Give the location to the larger plate.
                        self.hmap[k] = plate_map[j]
                        self.imap[k] = i
                        self.amap[k] = plate_age[j]

                    j += 1

        self.last_coll_count = (self.last_coll_count + 1) & - \
            (continental_collisions == 0)

        for i in range(self.num_plates):
            for j in range(len(self.subductions[i])):
                coll = self.subductions[i][j]

                # ifdef DEBUG
                if i == coll.index:
                    print("when subducting: SRC == DEST!")
                    exit(1)
                # endif

                # Do not apply friction to oceanic plates.
                # This is a very cheap way to emulate slab pull.
                # Just perform subduction and on our way we go!
                self.plates[i].add_crust_by_subduction(
                    coll.wx, coll.wy, coll.crust, self.iter_count, self.plates[coll.index].get_velocity_x(), self.plates[coll.index].get_velocity_y())
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
                if i == coll.index:
                    print("when colliding: SRC == DEST!")
                    exit(1)
                # endif

                # Collision causes friction. Apply it to both plates.
                self.plates[i].apply_friction(coll.crust)
                self.plates[coll.index].apply_friction(coll.crust)

                coll_count_i, coll_ratio_i = self.plates[i].get_collision_info(
                    coll.wx, coll.wy, coll_count_i, coll_ratio_i)
                coll_count_j, coll_ratio_j = self.plates[coll.index].get_collision_info(
                    coll.wx, coll.wy, coll_count_j, coll_ratio_j)

                # Find the minimum count of collisions between two
                # continents on different plates.
                # It's minimum because large plate will get collisions
                # from all over where as smaller plate will get just
                # a few. It's those few that matter between these two
                # plates, not what the big plate has with all the
                # other plates around it.
                coll_count = coll_count_i
                coll_count -= (coll_count - coll_count_j) & - \
                    (coll_count > coll_count_j)

                # Find maximum amount of collided surface area between
                # two continents on different plates.
                # Like earlier, it's the "experience" of the smaller
                # plate that matters here.
                coll_ratio = coll_ratio_i
                coll_ratio += (coll_ratio_j - coll_ratio) * \
                    (coll_ratio_j > coll_ratio)

                if (coll_count > self.aggr_overlap_abs) | (coll_ratio > self.aggr_overlap_rel):
                    amount = self.plates[i].aggregate_crust(
                        self.plates[coll.index], coll.wx, coll.wy)

                    # Calculate new direction and speed for the
                    # merged plate system, that is , for the
                    # receiving plate!
                    self.plates[coll.index].collide(
                        self.plates[i], coll.wx, coll.wy, amount)

            self.collisions[i] = []

        index_found = [0 for _ in self.num_plates]

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
                    self.amap[i] = self.iter_count
                    self.hmap[i] = OCEANIC_BASE * BUOYANCY_BONUS_X

                    self.plates[self.imap[i]].set_crust(
                        x, y, OCEANIC_BASE, self.iter_count)

                # TODO DOUBLE CHECK THIS OPErATER
                elif index_found[self.imap[i]] + 1 and self.hmap[i] <= 0:
                    index_found[self.imap[i]] += 1
                    print("Occupied point has no land mass!")
                    exit(1)

        # Remove empty plates from the system.
        i = 0
        while i < (self.num_plates):
            if (self.num_plates == 1):
                print("Only one plate left!")

            elif (index_found[i] == 0):
                self.plates[i] = self.plates[self.num_plates - 1]
                index_found[i] = index_found[self.num_plates - 1]

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
            crust_age = self.iter_count - self.amap[i]
            crust_age = MAX_BUOYANCY_AGE - crust_age
            crust_age &= -(crust_age <= MAX_BUOYANCY_AGE)

            self.hmap[i] += (self.hmap[i] < CONTINENTAL_BASE) * BUOYANCY_BONUS_X * \
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
                h0 = self.hmap[self.wd.index_of(x_mod, y_mod)]
                h1 = plate_map[j]
                a0 = self.amap[self.wd.index_of(x_mod, y_mod)]
                a1 = plate_age[j]

                self.amap[self.wd.index_of(x_mod, y_mod)] = (
                    h0 * a0 + h1 * a1) / (h0 + h1)
                self.hmap[self.wd.index_of(x_mod, y_mod)] += plate_map[j]
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

                        plate_age[j] = self.amap[self.wd.index_of(
                            x_mod, y_mod)]

            return

        # Add some "virginity buoyancy" to all pixels for a visual boost.
        for i in range((BUOYANCY_BONUS_X > 0) * map_area):
            crust_age = self.iter_count - self.amap[i]
            crust_age = MAX_BUOYANCY_AGE - crust_age
            crust_age &= -(crust_age <= MAX_BUOYANCY_AGE)

            self.hmap[i] += (self.hmap[i] < CONTINENTAL_BASE) * BUOYANCY_BONUS_X * \
                OCEANIC_BASE * crust_age * MULINV_MAX_BUOYANCY_AGE

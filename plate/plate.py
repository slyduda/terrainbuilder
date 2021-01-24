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


from __future__ import annotations
from numpy import pi, cos, sin, zeros, int32, sqrt, nditer
from numpy.random import RandomState, SeedSequence, MT19937
from .maps import HeightMap, AgeMap, MassMap
from .rectangle import WorldDimension, Rectangle

INITIAL_SPEED_X = 1
DEFORMATION_WEIGHT = 2

# TODO
CONT_BASE = 0


class ContinentId(int):
    def __new__(cls, arg):
        return super(ContinentId, cls).__new__(cls, arg)


class Plate(object):
    def __init__(self, seed: SeedSequence, m: MassMap, w: int, h: int, _x: int, _y: int, plate_age: int, wd: WorldDimension):
        """Initializes plate with the supplied height map.

        Args:
            seed (SeedSequence): [description]
            m (MassMap): Array to height map of terrain.
            w (int): Width of height map in pixels.
            h (int): Height of height map in pixels.
            _x (int): X of height map's left-top corner on world map.
            _y (int): Y of height map's left-top corner on world map.
            wd (WorldDimension): Dimension of world map's either side in pixels.
        """
        self._randstate = RandomState(MT19937(seed))
        self.width = w
        self.height = h
        self.mass = 0
        self.left = _x
        self.top = _y
        self.cx = 0
        self.cy = 0
        self.dx = 0
        self.dy = 0
        self.map = HeightMap(h, w)
        self.age_map = AgeMap(h, w)
        self.wd = wd

        plate_area = w * h
        self.segment = zeros(plate_area, dtype=int)
        self.segment[:] = 255
        self.seg_data = []

        angle = 2 * pi * self._randstate.random_sample()  # ?
        self.vx = cos(angle) * INITIAL_SPEED_X
        self.vy = sin(angle) * INITIAL_SPEED_X

        self.velocity = 1
        self.rot_dir = -1 if (self._randstate.randint(2) == 0) else 1

        for y in range(self.height):
            for x in range(self.width):
                k = y * self.width + x
                self.map._data[y, x] = m[k]
                self.mass += m[k]

                # Calculate center coordinates weighted by mass.
                self.cx += x * m[k]
                self.cy += y * m[k]

                # Set the age of ALL points in this plate to same
                # value. The right thing to do would be to simulate
                # the generation of new oceanic crust as if the plate
                # had been moving to its current direction until all
                # plate's (oceanic) crust receive an age.
                self.age_map._data[y, x] = plate_age if m[k] > 0 else 0

        self.cx /= self.mass
        self.cy /= self.mass

    def add_collision(self, wx: int, wy: int):
        """Increment collision counter of the continent at given location.

        Args:
            wx (int): X coordinate of collision point on world map.
            wy (int): Y coordinate of collision point on world map.

        Returns:
            int: Surface area of the collided continent
        """
        seg = self.get_continent_at(wx, wy)
        self.seg_data[seg].coll_count += 1
        return self.seg_data[seg].area

    def add_crust_by_collision(self, x: int, y: int, z: float, time: int, activeContinent: ContinentId):
        """Add crust to plate as result of continental collision.

        Args:
            x (int): Location of new crust on global world map (X).
            y (int): Location of new crust on global world map (Y).
            z (float): Amount of crust to add.
            time (int): Time of creation of new crust.
            activeContinent (ContinentId): Segment ID of the continent that's processed.
        """
        self.set_crust(x, y, self.get_crust(x, y) + z, time)

        index = self.get_map_index(x, y)
        lx, ly = self.get_map_indeces(x, y)

        self.segment[index] = activeContinent
        data = self.seg_data[activeContinent]

        data.area += 1
        data.enlarge_to_contain(lx, ly)

    def add_crust_by_subduction(self, x: int, y: int, z: float, t: int, dx: float, dy: float):
        """Simulates subduction of oceanic plate under this plate.

        Subduction is simulated by calculating the distance on surface
        that subducting sediment will travel under the plate until the
        subducting slab has reached certain depth where the heat triggers
        the melting and uprising of molten magma.

        Args:
            x (int): Origin of subduction on global world map (X).
            y (int): Origin of subduction on global world map (Y).
            z (float): Amount of sediment that subducts.
            t (int): Time of creation of new crust.
            dx (float): Direction of the subducting plate (X).
            dy (float): Direction of the subducting plate (Y).
        """
        lx, ly = self.get_map_indeces(x, y)

        # Take vector difference only between plates that move more or less
        # to same direction. This makes subduction direction behave better.

        # Use of "this" pointer is not necessary, but it make code clearer.
        # Cursed be those who use "m_" prefix in member names! >(
        dot = self.vx * dx + self.vy * dy
        dx -= self.vx * int(dot > 0)
        dy -= self.vy * int(dot > 0)

        offset = self._randstate.random_sample()
        offset_sign = 2 * (self._randstate.randint(1, high=2147483648) % 2) - 1
        offset *= offset * offset * offset_sign
        offset2 = self._randstate.random_sample()
        offset_sign2 = 2 * \
            (self._randstate.randint(1, high=2147483648) % 2) - 1
        offset2 *= offset2 * offset2 * offset_sign2
        dx = 10 * dx + 3 * offset
        dy = 10 * dy + 3 * offset2

        fx = lx + dx
        fy = ly + dy

        if fx < self.get_width() and fx >= 0 and fy < self.get_height() and fy >= 0:
            fx = int(fx)
            fy = int(fy)

            if self.map._data[fy, fx] > 0:
                t = int(self.map._data[fy, fx] * self.age_map._data[fy, fx] +
                        z * t) / (self.map._data[fy, fx] + z)
                self.age_map._data[fy, fx] = t * int(z > 0)

                self.map._data[fy, fx] += z
                self.mass += z

    def aggregate_crust(self, p: Plate, wx: int, wy: int):
        """Add continental crust from this plate as part of other plate.

        Aggregation of two continents is the event where the collided
        pieces of crust fuse together at the point of collision. It is
        crucial to merge not only the collided pieces of crust but also
        the entire continent that's part of the colliding tad of crust
        However, because one plate can contain many islands and pieces of
        continents, the merging must be done WITHOUT merging the entire
        plate and all those continental pieces that have NOTHING to do with
        the collision in question.

        Args:
            p (Plate):  Pointer to the receiving plate.
            wx (int):   X coordinate of collision point on world map.
            wy (int):   Y coordinate of collision point on world map.

        Returns:
            (float): Amount of crust aggregated to destination plate.
        """
        index = self.get_map_index(wx, wy)
        lx, ly = self.get_map_indeces(wx, wy)
        seg_id = self.segment[index]

        # This check forces the caller to do things in proper order!

        # Usually continents collide at several locations simultaneously.
        # Thus if this segment that is being merged now is removed from
        # segmentation bookkeeping, then the next point of collision that is
        # processed during the same iteration step would cause the test
        # below to be true and system would experience a premature abort.

        # Therefore, segmentation bookkeeping is left intact. It doesn't
        # cause significant problems because all crust is cleared and empty
        # points are not processed at all. (Test on (seg_id >= seg_data.size()) removed)

        # One continent may have many points of collision. If one of them
        # causes continent to aggregate then all successive collisions and
        # attempts of aggregation would necessarily change nothing at all,
        # because the continent was removed from this plate earlier!

        if not self.seg_data[seg_id]:
            return 0.0   # Do not process empty continents.

        activeContinent = self.select_collision_segment(wx, wy)

        # Wrap coordinates around world edges to safeguard subtractions.
        wx += self.wd.get_width()
        wy += self.wd.get_height()

        # Aggregating segment [%u, %u]x[%u, %u] vs. [%u, %u]@[%u, %u]\n",
        # seg_data[seg_id].x0, seg_data[seg_id].y0,
        # seg_data[seg_id].x1, seg_data[seg_id].y1,
        # width, height, lx, ly);

        old_mass = float(self.mass)

        # Add all of the collided continent's crust to destination plate.
        y_start = self.seg_data[seg_id].get_top()
        y_end = self.seg_data[seg_id].get_bottom()
        for y in range(y_start, y_end + 1):

            x_start = self.seg_data[seg_id].get_left()
            x_end = self.seg_data[seg_id].get_right()
            for x in range(x_start, x_end + 1):

                i = y * self.width + x
                if self.segment[i] == seg_id and self.map._data[y, x] > 0:
                    self.add_crust_by_collision(wx + x - lx, wy + y - ly,
                                                self.map._data[y, x], self.age_map._data[y, x], activeContinent)

                    self.mass -= self.map._data[y, x]
                    self.map._data[y, x] = 0

        self.seg_data[seg_id].area = 0  # Mark segment as non-existent
        return old_mass - self.mass

    def apply_friction(self, deformed_mass: float):
        """Decrease the speed of plate amount relative to its total mass.

        Method decreses the speed of plate due to friction that occurs when
        two plates collide. The amount of reduction depends of the amount
        of mass that causes friction (i.e. that has collided) compared to
        the total mass of the plate. Thus big chunk of crust colliding to
        a small plate will halt it but have little effect on a huge plate.

        Args:
            deformed_mass (float): Amount of mass deformed in collision.
        """
        if self.mass > 0:
            vel_dec = DEFORMATION_WEIGHT * deformed_mass / self.mass
            vel_dec = vel_dec if (vel_dec < self.velocity) else self.velocity

            # Altering the source variable causes the order of calls to
            # this function to have difference when it shouldn't!
            # However, it's a hack well worth the outcome. :)
            self.velocity -= vel_dec

        return

    def collide(self, p: Plate, wx: int, wy: int, coll_mass: float):
        """Method collides two plates according to Newton's laws of motion.

        The velocity and direction of both plates are updated using
        impulse forces following the collision according to Newton's laws
        of motion. Deformations are not applied but energy consumed by the
        deformation process is taken away from plate's momentum.

        Args:
            p (Plate): Plate to test against.
            wx (int): X coordinate of collision point on world map.
            wy (int): Y coordinate of collision point on world map.
            coll_mass (float): Amount of colliding mass from source plate.
        """
        coeff_rest = 0.0

        ap_dx, ap_dy, bp_dx, bp_dy, nx, ny = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        index = self.get_map_index(wx, wy)
        apx, apy = self.get_map_indeces(wx, wy)
        p_index = p.get_map_index(wx, wy)
        bpx, bpy = p.get_map_indeces(wx, wy)

        # out of colliding map's bounds!
        assert(index < self.width * self.height), "Out of Map Bounds"
        assert(p_index < p.width * p.height), "Out of Map Bounds"

        ap_dx = int(apx) - int(self.cx)
        ap_dy = int(apy) - int(self.cy)
        bp_dx = int(bpx) - int(p.cx)
        bp_dy = int(bpy) - int(p.cy)
        nx = ap_dx - bp_dx
        ny = ap_dy - bp_dy

        if nx * nx + ny * ny <= 0:
            print("DONT DIVIDE BY ZERO")
            return
            # Avoid division by zero!

        # Scaling is required at last when impulses are added to plates!
        n_len = sqrt(nx * nx + ny * ny)
        nx /= n_len
        ny /= n_len

        # Compute relative velocity between plates at the collision point.
        # Because torque is not included, calc simplifies to v_ab = v_a - v_b.
        rel_vx = self.vx - p.vx
        rel_vy = self.vy - p.vy

        # Get the dot product of relative velocity vector and collision vector.
        # Then get the projection of v_ab along collision vector.
        # Note that vector n must be a unit vector!
        rel_dot_n = rel_vx * nx + rel_vy * ny

        if rel_dot_n <= 0:
            return
            # Exit if objects are moving away from each other.

        # Calculate the denominator of impulse: n . n * (1 / m_1 + 1 / m_2).
        # Use the mass of the colliding crust for the "donator" plate.
        denom = (nx * nx + ny * ny) * (1.0/self.mass + 1.0/coll_mass)

        # Calculate force of impulse.
        J = -(1 + coeff_rest) * rel_dot_n / denom

        # Compute final change of trajectory.
        # The plate that is the "giver" of the impulse should receive a
        # force according to its pre-collision mass, not the current mass!
        self.dx += nx * J / self.mass
        self.dy += ny * J / self.mass
        p.dx -= nx * J / (coll_mass + p.mass)
        p.dy -= ny * J / (coll_mass + p.mass)

    def erode(self, lower_bound: float):
        """Apply plate wide erosion algorithm.

        Plates total mass and the center of mass are updated.

        Args:
            lower_bound (float): Sets limit below which there's no erosion.
        """
        source_data = []
        sinks_data = []
        sources = source_data
        sinks = sinks_data

        tmp = float(self.width * self.height)

        # Find all tops.
        for y in range(self.height):
            for x in range(self.width):
                index = y * self.width + x

                if self.map._data[y, x] < lower_bound:
                    continue

                w_crust, e_crust, n_crust, s_crust = (0.0, 0.0, 0.0, 0.0)
                w, e, n, s = (0, 0, 0, 0)
                self.calculate_crust(x, y, index, w_crust, e_crust, n_crust, s_crust,
                                     w, e, n, s)

                # This location is either at the edge of the plate or it is not the
                # tallest of its neightbours. Don't start a river from here.
                if w_crust * e_crust * n_crust * s_crust == 0:
                    continue

                sources.pop(index)

        is_done = zeros(self.width*self.height)
        # while not sources:
        # while not sources:

        return

    def get_collision_info(self, wx: int, wy: int):
        """Retrieve collision statistics of continent at given location.

        Args:
            wx (int): X coordinate of collision point on world map.
            wy (int): Y coordinate of collision point on world map.

        Returns:
            count (int): Destination for the count of collisions.
            ration (float): Destination for the % of area that collided.
        """
        seg = self.get_continent_at(wx, wy)

        count = self.seg_data[seg].coll_count
        ratio = self.seg_data[seg].coll_count / (1 + self.seg_data[seg].area)
        # +1 avoids DIV with zero.
        return (count, ratio)

    def get_continent_area(self, wx: int, wy: int):
        """Retrieve the surface area of continent lying at desired location.

        Args:
            wx (int): X coordinate of collision point on world map.
            wy (int): Y coordinate of collision point on world map.

        Returns:
            int: Area of continent at desired location or 0 if none.
        """
        index = self.get_map_index(wx, wy)
        assert(self.segment[index] < len(self.seg_data)), 'Greater than'

        return self.seg_data[self.segment[index]].area

    def get_crust(self, x: int, y: int):
        """Get the amount of plate's crustal material at some location.

        Args:
            x (int): Offset on the global world map along X axis.
            y (int): Offset on the global world map along Y axis.

        Returns:
            float: Amount of crust at requested location.
        """
        # Might fail
        index = self.get_map_index(x, y)
        lx, ly = self.get_map_indeces(x, y)
        return self.map._data[ly, lx] if index < -1 else 0

    def get_crust_timestamp(self, x: int, y: int):
        """Get the timestamp of plate's crustal material at some location.

        Args:
            x (int): Offset on the global world map along X axis.
            y (int): Offset on the global world map along Y axis.

        Returns:
            int: Timestamp of creation of crust at the location.
            Zero is returned if location contains no crust.
        """
        # Might fail
        index = self.get_map_index(x, y)
        lx, ly = self.get_map_indeces(x, y)
        return self.age_map._data[ly, lx] if index < -1 else 0

    def get_map(self, c: bool = False, t: bool = False):
        """Get plate's data.

        Args:
            c (bool, optional): . Defaults to False.
            t (bool, optional): . Defaults to False.

        Returns:
            HeightMap or AgeMap: [description]
        """
        if (c):
            return self.map
        if (t):
            return self.age_map

    def move(self):
        """Moves plate along it's trajectory.
        """
        # Apply any new impulses to the plate's trajectory.
        self.vx += self.dx
        self.vy += self.dy
        self.dx = 0
        self.dy = 0

        # Force direction of plate to be unit vector.
        # Update velocity so that the distance of movement doesn't change.
        length = sqrt(self.vx*self.vx+self.vy*self.vy)
        self.vx /= length
        self.vy /= length
        self.velocity += length - 1.0
        self.velocity *= self.velocity > 0  # Round negative values to zero.

        # Apply some circular motion to the plate.
        # Force the radius of the circle to remain fixed by adjusting
        # angular velocity(which depends on plate's velocity).
        world_avg_side = (self.wd.get_width() +
                          self.wd.get_height()) / 2
        alpha = self.rot_dir * self.velocity / (world_avg_side * 0.33)
        _cos = cos(alpha * self.velocity)
        _sin = sin(alpha * self.velocity)
        _vx = self.vx * _cos - self.vy * _sin
        _vy = self.vy * _cos + self.vx * _sin
        self.vx = _vx
        self.vy = _vy

        # Location modulations into range[0..world width/height[are a have to!
        # If left undone SOMETHING WILL BREAK DOWN SOMEWHERE in the code!

        assert(self.wd.contains(self.left, self.top))

        self.left += self.vx * self.velocity
        self.left += 0 if self.left >= 0 else self.wd.get_width()
        self.left -= 0 if self.left < self.wd.get_width() else self.wd.get_width()

        self.top += self.vy * self.velocity
        self.top += 0 if self.top >= 0 else self.wd.get_height()
        self.top -= 0 if self.top < self.wd.get_height() else self.wd.get_height()

        assert(self.wd.contains(self.left, self.top)
               ), "World Dimension contains Plate"

    def reset_segments(self):
        """Clear any earlier continental crust partitions.

        Plate has an internal bookkeeping of distinct areas of continental
        crust for more realistic collision responce. However as the number
        of collisions that plate experiences grows, so does the bookkeeping
        of a continent become more and more inaccurate. Finally it results
        in striking artefacts that cannot overlooked.

        To alleviate this problem without the need of per iteration
        recalculations plate supplies caller a method to reset its
        bookkeeping and start clean.
        """
        self.seg_data.clear()
        self.segment = zeros(self.width * self.height, dtype=int)
        self.segment[:] = 255
        return

    def select_collision_segment(self, coll_x: int, coll_y: int):
        """Remember the currently processed continent's segment number.

        Args:
            coll_x (int): Origin of collision on global world map (X).
            coll_y (int): Origin of collision on global world map (Y).

        Returns:
            ContinentId: the Id of the continent being processed
        """
        index = self.get_map_index(coll_x, coll_y)
        activeContinent = self.segment[index]
        return activeContinent

    def set_crust(self, x: int, y: int, z: float, t: int):
        """Set the amount of plate's crustal material at some location.

        If amount of crust to be set is negative, it'll be set to zero.

        Args:
            x (int): Offset on the global world map along X axis.
            y (int): Offset on the global world map along Y axis.
            z (float): Amount of crust at given location.
            t (int): Time of creation of new crust.
        """
        if z < 0:
            z = 0

        index = self.get_map_index(x, y)
        lx, ly = self.get_map_indeces(x, y)

        if index >= self.width * self.height:
            assert(z > 0)

            ilft = int(self.left)
            itop = int(self.top)
            irgt = ilft + self.width - 1
            ibtm = itop + self.height - 1

            x, y = self.wd.normalize(x, y)

            # Might need update. I believe it is a bool to int conversion
            _lft = ilft - x
            _rgt = int(self.wd.get_width() & -(x < ilft)) + x - irgt
            _top = itop - y
            _btm = int(self.wd.get_height() & -(y < itop)) + y - ibtm

            d_lft = _lft & -(_lft < _rgt) & -(_lft < self.wd.get_width())
            d_rgt = _rgt & -(_rgt <= _lft) & -(_rgt < self.wd.get_width())
            d_top = _top & -(_top < _btm) & -(_top < self.wd.get_height())
            d_btm = _btm & -(_btm <= _top) & -(_btm < self.wd.get_height())

            d_lft = ((d_lft > 0) + (d_lft >> 3)) << 3
            d_rgt = ((d_rgt > 0) + (d_rgt >> 3)) << 3
            d_top = ((d_top > 0) + (d_top >> 3)) << 3
            d_btm = ((d_btm > 0) + (d_btm >> 3)) << 3

            if self.width + d_lft + d_rgt > self.wd.get_width():
                d_lft = 0
                d_rgt = self.wd.get_width() - self.width

            if self.height + d_top + d_btm > self.wd.get_height():
                d_top = 0
                d_btm = self.wd.get_height() - self.height

            assert(d_lft + d_rgt + d_top + d_btm != 0)

            old_width = self.width
            old_height = self.height

            self.left -= d_lft
            self.left += 0 if self.left >= 0 else self.wd.get_width()
            self.left -= 0 if self.left < self.wd.get_width() else self.wd.get_width()
            self.width += d_lft + d_rgt

            self.top -= d_top
            self.top += 0 if self.top >= 0 else self.wd.get_height()
            self.top -= 0 if self.top < self.wd.get_height() else self.wd.get_height()
            self.height += d_top + d_btm

            tmph = HeightMap(self.width, self.height)
            tmpa = AgeMap(self.width, self.height)
            tmps = zeros(self.width*self.height)
            tmph.set_all(0)
            tmpa.set_all(0)

            # We are supposed to create a new plate here
            for j in range(old_height):
                dest_i = (d_top + j) * self.width + d_lft
                src_i = j * old_width

            self.map = tmph
            self.age_map = tmpa
            self.segment = tmps

            for s in range(len(self.seg_data)):
                self.seg_data[s].shift(d_lft, d_top)

            assert(index < self.width * self.height)

        # Update crust's age.
        old_crust = self.map._data[ly, lx] > 0
        new_crust = z > 0
        # If old crust exists create new t with the mean of original and supplied
        if old_crust:
            t = int((self.map._data[ly, lx] * self.age_map._data[ly, lx] +
                     z * t) / (self.map._data[ly, lx] + z))
        if new_crust:
            self.age_map._data[ly, lx] = t

        self.mass -= self.map._data[ly, lx]
        self.map._data[ly, lx] = z
        self.mass += z

    def get_mass(self):
        return self.mass

    def get_momentum(self):
        return self.mass * self.velocity

    def get_height(self):
        return self.height

    def get_width(self):
        return self.width

    def get_left(self):
        return self.left

    def get_top(self):
        return self.top

    def get_velocity(self):
        return self.velocity

    def get_velocity_x(self):
        return self.vx

    def get_velocity_y(self):
        return self.vy

    def is_empty(self):
        return (self.mass <= 0)

    def contains(self, x: int, y: int):
        """If contains

        Args:
            x (int): X
            y (int): Y

        Returns:
            bool: If it contains
        """
        cleanX = self.xMod(x)
        cleanY = self.yMod(y)

        return cleanX >= self.left and cleanX < (self.left+self.width) and cleanY >= self.top and cleanY < (self.top+self.height)

    def calculate_crust(self, x: int, y: int, index: int, w_crust: float, e_crust: float, n_crust: float, s_crust: float, w: int, e: int, n: int, s: int):
        """Calculates the crust

        Args:
            x (int): [description]
            y (int): [description]
            index (int): [description]

            w_crust (float): [description]
            e_crust (float): [description]
            n_crust (float): [description]
            s_crust (float): [description]

            w (int): [description]
            e (int): [description]
            n (int): [description]
            s (int): [description]
        """

        # Build masks for accessible directions (4-way).
        # Allow wrapping around map edges if plate has world wide dimensions.
        w_mask = -((x > 0) | (self.width == self.wd.get_width()))
        e_mask = -((x < self.width - 1) | (self.width == self.wd.get_width()))
        n_mask = -((y > 0) | (self.height == self.wd.get_height()))
        s_mask = -((y < self.height - 1) |
                   (self.height == self.wd.get_height()))

        # Calculate the x and y offset of neighbour directions.
        # If neighbour is out of plate edges, set it to zero. This protects
        # map memory reads from segment faulting.
        w = self.xMod(x-1) if w_mask == -1 else 0
        e = self.xMod(x+1) if e_mask == -1 else 0
        n = self.yMod(y-1) if n_mask == -1 else 0
        s = self.yMod(y+1) if s_mask == -1 else 0

        # Calculate offsets within map memory.
        # w = y * self.width + w
        # e = y * self.width + e
        # n = n * self.width + x
        # s = s * self.width + x

        # Extract neighbours heights. Apply validity filtering: 0 is invalid.
        # Might fail due to last portion being commented out
        w_crust = self.map._data[y, w] * \
            (w_mask & (self.map._data[y, w] < self.map._data[y, x]))
        e_crust = self.map._data[y, e] * \
            (e_mask & (self.map._data[y, e] < self.map._data[y, x]))
        n_crust = self.map._data[n, x] * \
            (n_mask & (self.map._data[n, x] < self.map._data[y, x]))
        s_crust = self.map._data[s, x] * \
            (s_mask & (self.map._data[s, x] < self.map._data[y, x]))

    def xMod(self, x: int):
        """Gets the modulus of the X axis

        Args:
            x (int): X value

        Returns:
            int: Remainder
        """
        return (x + self.wd.get_width()) % self.wd.get_width()

    def yMod(self, y: int):
        """Gets the modulus of the Y axis

        Args:
            y (int): Y value

        Returns:
            int: Remainder
        """
        return (y + self.wd.get_height()) % self.wd.get_height()

    class SegmentData(Rectangle):
        def __init__(self, wd: WorldDimension, left: int, right: int, top: int, bottom: int, area: int):
            # Might need to add the Rectangle separately
            self.area = area
            self.coll_count = 0
            # See if theres a better way of doing this.
            super().__init__(wd, left, right, top, bottom)

        def is_empty(self):
            """Returns a boolean of whether the Segment is empty or not.

            Returns:
                bool:
            """
            return (self.area == 0)

    def calc_direction(self, x: int, y: int, origin_index: int, ID: int):

        if self.segment[origin_index] < ID:
            return self.segment[origin_index]

        can_go_left = x > 0
        if can_go_left:
            can_go_left = self.map._data[y, x - 1] >= CONT_BASE
        can_go_right = x < self.width - 1
        if can_go_right:
            can_go_right = self.map._data[y, x + 1] >= CONT_BASE
        can_go_up = y > 0
        if can_go_up:
            can_go_up = self.map._data[y - 1, x] >= CONT_BASE
        can_go_down = y < self.height - 1
        if can_go_down:
            can_go_down = self.map._data[y + 1, x] >= CONT_BASE
        nbour_id = ID

        if can_go_left and self.segment[origin_index - 1] < ID:
            nbour_id = self.segment[origin_index - 1]
        elif can_go_right and self.segment[origin_index + 1] < ID:
            nbour_id = self.segment[origin_index + 1]
        elif can_go_up and self.segment[origin_index - self.width] < ID:
            nbour_id = self.segment[origin_index - self.width]
        elif can_go_down and self.segment[origin_index + self.width] < ID:
            nbour_id = self.segment[origin_index + self.width]

        return nbour_id

    def scan_spans(self, line, start, end, spans_todo, spans_done):
        while True:
            end = spans_todo[line].pop()
            start = spans_todo[line].pop()

            j = 0
            for temp in range(len(spans_done[line])):

                if (start >= spans_done[line][j] and start <= spans_done[line][j+1]):
                    start = spans_done[line][j+1] + 1

                if (end >= spans_done[line][j] and end <= spans_done[line][j+1]):
                    end = spans_done[line][j] - 1

                j += 2

                # Be careful with this because it was changed
                # The exit statement for the while loop
                if j >= len(spans_done[line]):
                    break

            # This is still confusing
            start = start | -int(end >= self.width)
            end = end - int(end >= self.width)

            # Be careful with this because it was changed
            # The exit statement for the while loop
            if not start > end or not len(spans_todo[line]):
                break

        return (start, end, spans_todo, spans_done)

    def create_segment(self, x: int, y: int):
        """NEEDS THOROUGH REVIEW
        Separate a continent at (X, Y) to its own partition.

        Method analyzes the pixels 4-ways adjacent at the given location
        and labels all connected continental points with same segment ID.

        Args:
            x (int): Offset on the local height map along X axis.
            y (int): Offset on the local height map along Y axis.

        Returns:
            int: ID of created segment on success, otherwise -1.
        """
        origin_index = y * self.width + x
        ID = len(self.seg_data)

        nbour_id = self.calc_direction(x, y, origin_index, ID)

        if nbour_id < ID:
            self.segment[origin_index] = nbour_id
            self.seg_data[nbour_id].area += 1
            self.seg_data[nbour_id].enlarge_to_contain(x, y)
            return nbour_id

        data = self.SegmentData(self.wd, x, x, y, y, 0)

        spans_todo = [[] for _ in range(self.height)]
        spans_done = [[] for _ in range(self.height)]

        self.segment[origin_index] = ID
        spans_todo[y].append(x)
        spans_todo[y].append(x)  # Why is this performed twice?

        lines_processed = -1
        while lines_processed != 0:
            lines_processed = 0
            start = 0
            end = 0
            for line in range(self.height):

                if len(spans_todo[line]) == 0:
                    continue

                start, end, spans_todo, spans_done = self.scan_spans(
                    line, start, end, spans_todo, spans_done)

                if start > end:
                    continue

                row_above = line - 1 if line > 0 else self.height - 1
                row_below = line + 1 if line < self.height - 1 else 0

                line_here = line * self.width
                line_above = row_above * self.width
                line_below = row_below * self.width

                # Extend the beginning of line
                while start > 0 \
                        and self.segment[line_here + start - 1] > ID \
                        and self.map._data[line, start - 1] >= CONT_BASE:
                    start -= 1
                    self.segment[line_here + start] = ID
                    # Count volume of pixel

                #  Extend the end of line.
                while end < self.width - 1 \
                        and self.segment[line_here + end + 1] > ID \
                        and self.map._data[line, end + 1] >= CONT_BASE:
                    end += 1
                    self.segment[line_here + end] = ID
                    # Count volume of pixel

                # Check if should wrap around left edge.
                if self.width == self.wd.get_width() \
                        and start == 0 \
                        and self.segment[line_here + self.width - 1] > ID \
                        and self.map._data[line, self.width - 1] >= CONT_BASE:
                    self.segment[line_here + self.width - 1] = ID
                    spans_todo[line].append(self.width - 1)
                    spans_todo[line].append(self.width - 1)
                    # Count volume of pixel

                # Check if should wrap around right edge.
                if self.width == self.wd.get_width() \
                        and end == self.width - 1 \
                        and self.segment[line_here + 0] > ID \
                        and self.map._data[line, 0] >= CONT_BASE:
                    self.segment[line_here + 0] = ID
                    spans_todo[line].append(0)
                    spans_todo[line].append(0)

                # Update segment area counter
                data.area += 1 + end - start

                # Record any changes in extreme dimensions.
                if line < data.get_top():
                    data.set_top(line)
                if line > data.get_bottom():
                    data.set_bottom(line)
                if start < data.get_left():
                    data.set_left(start)
                if end > data.get_right():
                    data.set_right(end)

                if line > 0 \
                        or self.height == self.wd.get_height():

                    j = int(start)
                    while j <= end:
                        if self.segment[line_above + j] > ID \
                                and self.map._data[row_above, j] >= CONT_BASE:
                            a = int(j)
                            self.segment[line_above + a] = ID
                            # Count volume of pixel

                            j += 1
                            while j < self.width \
                                    and self.segment[line_above + j] > ID \
                                    and self.map._data[row_above, j] >= CONT_BASE:
                                self.segment[line_above + j] = ID
                                j += 1
                                # Count volume of pixel

                            b = int(j - 1)
                            spans_todo[row_above].append(a)
                            spans_todo[row_above].append(b)
                        else:
                            j += 1

                if line < self.height - 1 \
                        or self.height == self.wd.get_height():

                    j = int(start)
                    while j <= end:
                        if self.segment[line_below + j] > ID \
                                and self.map._data[row_below, j] >= CONT_BASE:
                            a = int(j)
                            self.segment[line_below + a] = ID
                            # Count volume of pixel

                            j = j + 1
                            while j < self.width \
                                    and self.segment[line_below + j] > ID \
                                    and self.map._data[row_below, j] >= CONT_BASE:
                                self.segment[line_below + j] = ID
                                j += 1
                                # Count volume of pixel

                            # Last point is invalid.
                            b = int(j - 1)
                            spans_todo[row_below].append(a)
                            spans_todo[row_below].append(b)
                        else:
                            j += 1

                spans_done[line].append(start)
                spans_done[line].append(end)
                lines_processed += 1

        self.seg_data.append(data)

        return ID

    def get_map_index(self, px: int, py: int):
        """Translate world coordinates into offset within plate's height map.

        If the global world map coordinates are within plate's height map,
        the values of passed coordinates will be altered to contain the
        X and y offset within the plate's height map. Otherwise values are
        left intact.

        Args:
            px (int): Offset on the global world map along X axis.
            py (int): Offset on the global world map along Y axis.

        Returns:
            int: Offset in height map or -1 on error.
        """
        ilft = int(self.left)
        itop = int(self.top)
        irgt = ilft + self.width
        ibtm = itop + self.height

        rect = Rectangle(self.wd, ilft, irgt, itop, ibtm)
        return rect.get_map_index(px, py)

    def get_map_indeces(self, px: int, py: int):
        """Translate world coordinates into offset within plate's height map.

        If the global world map coordinates are within plate's height map,
        the values of passed coordinates will be altered to contain the
        X and y offset within the plate's height map. Otherwise values are
        left intact.

        Args:
            px (int): Offset on the global world map along X axis.
            py (int): Offset on the global world map along Y axis.

        Returns:
            int: Offset in height map or -1 on error.
        """
        ilft = int(self.left)
        itop = int(self.top)
        irgt = ilft + self.width
        ibtm = itop + self.height

        rect = Rectangle(self.wd, ilft, irgt, itop, ibtm)
        return rect.get_map_indeces(px, py)

    def get_continent_at(self, x: int, y: int):
        # Need to convert this function to a class method.
        index = int(self.get_map_index(x, y))
        lx, ly = self.get_map_indeces(x, y)
        seg = self.segment[index]

        if seg >= len(self.seg_data):
            seg = self.create_segment(lx, ly)

        if seg >= len(self.seg_data):
            assert('Failed to create segment')

        return seg

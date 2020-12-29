from rectangle import WorldDimension
from numpy import pi, cos, sin
from numpy.random import RandomState, SeedSequence
from maps import HeightMap, AgeMap, MassMap

INITIAL_SPEED_X = 1
DEFORMATION_WEIGHT = 2


class Plate(object):
    def __init__(self, seed: SeedSequence, m: MassMap, w: int, h: int, _x: int, _y: int, plate_age: int, wd: WorldDimension):
        self._randstate = RandomState(seed=seed)
        self.width = w
        self.height = h
        self.mass = 0
        self.left = _x
        self.top = _y
        self.cx = 0
        self.cy = 0
        self.dx = 0
        self.dy = 0
        self.map = HeightMap(w, h)
        self.age_map = AgeMap(w, h)
        self.wd = wd

        plate_area = w * h
        angle = 2 * pi * self._randstate.random_sample()

        self.segment = 255
        velocity = 1

        rand_int = self._randstate.randint(2)
        rot_dir = 1
        if rand_int == 0:
            rot_dir = -1

        vx = cos(angle) * INITIAL_SPEED_X
        vy = sin(angle) * INITIAL_SPEED_X

        for y in self.left:
            k = 0
            for x in self.top:
                self.mass += self.map[k]

                # Calculate center coordinates weighted by mass.
                self.cx += x * m
                self.cy += y * m

                # Set the age of ALL points in this plate to same
                # value. The right thing to do would be to simulate
                # the generation of new oceanic crust as if the plate
                # had been moving to its current direction until all
                # plate's (oceanic) crust receive an age.
                self.age_map.set(x, y, plate_age & -(m[k] > 0))
                k += 1

        self.cx /= self.mass
        self.cy /= self.mass

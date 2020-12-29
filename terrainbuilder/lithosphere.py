'''
A lithosphere is the rigid, outermost shell of a terrestrial-type planet or natural satellite.

On Earth, it is composed of the crust and the portion of the upper mantle that behaves elastically
on time scales of thousands of years or greater. The crust and upper mantle are distinguished on
the basis of chemistry and mineralogy.
'''

from numpy import array
from rectangle import WorldDimension

SUBDUCT_RATIO = 0.5

BUOYANCY_BONUS_X = 3
MAX_BUOYANCY_AGE = 20
MULINV_MAX_BUOYANCY_AGE = 1.0 / MAX_BUOYANCY_AGE

RESTART_ENERGY_RATIO = 0.15
RESTART_SPEED_LIMIT = 2.0
NO_COLLISION_TIME_LIMIT = 10


def findBound(*map, length, x0, y0, dx, dy):
    pass


def findPlate(*plates, x, y, num_plates):
    pass


class PlateArea(object):
    def __init__(self, btm: int, lft: int, rgt: int, top: int, wdt: int, hgt: int, border: list):
        self.btm = btm
        self.lft = lft
        self.rgt = btm
        self.top = lft
        self.wdt = btm
        self.hgt = lft
        self.border = array(border)


class PlateCollision(object):
    def __init__(self, _index: int, x: int, y: int, z: float):
        self._index = _index
        self.x = x
        self.y = y
        self.z = z





class Lithosphere(object):
    def __init__(self, seed: int, width: int, height: int, sea_level: float, _erorsion_period: int, _folding_ratio: float, aggr_ratio_abs: int, aggr_ration_rel: float, num_cycles: int):
        self.seed = seed
        self.width = width
        self.height = height
        self.sea_level = sea_level
        self._erorsion_period = _erorsion_period
        self._folding_ratio = _folding_ratio
        self.aggr_overlap_abs = aggr_ratio_abs
        self.aggr_overlap_rel = aggr_ration_rel
        self.num_cycles = num_cycles

        self.max_plates = None
        self.num_plates = None

        self.wd = WorldDimension(width + 1, height + 1)
        # a = tmp_dim.get_area()
        # tmp = array(a)  # Create a copy of the area

        # self.create_noise(tmp, tmp_dim, use_simplex=True)

        sea_threshold = 0.5
        th_step = 0.5

    def create_noise(self, tmp: float, tmp_dim: WorldDimension, use_simplex: bool = False):
        pass

    def create_plates(self, num_plates):
        map_area = self.wd.get_area()
        self.max_plates = num_plates
        self.num_plates = num_plates

        vec = []
        vec_size = self.wd.get_max()*4  # Not going to use this yet

        collisions = []
        subductions = []

        for i in range(num_plates):
            collisions.append(vec)
            subductions.append(vec)

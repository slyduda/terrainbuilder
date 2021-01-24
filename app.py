from plate.plateapi import platec_api_create, platec_api_step
from draw.plotter import plot_height
from numpy import int64, array
from numpy.random import randint
from numpy.random import SeedSequence, RandomState, MT19937
import time

from progress.bar import IncrementalBar

SEED = None


if not SEED:
    seed = SeedSequence(pool_size=6)
    SEED = seed.entropy

print(SEED)

seed = SeedSequence(SEED)
current_time = time.time()

X = 64
Y = 64
I = 24

p = platec_api_create(seed=seed, width=X, height=Y,
                      sea_level=0.25, erosion_period=60,
                      folding_ratio=0.02, aggr_overlap_abs=1000000,
                      aggr_overlap_rel=0.33, cycle_count=2, num_plates=6)
print(time.time() - current_time, "seconds")
current_time = time.time()
plot_height(p)

j_end = I
bar = IncrementalBar('Updating Plates', max=j_end)
for x in range(I):
    p = platec_api_step(p)
    bar.next()
bar.finish()

print(time.time() - current_time, "seconds")
current_time = time.time()
plot_height(p)

print(p.hmap._data.max(), p.hmap._data.min())

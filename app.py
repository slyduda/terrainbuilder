from plate.plateapi import platec_api_create
from numpy import int64
from numpy.random import randint
from numpy.random import SeedSequence, RandomState, MT19937

SEED = None

if not SEED:
    seed = SeedSequence(pool_size=6)
    SEED = seed.entropy

print(SEED)

seed = SeedSequence(SEED)

p = platec_api_create(seed=seed, width=32, height=32,
                      sea_level=0.65, erosion_period=60,
                      folding_ratio=0.02, aggr_overlap_abs=1000000,
                      aggr_overlap_rel=0.33, cycle_count=2, num_plates=10)

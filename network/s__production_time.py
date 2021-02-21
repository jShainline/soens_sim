import numpy as np

from _util import physical_constants

p = physical_constants()

#%%
num_light_sources = 1e6
time_per_light_source = 1 # seconds

print('total time to fabricate light sources = {:6.2f} days'.format(num_light_sources*time_per_light_source/(60*60)))
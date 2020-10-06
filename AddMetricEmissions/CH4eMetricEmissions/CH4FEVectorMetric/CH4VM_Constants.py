import numpy as np
import pandas as pd

# Test Emissions Timeseries
TIME_PERIOD = 300
LLCP_EMISSIONS = np.zeros(TIME_PERIOD).tolist()
LLCP_EMISSIONS[0] = 1
YEARS = np.arange(0, TIME_PERIOD).tolist()

LLCP_TIMESERIES = pd.DataFrame(
    data={'Year': YEARS,
          'SLCP': LLCP_EMISSIONS})

# CO2fe Calculator Constants

# CARBON
CARBON_BOXES = np.array([0.2173, 0.2240, 0.2824, 0.2763])
CARBON_LIFETIMES = np.array([1000000, 394.4, 36.54, 4.304])

# METHANE
METHANE_BOXES = np.array([1, 0, 0, 0])
METHANE_LIFETIMES = np.array([12.4, 1, 1, 1])

DECAY_CONSTANTS = pd.DataFrame(
    data={'CarbonBoxes': CARBON_BOXES,
          'CarbonLifetimes': CARBON_LIFETIMES,
          'MethaneBoxes': METHANE_BOXES,
          'MethaneLifetimes': METHANE_LIFETIMES})


# # CO2we Calculator Constants
# # DT = 20
# # H = 100
# # GWPH = 28
# # R = 0.75
# # S = 0.25


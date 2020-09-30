
from math import *
import mpmath as mp
mp.dps = 15; mp.pretty = True
import numpy as np
import pandas as pd
import AddMetricEmissions.CO2eMetricEmissions.CO2FEVectorMetric.CO2VM_Constants as Constants


# CO2fe and CO2we Emissions Calculator

def Dco2(k, t):
    Dco2 = Constants.CARBON_BOXES[k] * exp(-t / (Constants.CARBON_LIFETIMES[k]))
    return Dco2


def Dch4(k, t):
    Dch4 = Constants.METHANE_BOXES[k] * exp(-t / (Constants.METHANE_LIFETIMES[k]))
    return Dch4


def CalculateCO2ForcingEquivalent(slcpTimeSeries):
    CO2feTimeSeries = slcpTimeSeries.copy()

    # Create the output DataFrame and add Columns for CO2fe(t) and CO2we(t) emissions and R(t) and S(t) intermediaries.
    CO2feTimeSeries = slcpTimeSeries
    BlankEmissions = np.zeros(len(slcpTimeSeries))
    FourColumnBlankEmissions = pd.DataFrame(
        np.array([BlankEmissions, BlankEmissions, BlankEmissions, BlankEmissions]).transpose())

    #   Make Data-frames for R and S
    R = pd.concat([slcpTimeSeries[{'Year', 'SLCP'}], pd.DataFrame(FourColumnBlankEmissions)], axis=1)
    R['Total'] = BlankEmissions
    S = pd.concat([slcpTimeSeries[{'Year', 'SLCP'}], pd.DataFrame(FourColumnBlankEmissions)], axis=1)
    S['Total'] = BlankEmissions

    #   Compute the CO2fe values for each year
    for i in CO2feTimeSeries.index:

        for k in Constants.DECAY_CONSTANTS.index:
            #           Compute S(i) and enter that value for the current i (Year)
            S.loc[i, k] = mp.nsum(lambda j: CO2feTimeSeries.loc[int(j), 'SLCP'] * Dch4(k, i - j), [0, i])

            #           Compute R(i) given S(i) and existing values of R(i)
            R.loc[i, k] = mp.nsum(lambda j: (S.loc[int(j), 'Total'] - R.loc[int(j), 'Total']) * Dco2(k, i - j),
                                  [0, i - 1])

        S.loc[i, 'Total'] = mp.nsum(lambda k: S.loc[i, k], [0, len(Constants.DECAY_CONSTANTS.index) - 1])
        R.loc[i, 'Total'] = mp.nsum(lambda k: R.loc[i, k], [0, len(Constants.DECAY_CONSTANTS.index) - 1])

        #       Compute CO2fe(i) from these values of S(i) and R(i)
        CO2feTimeSeries.loc[i, 'CO2fe'] = (S.loc[i, 'Total'] - R.loc[i, 'Total'])

        # #       Compute CO2we
        # Et = CO2feTimeSeries.loc[i, 'SLCP']
        # EtMinusdt = 0
        # if i > 19:
        #     EtMinusdt = CO2feTimeSeries.loc[int(i - 20), 'SLCP']
        # CO2feTimeSeries.loc[i, 'CO2we'] = GWPH * ((r * H) / dt * (Et - EtMinusdt) + s * Et)

    return CO2feTimeSeries
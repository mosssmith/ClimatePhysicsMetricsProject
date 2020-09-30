from math import *
import numpy as np
import pandas as pd
import AddMetricEmissions.CO2eMetricEmissions.CO2FEVectorMetric.CO2VM_Constants as Constants
import mpmath as mp
mp.dps = 15
mp.pretty = True


# CO2fe and CO2we Emissions Calculator

def Dco2(k, t):
    return Constants.CARBON_BOXES[k] * exp(-t / (Constants.CARBON_LIFETIMES[k]))


def Dch4(k, t):
    return Constants.METHANE_BOXES[k] * exp(-t / (Constants.METHANE_LIFETIMES[k]))


def CalculateCH4ForcingEquivalent(llcp_time_series):
    ch4fe_time_series = llcp_time_series.copy()

    # Create the output DataFrame and add Columns for CO2fe(t) and CO2we(t) emissions and R(t) and S(t) intermediaries.
    ch4fe_time_series = llcp_time_series
    BlankEmissions = np.zeros(len(llcp_time_series))
    FourColumnBlankEmissions = pd.DataFrame(
        np.array([BlankEmissions, BlankEmissions, BlankEmissions, BlankEmissions]).transpose())

    #   Make Data-frames for R and S
    R = pd.concat([llcp_time_series[{'Year', 'SLCP'}], pd.DataFrame(FourColumnBlankEmissions)], axis=1)
    R['Total'] = BlankEmissions
    S = pd.concat([llcp_time_series[{'Year', 'SLCP'}], pd.DataFrame(FourColumnBlankEmissions)], axis=1)
    S['Total'] = BlankEmissions

    #   Compute the CH4fe values for each year
    for i in ch4fe_time_series.index:

        for k in Constants.DECAY_CONSTANTS.index:
            #           Compute S(i) and enter that value for the current i (Year)
            S.loc[i, k] = mp.nsum(lambda j: ch4fe_time_series.loc[int(j), 'SLCP'] * Dco2(k, i - j), [0, i])

            #           Compute R(i) given S(i) and existing values of R(i)
            R.loc[i, k] = mp.nsum(lambda j: (S.loc[int(j), 'Total'] - R.loc[int(j), 'Total']) * Dch4(k, i - j),
                                  [0, i - 1])

        S.loc[i, 'Total'] = mp.nsum(lambda k: S.loc[i, k], [0, len(Constants.DECAY_CONSTANTS.index) - 1])
        R.loc[i, 'Total'] = mp.nsum(lambda k: R.loc[i, k], [0, len(Constants.DECAY_CONSTANTS.index) - 1])

        #       Compute CO2fe(i) from these values of S(i) and R(i)
        ch4fe_time_series.loc[i, 'CH4fe'] = (S.loc[i, 'Total'] - R.loc[i, 'Total'])

    return ch4fe_time_series

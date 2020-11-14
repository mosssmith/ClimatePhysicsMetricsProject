import AddMetricEmissions.CO2eMetricEmissions.CO2FEVectorMetric.CO2VM_Constants as Constants
import AddMetricEmissions.CO2eMetricEmissions.CO2FEVectorMetric.CO2VM_Calculator as Calculator
import mpmath as mp
import numpy as np
mp.dps = 15
mp.pretty = True

CO2feVectorMetric = Calculator.CalculateCO2ForcingEquivalent(Constants.SLCP_TIMESERIES)

print(CO2feVectorMetric.head())

CO2feVectorMetric.to_csv('CO2feVectorMetricFactors.csv')

# Make the lower-diagonal tirpitz matrix for conversion between gases.
FullAGFP_CH4 = np.zeros((len(CO2feVectorMetric.index), len(CO2feVectorMetric.index)))

for i in np.arange(len(CO2feVectorMetric.index)):
    for j in np.arange(i+1):
        FullAGFP_CH4[i][j] = CO2feVectorMetric['CO2fe'].loc[i-j]


# CO2feVectorMetric.to_csv('CO2feVectorMetricMatrix.csv')
np.savetxt("CO2feVectorMetricMatrix.csv", FullAGFP_CH4, delimiter=",")

print(FullAGFP_CH4.size)

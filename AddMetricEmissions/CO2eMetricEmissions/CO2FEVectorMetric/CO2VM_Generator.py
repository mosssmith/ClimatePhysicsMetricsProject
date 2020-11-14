import AddMetricEmissions.CO2eMetricEmissions.CO2FEVectorMetric.CO2VM_Constants as Constants
import AddMetricEmissions.CO2eMetricEmissions.CO2FEVectorMetric.CO2VM_Calculator as Calculator
import mpmath as mp

mp.dps = 15
mp.pretty = True

CO2feVectorMetric = Calculator.CalculateCO2ForcingEquivalent(Constants.SLCP_TIMESERIES)

print(CO2feVectorMetric.head())

CO2feVectorMetric.to_csv('CO2feVectorMetricFactors.csv')

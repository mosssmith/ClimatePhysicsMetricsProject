import AddMetricEmissions.CH4eMetricEmissions.CH4FEVectorMetric.CH4VM_Constants as Constants
import AddMetricEmissions.CH4eMetricEmissions.CH4FEVectorMetric.CH4VM_Calculator as Calculator
import mpmath as mp
mp.dps = 15
mp.pretty = True

CH4feVectorMetric = Calculator.CalculateCH4ForcingEquivalent(Constants.LLCP_TIMESERIES)

print(CH4feVectorMetric.head())

CH4feVectorMetric.to_csv('CH4FEVectorMetricFactors.csv')

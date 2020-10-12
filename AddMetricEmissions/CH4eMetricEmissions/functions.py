import mpmath as mp
mp.dps = 15; mp.pretty = True
import pandas as pd
import numpy as np
from pathlib import Path

# Add VectorMetric into the functions directly
VECTOR_METRIC = pd.read_csv(Path(__file__).parent / './CH4FEVectorMetric/CH4FEVectorMetricFactors.csv')

METRIC_CONSTANTS = pd.DataFrame({
    'GWP100': [28],
    'GWP20': [84],
    'GTP100': [4],
    'GTP20': [67],
    'H': [100],
    'CGWP100': [4300],
    'CGTP75': [3700],
    'r': [0.75],
    's': [0.25],
    'dt': [20],
    'REch4': [0.000599],       #W/m2ppb AR5 Chapt.8 Appendix 8.A. 1.65 * 0.000363
    'REco2': [0.0000137],      #W/m2ppb AR5 Chapt.8 Appendix 8.A.
    'Convch4': [0.351828],  #ppb/MtCH4 GIR
    'Convco2': [0.1282496]   #ppb/MtCO2 [0.46895] #ppb/MtC GIR - where does this come from?
})


def addCH4eMetricEmissions(gir_emissions_series, methods):
    Scenarios = gir_emissions_series.columns.levels[0]

    gir_emissions_series_output = gir_emissions_series.copy()

    for scenario in Scenarios:
        scen_names = []
        gases_in = ['CO2', 'CH4', 'N2O']

        #       MAKE DATAFRAME LARGE ENOUGH TO HOLD METRIC EMISSIONS SCENARIOS
        for method in methods:
            # Add new scenario for each method with new scenarios named "Scenario" + " - " + "Method"
            ColumnName = scenario + " - " + method
            scen_names.append(ColumnName)

        # Make a correctly sized, shaped, and labelled dataframe to add to the emissions series.
        dfToAdd = pd.DataFrame(np.zeros((len(gir_emissions_series.index), len(gases_in) * len(methods))),
                               columns=pd.MultiIndex.from_product([scen_names, gases_in]),
                               index=gir_emissions_series_output.index)
        gir_emissions_series_output = gir_emissions_series_output.join(dfToAdd)

        #       MAKE LLCPEmissions DATAFRAME USING ORIGINAL SCENARIO EMISSIONS
        llcp_time_series = pd.DataFrame(data={'Year': gir_emissions_series[scenario, 'CO2'].index,
                                              'LLCP Emissions': gir_emissions_series[scenario, 'CO2'].to_list()})
        llcp_metric_time_series = addMetricEmissions(llcp_time_series)
        for method in methods:
            ColumnName = scenario + " - " + method
            gir_emissions_series_output[ColumnName, 'CH4'] = llcp_metric_time_series[['Year', method]].set_index('Year')
    return gir_emissions_series_output


def addMetricEmissions(llcp_emissions_series):
    # #   Key to Metric Constants
    llcp_emissions_series = addCH4feEmsColumn(VECTOR_METRIC, llcp_emissions_series)
    llcp_emissions_series = addGWPStarEmsColumn(llcp_emissions_series)
    llcp_emissions_series = addCGWPEmsColumn(llcp_emissions_series)
    llcp_emissions_series = addCGTPEmsColumn(llcp_emissions_series)
    llcp_emissions_series = addGWPEmsColumn(llcp_emissions_series)

    return llcp_emissions_series


def addGWPEmsColumn(llcp_emissions_series):
    llcp_emissions_series["GWP100"] = np.zeros(len(llcp_emissions_series.index)).tolist()
    llcp_emissions_series["GWP20"] = np.zeros(len(llcp_emissions_series.index)).tolist()
    llcp_emissions_series["GTP100"] = np.zeros(len(llcp_emissions_series.index)).tolist()
    llcp_emissions_series["GTP20"] = np.zeros(len(llcp_emissions_series.index)).tolist()

    GWP100 = METRIC_CONSTANTS['GWP100'][0]
    GWP20 = METRIC_CONSTANTS['GWP20'][0]
    GTP100 = METRIC_CONSTANTS['GTP100'][0]
    GTP20 = METRIC_CONSTANTS['GTP20'][0]

    for i in llcp_emissions_series.index:
        Et = llcp_emissions_series["LLCP Emissions"].loc[i]
        llcp_emissions_series["GWP100"].loc[i] = Et / GWP100
        llcp_emissions_series["GWP20"].loc[i] = Et / GWP20
        llcp_emissions_series["GTP100"].loc[i] = Et / GTP100
        llcp_emissions_series["GTP20"].loc[i] = Et / GTP20

    return llcp_emissions_series


def addCGWPEmsColumn(llcp_emissions_series):
    llcp_emissions_series["CGWP"] = np.zeros(len(llcp_emissions_series.index)).tolist()

    CGWP100 = METRIC_CONSTANTS['CGWP100'][0]

    for i in llcp_emissions_series.index:
        Et = llcp_emissions_series["LLCP Emissions"].loc[i]

        if i == 0:
            EtCH4dt = 0
            llcp_emissions_series["CGWP"].loc[i] = Et / CGWP100 + EtCH4dt
        else:
            EtCH4dt = llcp_emissions_series["CGWP"].loc[i - 1]
            llcp_emissions_series["CGWP"].loc[i] = Et / CGWP100 + EtCH4dt

    return llcp_emissions_series


def addCGTPEmsColumn(llcp_emissions_series):
    llcp_emissions_series["CGTP"] = np.zeros(len(llcp_emissions_series.index)).tolist()

    CGTP75 = METRIC_CONSTANTS['CGTP75'][0]

    for i in llcp_emissions_series.index:
        Et = llcp_emissions_series["LLCP Emissions"].loc[i]

        if i == 0:
            EtCH4dt = 0
            llcp_emissions_series["CGTP"].loc[i] = Et / CGTP75 + EtCH4dt
        else:
            EtCH4dt = llcp_emissions_series["CGTP"].loc[i - 1]
            llcp_emissions_series["CGTP"].loc[i] = Et / CGTP75 + EtCH4dt

    return llcp_emissions_series


def addGWPStarEmsColumn(llcp_emissions_series):
    llcp_emissions_series["GWP*"] = np.zeros(len(llcp_emissions_series.index)).tolist()

    # CO2we constants:
    r = METRIC_CONSTANTS['r'][0]
    s = METRIC_CONSTANTS['s'][0]
    H = METRIC_CONSTANTS['H'][0]
    GWP100 = METRIC_CONSTANTS['GWP100'][0]
    dt = METRIC_CONSTANTS['dt'][0]

    for i in llcp_emissions_series.index:
        #       Compute CH4we
        Et = llcp_emissions_series["LLCP Emissions"].loc[i]
        EtCH4weMinusdt = 0
        if i > 19:
            EtCH4weMinusdt = llcp_emissions_series["GWP*"].loc[int(i - 20)]

        llcp_emissions_series["GWP*"].loc[i] = 1 / (s + (r * H) / dt) * ((r * H) / dt * EtCH4weMinusdt + Et / GWP100)

    return llcp_emissions_series


# Add the "CO2feEms" column, Gamma is the ratio of radiative efficiencies
def addCH4feEmsColumn(vector_metric, llcp_emissions_series):
    REch4 = METRIC_CONSTANTS['REch4'][0]
    REco2 = METRIC_CONSTANTS['REco2'][0]
    Convch4 = METRIC_CONSTANTS['Convch4'][0]
    Convco2 = METRIC_CONSTANTS['Convco2'][0]

    Gamma = (Convco2 * REco2) / (Convch4 * REch4)

    llcp_emissions_series = addBValueColumn(vector_metric, llcp_emissions_series)

    llcp_emissions_series["CFE"] = Gamma * (llcp_emissions_series["LLCP Emissions"] + llcp_emissions_series["BValue"])

    return llcp_emissions_series.drop(['NormalisedYear', 'MultipliedEmissions', 'BValue'], axis=1)


# Add B value column and insert the b value for each year
def addBValueColumn(vector_metric, llcp_emissions_series):
    llcp_emissions_series["BValue"] = np.zeros(len(llcp_emissions_series.index)).tolist()

    for i in llcp_emissions_series.index:
        llcp_emissions_series["BValue"].loc[i] = generateBValue(vector_metric, llcp_emissions_series,
                                                                llcp_emissions_series["Year"].loc[i])

    return llcp_emissions_series


# Generate B value
def generateBValue(vector_metric, llcp_emissions_series, year_specified):
    #   llcpEmissionsSeries takes the form of a dataframe with columns "Year" and "LLCP Emissions"

    #   Create a "NormalisedYear" column
    llcp_emissions_series["NormalisedYear"] = year_specified - llcp_emissions_series["Year"]
    #   Create a new column for the MultipliedEmissions
    llcp_emissions_series["MultipliedEmissions"] = np.zeros(len(llcp_emissions_series.index)).tolist()

    #   Match indexes
    llcp_emissions_series = llcp_emissions_series.set_index('NormalisedYear')
    vector_metric = vector_metric.set_index('Year')

    #   Calculate the MultipliedEmissions
    llcp_emissions_series["MultipliedEmissions"] = llcp_emissions_series["LLCP Emissions"].mul(vector_metric["CH4fe"],
                                                                                               fill_value=0)

    #   Sum the Values in the "MultipliedEmissions" column to get to the CO2-fe emission for that year.
    return sum(llcp_emissions_series["MultipliedEmissions"]) - llcp_emissions_series["MultipliedEmissions"].loc[0]

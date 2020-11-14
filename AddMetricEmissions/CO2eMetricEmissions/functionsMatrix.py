import mpmath as mp
mp.dps = 15; mp.pretty = True
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm.autonotebook import tqdm

# Add VectorMetric into the functions directly
VECTOR_METRIC = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/CO2FEVectorMetricFactors.csv')

VECTOR_MATRIX = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/CO2FEVectorMetricMatrix.csv', header=None).to_numpy()

METRIC_CONSTANTS = pd.DataFrame({
    'GWP100': [28],
    'GWP20': [84],
    'GTP100': [4],
    'GTP20': [67],
    'GWP100N2O': [265],
    'GWP20N2O': [264],
    'GTP100N2O': [234],
    'GTP20N2O': [277],
    'H': [100],
    'CGWP100': [4300],
    'CGTP75': [3700],
    'r': [0.75],
    's': [0.25],
    'dt': [20],
    'REch4': [0.000599],        #W/m2ppb AR5 Chapt.8 Appendix 8.A. 1.65 * 0.000363
    'REco2': [0.0000137],       #W/m2ppb AR5 Chapt.8 Appendix 8.A.
    'Convch4': [0.351828],      #ppb/MtCH4 GIR
    'Convco2': [0.1282496]      #ppb/MtCO2 [0.46895] #ppb/MtC GIR - where does this come from?
})

def addCO2eMetricEmissions(gir_emissions_series, methods):
    Scenarios = gir_emissions_series.columns.levels[0]

    gir_emissions_series_output = gir_emissions_series.copy()

    for scenario in tqdm(Scenarios, desc="Scenarios Completed"):
        scen_names = []
        gases_in = ['CO2', 'CH4', 'N2O']

        # MAKE DATAFRAME LARGE ENOUGH TO HOLD METRIC EMISSIONS SCENARIOS
        for method in methods:
            # Add new scenario for each method with new scenarios named "Scenario" + " - " + "Method"
            ColumnName = scenario + " - " + method
            scen_names.append(ColumnName)

        # Make a correctly sized, shaped, and labelled dataframe to add to the emissions series.
        dfToAdd = pd.DataFrame(np.zeros((len(gir_emissions_series.index), len(gases_in) * len(methods))),
                               columns=pd.MultiIndex.from_product([scen_names, gases_in]),
                               index=gir_emissions_series_output.index)
        gir_emissions_series_output = gir_emissions_series_output.join(dfToAdd)

        # MAKE SLCPEmissions DATAFRAME USING ORIGINAL SCENARIO EMISSIONS
        if 'CH4' in gir_emissions_series.columns.levels[1]:
            SLCPTimeSeries = pd.DataFrame(data={'Year': gir_emissions_series[scenario, 'CH4'].index,
                                                'SLCP Emissions': gir_emissions_series[scenario, 'CH4'].to_list()})
            slcp_metric_time_series = addMetricEmissions(SLCPTimeSeries)

        for method in methods:
            ColumnName = scenario + " - " + method
            # Add CO2 equivalent emissions for CH4
            if 'CH4' in gir_emissions_series.columns.levels[1]:
                gir_emissions_series_output[ColumnName, 'CO2'] = slcp_metric_time_series[['Year', method]].set_index('Year')

            # Add CO2 equivalent emissions for CO2
            if 'CO2' in gir_emissions_series.columns.levels[1]:
                gir_emissions_series_output[ColumnName, 'CO2'] += gir_emissions_series[scenario, 'CO2']

            # # Add CO2 equivalent emissions for N20
            if 'N2O' in gir_emissions_series.columns.levels[1]:
                if method is 'GWP20':
                    CO2eValueN20 = METRIC_CONSTANTS['GWP20N2O'][0]
                elif method is 'GTP100':
                    CO2eValueN20 = METRIC_CONSTANTS['GTP100N2O'][0]
                elif method is 'GTP20':
                    CO2eValueN20 = METRIC_CONSTANTS['GTP20N2O'][0]
                else:
                    CO2eValueN20 = METRIC_CONSTANTS['GWP100N2O'][0]
                gir_emissions_series_output[ColumnName, 'CO2'] += CO2eValueN20*gir_emissions_series[scenario, 'N2O']
        pass
    return gir_emissions_series_output


def addMetricEmissions(slcp_emissions_series):
    # #   Key to Metric Constants
    slcp_emissions_series = addCFEEmsColumn(VECTOR_MATRIX, slcp_emissions_series)
    # slcp_emissions_series = addCO2feEmsColumn(VECTOR_METRIC, slcp_emissions_series)
    slcp_emissions_series = addGWPStarEmsColumn(slcp_emissions_series)
    slcp_emissions_series = addCGWPEmsColumn(slcp_emissions_series)
    slcp_emissions_series = addCGTPEmsColumn(slcp_emissions_series)
    slcp_emissions_series = addGWPEmsColumn(slcp_emissions_series)

    #     Remove unnecessary columns

    return slcp_emissions_series


def addGWPEmsColumn(slcp_emissions_series):
    slcp_emissions_series["GWP100"] = np.zeros(len(slcp_emissions_series.index)).tolist()
    slcp_emissions_series["GWP20"] = np.zeros(len(slcp_emissions_series.index)).tolist()
    slcp_emissions_series["GTP100"] = np.zeros(len(slcp_emissions_series.index)).tolist()
    slcp_emissions_series["GTP20"] = np.zeros(len(slcp_emissions_series.index)).tolist()

    GWP100 = METRIC_CONSTANTS['GWP100'][0]
    GWP20 = METRIC_CONSTANTS['GWP20'][0]
    GTP100 = METRIC_CONSTANTS['GTP100'][0]
    GTP20 = METRIC_CONSTANTS['GTP20'][0]

    for i in slcp_emissions_series.index:
        Et = slcp_emissions_series["SLCP Emissions"].loc[i]
        slcp_emissions_series["GWP100"].loc[i] = Et * GWP100
        slcp_emissions_series["GWP20"].loc[i] = Et * GWP20
        slcp_emissions_series["GTP100"].loc[i] = Et * GTP100
        slcp_emissions_series["GTP20"].loc[i] = Et * GTP20

    return slcp_emissions_series


def addCGWPEmsColumn(slcp_emissions_series):
    slcp_emissions_series["CGWP"] = np.zeros(len(slcp_emissions_series.index)).tolist()

    CGWP100 = METRIC_CONSTANTS['CGWP100'][0]

    for i in slcp_emissions_series.index:
        Et = slcp_emissions_series["SLCP Emissions"].loc[i]

        if i == 0:
            Etdt = 0
            slcp_emissions_series["CGWP"].loc[i] = (CGWP100 * (Et - Etdt))
        else:
            Etdt = slcp_emissions_series["SLCP Emissions"].loc[i - 1]
            slcp_emissions_series["CGWP"].loc[i] = (CGWP100 * (Et - Etdt))

    return slcp_emissions_series


def addCGTPEmsColumn(slcp_emissions_series):
    slcp_emissions_series["CGTP"] = np.zeros(len(slcp_emissions_series.index)).tolist()

    CGTP75 = METRIC_CONSTANTS['CGTP75'][0]

    for i in slcp_emissions_series.index:
        Et = slcp_emissions_series["SLCP Emissions"].loc[i]

        if i == 0:
            Etdt = 0
            slcp_emissions_series["CGTP"].loc[i] = (CGTP75 * (Et - Etdt))
        else:
            Etdt = slcp_emissions_series["SLCP Emissions"].loc[i - 1]
            slcp_emissions_series["CGTP"].loc[i] = (CGTP75 * (Et - Etdt))

    return slcp_emissions_series


def addGWPStarEmsColumn(slcp_emissions_series):
    slcp_emissions_series["GWP*"] = np.zeros(len(slcp_emissions_series.index)).tolist()

    # CO2we constants:
    r = METRIC_CONSTANTS['r'][0]
    s = METRIC_CONSTANTS['s'][0]
    H = METRIC_CONSTANTS['H'][0]
    GWP100 = METRIC_CONSTANTS['GWP100'][0]
    dt = METRIC_CONSTANTS['dt'][0]

    for i in slcp_emissions_series.index:
        Et = slcp_emissions_series["SLCP Emissions"].loc[i]

        if i > 19:
            Etdt = slcp_emissions_series["SLCP Emissions"].loc[i - int(dt)]
            slcp_emissions_series["GWP*"].loc[i] = GWP100 * ((r * (Et - Etdt) * H) / dt + s * Et)
        elif i < 20:
            slcp_emissions_series["GWP*"].loc[i] = GWP100 * ((r * Et * H) / dt + s * Et)

    return slcp_emissions_series


# Add CFE emissions column using the matrix method

def addCFEEmsColumn(vector_matrix, slcp_emissions_series):

    # resize VECTOR_MATRIX to mate with timeseries.
    AGFP_CH4 = vector_matrix.copy()

    pos1 = len(slcp_emissions_series.index)
    pos2 = AGFP_CH4.shape[0]

    AGFP_CH4 = np.delete(AGFP_CH4, np.arange(pos1, pos2), 0)
    AGFP_CH4 = np.delete(AGFP_CH4, np.arange(pos1, pos2), 1)

    # Calculate Gamma Value
    REch4 = METRIC_CONSTANTS['REch4'][0]
    REco2 = METRIC_CONSTANTS['REco2'][0]
    Convch4 = METRIC_CONSTANTS['Convch4'][0]
    Convco2 = METRIC_CONSTANTS['Convco2'][0]
    Gamma = (Convch4 * REch4) / (Convco2 * REco2)

    # Add CFE column with result of matrix multiple to a new column.
    slcp_emissions_series_to_multiply = slcp_emissions_series[['SLCP Emissions']].to_numpy()

    slcp_emissions_series['CFE'] = Gamma * np.dot(AGFP_CH4, slcp_emissions_series_to_multiply)

    return slcp_emissions_series

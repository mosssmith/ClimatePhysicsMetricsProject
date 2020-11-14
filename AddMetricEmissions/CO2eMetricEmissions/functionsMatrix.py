import mpmath as mp
mp.dps = 15; mp.pretty = True
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm.autonotebook import tqdm


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
            slcp_metric_time_series = addCH4MetricEmsColumns(SLCPTimeSeries)

        # MAKE N2OEmissions DATAFRAME USING ORIGINAL SCENARIO EMISSIONS
        if 'N2O' in gir_emissions_series.columns.levels[1]:
            N2OTimeSeries = pd.DataFrame(data={'Year': gir_emissions_series[scenario, 'N2O'].index,
                                               'N2O Emissions': gir_emissions_series[scenario, 'N2O'].to_list()})
            N2O_metric_time_series = addN2OMetricEmsColumns(N2OTimeSeries)

        for method in methods:
            ColumnName = scenario + " - " + method
            # Add CO2 equivalent emissions for CH4
            if 'CH4' in gir_emissions_series.columns.levels[1]:
                gir_emissions_series_output[ColumnName, 'CO2'] = slcp_metric_time_series[['Year', method]].set_index('Year')

            # Add CO2 equivalent emissions for CO2
            if 'CO2' in gir_emissions_series.columns.levels[1]:
                gir_emissions_series_output[ColumnName, 'CO2'] += gir_emissions_series[scenario, 'CO2']

            # Add CO2 equivalent emissions for N20
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


# Add and populate columns for each metric for N2O
def addN2OMetricEmsColumns(slcp_emissions_series):
    GWP_100_MATRIX = resizeMatrix(GWP_100_MATRIX_N2O, slcp_emissions_series)
    GWP_20_MATRIX = resizeMatrix(GWP_20_MATRIX_N2O, slcp_emissions_series)
    GTP_100_MATRIX = resizeMatrix(GTP_100_MATRIX_N2O, slcp_emissions_series)
    GTP_20_MATRIX = resizeMatrix(GTP_20_MATRIX_N2O, slcp_emissions_series)

    # Add CFE column with result of matrix multiple to a new column.
    slcp_emissions_series_to_multiply = slcp_emissions_series[['N2O Emissions']].to_numpy()

    slcp_emissions_series['GWP100'] = np.dot(GWP_100_MATRIX, slcp_emissions_series_to_multiply)
    slcp_emissions_series['GWP20'] = np.dot(GWP_20_MATRIX, slcp_emissions_series_to_multiply)
    slcp_emissions_series['GTP100'] = np.dot(GTP_100_MATRIX, slcp_emissions_series_to_multiply)
    slcp_emissions_series['GTP20'] = np.dot(GTP_20_MATRIX, slcp_emissions_series_to_multiply)

    return slcp_emissions_series


# Add and populate columns for each metric for CH4
def addCH4MetricEmsColumns(slcp_emissions_series):
    # Calculate Gamma Value
    REch4 = METRIC_CONSTANTS['REch4'][0]
    REco2 = METRIC_CONSTANTS['REco2'][0]
    Convch4 = METRIC_CONSTANTS['Convch4'][0]
    Convco2 = METRIC_CONSTANTS['Convco2'][0]
    Gamma = (Convch4 * REch4) / (Convco2 * REco2)

    # Resize metric matrices
    CFE_MATRIX = resizeMatrix(CFE_MATRIX_CH4, slcp_emissions_series)
    GWP_STAR_MATRIX = resizeMatrix(GWP_STAR_MATRIX_CH4, slcp_emissions_series)
    CGWP_MATRIX = resizeMatrix(CGWP_MATRIX_CH4, slcp_emissions_series)
    CGTP_MATRIX = resizeMatrix(CGTP_MATRIX_CH4, slcp_emissions_series)
    GWP_100_MATRIX = resizeMatrix(GWP_100_MATRIX_CH4, slcp_emissions_series)
    GWP_20_MATRIX = resizeMatrix(GWP_20_MATRIX_CH4, slcp_emissions_series)
    GTP_100_MATRIX = resizeMatrix(GTP_100_MATRIX_CH4, slcp_emissions_series)
    GTP_20_MATRIX = resizeMatrix(GTP_20_MATRIX_CH4, slcp_emissions_series)

    # Add CFE column with result of matrix multiple to a new column.
    slcp_emissions_series_to_multiply = slcp_emissions_series[['SLCP Emissions']].to_numpy()

    slcp_emissions_series['CFE'] = Gamma * np.dot(CFE_MATRIX, slcp_emissions_series_to_multiply)
    slcp_emissions_series['GWP*'] = np.dot(GWP_STAR_MATRIX, slcp_emissions_series_to_multiply)
    slcp_emissions_series['CGWP'] = np.dot(CGWP_MATRIX, slcp_emissions_series_to_multiply)
    slcp_emissions_series['CGTP'] = np.dot(CGTP_MATRIX, slcp_emissions_series_to_multiply)
    slcp_emissions_series['GWP100'] = np.dot(GWP_100_MATRIX, slcp_emissions_series_to_multiply)
    slcp_emissions_series['GWP20'] = np.dot(GWP_20_MATRIX, slcp_emissions_series_to_multiply)
    slcp_emissions_series['GTP100'] = np.dot(GTP_100_MATRIX, slcp_emissions_series_to_multiply)
    slcp_emissions_series['GTP20'] = np.dot(GTP_20_MATRIX, slcp_emissions_series_to_multiply)

    return slcp_emissions_series


# Resize the matrix to fit with the timeseries
def resizeMatrix(input_matrix, slcp_emissions_series):
    output_matrix = input_matrix.copy()

    pos1 = len(slcp_emissions_series.index)
    pos2 = output_matrix.shape[0]

    output_matrix = np.delete(output_matrix, np.arange(pos1, pos2), 0)
    output_matrix = np.delete(output_matrix, np.arange(pos1, pos2), 1)
    return output_matrix


# Import CH4 Matrices
CFE_MATRIX_CH4 = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/CFE_MATRIX_CH4.csv', header=None).to_numpy()
GWP_STAR_MATRIX_CH4 = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/GWP_STAR_MATRIX_CH4.csv', header=None).to_numpy()
CGWP_MATRIX_CH4 = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/CGWP_MATRIX_CH4.csv', header=None).to_numpy()
CGTP_MATRIX_CH4 = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/CGTP_MATRIX_CH4.csv', header=None).to_numpy()
GWP_100_MATRIX_CH4 = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/GWP_100_MATRIX_CH4.csv', header=None).to_numpy()
GWP_20_MATRIX_CH4 = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/GWP_20_MATRIX_CH4.csv', header=None).to_numpy()
GTP_100_MATRIX_CH4 = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/GTP_100_MATRIX_CH4.csv', header=None).to_numpy()
GTP_20_MATRIX_CH4 = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/GTP_20_MATRIX_CH4.csv', header=None).to_numpy()

# Import N2O Matrices
GWP_100_MATRIX_N2O = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/GWP_100_MATRIX_N2O.csv', header=None).to_numpy()
GWP_20_MATRIX_N2O = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/GWP_20_MATRIX_N2O.csv', header=None).to_numpy()
GTP_100_MATRIX_N2O = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/GTP_100_MATRIX_N2O.csv', header=None).to_numpy()
GTP_20_MATRIX_N2O = pd.read_csv(Path(__file__).parent / './CO2FEVectorMetric/Matrices/GTP_20_MATRIX_N2O.csv', header=None).to_numpy()


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

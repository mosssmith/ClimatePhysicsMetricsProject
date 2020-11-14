import numpy as np
import pandas as pd

# Import constants used to specify matrix elements
CO2feVectorMetric = pd.read_csv('CO2FEVectorMetricFactors.csv')
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

MatrixScale = len(CO2feVectorMetric.index)

# Generate CFE_MATRIX (lower-diagonal tirpitz matrix)
CFE_MATRIX_CH4 = np.zeros((MatrixScale, MatrixScale))
for i in np.arange(MatrixScale):
    for j in np.arange(i+1):
        CFE_MATRIX_CH4[i][j] = CO2feVectorMetric['CO2fe'].loc[i - j]

# Generate GWP_STAR_MATRIX
PosValue = METRIC_CONSTANTS['GWP100']*(METRIC_CONSTANTS['H'] * METRIC_CONSTANTS['r'] / METRIC_CONSTANTS['dt'] + METRIC_CONSTANTS['s'])
NegValue = METRIC_CONSTANTS['GWP100']*(METRIC_CONSTANTS['H'] * METRIC_CONSTANTS['r'] / METRIC_CONSTANTS['dt'])
GWP_STAR_MATRIX_CH4 = int(PosValue) * np.identity(MatrixScale)
for i in np.arange(int(METRIC_CONSTANTS['dt']), MatrixScale):
    GWP_STAR_MATRIX_CH4[i][i-int(METRIC_CONSTANTS['dt'])] = - NegValue

# Generate CGWP_MATRIX
Value = METRIC_CONSTANTS['CGWP100']
CGWP_MATRIX_CH4 = int(Value) * np.identity(MatrixScale)
for i in np.arange(1, MatrixScale):
    CGWP_MATRIX_CH4[i][i-1] = -Value

# Generate CGTP_MATRIX
Value = METRIC_CONSTANTS['CGTP75']
CGTP_MATRIX_CH4 = int(Value) * np.identity(MatrixScale)
for i in np.arange(1, MatrixScale):
    CGTP_MATRIX_CH4[i][i-1] = -Value

# Generate GWP100_MATRIX
GWP_100_MATRIX_CH4 = int(METRIC_CONSTANTS['GWP100']) * np.identity(MatrixScale)
GWP_100_MATRIX_N2O = int(METRIC_CONSTANTS['GWP100N2O']) * np.identity(MatrixScale)

# Generate GWP20_MATRIX
GWP_20_MATRIX_CH4 = int(METRIC_CONSTANTS['GWP20']) * np.identity(MatrixScale)
GWP_20_MATRIX_N2O = int(METRIC_CONSTANTS['GWP20N2O']) * np.identity(MatrixScale)

# Generate GTP100_MATRIX
GTP_100_MATRIX_CH4 = int(METRIC_CONSTANTS['GTP100']) * np.identity(MatrixScale)
GTP_100_MATRIX_N2O = int(METRIC_CONSTANTS['GTP100N2O']) * np.identity(MatrixScale)

# Generate GTP20_MATRIX
GTP_20_MATRIX_CH4 = int(METRIC_CONSTANTS['GTP20']) * np.identity(MatrixScale)
GTP_20_MATRIX_N2O = int(METRIC_CONSTANTS['GTP20N2O']) * np.identity(MatrixScale)

# Save Matrices
np.savetxt("Matrices/CFE_MATRIX_CH4.csv", CFE_MATRIX_CH4, delimiter=",")
np.savetxt("Matrices/GWP_STAR_MATRIX_CH4.csv", GWP_STAR_MATRIX_CH4, delimiter=",")
np.savetxt("Matrices/CGWP_MATRIX_CH4.csv", CGWP_MATRIX_CH4, delimiter=",")
np.savetxt("Matrices/CGTP_MATRIX_CH4.csv", CGTP_MATRIX_CH4, delimiter=",")
np.savetxt("Matrices/GWP_100_MATRIX_CH4.csv", GWP_100_MATRIX_CH4, delimiter=",")
np.savetxt("Matrices/GWP_20_MATRIX_CH4.csv", GWP_20_MATRIX_CH4, delimiter=",")
np.savetxt("Matrices/GTP_100_MATRIX_CH4.csv", GTP_100_MATRIX_CH4, delimiter=",")
np.savetxt("Matrices/GTP_20_MATRIX_CH4.csv", GTP_20_MATRIX_CH4, delimiter=",")
np.savetxt("Matrices/GWP_100_MATRIX_N2O.csv", GWP_100_MATRIX_N2O, delimiter=",")
np.savetxt("Matrices/GWP_20_MATRIX_N2O.csv", GWP_20_MATRIX_N2O, delimiter=",")
np.savetxt("Matrices/GTP_100_MATRIX_N2O.csv", GTP_100_MATRIX_N2O, delimiter=",")
np.savetxt("Matrices/GTP_20_MATRIX_N2O.csv", GTP_20_MATRIX_N2O, delimiter=",")
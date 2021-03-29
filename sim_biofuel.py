
"""
%
% Purpose:  The aim of this file is to simulate the biofuel model
%
% Inputs:
%   data_set_to_use     (integer) Set id for system constants needed for simulation
%   time_array              (array) An array of uniformly distributed time
%                           instances
%   init_bacteria_amount    The initial amount of bacteria
%   alpha_b                 Value of the parameter alpha_b which is the
%                           production rate of biofuel
%   alpha_p                 Value of the parameter alpha_p which is the
%                           production rate of efflux pumps
%
% Outputs:
%   bacteria_amount_array   (array) Amount of Bacteria
%   sensor_array            (array) Sensor output
%   pump_array              (array) Number of efflux pumps
%   biofuel_int_array       (array) Amount of biofuel inside bacteria
%   biofuel_ext_array       (array) Amount of biofuel outside of bacteria

"""

import biofuel_system_parameter_sets as bsps
import numpy as np


def sim_biofuel(data_set_to_use, time_array, init_bacteria_amount, alpha_b, alpha_p):

    # BEGIN - DO NOT REMOVE
    # Note: Please do not remove this
    sys_para = bsps.biofuel_system_parameter_sets(data_set_to_use)
    ALPHA_N = sys_para['ALPHA_N']  # Growth rate (1/h)
    ALPHA_R = sys_para['ALPHA_R']  # Basal repressor production rate (1/h)
    BETA_R = sys_para['BETA_R']    # Repressor degradation rate (1/h)
    BETA_P = sys_para['BETA_P']    # Pump degradation rate (1/h)
    DELTA_N = sys_para['DELTA_N']   # Biofuel toxicity coefficient (1/(Mh))
    DELTA_B = sys_para['DELTA_B']   # Biofuel export rate per pump (1/(Mh))
    GAMMA_P = sys_para['GAMMA_P']   # Pump toxicity threshold
    GAMMA_I = sys_para['GAMMA_I']   # Inducer saturation threshold (M)
    GAMMA_R = sys_para['GAMMA_R']   # Repressor saturation threshold
    K_R = sys_para['K_R']       # Repressor activation constant (h)
    K_P = sys_para['K_P']       # Pump activation constant (1/h)
    K_B = sys_para['K_B']       # Repressor deactivation constant (1/M)
    V = sys_para['V']         # Ratio of intra to extracellular volume
    I = sys_para['I']         # Amount of inducer
    # The above lines set the following constants:
    # ALPHA_N ALPHA_R BETA_R  BETA_P  DELTA_N DELTA_B GAMMA_P GAMMA_I
    # GAMMA_R K_R K_P K_B V I
    # END - DO NOT REMOVE

    # You should put your work below this line

    # create output arrays that have the same length as time_array
    bacteria_amount_array = np.zeros_like(time_array)
    sensor_array = np.zeros_like(time_array)
    pump_array = np.zeros_like(time_array)
    biofuel_int_array = np.zeros_like(time_array)
    biofuel_ext_array = np.zeros_like(time_array)

    # set first element of bacteria_amount_array to start num of bacteria
    bacteria_amount_array[0] = init_bacteria_amount

    # we need time increment for equation-->  represented by delta triangle
    delta = time_array[1]-time_array[0]

    # for loop over each array --> use len(time_array)-1 due to 0 index
    for i in range(len(time_array)-1):

        # bacteria amount --> use equation 1 (Euler's forward method)
        bacteria_amount_array[i+1] = bacteria_amount_array[i] + \
            (ALPHA_N * (1-bacteria_amount_array[i]) - DELTA_N * biofuel_int_array[i]
             - (ALPHA_N * pump_array[i])/(pump_array[i] + GAMMA_P)) * \
            bacteria_amount_array[i] * delta

        # sensor array --> must be equation number 2
        sensor_array[i + 1] = sensor_array[i] + (ALPHA_R + K_R*(I/(I+GAMMA_I)) -
                                                 BETA_R * sensor_array[i]) * delta

        # pump array --> equation 3
        pump_array[i+1] = pump_array[i] + (alpha_p + K_P * (1/((sensor_array[i] /
                                                                (1 + K_B * biofuel_int_array[i])) + GAMMA_R))
                                           - BETA_P * pump_array[i]) * delta

        # biofuel internal array --> equation 4
        biofuel_int_array[i+1] = biofuel_int_array[i] + (alpha_b * bacteria_amount_array[i]
                                                         - DELTA_B * pump_array[i] * biofuel_int_array[i]) * delta

        # biofuel external array --> equation 5
        biofuel_ext_array[i+1] = biofuel_ext_array[i] + (V * DELTA_B * pump_array[i] *
                                                         biofuel_int_array[i] *
                                                         bacteria_amount_array[i]
                                                         * delta)

    return bacteria_amount_array, sensor_array, pump_array, biofuel_int_array, biofuel_ext_array


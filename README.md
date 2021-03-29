# Biofuel-Simulation
Model of Biofuel Production (SImulation)


Model aims to create a simulation of a biofuel production system. The simulation is created using different design parameters. Once the simulation is complete, the performance of the system is evaluated and the design parameters are chosen that maximize the amount of biofuel created. Completed as coursework for Masters of Health Data Science at University of New South Wales.

Scripts included :

sim_biofuel.py : Define a function that can be used to simulate the biofuel production process.

find_max_and_oscillation.py : Define a function that finds the maximum amount of internal biofuel and the amount of oscillation in the level of internal biofuel.

generate.py : Define a function that uses many different pairs of alpha_b and alpha_p to simulate the biofuel production

Design.py : Define a function that chooses the best designs of alpha_b and alpha_p parameters based on the most external biofuel product. The function chooses a best design based on threshold criteria and a poor design based only on output measures.

biofuel_simulation_design_parameter_sets.py : Defines constants and starting parameter values

biofuel_system_parameter_sets.py :  Defines constants and starting parameter values

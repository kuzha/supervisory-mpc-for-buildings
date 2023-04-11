# Supervisory Model Predictive Control for Buildings
Model predictive control (MPC) is a very promising control approach for energy
efficiency and demand response for building energy systems. This repository
includes example files for co-simulation between TRNSYS and Matlab.

If you are a building energy modeller who is not yet familiar with MPC, this
repository is a good start for you to try out MPC strategies applied to
buildings. A simple yet general supervisory control scheme is employed in the
controller in Matlab.

In this repository, a detailed building (white box) model is built in TRNSYS
which acts as the virtual building, while a Resistance-Capacitance (grey box)
model is realized in Matlab. For the virtual building model, other building
performance simulation tools such as EnergyPlus and ESP-r can also be used for
the co-simulation. A python-based controller will be added later.


# Citation
If you find the scripts useful, please kindly cite the following work:

Zhang, K., Kummert, M. Evaluating the impact of thermostat control strategies on the energy flexibility of residential buildings for space heating. Build. Simul. (2021). https://doi.org/10.1007/s12273-020-0751-x

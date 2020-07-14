# EKF_MagBias
This repository contains code for calibration of Magnetometers in Inertial Measurement Units (IMUs), with applications to Underwater vehicles. Magnetometer calibration is essential for accurate navigation in GPS denied environments, such as underwater. A Extended Kalman Filter framework is used to estimate the biases in magnetometer signals as well as in angular-rate sensor signals. These values are then used for heading estimation, resulting in significantly more accurate vehicle navigation. This code can be used with synthetic simulation data as well as with experimental data for a sensor on a vehicle. 

Some results of experimental evaluation with sensor data are shown in the figures below. The first figure shows that the biases converge to roughly constant values in about 15 minutes time. 

![Biases](/Figures/exp1_biases_18.png)

The second figure shows that the estimated heading, pitch and roll of the vehicle are very close to the true heading, pitch and roll.
![Heading](/Figures/exp1_att_err_18.png)

## To Run this code: 

- To run the code with simulated data, modify gen_samples.m as required. 

- To run the code with experimental data, call the data file using the read_microstrain or read_kvh functions. 

- The folder EKF_data contains all the data used for experimental evaluation as well as some notes and helper functions for reading the data.

- The folder noise_plots contains analysis of the noise statistics of the LORD Microstrain 3dm-gx5-25 used.

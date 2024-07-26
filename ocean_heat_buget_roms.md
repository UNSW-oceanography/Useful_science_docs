# Ocean Heat Budget in ROMS
### Author: Fernando Sobral
### Date: 26/July/2024
### Summary
This document aims to provide useful information about the heat budget in ROMS and best practices after a long process of analysis of the model outputs.

---

## Table of Contents
1. [Introduction](#introduction)
2. [The ROMS diagnostic and average outputs](#the-roms-diagnostic-and-average-outputs)
3. [The heat budget equation](#the-heat-budget-equation)
4. [Huon_temp and Hvom_temp and temp_hadv relationship](#huon_temp-and-hvom_temp-and-temp_hadv-relationship)
5. [References](#references)



---

## Introduction
The Ocean Heat Budget (OHB) is a usefull approach to get information about the ocean heat drivers. It is meant to 
give you an accurate decomposition of the terms controlling the temperature change over time in a specific
area. By design, the ROMS model satisfy the heat conservation equation at every grid poing. The time evolution of temperature in the ocean is given by the sum of net heat exchange with the atmosphere, divergence of advective heat transport by horizontal and vertical velocities, and three-dimensional diffusive processes:

$$
\frac{\partial T}{\partial t} = \frac{\partial Q}{\partial z} - \rho_0 c_p \Big[ \textbf{u} \cdot \nabla T - \Big( \kappa_H \nabla^2_H T + \kappa_z \frac{\partial^2 T}{\partial z^2} + K^{turb}_T \Big)\Big ]
$$



## The ROMS diagnostic and average outputs
In terms of ROMS, this equation is in the following format:

temp\_rate = temp\_xadv + temp\_yadv + temp\_vadv + temp\_xdiff + temp\_ydiff + \int temp\_vdiff dz





## The heat budget equation


## Huon_temp and Hvom_temp and temp_hadv relationship
[Summarize the main points of the document and provide any concluding remarks.]

## References
[1] [Reference 1]
[2] [Reference 2]

---

### Appendix (Optional)
[Include any additional information or appendices here.]

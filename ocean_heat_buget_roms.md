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
5. [The diffusion term](#the-diffusion-term)
6. [Good practices when analysing the data](#good-practices-when-analysing-the-data)
<!-- 5. [References](#references) -->



---

## Introduction
The Ocean Heat Budget (OHB) is a usefull approach to provide information about the ocean heat drivers. It is meant to give you an accurate decomposition of the terms controlling the temperature change over time in a specific area. 

<br>
<br>

## The heat budget equation

By design, the ROMS model satisfy the heat conservation equation at every grid poing. The time evolution of temperature in the ocean is given by the sum of net heat exchange with the atmosphere, divergence of advective heat transport by horizontal and vertical velocities, and three-dimensional diffusive processes:

$$
\frac{\partial T}{\partial t} = \frac{\partial Q}{\partial z} - \rho_0 c_p \Big[ \textbf{u} \cdot \nabla T - \Big( \kappa_H \nabla^2_H T + \kappa_z \frac{\partial^2 T}{\partial z^2} + K^{turb}_T \Big)\Big ]
$$


<br>
<br>

## The ROMS diagnostic and average outputs
To obtain the OHB terms in ROMS you must use the diagnostic output. This output provides you all the necessary terms to close the heat budget.
In terms of ROMS diagnostic outputs, the temperature rate of change is represented by the following equation:

temp_rate = temp_xadv + temp_yadv + temp_vadv + temp_xdiff + temp_ydiff + temp_vdiff

Using these terms, you can close the budget in any selected area within your domain.
You can replace temp_xadv + temp_yadv by temp_hadv and the same for diffusion, temp_xdiff + temp_ydiff by temp_hdiff. Following is a image showing the equivalency of temp_xadv+temp_yadv = temp_hadv.

![Reconstructing temp_hadv with temp_xadv and temp_yadv](images/ohb_images/hadv_xadv_yadv.png)
*Fig: Calculating temp_hadv based on the x and y terms.*

The temp_rate variable is provided, but you can also calculate it using all the variables above. 

```python
volume = dx * dy * dz

# The order is: horizontal advection (x and y), vertical adv, horizontal (x and y) and vertical diffusitity (z)
# Multiplying per volume and summing the 3 dimensions: integral
mine_heat_budget = (diag.temp_xadv * volume).sum(['s_rho', 'xi_rho', 'eta_rho']) + (diag.temp_yadv * volume).sum(['s_rho', 'xi_rho', 'eta_rho']) + \
                    (diag.temp_vadv * volume).sum(['s_rho', 'xi_rho', 'eta_rho']) + \
                    (diag.temp_xdiff * volume).sum(['s_rho', 'xi_rho', 'eta_rho']) + (diag.temp_ydiff * volume).sum(['s_rho', 'xi_rho', 'eta_rho']) + \
                    (diag.temp_vdiff * volume).sum(['s_rho', 'xi_rho', 'eta_rho'])


# Temp rate is the temperature tendency output from the model. 
# It's a daily 3D field and to compare with my calculation I also have to integrate it
output_temp_tendency = (diag.temp_rate * volume).sum(['s_rho', 'xi_rho', 'eta_rho'])

```

<br>
<br>

![Comparing temp_rate](images/ohb_images/temp_rate_comparison.png)
*Fig: Comparison between output temp_rate and calcualted by hand using the diagnostic terms.*



<br>
<br>

About the air-sea heat flux in [ROMS forum](https://www.myroms.org/forum/viewtopic.php?t=2420).
Because the air-sea heat flux is applied as the surface boundary condition to temp_vdiff and so is already included in the vertical divergence.
> If you vertically integrate temp_vdiff = d/dz*(K_v*dT/dz) between the limits z = -h and z = zeta then you should simply get K_v*dT/dz |z = zeta minus K_v*dT/dz |z = -h which is shflux/(Cp*rho0) minus 0. 

```python
Cp = 4181.3

# Model air-sea heat flux
((ds.shflux / (ds.rho0 * Cp)) * dx * dy).sum(['xi_rho', 'eta_rho']).plot(label='output')

# My calculation of air-sea flux
(diag.temp_vdiff * volume).sum(['xi_rho', 'eta_rho', 's_rho']).plot(label='mine')

plt.title('Comparisson between model output air-sea heat flux and calculated air-sea heat flux')
plt.legend()

```

<br>
<br>

![air-sea-flux](images/ohb_images/air_sea_flux_comparison.png)
*Fig: Air-Sea heat flux Volume integrated for the whole domain over time.*

<br>
<br>

```python
# Comparing the vertical diffusivity with the air sea heat flux
# Firstly over time, So I have to integrate horizontally the model output
Cp = 4181.3

fig, ax = plt.subplots(nrows=2)
# Model air-sea heat flux
(ds.shflux/(Cp * ds.rho0)).mean('ocean_time').plot(ax=ax[0])
ax[0].set_title('Model shflux')
# My calculation of air-sea flux
(diag.temp_vdiff * dz).mean('ocean_time').sum(['s_rho']).plot(ax=ax[1])
ax[1].set_title('My calculation')

plt.suptitle('Air-sea heat flux')
fig.tight_layout()

```

<br>
<br>

![air-sea-flux](images/ohb_images/air_sea_flux_comp_spatial.png)
*Fig: Time-mean air-Sea heat flux vertical integrated for the whole domain.*


<br>
<br>

## Huon_temp and Hvom_temp and temp_hadv relationship
One of the big lessons while trying to compare variables calculated in AVG (in this case Huon_temp and Hvom_temp) with variables calculated in DIA (like temp_xadv and temp_yadv) is that they are not comparable.

There is a good topic in [ROMS' forum](https://www.myroms.org/forum/viewtopic.php?t=5481) about it.

There are important differences while theses variables are being calculated when runing the model.
>The _xadv and _yadv terms (and their divergence which is already saved for you in the companion _hadv diagnostics) are the fluxes through the faces (time averaged) exactly as ROMS computed them according to the selected advection scheme. In the case of the high order Akima and weighted-upwind schemes, these fluxes are computed over a 3 or 4 grid cell stencil so their divergence is not a simple difference of the u*temp terms on the faces of a single cell. 
>Moreover, the ROMS time varying vertical s-coordinate means that the cell thickness, H, and hence cell face area itself (H/n) varies with time on every time step. Thus the time average <u*temp> multiplied by the time average layer area <H/n> is not exactly equal to <H/n*u*temp> because the triple nonlinearity of the perturbations <H'u'temp'> is not zero.

Before knowing that, some tests were performed and even the result is quite the same, the two outputs have similarities.

```python



```

<br>
<br>


## The diffusion term


<br>
<br>

## Good practices when analysing the data
- Masking land with NaN, to make sure it won't influentiate when integrating.
- Calculate your own temp_rate and comparing with the model output
- Calculate your own air-sea flux as the vertical integral of the temp_vdiff and compare with the air-sea flux provided in the average data (to be confirmed...)
- Because the DIA files can't be read with xroms, since it does not contain the zeta variable (for the version used as reference while writing this document), to easily have the metrics calculated, read the AVG file with xroms and assign the metrics to the DIA files read with xarray.



<!-- ## References
[1] [Reference 1]
[2] [Reference 2] -->

---
<!-- 
### Appendix (Optional)
[Include any additional information or appendices here.] -->

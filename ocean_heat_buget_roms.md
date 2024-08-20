# Ocean Heat Budget in ROMS
### Author: Fernando Sobral
### Date: 26/July/2024
### Summary

This document aims to provide useful information about the heat budget in ROMS and best practices following a lengthy process of analyzing the model outputs. It not only shows the correspondence of Ocean Heat Budget (OHB) terms, but also what can’t be correlated. The key takeaway is that there isn’t a straightforward correlation between the diagnostic terms (DIA) and the averaged terms (AVG) due to differences in the advection and time-stepping schemes, so they won’t be proportional as provided. Therefore, it is easier to use them separately. 


---

## Table of Contents
0. [The problem being addressed](#the-problem-being-addressed)
1. [Introduction](#introduction)
2. [The heat budget equation](#the-heat-budget-equation)
3. [The ROMS diagnostic](#the-roms-diagnostic)
    - [The Air-sea heat flux](#the-air-sea-heat-flux)
4. [Huon_temp and Hvom_temp and temp_hadv relationship](#huon_temp-and-hvom_temp-and-temp_hadv-relationship)
5. [ The horizontal diffusion term calculation](#the-horizontal-diffusion-term-calculation)
6. [Good practices when analysing the data](#good-practices-when-analysing-the-data)


---

## The problem being addressed

The main challenge that initiated the analysis was to find the relationship between the divergent diagnostic terms and the flux terms (Huon_temp and Hvom_temp) from the averaged output. If this correlation could be established, it would be possible to use the flux terms on the grid faces with a proportional correspondence to the change in temperature at the grid cell center. The application of this correlation is that it could provide confidence that a cross-contour heat transport amount is correct because it conserves the heat within the system.

<br>


## Introduction
The Ocean Heat Budget (OHB) is a useful approach to provide information about the drivers of ocean heat. It is designed to offer an accurate decomposition of the terms controlling the temperature change over time in a specific area.

This document contains relevant discussions that I have had with Neil Malan, Fabio Dias, Ryan Holmes, and John Wilkin.

<br>

## The heat budget equation

By design, the ROMS model satisfies the heat conservation equation at every grid point. The time evolution of temperature in the ocean is given by the sum of net heat exchange with the atmosphere, divergence of advective heat transport by horizontal and vertical velocities, and three-dimensional diffusive processes:

$$
\frac{\partial T}{\partial t} = \frac{\partial Q}{\partial z} - \rho_0 c_p \Big[ \textbf{u} \cdot \nabla T - \Big( \kappa_H \nabla^2_H T + \kappa_z \frac{\partial^2 T}{\partial z^2} + K^{turb}_T \Big)\Big ]
$$


<br>

## The ROMS diagnostic
To obtain the OHB terms in ROMS, you must use the diagnostic output. This output provides all the necessary terms to close the heat budget. When considering the diagnostic outputs provided by ROMS, the rate of temperature change is represented by the following equation:

temp_rate = temp_xadv + temp_yadv + temp_vadv + temp_xdiff + temp_ydiff + temp_vdiff

Using these terms, you can close the budget in any selected area within your domain. You can replace temp_xadv + temp_yadv with temp_hadv and similarly for diffusion, replace temp_xdiff + temp_ydiff with temp_hdiff. Following is an image showing the equivalency of temp_xadv + temp_yadv = temp_hadv.

<br>

![Reconstructing temp_hadv with temp_xadv and temp_yadv](images/ohb_images/hadv_xadv_yadv.png)
*Fig: Calculating temp_hadv based on the x and y terms.*

<br>

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
*Fig: Comparison between the volume integrated output temp_rate and calculated by hand using the diagnostic terms.*



<br>
<br>

### The Air-sea heat flux

About the air-sea heat flux in [ROMS forum](https://www.myroms.org/forum/viewtopic.php?t=2420).
Because the air-sea heat flux is applied as the surface boundary condition to temp_vdiff and so is already included in the vertical divergence.
> If you vertically integrate temp_vdiff = d/dz*(K_v*dT/dz) between the limits z = -h and z = zeta then you should simply get K_v*dT/dz |z = zeta minus K_v*dT/dz |z = -h which is shflux/(Cp*rho0) minus 0. 

The small differences is likelly due to the wrong Cp value used here, that was an approximation. The proper value used in the model can be found in [Good practices when analysing the data](#specific-heat-capacity).

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

Before knowing that, some tests were performed and even that the result are not quite the same, the two outputs have similarities. I did some analysis to compare them and the following plot shows the area used for that. It was chosen a shelf water area.

<br>
<br>

![area analysed](images/ohb_images/area_analysed.png)
*Fig: First pannel to the left to show that it was chosen columns and the other pannels shows the area used for the analysis.*

<br>
<br>

I made sure that I was taking the right grid cells for the analysis.

<br>
<br>

![grid selected](images/ohb_images/grid_selected_points.png)

*Fig: Grid cells selected for both temp_yadv and Hvom_temp.*

<br>
<br>

And the next plots are the comparison calculated as the vertical integral and then calculating the divergent and the difference with the diagnostic divergent. 
The red circles are from the temp_yadv and black circles from Hvom_temp. From the difference, the error is of less than 1% for these examples. 


```python
# All dots selected
def new_fig():
    fig, ax = plt.subplots(ncols=3, figsize=(15, 5))
    ax[0].scatter(lon_v, lat_v, label='Hvom')
    ax[0].scatter(lon_rho[1:, 1:], lat_rho[1:, 1:], label='Yadv')
    ax[0].set_title('Both SLICED together');
    ax[0].legend();
    return fig, ax


# Selected area to be used in the calculation
hvom0 = hvom_temp.isel(xi_rho=slice(89, 95)).isel(eta_v=slice(0, None)).sum('s_rho')
yadv0 = (temp_yadv * (dx * dy * dz)).isel(xi_rho=slice(89, 95), eta_rho=slice(1, None)).sum('s_rho')

for down, up in zip(range(0, hvom0.eta_v.shape[0]), range(1, hvom0.eta_v.shape[0])):
    print(down, up)
    fig, ax = new_fig()
    hvom00 = hvom0.isel(eta_v=down)
    hvom01 = hvom0.isel(eta_v=up)
    yadv00 = yadv0.isel(eta_rho=down)
    ax[0].plot(hvom00.lon_v, hvom00.lat_v, 'ok')
    ax[0].plot(hvom01.lon_v, hvom01.lat_v, 'ok')
    ax[0].plot(yadv00.lon_rho, yadv00.lat_rho, 'or')

    # Calculations
    dvg = hvom00 - hvom01

    # Comparing plot
    ax[1].plot(dvg, 'o-k', label='Div Hvom')
    ax[1].plot(yadv00, '--or', label='yadv')
    # ax[0].plot(yadv00.lon_rho, yadv00.lat_rho, 'or')
    ax[1].set_title('Comparing')
    ax[1].legend()

    # Difference
    diff_dvg = yadv00 - dvg
    ax[2].plot(diff_dvg, '-ok')
    ax[2].set_title('Difference yadv - div Hvom')
    fig.tight_layout()
```


<br>
<br>

![comp1](images/ohb_images/comp1.png)
![comp2](images/ohb_images/comp2.png)
*Fig: comparison between temp_yadv and Hvom_temp. The x-axis for the subplot in the center and to the right represents the six grid cells being compared.*
<br>
<br>

And the comparison over time to see if there was an accumulative factor. Which can't be seen.

<br>
<br>

![comp_overtime](images/ohb_images/comp_overtime.png)
*Fig: comparison between temp_yadv and Hvom_temp over time.*

<br>
<br>

## The horizontal diffusion term calculation
Based on all the information provided, how would be possible to calculate the heat transport across contours? In this case, only Huon/Hvom_temp aren't enough. We also need the diffusion terms. And for those, we need to calculate them. 

First of all, you need to figure out, which scheme you have used in the model setup. And from that, calculate [manually the diffusion](https://www.myroms.org/wiki/Horizontal_Mixing#Horizontal_Diffusion). 

In my case, I have to use the Laplacian horizontal diffusion. A way to find this is looking the cpp defs. My model has TS_DIF2 (harmonic mixing tracers) and if biharmonic, it would have defined TS_DIF4 ([have a look into the CPP defs](https://www.myroms.org/wiki/cppdefs.h#Options_for_horizontal_mixing_of_tracers)).


$$\frac{\partial}{\partial \xi}\Big(\frac{\nu_2H_zm}{n}\frac{\partial C}{\partial \xi}\Big) + \frac{\partial}{\partial \eta}\Big(\frac{\nu_2H_zn}{m}\frac{\partial C}{\partial \eta}\Big)$$


Important concepts to have in mind before you start the calculation. 
1. You need to consider the s-layers when doing a horizontal derivative in ROMS. You need to calculate the gradient for the same z-layer. The easiest way of doing it is to use the already made tools. In my case, I used xroms ddxi and ddeta.

2. The notation $m$ and $n$ here is equivalent to $pm$ and $pn$, respectivelly. Given:

$$
dx = \frac{1}{pm} \ \mathrm{or} \ \frac{1}{m} 
$$

and

$$
dy = \frac{1}{pn} \ \mathrm{or} \ \frac{1}{n}
$$

3. $\xi$ and $\eta$ are non-dimension in the equation and what gives them dimension is the multiplication of:

$$
m \frac{\partial C}{\partial \xi}
$$

and 

$$
n \frac{\partial C}{\partial \eta}
$$

or using dx and dy:

$$
\frac{1}{dx} \frac{\partial C}{\partial \xi}
$$

and

$$
\frac{1}{dy} \frac{\partial C}{\partial \eta}
$$

4. Always check the units! If your goal it to calculate the diffusive flux, them it need $m^3s^{-1} \degree C$. For the divergence and therefore to be able to do your sanity check, you'll need units like $\degree C s^{-1}$.

You will need the diffusive coefficient $\nu2$ value, and can find it in the *ocean.in* file, and it is equal to 55 in my case.


Here is the code I use to calculate the diffusive flux and the divergence.

```python

# Value for the viscosity horizontal diffusion coefficient (from ocean.in)
nu2 = 55


'''-------------------------------------
First I calculate the gradient. And for this, I am already considering dx, dy for the x- and y-derivatives
-------------------------------------'''

# Partial derivatives in XI and ETA.

dtemp_dxi_ = xroms.ddxi(ds_avg.temp, xgrid, scoord='s_rho')
dtemp_deta_ = xroms.ddeta(ds_avg.temp, xgrid, scoord='s_rho')


'''-------------------------------------
Second multiplying by the diffusive coefficient and the respective areas: dz*dy for x-direction and dz*dx for y-direction. Before, I was also dividing by dx;dy here, because I didn't consider that xroms has already consider dx;dy areas. At this step, my units are m3 deg C/s
-------------------------------------'''

# This is the calculation of the horizontal diffusive flux

term1_ = ((nu2 * ds_avg.dz_u * ds_avg.dy_u)) * dtemp_dxi_
term2_ = ((nu2 * ds_avg.dz_v * ds_avg.dx_v)) * dtemp_deta_


'''-------------------------------------
And then calculating the divergence of the diffusive flux. Again, it already consider dx;dy when doing the horizontal derivatives. At this step, I have units like m2 deg C/s
----------------------------------------'''

# Divergence of the diffusive flux
dterm1_dxi_ = xroms.ddxi(term1_, xgrid, scoord='s_rho', attrs=term1_.attrs)
dterm2_deta_ = xroms.ddeta(term2_, xgrid, scoord='s_rho', attrs=term2_.attrs)


'''-------------------------------------
And finally dividing by the respective areas to get deg C/s. Before, I thought that I would also needed to divide by dx the xdiff_temp_ and by dy the ydiff_temp_. But that is wrong since, again, xroms has done this already.
-------------------------------------'''

# Dividing by the metrics that havent participated in the divergence to get deg C/s units and to be comparable with temp_hdiff.

xdiff_temp_ = dterm1_dxi_ / (ds_avg.dz * ds_avg.dy)
ydiff_temp_ = dterm2_deta_ / (ds_avg.dz * ds_avg.dx)


# Summing xdiff and ydiff to get hdiff_temp
hdiff_temp_ = xdiff_temp_ + ydiff_temp_

```
<br>
<br>

![horizontal diff](images/ohb_images/comparing_hdiff.png)

*Fig: Rate of change comparisson between the temp_hdiff from diagnostic and my manual calculation shown here.*

John Wilkin words about the topic can be found on the [ROMS forum](https://www.myroms.org/forum/viewtopic.php?p=25568#p25568) and in private messages:


![horizontal diff](images/ohb_images/wilkin_PM_01.png)

<br>
<br>

![horizontal diff](images/ohb_images/wilkin_PM_02.png)

<br>
<br>

![horizontal diff](images/ohb_images/wilkin_PM_03.png)

<br>
<br>

![horizontal diff](images/ohb_images/wilkin_PM_04.png)




## Good practices when analysing the data
- Masking land with NaN, to make sure it won't influentiate when integrating.
- Calculate your own temp_rate and comparing with the model output
- Calculate your own air-sea flux as the vertical integral of the temp_vdiff and compare with the air-sea flux provided in the average data (to be confirmed...)
- Because the DIA files can't be read with xroms, since it does not contain the zeta variable (for the version used as reference while writing this document), to easily have the metrics calculated, read the AVG file with xroms and assign the metrics to the DIA files read with xarray.
<a name="specific-heat-capacity"></a>
- When calculating the heat you'll need the $\rho_0$ and $C_p$. These values have to come from the model, and you can find them in the file called /g/data/fu5/trunk/ROMS/Modules/mod_scalars.F on Gadi. Sometimes I used approximate values showed in this file.


---


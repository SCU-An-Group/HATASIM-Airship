# HATASIM-Airship

HATASIM-Airship is a simulation tool for predicting the quasi-static ascent and equilibrium attitude of tethered airships under wind loading. It couples the airship's attitude-dependent aerodynamic forces and pitching moment with the mechanical response of the tether. MATLAB and browser-based implementations are provided.

## What the tool calculates

For specified airship, tether, atmospheric, and wind parameters, the tool calculates equilibrium configurations at a sequence of airship altitudes. At each altitude it:

1. Evaluates the atmospheric conditions and the airship's buoyancy and weight.
2. Solves the airship moment balance for its equilibrium angle of attack, then obtains the tether tension and inclination at the confluence point from force balance.
3. Solves the equilibrium of successive elastic tether segments from the confluence point toward the ground.

The results include the airship attitude and horizontal drift, tether profile and tension, and aerodynamic and environmental quantities. The implementations also calculate local moment derivatives with respect to angle of attack and wind speed.

## Files

| File | Description |
| --- | --- |
| `Lotte_EN.html` | English-language browser interface. |
| `Lotte.html` | Chinese-language browser interface. |
| `Lotte.m` | MATLAB implementation with LOTTE airship parameters. |

## Getting started

### Browser interface

1. Open `Lotte_EN.html` for the English interface, or `Lotte.html` for the Chinese interface, in a modern web browser.
2. Enter the airship, tether, wind, and atmospheric parameters and the aerodynamic coefficient functions.
3. Select **Run Simulation**. The interface displays available initial solutions and the calculated ascent results. Use **Query** to examine a selected altitude.

The HTML files load Plotly from an online content delivery network. An internet connection is therefore needed to load the plotting library unless a local copy is configured.

### MATLAB

Open `Lotte.m` in MATLAB, review the parameters at the beginning of the file, and run the script. The code uses `fsolve` from MATLAB's Optimization Toolbox.

## Inputs and conventions

- Enter geometry, masses, tether properties, wind conditions, and aerodynamic data for the configuration being studied. The supplied defaults provide an example configuration and should be checked before use for another airship.
- Buoyancy-center and confluence-point coordinates are body-fixed offsets from the center of gravity. The positive body-frame x-axis points toward the nose, and the positive y-axis points upward.
- Aerodynamic coefficient inputs are polynomial coefficients in descending order of power. The polynomials are evaluated using the angle of attack in radians; angles displayed in the interface are in degrees.
- The lift, drag, and pitching-moment data must use reference dimensions consistent with the inputs. The pitching moment is represented about the center of gravity; coefficients referenced to another point require the corresponding moment transfer.
- The model assumes planar quasi-static equilibrium, a rigid airship of constant volume, a fixed confluence-point position relative to the airship, and a segmented elastic main tether. The computed moment derivatives describe local pitch response and wind sensitivity.

## Airship manuscript

The accompanying study, **“Predicting the ascent and attitude of tethered airships under wind loading,”** is in preparation. Its citation and the link to the HATASIM-Airship repository can be added when available.

## Credits and contact

HATASIM-Airship is developed by Zhongwei Ni and collaborators under the supervision of Ning An, Sichuan University. For correspondence, contact Ning An at anning@scu.edu.cn.

## Related project: HATASIM for spherical balloons

HATASIM-Airship extends the framework developed for [HATASIM](https://github.com/SCU-AN-Group/HATASIM), a separate program for the quasi-static ascent of tethered spherical balloons. The balloon framework is described in:

Ni, Z., Zhu, S., Dou, S., Song, J., and An, N. “Predicting the ascent of high-altitude tethered spherical balloons.” *Aerospace Science and Technology* 168 (2026), 111087. [Read the paper](https://www.sciencedirect.com/science/article/pii/S1270963825011502).

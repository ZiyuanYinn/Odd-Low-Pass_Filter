# Odd-Low-Pass-Filter
The Odd-Low-Pass Filter is a tool to remove smaller, odd-parity structures like the vertical waves, and use the filtered data to examine the location of the Galaxy's mid-plane.

## R_cuts_and_CMD_trim.py
A Python file that cuts the sample's R range and the color and magnitude ranges, producing trimmed.csv files for each bin of R = 0.1 kpc. Apply this code first before applying the errorbar.py

## OLPF.py (or errorbar.py)
A Python file that contains all OLPF codes, the average vertical height (for each radial bin) calculated with vertical waves removed, and the error bars of that calculation. The file also contains a way to investigate the effects of OLPF by subtracting the OLPF'd data from the raw data.

## Verlet Integration Simulation.py
A simulation using Verlet integration that plots the motion of stars in the vertical direction after every timestep.

## analysis without OLPF.py
A mid-plane analysis without OLPF.

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PHYS20161 - First Assignment: Forced Oscillations

This code describes the motion of a system undergoing forced oscillations.

Input parameters are: a fixed constant A_0 and a parameter a_1, which describe the initial
amplitude and the attenuation coefficient of the system, frequency of the oscillations,
and a minimum fractional intensity intensity_frac_min, where the fractional intensity is
the intensity divided by the maximum intensity.

Amplitude A as a function of time:  A(time) = 1 / (A0 + a1 * time^2) * cos(a2 * time),
where a2 = 2 * pi * f - angular frequency.

Fractional intensity at a given time is equal to (A(time))^2 / (A(0))^2.
Hence, formula for fractional intensity:
intensity_frac(time) = (1 / (A_0 + a_1 * time^2) * cos(a_2 * time))^2 / (1 / A_0)^2

The code gives out the number of oscillations, where the fractional intensity is
larger than intensity_frac_min, and the time time_of_osciallations it takes to complete
the oscillations. It also plots the Fractional Intensity vs Time plot for the given
system.

Fedir Boreiko, 21/10/2023
"""

import time as t
from math import pi, cos, tan, log
import numpy as np
import matplotlib.pyplot as plt

A_0 = 2.0

def input_validated(param_name:str, limits=str, param_units=None):
    """
    Takes an input value of the needed parameter (float by default) and checks
    wether the entered value satisfies the limits, if any given.
    """
    parameter_not_valid = True

    while parameter_not_valid:

        if param_units is not None:

            parameter = input(f'\nEnter the value of {param_name} in {param_units}: ')

        else:

            parameter = input(f'\nEnter the value of {param_name}: ')

        parameter = float(parameter)

        if limits is not None:

            limit_values = list(map(float, limits[1:-1].split(', ')))

            limit_1, limit_2 = limit_values

            condition = {
                '[]': limit_1 <= parameter <= limit_2,
                '()': limit_1 < parameter < limit_2,
                '[)': limit_1 <= parameter < limit_2,
                '(]': limit_1 < parameter <= limit_2
            }[limits[0] + limits[-1]]

            if condition:

                parameter_not_valid = False

            else:

                if param_units is not None:

                    print(f"\nValue of {param_name} should lay in the range "
                          f"{limits[0]}{limit_1}, {limit_2}{limits[-1]} {param_units}."
                          "\nPlease enter the value again.")

                else:

                    print(f"\nValue of {param_name} should lay in the range "
                          f"{limits[0]}{limit_1}, {limit_2}{limits[-1]}."
                          "\nPlease enter the value again.")

    return parameter


def intensity_fractional(time:float):
    """
    Evaluated fractional intensity at a given time.
    """

    i_frac = (1 / (A_0 + a_1 * time**2) * cos(a_2 * time))**2 / (1/A_0)**2

    return i_frac

def intensity_fractional_numpy(time:tuple):
    """
    Returns np.array of fractional intensities corresponding to the given time array.
    """

    i_frac = (1 / (A_0 + a_1 * np.power(time, 2)) * np.cos(a_2 * time))**2 / (1/A_0)**2

    return i_frac

def intensity_derivative(time:float):
    """
    Evaluated part of the derivative of the fractional intensity function in order to find points
    of local maxima.
    """

    derivative = 2 * a_1 * time / (a_2 * (2 + a_1 * time**2)) + tan(a_2 * time)

    return derivative

def bisection_iteration(endpoint_1:float, endpoint_2:float, function):
    """
    Given a function and two values bracketing the root this function, evaluates
    the mid-point and chooses the next two values bracketing the root.
    """

    mid_point = (endpoint_1 + endpoint_2) / 2

    if function(mid_point) * function(endpoint_1) < 0:

        endpoint_1_update = endpoint_1
        endpoint_2_update = mid_point

    else:

        endpoint_1_update = mid_point
        endpoint_2_update = endpoint_2

    return endpoint_1_update, endpoint_2_update

a_1 = input_validated(param_name='a1',
                      limits="(0.1, 50)",
                      param_units='m^(-1)s^(-2)')
frequency = input_validated(param_name='frequency of oscillations',
                            limits="[1, 200]",
                            param_units='Hz')
intensity_frac_min = input_validated(param_name='minimum fractional intensity',
                                     limits="(0, 1]")

start_execution = t.time()

a_2 = 2 * pi * frequency     #angular frequency
period = 1 / frequency       #period of oscillations

"""
We compute the time derivative of the fractional intensity function in order to find the points of
local minima and local maxima. Hence:

time = 1/a_2 * (num + 1/2) * pi,  where num = 0, 1, 2, ...  - points of local minima.

To obtain points of the local maxima the roots of following equation need to be found:

derivative(time) = 2 * a_1 * time / (a_2 * (2 + a_1 * time**2)) + tan(a_2 * time)

The roots will be found via bisection method.
"""

intensity_frac = 1
num = 0                      #current oscillation number

#number of bisection iterations needed for obtaining error 10^(-10) of local maxima value:
num_of_iterations = int((log(pi / a_2) - log(10**(-10))) / log(2))

print(f'\nNumber of bisection iterations: {num_of_iterations}')

if intensity_frac_min == 1:

    num_of_oscillations = 0
    time_of_osciallations = 0

else:

    while intensity_frac > intensity_frac_min:

        #checking whether the current local maxima of i_fractional function is larger than
        #the i_fractional_min and if so, proceeding to the next local maxima

        num += 1

        #for oscillations, each local maxima is located between two subsequent local minima
        time_of_current_min = (num - 1/2) * pi / a_2
        time_of_next_min = (num + 1/2) * pi / a_2

        time_narrowing = (time_of_next_min - time_of_current_min) / 100

        #since derivative(time) is undefined in local minimas of fractional intensity function,
        #decrease the derivative's root search area by 1% from each endpoint

        time_endpoint_1 = time_of_current_min + time_narrowing
        time_endpoint_2 = time_of_next_min - time_narrowing

        for iteration in range(num_of_iterations):
            #updating the time coordinates bracketing a derivative function's root

            time_endpoint_1, time_endpoint_2 = bisection_iteration(time_endpoint_1,
                                                                   time_endpoint_2,
                                                                   intensity_derivative)

        time_of_current_max = (time_endpoint_1 + time_endpoint_2) / 2

        intensity_frac = intensity_fractional(time_of_current_max)

    num_of_oscillations = num - 1                  #returning to the last valid oscillation's number
    time_of_osciallations = time_of_current_min    #obtaining the time coordinate of current minima

print("\nNumber of oscillations, where the fractional intensity is larger than minimum "
      f"fractional intensity: {num_of_oscillations} \nTime it takes to complete the "
      f"oscillations: {time_of_osciallations:.3f} sec")

#Plotting the Fractional Intensity vs Time plot for the system with given parameters,
#only if intensity_frac_min is not 1.

if intensity_frac_min != 1:

    time_frames = np.arange(0, time_of_osciallations + 2 * period, period / 100)

    intensity_frac_frames = intensity_fractional_numpy(time_frames)

    fix, ax = plt.subplots()

    plt.plot(time_frames, intensity_frac_frames, color='blue', label='I_fractional(t)')

    ax.axhline(y=intensity_frac_min, color='green', linestyle='dashed', label='I_fractional_min')

    ax.axvline(x=time_of_osciallations, color='red', linestyle='dashed', label='t_oscillations')

    ax.set_xlabel('Time (secs)')
    ax.set_ylabel('Fractional Intensity')

    ax.set_xlim([0, time_of_osciallations + 2 * period])
    ax.set_ylim([0, 1])

    ax.legend(loc=1, fontsize=8)

    plt.title('Fractional Intensity vs Time')

    plt.show()

end_execution = t.time()

execution_time = end_execution - start_execution

print(f"\nExecution time: {execution_time} sec")

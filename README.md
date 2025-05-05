# Gravitational Orbit Simulator

## Description

This MATLAB program simulates planetary orbits in universes with modified gravitational laws, where the gravitational force follows the formula F = G*m₁*m₂/rⁿ. The simulator allows exploration of how different values of the exponent 'n' affect orbital dynamics and stability.

Unlike our universe where n=2 (Newton's inverse square law), this simulator lets you investigate alternative gravitational physics and visualize the resulting orbits. The program includes:

- Interactive parameter input for gravitational exponent (n)
- Dynamic visualization of orbital trajectories
- Real-time calculation of orbital parameters (perihelion, aphelion, period)
- Energy diagram showing the effective potential for the selected parameters
- Bertrand's theorem analysis for orbital stability

## Scientific Background

The simulator is based on Bertrand's theorem, which proves that only specific values of 'n' produce stable, closed orbits, with n=2 (Newton's inverse square law) being the one we observe in our universe.

For other values of 'n', orbits will generally precess and not close upon themselves. This has profound implications for the possibility of life in universes with different gravitational laws.

## Features

- Numerically solves the equations of motion for arbitrary gravitational laws
- Visualizes trajectories with dynamic zooming as orbits evolve
- Displays real-time orbital parameters (distance, velocity, elapsed time)
- Color-coded indication of orbital stability based on Bertrand's theorem
- Energy diagram to visualize the effective potential well

## Usage

1. Run the program through MATLAB
2. Enter a value for the gravitational exponent 'n' between 1.904 and 2.764
3. The simulation will calculate and display the orbital parameters
4. An animated visualization will show the resulting orbit with a dynamic trail

The program provides a valuable educational tool for understanding how fundamental physical laws shape the possibility of stable planetary systems.

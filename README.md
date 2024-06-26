# Langevin Dynamics Code for Brownian Particles



## Table of Contents
- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [Acknowledgements](#acknowledgements)
- [Analysis](#analysis)

## Introduction

This C code uses Langevin Dynamics (under-damped Brownian Dynamics) to explore the spatial-temporal correlations of Brownian particles interacting via competitive short-ranged attractive and long-ranged repulsive potentials. Specifically, the interaction potential comprises the 100-50 (Lennard-Jones-like) Mie potential and the relatively long-ranged Yukawa (screened Coulomb) potential.

## Features

- **Parallelization**: The code is parallelized using OpenMP, providing decent computing speed for simulations involving tens of thousands of Brownian particles.
- **Static Structure Factor (S(q))**: Direct computation within the code.
- **Radial Distribution Function (g(r))**: Direct computation within the code.
- **Dynamical Trajectory Analysis**: Compute mean-squared displacement functions directly in the code.
- **Trajectories Analysis**: Obtain additional properties by analyzing the particle trajectories.

## Installation

### Prerequisites

Ensure you have a C compiler such as ICC or GCC installed on your system.

### Steps

1. Clone the repository using the permalink:
   ```bash
   git clone https://github.com/ziihan/Langevin_Dynamics_SALR/tree/1af7a3bcb3cd78005e03ccb661c0672046f7fd62/LD_code_SALR

2. Navigate to the project directory:
   cd Langevin_Dynamics_SALR/LD_code_SALR
   
3. Compile the code:
Using ICC:
make -f Makefile

Using GCC:
make -f Makefile_gcc

## Usage
It is a simple LD code and it is free to use. But please cite me if it contributes/inspires your research.

1. Run the executable generated after compilation
   ./langevin_dynamics
   
2. Input your simulation parameters as required by the program prompts.

## Acknowledgements
Zihan Tan developed the code. 
Zihan Tan, along with his collaborators Vania Calandrini, Jan K.G. Dhont, Gerhard Nägele, and Forschungszentrum Jülich (where he was employed), owns this code.

## Analysis
Inside the "LD_code_SALR" folder, there are the scripts for cluster analysis to compute: cluster size distribution (CSD), coordination number (z_b) and its distribution, local hexagonal order parameter |q_6^2| and its distribution and average. Those are developed based partially on my former colleagues' knowledge.

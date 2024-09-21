# Wildfire Simulator

### Author: Alberto Cuoci  
**CRECK Modeling Lab, Politecnico di Milano**  
Email: [alberto.cuoci@polimi.it](mailto:alberto.cuoci@polimi.it)  
Website: [CRECK Modeling Lab](http://creckmodeling.chem.polimi.it/)

---

## Overview

*Wildfire Simulator* is a computational tool designed to simulate the spread of wildfires in forest environments based on a **cellular automata (CA)** model, as proposed by Karafyllidis and Thanailakis in their 1997 paper: *"A model for predicting forest fire spreading using cellular automata"*.

### Key Features
- Simulation of fire spread in **homogeneous** and **inhomogeneous** forests.
- Inclusion of **weather conditions** (e.g., wind direction and speed) to affect fire propagation.
- Incorporation of **terrain topography** to model fire behavior on slopes.
- Discrete **cellular automata** approach for efficient simulations, suitable for real-time prediction and decision support.
- Designed for parallel computation and scalability on modern computing architectures.

## Cellular Automata Model

The wildfire spreading model uses a 2D grid of cells to represent the forest. Each cell has a state indicating the proportion of the cell that is burned. The state of each cell is influenced by:
1. **Rate of fire spread** in its own cell.
2. The state of **neighboring cells** (both adjacent and diagonal).
3. **Wind speed and direction**, which can accelerate or decelerate the spread in certain directions.
4. **Topography**, where uphill spread is faster and downhill spread is slower.

### Model Parameters

The CA model is defined by:
- A **square grid** of forest cells.
- **Local state** for each cell, representing the fraction of the area burned.
- **Fire spread rate (R)**, which depends on factors like vegetation type and terrain.
- **Transition rules** for updating the state of each cell based on its neighbors.

### Fire Spread Dynamics

At each timestep, the fire spreads from burning cells to neighboring cells based on the fire spread rate (R). The model differentiates between:
- **Adjacency burning**: Spread to neighboring cells sharing a side.
- **Diagonal burning**: Spread to cells sharing a corner, which takes longer due to distance.
- **Wind effects**: Fire spreads faster in the direction of the wind and slower against it.
- **Topographical effects**: Fire spreads faster uphill and slower downhill.

## Usage

The *Wildfire Simulator* uses input files to initialize the forest environment, weather conditions, and simulation parameters. The user can configure:
- Forest map (vegetation types, terrain elevation, wind speed and direction).
- Simulation end time and output intervals for visualization.
- Output folder for saving simulation results.

### Command Line Arguments

- `--map <path_to_xml>`: Specifies the XML file defining the forest map.
- `--tEnd <time>`: Sets the simulation end time (default is 1000s).
- `--dtMap <interval>`: Time interval for saving fire spread maps (default is 10s).
- `--dtLog <interval>`: Time interval for writing log files (default is 10s).
- `--dtScreen <interval>`: Time interval for printing to the screen (default is 1s).
- `--Co <Courant_number>`: Sets the Courant number for numerical stability (default is 0.25).
- `--outputFolder <path>`: Sets the folder where results will be saved.

### Example

```bash
./wildfire_simulator --map forest_map.xml --tEnd 2000 --dtMap 50 --dtLog 50 --dtScreen 10 --Co 0.5 --outputFolder results/
```

This will run the simulation with the specified settings, using the forest map provided in `forest_map.xml` and saving results to the `results/` folder.

## Output

The simulator produces various outputs including:

- **Fire spread maps**: Saved in VTK or Tecplot formats for easy visualization of fire progression.
- **Log files**: Containing detailed information about the simulation steps.
- **Screen output**: Real-time updates on the simulation progress.

---

## Installation

The simulator is built using C++ with dependencies on the Boost library and Eigen for matrix operations.

To compile:

```bash
mkdir build
cd build
cmake ..
make
```

Ensure that you have the necessary dependencies installed:

- **Boost**: For file handling, property trees, and command-line options.
- **Eigen**: For matrix operations.

---

## References

- **Karafyllidis, I., & Thanailakis, A.** (1997). *A model for predicting forest fire spreading using cellular automata*. Ecological Modelling, 99(1), 87-97.

---

## License

This project is licensed under the MIT License - see the LICENSE file for details.
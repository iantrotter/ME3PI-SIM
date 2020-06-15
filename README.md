# ME3PI-SIM: Macroeconomic-Epidemiological Emergency Policy Intervention Simulator

The ME3PI-SIM simulates **emergency policy interventions** in an **economy** under the impact of an **epidemic**.

## Introduction
This program runs the simulation model described in the paper [COVID-19 and Global Economic Growth: Policy Simulations with a Pandemic-Enabled Neoclassical Growth Model](https://arxiv.org/abs/2005.13722).

The paper integrates a variant of the SIR (Susceptible-Infected-Recovered) epidemiological model into a standard neoclassical model of economic growth, then runs a series of numerical experiments with "emergency policies", which simultaneously reduce the infection rate and the production of economic output. For each of the policies, we simulate the trajectory of the pandemic, and of the economic variables.

## Usage
The configuration for a simulation is stored in a JSON file. All the parameters are explained in detail in the paper, and we provide several examples in this repository.

The program is a Python script, so to run the simulation, navigate to the base repository directory in a terminal, and run the script using the Python command, with the JSON configuration file name as the argument, for example:
```
$ python ME3PI-SIM.py configuration_file.json
```

The simulation might take several minutes to run, depending on the parameters of the simulation. The output will be written to the output file (a Pickle of a Pandas dataframe, and an Excel spreadsheet with the same data) specified in the configuration file.

## Disclaimer
This is to be considered prototype code that has escaped from the laboratory.
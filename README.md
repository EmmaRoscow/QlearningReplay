README

Contains code and data for reproducing the analyses presented in Roscow et al. 2025. All code written in MATLAB, developed and tested in MATLAB R2016a.

Repo structure is as follows.

QlearningReplay/
│
├─ code/
│  ├─ ephys/  # Code for analysing electrophysiological data
│  │  ├─ tests/	 # Unit tests for the ephys code
│  │  └─ main_ephys_script.m  # Script for analysing ephys data from scratch
│  │
│  └─ q_learning_optimisation/  # Code for fitting Q-learning parameters to behavioural data
│
├─ code_for_reproducing_figures/
│  └─ utils/
│
├─ data/
│  ├─ behavioural_data/
│  ├─ ephys_data/
│  ├─ q_learning_raw_outputs/
│  ├─ q_learning_results/
│  └─ explainedVarianceReactivationAnalysis.mat	# Custom object containing the processed, analysed ephys data
│
├─ .gitignore
│
└─ README.md

Recommended use:
* Run all scripts from the base directory (rpe-replay), and add all folders and subfolders to the path
* To reproduce figures from the paper, execute the scripts in 'QlearningReplay/code_for_reproducing_figures/'
* To explore the Q-learning modelling, start from 'QlearningReplay/code/q_learning_optimisation/dynaQ.m'
* To explore the electrophysiological data analysis, load 'QlearningReplay/data/explainedVarianceReactivationAnalysis.mat'
* To run the Q-learning modelling from scratch, run the scripts with names starting 'runOptimisation...' in 'QlearningReplay/code/q_learning_optimisation'. Optimisation requires the MATLAB package bads (https://github.com/acerbilab/bads)
* To run the electrophysiological data analysis from scratch, run 'QlearningReplay/code/ephys/main_ephys_script.m'

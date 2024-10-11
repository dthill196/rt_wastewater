#!/bin/bash


# Set up the Conda environment in the current Bash session
eval "$(/home/$(whoami)/miniconda3/bin/conda shell.bash hook)"

# Extract argument
county_name=$1

# Specify the R script to run, provide with county name
Rscript rt_episewer_process.R "$county_name"
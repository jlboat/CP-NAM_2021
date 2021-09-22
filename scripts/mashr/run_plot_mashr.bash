#!/usr/bin/env bash

# singularity run -B /zfs ~/singularity_containers/sap_mashr_v0.1.sif Rscript MashRScript.R

singularity run -B /zfs ~/singularity_containers/mashr.sif Rscript plot_mash.R

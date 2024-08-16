#!/bin/bash

# Load the configuration file
source config.nf

# Execute the Nextflow pipeline
nextflow run ngs_workflow.nf -c config.nf


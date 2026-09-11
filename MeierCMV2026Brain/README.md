In the `Code` directory, you'll find various R scripts for final analyses.

# MRI processing
Data processing was done using Singularity (now Apptainer) containers with the following software versions:
- [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki) 7.2 (freesurfer-linux-centos6_x86_64-7.2.0-20210720-aa8f76b)
- [HALFpipe](https://github.com/HALFpipe/HALFpipe) 1.2.2 (1.2.2.post1.dev189+g5579a7f)
- [QSIPrep](https://github.com/PennLINC/qsiprep) 1.0.2 (v1.0.2.dev0+gfc89945.d20250405)
- [QSIRecon](https://github.com/PennLINC/qsirecon) 1.1.1 (v1.1.1.dev0+gaf43da9.d20250414)

<!--- could add an example spec.json for halfpipe; do we include basic info on the R scripts in Code? --->
<!--- do we include the networks script for connectivity stuff? --->

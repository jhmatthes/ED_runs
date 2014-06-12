These files contain pre- and post-processing code for running the ED2 model.

The SAS_spinup/ directory contains code for running the semi-analytical solution for ED:
1. make_spinup_subdirs.sh constructs a directory and set of sub-directories, copies, and edits
   an example ED2IN file in order to run a spin-up from bare-ground at multiple sites.
2. restart_spinup_histo_lite.sh restarts ED2 runs that have stopped by finding the last
   year and month written to the histo/ directory, editing the ED2IN file appropriately,
   and then re-submitting the job to the cluster (paleon_ed2_geo.sh in each site-level folder 
   where ED2IN is).
3. compile_SAS_runs.R contains the code to loop over the ED2 output analy and histo files 
   to compile the relevant pools and fluxes to construct the .css and .pss files for the semi-
   analytical solution. 
4. phase1a_spinup/ directory contains the output (SAS .css and .pss files) from the 
   compile_SAS_runs.R code.


The ED_process_output/ directory contains some random bits and pieces of plotting code.

The met_drivers/ directory contains the PalEON met drivers for ED. Within this directory the script
process_paleon_met.R will format the PalEON Phase 1a met drivers for input into ED2. 

The phase1a_spininitial/ directory contains the Phase 1a 150-year initial spin-up results (before
the semi-analytical solution is applied).

The phase1a_spinfinish/ directory contains the Phase 1a 150-year settling period run results (after
the semi-analytical solution is applied, using .css and .pss files from SAS_spinup/phase1a_spinup/).

The phase1a_runs/ directory contains the Phase 1a PalEON runs post-spin-up period (where the years
850-2010 are represented as 1850-3010 due to problems running with 3-digit years).

paleon_reg/ contains some preliminary output for regional runs. Mostly just saving these because Mike 
wanted to give some regional data to Paul Duffy for preliminary script-making.

phase1a_env_drivers/ contains the PalEON Phase 1a site-level environmental drivers for the ED runs.

The archive PalEON_old_runz.tar.gz contains old runs (initial tests and things that didn't work out -
mostly junk) that Jackie doesn't want to throw away quite yet.







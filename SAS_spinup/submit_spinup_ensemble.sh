#!bin/bash
#Loop through the set of subdirectories made by make_spinup_subdirs.sh and submit the ensemble 
#of spin-up runs to the geo cluster
#Jaclyn Hatala Matthes, 2/5/14
#jaclyn.hatala.matthes@gmail.com

sites=(PBL PDL PHA PHO PMB PUN)
outdir=/projectnb/cheas/paleon/ED_runs/phase1a_spinup/

if [ ! -d ${outdir} ]
then
    echo "ERROR: you must run make_spinup_subdirs.sh before this!"
fi

for SITE in ${sites[@]}
do
    for n in {1..15}
    do
	rep=$(printf "%02d" ${n})
	pushd ${outdir}$SITE/spin$rep/
	qsub paleon_ed2_geo.sh
	popd
    done
done
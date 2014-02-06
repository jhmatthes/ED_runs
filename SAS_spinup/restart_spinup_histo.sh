#!bin/bash
#This file finds the most recent history file and then 

USER=`whoami`
sites=(PBL PDL PHA PHO PMB PUN)
finalyear=3750
runtype="'INITIAL'"
outdir=/projectnb/cheas/paleon/ED_runs/phase1a_spinup/

if [ ! -d ${outdir} ]
then
    echo "ERROR: you must run make_spinup_subdirs.sh before this!"
fi

#while true
#do
    for SITE in ${sites[@]}
    do
	for n in {1..15}
	do
	    rep=$(printf "%02d" ${n})
	    path=${outdir}${SITE}/spin${rep}/
	    runstat=$(qstat -u $USER | grep ${SITE}spin${rep} | wc -l)

	    lastday=`ls -l -rt ${path}/histo| tail -1 | rev | cut -c15-16 | rev`
	    lastmonth=`ls -l -rt ${path}/histo| tail -1 | rev | cut -c18-19 | rev`
	    lastyear=`ls -l -rt ${path}/histo| tail -1 | rev | cut -c21-24 | rev`
	    echo $lastyear
	    echo $runstat

	    #if run stopped and histo file last year is less than final year
	    if [[ "${runstat}" -eq 0 && "${lastyear}" -lt "${finalyear}" ]]
	    then
		pushd ${path} 
		echo "got here!"

		#edit ED2IN dates for restart
		

		#edit ED2IN params for restart


#	        qsub paleon_ed2_geo.sh
		popd
	    fi
	done
    done
#done

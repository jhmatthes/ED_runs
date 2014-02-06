#!bin/bash
#This file queries the queue to make sure the runs are still going. If not, it edits the ED2IN file
#based on the last histo file printed and restarts a HISTORY run.
#Jaclyn Hatala Matthes, 2/6/14
#jaclyn.hatala.matthes@gmail.com

USER=`whoami`
sites=(PBL PDL PHA PHO PMB PUN)
finalyear=3750
outdir=/projectnb/cheas/paleon/ED_runs/phase1a_spinup/

if [ ! -d ${outdir} ]
then
    echo "ERROR: you must run make_spinup_subdirs.sh before this!"
    echo "ERROR: you must run submit_spinup_ensemble.sh before this to start some runs!"
fi

while true
do
    sleep 300 #only run every 5 minutes

    for SITE in ${sites[@]}
    do
	for n in {1..15}
	do
	    #check if site rep is running
	    rep=$(printf "%02d" ${n})
	    path=${outdir}${SITE}/spin${rep}/
	    runstat=$(qstat -u ${USER} | grep ${SITE}spin${rep} | wc -l)

	    #get dates of last histo file
	    lastday=`ls -l -rt ${path}/histo| tail -1 | rev | cut -c15-16 | rev`
	    lastmonth=`ls -l -rt ${path}/histo| tail -1 | rev | cut -c18-19 | rev`
	    lastyear=`ls -l -rt ${path}/histo| tail -1 | rev | cut -c21-24 | rev`

	    #if run has stopped and last year in histo directory is less than the final run year
	    if [[ "${runstat}" -eq 0 && "${lastyear}" -lt "${finalyear}" ]]
	    then
		pushd ${path} 

		#edit ED2IN dates for restart
		sed -i "s/IYEARA   = .*/IYEARA   = ${lastyear}   ! Year/" ED2IN 
		sed -i "s/IDATEA   = .*/IDATEA   = ${lastday}     ! Day/" ED2IN 
		sed -i "s/IMONTHA  = .*/IMONTHA  = ${lastmonth}     ! Month/" ED2IN 
		sed -i "s/IYEARH   = .*/IYEARH   = ${lastyear}   ! Year/" ED2IN 
		sed -i "s/IDATEH   = .*/IDATEH   = ${lastday}     ! Day/" ED2IN 
		sed -i "s/IMONTHH  = .*/IMONTHH  = ${lastmonth}     ! Month/" ED2IN 

		#edit ED2IN params for history restart
		sed -i 's/IED_INIT_MODE   = .*/IED_INIT_MODE   = 5/' ED2IN
		sed -i "s/RUNTYPE  = .*/RUNTYPE  = 'HISTORY'/" ED2IN

	        qsub paleon_ed2_geo.sh
		popd
	    fi
	done
    done
done

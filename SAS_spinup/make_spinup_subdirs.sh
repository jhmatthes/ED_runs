#!bin/bash
#This file sets up the directories for the ensemble SAS spin-up 
#Jaclyn Hatala Matthes, 2/5/14
#jaclyn.hatala.matthes@gmail.com

sites=(PBL PDL PHA PHO PMB PUN)
startyear=1850
runtype="'INITIAL'"
outdir=/projectnb/cheas/paleon/ED_runs/phase1a_spinup/
cpdir=/projectnb/cheas/paleon/ED_runs/phase1a_runs/

if [ ! -d ${outdir} ]
then
    mkdir -p ${outdir}
fi

for SITE in ${sites[@]}
do
    #make site output directory
    if [ ! -d ${outdir}$SITE/ ]
    then
        mkdir -p ${outdir}/$SITE/
    fi

    #make directories for ensemble runs
    for n in {1..15}
    do
	rep=$(printf "%02d" ${n})
	echo $rep
	if [ ! -d ${outdir}/$SITE/spin$rep/ ]
	then
            mkdir -p ${outdir}/$SITE/spin$rep/
	fi

	#copy files and make directories for each run
	cp ${cpdir}$SITE/ED2IN ${outdir}/$SITE/spin$rep/
	cp ${cpdir}$SITE/paleon_ed2_geo.sh ${outdir}/$SITE/spin$rep/
	cp ${cpdir}$SITE/PL_MET_HEADER ${outdir}/$SITE/spin$rep/
	cp ${cpdir}$SITE/config.xml ${outdir}/$SITE/spin$rep/
	mkdir ${outdir}$SITE/spin$rep/histo
	mkdir ${outdir}$SITE/spin$rep/analy
	
	#edit ED2IN file with correct params & path
	pushd ${outdir}$SITE/spin${rep}/
	ln -s /projectnb/cheas/paleon/ED_runs/T_ED2/ED/build/ed_2.1-opt .
	newbase=${outdir}$SITE/spin${rep}
	oldbase=${cpdir}$SITE
	newpath1="'${outdir}${SITE}/spin${rep}/analy/${SITE}${rep}spin'"
        newpath2="'${outdir}${SITE}/spin${rep}/histo/${SITE}${rep}spin'"
        oldpath1="'${cpdir}$SITE/analy/${SITE}spin'"
	oldpath2="'${cpdir}$SITE/histo/${SITE}spin'"

	sed -i "s/IYEARA   = [1-9][0-9][0-9][0-9]/IYEARA   = $startyear/" ED2IN #change start year value
 	sed -i "s/RUNTYPE  = 'HISTORY'/RUNTYPE  = $runtype/" ED2IN #change run type (INITIAL or HISTORY)
	sed -i "s,$oldpath1,$newpath1,g" ED2IN #change output paths
	sed -i "s,$oldpath2,$newpath2,g" ED2IN #change output paths
	sed -i 's/IED_INIT_MODE   = 5/IED_INIT_MODE   = 0/' ED2IN #change init mode from history to bare ground
	sed -i "s,$oldbase.*,$newbase,g" paleon_ed2_geo.sh #change path in submit script
	sed -i "s,${SITE}spin,${SITE}spin${rep},g" paleon_ed2_geo.sh #change job name
	popd
    done
done







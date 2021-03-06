#!bin/bash
#This file sets up the directories for the runs with the .css and .pss SAS output 
#Jaclyn Hatala Matthes, 2/20/14
#jaclyn.hatala.matthes@gmail.com

sites=(PBL PDL PHA PHO PMB PUN)
startyear=1850
runtype="'INITIAL'"
outdir=/projectnb/cheas/paleon/ED_runs/phase1a_spinfinish/
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
    
	#copy files and make directories for each run
    cp ${cpdir}$SITE/ED2IN ${outdir}/$SITE/
    cp ${cpdir}$SITE/paleon_ed2_geo.sh ${outdir}/$SITE/
    cp ${cpdir}$SITE/PL_MET_HEADER ${outdir}/$SITE/
    cp ${cpdir}$SITE/PalEON_Phase1a.xml ${outdir}/$SITE/
    mkdir ${outdir}$SITE/histo
    mkdir ${outdir}$SITE/analy
    
	#edit ED2IN file with correct params & path
    pushd ${outdir}$SITE/
    ln -s /usr4/spclpgm/jmatthes/ED.r82/ED/build/ed_2.1-opt .
    newbase=${outdir}$SITE
    oldbase=${cpdir}$SITE
    newpath1="'${outdir}${SITE}/analy/${SITE}spin'"
    newpath2="'${outdir}${SITE}/histo/${SITE}spin'"
    oldpath1="'${cpdir}$SITE/analy/${SITE}spin'"
    oldpath2="'${cpdir}$SITE/histo/${SITE}spin'"

    sed -i "s/IYEARA   = [1-9][0-9][0-9][0-9]/IYEARA   = $startyear/" ED2IN #change start year value
    sed -i "s/RUNTYPE  = 'HISTORY'/RUNTYPE  = $runtype/" ED2IN #change run type (INITIAL or HISTORY)
    sed -i "s,$oldpath1,$newpath1,g" ED2IN #change output paths
    sed -i "s,$oldpath2,$newpath2,g" ED2IN #change output paths
    sed -i 's/IED_INIT_MODE   = 0/IED_INIT_MODE   = 3/' ED2IN #change init mode from history to css/pss run
    sed -i "s,$oldbase.*,$newbase,g" paleon_ed2_geo.sh #change path in submit script
    sed -i "s,${SITE}spin,${SITE}spin${rep},g" paleon_ed2_geo.sh #change job name
    popd
done







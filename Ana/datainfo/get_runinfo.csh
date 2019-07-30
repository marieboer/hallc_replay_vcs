#!/bin/csh

set org="/home/cdaq/mboer/hallc_replay_vcs"

echo "arguments: 1) run number 2) events" 
set dir="${org}/REPORT_OUTPUT/COIN/PRODUCTION/"

set file="replay_coin_production"
set loc="${org}/Ana/datainfo/singleruninfo"
set outfile1="${loc}/run_info_${1}.txt"
set outfile2="${loc}/eff_info_${1}.txt"
set outfile3="${loc}/trigger_info_${1}.txt"
set outfile4="${loc}/missref_info_${1}.txt"
rm -f $outfile1 $outfile2 $outfile3 $outfile4
touch $outfile1 $outfile2 $outfile3 $outfile4

echo "run number :"
grep 'Run #' ${dir}/${file}_${1}_${2}.report | awk '{printf $3 " "}' >> $outfile1
grep 'Run #' ${dir}/${file}_${1}_${2}.report | awk '{printf $3 " "}' >> $outfile2
grep 'Run #' ${dir}/${file}_${1}_${2}.report | awk '{printf $3 " "}' >> $outfile3
grep 'Beam energy' ${dir}/${file}_${1}_${2}.report | awk '{printf $4 " "}' >> $outfile1
grep 'Target mass (amu)' ${dir}/${file}_${1}_${2}.report | awk '{printf $5 " "}' >> $outfile1
grep 'HMS P Central' ${dir}/${file}_${1}_${2}.report | awk '{printf $5 " "}' >> $outfile1
grep 'HMS Angle' ${dir}/${file}_${1}_${2}.report | awk '{printf $4 " "}' >> $outfile1
grep 'HMS  Run Length'  ${dir}/${file}_${1}_${2}.report | awk '{printf $5 " "}' >> $outfile1
grep 'SHMS Run Length'  ${dir}/${file}_${1}_${2}.report | awk '{printf $5 " "}' >> $outfile1
grep 'HMS BCM2 Current:' ${dir}/${file}_${1}_${2}.report | awk '{printf $4 " "}' >> $outfile1
grep 'Cut on BCM4 current of' ${dir}/${file}_${1}_${2}.report | awk '{printf $6 " " }' >> $outfile1
grep 'Cut on BCM4 current of' ${dir}/${file}_${1}_${2}.report | awk '{printf $12 " "}' >> $outfile1
grep 'HMS BCM2 Beam Cut Current:' ${dir}/${file}_${1}_${2}.report | awk '{printf $6 " "}' >> $outfile1
grep 'HMS Computer Live Time' ${dir}/${file}_${1}_${2}.report | awk '{printf $6 " "}' >> $outfile1
grep 'HMS 3/4 Trigger Rate' ${dir}/${file}_${1}_${2}.report | awk '{printf $6 " "}' >> $outfile3
grep 'HMS Singles Triggers' ${dir}/${file}_${1}_${2}.report | awk '{printf $5 " "}' >> $outfile3
grep 'Coincidence Triggers' ${dir}/${file}_${1}_${2}.report | awk '{printf $4 " "}' >> $outfile3
grep 'All Triggers' ${dir}/${file}_${1}_${2}.report | awk '{printf $4 " "}' >> $outfile3
grep 'E SING FID TRACK EFFIC' ${dir}/${file}_${1}_${2}.report | awk '{printf $7 " "}' >> $outfile2
grep 'HADRON SING FID TRACK EFFIC' ${dir}/${file}_${1}_${2}.report | awk '{printf $7 " "}' >> $outfile2

grep 'Run #' ${dir}/${file}_${1}_${2}.report | awk '{printf $3 " "}' >> $outfile4
grep 'Missing Ref times:' ${dir}/printed_coin_production_${1}_${2}.report | awk '{printf $5 " " $6 " "}' >> $outfile4


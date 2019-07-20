#!/bin/csh

set dir=`pwd`
cd ..
set where=`pwd`
cd $dir
ln -s UTIL_OL/SCRIPTS/SHELL/run_shms.sh run_shms.sh
ln -s UTIL_OL/SCRIPTS/SHELL/run_hms.sh run_hms.sh
ln -s UTIL_OL/SCRIPTS/SHELL/run_coin_shms.sh run_coin_shms.sh 
ln -s UTIL_OL/SCRIPTS/SHELL/run_coin_shms.sh run_coin_hms.sh 
ln -s UTIL_OL/CHRG_MON/run_charge_counter.csh run_charge_counter.csh
ln -s ${where}/hcana/hcana hcana 
ln -s /cache/hallc/spring17/raw cache 
ln -s /net/cdaq/cdaql1data/cdaq/vcs2019/rootfiles ROOTfiles 
ln -s /net/cdaq/cdaql1data/cdaq/vcs2019/REPORT_OUTPUT REPORT_OUTPUT  
ln -s /net/cdaq/cdaql1data/coda/data/raw raw
ln -s /net/cdaq/cdaql1data/coda/data/raw.copiedtotape raw.copiedtotape
ln -s /net/cdaq/cdaql1data/cdaq/vcs2019/HISTOGRAMS/ HISTOGRAMS
ln -s /net/cdaq/cdaql1data/cdaq/vcs2019/MON_OUTPUT MON_OUTPUT


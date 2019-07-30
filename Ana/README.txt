Code for data analysis from coin replay file and creation of a new root tree

-------------------------------------------------------
Steps before running this code:

1) Make sure to create first the coin_reco...root file and get the report output from the CODA data

2) Run get_runinfo.csh in the datainfo repository: it will create the text files for input data needed by this code

3) Change cuts in ANA_CUTS.input if needed (see description in the file)

4) Then run the code as decribed below

> The will create:
- 1 monitoring pdf file ana_monitor_(index).pdf in the "Results"
- 5 pdf files with deeper monitoring figures in the "Results" directory (if "ana" option is called). files are:
  exclusivity_(index).pdf, cointime_(index).pdf, kinematics_(index).pdf, kinematics2D_(index).pdf, positions_(index).pdf
- a root file of reduced data with 4-vectors in the "files" directory: files/HallCData_(index).root,
- 3 text files timing_(index).dat, missmass2_(index).dat, lumi_(index).dat containing information to build database in the "files" directory.

> you can find all these files in Results/ and sometimes older files moved to Results1 (on /net)

---------------------------------------------------------

- if first time running the code over one set of data, call the option "ana" and check all monitoring plots

- check if the timing cut is good, if not modify it in the header file before using the created root tree

- if you run this code just to produce a new root file (reduce the data), just call the option "reduce". Figures will not be created again 

-------------------------------------------------------

Running HallCdata code with arguments:

example: ./HallCana vcs ana LH2 8585 -f /home/cdaq/vcs2019/hallc_replay_vcs/ROOTfiles/coin_replay_production_8585_10000000.root

Arguments options are:

1) Process:
	  elastic or vcs or pi0 or pi+ or K+ = run setup with electron in SHMS
	  elasticLT or vcsLT or pi0LT or pi+LT or K+LT = run setup with electron in HMS

2) What to do:
   ana (create ana histos and root files), reduce (only root files)

3) which target: 
LH2 or dummy

4) run number: 

5) -f File name(s) or file list: 
   -f or -files = single file or root file list given in header
   -l or -list = name of text file containing the list of root files
   provide only one name if the run is in a single root file, 
   provide a serie or list of file names if the run is split into different root files
   example: -f file1.root (file2.root file3.root) or -l list.txt (list.txt = file1.root file2.root...)
   then provide the list of files txt or the root file(s) names


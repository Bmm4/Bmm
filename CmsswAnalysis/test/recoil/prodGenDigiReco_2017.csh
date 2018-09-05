#! /bin/csh -f


# ----------------------------------------------------------------------
# example submission:
# -------------------
# t3ui01>wpd                                                                                                                                                    # /t3home/ursl/recoil/mcprod/CMSSW_9_3_1/jobs/recoil0
# t3ui01>$BASE/perl/run -t ../../gendigireco-180904.tar.gz -m local -c $BASE/CmsswAnalysis/test/recoil/prodGenDigiReco_2017.csh -r 'PFNS srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat%STORAGE3 /store/user/ursl/recoil/mcprod2017/aodsim/recoil0%SITE T3_CH_PSI' recoil0_gensim-40000.py
# t3ui01>$BASE/perl/run -t ../../gendigireco-180904.tar.gz -q all.q -c $BASE/CmsswAnalysis/test/recoil/prodGenDigiReco_2017.csh -r 'PFNS srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat%STORAGE3 /store/user/ursl/recoil/mcprod2017/aodsim/recoil0%SITE T3_CH_PSI' recoil0_gensim-40000.py
#
# Note: this script uses the py file with which it is submitted for GENERATION,
#       it uses step1.py and step2.py which should be in the digireco.tar.gz file.
# ----------------------------------------------------------------------
setenv GENRELEASE CMSSW_9_3_1
setenv RECRELEASE CMSSW_9_4_0_patch1

setenv SCRAM_ARCH
setenv SRMCP

setenv JOB
setenv STORAGE3
setenv FILE1    file:./$JOB.root
setenv FILE2    file:./$JOB:s/gensim/step1/.root
setenv FILE3    $JOB:s/gensim/step2/.root
setenv PFNS
setenv SITE



echo "========================"
echo "====> SGE  wrapper <===="
echo "========================"

echo "--> Running SGE 2017 gen-digi-reco job wrapper"
echo $JOB
echo $FILE1
echo $FILE2
echo $FILE3

# ----------------------------------------------------------------------
# -- The Basics
# ----------------------------------------------------------------------
echo "--> Environment"
date
hostname
uname -a
#df -kl
limit coredumpsize 0

echo $VO_CMS_SW_DIR
ls -l $VO_CMS_SW_DIR
source $VO_CMS_SW_DIR/cmsset_default.csh
echo "-> which edg-gridftp-ls"
which edg-gridftp-ls
echo "-> which globus-url-copy"
which globus-url-copy
echo "-> which srmcp"
which srmcp

pwd
echo "--> End of env testing"

# BATCH START


# ----------------------------------------------------------------------
# -- Setup CMSSW for GENSIM
# ----------------------------------------------------------------------
echo "--> Setup CMSSW for GENSIM"
pwd
date
cmsrel $GENRELEASE
cd $GENRELEASE
eval `scramv1 runtime -csh`
pwd

#echo "--> Extract tar file"
#date
tar zxvf ../$JOB.tar.gz
cd src
mv ../../$JOB.py .


# ----------------------------------------------------------------------
# -- Run cmsRun for GENSIM
# ----------------------------------------------------------------------
echo "--> Run cmsRun for GENSIM"
pwd
date
which cmsRun
echo "cmsRun $JOB.py "
cmsRun $JOB.py |& tee $JOB.log
date
pwd
ls -rtl


# ----------------------------------------------------------------------
# -- Setup CMSSW for DIGIRECO
# ----------------------------------------------------------------------
echo "--> Setup CMSSW for DIGIRECO"
# -- note the two spaces to avoid 'run' inserting its value for the variable
setenv  SCRAM_ARCH  slc6_amd64_gcc493
cd ../../
pwd
date
cmsrel $RECRELEASE
cd $RECRELEASE
eval `scramv1 runtime -csh`
pwd

#echo "--> Extract tar file"
#date
tar zxvf ../$JOB.tar.gz step1.py
tar zxvf ../$JOB.tar.gz step2.py
cd src
cp ../step1.py .
cp ../step2.py .
pwd
ls -rtl
mv ../../$GENRELEASE/src/$JOB.root .
ls -rtl

# ----------------------------------------------------------------------
# -- Run cmsRun for DIGIRECO
# ----------------------------------------------------------------------
echo "--> Run cmsRun for DIGIRECO"
pwd
date
which cmsRun
echo "cmsRun step1.py "
cmsRun step1.py |& tee step1.log
date
pwd
cp ../step2.py .
ls -rtl
echo "cmsRun step2.py "
cmsRun step2.py |& tee step2.log
date
pwd
ls -rtl
echo "check parent directory"
ls -rtl ..

# ----------------------------------------------------------------------
# -- Save Output to SE
# ----------------------------------------------------------------------
echo "--> Save output to SE: $PFNS/$STORAGE3/$FILE3"
echo " job   rootfiles: $FILE3"

echo lcg-del -b -D srmv2 -l  "$PFNS/$STORAGE3/$FILE3"
lcg-del -b -D srmv2 -l "$PFNS/$STORAGE3/$FILE3"
# -- switch to data_replica.py
#ls `pwd`/$FILE3 > dr.list
#echo "--> cat dr.list: "
#cat dr.list
#echo "--> AM running data_replica.py: /mnt/t3nfs01/data01/swshare/psit3/bin/data_replica.py --from-site LOCAL --to-site $SITE dr.list $STORAGE3 "
#/mnt/t3nfs01/data01/swshare/psit3/bin/data_replica.py --from-site LOCAL --to-site $SITE dr.list "$STORAGE3"

# -- switch to xrdcp
xrdcp $FILE3 root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/$STORAGE3

echo "--> lcg-ls : $PFNS/$STORAGE3/$FILE3"
echo lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE3/$FILE3"
lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE3/$FILE3"

date

# BATCH END

echo "run: This is the end, my friend"

#! /bin/csh -f


# ----------------------------------------------------------------------
# example submission:
# -------------------
#
# $BMMBASE/perl/run -t $BMMBASE/../../160707.tar.gz -c $BMMBASE/CmsswAnalysis/test/bmm4/gridNoComp.csh
# -r 'STORAGE1 /store/user/ursl/bmm4/cmsRun/v01/bmmSingleMuon2016C%XRD root://t3se01.psi.ch:1094%SITE T3_CH_PSI'
# -m grid -D cern.ch bmm-prompt-Run2016C-SingleMuon_Run2016C-PromptReco-v01-027*.py >& submit.log
# ----------------------------------------------------------------------

setenv CMSSW
setenv SCRAM_ARCH
setenv SRMCP

setenv JOB
setenv FILE1    $JOB.root
setenv STORAGE1
setenv XRD
setenv SITE


echo "========================"
echo "====> grid wrapper <===="
echo "========================"

echo "--> Running grid job wrapper"

# ----------------------------------------------------------------------
# -- The Basics
# ----------------------------------------------------------------------
echo "--> Environment"
date
echo "--> hostname"
hostname
echo "--> cat /proc/cpuinfo"
cat /proc/cpuinfo
echo "--> uname -a"
uname -a
echo "--> df -kl"
df -kl
echo "--> printenv"
printenv
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
# -- Setup CMSSW
# ----------------------------------------------------------------------
echo "--> Setup CMSSW"
pwd
date
cmsrel $CMSSW
cd $CMSSW
eval `scramv1 runtime -csh`
pwd

echo "--> Extract tar file"
date
tar zxf ../$JOB.tar.gz
cd src
mv ../../$JOB.py .
mv ../../data_replica.py .
chmod 755 data_replica.py


# ----------------------------------------------------------------------
# -- Run cmsRun
# ----------------------------------------------------------------------
echo "--> Run cmsRun"
pwd
date
which cmsRun
echo "cmsRun $JOB.py "
cmsRun $JOB.py |& tee $JOB.log
date
pwd
ls -rtl

setenv LD_LIBRARY_PATH /lib64:${LD_LIBRARY_PATH}

setenv ROOTFILE `ls *.root`

# ----------------------------------------------------------------------
# -- Save Output to SE
# ----------------------------------------------------------------------
echo "--> Save output to SE: $XRD/$STORAGE1/$FILE1"
echo " local rootfile: $ROOTFILE"
echo " job   rootfile: $FILE1"

echo "--> AM running xrdcp: " xrdcp -f -d 3 `pwd`/$FILE1 $XRD/$STORAGE1/$FILE1
xrdcp -f -d 1 `pwd`/$FILE1 $XRD/$STORAGE1/$FILE1

echo xrdfs root://t3se01.psi.ch:1094 stat $STORAGE1/$FILE1
xrdfs root://t3se01.psi.ch:1094 stat $STORAGE1/$FILE1

date

# BATCH END

echo "run: This is the end, my friend"

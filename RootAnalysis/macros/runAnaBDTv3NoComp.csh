#! /bin/csh -f

setenv CMSSW
setenv SCRAM_ARCH
setenv SRMCP

setenv JOB
setenv EXECUTABLE
setenv FILE1 $JOB.root
setenv STORAGE1
setenv PFNS
setenv SITE

echo "========================"
echo "====> SGE wrapper <===="
echo "========================"

echo "--> Running grid job wrapper"

# ----------------------------------------------------------------------
# -- The Basics
# ----------------------------------------------------------------------
echo "--> Environment:"
date
hostname
#cat /proc/cpuinfo
uname -a
echo "--> df -kl:"
df -kl

echo $VO_CMS_SW_DIR
#ls -l $VO_CMS_SW_DIR
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
# -- Setup TreeReader
# ----------------------------------------------------------------------
echo "--> Setup TreeReader"
pwd
date
cmsrel $CMSSW

cd $CMSSW/
eval `scramv1 runtime -csh`
pwd

echo "--> Extract tar file"
date
tar zxf ../$JOB.tar.gz
mv ../$JOB src/Bmm/RootAnalysis/macros

# ----------------------------------------------------------------------
# -- Run Treereader
# ----------------------------------------------------------------------

echo "--> Run Treereader"
date
cd src/Bmm/RootAnalysis/macros/
pwd
ls -rtl

echo "OUTPUTFILE: $FILE1"
echo "EXECUTABLE: $EXECUTABLE $JOB  |& tee $JOB.log"
$EXECUTABLE  $JOB  |& tee $JOB.log
date

touch /scratch/ursl/bmm4/bdt/*.root
/bin/ls -l /scratch/ursl/bmm4/bdt/*.root

date

# BATCH END


# -- cleanup
#/bin/rm -rf /tmp/

echo "run: This is the end, my friend"

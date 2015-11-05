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
# -- use special root version
#setenv ROOTSYS /shome/naegelic/root
#setenv LD_LIBRARY_PATH ${ROOTSYS}/lib:${LD_LIBRARY_PATH}
#setenv PATH ${ROOTSYS}/bin:${PATH}
#echo "--> Unsetting SCRAM_ARC"
#unsetenv SCRAM_ARCH
#echo "--> Sourcing Christophs script"
#source /shome/naegelic/root/bin/thisroot.csh
#echo "--> Sourced Christophs script"

cd $CMSSW/
eval `scramv1 runtime -csh`
#?? setenv LD_LIBRARY_PATH /swshare/glite/d-cache/dcap/lib/:${LD_LIBRARY_PATH}
pwd

echo "--> Extract tar file"
date
tar zxf ../$JOB.tar.gz
mv ../$JOB src/Bmm/RootAnalysis/inclb

# ----------------------------------------------------------------------
# -- Run Treereader
# ----------------------------------------------------------------------

echo "--> Run Treereader"
date
cd src/Bmm/RootAnalysis/inclb
pwd

echo "$EXECUTABLE -c $JOB  -o $FILE1 |& tee $JOB.log"
$EXECUTABLE -c $JOB -o $FILE1 |& tee $JOB.log
date

# ----------------------------------------------------------------------
# -- Save Output to NFS, not the SE
# ----------------------------------------------------------------------

echo "--> Save output to SE: $PFNS/$STORAGE1/$FILE1"
echo " job   rootfile: $FILE1"

echo lcg-del -b -D srmv2 -l  "$PFNS/$STORAGE1/$FILE1"
lcg-del -b -D srmv2 -l "$PFNS/$STORAGE1/$FILE1"
# -- switch to data_replica.py
ls `pwd`/$FILE1 > dr.list
echo "--> cat dr.list: " 
cat dr.list
echo "--> AM running data_replica.py: " 
/swshare/psit3/bin/data_replica.py --from-site LOCAL --to-site $SITE dr.list "$STORAGE1"

echo "--> lcg-ls : $PFNS/$STORAGE1/$FILE1" 
echo lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE1/$FILE1"
lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE1/$FILE1"

date

# BATCH END


# -- cleanup

echo "run: This is the end, my friend"

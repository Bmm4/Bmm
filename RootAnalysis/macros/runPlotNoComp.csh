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
ls -rtl
ls -rtl results


# ----------------------------------------------------------------------
# -- Save Output to SE
# ----------------------------------------------------------------------
echo "--> Save output to SE: $PFNS/$STORAGE1/$FILE1"

# -- multi file output
cd results
foreach x ( *.root )
  echo $x
  echo "rootfile: $x"
  echo lcg-del -b -D srmv2 -l  "$PFNS/$STORAGE1/$x"
  lcg-del -b -D srmv2 -l "$PFNS/$STORAGE1/$x"
  echo "--> AM running xrdcp $x root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/$STORAGE1/$x "
  xrdcp $x root://t3dcachedb.psi.ch:1094///pnfs/psi.ch/cms/trivcat/$STORAGE1/$x
  echo "--> lcg-ls : $PFNS/$STORAGE1/$x"
  echo lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE1/$x"
  lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE1/$x"
  xrdfs t3dcachedb03.psi.ch ls -l -u //pnfs/psi.ch/cms/trivcat/$STORAGE1/ | /bin/grep $FILE1
end

date

# BATCH END


# -- cleanup
#/bin/rm -rf /tmp/

echo "run: This is the end, my friend"

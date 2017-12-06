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
ls -rtl dataset/weights
cp dataset/weights/TMVA-$JOB-Events0_BDT.weights.xml weights/TMVA-$JOB-Events0_BDT.weights.xml
cp dataset/weights/TMVA-$JOB-Events1_BDT.weights.xml weights/TMVA-$JOB-Events1_BDT.weights.xml
cp dataset/weights/TMVA-$JOB-Events2_BDT.weights.xml weights/TMVA-$JOB-Events2_BDT.weights.xml
date
rm -f results/baseCuts-$JOB.cuts
cat cuts/baseCuts.nobdt.cuts append-basecuts.txt >> results/baseCuts-$JOB.cuts
ls -rtl results
bin/runPlot -p results  -m bdtopt -c baseCuts-$JOB.cuts -y 2016BF -f plotResults.2016GHse.files |& tee -a $JOB.log
date
ls -rtl
bin/runPlot -p results  -m bdtopt -c baseCuts-$JOB.cuts -y 2016BF -f plotResults.2016BFse.files |& tee -a $JOB.log
date
ls -rtl
bin/runPlot -p overlays -m bdtopt -c baseCuts-$JOB.cuts -y 2016BF -f plotResults.2016BFse.files |& tee -a $JOB.log
date
ls -rtl
bin/runPlot -p overlays -m bdtopt -c baseCuts-$JOB.cuts -y 2016BF -f plotResults.2016GHse.files |& tee -a $JOB.log
date
ls -rtl

# ----------------------------------------------------------------------
# -- Save Output to SE
# ----------------------------------------------------------------------
echo "--> Save output to SE: $PFNS/$STORAGE1/$FILE1"

# -- multi file output
cd results
foreach x ( *.root *.tex)
  echo $x
  echo "rootfile: $x"
  echo lcg-del -b -D srmv2 -l  "$PFNS/$STORAGE1/$x"
  lcg-del -b -D srmv2 -l "$PFNS/$STORAGE1/$x"
  echo "--> AM running xrdcp $x root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/$STORAGE1/$x "
  xrdcp $x root://t3se01.psi.ch:1094/pnfs/psi.ch/cms/trivcat/$STORAGE1/$x
  echo "--> lcg-ls : $PFNS/$STORAGE1/$x"
  echo lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE1/$x"
  lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE1/$x"
end

date

# BATCH END


# -- cleanup
#/bin/rm -rf /tmp/

echo "run: This is the end, my friend"

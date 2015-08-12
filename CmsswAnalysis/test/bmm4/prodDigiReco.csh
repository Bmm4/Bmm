#! /bin/csh -f


# ----------------------------------------------------------------------
# example submission: 
# -------------------
# $BMMBASE/perl/run -t ../../../digireco.tar.gz -m batch -c ../../prodDigiReco.csh -r 'PFNS srm://t3se01.psi.ch:8443/srm/managerv2\?SFN=/pnfs/psi.ch/cms/trivcat%STORAGE1 /store/user/ursl/bmm4/prod/gen/Bs2JpsiPhi_EtaPtFilter%STORAGE2 /store/user/ursl/bmm4/prod/aodsim/Bs2JpsiPhi_EtaPtFilter%SITE T3_CH_PSI'  PYTHIA8_Bs2JpsiPhi_EtaPtFilter_CUEP8M1_13TeV_step1-70000
# ----------------------------------------------------------------------

setenv CMSSW       
setenv SCRAM_ARCH  
setenv SRMCP       

setenv JOB      
setenv STORAGE1 
setenv STORAGE2 
setenv FILE1    $STORAGE1/$JOB.root
setenv FILE2    $JOB:s/step1/step2/.root
setenv FILE3    $JOB:s/step1/step3/.root
setenv PFNS     
setenv SITE     

echo "========================"
echo "====> SGE  wrapper <===="
echo "========================"

echo "--> Running SGE digi-reco job wrapper"
echo $JOB
echo $FILE1 
echo $FILE2
echo $FILE3 
echo $STORAGE1
echo $STORAGE2
# ----------------------------------------------------------------------
# -- The Basics
# ----------------------------------------------------------------------
echo "--> Environment"
date
hostname
uname -a
df -kl 

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

#echo "--> Extract tar file"
#date
tar zxf ../$JOB.tar.gz
cd src
mv ../step2.py .
mv ../step3.py .


# ----------------------------------------------------------------------
# -- Run cmsRun
# ----------------------------------------------------------------------
echo "--> Run cmsRun"
pwd
date
which cmsRun
echo "cmsRun step2.py "
cmsRun step2.py |& tee step2.log
date
pwd
ls -rtl 
echo "cmsRun step3.py "
cmsRun step3.py |& tee step3.log
date
pwd
ls -rtl 

exit(0)

#setenv ROOTFILE `ls *.root | /bin/grep step3`
setenv ROOTFILE FILE3


# ----------------------------------------------------------------------
# -- Save Output to SE
# ----------------------------------------------------------------------
echo "--> Save output to SE: $PFNS/$STORAGE2/$FILE3"
echo " local rootfile: $ROOTFILE"
echo " job   rootfiles: $FILE1, $FILE2, $FILE3"

echo lcg-del -b -D srmv2 -l  "$PFNS/$STORAGE2/$FILE3"
lcg-del -b -D srmv2 -l "$PFNS/$STORAGE2/$FILE3"
# -- switch to data_replica.py
ls `pwd`/$FILE3 > dr.list
echo "--> cat dr.list: " 
cat dr.list
echo "--> AM running data_replica.py: " 
/swshare/psit3/bin/data_replica.py --from-site LOCAL --to-site $SITE dr.list "$STORAGE2"

echo "--> lcg-ls : $PFNS/$STORAGE2/$FILE3" 
echo lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE2/$FILE3"
lcg-ls -b -D srmv2 -l  "$PFNS/$STORAGE2/$FILE3"

date

# BATCH END

echo "run: This is the end, my friend"

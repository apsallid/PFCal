#!/bin/bash

echo "The script starts now."

echo "System: "
uname -a

source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh

source /afs/cern.ch/work/a/apsallid/CMS/Geant4/PFCal/PFCalEE/g4env.sh

export EOSMAINPATH=/eos/cms/store/group/phys_b2g/apsallid/PFCal

export PWD=`pwd`

cp /afs/cern.ch/work/a/apsallid/CMS/Geant4/PFCal/PFCalEE/runparallel.mac .

sed -e "s/PARTICLEGUN/TYPEGUN/g" runparallel.mac > voodoo
sed -e "s/ENERGYOFGUN/TYPEENE/g" voodoo > voodoo1
sed -e "s/SEEDNUM/THESEEDNUM/g" voodoo1 > voodoo2

mv voodoo2 run.mac 

PFCalEE run.mac THEDETECTORCONFIG 0

mv PFcal.root OUTPUT  

#cp OUTPUT /afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/Parallel/DetectorConfigurations/Config_THEDETECTORCONFIG/TYPEGUN/TYPEENE/results 

#/afs/cern.ch/project/eos/installation/cms/bin/eos.select cp $PWD/OUTPUT $EOSMAINPATH/DetectorConfigurations/Config_THEDETECTORCONFIG/TYPEGUN/TYPEENE/results/OUTPUT 

/afs/cern.ch/project/eos/installation/cms/bin/eos.select cp $PWD/OUTPUT $EOSMAINPATH/DetectorConfigurations/NospreadNoupstream/Config_THEDETECTORCONFIG/TYPEGUN/TYPEENE/results/OUTPUT



 

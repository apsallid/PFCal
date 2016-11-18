#!/bin/tcsh

# Detector Configurations
# For details: https://twiki.cern.ch/twiki/bin/view/Main/SiWECALAnalysisDataStructure
#setenv detectorconfigurations "1 2 3 4 5 6 7 8 9 10 11"
setenv detectorconfigurations "12"
# Particle 
#setenv suddirineoslistguns "e+ pi+ mu+"
#setenv suddirineoslistguns "mu-"
setenv suddirineoslistguns "e-"
# Energies in GeV
setenv suddirineoslistenergies "2 3 5 8 10 15 30 50 80 100 120 150 180 200 250 300 400 500"
#setenv suddirineoslistenergies "15 30 50"
#setenv suddirineoslistenergies "15"
#setenv suddirineoslistenergies "15 30 50 80 100 150"
setenv eospath "/eos/cms/store/group/phys_b2g/apsallid/PFCal/DetectorConfigurations/NospreadNoupstream"
# Starting the loop through all test beam setups, particles and energies
foreach detconf ($detectorconfigurations)
echo "===================================================================================="
echo "Detector configuration $detconf"
foreach subdirineosguns  ($suddirineoslistguns)
echo "===================================================================================="
foreach subdirineosenergies  ($suddirineoslistenergies)
echo "------------------------"
echo "${subdirineosguns} ${subdirineosenergies} GeV"

#setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/Parallel/DetectorConfigurations/Config_${detconf}/$subdirineosguns/$subdirineosenergies"
setenv workpath "$eospath/Config_${detconf}/$subdirineosguns/$subdirineosenergies"

setenv output PFCal_${detconf}_${subdirineosguns}_${subdirineosenergies}GeV.root

#mkdir -p /tmp/apsallid/Config_${detconf}/$subdirineosguns/$subdirineosenergies

#/afs/cern.ch/project/eos/installation/cms/bin/eos.select cp -r ${workpath}/ /tmp/apsallid/Config_${detconf}/$subdirineosguns/$subdirineosenergies/

/afs/cern.ch/project/eos/installation/cms/bin/eos.select cp -r ${workpath}/ /tmp/apsallid/Configs/


#hadd $output ${workpath}/SiWEcal_${detconf}_${subdirineosguns}_${subdirineosenergies}GeV_*.root
#echo "/tmp/apsallid/Config_${detconf}/$subdirineosguns/$subdirineosenergies/SiWEcal_${detconf}_${subdirineosguns}_${subdirineosenergies}GeV_*.root"
 
hadd /tmp/apsallid/$output /tmp/apsallid/Configs/results/PFCal_${detconf}_${subdirineosguns}_${subdirineosenergies}GeV_*.root

/afs/cern.ch/project/eos/installation/cms/bin/eos.select rm $eospath/Config_${detconf}/$output

/afs/cern.ch/project/eos/installation/cms/bin/eos.select cp /tmp/apsallid/$output $eospath/Config_${detconf}/$output

end 

end

end


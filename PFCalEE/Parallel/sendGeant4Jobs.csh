#!/bin/tcsh

# Detector Configurations
# For details: https://twiki.cern.ch/twiki/bin/view/Main/SiWECALAnalysisDataStructure
#setenv detectorconfigurations "1 2 3 4 5 6 7 8 9 10 11"
setenv detectorconfigurations "12"
# Particle 
#setenv suddirineoslistguns "e+ pi+ mu+"
setenv suddirineoslistguns "e-"
# Energies in GeV
setenv suddirineoslistenergies "2 3 5 8 10 15 30 50 80 100 120 150 180 200 250 300 400 500"
# Starting the loop through all test beam setups, particles and energies
foreach detconf ($detectorconfigurations)
echo "===================================================================================="
echo "Detector configuration $detconf"
foreach subdirineosguns  ($suddirineoslistguns)
echo "===================================================================================="
foreach subdirineosenergies  ($suddirineoslistenergies)
echo "------------------------"
echo "${subdirineosguns} ${subdirineosenergies} GeV"

setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/Geant4/PFCal/PFCalEE/Parallel/DetectorConfigurations/NospreadNoupstream/Config_${detconf}/$subdirineosguns/$subdirineosenergies/jobs"

setenv runumberslist ` ls -ltr ${workpath} | grep .job |  awk '{print $9}' `

foreach run  ($runumberslist)

#echo ${run}
chmod 755 ${workpath}/${run}

echo "Sending ${run}"
bsub -q 1nd -o /tmp/junk ${workpath}/${run}

end

end 

end

end




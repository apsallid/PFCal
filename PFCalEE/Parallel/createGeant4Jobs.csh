#!/bin/tcsh

# Detector Configurations
# For details look in DetectorConstruction.hh
# v_HGCALEE_v5=12
setenv detectorconfigurations "12"
# Particle 
#setenv suddirineoslistguns "e+ pi+ mu+"
#setenv suddirineoslistguns "mu- e-"
setenv suddirineoslistguns "e-"
# Energies in GeV
setenv suddirineoslistenergies "2 3 5 8 10 15 30 50 80 100 120 150 180 200 250 300 400 500"
# Delete some auxilliary files just in case
rm voodoo voodoo1 voodoo2 voodoo3 voodoo4 voodoo5
# Clear eos space for new files
setenv eospath "/eos/cms/store/group/phys_b2g/apsallid/PFCal/DetectorConfigurations/NospreadNoupstream"
/afs/cern.ch/project/eos/installation/cms/bin/eos.select rm -r $eospath
# Starting the loop through all test beam setups, particles and energies
foreach detconf ($detectorconfigurations)
echo "===================================================================================="
echo "Detector configuration $detconf" 
foreach subdirineosguns  ($suddirineoslistguns)
foreach subdirineosenergies  ($suddirineoslistenergies)
echo "------------------------"
echo "${subdirineosguns} ${subdirineosenergies} GeV"

foreach num (`seq 1 100`)

set num1=`expr ${num} \* 1257 + 158725983 `

cat setupGeant4.sh > voodoo

setenv file PFCal_${detconf}_${subdirineosguns}_${subdirineosenergies}GeV_${num}.root  
setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/Geant4/PFCal/PFCalEE/Parallel/DetectorConfigurations/NospreadNoupstream/Config_${detconf}/$subdirineosguns/$subdirineosenergies"

if ( $num == 1 ) then
rm -rf $workpath/jobs $workpath/logfiles $workpath/results
mkdir -p $workpath/jobs $workpath/logfiles $workpath/results
chmod 755 $workpath/jobs $workpath/logfiles $workpath/results

/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir $eospath 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir $eospath/Config_${detconf} 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir $eospath/Config_${detconf}/$subdirineosguns 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir $eospath/Config_${detconf}/$subdirineosguns/$subdirineosenergies 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir $eospath/Config_${detconf}/$subdirineosguns/$subdirineosenergies/results 

/afs/cern.ch/project/eos/installation/cms/bin/eos.select chmod 755 $eospath 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select chmod 755 $eospath/Config_${detconf} 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select chmod 755 $eospath/Config_${detconf}/$subdirineosguns 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select chmod 755 $eospath/Config_${detconf}/$subdirineosguns/$subdirineosenergies 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select chmod 755 $eospath/Config_${detconf}/$subdirineosguns/$subdirineosenergies/results 
endif

sed -e "s/OUTPUT/$file/g" voodoo > voodoo1
sed -e "s/TYPEGUN/$subdirineosguns/g" voodoo1 > voodoo2
sed -e "s/TYPEENE/$subdirineosenergies/g" voodoo2 > voodoo3
sed -e "s/THESEEDNUM/${num1}/g" voodoo3 > voodoo4
sed -e "s/THEDETECTORCONFIG/${detconf}/g" voodoo4 > voodoo5

mv voodoo5 $workpath/jobs/geant_${num}.job
chmod 755 $workpath/jobs/geant_${num}.job

echo geant_${num}.job

rm voodoo 
rm voodoo1
rm voodoo2
rm voodoo3
rm voodoo4

end

end

end

end
~


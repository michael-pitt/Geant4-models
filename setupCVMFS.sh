#The following lines will check your OS and setup correct LCG enviroment
if [[ `lsb_release -d | awk -F"\t" '{print $2}'` =~ "CentOS" ]]; then
  echo running with `lsb_release -d | awk -F"\t" '{print $2}'`, SETUP LCG_96/x86_64-centos7-gcc8-opt
  source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-centos7-gcc8-opt/setup.sh
else
  echo running with `lsb_release -d | awk -F"\t" '{print $2}'`, SETUP LCG_96/x86_64-slc6-gcc8-opt
  source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh
fi


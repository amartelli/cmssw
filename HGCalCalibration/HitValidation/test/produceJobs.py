#!/usr/bin/python                                                                                                                                                                                                                 

#nRUN = 245192
#nRUN = "MinBias"
## change run number

import os
currentDir = os.getcwd()
print 'currentDir = %s' %(currentDir)

##settings
nRUN = ["Pho22_E60"]
nRUN.append("Pho22_E5")
nRUN.append("Pho22_E10")
nRUN.append("Pho22_E30")
nRUN.append("Pho22_E100")
nRUN.append("Pho22_E200")
nRUN.append("Pion211_E5")
nRUN.append("Pion211_E10")
nRUN.append("Pion211_E30")
nRUN.append("Pion211_E60")
nRUN.append("Pion211_E100")
nRUN.append("Pion211_E200")

print 'number of jobs = %d' %(len(nRUN))

outLancia = open('%s/lancia.sh' %(currentDir),'w');

cfgTemp = "partGun_NTUP_template.py"

for num in nRUN:
    print num;
    subDir = currentDir+"/JOB_%s" %(num);
    os.system('mkdir %s' %(subDir));
    print 'subDir = %s  >>>> DID YOU CREATE THIS?' %(subDir);

    inFileList = currentDir+"/fileList/list%s.txt" %(num);
    print 'inFileList = %s' %(inFileList);

    cfgFile = "partGun_NTUP_fromtemp_%s.py" %(num);
    os.system("cp %s %s/%s" %(cfgTemp, subDir, cfgFile));
    os.system("sed -i s~INPUTFILELIST~%s~g %s" %(currentDir+"/fileList/list"+num+".txt", subDir+"/"+cfgFile));
    os.system("sed -i s~OUTFILE~%s~g %s" %("Calib_"+num, subDir+"/"+cfgFile));

    outScript = open('%s/bjob.sh' %(subDir),'w');
    outScript.write('#!/bin/bash \n');
    outScript.write('cd %s \n' %(subDir));
    outScript.write('export SCRAM_ARCH=slc6_amd64_gcc530 \n');
    outScript.write('eval `scramv1 ru -sh` \n');
    #outScript.write('cd - \n');
    #outScript.write('cmsMkdir %s \n' %(outFolder));
    outScript.write('pwd \n ')
    outScript.write('cmsRun  %s \n' %(cfgFile) );
    #outScript.write('cmsStage coll_timing_%d_%d.root %s/' %(num, num+plusNum, outFolder));
    outScript.write( '\n ' );
    os.system('chmod 777 %s/bjob.sh' %(subDir));
    outLancia.write(' bsub -cwd %s -q cmscaf1nd %s/bjob.sh \n' %(subDir, subDir));


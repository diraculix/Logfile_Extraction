# -*- coding: utf-8 -*-
import os

raw = r'N:/fs4-HPRT/HPRT-Data/ONGOING_PROJECTS/AutoPatSpecQA/Logfiles_Spotshape_QA/raw'
conv = r'N:/fs4-HPRT/HPRT-Data/ONGOING_PROJECTS/AutoPatSpecQA/Logfiles_Spotshape_QA/converted'

if not os.path.isdir(raw) or not os.path.isdir(conv):
    raw = r'C:\Users\lukas\Documents\OncoRay HPRT\Logfile_Extraction_mobile\raw'
    conv = r'C:\Users\lukas\Documents\OncoRay HPRT\Logfile_Extraction_mobile\converted'

for logfilezip in os.listdir(raw):
    date = logfilezip.split('_')[0]
    beam = logfilezip.split('_')[1]
    inputdir = os.path.join(raw, logfilezip)
    outputdir = os.path.join(conv, date, beam)
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)

    command = ('java -classpath sites.configs.LLN.LAB-continuous-SNAPSHOT-deploy.jar;//C:/Users/lukas/Documents/OncoRay HPRT/Logfile_Extraction_mobile/data-recorder-proc-8.5.1.2-deploy.jar com.iba.icomp.core.launch.Launcher config/bms/dataRecorder/container.xml -Dconfig.recorder.mode=converter -Dbms.datarecorder.converter.file='+inputdir+' -Dbms.datarecorder.converter.outputdir='+outputdir)
    os.system(command)

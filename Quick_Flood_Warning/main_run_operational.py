# -*- coding: utf-8 -*-
"""
Created on Sat Jun 12 18:16:18 2021

@author: antonioh
"""

import datetime as dt
from datetime import timedelta
import os
import shutil
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from bs4 import BeautifulSoup
import sys
#sys.path.insert(1,"C:\Users\moritzw\OneDrive - SPC\Documents\Projects\ECIKS\Quick_Flood_Warning\codes")
import step1_download_NCEP as step1
import step2_download_CMEMS as step2
import step3_gen_tide_TPOX8 as step3
import step4_make_wave_forcing as step4
import step5_ingest as step5
#####################################


def list_available_runs(url):
    session = requests.Session()
    retry = Retry(connect=5, backoff_factor=1)
    adapter = HTTPAdapter(max_retries=retry)
    session.mount(url, adapter)

    try:
        req = session.get(url).text
        session.close() 
        soup = BeautifulSoup(req, 'html.parser')     
        x = (soup.find_all('a'))
        runs = []
        for i in x:
            file_name = i.extract().get_text()
            runs.append(int(file_name))
    except:
        
        runs = []
        print('Keep working on making the downloading process more robust')
      
    return(runs)

###############################################################################
if __name__ == "__main__":
    
    now = dt.datetime.now()
    current_time = now.strftime("%H:%M:%S")
    current_date = now.strftime("%Y/%m/%d")
    print("Current Date = ",current_date)
    print("Current Time = ", current_time)
    file = open("last_run_start.txt", "w")
    a = file.write("Last run started on "+ current_date +" at "+ current_time)
    file.close()
	
    #- Find the latest available run in nomads.ncep.noaa.gov
    now = dt.datetime.utcnow()
    url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gefs_wave_0p25.pl?dir=%2Fgefs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
    runs = list_available_runs(url)
    if len(runs)==0:
        now = dt.datetime.utcnow()-timedelta(1)
        url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gefs_wave_0p25.pl?dir=%2Fgefs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
        runs = list_available_runs(url)
    
    #- Define the run to be used
    runs=sorted(runs)
    now = now.replace(hour=runs[-1],minute=0,second=0,microsecond=0)
    #now = dt.datetime(2023,4,29,0)

    step1.download_NCEP(now)
    step2.download_CNEMS(now)
    step3.gen_tide(now)
    step4.make_waves(now)
    step5.ingest_products()
	
    now = dt.datetime.now()
    current_time = now.strftime("%H:%M:%S")
    current_date = now.strftime("%Y/%m/%d")
    print("Current Date = ",current_date)
    print("Current Time = ", current_time)
    file = open("last_run_end.txt", "w")
    a = file.write("Last run completed on "+ current_date +" at "+ current_time)
    file.close()
	
    #import pysftp
    #rootDir='D:\\ECIKS\\Quick_Flood_Warning\\results\\'
    #Hostname = "192.168.53.43"
    #Username = "ftpanuj"
    #Password = "0cean2022"
    #cnopts = pysftp.CnOpts()
    #cnopts.hostkeys = None
    #for file in os.listdir(rootDir):
    #    with pysftp.Connection(host=Hostname, username=Username, password=Password,cnopts=cnopts) as sftp:
    #        print("Connection successfully established ... ")
    #        sftp.cwd('/home/ftpanuj')
    #        sftp.put(rootDir+file, '/home/ftpanuj/'+file)
	
    # step5.make_winds(now)
    # step6.par_run(now)
    # step7.postprocess_SWAN(now,1)

    # #step 8 making flood risk using multiprocessing
    # # creating processes
    # mp.freeze_support()
    # p1 = mp.Process(target=step8_Niutao.inundation_forecast, args = (now,))
    # p2 = mp.Process(target=step8_Niulakita.inundation_forecast, args = (now,))
    # p3 = mp.Process(target=step8_Nanumaga.inundation_forecast, args = (now,))
    # p4 = mp.Process(target=step8_Nanumea.inundation_forecast, args = (now,))
    # p5 = mp.Process(target=step8_Funafuti.inundation_forecast, args = (now,))
    # p6 = mp.Process(target=step8_Nukufetau.inundation_forecast, args = (now,))
    # p7 = mp.Process(target=step8_Nukulaelae.inundation_forecast, args = (now,))
    # p8 = mp.Process(target=step8_Nui.inundation_forecast, args = (now,))
    # p9 = mp.Process(target=step8_Vaitupu.inundation_forecast, args = (now,))
    # p10 = mp.Process(target=step8_Fun_lagoon, args = (now,lagoon_path,))

    # # starting processes
    # p1.start()
    # p2.start()
    # p3.start()
    # p4.start()
    # p5.start()
    # p6.start()
    # p7.start()
    # p8.start()
    # p9.start()
    # p10.start()

    # # wait until processes are finished
    # p1.join()
    # p2.join()
    # p3.join()
    # p4.join()
    # p5.join()
    # p6.join()
    # p7.join()
    # p8.join()
    # p9.join()
    # p10.join()

    # # delete previous archived swan outputs older than 3 months
    # try:
    #     delete_nmonthsbefore(now,3)
    # except Exception as e:
    #     print(e)
        
    # step9.archive_output(now)
    # step10.ingest2GUI(now)


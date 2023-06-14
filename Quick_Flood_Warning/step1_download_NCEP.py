import os, glob
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import subprocess
from datetime import timedelta
from bs4 import BeautifulSoup
import datetime as dt
from datetime import timedelta
import shutil
import step2_download_CMEMS as step2
import step3_gen_tide_TPOX8 as step3

def convert_grib_2_nc(grb_fl_name,nc_fl_name):
# function to convert grib2 files to netCDF, it uses the java script "toolsUI-5.4.1.jar" that needs to be placed in the same folder where this script is being called    
    subprocess.call(["java", "-Xmx512m", "-cp", "toolsUI-5.4.1.jar", "ucar.nc2.dataset.NetcdfDataset", "-in", grb_fl_name, "-out", nc_fl_name], shell=True,)
    print(grb_fl_name + ' converted to ' + nc_fl_name)
    return()


def download_all_wind_grb_and_convert_2_nc(mydate,Tcycle,dt,end_tt,leftlon,rightlon,toplat,bottomlat,grb_out,nc_out):

# function to download the winds
# mydate= current date
# Tcycle= forecast run, 00,06,12,18
# dt= time step, usually 3 hours
# end_dt= end of the forecast period
# leftlon, rightlon,toplat, bottomlat= limits of the area to be downloaded
# grb_out= root name of grib2 files
# nc_out= root name of netcdf file
    
    # define parameters for the connection with the server
    session = requests.Session()
    retry = Retry(connect=50, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
      
    # there are independent files for each time step, we need to go through all the times in the forecast period
    for i in range(0,end_tt,dt):
        
        # url to define the time and region to download
        url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?' + \
        'file=gfs.t' + Tcycle + 'z' + \
        '.pgrb2.0p25.f' + "{0:0>3}".format(i) + \
        '&lev_10_m_above_ground=on&var_UGRD=on&var_VGRD=on&subregion=&leftlon=' + leftlon + \
        '&rightlon=' + rightlon + \
        '&toplat=' + toplat + \
        '&bottomlat=' + bottomlat + \
        '&dir=%2Fgfs.' + mydate + '%2F' + Tcycle + '%2Fatmos'
        
        session.mount(url, adapter)  
        r=session.get(url)  
        

        grb_outX = grb_out + '_' + "{0:0>3}".format(i) +'.grib2'
        nc_outX = nc_out + '_' + "{0:0>3}".format(i) +'.nc'

        # write grib2 on disk
        open(grb_outX, 'wb').write(r.content)
        print('Grib file downloaded and stored as ' + grb_outX) 
	# convert to netcdf       
        convert_grib_2_nc(grb_outX,nc_outX)
    session.close()    
    return()

def download_all_wave_grb_and_convert_2_nc(mydate,Tcycle,dt,end_tt,leftlon,rightlon,toplat,bottomlat,grb_out,nc_out):

# function to download wave partitions
# mydate= current date
# Tcycle= forecast run, 00,06,12,18
# dt= time step, usually 3 hours
# end_dt= end of the forecast period
# leftlon, rightlon,toplat, bottomlat= limits of the area to be downloaded
# grb_out= root name of grib2 files
# nc_out= root name of netcdf file
    
    session = requests.Session()
    retry = Retry(connect=10, backoff_factor=0.5)
    adapter = HTTPAdapter(max_retries=retry)
    # session.mount('http://', adapter)
    
    for i in range(0,end_tt,dt):
        # url to define the time and region to download
        url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gefs_wave_0p25.pl?' + \
        'file=gefs.wave.t' + Tcycle + 'z.c00' + \
        '.global.0p25.f' + "{0:0>3}".format(i) + \
        '.grib2&all_lev=on&all_var=on&subregion=&leftlon=' + leftlon + \
        '&rightlon=' + rightlon + \
        '&toplat=' + toplat + \
        '&bottomlat=' + bottomlat + \
        '&dir=%2Fgefs.' + mydate + '%2F' + Tcycle + '%2Fwave%2Fgridded'
        print(url)
        session.mount(url, adapter)  
        r=session.get(url)

        grb_outX = grb_out + '_' + "{0:0>3}".format(i) +'.grib2'
        nc_outX = nc_out + '_' + "{0:0>3}".format(i) +'.nc'

        # write grib2 on disk
        open(grb_outX, 'wb').write(r.content)
        print('Grib file downloaded and stored as ' + grb_outX)
        # convert to netcdf 
        convert_grib_2_nc(grb_outX,nc_outX)
    session.close()    
    return()
###############################################################################

def download_NCEP(now):

    # now= date and hour of the NOAA/NCEP run to download, it is a datetime object but this can be changed to our needs. I define this date in a previous script that checks the latest available run on the server

    Tcycle = str(now.hour).zfill(2)
    # time interval of the forecast, usually 3 hours
    wave_dt = 3
    wind_dt = 3

    # from the current time downlad the next 180 (+1) hours (7.5 days)
    time_length = 181
    # I also dowload the first 21 hours of the previous day to spin up the model (1 day), if you use hotstart you can negret this step
    htime_length = 22
    # coordinates of the region to be downloades
    leftlon = '191.0'
    rightlon = '206.0'
    toplat = '-5.0'
    bottomlat = '-26.0'
    
    # define folder and root names to download the date, I put them all in a temporal folder 
    wave_grb_out = 'tmp/wave_tmp'
    wave_nc_out = 'tmp/wave_tmp'
    wnd_grb_out = 'tmp/wind_tmp'
    wnd_nc_out = 'tmp/wind_tmp'
    
    # wave_grb_outh = 'tmp/h_wave_tmp'
    # wave_nc_outh = 'tmp/h_wave_tmp'
    # wnd_grb_outh = 'tmp/h_wind_tmp'
    # wnd_nc_outh = 'tmp/h_wind_tmp'
    
    # wave_grb_outhh = 'tmp/hh_wave_tmp'
    # wave_nc_outhh = 'tmp/hh_wave_tmp'
    # wnd_grb_outhh = 'tmp/hh_wind_tmp'
    # wnd_nc_outhh = 'tmp/hh_wind_tmp'
    
    ###############################################################################
    # mydate is a string 20210902
    mydate = now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
    # yesterday is the day before to do the same for myhdate
  #  yesterday = now - timedelta(1)
   # myhdate = yesterday.strftime("%Y") + yesterday.strftime("%m") + yesterday.strftime("%d")
    
    #byesterday = yesterday - timedelta(1)
   # myhhdate = byesterday.strftime("%Y") + byesterday.strftime("%m") + byesterday.strftime("%d") 

    # download wind and waves for the forecast period (7.5 days)
    #download_all_wind_grb_and_convert_2_nc(mydate,Tcycle,wind_dt,time_length,leftlon,rightlon,toplat,bottomlat,wnd_grb_out,wnd_nc_out)
    download_all_wave_grb_and_convert_2_nc(mydate,Tcycle,wave_dt,time_length,leftlon,rightlon,toplat,bottomlat,wave_grb_out,wave_nc_out)

    # # download wind and waves for yesterday)
    # download_all_wind_grb_and_convert_2_nc(myhdate,Tcycle,wind_dt,htime_length,leftlon,rightlon,toplat,bottomlat,wnd_grb_outh,wnd_nc_outh)
    # download_all_wave_grb_and_convert_2_nc(myhdate,Tcycle,wave_dt,htime_length,leftlon,rightlon,toplat,bottomlat,wave_grb_outh,wave_nc_outh)    
    
    # # download wind and waves for byesterday)
    # download_all_wind_grb_and_convert_2_nc(myhhdate,Tcycle,wind_dt,htime_length,leftlon,rightlon,toplat,bottomlat,wnd_grb_outhh,wnd_nc_outhh)
    # download_all_wave_grb_and_convert_2_nc(myhhdate,Tcycle,wave_dt,htime_length,leftlon,rightlon,toplat,bottomlat,wave_grb_outhh,wave_nc_outhh)    
    
    
    # remove all the grib2files
    for filename in glob.glob("../tmp/*.grib2"):
        os.remove(filename) 
    

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



# now = dt.datetime.utcnow()
# url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gefs_wave_0p25.pl?dir=%2Fgefs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
# runs = list_available_runs(url)
# if len(runs)==0:
#     now = dt.datetime.utcnow()-timedelta(1)
#     url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gefs_wave_0p25.pl?dir=%2Fgefs.' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d")
#     runs = list_available_runs(url)
   
# #- Define the run to be used
# runs=sorted(runs)
# now = now.replace(hour=runs[-1],minute=0,second=0,microsecond=0)
# #now = dt.datetime(2022,7,12,0)
   
   
# # delete previous runs older than 14 days
# try:
#     delete_ndaysbefore(now,14)
# except Exception as e:
#     print(e)



# download_NCEP(now)
# step2.download_CNEMS(now)
# step3.gen_tide(now)










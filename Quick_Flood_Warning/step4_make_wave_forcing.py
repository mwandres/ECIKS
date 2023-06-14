# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 09:40:26 2021

@author: 
"""
import os
import netCDF4
import numpy as np
from netCDF4 import num2date
import pandas as pd
#import create_spec_from_partitions
import netCDF4 as nc
from scipy.interpolate import griddata as griddata
import datetime as dt
from scipy.interpolate import interp1d
import matplotlib.dates as mdates
from matplotlib import pyplot as plt
# Functions
def read_netcdf(nc_fname,dt,length_dt):
# Function to read wave partitions contained in the netCDF files from step1_download_NCEP.py, gives back wave partitions for the whole tile. 
# This function also corrects wave height partitions to make their summation to match with the total Hs produced by WWIII before partitioning the spectra

    # read wave partition for hour=0
    nc_fnameX = nc_fname + '_' + "{0:0>3}".format(0) +'.nc'
    nc = netCDF4.Dataset(nc_fnameX)
    # kk = np.array(nc['ordered_sequence_of_data'])
    ln = np.array(nc['lon'])
    lt = np.array(nc['lat'])
    tt = np.array(nc['time'])
    nc_unit = nc['time'].units
    nc_calendar = nc['time'].calendar
    Hs_ST = np.nan_to_num(np.array(nc['Significant_height_of_combined_wind_waves_and_swell_surface'])[:,:,:])
    Tm_ST = np.nan_to_num(np.array(nc['Primary_wave_mean_period_surface'])[:,:,:])
    Dir_ST = np.nan_to_num(np.array(nc['Primary_wave_direction_surface'])[:,:,:])
    Hs_S1 = np.nan_to_num(np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:])
    Hs_S2 = np.nan_to_num(np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:])
    Hs_S3 = np.nan_to_num(np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:])
    Hs_W = np.nan_to_num(np.array(nc['Significant_height_of_wind_waves_surface'])[:,:,:])   
    Tm_S1 = np.nan_to_num(np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:])
    Tm_S2 = np.nan_to_num(np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:])
    Tm_S3 = np.nan_to_num(np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:])
    Tm_W = np.nan_to_num(np.array(nc['Mean_period_of_wind_waves_surface'])[:,:,:])
    Dir_S1 = np.nan_to_num(np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:])
    Dir_S2 = np.nan_to_num(np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:])
    Dir_S3 = np.nan_to_num(np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:])
    Dir_W = np.nan_to_num(np.array(nc['Direction_of_wind_waves_surface'])[:,:,:])    
   # read wave partition from hour=3 to hour=180 every 3 hours
    for i in range(dt,length_dt,dt):
        nc_fnameX = nc_fname + '_' +  "{0:0>3}".format(i) +'.nc'
        nc = netCDF4.Dataset(nc_fnameX)
        tt = np.append(tt,np.array(nc['time']),axis=0)
        Hs_ST = np.nan_to_num(np.append(Hs_ST,np.array(nc['Significant_height_of_combined_wind_waves_and_swell_surface'])[:,:,:],axis=0))
        Tm_ST = np.nan_to_num(np.append(Tm_ST,np.array(nc['Primary_wave_mean_period_surface'])[:,:,:],axis=0))
        Dir_ST = np.nan_to_num(np.append(Dir_ST,np.array(nc['Primary_wave_direction_surface'])[:,:,:],axis=0))
        
        Hs_S1 = np.nan_to_num(np.append(Hs_S1,np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:],axis=0))
        Hs_S2 = np.nan_to_num(np.append(Hs_S2,np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:],axis=0))
        Hs_S3 = np.nan_to_num(np.append(Hs_S3,np.array(nc['Significant_height_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:],axis=0))
        Hs_W = np.nan_to_num(np.append(Hs_W,np.array(nc['Significant_height_of_wind_waves_surface'])[:,:,:],axis=0))
        Tm_S1 = np.nan_to_num(np.append(Tm_S1,np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:],axis=0))
        Tm_S2 = np.nan_to_num(np.append(Tm_S2,np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:],axis=0))
        Tm_S3 = np.nan_to_num(np.append(Tm_S3,np.array(nc['Mean_period_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:],axis=0))
        Tm_W = np.nan_to_num(np.append(Tm_W,np.array(nc['Mean_period_of_wind_waves_surface'])[:,:,:],axis=0))
        Dir_S1 = np.nan_to_num(np.append(Dir_S1,np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,0,:,:],axis=0))
        Dir_S2 = np.nan_to_num(np.append(Dir_S2,np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,1,:,:],axis=0))
        Dir_S3 = np.nan_to_num(np.append(Dir_S3,np.array(nc['Direction_of_swell_waves_ordered_sequence_of_data'])[:,2,:,:],axis=0))
        Dir_W = np.nan_to_num(np.append(Dir_W,np.array(nc['Direction_of_wind_waves_surface'])[:,:,:],axis=0))
        print(nc_fnameX)
        
    # determine calbration coeficients
    factor= Hs_ST/np.sqrt(Hs_S1**2+Hs_S2**2+Hs_S3**2+Hs_W**2)
    # apply calbration coeficients
    Hs_S1=Hs_S1*factor
    Hs_S2=Hs_S2*factor
    Hs_S3=Hs_S3*factor
    Hs_W=Hs_W*factor
    
    
    dates = num2date(tt[:], units=nc_unit, calendar = nc_calendar)
    return(ln,lt,dates,Hs_ST,Tm_ST,Dir_ST,Hs_S1,Hs_S2,Hs_S3,Hs_W,Tm_S1,Tm_S2,Tm_S3,Tm_W,Dir_S1,Dir_S2,Dir_S3,Dir_W)

def SWAN2nc(ln,lt,tt,Hs_ST,Tm_ST,Dir_ST,xll,yll,xur,yur,inc,namenc):
    
    x = np.arange(xll,xur,inc) 
    y = np.arange(yll,yur,inc)
    lgrid_x, lgrid_y = np.meshgrid(x,y)
    ln_x, lt_y = np.meshgrid(ln,lt)

    Hsp = np.ones((len(tt),lgrid_x.shape[0],lgrid_x.shape[1]))*-999
    Tmp = np.ones((len(tt),lgrid_x.shape[0],lgrid_x.shape[1]))*-999
    Dirp = np.ones((len(tt),lgrid_x.shape[0],lgrid_x.shape[1]))*-999 

    for ts in range(len(tt)): # change the starting hour into range to the lenght of the spin-up 
        
        Hsp[ts,:,:] = griddata((ln_x.flatten(),lt_y.flatten()),Hs_ST[ts,:,:].flatten(),(lgrid_x,lgrid_y))
        Tmp[ts,:,:] = griddata((ln_x.flatten(),lt_y.flatten()),Tm_ST[ts,:,:].flatten(),(lgrid_x,lgrid_y))
        Dirp[ts,:,:] = griddata((ln_x.flatten(),lt_y.flatten()),Dir_ST[ts,:,:].flatten(),(lgrid_x,lgrid_y))
        
    ds = nc.Dataset(namenc, 'w', format='NETCDF4')
    time = ds.createDimension('time', None)
    lon = ds.createDimension('lon', None)
    lat = ds.createDimension('lat', None)
    times = ds.createVariable('time', 'f4', ('time',))
    times.units='hours since 1950-01-01 00:00:00'
    times.calendar='gregorian'
    times[:] = [nc.date2num(ti,units=times.units,calendar=times.calendar) for ti in tt]
    lonnc =ds.createVariable('lon', 'f4', ('lon',))
    lonnc.units ='degrees_east'
    lonnc[:] = x
    latnc = ds.createVariable('lat', 'f4', ('lat',))
    latnc.units ='degrees_north'
    latnc[:] = y
    Hsnc= ds.createVariable('Hs', 'f4', ('time','lat','lon'))
    Hsnc.units = 'm'
    Hsnc[:] =  Hsp
    Tmnc= ds.createVariable('Tm', 'f4', ('time','lat','lon'))
    Tmnc.units = 's'
    Tmnc[:] =  Tmp
    Dirnc= ds.createVariable('Dir', 'f4', ('time','lat','lon'))
    Dirnc.units = 'degress from north (north=0, east=90)'
    Dirnc[:] =  Dirp
    ds.close()
    return(Hsp,Tmp,Dirp,tt,x,y)
    
def read_loc_csv(fpath):
    df = pd.read_csv(fpath)
    point_no = df['point_no'].tolist()
    atoll = df['atoll'].tolist()
    lon = df['lon'].tolist()
    lat = df['lat'].tolist()
    sla_offset = df['sla_offset'].tolist()
    shore_normal_angle = df['shore_normal_angle'].tolist()
    threshold_1 = df['threshold_1'].tolist()
    threshold_2 = df['threshold_2'].tolist()
    return(point_no,atoll,lon,lat,sla_offset,shore_normal_angle,threshold_1,threshold_2)

def find_nearest(array,target_value):
    idx = (np.abs(array-target_value)).argmin()
    return(array[idx], idx)

def get_closest_wave_data(Hsp,Tmp,Dirp,wave_dates,x,y,target_lon,target_lat):
    # Function to extract time series of wave partitions at the required location,
    # determine directional spreading accoring to a relation obtained in the area
    # returns the 2D spectra for all times at one location
    # x,y are the indices of the point
    # freqs frequency vector for the 2d spectra (in Hz)
    # theta, directional vector in Deg
    # Hs_S... wave partition parameters
    lonX, lon_ix = find_nearest(x,target_lon)
    latX, lat_ix = find_nearest(y,target_lat)
    pHsp = Hsp[:,lat_ix,lon_ix]
    pTmp = Tmp[:,lat_ix,lon_ix]
    pDirp = Dirp[:,lat_ix,lon_ix]                
    return(pHsp,pTmp,pDirp)

def readvars_nc(nc_fname,varname):
    ncc = nc.Dataset(nc_fname,'r') 
    if varname == 'time':
        time_var = ncc.variables[varname] 
        time_or = nc.num2date(ncc[varname],time_var.units, time_var.calendar)
        time_str=[time_or[i].strftime('%Y-%m-%d %H:%M') for i in range(0,len(time_or))]
        var = [dt.datetime.strptime(time_str[i],'%Y-%m-%d %H:%M') for i in range(0,len(time_str))]
    else:
        var = np.array(ncc[varname])
    ncc.close()
    return(var)

def calculate_runup_from_Merrifield(theta_N,Hs,Tp,Dir):
    gamma = 1.0
    b1 = 0.33
    b0 = -0.1
    
    angles = np.cos(np.deg2rad(Dir - theta_N))
    if angles < 0:
        angles = 0
        
    Hb = ((Hs**2) * Tp * (4*np.pi)**(-1) * angles * np.sqrt(gamma*9.81))**(2/5)
    RU2 = b1 * Hb + b0
    
    return(RU2)




def make_ts_fig(time_waves,twl_nearshore,tide,time_tide_minute,tide_minute,mx_lvl1,mx_lvl2,outfile_name):
    time_waves_pd = pd.to_datetime(time_waves)
    
    xor = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_waves]
    xnew = [nc.date2num(i,'minutes since 1950-01-01 00:00:00') for i in time_tide_minute]
    f = interp1d(xor,twl_nearshore,kind='cubic',fill_value="extrapolate")
    twl_min=f(xnew)
    
    limits = [-1,np.max(mx_lvl2)+0.25]
    dayFmt = mdates.DateFormatter('%d-%m-%Y %H:%M')
    
    #-- differentiate to calculate high and low tides
    diff = np.zeros_like(time_tide_minute, dtype=np.float64)
    #-- forward differentiation for starting point
    diff[0] = tide_minute[1] - tide_minute[0]
    #-- backward differentiation for end point
    diff[-1] = tide_minute[-1] - tide_minute[-2]
    #-- centered differentiation for all others
    diff[1:-1] = (tide_minute[2:] - tide_minute[0:-2])/2.0
    htindex, = np.nonzero((np.sign(diff[0:-1]) >= 0) & (np.sign(diff[1:]) < 0))
    # ltindex, = np.nonzero((np.sign(diff[0:-1]) <= 0) & (np.sign(diff[1:]) > 0))

    #plt.figure(figsize=(10,6))
    f, ax = plt.subplots(figsize=(12,6))
    
    #ax.axhspan(mx_lvl1,mx_lvl2,facecolor='yellow',alpha=0.5)
    #ax.axhspan(mx_lvl2,8,color='r',alpha=0.5)
    
    ax.plot(time_waves_pd[0:],np.matlib.repmat(mx_lvl2,len(time_waves_pd[0:]),1),color='firebrick')
    ax.text(time_waves_pd[-5600], mx_lvl2+0.05, 'Moderate flood threshold',color='firebrick',size=11)
    ax.plot(time_waves_pd[0:],np.matlib.repmat(mx_lvl1,len(time_waves_pd[0:]),1),color='darkorange')
    plt.text(time_waves_pd[-5600], mx_lvl1+0.05, 'Minor flood threshold',color='darkorange',size=11)
    
    ax.fill_between(time_waves_pd[0:],twl_nearshore[0:],-2,color='lightsteelblue',label='Nearshore TWL')
    ax.fill_between(time_waves_pd[0:],tide[0:],-2,color='royalblue',label='Offshore Tide')
    for h in range(0,len(htindex)-4):
            text=time_tide_minute[htindex[h]].strftime("%H:%M")
            plt.plot(time_tide_minute[htindex[h]],twl_min[htindex[h]],"v",color= 'gray',markersize=5)
            plt.text(time_tide_minute[htindex[h]]-dt.timedelta(hours=3),twl_min[htindex[h]]+0.08,text,color= 'gray',size=15) 
    ax.set_xlim(time_waves_pd[0],time_waves_pd[-2940])
    ax.set_ylim(limits)
    ax.set_ylabel('Water level (m)')
    ax.grid(True,alpha = 0.7)
    ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.xaxis.set_minor_locator(mdates.HourLocator(byhour = [6,12,18]))
    ax.xaxis.set_major_formatter(dayFmt)
    f.autofmt_xdate()
    ax.legend()
    f.savefig(outfile_name)
    plt.close(f)
    return()


def make_risk_csv_file(output_twl_nearshore,mx_lvl1,mx_lvl2,lon,lat,risk_outfile_name):
    if output_twl_nearshore.max() >= mx_lvl2:
        risk_cat = 2
    elif output_twl_nearshore.max() >= mx_lvl1:
        risk_cat = 1
    else:
        risk_cat = 0
    
    risk_output = np.array((lon,lat,risk_cat))
    risk_output = np.reshape(risk_output,(len(risk_output),1)).T
    np.savetxt(risk_outfile_name, risk_output, delimiter=",")
    return()   

###############################################################################
def make_waves(now):
    print('Generating spectral wave forcing')  
    # define the output folder were SWAN will be run
    out_name = 'runs/' + now.strftime("%Y") + now.strftime("%m") + now.strftime("%d") +  now.strftime("%H") 

    ###############################################################################
    #freqs = np.arange(0.03,1.,.01)
    # logarithmic scale definition inf frequencise as it is in SWAN manual,  increases resolution on the long periods, larger bins on the short waves 
    # freqs = 0.0373*np.power(1.1,np.arange(0,32,1))
    # theta = np.arange(0,360,10)
    ###############################################################################
   
    # directory of the nc files
    wave_nc = 'tmp/wave_tmp'
    # h_wave_nc = '../tmp/h_wave_tmp'
    # hh_wave_nc = '../tmp/hh_wave_tmp'
     
    wave_dt = 3
    time_length = 181
    # h_time_length = 22
     
    # read wave nc files
    ln,lt,wave_dates,Hs_ST,Tm_ST,Dir_ST,Hs_S1,Hs_S2,Hs_S3,Hs_W,Tm_S1,Tm_S2,Tm_S3,Tm_W,Dir_S1,Dir_S2,Dir_S3,Dir_W = read_netcdf(wave_nc,wave_dt,time_length) 
    #ln,lt,h_tt,h_Hs_S1,h_Hs_S2,h_Hs_S3,h_Hs_W,h_Tm_S1,h_Tm_S2,h_Tm_S3,h_Tm_W,h_Dir_S1,h_Dir_S2,h_Dir_S3,h_Dir_W = read_netcdf(h_wave_nc,wave_dt,h_time_length) 
    #ln,lt,hh_tt,hh_Hs_S1,hh_Hs_S2,hh_Hs_S3,hh_Hs_W,hh_Tm_S1,hh_Tm_S2,hh_Tm_S3,hh_Tm_W,hh_Dir_S1,hh_Dir_S2,hh_Dir_S3,hh_Dir_W = read_netcdf(hh_wave_nc,wave_dt,h_time_length) 
     
    xll = 191.0
    yll = -26.0
    xur = 206.0
    yur = -5.0
    inc =0.01# about1km
    namenc = 'results/CookIslands.nc'
    output_table = 'extras/output_points.csv'
    Hsp,Tmp,Dirp,wave_dates,x,y = SWAN2nc(ln,lt,wave_dates,Hs_ST,Tm_ST,Dir_ST,xll,yll,xur,yur,inc,namenc)
    
    point_no,atoll,point_lon,point_lat,sla_offset,shore_normal_angle,threshold_1,threshold_2 = read_loc_csv(output_table)
    for i in range(len(atoll)):
        print("\n",atoll[i])
        target_lon = point_lon[i]
        target_lat = point_lat[i]
        pHsp,pTmp,pDirp = get_closest_wave_data(Hsp,Tmp,Dirp,wave_dates,x,y,target_lon,target_lat)
        
        sla_name = 'tmp/sla_hourly_' + atoll[i] + '.nc'
        tide_name = 'tmp/tide_minute_' + atoll[i] + '.nc'
        
        time_sla = readvars_nc(sla_name,'time')
        sla = readvars_nc(sla_name,'SLA')/100
        time_tide = readvars_nc(tide_name,'time')
        tide = readvars_nc(tide_name,'tide')/100
        
        theta_N = shore_normal_angle[i]
        
        RU2 = np.ones((len(pHsp)))*-999
        for jj in range(len(RU2)):
            RU2[jj] = calculate_runup_from_Merrifield(theta_N,pHsp[jj],pTmp[jj],pDirp[jj])
        
        xor = [nc.date2num(jj,'minutes since 1950-01-01 00:00:00') for jj in wave_dates]
        xnew = [nc.date2num(jj,'minutes since 1950-01-01 00:00:00') for jj in time_tide]
        f = interp1d(xor,RU2,kind='linear',fill_value="extrapolate")
        RU2_min=f(xnew)
        
        xor = [nc.date2num(jj,'minutes since 1950-01-01 00:00:00') for jj in time_sla]
        xnew = [nc.date2num(jj,'minutes since 1950-01-01 00:00:00') for jj in time_tide]
        f = interp1d(xor,sla,kind='linear',fill_value="extrapolate")
        sla_min = f(xnew)
        
        twl_nearshore = RU2_min + sla_min + tide
        
        mx_lvl1 = threshold_1[i]
        mx_lvl2 = threshold_2[i]
        outfile_name = 'results/' + atoll[i] + '_' + str(point_no[i]) + '.png'
        risk_outfile_name = 'results/' + atoll[i] + '_' + str(point_no[i]) + '_risk.csv'
        make_ts_fig(time_tide,twl_nearshore,tide,time_tide,tide,mx_lvl1,mx_lvl2,outfile_name)
        
        make_risk_csv_file(twl_nearshore[0:-2940],mx_lvl1,mx_lvl2,point_lon[i],point_lat[i],risk_outfile_name)
    return()
    
    
    
    
    
    
    
    
    # #%%
      
    # try:
    #     os.mkdir(out_name)
    # except OSError as error:
    #     print(error)    
  
    #  # read the location of the boundary points of the SWAN unestructured mesh
    # f14bnds = pd.read_csv("../extras/f14_boundary_points.csv",sep=",")
    # f14_xb = np.array(f14bnds["x_f14_bnd"])
    # f14_yb = np.array(f14bnds["y_f14_bnd"])
    # f14_id = np.array(f14bnds["adcirc_index"])
  
    
    # # find the indices of the closest point in the wave grid to each boundary point of SWAN mesh
    # yv, xv = np.meshgrid(lt,ln, sparse=False, indexing='ij')   
    # lon=np.reshape(xv, (np.size(xv), 1)) 
    # lat=np.reshape(yv, (np.size(yv), 1))
     
    # Wf_index=[]
    # for k in range(len(f14_xb)):
    #     dist=np.squeeze(np.sqrt((lon-f14_xb[k])**2+(lat-f14_yb[k])**2))
    #     kk=np.argsort(dist)[0]
    #     Wf_index.append(kk) 
       
    #  # this part smooths the wave spectra along the boundary, needed when you work with unstructured SWAN meshes where boundary points are  defined independently  
    # shape=np.shape(xv)
    # Wf_index=np.asarray( Wf_index)
     
    # for l in range(len(f14_xb)):
        
    #     pos=np.arange(l-1,l+2)
    #     pos[pos<0]=pos[pos<0]+len(f14_xb)
    #     pos[pos>len(f14_xb)-1]=pos[pos>len(f14_xb)-1]-len(f14_xb)
        
    #     out_nameXX=out_name+'/Pto_' + str(f14_id[l]) + '.sp2'   
    #     [i, j] = np.unravel_index(Wf_index[pos], shape) 
         
    #     indices=np.vstack([i,j]).T
    #     uind=np.unique(indices,axis=0)
        
    #     MATT_hh=np.zeros((np.size(hh_tt),np.size(freqs),np.size(theta)))
    #     MATT_h=np.zeros((np.size(h_tt),np.size(freqs),np.size(theta)))
    #     MATT_f=np.zeros((np.size(tt),np.size(freqs),np.size(theta)))
        
  
        
    #     for m in range(len(uind)):
            
    #         iii=uind[m,0]
    #         jjj=uind[m,1]
           
    #         MAT_hh = extract_entire_time_series(hh_tt,iii,jjj,freqs,theta,hh_Hs_S1,hh_Hs_S2,hh_Hs_S3,hh_Hs_W,hh_Tm_S1,hh_Tm_S2,hh_Tm_S3,hh_Tm_W,hh_Dir_S1,hh_Dir_S2,hh_Dir_S3,hh_Dir_W)
    #         MATT_hh=MATT_hh+MAT_hh
    #         MAT_h = extract_entire_time_series(h_tt,iii,jjj,freqs,theta,h_Hs_S1,h_Hs_S2,h_Hs_S3,h_Hs_W,h_Tm_S1,h_Tm_S2,h_Tm_S3,h_Tm_W,h_Dir_S1,h_Dir_S2,h_Dir_S3,h_Dir_W)
    #         MATT_h=MATT_h+MAT_h
    #         MAT_f = extract_entire_time_series(tt,iii,jjj,freqs,theta,Hs_S1,Hs_S2,Hs_S3,Hs_W,Tm_S1,Tm_S2,Tm_S3,Tm_W,Dir_S1,Dir_S2,Dir_S3,Dir_W)
    #         MATT_f=MATT_f+MAT_f
            
    #     TT=np.concatenate((hh_tt,h_tt,tt), axis=0) 
    #     MATT=np.concatenate((MATT_hh,MATT_h,MATT_f), axis=0)/len(uind)
        
    #     # create SWAN boundary file for each point in the boundary, this can be used to generate a file with all the points for regular SWAN grids
    #     #makeSWANspec(out_nameXX,TT,freqs,theta,MATT,f14_xb[l],f14_yb[l],1)
    #     # time_str=[TT[i].strftime('%Y-%m-%d %H:%M') for i in range(0,len(tt))]
    # print('Wave forcing generated in: ' + out_name)
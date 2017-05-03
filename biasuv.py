#!/usr/bin/env python
from mpi4py import MPI
#me = MPI.COMM_WORLD.Get_rank()
#nproc = MPI.COMM_WORLD.Get_size()
#print me, nproc

from mpl_toolkits.basemap import Basemap, cm, shiftgrid,interp
from netCDF4 import Dataset as NetCDFFile
import numpy as N
import matplotlib
import glob
import numpy.ma as ma
from scipy.interpolate import griddata
from datetime import datetime, timedelta
from netCDF4 import num2date, date2num



tk=124+112+124+120+124+120+124+124+120+124+120+124
#print tk
clmday=(124+112+124+120+124+120+124+124+120+124+120+124)*2
#print clmday
tstep=range(1,1461)
#print tstep
uall= N.zeros((11,tk, 360, 720))
vall= N.zeros((11,tk, 360, 720))

clm85uall= N.zeros((11,clmday, 360, 720))
clm85vall= N.zeros((11,clmday, 360, 720))

clm85allu= N.zeros((11,tk, 360, 720))
clm85allv= N.zeros((11,tk, 360, 720))


crudevu=N.zeros((tk, 360, 720))
clmdevu=N.zeros((tk, 360, 720))
crudevv=N.zeros((tk, 360, 720))
clmdevv=N.zeros((tk, 360, 720))






years = range(2006, 2017)
months= range(1,13)
q=0
c=0
kx=0

for i, year in enumerate(years):
        q=0
        q1=0
        for m, month in enumerate(months):
                if month < 10 :
                        month = '0{}'.format(month)
                print "/scratch2/scratchdirs/tslin2/bias/data/cru/{0}-{1}.nc".format(year,month)
                base = NetCDFFile ("/scratch2/scratchdirs/tslin2/bias/data/cru/{0}-{1}.nc".format(year,month), mode='r')
                lon = base.variables["longitude"][:]
                lat = base.variables["latitude"][:]
                time  = base.variables["Time"][:]
#                print time,time.shape 
                for x, tt in enumerate(time):

                   c=x+q
                   tempa = base.variables["UWIND"][x,:,:]
                   tempb = base.variables["VWIND"][x,:,:]

                   uall[i,c, :, :] = N.sqrt((tempa*tempa)+(tempb*tempb))


                q=q+tt
                #print q

                print "/scratch2/scratchdirs/tslin2/bias/data/clm45/{0}_{1}.nc".format(year,month)
                base1 = NetCDFFile ("/scratch2/scratchdirs/tslin2/bias/data/clm45/{0}_{1}.nc".format(year,month), mode='r')
                trclm=N.max(time)*2
                ab=0
                for xy in range(0,300):
                    if  ab < trclm:
                         #print tclm,ab
                         kx=ab+q1
                         temp2= base1.variables["UWIND"][ab,:,:]
                         temp3= base1.variables["VWIND"][ab,:,:]
                         xyz=(temp2*temp2)+(temp3*temp3)
                         clm85uall[i,kx, :, :] = N.sqrt(xyz)


                         ab=ab+1
                q1=q1+trclm
                #print q1
                base.close()
                base1.close()

#retriece 6 hourly
        for dt in range(0,tk):
            dd=dt*2
            clm85allu[i,dt,:,:]=clm85uall[i,dd, :, :]



uempymean=N.average(uall,axis=0)
umean=N.average(uempymean,axis=0)





uclmy=N.average(clm85allu,axis=0)
uclm=N.average(uclmy,axis=0)

crudevu=N.std(uall,axis=0)
clmdevu=N.std(clm85allu,axis=0)



ncfile=NetCDFFile('meanws45.nc','w',format='NETCDF4')
ncfile.description='clm and cruncep wind speed no direction 2006-2016 mean 6hourly'

# dimensions

timek=ncfile.createDimension('time',tk)
latitude=ncfile.createDimension('lat', 360)
longitude=ncfile.createDimension('lon', 720)


# variables


latitudes = ncfile.createVariable('lat', 'f8', ('lat',))
longitudes = ncfile.createVariable('lon', 'f8', ('lon',))
times = ncfile.createVariable('time', 'f8', ('time',))




au= ncfile.createVariable('cruu', 'f8', ('time','lat','lon'),fill_value=-9999.)
bu= ncfile.createVariable('clm85u', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevau= ncfile.createVariable('crudevu', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevbu= ncfile.createVariable('clmdevu', 'f8', ('time','lat','lon'),fill_value=-9999.)

#print fractions

# data


latitudes[:] = lat
longitudes[:] = lon
#dates = [datetime(2006,1,1)+n*timedelta(hours=06) for n in range(tstep.shape[0])]
#times[:]=date2num(dates,units=times.units,calendar=times.calendar)
times[:]=tstep



au[:,:,:]=uempymean
bu[:,:,:]=uclmy
odevau[:,:,:]=crudevu
odevbu[:,:,:]=clmdevu



latitudes.units = 'degrees_north'
longitudes.units = 'degrees_east'


ncfile.close()

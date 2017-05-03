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
tall= N.zeros((11,tk, 360, 720))
psall= N.zeros((11,tk, 360, 720))

clm85psall= N.zeros((11,clmday, 360, 720))

clm85allps= N.zeros((11,tk, 360, 720))

crudevps=N.zeros((tk, 360, 720))
clmdevps=N.zeros((tk, 360, 720))


clm85tall= N.zeros((11,clmday, 360, 720))

clm85allt= N.zeros((11,tk, 360, 720))

crudevt=N.zeros((tk, 360, 720))
clmdevt=N.zeros((tk, 360, 720))


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
                   tempe = base.variables["PSRF"][x,:,:]
                   psall[i,c, :, :] = tempe


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
                         temp7= base1.variables["PRSF"][ab,:,:]
                         clm85psall[i,kx, :, :] = temp7

                         ab=ab+1
                q1=q1+trclm
                #print q1
                base.close()
                base1.close()

#retriece 6 hourly
        for dt in range(0,tk):
            dd=dt*2
            clm85allps[i,dt,:,:]=clm85psall[i,dd, :, :]




psempymean=N.average(psall,axis=0)
psmean=N.average(psempymean,axis=0)



psclmy=N.average(clm85allps,axis=0)
psclm=N.average(psclmy,axis=0)

crudevps=N.std(psall,axis=0)
clmdevps=N.std(clm85allps,axis=0)


ncfile=NetCDFFile('meanps45.nc','w',format='NETCDF4')
ncfile.description='clm and crubecp t and pressure 2006-2016 mean 6hourly'

# dimensions

timek=ncfile.createDimension('time',tk)
latitude=ncfile.createDimension('lat', 360)
longitude=ncfile.createDimension('lon', 720)


# variables


latitudes = ncfile.createVariable('lat', 'f8', ('lat',))
longitudes = ncfile.createVariable('lon', 'f8', ('lon',))
times = ncfile.createVariable('time', 'f8', ('time',))


aps= ncfile.createVariable('crups', 'f8', ('time','lat','lon'),fill_value=-9999.)
bps= ncfile.createVariable('clm85ps', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevaps= ncfile.createVariable('crudevps', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevbps= ncfile.createVariable('clmdevps', 'f8', ('time','lat','lon'),fill_value=-9999.)

latitudes[:] = lat
longitudes[:] = lon
#dates = [datetime(2006,1,1)+n*timedelta(hours=06) for n in range(tstep.shape[0])]
#times[:]=date2num(dates,units=times.units,calendar=times.calendar)
times[:]=tstep

aps[:,:,:]=psempymean
bps[:,:,:]=psclmy
odevaps[:,:,:]=crudevps
odevbps[:,:,:]=clmdevps


latitudes.units = 'degrees_north'
longitudes.units = 'degrees_east'


ncfile.close()

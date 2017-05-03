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
                   tempes = base.variables["TBOT"][x,:,:]
                   tall[i,c, :, :] = tempes


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
                         temp1= base1.variables["TBOT"][ab,:,:]
                         clm85tall[i,kx, :, :] = temp1

                         ab=ab+1
                q1=q1+trclm
                #print q1
                base.close()
                base1.close()

#retriece 6 hourly
        for dt in range(0,tk):
            dd=dt*2
            clm85allt[i,dt,:,:]=clm85tall[i,dd, :, :]

tempymean=N.average(tall,axis=0)
#print tempymean.shape
tmean=N.average(tempymean,axis=0)
#print tmean.shape


tclmy=N.average(clm85allt,axis=0)
tclm=N.average(tclmy,axis=0)



crudevt=N.std(tall,axis=0)
clmdevt=N.std(clm85allt,axis=0)


ncfile=NetCDFFile('meant45.nc','w',format='NETCDF4')
ncfile.description='clm and crubecp t  2006-2016 mean 6hourly'

# dimensions

timek=ncfile.createDimension('time',tk)
latitude=ncfile.createDimension('lat', 360)
longitude=ncfile.createDimension('lon', 720)


# variables


latitudes = ncfile.createVariable('lat', 'f8', ('lat',))
longitudes = ncfile.createVariable('lon', 'f8', ('lon',))
times = ncfile.createVariable('time', 'f8', ('time',))

aps1= ncfile.createVariable('cru', 'f8', ('time','lat','lon'),fill_value=-9999.)
bps1= ncfile.createVariable('clm85', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevaps1= ncfile.createVariable('crudev', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevbps1= ncfile.createVariable('clmdev', 'f8', ('time','lat','lon'),fill_value=-9999.)



latitudes[:] = lat
longitudes[:] = lon
#dates = [datetime(2006,1,1)+n*timedelta(hours=06) for n in range(tstep.shape[0])]
#times[:]=date2num(dates,units=times.units,calendar=times.calendar)
times[:]=tstep


aps1[:,:,:]=tempymean
bps1[:,:,:]=tclmy
odevaps1[:,:,:]=crudevt
odevbps1[:,:,:]=clmdevt

latitudes.units = 'degrees_north'
longitudes.units = 'degrees_east'


ncfile.close()



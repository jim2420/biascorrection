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
sall= N.zeros((11,tk, 360, 720))
lall= N.zeros((11,tk, 360, 720))

clm85sall= N.zeros((11,clmday, 360, 720))
clm85lall= N.zeros((11,clmday, 360, 720))

clm85alls= N.zeros((11,tk, 360, 720))
clm85alll= N.zeros((11,tk, 360, 720))

crudevs=N.zeros((tk, 360, 720))
clmdevs=N.zeros((tk, 360, 720))
crudevl=N.zeros((tk, 360, 720))
clmdevl=N.zeros((tk, 360, 720))



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
                   tempf = base.variables["FSDS"][x,:,:]
                   sall[i,c, :, :] = tempf
                   tempg = base.variables["LWDOWN"][x,:,:]
                   lall[i,c, :, :] = tempg


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
                         temp6= base1.variables["FSDS"][ab,:,:]
                         clm85sall[i,kx, :, :] = temp6
                         temp8= base1.variables["LWDOWN"][ab,:,:]
                         clm85lall[i,kx, :, :] = temp8

                         ab=ab+1
                q1=q1+trclm
                #print q1
                base.close()
                base1.close()

#retriece 6 hourly
        for dt in range(0,tk):
            dd=dt*2
            clm85alls[i,dt,:,:]=clm85sall[i,dd, :, :]
            clm85alll[i,dt,:,:]=clm85lall[i,dd, :, :]

sempymean=N.average(sall,axis=0)
smean=N.average(sempymean,axis=0)

lempymean=N.average(lall,axis=0)
lmean=N.average(lempymean,axis=0)







sclmy=N.average(clm85alls,axis=0)
sclm=N.average(sclmy,axis=0)

lclmy=N.average(clm85alll,axis=0)
lclm=N.average(lclmy,axis=0)



crudevs=N.std(sall,axis=0)
clmdevs=N.std(clm85alls,axis=0)


crudevl=N.std(lall,axis=0)
clmdevl=N.std(clm85alll,axis=0)

ncfile=NetCDFFile('meanra45.nc','w',format='NETCDF4')
ncfile.description='clm and cruncep 2006-2016 mean 6hourly'


# dimensions

timek=ncfile.createDimension('time',tk)
latitude=ncfile.createDimension('lat', 360)
longitude=ncfile.createDimension('lon', 720)


# variables


latitudes = ncfile.createVariable('lat', 'f8', ('lat',))
longitudes = ncfile.createVariable('lon', 'f8', ('lon',))
times = ncfile.createVariable('time', 'f8', ('time',))



ass= ncfile.createVariable('crus', 'f8', ('time','lat','lon'),fill_value=-9999.)
bss= ncfile.createVariable('clm85s', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevass= ncfile.createVariable('crudevs', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevbss= ncfile.createVariable('clmdevs', 'f8', ('time','lat','lon'),fill_value=-9999.)


al= ncfile.createVariable('crul', 'f8', ('time','lat','lon'),fill_value=-9999.)
bl= ncfile.createVariable('clm85l', 'f8', ('time','lat','lon'),fill_value=-9999.)
odeval= ncfile.createVariable('crudevl', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevbl= ncfile.createVariable('clmdevl', 'f8', ('time','lat','lon'),fill_value=-9999.)




latitudes[:] = lat
longitudes[:] = lon
#dates = [datetime(2006,1,1)+n*timedelta(hours=06) for n in range(tstep.shape[0])]


#times[:]=date2num(dates,units=times.units,calendar=times.calendar)
times[:]=tstep


al[:,:,:]=lempymean
bl[:,:,:]=lclmy
odeval[:,:,:]=crudevl
odevbl[:,:,:]=clmdevl

ass[:,:,:]=sempymean
bss[:,:,:]=sclmy
odevass[:,:,:]=crudevs
odevbss[:,:,:]=clmdevs



latitudes.units = 'degrees_north'
longitudes.units = 'degrees_east'


ncfile.close()

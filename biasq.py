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


cruday=[124,112,124,120,124,120,124,124,120,124,120,124]
tk=124+112+124+120+124+120+124+124+120+124+120+124
#print tk
clmday=(124+112+124+120+124+120+124+124+120+124+120+124)*2
#print clmday
tstep=range(1,1461)
#print tstep
qall= N.zeros((11,tk, 360, 720))
pall= N.zeros((11,tk, 360, 720))

clm85qall= N.zeros((11,clmday, 360, 720))
clm85pall= N.zeros((11,clmday, 360, 720))

clm85allq= N.zeros((11,tk, 360, 720))
clm85allp= N.zeros((11,tk, 360, 720))

crudevq=N.zeros((tk, 360, 720))
clmdevq=N.zeros((tk, 360, 720))
crudevp=N.zeros((tk, 360, 720))
clmdevp=N.zeros((tk, 360, 720))






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
                   tempc = base.variables["QBOT"][x,:,:]
                   qall[i,c, :, :] = tempc
                   tempd = base.variables["PRECTmms"][x,:,:]
                   pall[i,c, :, :] = tempd


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
                         temp4= base1.variables["QBOT"][ab,:,:]
                         clm85qall[i,kx, :, :] = temp4
                         temp5= base1.variables["PRECIP"][ab,:,:]
                         clm85pall[i,kx, :, :] = temp5

                         ab=ab+1
                q1=q1+trclm
                #print q1
                base.close()
                base1.close()

#retriece 6 hourly
        for dt in range(0,tk):
            dd=dt*2
            clm85allq[i,dt,:,:]=clm85qall[i,dd, :, :]
            clm85allp[i,dt,:,:]=clm85pall[i,dd, :, :]

qempymean=N.average(qall,axis=0)
qmean=N.average(qempymean,axis=0)

pempymean=N.average(pall,axis=0)
pmean=N.average(pempymean,axis=0)


qclmy=N.average(clm85allq,axis=0)
qclm=N.average(qclmy,axis=0)

pclmy=N.average(clm85allp,axis=0)
pclm=N.average(pclmy,axis=0)

pficlm=N.zeros((tk,360,720))
pfcru=N.zeros((tk,360,720))
qc=0
de=0
for a in range(0,12):
        jj=cruday[a]
        de1=N.sum(cruday[0:a+1])
        de2=N.sum(cruday[0:a])
        print de1,de2

        #print a,jj
        for kk in range(0,jj):
                #print kk,qc
                qc=kk+de
                for x1 in range(0,360):
                                for y1 in range(0,720):
                                     pficlm[qc,x1,y1]=N.sum(pclmy[de2:de1,x1,y1])
                                     pfcru[qc,x1,y1]=N.sum(pempymean[de2:de1,x1,y1])
        de=N.sum(cruday[0:a+1])
crudevp=N.std(pall,axis=0)
clmdevp=N.std(clm85allp,axis=0)

crudevq=N.std(qall,axis=0)
clmdevq=N.std(clm85allq,axis=0)




ncfile=NetCDFFile('meanq45.nc','w',format='NETCDF4')
ncfile.description='clm and cruncep 2006-2016 mean 6hourly'

# dimensions

timek=ncfile.createDimension('time',tk)
latitude=ncfile.createDimension('lat', 360)
longitude=ncfile.createDimension('lon', 720)


# variables


latitudes = ncfile.createVariable('lat', 'f8', ('lat',))
longitudes = ncfile.createVariable('lon', 'f8', ('lon',))
times = ncfile.createVariable('time', 'f8', ('time',))


aq= ncfile.createVariable('cruq', 'f8', ('time','lat','lon'),fill_value=-9999.)
bq= ncfile.createVariable('clm85q', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevaq= ncfile.createVariable('crudevq', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevbq= ncfile.createVariable('clmdevq', 'f8', ('time','lat','lon'),fill_value=-9999.)

ap= ncfile.createVariable('crup', 'f8', ('time','lat','lon'),fill_value=-9999.)
bp= ncfile.createVariable('clm85p', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevap= ncfile.createVariable('crudevp', 'f8', ('time','lat','lon'),fill_value=-9999.)
odevbp= ncfile.createVariable('clmdevp', 'f8', ('time','lat','lon'),fill_value=-9999.)

msumpclm= ncfile.createVariable('pficlm', 'f8', ('time','lat','lon'),fill_value=-9999.)
msumpcru= ncfile.createVariable('pfcru', 'f8', ('time','lat','lon'),fill_value=-9999.)



latitudes[:] = lat
longitudes[:] = lon
#dates = [datetime(2006,1,1)+n*timedelta(hours=06) for n in range(tstep.shape[0])]
#times[:]=date2num(dates,units=times.units,calendar=times.calendar)
times[:]=tstep

msumpcru[:,:,:]=pfcru
msumpclm[:,:,:]=pficlm
aq[:,:,:]=qempymean
bq[:,:,:]=qclmy
odevaq[:,:,:]=crudevq
odevbq[:,:,:]=clmdevq

ap[:,:,:]=pempymean
bp[:,:,:]=pclmy
odevap[:,:,:]=crudevp
odevbp[:,:,:]=clmdevp



latitudes.units = 'degrees_north'
longitudes.units = 'degrees_east'


ncfile.close()






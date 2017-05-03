#!/usr/bin/env python
#from mpi4py import MPI
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
cruday=[124,112,124,120,124,120,124,124,120,124,120,124]
clmday=(124+112+124+120+124+120+124+124+120+124+120+124)*2
#print clmday
tstep=range(1,1461)
#print tstep

base = NetCDFFile ("/scratch2/scratchdirs/tslin2/bias/meant45.nc", mode='r')
baseu = NetCDFFile ("/scratch2/scratchdirs/tslin2/bias/meanws45.nc", mode='r')
basera = NetCDFFile ("/scratch2/scratchdirs/tslin2/bias/meanra45.nc", mode='r')
baseps = NetCDFFile ("/scratch2/scratchdirs/tslin2/bias/meanps45.nc", mode='r')
baseq = NetCDFFile ("/scratch2/scratchdirs/tslin2/bias/meanq45.nc", mode='r')

lon = base.variables["lon"][:]
lat = base.variables["lat"][:]
time  = base.variables["time"][:]

tall= base.variables["cru"][:]
uall= baseu.variables["cruu"][:]
qall= baseq.variables["cruq"][:]
pall= baseq.variables["pfcru"][:]
sall= basera.variables["crus"][:]
lall= basera.variables["crul"][:]
psall= baseps.variables["crups"][:]

clm85tall= base.variables["clm85"][:]
clm85uall= baseu.variables["clm85u"][:]
clm85qall= baseq.variables["clm85q"][:]
clm85pall= baseq.variables["pficlm"][:]
clm85psall= baseps.variables["clm85ps"][:]
clm85sall= basera.variables["clm85s"][:]
clm85lall= basera.variables["clm85l"][:]

crudevt=base.variables["crudev"][:]
clmdevt=base.variables["clmdev"][:]
crudevu=baseu.variables["crudevu"][:]
clmdevu=baseu.variables["clmdevu"][:]
crudevq=baseq.variables["crudevq"][:]
clmdevq=baseq.variables["clmdevq"][:]
crudevps=baseps.variables["crudevps"][:]
clmdevps=baseps.variables["clmdevps"][:]
crudevs=basera.variables["crudevs"][:]
clmdevs=basera.variables["clmdevs"][:]
crudevl=basera.variables["crudevl"][:]
clmdevl=basera.variables["clmdevl"][:]



years = range(2009, 2011)
months= range(1,13)
q=0
c=0
kx=0

for i, year in enumerate(years):
        q=0
        q1=0
        kc=0
        num=0
        for m, month in enumerate(months):
                if month < 10 :
                        month = '0{}'.format(month)

                #print "/scratch2/scratchdirs/tslin2/bias/data/clm85/{0}_{1}.nc".format(year,month)
                base1 = NetCDFFile ("/scratch2/scratchdirs/tslin2/bias/data/clm45/{0}_{1}.nc".format(year,month), mode='r')
                clmtft= base1.variables["TBOT"][:]
                clmtfu= base1.variables["UWIND"][:]
                clmtfq= base1.variables["QBOT"][:]
                clmtfps= base1.variables["PRSF"][:]
                clmtfp= base1.variables["PRECIP"][:]
                clmtfl= base1.variables["LWDOWN"][:]
                clmtfs= base1.variables["FSDS"][:]
                cct=cruday[num]
                clm85t=N.zeros((cct,360,720))
                clm85u=N.zeros((cct,360,720))
                clm85p=N.zeros((cct,360,720))
                clm85ps=N.zeros((cct,360,720))
                clm85s=N.zeros((cct,360,720))
                clm85l=N.zeros((cct,360,720))
                clm85q=N.zeros((cct,360,720))


                for dt in range(0,cct):
                    dd=dt*2
                    gg=dt+kc
                    #print gg
                    clm85t[dt,:,:]=tall[gg,:,:]+((crudevt[gg,:,:]/clmdevt[gg,:,:])*(clmtft[dd, :, :]-clm85tall[gg,:,:]))
                    clm85q[dt,:,:]=qall[gg,:,:]+((crudevq[gg,:,:]/clmdevq[gg,:,:])*(clmtfq[dd, :, :]-clm85qall[gg,:,:]))
                    clm85ps[dt,:,:]=psall[gg,:,:]+((crudevps[gg,:,:]/clmdevps[gg,:,:])*(clmtfps[dd, :, :]-clm85psall[gg,:,:]))
                    clm85s[dt,:,:]=sall[gg,:,:]+((crudevs[gg,:,:]/clmdevs[gg,:,:])*(clmtfs[dd, :, :]-clm85sall[gg,:,:]))
                    clm85l[dt,:,:]=lall[gg,:,:]+((crudevl[gg,:,:]/clmdevl[gg,:,:])*(clmtfl[dd, :, :]-clm85lall[gg,:,:]))

                    clm85u[dt,:,:]=uall[gg,:,:]+(((crudevu[gg,:,:]/clmdevu[gg,:,:])*(N.sqrt(2*clmtfu[dd, :, :]*clmtfu[dd, :, :])-clm85uall[gg,:,:])))
                    clm85p[dt,:,:]=clmtfp[dd, :, :]*pall[gg,:,:]/clm85pall[gg,:,:]
                kc=N.sum(cruday[0:num+1])
                base1.close()
                num=num+1

                ncfile=NetCDFFile("/scratch2/scratchdirs/tslin2/bias/output/clm45/{0}_{1}.nc".format(year,month),'w',format='NETCDF3_64BIT_OFFSET')
                ncfile.description='clm rcp 4.5 6hourly'
# dimensions
                timek=ncfile.createDimension('time',cct)
                latitude=ncfile.createDimension('latitude', 360)
                longitude=ncfile.createDimension('longitude', 720)
# variables



                latitudes = ncfile.createVariable('latitude', 'f8', ('latitude',))
                longitudes = ncfile.createVariable('longitude', 'f8', ('longitude',))
                times = ncfile.createVariable('time', 'f8', ('time',))

                at= ncfile.createVariable('TBOT', 'f8', ('time','latitude','longitude'),fill_value=-9999.)
                bt= ncfile.createVariable('UWIND', 'f8', ('time','latitude','longitude'),fill_value=-9999.)
                odevat= ncfile.createVariable('VWIND', 'f8', ('time','latitude','longitude'),fill_value=-9999.)
                odevbt= ncfile.createVariable('QBOT', 'f8', ('time','latitude','longitude'),fill_value=-9999.)

                aps= ncfile.createVariable('PRECIP', 'f8', ('time','latitude','longitude'),fill_value=-9999.)
                bps= ncfile.createVariable('FSDS', 'f8', ('time','latitude','longitude'),fill_value=-9999.)
                odevaps= ncfile.createVariable('PRSF', 'f8', ('time','latitude','longitude'),fill_value=-9999.)
                odevbps= ncfile.createVariable('LWDOWN', 'f8', ('time','latitude','longitude'),fill_value=-9999.)

                latitudes[:] = lat
                longitudes[:] = lon
                times[:]=range(0,cct)
                at[:,:,:]=clm85t
                bt[:,:,:]=clm85u/(N.sqrt(2))
                odevat[:,:,:]=clm85u/(N.sqrt(2))
                odevbt[:,:,:]=clm85q
                aps[:,:,:]=clm85p
                bps[:,:,:]=clm85s
                odevaps[:,:,:]=clm85ps
                odevbps[:,:,:]=clm85l


                latitudes.long_name = "latitude"
                latitudes.units = 'degrees_north'
                longitudes.long_name = "longitude"
                longitudes.units = 'degrees_east'
                at.long_name = "averaged temperature"
                at.units = 'K'
                bt.long_name = "u wind"
                bt.units = 'm/s'
                odevat.long_name = "v wind"
                odevat.units = 'm/s'
                odevbt.long_name = "humidity"
                odevbt.units = 'g/g'
                aps.long_name = "precipitation rate"
                aps.units = 'mm/s'
                bps.long_name = "shortwave radiation"
                bps.units = 'W/m2'
                odevaps.long_name = "surface atmosphere pressure"
                odevaps.units = 'Pa'
                odevbps.long_name = "incoming longwave radiation"
                odevbps.units = 'W/m2'
                odevbps.missing_value = -9999.
                odevaps.missing_value = -9999.
                aps.missing_value = -9999.
                bps.missing_value = -9999.
                bt.missing_value = -9999.
                at.missing_value = -9999.
                odevat.missing_value = -9999.
                odevbt.missing_value = -9999.

                ncfile.close()

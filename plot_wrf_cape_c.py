#!/usr/bin/python


import Nio as netcdf
import numpy as np
import math
import matplotlib
matplotlib.use('agg')
import pylab
import os
import sys,getopt

from compute_stability import SBStability

from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap

# Set the default domain to be d01
dom = 'd01'
var = 'all'
export_flag = 0



# Set up a command-line argument structure to allow
# for command-line changes of variables.
# f --> the name of the domain we want to use
(opts,args)=getopt.getopt(sys.argv[1:],'f:v:e')
for o,a in opts:
	if o=="-f":
		dom = a
	if o=="-v":
		var = str(a)
	if o=="-e":
		export_flag = 1	

	

# Skip is the length between outputs
# skip = 1
if dom == 'd02':
	skip = 1
	DP_CLEVS = range(55,90,5)
else:
	skip = 3
	DP_CLEVS = range(55,90,5)


#filename = '../wrfout_' + dom + '.nc'
filename = '/home/wrfuser/DOMAINS/TSupinie/wrfout_d01_2007-05-04_12:00:00.nc'
nc = netcdf.open_file(filename)


PHB = nc.variables['PHB']
PH  = nc.variables['PH']
PB  = nc.variables['PB']
P   = nc.variables['P']
T  = nc.variables['T']
TB = 300
QVAP = nc.variables['QVAPOR']
PSFC = nc.variables['PSFC']
HGT = nc.variables['HGT']
Q2 = nc.variables['Q2']
T2 = nc.variables['T2']
times = nc.variables['Times']

max_metgrid_level = 5


x_dim  = len(nc.variables['XLONG'][0,0,:])
y_dim  = len(nc.variables['XLAT'][0,:,0])
z_dim = 39
time = 0
# Set grid spacing
dx = 12000
dy = 12000


# Thin factor is used for thinning out wind barbs
thin = 10

# x_dim and y_dim are the x and y dimensions of the model
# domain in gridpoints
x_dim = len(nc.variables['XLAT'][0,0,:])
y_dim = len(nc.variables['XLONG'][0,:,0])


# Central latitude and longitude are found by taking the
# XLAT and XLON at the points halfway across the x and y dims
lon_ctr = nc.variables['XLONG'][0,0,int(x_dim*.5)]
lat_ctr = nc.variables['XLAT'][0,int(y_dim*.5),0]



# Set up MANUAL dimensions for the image size
# man_x_ll = 0
# man_y_ll = 0
# man_x_ur = x_dim
# man_y_ur = y_dim

man_x_ll = 0
man_y_ll = 75
man_x_ur = x_dim - 54
man_y_ur = y_dim - 365


# Lower left and upper right lats and lons are found by 
# getting the lat and lon at the lower left corner (0,0)
# of the domain and the upper right (x_dim, y_dim) corner
ll_lat = np.min(nc.variables['XLAT'][0,:,man_y_ll])
ll_lon = np.min(nc.variables['XLONG'][0,man_x_ll,:])
ur_lat = np.max(nc.variables['XLAT'][0,:,man_y_ur])
ur_lon = np.max(nc.variables['XLONG'][0,man_x_ur,:])




# Draw the base map behind it with the lats and
# lons calculated earlier
map = Basemap(resolution='i',projection='lcc',\
	llcrnrlon= ll_lon, llcrnrlat=ll_lat,\
	urcrnrlon= ur_lon, urcrnrlat= ur_lat,\
	lat_0=lat_ctr,lon_0=lon_ctr,lat_1=38)

# This sets the standard grid point structure at full resolution
x,y = map(nc.variables['XLONG'][0],nc.variables['XLAT'][0])

# This sets a thinn-ed out grid point structure for plotting
# wind barbs at the interval specified in "thin"
x_th,y_th = map(nc.variables['XLONG'][0,::thin,::thin],\
	nc.variables['XLAT'][0,::thin,::thin])


# This sets a subgrid for voricity and other dx/dy calculations
x_1,y_1 = map(nc.variables['XLONG'][0,1:(y_dim-1),1:(x_dim-1)],\
	nc.variables['XLAT'][0,1:(y_dim-1),1:(x_dim-1)])



# Set universal figure margins
width = 10
height = 8

pylab.figure(figsize=(width,height))
pylab.rc("figure.subplot", left = .001)
pylab.rc("figure.subplot", right = .999)
pylab.rc("figure.subplot", bottom = .001)
pylab.rc("figure.subplot", top = .999)





def timestring(wrftime,curtime):
	curtime_str = '%02d' % curtime
	year  = str(wrftime[2])  + str(wrftime[3])
	month = str(wrftime[5])  + str(wrftime[6])
	day   = str(wrftime[8])  + str(wrftime[9])
	hour  = str(wrftime[11]) + str(wrftime[12])
	outtime = year + month + day + '/' + hour + '00Z F'+ curtime_str
	return outtime

def drawmap(DATA,TITLESTRING,PROD):
	# rect = (.25,.05,.5,.05)
	# pylab.add_axis(rect,axisbg='white',alpha=.5)
	pylab.colorbar(DATA,orientation='horizontal',\
		extend='both',aspect=65,\
		shrink=.875,pad=0)
	map.drawstates(color='k')
	map.drawcoastlines(color='k')
	map.drawcountries(color='k')
	pylab.title('OWL/KAHOOLAWE WRF-ARW %s   Valid: %s' % (TITLESTRING, curtimestring), \
		fontsize=11,bbox=dict(facecolor='white', alpha=0.6),\
		x=0.5,y=.95,weight = 'demibold',style='oblique', \
		stretch='normal', family='sans-serif')
	file_id = '%s_%s_f%02d' % (dom, PROD, time*skip)
	filename = '%s.png' % (file_id)
	
	pylab.savefig(filename)
	pylab.close()

	if export_flag == 1:
		os.system('convert %s %s.gif' % (filename, file_id))
		os.system('rm -f %s' % filename)





# Constants

grav = 9.81

# Epsilons for moisture
ezero = 6.112

eps = .622

cp = 1004
rgas = 287.04

gamma = (rgas/cp)

# cp_moist = cp*(1 + cpmd*qvp)
cpmd = .887

# rgas_moist = rgas * (1+rgasmd*qvp)
rgasmd = .608

g = 9.81

def calc_T(theta_prime,press):
	theta = np.add(theta_prime,300)
	temp = np.multiply(theta,np.power(np.divide(press,100000.),(287.04/1004)))
	return temp

def plot_cape():
	hgt = np.divide(np.add(PHB[time],PH[time]),g)
	th   = T[time]
	p = PB[time]+P[time]
	t = calc_T(th,p)
	qvp = QVAP[time]

	cape = SBStability(t, qvp, p / 100)

	print("    SBCAPE")
	print np.shape(cape)
#	print np.shape(li)
	#Set Figure Size (1000 x 800)
	pylab.figure(figsize=(width,height),frameon=False)
	cdict_vort ={'red':	((0.00, 1.00, 1.00),
				(0.10, 1.00, 1.00),
				(0.25, 1.00, 1.00),
				(0.35, 0.82, 0.82),
				(0.40, 0.63, 0.63),
				(0.45, 0.00, 0.00),
				(0.50, 0.12, 0.12),
				(0.60, 0.00, 0.00),
				(0.70, 0.25, 0.25),
				(0.85, 0.50, 0.50),
				(1.00, 0.65, 0.65)),

		'green':	((0.00, 0.43, 0.43),
				(0.10, 0.88, 0.88),
				(0.25, 1.00, 1.00),
				(0.35, 0.96, 0.96),
				(0.40, 0.82, 0.82),
				(0.45, 0.75, 0.75),
				(0.50, 0.56, 0.56),
				(0.60, 0.41, 0.41),
				(0.70, 0.00, 0.00),
				(0.85, 0.00, 0.00),
				(1.00, 0.13, 0.13)),

		'blue':		((0.00, 0.00, 0.00),
				(0.10, 0.20, 0.20),
				(0.25, 1.00, 1.00),
				(0.35, 1.00, 1.00),
				(0.40, 1.00, 1.00),
				(0.45, 1.00, 1.00),
				(0.50, 1.00, 1.00),
				(0.60, 0.88, 0.88),
				(0.70, 0.80, 0.80),
				(0.85, 0.59, 0.59),
				(1.00, 0.94, 0.94))}
	
	vort_coltbl = LinearSegmentedColormap('VORT_COLTBL',cdict_vort)
	
	CAPE=pylab.contourf(x,y,cape[0],cape_clevs)

	title = 'CAPE (J/kg)'
	prodid = 'cape'

	drawmap(CAPE, title, prodid)

	LI=pylab.contourf(x,y,cape[1],li_clevs)

	title = '500mb LI (K)'
	prodid = 'li'

	drawmap(LI, title, prodid)



		
# Set contour levels
cape_clevs = range(500,5500,500)
li_clevs = range(-12,0,2)


# Check to see if we are exporting
if export_flag == 1:
	dom = 'wrf'
# Begin looping through times
for time in range(4,len(nc.variables['XLAT'])):
	print 'Plotting time ',time*skip
	curtimestring = timestring(times[time],time*skip)
	if var == 'vort':
		for level in range(1,3):
			plot_ua_avort()
	if var == 'wind':
		for level in range(0,5):
			plot_ua_winds()
	if var == 'vvel':
		for level in range(0,3):
			plot_ua_vvel()
	else:	
		plot_cape()

	os.system('scp d01_cape_f12.png tsupinie@somclass01.som.nor.ouint:~/http/.')
	os.system('scp d01_li_f12.png tsupinie@somclass01.som.nor.ouint:~/http/.')
	break

# Copy the files over to the appropriate locations on HOOT

#if export_flag == 1:
#	os.system('scp wrf_*mb_wind*.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')
#	os.system('scp wrf_*mb_temp*.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')
#	os.system('scp wrf_*mb_vert_vel*.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')
#	os.system('scp wrf_*mb_vort_*.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')
#
#
#else:
#	os.system('scp d01_li_* hoot@10.197.1.220:/usr/home/hoot/http/wrf_data/SPlains_d01/.')



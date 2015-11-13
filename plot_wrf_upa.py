#!/usr/bin/python

# Working script to generate maps from wrfout netCDF files
# using matplot lib with basemap
# Basemap coding from David John Gagne II
# Written by Luke Madaus for use with operational WRF domains

import Nio as netcdf
import matplotlib
matplotlib.use('agg')
import pylab
import numpy as np
import os
import sys,getopt
import math
import datetime
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap

# Set the default domain to be d01
dom = 'd01'
var = 'all'
export_flag = 0

restart_time = 0

# Set up a command-line argument structure to allow
# for command-line changes of variables.
# f --> the name of the domain we want to use
(opts,args)=getopt.getopt(sys.argv[1:],'f:v:r:e')
for o,a in opts:
	if o=="-f":
		dom = a
	if o=="-v":
		var = str(a)
	if o=="-e":
		export_flag = 1	
	if o=="-r":
		restart_time = int(a)
	

# Skip is the length between outputs
# skip = 1
if dom == 'd02':
	skip = 1
	DP_CLEVS = range(55,90,5)
else:
	skip = 3
	DP_CLEVS = range(55,90,5)


filename = '../wrfout_' + dom + '_PLEV.nc'
nc = netcdf.open_file(filename)


# Grab three variables for now
#temps_ua = nc.variables['TT']
temps_ua =  nc.variables['T']
#temps_base_ua = nc.variables['TB']
temps_base_ua = 300
u_wind_ms_ua = nc.variables['UU'][:,:,:,:-1]
v_wind_ms_ua = nc.variables['VV'][:,:,:-1,:]
w_wind_ms_ua = nc.variables['W']
#rhum_ua = nc.variables['RH']
qvap = nc.variables['QVAPOR']
#ght = nc.variables['GHT']
phb = nc.variables['PHB']
ph = nc.variables['PH']
w_wind_ua = nc.variables['W']
f=nc.variables['F']
times = nc.variables['Times']

# Set pressure levels
plev = [850, 700, 500, 300, 250]


# For just those levels in wrfout_PLEV file
#levid = {'850' : 0,
#	'700'  : 1,
#	'500'  : 2,
#	'300'  : 3,
#	'250'  : 4}

# For all levels in wrfout_PLEV file
levid = {'850' : 9,
	'700'  : 13,
	'500'  : 17,
	'300'  : 21,
	'250'  : 22}


contlevid = {'850' : 0,
	'700'  : 1,
	'500'  : 2,
	'300'  : 3,
	'250'  : 4}


dayofweek = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']

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





# Set contouring ranges

hght_clevs=[]
wspd_clevs=[]
temp_clevs=[]
rhum_clevs = [60,70,80,90]

hghts_850=range(1200,1620,15)
wspds_850=range(20,85,5)
temps_850=range(-45,45,5)
hght_clevs.append(hghts_850)
wspd_clevs.append(wspds_850)
temp_clevs.append(temps_850)

hghts_700=range(2700,3500,15)
wspds_700=range(20,85,5)
temps_700=range(-45,45,5)
hght_clevs.append(hghts_700)
wspd_clevs.append(wspds_700)
temp_clevs.append(temps_700)

hghts_500=range(5100,6000,30)
wspds_500=range(50,125,5)
temps_500=range(-55,25,5)
hght_clevs.append(hghts_500)
wspd_clevs.append(wspds_500)
temp_clevs.append(temps_500)


hghts_300=range(8520,9700,30)
wspds_300=range(60,220,10)
temps_300=range(-45,45,5)
hght_clevs.append(hghts_300)
wspd_clevs.append(wspds_300)
temp_clevs.append(temps_300)


hghts_250=range(9900,11000,30)
wspds_250=range(60,220,10)
temps_250=range(-45,45,5)
hght_clevs.append(hghts_250)
wspd_clevs.append(wspds_250)
temp_clevs.append(temps_250)

# For Vorticity
#neg_vorts = range(-5,1,1)
vort_clevs = range(-13,45,2)

# For Vertical Velocity
vvel_clevs = range(-22,36,2)


def timestring(wrftime,curtime):
	curtime_str = '%02d' % curtime
	year  = str(wrftime[2])  + str(wrftime[3])
	month = str(wrftime[5])  + str(wrftime[6])
	day   = str(wrftime[8])  + str(wrftime[9])
	hour  = str(wrftime[11]) + str(wrftime[12])
	numdow = datetime.date.weekday(datetime.date(int(year), int(month), int(day)))
	dow  = dayofweek[numdow]
	outtime = dow + ' ' + year + month + day + '/' + hour + '00Z F'+ curtime_str
	return outtime

def drawmap(DATA,TITLESTRING,PROD,UNITS):
	# rect = (.25,.05,.5,.05)
	# pylab.add_axis(rect,axisbg='white',alpha=.5)
	pylab.colorbar(DATA,orientation='horizontal',\
		extend='both',aspect=65,\
		shrink=.875,pad=0)
	map.drawstates(color='k')
	map.drawcoastlines(color='k')
	map.drawcountries(color='k')
	#map.bluemarble()
	pylab.suptitle('%s' % UNITS, fontsize = 11, x = 0.08, y = 0.105)
	pylab.title('OWL/KHLWE WRF-ARW %s   Valid: %s' % (TITLESTRING, curtimestring), \
		fontsize=11,bbox=dict(facecolor='white', alpha=0.6),\
		x=0.5,y=.95,weight = 'demibold',style='oblique', \
		stretch='normal', family='sans-serif')
	file_id = '%s_%s_f%02d' % (dom, PROD, time*skip+restart_time)
	filename = '%s.png' % (file_id)
	
	pylab.savefig(filename)
	pylab.close()

	if export_flag == 1:
		os.system('convert -render -flatten %s %s.gif' % (filename, file_id))
		os.system('rm -f %s' % filename)
		# os.system('scp %s.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.' % file_id)

def plot_ua_winds(press):
	from scipy.ndimage.filters import gaussian_filter
	print("    %smb WINDS" % press)
	level = levid[str(press)]

	#Set Figure Size (1000 x 800)
	pylab.figure(figsize=(width,height),frameon=False)
	heights=hght_clevs[contlevid[str(press)]]
	
	phtotal = np.add(phb[time,level],ph[time,level])
	gheight = gaussian_filter(np.divide(phtotal, 9.81),1)
	#gheight = np.divide(phtotal, 9.81)
	#gheight = ght[time,level]

	# Convert winds from m/s to kts and then draw barbs	
	u_wind_kts = u_wind_ms_ua[time,level] * 1.94384449
	v_wind_kts = v_wind_ms_ua[time,level] * 1.94384449
	
	u_wind_kts_sq = np.power(u_wind_kts, 2)
	v_wind_kts_sq = np.power(v_wind_kts, 2)
	wind_sum = np.add(u_wind_kts_sq,v_wind_kts_sq)
	# Find actual wind speed and contour
	windspd = np.sqrt(wind_sum)
	WSPD=pylab.contourf(x,y,windspd,wspd_clevs[contlevid[str(press)]],extend='max')

	# Contour the heights

	HGHT=pylab.contour(x,y,gheight,heights,linewidths=1.5, colors = 'k')
	pylab.clabel(HGHT,inline=1,fontsize=8,fmt='%1.0f',inline_spacing=1)

	#pylab.clabel(T,inline=1,fontsize=10)


	pylab.barbs(x_th,y_th,u_wind_kts[::thin,::thin],\
		v_wind_kts[::thin,::thin], length=5, sizes={'spacing':0.2},pivot='middle')

	title = '%s mb Height (m), Wind' % press
	prodid = '%smb_wind' % press
	units = 'kts'

	drawmap(WSPD, title, prodid, units)

def plot_ua_temps(press):
	print("    %smb TEMPERATURE" % press)
	level = levid[str(press)]
	#Set Figure Size (1000 x 800)
	pylab.figure(figsize=(width,height),frameon=False)
	heights=hght_clevs[contlevid[str(press)]]


	# Temperature calculation here
	
	#TC = temps_ua[time,level] - 273
	TH_K = np.add(temps_ua[time,level],temps_base_ua)
	TK = np.multiply(TH_K,math.pow((float(press)/1000.),(287.04/1004.))) 
	TC = TK - 273
	TEMP=pylab.contourf(x,y,TC,temp_clevs[contlevid[str(press)]],extend='both')

	# Contour the heights

	phtotal = np.add(phb[time,level],ph[time,level])
	gheight = np.divide(phtotal, 9.81)
	#gheight = ght[time,level]
	HGHT=pylab.contour(x,y,gheight,heights,colors='k',linewidths=1.5)
	pylab.clabel(HGHT,inline=1,fontsize=8,fmt='%1.0f',inline_spacing=1)


	# Convert winds from m/s to kts and then draw barbs	
	u_wind_kts = u_wind_ms_ua[time,level] * 1.94384449
	v_wind_kts = v_wind_ms_ua[time,level] * 1.94384449
	

	pylab.barbs(x_th,y_th,u_wind_kts[::thin,::thin],\
		v_wind_kts[::thin,::thin], length=5, sizes={'spacing':0.2},pivot='middle')


	title = '%s mb Temp, Height (m), Wind (kts)' % press
	prodid = '%smb_temp' % press
	units = u"\u00B0" + "C"

	drawmap(TEMP, title, prodid, units)

def plot_ua_rhum(press):
	print("    %smb RHUM" % press)
	level = levid[str(press)]

	#Set Figure Size (1000 x 800)
	pylab.figure(figsize=(width,height),frameon=False)
	heights=hght_clevs[contlevid[str(press)]]


	cdict =	{'red':		((0.00, 0.76, 0.76),
				(0.25, 0.64, 0.64),
				(0.50, 0.52, 0.52),
				(0.75, 0.42, 0.42),
				(1.00, 0.32, 0.32)),

		'green':	((0.00, 0.90, 0.90),
				(0.25, 0.80, 0.80),
				(0.50, 0.70, 0.70),
				(0.75, 0.60, 0.60),
				(1.00, 0.50, 0.50)),

		'blue':		((0.00, 0.49, 0.49),
				(0.25, 0.32, 0.32),
				(0.50, 0.17, 0.17),
				(0.75, 0.06, 0.06),
				(1.00, 0.05, 0.05))}

	 
	rhum_coltbl = LinearSegmentedColormap('RHUM_COLTBL',cdict)

	# Temperature calculation here
	
	#TC = temps_ua[time,level] - 273
	TH_K = np.add(temps_ua[time,level],temps_base_ua)
	TK = np.multiply(TH_K,math.pow((float(press)/1000.),(287.04/1004.))) 
	TC = TK - 273

	# RHUM Calculation here
	es = 6.112*np.exp(17.67*TC/(TC+243.5))
	w = np.divide(qvap[time,level],(1-qvap[time,level]))
	e = (w*float(press))/(0.622+w)
	relh = e / es * 100.

	RHUM=pylab.contourf(x,y,relh,rhum_clevs,extend='max',cmap=rhum_coltbl)
	
	# Contour the heights

	phtotal = np.add(phb[time,level],ph[time,level])
	gheight = np.divide(phtotal, 9.81)
	#gheight = ght[time,level]
	HGHT=pylab.contour(x,y,gheight,heights,colors='k',linewidths=1.5)
	pylab.clabel(HGHT,inline=1,fontsize=8,fmt='%1.0f',inline_spacing=1)


	# Convert winds from m/s to kts and then draw barbs	
	u_wind_kts = u_wind_ms_ua[time,level] * 1.94384449
	v_wind_kts = v_wind_ms_ua[time,level] * 1.94384449
	

	pylab.barbs(x_th,y_th,u_wind_kts[::thin,::thin],\
		v_wind_kts[::thin,::thin], length=5, sizes={'spacing':0.2},pivot='middle')

	title = '%s mb RELH, Height (m), Wind (kts)' % press
	prodid = '%smb_rel_hum' % press
	units = "percent"

	drawmap(RHUM, title, prodid, units)


def plot_ua_vvel(press):
	print("    %smb VVEL" % press)
	level = levid[str(press)]

	#Set Figure Size (1000 x 800)
	pylab.figure(figsize=(width,height),frameon=False)
	
	cdict_vvel ={'red':	((0.00, 0.50, 0.50),
				(0.10, 1.00, 1.00),
				(0.20, 1.00, 1.00),
				(0.35, 1.00, 1.00),
				(0.40, 1.00, 1.00),
				(0.45, 0.80, 0.80),
				(0.50, 0.72, 0.72),
				(0.55, 0.68, 0.68),
				(0.60, 0.64, 0.64),
				(0.65, 0.60, 0.70),
				(0.70, 0.40, 0.45),
				(0.80, 0.20, 0.20),
				(0.90, 0.00, 0.00),
				(1.00, 0.00, 0.00)),

		'green':	((0.00, 0.35, 0.35),
				(0.10, 0.45, 0.45),
				(0.20, 0.55, 0.55),
				(0.35, 1.00, 1.00),
				(0.40, 1.00, 1.00),
				(0.45, 1.00, 1.00),
				(0.50, 1.00, 1.00),
				(0.55, 1.00, 1.00),
				(0.60, 1.00, 1.00),
				(0.65, 1.00, 1.00),
				(0.70, 1.00, 1.00),
				(0.80, 1.00, 1.00),
				(0.90, 1.00, 1.00),
				(1.00, 0.80, 0.80)),

		'blue':		((0.00, 0.00, 0.00),
				(0.10, 0.00, 0.00),
				(0.20, 0.20, 0.20),
				(0.35, 1.00, 1.00),
				(0.40, 1.00, 1.00),
				(0.45, 1.00, 1.00),
				(0.50, 1.00, 1.00),
				(0.55, 1.00, 1.00),
				(0.60, 1.00, 1.00),
				(0.65, 1.00, 1.00),
				(0.70, 1.00, 1.00),
				(0.80, 1.00, 1.00),
				(0.90, 1.00, 1.00),
				(1.00, 0.80, 0.80))}

	vvel_coltbl = LinearSegmentedColormap('VVEL_COLTBL',cdict_vvel)

	vert_vel = np.multiply(w_wind_ua[time,level],100)	
	vert_vel = np.nan_to_num(vert_vel)

	VVEL=pylab.contourf(x,y,vert_vel,vvel_clevs,cmap=vvel_coltbl,extend='both')

	# Contour the heights

	heights=hght_clevs[contlevid[str(press)]]
	phtotal = np.add(phb[time,level],ph[time,level])
	gheight = np.divide(phtotal, 9.81)
	#gheight = ght[time,level]

	HGHT=pylab.contour(x,y,gheight,heights,colors='k',linewidths=1.5)
	pylab.clabel(HGHT,inline=1,fontsize=8,fmt='%1.0f',inline_spacing=1)


	# Convert winds from m/s to kts and then draw barbs	
	u_wind_kts = u_wind_ms_ua[time,level] * 1.94384449
	v_wind_kts = v_wind_ms_ua[time,level] * 1.94384449


	pylab.barbs(x_th,y_th,u_wind_kts[::thin,::thin],\
		v_wind_kts[::thin,::thin], length=5, sizes={'spacing':0.2},pivot='middle')


	title = '%s mb VVel, Height (m), Wind (kts)' % press
	prodid = '%smb_vert_vel' % press
	units = "cm/s"

	drawmap(VVEL, title, prodid, units)

def plot_ua_avort(press):
	print("    %smb VORICITY" % press)
	level = levid[str(press)]
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

	# Vorticity calculation goes here
	dx_2 = 2 * dx
	dy_2 = 2 * dy
	for xs in range(x_dim - 2):
		cur_column = []
		for ys in range(y_dim - 2):
			du = np.subtract(u_wind_ms_ua[time,level,(ys+1),(xs+2)],u_wind_ms_ua[time,level,(ys+1),xs]) 
			dv = np.subtract(v_wind_ms_ua[time,level,(ys+2),(xs+1)],v_wind_ms_ua[time,level,ys,(xs+1)])
			f_val = f[time,(ys+1),(xs+1)]
			cur_avort = (dv/dx_2) - (du/dy_2) + f_val
			cur_column.append(cur_avort)
			np_column=np.array(cur_column)

		if xs == 0:
			vort_grid=np_column
		else:
			vort_grid = np.column_stack((vort_grid,np_column)) 
			
	avort = np.multiply(vort_grid, 100000)
	#np.reshape(avort,-1)
	avort = np.nan_to_num(avort)
	VORT=pylab.contourf(x_1,y_1,avort,vort_clevs,cmap=vort_coltbl, extend='both')

	# Contour the heights

	heights=hght_clevs[contlevid[str(press)]]
	phtotal = np.add(phb[time,level],ph[time,level])
	gheight = np.divide(phtotal, 9.81)
	#gheight = ght[time,level]

	HGHT=pylab.contour(x,y,gheight,heights,colors='k',linewidths=1.5)
	pylab.clabel(HGHT,inline=1,fontsize=8,fmt='%1.0f',inline_spacing=1)


	# Convert winds from m/s to kts and then draw barbs	
	u_wind_kts = u_wind_ms_ua[time,level] * 1.94384449
	v_wind_kts = v_wind_ms_ua[time,level] * 1.94384449


	pylab.barbs(x_th,y_th,u_wind_kts[::thin,::thin],\
		v_wind_kts[::thin,::thin], length=5, sizes={'spacing':0.2},pivot='middle')


	title = '%s mb avort (10^-5 s-1) Height, Wind (kts)' % press
	prodid = '%smb_vort' % press
	units = "10^-5 / s"

	drawmap(VORT, title, prodid, units)

def plot_LI():
	x = 1	


# Check to see if we are exporting
if export_flag == 1:
	dom = 'wrf'
# Begin looping through times
for time in range(0,len(nc.variables['T'])):
	print 'Plotting time ',time*skip+restart_time
	curtimestring = timestring(times[time],time*skip+restart_time)
	if var == 'vort':
		for level in range(1,3):
			try:
				plot_ua_avort(prlev)
			except:
				x=1
	if var == 'wind':
		for level in range(0,5):
			plot_ua_winds(prlev)
	if var == 'vvel':
		for level in range(0,3):
			plot_ua_vvel(prlev)
	if var == 'temp':
		for level in range(0,3):
			plot_ua_temps(prlev)
	if var == 'rhum':
		for level in range(0,3):
			plot_ua_rhum(prlev)
	if var == 'all':	
		for prlev in plev:
			plot_ua_winds(prlev)
			if (levid[str(prlev)] <= levid['500']):
				plot_ua_temps(prlev) 
				plot_ua_vvel(prlev)
				plot_ua_rhum(prlev)
				if (levid[str(prlev)] >= levid['700']):
					plot_ua_avort(prlev)

# Copy the files over to the appropriate locations on HOOT
if export_flag == 1:
	os.system('scp wrf_*mb_wind*.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')
	os.system('scp wrf_*mb_temp*.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')
	os.system('scp wrf_*mb_rel_hum*.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')
	os.system('scp wrf_*mb_vert_vel*.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')
	os.system('scp wrf_*mb_vort_*.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')


else:
	os.system('scp d01_*_wind_* hoot@10.197.1.220:/usr/home/hoot/http/wrf_data/SPlains_d01/500mb/.')
	os.system('scp d01_*_temp_* hoot@10.197.1.220:/usr/home/hoot/http/wrf_data/SPlains_d01/500mb/.')
	os.system('scp d01_*_vert_vel_* hoot@10.197.1.220:/usr/home/hoot/http/wrf_data/SPlains_d01/500mb/.')
	os.system('scp d01_*_vort_* hoot@10.197.1.220:/usr/home/hoot/http/wrf_data/SPlains_d01/500mb/.')

#!/usr/bin/python

# Working script to generate maps from wrfout netCDF files
# using matplot lib with basemap
# Original code by David John Gagne II
# Additions with mesonet by Luke Madaus
# Some mesonet acquisition code from Dr. Brian Fiedler
import sys
import urllib2
import Nio
import matplotlib
matplotlib.use('agg')
import pylab
import math
import numpy as np
import os, getopt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

cursite = 'OUN'

(opts,args) = getopt.getopt(sys.argv[1:],'s:')
for o,a in opts:
	if o == '-s':
		cursite = str(a)



filename = '../wrfout_d01.nc'
nc = Nio.open_file(filename)
#print nc.variables.keys()


sfcT = nc.variables['T2']
sfcQ = nc.variables['Q2']
sfcP = nc.variables['PSFC']
sfcSwdown = nc.variables['SWDOWN']
sfcRain = np.add(nc.variables['RAINNC'],nc.variables['RAINC'])
sfcU = nc.variables['U10']
sfcV = nc.variables['V10']



times = nc.variables['Times']

# Skip is the length between outputs
skip = 3

x_stations = 	{'nrmn' : 144,
		 'OUN'  : 144,
		 'GUY'  : 114,
		 'GAG'  : 129,
		 'LTS'  : 130,
		 'PNC'  : 145,
		 'ARD'  : 147,
		 'TUL'  : 155,
		 'MCL'  : 156}

y_stations = 	{'nrmn' : 76,
		 'OUN'  : 76,
		 'GUY'  : 89,
		 'GAG'  : 86,
		 'LTS'  : 70,
		 'PNC'  : 89,
		 'ARD'  : 66,
		 'TUL'  : 82,
		 'MCL'  : 73}

meso_ids = 	{'nrmn' : 'nrmn',
		 'OUN'  : 'nrmn',
		 'GUY'  : 'good',
		 'GAG'  : 'wood',
		 'LTS'  : 'altu',
		 'PNC'  : 'blac',
		 'ARD'  : 'ard2',
		 'TUL'  : 'bixb',
		 'MCL'  : 'mcal'}



#press_levels = [1000.,987.5,975.,962.5,950.,937.5,925., 
#                 900.,875.,850.,825.,800.,750.,700.,650.,  
#                 600.,550.,500.,450.,400.,350.,300.,250., 
#                 225.,200.,175.,150.,137.5,125.,112.5,100., 
#                 87.5,75.,62.5,50.,37.5,25.,12.5]

press_levels = [1000.,987.5,975.,962.5,950.,937.5,925., 
                 900.,875.,850.,825.,800.,750.,700.,650.,  
                 600.,550.,500.,450.,400.,350.,300.,250., 
                 225.,200.,175.,150.,137.5,125.,112.5,100., 
                 87.5,75.,62.5]
wind_levs =  [1000.,975.,950.,900.,850.,800.,750.,700.,650.,  
                 600.,550.,500.,450.,400.,350.,300.,250., 
                 225.,200.,175.,150.,137.5,125.,112.5,100., 
                 87.5,75.,62.5]


point_y = 76
point_x = 144
siteid = 'nrmn'
#siteids = ['Guymon','Gage','Altus','Ponca City','Norman','Ardmore','Tulsa']
#sitecalls = ['GUY','GAG','LTS','PNC','OUN','ARD','TUL']
#y_points = [89,85,71,89,76,66,85]
#x_points = [114,127,130,146,144,146,155]
siteids = ['Norman']
y_points = [76]
x_points = [144]
sitecalls = ['OUN']
skewness = 75


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





def timestring(wrftime):
	year  = str(wrftime[2])  + str(wrftime[3])
	month = str(wrftime[5])  + str(wrftime[6])
	day   = str(wrftime[8])  + str(wrftime[9])
	hour  = str(wrftime[11]) + str(wrftime[12])
	outtime = year + month + day + '/' + hour + '00Z'
	return outtime

def hourstring(wrftimes):
	hours = []
	dates = []
	dateticks = []
	ts = wrftimes[0]
	sdstring = str(ts[5]) + str(ts[6]) + '/' + str(ts[8]) + str(ts[9])
	dates.append(sdstring)
	dateticks.append(0)
	for t in range(len(wrftimes)):
		curhour = str(wrftimes[t][11]) + str(wrftimes[t][12]) + 'Z'
		hours.append(curhour)
		if float(curhour[0:2]) == 0.:
			ts = wrftimes[t]
			sdstring = str(ts[5]) + str(ts[6]) + '/' + str(ts[8]) + str(ts[9])
			dates.append(sdstring)
			dateticks.append(t)
	hourticks = range(0,len(wrftimes) * 3 , 3)
	return hours, hourticks, dates, dateticks
			

def meso_comp(ID):
	import time
	import datetime
	import matplotlib.numerix.ma as M
	from matplotlib.dates import YearLocator, MonthLocator, DayLocator, DateFormatter, HourLocator, date2num
	import datetime
	####################################
	file="http://www.mesonet.org/data/public/mesonet/mts/YYYY/MM/DD/YYYYMMDDSTID.mts"
	#######################################################
	# process parameters passed to this script from the command line:
	# Default date is today, length is 6 hours, wrf file is wrfout_d01
	# and mesonetfile is latest_mts for Norman 
	wrf_in = "wrfout_d01"
	meso_file = "latest_mts"
	BT = 0
	degrees = 'fahrenheit'
	fulllength = 72
	skip = 3

	length_flag = 0
	# The length option at the moment will be overridden by the auto-
	# calculation based on the current time

	NS = y_stations[ID]
	WE = x_stations[ID]

	# If there are no matches, we default to Norman (nrmn, NS=47, WE=99)


	# Open the netCDF file and get the starting time

	filetime = times

	str_startyear = str(filetime[0,0])+str(filetime[0,1])+str(filetime[0,2])+str(filetime[0,3])
	str_startmonth = str(filetime[0,5])+str(filetime[0,6])
	str_startday = str(filetime[0,8])+str(filetime[0,9])
	str_starttime = str(filetime[0,11])+str(filetime[0,12])

	str_endyear = str(filetime[-1,0])+str(filetime[-1,1])+str(filetime[-1,2])+str(filetime[-1,3])
	str_endmonth = str(filetime[-1,5])+str(filetime[-1,6])
	str_endday = str(filetime[-1,8])+str(filetime[-1,9])
	str_endtime = str(filetime[-1,11])+str(filetime[-1,12])

	# Produce a string of YYYYMMDD for the start time
	date = str_startyear + str_startmonth + str_startday


	# Make integers of the time values for computational purposes
	startyear = int(str_startyear)
	startmonth = int(str_startmonth)
	startday = int(str_startday)
	starttime = int(str_starttime)

	# Make a time tuple of the start time, also for computational
	# purposes
	startdate = []
	startdate.append(startyear)
	startdate.append(startmonth)
	startdate.append(startday)
	startdate.append(starttime)
	startdate.append(0)
	startdate.append(0)
	startdate.append(0)
	startdate.append(0)
	startdate.append(0)
	#print starttime


	# This is in case we need it for future work, makes a
	# YYYYMMDD string for the end date-time
	enddate = str_endyear + str_endmonth + str_endday + str_endtime
	#print enddate


	# Now, check to see if we have manually specified a length.
	# Otherwise, calculate how many hours it has been since the
	# model's start
	if length_flag == 0:
		from time import localtime
		from time import mktime
		endtime = []
		# This generates a time tuple of the local time
		nowdate = localtime()
		for k in range(len(nowdate)):
			if k<=3:
				endtime.append(nowdate[k])
			else:
				# We don't care about minutes and seconds, so
				# we just make them all 0
				endtime.append(0)
		# See if we are within ten minutes of the new hour and use
		# the previous hour's data if this is so.
		if nowdate[4] < 10:
			endtime[3] = endtime[3] - 1
	
		#print startdate
		#print endtime

		# Get the unix times in seconds for both
		unix_start = mktime(startdate)
		unix_end_lcl = mktime(endtime)
		# Convert the end time in unix seconds to a time
		# tuple in GMT, then back to unix seconds in
		# GMT
		unix_end_gmt = time.gmtime(unix_end_lcl)
		unix_end = mktime(unix_end_gmt)
		#print unix_end
		# Our length is simply the change in seconds
		# divided by 3600 seconds per hour
		length_sec = unix_end - unix_start
		length = int(length_sec / 3600) - 1
		# If the length is longer than the end of the model,
		# make the length the length of the model run
		if length > fulllength:
			length = fulllength
	
		print "Length in seconds is: ", length_sec
		print "Length in hours is: ", length




	# Figure out if more than one day has been selected, and
	# if so, gather the appropriate number of mesonet files
	if ((24-starttime) < length):
		numdays = math.ceil((length-(24-starttime))/24) + 1
		#print "Numdays: ", numdays
	else:
		numdays = 0


	datelist=[]
	### Now we take the date and get the correct Mesonet file
	site=meso_ids[ID]
	yyyymmdd=date
	datelist.append(date)
	numdate = int(date)
	#print "Numdate: ",numdate
	for n in range(int(numdays)):
		numdate = numdate + 1
		#print "Numdate: ",numdate
		datelist.append(str(numdate))
	#print "Datelist len: ", len(datelist)
	#print datelist
	tairs=[]
	wspds=[]
	relhs=[]
	press=[]
	srads=[]
	rains=[]
	wdirs=[]
	dates=[]

	#####################################################
	# construct mesonet file name, and retrieve file
	for j in range(len(datelist)):
		try:
			file="http://www.mesonet.org/data/public/mesonet/mts/YYYY/MM/DD/YYYYMMDDSTID.mts"
			yyyy=datelist[j][0:4]
			mm=datelist[j][4:6]
			dd=datelist[j][6:8]
			file=file.replace('YYYY',yyyy)
			file=file.replace('MM',mm)
			file=file.replace('DD',dd)
			file=file.replace('STID',site.lower())
			#print datelist[j]
			print 'getting file: %s' % file
			content=urllib2.urlopen(file).readlines() # this reads the text at the URL "file"

		except:
			try:   # Maybe it's a different month
				
				file="http://www.mesonet.org/data/public/mesonet/mts/YYYY/MM/DD/YYYYMMDDSTID.mts"
				yyyy=datelist[j][0:4]
				mm=datelist[j][4:6]
				dd=datelist[j][6:8]
				# Add one to the month, make the date the first 
				# and check for new year
				num_month = int(mm)
				num_year = int(yyyy)
	
				num_month = num_month + 1
				if num_month == 13:
					num_year = num_year + 1
					num_month = 1
				num_day = int(dd)
				num_day = 1
				# format as strings
				yyyy = str(num_year)			
				mm = "%(#)02d" % num_month
				dd = "%(#)02d" % num_day
	 
				file=file.replace('YYYY',yyyy)
				file=file.replace('MM',mm)
				file=file.replace('DD',dd)
				file=file.replace('STID',site.lower())
				#print datelist[j]
				print 'getting file: %s' % file
				content=urllib2.urlopen(file).readlines() # this reads the text at the URL "file"

			except:
				print "<p>shucks, we cannot open the above mesonet data file"
				sys.exit()
		##############################################
		# process the data in the file
		nl=0
		nummissing=0
		for line in content:
			nl+=1
			if nl==1 : continue #copyright line, ignore
			if nl==2 : # time header
				try:
					line=line.strip()
					year,month,day=line.split()[1:4]
					yr=int(year)
					mn=int(month)
					dy=int(day)
				except:
					print "choked on parsing time time info"
					sys.exit()
				continue
			if nl==3 : continue # column labels, ignore  
			line=line.strip()
			items=line.split()
			try:		
				stid,stnm,minutes,relh,tair,wspd,wvec,wdir,wdsd,wssd,wmax,rain,pres,srad=items[0:14]
			except:
				print "It crashed! The line did not split into items correctly!"
				sys.exit()
			
			thedate = datetime.date( yr, mn, dy )
			dnum=date2num(thedate)+float(minutes)/1440. # add partial day to matplotlib's date
			dates.append(dnum)
			tairfloat=float(tair)
			if tairfloat < -990.:
				nummissing+=1
			elif  degrees=='fahrenheit':
				tairfloat=tairfloat*1.8+32.
			tairs.append(tairfloat)
			wspdfloat=float(wspd)
			if wspdfloat < -990.:
				nummissing+=1
			wspds.append(wspdfloat)
			relhfloat=float(relh)
			if relhfloat < -990.:
				nummissing+=1
			relhs.append(relhfloat)
			presfloat = float(pres)
			press.append(presfloat)
			rainfloat = float(rain)
			rains.append(rainfloat)
			sradfloat = float(srad)
			srads.append(sradfloat)
			wdirfloat = float(wdir)
			wdirs.append(wdirfloat)
		#print "number of missing values =",nummissing


	#print "Tairs length: ",len(tairs)
	#print "Wspds length: ",len(wspds)
	#print "Relhs length: ",len(relhs)

	####################################################
	# Since we now have a tairs list, lets further parse
	# it down to just the data every hour for the first
	# 6 hours
	###################################################
	obsperhour = 12
	if starttime > 0:
		startob = obsperhour * starttime - 1
	else:
		startob = 0
	totalobs = (obsperhour * length) + startob + 1
	
	hourly_t = tairs[startob:totalobs:obsperhour]
	hourly_spd = wspds[startob:totalobs:obsperhour]
	hourly_relh = relhs[startob:totalobs:obsperhour]
	hourly_pres = press[startob:totalobs:obsperhour]
	hourly_srad = srads[startob:totalobs:obsperhour]
	hourly_rain = rains[startob:totalobs:obsperhour]
	hourly_wdir = wdirs[startob:totalobs:obsperhour]

	# For 3 hour output
	#hourly_t = hourly_t[::3]
	#hourly_spd = hourly_spd[::3]
	#hourly_relh = hourly_relh[::3]
	#hourly_pres = hourly_pres[::3]
	#length = int(length/3) + 1
	print "Length hourly_t: ", len(hourly_t)

	return hourly_t, hourly_relh, hourly_spd, hourly_pres, hourly_rain, hourly_srad, hourly_wdir



def draw_meteo(ID):

	x = x_stations[ID]
	y = y_stations[ID]

	# Compile time list
	# First, handle temperature and dewpoint
	ftime = []
	temp = []
	q2 = []
	uwind = []
	vwind = []
	psfc = []
	rain = []
	swdown = []
	for t in range(len(times)):
		temp.append(sfcT[t,y,x] - 273.)
		q2.append(sfcQ[t,y,x])
		uwind.append(sfcU[t,y,x])
		vwind.append(sfcV[t,y,x])
		psfc.append(sfcP[t,y,x] / 100.)
		rain.append(sfcRain[t,y,x] * 0.03937)
		swdown.append(sfcSwdown[t,y,x])
		ftime.append(t * skip)

	# Get rid of extremely low precip values
	for r in range(len(rain)):
		if rain[r] < 0.01:
			rain[r] = 0.0


	# Dewpoint Calculation
	w = np.divide(q2, np.subtract(1,q2))
	e = np.divide(np.multiply(w,psfc), np.add(.622, w))
	numer = np.multiply(243.5,np.log(np.divide(e, 6.112)))
	denom = np.subtract(17.67, np.log(np.divide(e,6.112)))
	dewp = np.divide(numer,denom)
	del numer, denom, e

	# Wind speed calc
	wind_v = np.power(np.add(np.power(uwind,2),np.power(vwind,2)),0.5)	
	wind_v = np.multiply(wind_v, 2.2369)

	# Calculate a wind direction from u and v
	dirs = []
	for v in range(len(uwind)):
		rawdir = math.atan2(-1 * vwind[v], -1 * uwind[v]) * 180. / 3.14159
		rawdir = 90. - rawdir
		if rawdir < 0:
			rawdir = 360 + rawdir
		dirs.append(rawdir)


	# Convert to Fahrenheit
	temp_F = np.add(np.multiply(temp, (9./5.)), 32)
	dewp_F = np.add(np.multiply(dewp, (9./5.)), 32)
		

	# Grab mesonet data
	mesodata = meso_comp(ID)
	meso_temp = mesodata[0]
	meso_relh = mesodata[1]
	meso_wspd = np.multiply(mesodata[2], 2.2369)
	meso_press = mesodata[3]
	meso_rain = mesodata[4]
	meso_srad = mesodata[5]
	meso_wdir = mesodata[6]

	# Convert RELH to dewpoint
	meso_temp_C = np.multiply(np.subtract(meso_temp, 32), 5./9.)
	es_num = np.multiply(17.67,meso_temp_C)
	es_den = np.add(meso_temp_C, 243.5)
	es = np.multiply(6.112, np.exp(np.divide(es_num, es_den)))
	e = np.divide(np.multiply(meso_relh,es), 100.)
	numer = np.multiply(243.5,np.log(np.divide(e, 6.112)))
	denom = np.subtract(17.67, np.log(np.divide(e,6.112)))
	meso_dewp = np.divide(numer,denom)
	meso_dewp_F = np.add(np.multiply(meso_dewp, 9./5.), 32)
	
	meso_rain = np.multiply(meso_rain, 0.03937)

	meso_len = len(meso_temp)
	obtime = []
	for n in range(meso_len):
		obtime.append(n)

	# Begin plotting
	pylab.figure(figsize=(10,8), frameon = False)

	# Get our dates and times
	temporalticks = hourstring(times)
	hours = temporalticks[0]
	hourticks = temporalticks[1]
	dates = temporalticks[2]
	dateticks = temporalticks[3]
	#print hours, hourticks

	# Now make temperature plot
	pylab.subplot(511)
	#print temp_F
	#print meso_temp
	tempplot = pylab.plot(ftime,temp_F,'r-')
	dewpplot = pylab.plot(ftime,dewp_F,'g-')
	mesodewpplot1 = pylab.plot(obtime,meso_dewp_F,'go-')
	ymin = pylab.axis()[2]
	tempfplot = pylab.fill_between(ftime,ymin,temp_F, facecolor = 'pink')
	dewpfplot = pylab.fill_between(ftime,ymin,dewp_F, facecolor = 'palegreen')
	mesotempplot = pylab.plot(obtime,meso_temp,'ro-')
	mesodewpplot = pylab.plot(obtime,meso_dewp_F,'go-')
	pylab.grid(True)
	pylab.title('Forecast Meteogram for %s   WRF init %s' % (ID, timestring(times[0])))
	ax = pylab.gca()
	ax.set_xlim([0,ftime[-1]])
	pylab.xticks(hourticks, hours)
	for label in ax.get_yticklabels():
		label.set_fontsize(8)
	for label in ax.get_xticklabels():
		label.set_fontsize(8)
	pylab.text(-0.07,0.5,'Air Temperature [F]', rotation='vertical',fontsize = 8, color = 'r', verticalalignment='center', transform = ax.transAxes)
	pylab.ylabel('Dewpoint Temp [F]', fontsize = 8, color = 'g')

	# Now for the wind speed / direction plot
	ax1 = pylab.subplot(512)
	windvplot = ax1.plot(ftime,wind_v, 'b-')
	mesowindvplot1 = ax1.plot(obtime,meso_wspd, 'bo-')
	ymin = ax1.axis()[2]
	windvfplot = ax1.fill_between(ftime,ymin,wind_v, facecolor = 'lightskyblue')
	mesowindvplot = ax1.plot(obtime,meso_wspd, 'bo-')
	ax1.grid(True)
	for label in ax1.get_yticklabels():
		label.set_fontsize(8)
	ax1.set_ylabel('Wind Speed [mph]', fontsize = 8, color = 'b')

	ax2 = ax1.twinx()
	winddplot = ax2.plot(ftime, dirs, marker = '^', mec = 'darkslateblue', mfc = 'darkseagreen', linestyle = 'none')
	mesowinddplot = ax2.plot(obtime, meso_wdir, marker = 'o', mec='darkslategray', mfc = 'darkgray', linestyle = 'none')
	ax2.set_ylim([0,360])
	ax2.set_xlim([0,ftime[-1]])
	pylab.xticks(hourticks, hours)
	for label in ax2.get_yticklabels():
		label.set_fontsize(8)
	for label in ax2.get_xticklabels():
		label.set_fontsize(8)
	for label in ax1.get_xticklabels():
		label.set_fontsize(8)
	pylab.yticks([0.,45.,90.,135.,180.,225.,270.,315.,360], ['N','NE','E','SE','S','SW','W','NW','N'])
	ax2.set_ylabel('Wind Direct. [deg]', fontsize = 8, color = 'brown')

	# Pressure Plot
	pylab.subplot(513)
	print psfc
	print meso_press
	press_plot = pylab.plot(ftime, psfc, color='brown', linestyle='solid')
	mesopressplot1 = pylab.plot(obtime, meso_press, color = 'brown', linestyle = '-', marker = 'o')
	ymin = pylab.axis()[2]
	pressfplot = pylab.fill_between(ftime, ymin, psfc, facecolor = 'rosybrown')
	mesopressplot = pylab.plot(obtime, meso_press, color = 'brown', linestyle = '-', marker = 'o')
	pylab.grid(True)
	ax = pylab.gca()
	ax.set_xlim([0,ftime[-1]])
	pylab.xticks(hourticks, hours)
	for label in ax.get_yticklabels():
		label.set_fontsize(8)
	for label in ax.get_xticklabels():
		label.set_fontsize(8)
	pylab.ylabel('Pressure [mb]', fontsize = 8, color = 'brown')

	# Rainfall Plot
	pylab.subplot(514)
	rain_plot = pylab.plot(ftime, rain, 'g-')
	pylab.grid(True)
	ax = pylab.gca()
	ax.set_ylim([0,1.])
	ymin = pylab.axis()[2]
	rainfplot = pylab.fill_between(ftime,ymin,rain, facecolor = 'mediumaquamarine')
	mesorainplot = pylab.plot(obtime, meso_rain, 'go-')
	ax.set_xlim([0,ftime[-1]])
	pylab.xticks(hourticks, hours)
	for label in ax.get_yticklabels():
		label.set_fontsize(8)
	for label in ax.get_xticklabels():
		label.set_fontsize(8)
	pylab.ylabel('Accum. Rainfal [in]', fontsize = 8, color = 'g')

	# SRad Plot
	pylab.subplot(515)
	srad_plot = pylab.plot(ftime, swdown, 'y-')
	pylab.grid(True)
	ax = pylab.gca()
	ax.set_ylim([0,1000])
	ymin = pylab.axis()[2]
	sradfplot = pylab.fill_between(ftime,ymin,swdown, facecolor = 'khaki')
	mesorainplot = pylab.plot(obtime, meso_srad,color='goldenrod', marker = 'o', linestyle='solid')
	ax.set_xlim([0,ftime[-1]])
	pylab.xticks(hourticks, hours)
	for label in ax.get_yticklabels():
		label.set_fontsize(8)
	for label in ax.get_xticklabels():
		label.set_fontsize(8)
	for d in range(len(dates)):
		pylab.text(dateticks[d] * 3, -325., dates[d],rotation='horizontal',fontsize = 10, color = 'k', horizontalalignment = 'center')
	pylab.ylabel('Solar Radiation [W m-2]', fontsize = 8, color = 'y')
	

	dom = 'wrf'
	filename = '%s_%s_meteo' % (dom, ID)
	pylab.savefig(filename)
       	pylab.close()
	os.system('convert -render -flatten %s.png %s.gif' % (filename,filename))
	os.system('rm %s.png' % filename)
			


print "Plotting"
draw_meteo(cursite)


os.system('scp *meteo.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')
print "Done."

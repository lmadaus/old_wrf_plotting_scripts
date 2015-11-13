#!/usr/bin/python

# Working script to generate maps from wrfout netCDF files
# using matplot lib with basemap
# Written by David John Gagne II
# Amended by Luke Madaus

import Nio
import matplotlib
matplotlib.use('agg')
import pylab
import math
import numpy as np
import datetime
import os
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import calc_wrf_severe as severe


filename = '../wrfout_d01_PLEV.nc'
filename_sfc = '../wrfout_d01.nc'
nc = Nio.open_file(filename)
ncsfc = Nio.open_file(filename_sfc)
#print nc.variables.keys()
thetas_in = nc.variables['T']
qhum = nc.variables['QVAPOR']
uwind = nc.variables['UU']
vwind = nc.variables['VV']
PH = nc.variables['PH']
PHB = nc.variables['PHB']
zheight = np.divide(np.add(PH,PHB),9.81)

elev = ncsfc.variables['HGT']
sfcT = ncsfc.variables['T2']
sfcQ = ncsfc.variables['Q2']
sfcP = ncsfc.variables['PSFC']

times = ncsfc.variables['Times']

# Skip is the length between outputs
skip = 3

restart_time = 0

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
siteid = 'Norman'
#siteids = ['Guymon','Gage','Altus','Ponca City','Norman','Ardmore','Tulsa']
#sitecalls = ['GUY','GAG','LTS','PNC','OUN','ARD','TUL']
#y_points = [89,85,71,89,76,66,85]
#x_points = [114,127,130,146,144,146,155]
siteids = ['Norman']
y_points = [76]
x_points = [144]
sitecalls = ['OUN']
skewness = 75


dayofweek = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']
# Constants

grav = 9.81

# Epsilons for moisture
ezero = 6.112

eps = .622

cp = 1004.
rgas = 287.04

gamma = (rgas/cp)

# cp_moist = cp*(1 + cpmd*qvp)
cpmd = .887

# rgas_moist = rgas * (1+rgasmd*qvp)
rgasmd = .608

g = 9.81


def draw_isotherms():
	itherms = [[]]
	itherms.append([])
	for n in [-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40]:
		itherms[0].append(n)
		itherms[1].append(n + skewness)
	endvals = [1000.,100.]
	for n in range(len(itherms[1])):
		endpoints = []
		endpoints.append(itherms[0][n])
		endpoints.append(itherms[1][n])
		if itherms[0][n] < 0:
			pylab.plot(endpoints,endvals,color = 'blue', \
				linestyle='dashed', linewidth = .5)
		elif itherms[0][n] > 0:
			pylab.plot(endpoints,endvals,color = 'red', \
				linestyle='dashed', linewidth = .5)
		else:
			pylab.plot(endpoints,endvals,color = 'blue', \
				linestyle='solid', linewidth = .5)


def draw_dry_adiabats():
	dadiabats = []
	for n in [-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150]:
		theta = n + 273
		curadiabat = []
		for m in range(len(press_levels)):		
			curadiabat.append((theta * np.power(press_levels[m]/1000.,(287.04/1004.)) - 273) + (skewness * math.log(1000./press_levels[m],10)))
		dadiabats.append(curadiabat)
	endvals = press_levels
	for n in range(len(dadiabats)):
		pylab.plot(dadiabats[n],endvals,color = 'brown', \
			linestyle='dashed', linewidth = .5)



def draw_isobars():
	endxs = [-40,50]
	for n in range(100,1000,100):
		endys = [n,n]
		pylab.plot(endxs, endys, color = 'k',\
			linewidth = .5)
	

#def draw_moist_adiabats(skewness):
#	infile = open('./capeline.txt','r')
#	Tp = []
#	Tp_trans = []
#	press_p = []
#	cape = 0
#	for line in infile:
#		linelist = line.split()
#		if len(linelist) > 1:
#			Tp_curr = float(linelist[1])
#			press_curr = float(linelist[0])
#			Tp.append(Tp_curr)
#			press_p.append(press_curr)
#			Tp_trans.append(Tp_curr + skewness * math.log(1000./press_curr,10))
#		else:
#			cape = float(linelist[0])		
#	infile.close()
#	return Tp_trans, press_p, cape, Tp
def calc_es(T_K):
	es = 6.112*np.exp(17.67*(T_K-273)/((T_K-273)+243.5))
	return es 

def draw_moist_adiabats(skewness,parcel,parcel_T,parcel_P):
	def calc_gamma_s(T,p):
		p = p/1000.
		es = .611 * math.exp(5423*((1/273.)-(1./(T))))
		r = eps * es / (p - es)
		gamma_d = .0098
		a = 0.28571
		b = 1.35e7
		c = 2488.4
		numer = (a * T) + (c * r)
		denom = p * (1 + (b * r / (T*T))) 
		gamma_s = numer/denom
		# print "P: %f T: %f es: %f Gs: %f" % (p, T, es, gamma_s)
		return gamma_s	

	def RK_4(p,T,dp):
		k1 = calc_gamma_s(T,p)
		k2 = calc_gamma_s((T + .5 * dp * k1),(p + .5 * dp))
		k3 = calc_gamma_s((T + .5 * dp * k2),(p + .5 * dp))
		k4 = calc_gamma_s((T + dp * k3),(p + dp))
	
		Tp = T + dp * (1/6.) * (k1 + 2. * k2 + 2. * k3 + k4)
		return Tp

	def RK_45(p,T,dp):
		k1 = calc_gamma_s(T,p)
		k2 = calc_gamma_s((T + .25 * dp * k1),(p + .25 * dp))
		k3 = calc_gamma_s((T + 3/32. * dp * k1 + 9/32. * dp * k2),(p + 3/8. * dp))
		k4 = calc_gamma_s((T + dp * 1932/2197. * k1 - 7200/2197. * dp * k2 + 7296/2197. * dp * k3),(p + 12/13. * dp))
		k5 = calc_gamma_s((T + dp * 439/216. * k1 - 8. * dp * k2 + 3680/513. * dp * k3 - 845/4104 * dp * k4),(p + dp))
		k6 = calc_gamma_s((T - dp * 8/27. * k1 + 2. * dp * k2 - 3544/2565. * dp * k3 + 1859/4104 * dp * k4 - 11/40. * dp * k5),(p + .5 * dp))
	
		Tp = T + dp * (16/135. * k1 + 6656/12825. * k3 + 28561/56430. * k4 - 9/50. * k5 + 2/55. * k6)
		return Tp

	if parcel == 0:
		for base_T in [-10,0,10,14,18,22,26,31,35]:
			madiabat = []
			madiabat_trans = []
			Tp = base_T + 273.
			
			for z in range(len(press_levels)):	
				p = press_levels[z]
				if z == 0:
					dp = 0
				else:
					dp = p - press_levels[z-1]		
				dt = calc_gamma_s(Tp,p) * (dp/1000.)
				#print "dt: ", dt
				#new_Tp = Tp + dt
				new_Tp = RK_4(p*100.,Tp,(dp/10.)) 
				Tp = new_Tp
				madiabat.append(Tp)
				Tp_C = Tp - 273.
				Tp_trans = Tp_C + skewness * math.log(1000./p,10)

				#print p, Tp-273., Tp_trans
				madiabat_trans.append(Tp_trans)	
			pylab.semilogy(madiabat_trans, press_levels, color = 'green', basey = 10, linestyle = 'dotted', linewidth = .5)

	if parcel == 1:
		madiabat = []
		madiabat_trans = []
		Tp = parcel_T + 273.
		madiabat_press_levels = []
		for z in range(len(press_levels[parcel_P:])):	
			p = press_levels[z+parcel_P]
			if z == 0:
				dp = 0
			else:
				dp = p - press_levels[z-1 + parcel_P]		
			dt = calc_gamma_s(Tp,p) * (dp/1000.)
			#print "dt: ", dt
			#new_Tp = Tp + dt
			new_Tp = RK_45(p*100.,Tp,(dp/10.)) 
			Tp = new_Tp
			madiabat.append(Tp)
			Tp_C = Tp - 273.
			Tp_trans = Tp_C + skewness * math.log(1000./p,10)

			#print p, Tp-273., Tp_trans
			madiabat_trans.append(Tp_trans)	
			madiabat_press_levels.append(p)
		return np.subtract(madiabat, 273)
		#print "Drawing parcel moist adiabat..."
		#print madiabat
		#print madiabat_press_levels
		#raw_input()
		#pylab.semilogy(madiabat_trans, madiabat_press_levels, color = 'brown', basey = 10, linestyle = 'dashed', linewidth = 1)



def draw_parcel_trace(Tparc, Press):

	# Convert Pressures to log scale
	Pfact = np.multiply(skewness,np.log10(np.divide(1000., Press)))

	minlen = min(len(Tparc), len(Pfact))
	
	dry_parcel_trace_trans = np.add(Tparc[:minlen],Pfact[2:minlen+2])	


	pylab.semilogy(dry_parcel_trace_trans,Press[2:minlen+2],\
		basey=10, color = 'brown', linestyle = 'dashed',\
		linewidth = 1.5)

def new_draw_parcel_trace(Tb, PLCL, Press):

	# Convert Pressures to log scale
	Pfact = np.multiply(skewness,np.log10(np.divide(1000., Press)))

	parcelT = []
	flag = 1

	for p in range(len(Press)):
		if Press[p] >= PLCL:
			newTB = ((Tb + 273.) * (Press[p]/Press[0]) ** (287.04/1004.)) - 273.
			parcelT.append(newTB)
		else:
			if flag:
				if p == 0:
					moists = draw_moist_adiabats(0, 1, Tb, 0)
				else:
					moists = draw_moist_adiabats(0,1,parcelT[p-1], (p - 1 + len(press_levels) - len(Press)))
				for m in moists:
					parcelT.append(m)
				flag = 0


	minlen = min(len(parcelT), len(Pfact))
	
	dry_parcel_trace = np.add(parcelT[:minlen], Pfact[:minlen])



	pylab.semilogy(dry_parcel_trace,Press[:minlen],\
		basey=10, color = 'brown', linestyle = 'dotted',\
		linewidth = 1.5)
	


def get_severes():
	# All we need is pressures, temperatures and humidities
	#W = np.divide(qvap_list, np.subtract(1, qvap_list))
	W = qvap_list
	PR_h = press_levels
	T_Ks = np.add(temp_list, 273.)
	minlen = min(len(PR_h), len(W), len(T_Ks))
	PR_h = PR_h[:minlen]
	T_Ks = T_Ks[:minlen]
	W = W[:minlen]
	severes = severe.CAPESOUND(PR_h, T_Ks, W)
	scape = severes[1]
	mucape = severes[0]
	cinh = severes[4]
	parcelT = severes[5]
	PLCL = severes[2]
	if scape < 0:
		scape = 0.0
		
	return scape,mucape,cinh,parcelT,PLCL


def timestring(wrftime,curtime):
	curtime_str = '%02d' % curtime
	year  = str(wrftime[2])  + str(wrftime[3])
	month = str(wrftime[5])  + str(wrftime[6])
	day   = str(wrftime[8])  + str(wrftime[9])
	hour  = str(wrftime[11]) + str(wrftime[12])
	numdow = datetime.date.weekday(datetime.date(int(year), int(month), int(day)))
	dow = dayofweek[numdow]
	outtime = dow + ' ' + year + month + day + '/' + hour + '00Z F'+ curtime_str
	return outtime



for k in range(len(siteids)):
	siteid = siteids[k]
	point_y = y_points[k]
	point_x = x_points[k]
	sitecall = sitecalls[k]

	# Begin looping through times
	for time in range(len(thetas_in)):
	#for time in [1]:
		print 'At ',time*skip
		temp_list = []
		dewp_list = []
		qvap_list = []
		temp_trans = []
		dewp_trans = []
		actual_press_levels = []
		u_wind_kts = []
		v_wind_kts = []
		flag = 0
		actual_wind_levs = []
		sfctemp_C = sfcT[time,point_y,point_x] - 273.
		temp_list.append(sfctemp_C)
		sfc_press_hpa = sfcP[time,point_y,point_x] * 0.01
		actual_press_levels.append(sfc_press_hpa)
		sfcqvap = sfcQ[time,point_y,point_x]
		w = sfcqvap / (1 - sfcqvap)
		e = (w * sfc_press_hpa/(0.622+w))
		#es = 6.112*np.exp(17.67*(T_K-273)/((T_K-273)+243.5))
		#e = RH *.01 * es
		sfcTd_C = (243.5*np.log(e/6.112))/(17.67-np.log(e/6.112))
		dewp_list.append(sfcTd_C)

		temp_trans.append(sfctemp_C + skewness * math.log(1000./sfc_press_hpa,10))
		dewp_trans.append(sfcTd_C + skewness * math.log(1000./sfc_press_hpa,10))
		
		for height in range(len(press_levels)):
			try:
				TH_K = thetas_in[time,height,point_y,point_x] + 300.0
				T_K = TH_K * math.pow((press_levels[height]/1000.),\
					(287.04/1004))
				temp_list.append(T_K - 273.)
				
				qvap = qhum[time,height,point_y,point_x]
				qvap_list.append(qvap)
				w = qvap / (1 - qvap)
				e = (w * press_levels[height]/(0.622+w))
				es = 6.112*np.exp(17.67*(T_K-273)/((T_K-273)+243.5))
				#e = RH *.01 * es
				Td_C = (243.5*np.log(e/6.112))/(17.67-np.log(e/6.112))
				dewp_list.append(Td_C)
				goodwind = 0
				curr_uwind = uwind[time,height,point_y,point_x] * 1.9438
				curr_vwind = vwind[time,height,point_y,point_x] * 1.9438
				flag = 1
				print press_levels[height], T_K-273, Td_C
			
			except:
				flag = 0
			if flag:
				actual_press_levels.append(press_levels[height])
				dewp_trans.append(Td_C + skewness * math.log(1000./press_levels[height],10))
				temp_trans.append((T_K - 273) + skewness * math.log(1000./press_levels[height],10))
				for g in wind_levs:
					if press_levels[height] == g:
						u_wind_kts.append(curr_uwind)
						v_wind_kts.append(curr_vwind)
						actual_wind_levs.append(press_levels[height])





		# Set Figure Size (1000 x 800)
		pylab.figure(figsize=(10,8), frameon=False)
		draw_isotherms()
		draw_dry_adiabats()
		draw_isobars()
		#Tp_trans = draw_moist_adiabats(skewness)[0]
		#press_p = draw_moist_adiabats(skewness)[1]
		#cape = draw_moist_adiabats(skewness)[2]
		#Tp = draw_moist_adiabats(skewness)[3]
		#print len(dewp_trans), len(actual_press_levels)
	

		# Get CAPE
		severes = get_severes()
		sbcape = severes[0]
		mucape = severes[1]
		cinh = severes[2]
		Tparc = severes[3]
		PLCL = severes[4]

		#draw_parcel_trace(Tparc,press_levels)
		#new_draw_parcel_trace((Tparc[0]), PLCL, actual_press_levels)

		#Kitbash parcel trace


		draw_moist_adiabats(skewness,0,0,0)
		#draw_parcel_trace(temp_list[0] + 273,qvap_list[0],actual_press_levels[0])
		pylab.semilogy(dewp_trans,actual_press_levels, \
			color = 'green',basey=10, linewidth = 2)
		pylab.semilogy(temp_trans,actual_press_levels, \
			color = 'red',basey=10, linewidth = 2)
	



		#print Tp, press_p
		#raw_input()
		#pylab.semilogy(Tp_trans,press_p,color = 'brown',\
		#	basey = 10, linewidth = 2, linestyle = 'dashed')
		baraxis = []
		# Need this -3 to the list to stop them from going off the top
		for n in range(len(u_wind_kts[:-5])):
			baraxis.append(45.)
		pylab.plot([45,45],[100,1000],linewidth = .75, color = 'k')
		bb = pylab.barbs(baraxis,actual_wind_levs[:-5],u_wind_kts[:-5],v_wind_kts[:-5],\
			linewidth = .75)
		bb.set_clip_box(None)
		ax = pylab.gca()
		#ax.set_ylim(ax.get_ylim()[::-1])
		ax.set_ylim([1000,100])
		ax.set_xlim([-40,50])
		majorLocator = MultipleLocator(100.)
		majorFormatter = FormatStrFormatter('%4.0f')
		ax.yaxis.set_major_locator(majorLocator)
		ax.yaxis.set_major_formatter(majorFormatter)

		stemps = sfcT[time,point_y,point_x]+6.5*ncsfc.variables['HGT'][time,point_y,point_x]/1000.
		mslp = sfcP[time,point_y,point_x]*np.exp(9.81/(287.0*stemps)*ncsfc.variables['HGT'][time,point_y,point_x])*0.01 + (6.7 * ncsfc.variables['HGT'][time,point_y,point_x] / 1000)

		# Convert Celsius Temps to Fahrenheit
		ftemp = (9./5.)*(sfcT[time,point_y,point_x]-273) + 32
	
		# Find saturation vapor pressure
		es = 6.112 * np.exp(17.67 * sfcT[time,point_y,point_x]/(sfcT[time,point_y,point_x] + 243.5))
		w = sfcQ[time,point_y,point_x]/(1-sfcQ[time,point_y,point_x])
		e = (w * sfcP[time,point_y,point_x] / (.622 + w)) / 100
		Td_C = (243.5 * np.log(e/6.112))/(17.67-np.log(e/6.112))
		Td_F = (Td_C * 9 / 5) + 32


		print 'Plotting time ',time*skip + restart_time
		curtimestring = timestring(times[time],time*skip + restart_time)
	
	
		TITLESTRING = 'Skew-T at %s' % (siteid) 
		dom = 'wrf'
		pylab.title('OWL/KHLWE WRF-ARW %s   Valid: %s' % (TITLESTRING, curtimestring), \
			fontsize=11,bbox=dict(facecolor='white', alpha=0.65),\
			x=0.5,y=.95,weight = 'demibold',style='oblique', \
			stretch='normal', family='sans-serif')
		filename = '%s_%s_fsound_f%02d' % (dom,sitecall,time*skip + restart_time )

		pylab.suptitle('Surface: T(F): %4.1f    Td(F): %3.1f     MSLP(mb): %5.1f' % (ftemp,Td_F,mslp), \
			fontsize=11, x = 0.5, y = 0.07)		
		#pylab.suptitle('CAPES (J/kg):  SB: %4.0f    MU: %4.0f     SCINH: %3.0f     PLCL(mb): %3.0f' % (sbcape,mucape,cinh,PLCL), \
		#	fontsize=11, x = 0.5, y = 0.045)		
		pylab.suptitle('PLCL(mb): %3.0f' % (PLCL), \
			fontsize=11, x = 0.5, y = 0.045)		

		pylab.suptitle(u'\u00B0' + 'C', fontsize=11, x = 0.13, y = 0.07)		
		pylab.suptitle('mb', fontsize=11, x = 0.109, y = 0.89)		


		#filename = 'skewt_%02d.png' % (time)
		#filename = 'skewt.png'
		pylab.savefig(filename)
        	pylab.close()
		os.system('convert -render -flatten %s.png %s.gif' % (filename,filename))
		os.system('rm %s.png' % filename)


os.system('scp *fsound*.gif hoot@10.197.1.220:/usr/home/hoot/http/models_data/wrf/.')
print "Done."

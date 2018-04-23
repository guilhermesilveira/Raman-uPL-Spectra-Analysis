import numpy as np
import csv
import math
import matplotlib.pyplot as plt 
import matplotlib 
import scipy.optimize as spop
import scipy.interpolate as spip
import sys
import os
import datetime
from os.path import exists


h 		= 4.136e-15 #[eV-s]
c 		= 3e10 #[cm/s]
nm2cm 	= 1e-7

def FitGaussian(x, y, xmin=430, xmax=460, DEBUGFILE=False):
	'''
	Fit a single Gaussian peak plus a linear background to dataset (x,y), where x is in [eV] and y is counts
	
	PARAMETERS:
	x: wavelength vector [numpy vector]
	y: intensity vector [numpy vector]
	(Optional) xmin: low energy cut-off for fit [float/eV]
	(Optional) xmax: high energy cut-off for fit [float/eV]
	(Optional) DEBUGFILE: If this is not false, then raw data and fit will be plotted and saved to 
	the file name <DEBUGFILE>
	
	RETURNS:
	A list, where element #<x> is:
	0: The Gaussian function amplitude  [float/arb]
	1: The Gaussian function center [float/eV]
	2: y(x)-gausienne
	'''
	ii 		= np.nonzero(np.logical_and(x >= xmin, x <= xmax))[0]
	xm		= x
	ym		= y
	x 		= x[ii]
	y 		= y[ii]
	a0 		= np.max(y)
	imx 	= np.argmax(y)
	x0 		= x[imx]
	ilo 	= imx - 1 
	ihi 	= imx + 1
	#Estimate the FWHM 
	while ilo > 1 and y[ilo] > 0.5*a0:
		ilo 	= ilo - 1
	while ihi < ii.shape[0] - 2 and y[ihi] > 0.5*a0:
		ihi 	= ihi + 1	
	sig0 	= np.abs(x[ihi]-x[ilo]) / (2*np.sqrt(2*np.log(2)))
	#Estimate the initial background (form <y = m*x + b>)
	m 		= 0 
	b 		= 0 
	#Define the function to minimize via the least-squares method 
	def ObjFxn(xf):
		A 		= np.abs(xf[0])
		x0 		= xf[1]
		sig 	= xf[2]
		m 		= xf[3]
		b 		= xf[4]
		error 	=  y - A*np.exp( -np.power((x-x0)/(math.sqrt(2)*sig), 2)) - m*x - b
		return sum(error**2)
	#Do the fit
	x0 		= np.array([a0,x0,sig0,m,b])
	xbound = [(1000,100000),(440,450),(4,50),(None,None),(0,3000)]
	xfit 	= spop.fmin_l_bfgs_b(ObjFxn, x0, bounds=xbound, approx_grad=True)[0]

	#Define the fit function
	def FitFxn(xf):
		A 		= np.abs(xf[0])
		x0 		= xf[1]
		sig 	= xf[2]
		m 		= xf[3]
		b 		= xf[4]
		return lambda x: A*np.exp( -np.power((x-x0)/(math.sqrt(2)*sig), 2)) + m*x + b		
	xf  	= np.linspace(xmin,xmax,1000)
	fit 	= FitFxn(xfit)
	diff	= ym-fit(xm)
	
	if DEBUGFILE:	
		plt.plot(xm,ym,'ko',ms=4,label='Data')
		plt.plot(xm,fit(xm),'r-',lw=1.5,label='Fit')
		plt.plot(xm,diff,'r-',lw=1.5,label='Diff')
		plt.legend(loc='upper left')
		plt.xlabel('wavelength [nm]')
		plt.ylabel('Count rate [1/s]')
		plt.savefig('DEBUGFILE/%s'% (DEBUGFILE))
		plt.clf()
	
	
	#Return the fit parameters
	return [xfit[0],xfit[1],diff]

def FitLorrentz1(x, y, xmin=445, xmax=447, DEBUGFILE=False):
	'''
	Fit a single Gaussian peak plus a linear background to dataset (x,y), where x is in [eV] and y is counts
	
	PARAMETERS:
	x: wavelength vector [numpy vector]
	y: intensity vector [numpy vector]
	(Optional) xmin: low energy cut-off for fit [float/eV]
	(Optional) xmax: high energy cut-off for fit [float/eV]
	(Optional) DEBUGFILE: If this is not false, then raw data and fit will be plotted and saved to 
	the file name <DEBUGFILE>
	
	RETURNS:
	A list, where element #<x> is:
	0: A
	1: x0
	2: gam = fwhm
	'''
	ii 		= np.nonzero(np.logical_and(x >= xmin, x <= xmax))[0]
	xm		= x
	ym		= y
	x 		= x[ii]
	y 		= y[ii]
	a0 		= np.max(y)
	imx 	= np.argmax(y)
	x0 		= x[imx]
	ilo 	= imx - 0 
	ihi 	= imx + 0
	#Estimate the FWHM 
	while ilo > 1 and y[ilo] > 0.5*a0:
		ilo 	= ilo - 1
	while ihi < ii.shape[0] -1 and y[ihi] > 0.5*a0:
		ihi 	= ihi + 1
	gam0 	= np.abs(x[ihi]-x[ilo])
	#Estimate the initial background (form <y = m*x + b>)
	m 		= 0 
	b 		= 0 
	#Define the function to minimize via the least-squares method 
	def ObjFxn(xf):
		A 		= np.abs(xf[0])
		x0 		= xf[1]
		gam 	= xf[2]
		m 		= xf[3]
		b 		= xf[4]
		return y - A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2)))
		#ne pas mettre A au mumerateur, ca compliaue les calculs
		
	#Do the fit
	x0 		= np.array([a0,x0,gam0,m,b])
	xfit 	= spop.leastsq(ObjFxn, x0)[0]

	#For debugging purposes, plot each fit individually
	if DEBUGFILE:
		#Define the fit function
		def FitFxn(xf):
			A 		= np.abs(xf[0])
			x0 		= xf[1]
			gam 	= xf[2]
			m 		= xf[3]
			b 		= xf[4]
			return lambda x: A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2)))	
		xf  	= np.linspace(xmin,xmax,1000)
		fit 	= FitFxn(xfit)		
		plt.plot(xm,ym,'ko',ms=4,label='Data')
		plt.plot(xm,fit(xm),'r-',lw=1.5,label='Fit')
		plt.legend(loc='upper left')
		plt.xlabel('wavelength [nm]')
		plt.ylabel('Count rate [1/s]')
		plt.savefig('DEBUGFILE/%s'% (DEBUGFILE))
		plt.clf()
		
	#Return the fit parameters
	return [xfit[0],xfit[1],fit(xm)]

def FitLorrentz2(x, y, xmin=453.5, xmax=456.5, DEBUGFILE=False):
	'''
	Fit a single Gaussian peak plus a linear background to dataset (x,y), where x is in [eV] and y is counts
	
	PARAMETERS:
	x: wavelength vector [numpy vector]
	y: intensity vector [numpy vector]
	(Optional) xmin: low energy cut-off for fit [float/eV]
	(Optional) xmax: high energy cut-off for fit [float/eV]
	(Optional) DEBUGFILE: If this is not false, then raw data and fit will be plotted and saved to 
	the file name <DEBUGFILE>
	
	RETURNS:
	A list, where element #<x> is:
	0: A
	1: x0
	2: gam = fwhm
	'''
	ii 		= np.nonzero(np.logical_and(x >= xmin, x <= xmax))[0]
	xm		= x
	ym		= y
	x 		= x[ii]
	y 		= y[ii]
	a0 		= np.max(y)
	imx 	= np.argmax(y)
	x0 		= x[imx]
	ilo 	= imx - 1 
	ihi 	= imx + 1
	#Estimate the FWHM 
	while ilo > 1 and y[ilo] > 0.5*a0:
		ilo 	= ilo - 1
	while ihi < ii.shape[0] -1 and y[ihi] > 0.5*a0:
		ihi 	= ihi + 1
	gam0 	= np.abs(x[ihi]-x[ilo])
	#Estimate the initial background (form <y = m*x + b>)
	m 		= 0 
	b 		= 0 
	#Define the function to minimize via the least-squares method 
	def ObjFxn(xf):
		A 		= np.abs(xf[0])
		x0 		= xf[1]
		gam 	= xf[2]
		m 		= xf[3]
		b 		= xf[4]
		return y - A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2)))
		#ne pas mettre A au mumerateur, ca compliaue les calculs
		
	#Do the fit
	x0 		= np.array([a0,x0,gam0,m,b])
	xfit 	= spop.leastsq(ObjFxn, x0)[0]

	#For debugging purposes, plot each fit individually
	if DEBUGFILE:
		#Define the fit function
		def FitFxn(xf):
			A 		= np.abs(xf[0])
			x0 		= xf[1]
			gam 	= xf[2]
			m 		= xf[3]
			b 		= xf[4]
			return lambda x: A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2)))	
		xf  	= np.linspace(xmin,xmax,1000)
		fit 	= FitFxn(xfit)		
		plt.plot(xm,ym,'ko',ms=4,label='Data')
		plt.plot(xm,fit(xm),'r-',lw=1.5,label='Fit')
		plt.legend(loc='upper left')
		plt.xlabel('wavelength [nm]')
		plt.ylabel('Count rate [1/s]')
		plt.savefig('DEBUGFILE/%s'% (DEBUGFILE))
		plt.clf()
		
	#Return the fit parameters
	return [xfit[0],xfit[1],fit(xm)]	
	
def ParseFile(f):
	with open(f, 'rb') as csvfile:
		rd = csv.reader(csvfile, delimiter='	')
		wl  		= []
		dat 		= []
		for row in rd:
			wl.append(float(row[0]))
			dat.append(np.array(map(lambda x: float(x), row[1:])))
		dat 		= np.array(dat)
		wl 			= np.array(wl)
	d			= os.path.getmtime(f)
	t		= datetime.datetime.fromtimestamp(d)
	return [t,wl, dat]

if __name__=="__main__":
	#Increase the matplotlib font size 
	font = {'family' : 'normal', 'size'   : 16}
	matplotlib.rc('font', **font)
	
	#Read EPFL's data from the file	
	n = 13 #nombre de fichiers dans le dossier
	dat = np.zeros((2048,n))
	wl = np.zeros((2048,n))
	i = 0
	folder_path = "measurements"
	f  = []
	t  = []
	
	for root, dirs, files in os.walk(folder_path):
		for filename in files:
			[tf, wlf, datf] 	= ParseFile('%s'% (os.path.join(root, filename)))
	
			for a in range (datf.shape[0]): #on rempli dat
				dat[a,i]= datf[a]
				wl[a,i]= wlf[a]
			f.append(filename)
			t.append(tf)
			i=i+1
	
	#Plot the EPFL spectra
	plt.plot(wl[:,0], dat,lw=1)
	plt.xlabel('wavelenght [nm]')
	plt.ylabel('Count rate [1/s]')
	plt.grid()
	plt.savefig('mPL.png')
	plt.show()
	plt.clf()
	
	
	#Fit each EPFL spectrum (one for each temperature) to a single Gaussian function plus a linear background
	
	if 	exists("DEBUGFILE") == False :
		os.makedirs('DEBUGFILE')
		
	E0  	= np.zeros(dat.shape[1])
	A 		= np.zeros(dat.shape[1])
	fwhm 	= np.zeros(dat.shape[1])
	myfile = open ("fit_results.csv","w")
	myfile.write ("Fitting data: Tension, Time, <QW>, Lorentz1, Lorentz2")
	myfile.close()
	for i in range(dat.shape[1]):	
		[Ai,e0i,y0i] 	= FitGaussian(wl[:,0], dat[:,i], DEBUGFILE='%s_%s_G.png'%(f[i][0:2],f[i][30:-4]))
		print "Fitting data Gaussian: %0.1e 1/s\tPeak energy: %0.2f nm" % (Ai, e0i)
		A[i] 	= Ai
		E0[i] 	= e0i
		
		[Ai1,e0i1,y1i] 	= FitLorrentz1(wl[:,0], y0i, DEBUGFILE='%s_%s_L1.png'%(f[i][0:2],f[i][30:-4]))
		print "Fitting data Lorrentz1: %0.1e 1/s\tPeak energy: %0.2f nm" % (Ai1, e0i1)
		[Ai2,e0i2,y12] 	= FitLorrentz2(wl[:,0], y0i, DEBUGFILE='%s_%s_L2.png'%(f[i][0:2],f[i][30:-4]))
		print "Fitting data Lorrentz2: %0.1e 1/s\tPeak energy: %0.2f nm" % (Ai2, e0i2)
		ti='%s'%(t[i])
	
		myfile = open ("fit_results.csv","a")
		myfile.write ("\n%s_%s;%s;%0.3f;%0.3f" % (f[i][0:2],f[i][30:-4],ti[-8:], e0i1, e0i2))
		myfile.close()
	
	
	
	
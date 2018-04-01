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


def FitGaussian(x, y, xmin=430, xmax=470, DEBUGFILE=False):

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
		error 	=  y - A*np.exp( -np.power((x-x0)/(math.sqrt(2)*sig), 2)) - b
		return sum(error**2)
	#Do the fit
	x0 		= np.array([a0,x0,sig0,m,b])
	xbound = [(0,100000),(440,450),(4,50),(None,None),(None,None)]
	xfit 	= spop.fmin_l_bfgs_b(ObjFxn, x0, bounds=xbound, approx_grad=True)[0]

	#Define the fit function
	def FitFxn(xf):
		A 		= np.abs(xf[0])
		x0 		= xf[1]
		sig 	= xf[2]
		m 		= xf[3]
		b 		= xf[4]
		return lambda x: A*np.exp( -np.power((x-x0)/(math.sqrt(2)*sig), 2)) + b		
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
	return [xfit[0],xfit[1],fit(xm)]

def FitLorrentz1(x, y, xmin=445, xmax=447, DEBUGFILE=False):

	ii 		= np.nonzero(np.logical_and(x >= xmin, x <= xmax))[0]
	xm		= x
	ym		= y
	x 		= x[ii]
	y 		= y[ii]
	a0 		= np.max(y)
	imx 	= np.argmax(y)
	x0 		= x[imx]
	ilo 	= imx - 0.5 
	ihi 	= imx + 0.5
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
		return y - A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) - b
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
			return lambda x: A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) + b	
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

	ii 		= np.nonzero(np.logical_and(x >= xmin, x <= xmax))[0]
	xm		= x
	ym		= y
	x 		= x[ii]
	y 		= y[ii]
	a0 		= np.max(y)
	imx 	= np.argmax(y)
	x0 		= x[imx]
	ilo 	= imx - 0.5 
	ihi 	= imx + 0.5
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
		return y - A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) - b
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
			return lambda x: A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) + b	
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

def FitLorrentz3(x, y, xmin=461, xmax=463, DEBUGFILE=False):

	ii 		= np.nonzero(np.logical_and(x >= xmin, x <= xmax))[0]
	xm		= x
	ym		= y
	x 		= x[ii]
	y 		= y[ii]
	a0 		= np.max(y)
	imx 	= np.argmax(y)
	x0 		= x[imx]
	ilo 	= imx - 0.5 
	ihi 	= imx + 0.5
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
		return y - A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) - b
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
			return lambda x: A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) + b	
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

def FitLorrentz4(x, y, xmin=467, xmax=469, DEBUGFILE=False):

	ii 		= np.nonzero(np.logical_and(x >= xmin, x <= xmax))[0]
	xm		= x
	ym		= y
	x 		= x[ii]
	y 		= y[ii]
	a0 		= np.max(y)
	imx 	= np.argmax(y)
	x0 		= x[imx]
	ilo 	= imx - 0.5 
	ihi 	= imx + 0.5
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
		return y - A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) - b
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
			return lambda x: A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) + b	
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
		w			= []
		dat			= []
		
		for row in rd:
			
			w.append(float(row[0]))
			dat.append(float(row[1]))
			
		w 		= np.array(w)
		dat 	= np.array(dat)

		return [w, dat]

if __name__=="__main__":
	#Increase the matplotlib font size 
	font = {'family' : 'normal', 'size'   : 16}
	matplotlib.rc('font', **font)
	
	
	#Read EPFL's data from the file	
	[w,dat] 	= ParseFile('data.dat')
	
	plt.plot (w, dat, 'ko',ms=2, color='black')
	plt.xlabel("Longueur d'onde [nm]")
	plt.ylabel("Intensite [1/s]")
	plt.savefig('plot_init.png')
	plt.show()
	plt.clf()
	

	#Fit each spectrum		
	[A1,e01,y1] 	= FitLorrentz1(w, dat, DEBUGFILE='a.png')
	print "Fitting data Lorrentz1: %0.1e 1/s\tPeak energy: %0.2f nm" % (A1, e01)
	[A2,e02,y2] 	= FitLorrentz2(w, dat, DEBUGFILE='b.png')
	print "Fitting data Lorrentz2: %0.1e 1/s\tPeak energy: %0.2f nm" % (A2, e02)	
	[A,e0,y0] 	= FitGaussian(w, dat-y1-y2, DEBUGFILE='c.png')
	print "Fitting data Gaussian: %0.1e 1/s\tPeak energy: %0.2f nm" % (A, e0)
	[A3,e03,y3] 	= FitLorrentz3(w, dat, DEBUGFILE='d.png')
	print "Fitting data Lorrentz1: %0.1e 1/s\tPeak energy: %0.2f nm" % (A3, e03)
	[A4,e04,y4] 	= FitLorrentz4(w, dat, DEBUGFILE='e.png')
	print "Fitting data Lorrentz1: %0.1e 1/s\tPeak energy: %0.2f nm" % (A4, e04)

	ii1 		= np.nonzero(np.logical_and(w >= 444.8, w <= 446))[0]
	ii2 		= np.nonzero(np.logical_and(w >= 454, w <= 456))[0]
	#ii			= np.nonzero(np.logical_and(w >= xmin, w <= xmax))[0]
	ii3 		= np.nonzero(np.logical_and(w >= 461.5, w <= 462.2))[0]
	ii4 		= np.nonzero(np.logical_and(w >= 467.5, w <= 468))[0]
	
	w1 			= w[ii1]
	dat1 		= dat[ii1]
	
	w2 			= w[ii2]
	dat2 		= dat[ii2]	
	
	w3 			= w[ii3]
	dat3 		= dat[ii3]

	w4 			= w[ii4]
	dat4 		= dat[ii4]

	
	#Plot
	plt.plot (w, dat, 'ko',ms=2, color='black')
	plt.plot (w1,dat1, lw=2, color='blue')
	plt.plot (w2,dat2, lw=2, color='red')
	plt.plot (w3,dat3, lw=2, color='green')
	plt.plot (w4,dat4, lw=2, color='purple')
	plt.xlabel("Longueur d'onde [nm]")
	plt.ylabel("Intensite [1/s]")
	plt.savefig('plot_fit.png')
	plt.show()
	plt.clf()
	
	
	
	
import numpy as np
import csv
import math
import matplotlib.pyplot as plt 
import matplotlib 
import scipy.optimize as spop
import scipy.interpolate as spip
import sys
import os
from os.path import exists

h 		= 4.136e-15 #[eV-s]
c 		= 3e10 #[cm/s]
nm2cm 	= 1e-7

def FitLorrentz1(x, y, xmin=720, xmax=750, DEBUGFILE=False):
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
		return y - A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) - m*x - b
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
			return lambda x: A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) + m*x + b
		xf  	= np.linspace(xmin,xmax,1000)
		fit 	= FitFxn(xfit)		
		plt.plot(xm,ym,'ko',ms=4,label='Data')
		plt.plot(xf,fit(xf),'r-',lw=1.5,label='Fit')
		plt.legend(loc='upper left')
		plt.xlabel('Energy [eV]')
		plt.ylabel('Count rate [1/s]')
		plt.xlim(xmin,xmax)
		plt.ylim(0,a0*1.2)	
		plt.savefig('DEBUGFILE/%s'% (DEBUGFILE))		
		plt.clf()
		
	#Return the fit parameters
	return [xfit[0],xfit[1],fit(xm)]

def FitLorrentz2(x, y, xmin=555, xmax=572, DEBUGFILE=False):

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
		return y - A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2))) - m*x - b
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
			return lambda x: A*np.power(gam,2)/((np.power((x-x0), 2) + np.power(gam, 2)))+ m*x + b
		xf  	= np.linspace(xmin,xmax,1000)
		fit 	= FitFxn(xfit)		
		plt.plot(xm,ym,'ko',ms=4,label='Data')
		plt.plot(xf,fit(xf),'r-',lw=1.5,label='Fit')
		plt.legend(loc='upper left')
		plt.xlabel('Energy [eV]')
		plt.ylabel('Count rate [1/s]')
		plt.xlim(xmin,xmax)
		plt.ylim(0,a0*1.2)	
		plt.savefig('DEBUGFILE/%s'% (DEBUGFILE))
		plt.clf()
		
	#Return the fit parameters
	return [xfit[0],xfit[1],fit(xm)]	

def ParseFile(f):
	with open(f, 'r') as csvfile:
		rd = csv.reader(csvfile, delimiter='	')
		row 		= next(rd)
		x			= []
		y			= []
		eV  		= []
		dati 		= []
		for row in rd:
			x.append(float(row[0]))
			y.append(float(row[1]))
			eV.append(float(row[2]))
			dati.append(float(row[3]))
		x 			= np.array(x)
		y 			= np.array(y)
		eV 			= np.array(eV)
		dati 		= np.array(dati)
		
		last_i = 0
		j	= 0

		all_data = np.column_stack((x,y))
		m = len({(x,y) for x,y in all_data})
		
		for i in range (dati.shape[0]):
			
			is_last_reading = i == dati.shape[0]-1
			if x[i]!=x[last_i] or y[i]!=y[last_i] or is_last_reading:
				
				if last_i == 0:	#on selectionne le raman shift mais la premiere fois est suffisante
					dat = np.zeros((i+2, m)) #dat est de dimension 2, contrairement a dati
					eVj = eV[:i+2]

				dat[:-2, j] = dati[last_i:i + is_last_reading]
				dat[0,j]	= x[last_i]
				dat[1,j]	= y[last_i]

				last_i = i
				j += 1
			
	
		return [m, eVj, dat]
 
		
if __name__=="__main__":
	#Increase the matplotlib font size 
	font = {'family' : 'normal', 'size'   : 16}
	matplotlib.rc('font', **font)
	
	#Read EPFL's data from the file	
	[m, eV, dat] 	= ParseFile("B1_12_488_x22_y2_1um.txt")
	
	
	#Plot the EPFL spectra
	plt.plot(eV, dat)
	plt.xlabel('Raman shift [1/cm]')
	plt.ylabel('Count rate [1/s]')
	plt.xlim(550,800)
	plt.ylim(0,18000)
	plt.grid()
	plt.savefig('raman.png')
	plt.show()
	plt.clf()
	
	if 	exists("DEBUGFILE") == False :
		os.makedirs('DEBUGFILE')
	
	#Fit each EPFL spectrum (one for each temperature) to a single Gaussian function plus a linear background
	myfile = open ("fitresultB1_12_xy.csv","w")
	myfile.write ("fit result\nFitting data: x, y, aA1, x0A1, aE2, x0E2")
	myfile.close()
	for i in range(m):	
		xi = dat[0,i]
		yi = dat[1,i]
		[A1i,x10i,gam1i] 	= FitLorrentz1(eV, dat[:,i] , DEBUGFILE='mes1_%d.png'%(i))
		[A2i,x20i,gam2i] 	= FitLorrentz2(eV, dat[:,i] , DEBUGFILE='mes2_%d.png'%(i))
		print("Fitting data: x=%0.1f, y=%0.1f, %0.3e,%0.3e,%0.3e,%0.3e" % (xi, yi, A1i, x10i, A2i, x20i))
		myfile = open ("fitresultB1_12_xy.csv","a")
		myfile.write ("\n%0.1f,%0.1f,%0.3e,%0.3f,%0.3e,%0.3f" % (xi, yi, A1i, x10i, A2i, x20i))
		myfile.close()
	

	
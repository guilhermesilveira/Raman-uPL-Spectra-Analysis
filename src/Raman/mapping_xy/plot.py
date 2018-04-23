import numpy as np
import csv
import math
import matplotlib.pyplot as plt 
import matplotlib 
import scipy.optimize as spop
import scipy.interpolate as spip
import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import os
from os.path import exists


def ParseFile(f):
	with open(f, 'rb') as csvfile:
		rd = csv.reader(csvfile, delimiter=',')
		row 		= rd.next()
		row 		= rd.next()
		x  			= []
		y  			= []
		A1			= []
		E2			= []
		for row in rd:
			x.append(float(row[0]))
			y.append(float(row[1]))
			A1.append(float(row[3]))
			E2.append(float(row[5]))
		x 		= np.array(x)
		y 		= np.array(y)
		A1 		= np.array(A1)
		E2 		= np.array(E2)
		return [x, y, A1, E2]
	
def plot(f):
	#Increase the matplotlib font size 
	font = {'family' : 'normal', 'size'   : 16}
	matplotlib.rc('font', **font)
	
	#Read  from the file	
	[x, y, A1, E2] 	= ParseFile(f)
	
	#construire Xf et Yf
	
	x0 	= x[0]
	y0 	= y[0]
	xf = [x0]
	yf = [y0]


	for i in range(x.shape[0]):
		if x[i]!=x0 and y[i]==y[0]:
			xf.append(x[i])
			x0=x[i]
			
		if y[i]!=y0:
			yf.append(y[i])
			y0=y[i]	
		
				
	xf 		= np.array(xf)
	yf 		= np.array(yf)

	'''	
	xf 		= np.arange(-5,5.5,0.5)
	yf 		= xf
	'''
	
	zdataA1p = np.zeros((xf.shape[0],yf.shape[0]))
	zdataE2p = np.zeros((xf.shape[0],yf.shape[0]))

	zdataA1 = np.zeros((yf.shape[0],xf.shape[0]))
	zdataE2 = np.zeros((yf.shape[0],xf.shape[0]))
	
	#construire z (2D array)
	for j in range(xf.shape[0]):
		zdataA1p[j,0] = A1[j]
		zdataA1p[j,1] = A1[j]
		zdataE2p[j,0] = E2[j]
		zdataE2p[j,1] = E2[j]
	
	#inverser ligne et colone
	for k in range(xf.shape[0]):
		for l in range(yf.shape[0]):
			zdataA1[l,k] = zdataA1p[k,l]
			zdataE2[l,k] = zdataE2p[k,l]
	
	#plot the graph
	levelsA1 = np.linspace(730, 740, 40)
	cs1 = plt.contourf(xf, yf, zdataA1, levels=levelsA1)
	plt.colorbar(cs1, format="%.2f")
	plt.title('%sA1.png'% (f[:-7]))
	plt.xlim(-20,20)
	plt.ylim(-20,20)	
	plt.savefig('%sA1.png'% (f[:-7]))
	plt.show()
	plt.clf()
	
	levelsE2 = np.linspace(565, 569, 40)
	cs2 = plt.contourf(xf, yf, zdataE2, levels=levelsE2)
	plt.colorbar(cs2, format="%.2f")
	plt.title('%sE2.png'% (f[:-7]))
	plt.xlim(-20,20)
	plt.ylim(-20,20)
	plt.savefig('%sE2.png'% (f[:-7]))
	plt.show()
	plt.clf()	

if __name__=="__main__":
		plot('fitresultB1_12_xy.csv')


# -*- coding: utf-8 -*-
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

	
def ParseFile(f):
	with open(f, 'rb') as csvfile:
		rd = csv.reader(csvfile, delimiter='	')
		row 		= rd.next()
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
	[w,dat] 	= ParseFile('data.txt')



		
#Plot
plt.semilogy(w, dat, lw=2, color='b')
plt.xlabel(u'Déplacement Raman [1/cm]')
plt.ylabel(u'Intensité [1/s]')	
plt.xlim(450,750)
plt.savefig('raman.png')
plt.show()
plt.clf()
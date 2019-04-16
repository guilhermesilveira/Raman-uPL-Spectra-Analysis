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
		rd = csv.reader(csvfile, delimiter=';')
		row 		= rd.next()
		V			= []
		wl1  		= []
		wl2 		= []
		t			= []
		
		
		myfile = open ("test.csv","w")
		myfile.write ("test")
		myfile.close()
		
		for row in rd:
			Vr		= row[0]
			if Vr[-1:]	== 'V':
				Vr2=Vr[:-1]
			if Vr[-1:]	!= 'V':
				Vr2=Vr[:-3]
			V.append(Vr2)
			
			tl		= row[1]
			tr		= float(tl[-2:]) + float(tl[-5:-3])*60 + float(tl[:-6])*3600
			
			t.append(float(tr))
			
			wl1.append(float(row[2]))
			wl2.append(float(row[3]))
			
		wl1 		= np.array(wl1)
		wl2 		= np.array(wl2)
		t 			= np.array(t)
		
		Va			= []
		ta			= []
		wl1a 		= []
		wl2a		= []
		Vb			= []
		tb			= []
		wl1b  		= []
		wl2b 		= []
		
		for i in range(wl1.shape[0]):
			if V[i][0:1] == 'A':
				Va.append(float(V[i][3:]))
				wl1a.append(float(wl1[i]))
				wl2a.append(float(wl2[i]))
				ta.append(float(t[i]))
			if V[i][0:1] == 'B':
				Vb.append(float(V[i][3:]))
				wl1b.append(float(wl1[i]))
				wl2b.append(float(wl2[i]))
				tb.append(float(t[i]))
				
		wl1a 		= np.array(wl1a)
		wl2a 		= np.array(wl2a)
		wl1b 		= np.array(wl1b)
		wl2b 		= np.array(wl2b)
		Va			= np.array(Va)
		Vb			= np.array(Vb)
		
		return [Va, ta-np.min(ta), wl1a, wl2a, Vb, tb-np.min(tb), wl1b, wl2b]

def Plot(x,y,xlim,ylim,name,color):	
	return 1
		
if __name__=="__main__":
	#Increase the matplotlib font size 
	font = {'family' : 'normal', 'size'   : 16}
	matplotlib.rc('font', **font)
	
	#Read EPFL's data from the file	
	[Va, ta, wl1a, wl2a, Vb, tb, wl1b, wl2b] 	= ParseFile('fit_results.csv')
	
	i=np.argsort(ta)
	ta=ta[i]
	wl1a=wl1a[i]
	wl2a=wl2a[i]
	Va=Va[i]
	
	
	j=np.argsort(tb)
	tb=tb[j]
	wl1b=wl1b[j]
	wl2b=wl2b[j]
	Vb=Vb[j]
	
	
	COLORS = ['b','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'w', 'BlueViolet', 'Crimson', 'ForestGreen','ForestGreen', 'Indigo', 'Tomato', 'Maroon', 'b', 'g', 'r', 'c', 'm', 'y','ForestGreen', 'k', 'w', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
	

		
#Plot Volt et Time


f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(ta, Va, "b:o", color='black')
axarr[0].set_ylim(-230,230)
axarr[0].set_xlim(-70,2300)
axarr[0].set_ylabel('tension [V]')
axarr[0].set_title('A_Wl1')
axarr[1].plot(ta, wl1a-wl1a[0], "b:o", color='b')
axarr[1].plot(ta, wl2a-wl2a[0], "b:o", color='r')
axarr[1].set_ylabel('shift [nm]')
axarr[1].set_xlabel('time [s]')
plt.savefig('A1.png')
plt.show()
plt.clf()


f, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(tb, Vb, "b:o", color='black')
axarr[0].set_ylim(-230,230)
axarr[0].set_xlim(-70,4250)
axarr[0].set_ylabel('tension [V]')
axarr[0].set_title('B_Wl1')
axarr[1].plot(tb, wl1b-wl1b[0], "b:o", color='b')
axarr[1].plot(tb, wl2b-wl2b[0], "b:o", color='r')
axarr[1].set_ylabel('shift [nm]')
axarr[1].set_xlabel('time [s]')
plt.savefig('B1.png')
plt.show()
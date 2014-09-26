import coro as ok
import csv
import numpy
import scipy
import scipy.interpolate
import scipy.optimize
import pygame
import matplotlib.pyplot as plt
import mdp
import pylab as pl
import os
import sys
import scipy.cluster.hierarchy as cluster
import main as mk

files = [f for f in os.listdir('.') if os.path.isfile(f)]
print '1.1 - ReadingFiles:'
Data={}
Twep= ''
for f in files:
	if f[0]=='_':
		Data[f]=mk.boolienize(mk.clip(mk.readFileC(f)))		# Mark files for input with an | as the first charactor in their name
		print f
		if f[1]=='_':
			Twep=f					# Mark the TableWithEndPoints with ||

DataT= [mk.getTypes(mk.getRowNames(Data[k])) for k in files types[k]]

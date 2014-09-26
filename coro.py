#########################################################
##   Author : Kevin Nelson                             ## 
#	Code is property of Nelson Lab @ Uconn          #
##	Analysis of Gene Expression in Re-programming  ##
#########################################################

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
#data={}
#Data={}

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#fileIN=sys.argv[1]
#Data={}
#with open(fileIN) as csvfile:
#	tableMain=csv.reader(csvfile,delimiter='\t')
#	k=0
#	for row in tableMain:
#		Data[k]=row
#		k+=1							#Read the Data


def clip(D,booler):
	BM96={}
	k = 12
	while(is_number(D[k][2])):
		k+=1
	b=k
	n=0
	for k in range(11,b):					#SHEET SPECIFIC
		row=[]
		for p in range (1,98,1):			#98 well plate
			if(is_number(D[k][p])):
				if booler:
					if(float(D[k][p])>0 and float(D[k][p])!=999):
						row.append(1)
					else:
						row.append(0)
				else:
					if(float(D[k][p])!=999):
						row.append(D[k][p])
					else:
						row.append(0)
			else:
				row.append(D[k][p])
		BM96[k-11]=row
	#print '\#'
	#print BM96
	return BM96
#data={}
#nata={}
#data=clip(Data,True)						#data -> 0,1
#nata=clip(Data,False)						#Data 25 Delta CT *
#print data



def getColNames(data):
	names=[""]*(len(data[0]))
	for k in range(0,len(data[0]),1):
		names[k]=data[0][k]
	#print names
	return(names)
#Genes=getColNames(data)

def getRowNames(data):
	names=[""]*(len(data))
	for k in range(0,len(data),1):
		names[k]=data[k][0]
	#print names
	return(names)
#samples=getRowNames(data)

def npClip(data):						#Remove gene names => numpy.array
	so=numpy.zeros((len(data)-1,len(data[0])-1))	
	for k in range(1,len(data)):
		for k2 in range(1,len(data[1])):
			so[k-1][k2-1]=float(data[k][k2])
	return(so)
#d=npClip(data)
#n=npClip(nata)

def npClip2(data):						#Remove gene names => numpy.array
	so=numpy.zeros((len(data)-1,len(data[0])-3))	
	for k in range(1,len(data)):
		for k2 in range(1,len(data[1])-2):
			so[k-1][k2-1]=float(data[k][k2])
	return(so)

def npClip3(data):						#Remove gene names => numpy.array
	#print data
	#print len(data[1])-3
	so=numpy.zeros((len(data),len(data[0])-1))	
	for k in range(0,len(data)):
		for k2 in range(1,len(data[0])):
			if data[k][k2]=='':
				so[k][k2-1]=0.
			else:
				so[k][k2-1]=float(data[k][k2])
	return(so)
#print n

def eucDistMat(n,condensed):
	Sample=numpy.zeros((len(n),len(n)))
	Gene=numpy.zeros((len(n[0]),len(n[0])))
	for k in range(0,len(n)):
		for k2 in range(0,len(n)):
			if k2 > k or condensed == False:
				dist=0
				for k3 in range(0,len(n[0])):
					dist+=(n[k][k3]-n[k2][k3])**2
				dist=dist**0.5
				Sample[k][k2]=dist
					
	for k in range(0,len(n[0])):
		for k2 in range(0,len(n[0])):
			dist=0
			if k2 > k or condensed == False:
				for k3 in range(0,len(n)):
					dist+=(n[k3][k]-n[k3][k2])**2
				dist=dist**0.5
			Gene[k][k2]=dist
	return (Sample,Gene)

def StdMat(n):
	Std=numpy.zeros(len(n))
	Gtd=numpy.zeros(len(n[0]))
	for k in range(0,len(n)):
		dev=0
		m=sum(n[k])/len(n[0])
		for k2 in range(0,len(n[0])):
			dev+=(n[k][k2]-m)**2
		Std[k]=dev**0.5
					
	for k in range(0,len(n[0])):
		dev=0.0
		m=0.0
		for k2 in range(0,len(n)):
			m+=n[k2][k]
		m=m/len(n[0])
		for k2 in range(0,len(n)):
			dev+=(n[k2][k]-m)**2
		Gtd[k]=dev**0.5
		
	return (Std,Gtd)
#[std,gtd]=StdMat(n)

#def coroMat(n,d):
#	[stdN,gtdN]=StdMat(n)
#	[stdB,gtdB]=StdMat(d)
#	[sdmB,gdmB]=eucDistMat(d,True)
#	[sdmN,gdmN]=eucDistMat(n,True)

#	gcmN=numpy.zeros((len(n[0]),len(n[0])))
#	gcmB=numpy.zeros((len(n[0]),len(n[0])))
#	scmN=numpy.zeros((len(n),len(n)))
#	scmB=numpy.zeros((len(n),len(n)))
#	
#	gcmN0=numpy.zeros((len(n[0]),len(n[0])))
#	gcmB0=numpy.zeros((len(n[0]),len(n[0])))
#	scmN0=numpy.zeros((len(n),len(n)))
#	scmB0=numpy.zeros((len(n),len(n)))
#	def rxy(n):
def CORR(n):
	[r,c]=StdMat(n)
	#n0=numpy.mean(n,axis=0)
	nm=numpy.mean(n,axis=1)
	mn=numpy.mean(n,axis=0)
	T=N=numpy.zeros((len(n),len(n)))
	#print len(n)
	#print len(n[0])
	M=N=numpy.zeros((len(n[0]),len(n[0])))
	t=0.
	for k in range(0,len(n),1):
		for k2 in range(0,len(n),1):
			t=0.					
			if k>k2:
				for k3 in range(0,len(n[0])):
					t+=(n[k][k3]-nm[k])*(n[k2][k3]-nm[k2])
				T[k][k2]=t/(r[k]*r[k2])
#	for k in range(0,len(n),1):
#		for k2 in range(0,len(n),1):
#			T[k][k2]=T[k][k2]/(r[k]*r[k2])
	for k in range(0,len(n[0]),1):
		for k2 in range(0,len(n[0]),1):
			t=0.					
			if k>k2:
				for k3 in range(0,len(n)):
					t+=(n[k3][k]-mn[k])*(n[k3][k2]-mn[k2])
				M[k][k2]=t/(c[k]*c[k2])
#	for k in range(0,len(n[0]),1):
#		for k2 in range(0,len([0]),1):
#			M[k][k2]=M[k][k2]/(c[k]*r[k2])
	return (T,M)

#		N=numpy.zeros((len(n),len(n)))	
#		for k in range(0,len(n),1):
#			for k2 in range(0,len(n),1):
#			N[k][k2]=
#			
#			
#		
#	n0=numpy.mean(n,axis=0)
#	n1=numpy.mean(n,axis=1)
#	
#	d0=numpy.mean(d,axis=0)
#	d1=numpy.mean(d,axis=1)
#	

#	for k in range(0,len(n)):
#		for k2 in range(0,len(n)):
#			if ((stdN[k]-stdN[k2])!= 0.):
#				scmN[k][k2]=len(n)*sdmN n[] n0[k]*[k2]/((stdN[k]-stdN[k2])**2)
#			else:
#				scmN[k][k2]=sdmN[k][k2]
#				scmN0[k][k2]=sdmN[k][k2]
#			if ((stdB[k]-stdB[k2])!=0.):			
#				scmB[k][k2]=sdmB[k][k2]/((stdB[k]-stdB[k2])**2)
#			else:
#				scmB[k][k2]=sdmB[k][k2]
#				scmB0[k][k2]=sdmB[k][k2]
#					
#	for k in range(0,len(n[0])):
#		for k2 in range(0,len(n[0])):
#			if (gtdN[k2]-gtdN[k2])!=0:
#				gcmN[k][k2]=gdmN[k][k2]/((gtdN[k2]-gtdN[k2])**2)
#			else:
#				gcmN[k][k2]=gdmN[k][k2]/((gtdN[k2]-gtdN[k2])**2)
#				gcmN0[k][k2]=gdmN[k][k2]/((gtdN[k2]-gtdN[k2])**2)
#			if (gtdB[k2]-gtdB[k2])!=0:
#				gcmB[k][k2]=gdmB[k][k2]/((gtdB[k2]-gtdB[k2])**2)
#			else:
#				gcmB[k][k2]=gdmB[k][k2]/((gtdB[k2]-gtdB[k2])**2)
#				gcmB0[k][k2]=gdmB[k][k2]/((gtdB[k2]-gtdB[k2])**2)			
#	return (scmN,scmB,gcmN,gcmB,scmN0,scmB0,gcmN0,gcmB0)



def octPrt(sdmB,name):
	op=name+'=['
	for k in range(0,len(sdmB)):
		for k2 in range(0,len(sdmB[0])):
			op += str(sdmB[k][k2])
			if k2+1 == len(sdmB[0]):
				op+=';'
			else:
				op+=','
	op+=']\n'
	return(op)
	
def BCM(d,n):
	sdmB,gdmB = CORR(d)
	sdmN,gdmN = CORR(n)
	#[scmN,scmB,gcmN,gcmB,scmN0,scmB0,gcmN0,gcmB0]=coroMat(n,d)
	op=''
	op+=octPrt(sdmB,'scmB')
	op+=octPrt(gdmB,'gcmB')
	op+=octPrt(sdmN,'scmN')
	op+=octPrt(gdmN,'gcmN')

#	op=''
#	op+=octPrt(smB0,'scmB0')
#	op+=octPrt(gmB0,'gcmB0')
#	op+=octPrt(smN0,'scmN0')
#	op+=octPrt(gmN0,'gcmN0')

	#output=open('c0r0.m',"w")
	#output.write(op)


	#print sdmB
	#print gdmB
	#print sdmN
	#print gdmN
	return(op)





#output=open('c0r0.m',"w")
#op=''
#op+='\nfunction m()\nmesh(1:length(n),1:length(n),n)\nendfunction'
#output.write(BCM(d,n)+op)

def CorrScat(a,b):
	#print 'a'
	#print a
	#print 'b'
	#print b
	A1={}
	B1={}
	n=0;
	for k in range(1,len(a[0])):
		for k2 in range(1,len(b[0])):
			if a[0][k]==b[0][k2]:
				row=[]
				for k3 in range(0,len(a[1])):
					row.append(a[k3][k])
				A1[n]=row
				row=[]
				for k3 in range(0,len(b[1])):
					row.append(b[k3][k2])
				B1[n]=row
				n+=1
				#print n
	#print A1
	#print B1
	A2,A3=CORR(npClip3(A1))
	B2,B3=CORR(npClip3(B1))
	A=numpy.zeros((len(A2)*(len(A2)-1))/2)
	B=numpy.zeros((len(A2)*(len(A2)-1))/2)
	n=0
	for k in range(0,len(A2)):
		for k2 in range(0,len(A2[0])):
			if(k>k2):
				A[n]=A2[k][k2]
				B[n]=B2[k][k2]
				n+=1
	return A,B


					
					





















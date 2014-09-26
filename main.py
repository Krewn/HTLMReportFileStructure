#########################################################
##   Author : Kevin Nelson                             ## 
#					                #
##	Analysis of Gene Expression in Re-programming  ##
#########################################################
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
#data={}
#Data={}

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

def is_bool(s):
	try:
		bool(s)
		return True
	except ValueError:
		return False

#def is_not(s):
#	try:
#		str(s)
#		return True
#	except ValueError:
#		return False


def readFileT(n):
	fileIN=sys.argv[n]
	Data={}
	with open(fileIN) as csvfile:
		tableMain=csv.reader(csvfile,delimiter='\t')
		k=0
		for row in tableMain:
			Data[k]=row
			k+=1							#Read the Data into dict Data
	return(Data)
def readFileC(n):
	try:
		fileIN=sys.argv[n]
	except(TypeError):
		fileIN=n
	Data={}
	with open(fileIN) as csvfile:
		tableMain=csv.reader(csvfile,delimiter=',')
		k=0
		for row in tableMain:
			Data[k]=row
			k+=1							#Read the Data into dict Data
	return(Data)
#Data=readFile(1)
def prDictStr(Dists):
	once = True
	Head="\t"
	tail=""
	n=0
	for k in Dists:
		if once == True and n>0:
			try:
				for k2 in Dists[k]:
					Head+=str(k2)+'\t'
			except(TypeError):
				Head=str(k)
			once =False
		#Head+=str(k)+"\t"
		tail+=str(k)+"\t"
		n2 = 0
		try:
			for k2 in Dists[k]:
				try:
					tail+=str(Dists[k][k2])+"\t"
				except(TypeError):
					tail+=str(Dists[k][n2])+"\t"
				n2+=1
			tail+="\n"
		except(TypeError):
			pass
		n+=1
	tail=tail[0:len(tail)-1]
	Head+="\n"
	return(Head+tail)
def prDict(Dists):
	once = True
	Head="\t"
	tail=""
	n=0
	for k in Dists:
		if once == True and n>0:
			try:
				for k2 in Dists[k]:
					Head+=str(k2)+'\t'
			except(TypeError):
				Head=str(k)
			once =False
		#Head+=str(k)+"\t"
		tail+=str(k)+"\t"
		n2 = 0
		try:
			for k2 in Dists[k]:
				try:
					tail+=str(Dists[k][k2])+"\t"
				except(TypeError):
					tail+=str(Dists[k][n2])+"\t"
				n2+=1
			tail+="\n"
		except(TypeError):
			pass
		n+=1
	tail=tail[0:len(tail)-1]
	Head+="\n"
	print Head+tail
	pass

def deBlank(n):
	for k in range(0,len(n)):
		for k2 in range(0,len[0]):
			if is_not(n[k][k2]):
				n[k,k2]=0.


def clip(D):
	#print '\n.\n..\n...\n....\n...\n..\n.\n'
	#prDict2(D)
	BM96={}
	n=0
	b=12
	while(is_number(D[b][2])):
		b+=1
	#print b
	for k in range(11,b):
		row=[]
		for p in range (1,98,1):
			if(is_number(D[k][p])):
				if(float(D[k][p])>0 and float(D[k][p])!=999):
					row.append(float(D[k][p]))
				else:
					row.append(0)
			else:
				row.append(D[k][p])
		BM96[k-11]=row
	#print '\#'
	#print BM96
	return BM96
#data={}
#data=clip(Data)
#print data



def getColNames(data):
	try:
		names=[""]*(len(data[0]))
		for k in range(0,len(data[0]),1):
			names[k]=data[0][k]
		#print names
		return(names)
	except(KeyError):
		ky = data.keys()
		names=[""]*(len(data[ky[0]]))
		for k in range(0,len(data[k2]),1):
			names[k]=data[k2][k]
		#print names
		return(names)
#Genes=getColNames(data)

def getRowNames(data):
	try:
		names=[""]*(len(data))
		for k in range(0,len(data),1):
			names[k]=data[k][0]
		#print names
		return(names)
	except(KeyError):
		names=[""]*(len(data))
		n=0
		for k in data:
			for k2 in data[k]:
				break
			names[n]=k2
			n+=1
		#print names
		return(names)
		
#samples=getRowNames(data)

#samples=getRowNames(data)

#def npClip(data):						#Remove gene names => numpy.array
#	so=numpy.zeros((len(data)-1,len(data[1])-1))	
#	for k in range(1,len(data)):
#		for k2 in range(1,len(data[1])):
#			so[k-1][k2-1]=float(data[k][k2])
#	return(so)
#nata=npClip(data)

def boolienize(data):
	BData={}
	n=0
	for k in range(0,len(data),1):					#boolienize Data
		row=[]
		for p in data[k]:
			if(is_number(p)):
				if(float(p)>0):
					row.append(1)
				else:
					row.append(0)
			else:
				row.append(p)
		BData[k]=row
	return BData

def Boolienize(data):
	BData={}
	n=0
	try:
		for k in range(0,len(data),1):					#boolienize Data
			row=[]
			for p in data[k]:
				if(is_number(p)):
					if(float(p)>0 and float(p)!=999):
						row.append(1)
					else:
						row.append(0)
				else:
					row.append(p)
			BData[k]=row
			#print BData
		return BData
	except(KeyError):
		for k in data:					#boolienize Data
			row=[]
			for p in data[k]:
				if(is_number(p)):
					if(float(p)>0 and float(p)!=999):
						row.append(1)
					else:
						row.append(0)
				else:
					row.append(p)
			BData[k]=row
			#print BData
		return BData

#BData = boolienize(data)
#print BData
#day=[0]*len(data)
#samples=getRowNames(data)
#print samples

def getdays(samples):								#Parsing based on predetermined naming convention.
	day=[0]*len(samples)							#
	for k in range(1,len(samples)-1,1):
		for k2 in range(0,len(samples[k])):
			if samples[k][k2]=='D':
				#print samples[k]
				if is_number(samples[k][k2+1]+samples[k][k2+2]):			
					day[k]=float(samples[k][k2+1]+samples[k][k2+2])
				elif is_number(samples[k][k2+1]):
					day[k]=flaot(samples[k][1])
	#print '##'	
	#print (day)	
	return(day)
#days=getdays(samples)

def getTypes(samples):
	Type=['']*len(samples)
	for k in range(0,len(samples),1):
		for k2 in range(0,len(samples[k]),1):
			if samples[k][k2]=='B' and samples[k][k2+1]=='J':
				Type[k]='Fib'
			if samples[k][k2]=='H' and (samples[k][k2+1]=='1' or samples[k][k2+1]=='9'):
				Type[k]='hESC'
			if samples[k][k2]=='S' and samples[k][k2+5]=='-' and samples[k][k2+14]=='-':
				Type[k]='--'
			if samples[k][k2]=='S' and samples[k][k2+5]=='+' and samples[k][k2+14]=='-':
				Type[k]='+-'
			if samples[k][k2]=='S' and samples[k][k2+5]=='+' and samples[k][k2+14]=='+':
				Type[k]='++'
	return (Type)
def ltd(l):
		n=0
		d={}
		for k in l:
			d[n]=k
			n+=1
		return(d)
def dtl(d):
		n=0
		l=['']*len(d)
		for k in d:
			l[n]=k
			n+=1
		return(l)
def getST(samples):			
	S=['']*len(samples)
	T=['']*len(samples)
	hesc=['']*len(samples)
	fib = ['']*len(samples)
	for k in range(0,len(samples),1):
		for k2 in range(0,len(samples[k]),1):
			if samples[k][k2]=='B' and samples[k][k2+1]=='J':
				fib[k]='+'
			if samples[k][k2]=='H' and (samples[k][k2+1]=='1' or samples[k][k2+1]=='9'):
				hesc[k]='+'
			if samples[k][k2]=='S' and samples[k][k2+5]=='-' and samples[k][k2+14]=='-':
				S[k]='-'
				T[k]='-'
			if samples[k][k2]=='S' and samples[k][k2+5]=='+' and samples[k][k2+14]=='-':
				S[k]='+'
				T[k]='-'
			if samples[k][k2]=='S' and samples[k][k2+5]=='+' and samples[k][k2+14]=='+':
				S[k]='+'
				T[k]='+'
	#print S, T, hesc, fib
	a = {'SEA4':ltd(S),'TRA160':ltd(T),'hesc':ltd(hesc),'fib':ltd(fib)}
	return (a)

def getTypesMP(samples):
	Type=getTypes(samples)
	TypeMP= ['']*len(samples)
	for k in range(0,len(samples),1):
		for k2 in range(0,len(samples[k]),1):
			if samples[k][k2]=='M' and samples[k][k2+1]=='O' and samples[k][k2+2]=='N' and samples[k][k2+3]=='O':
				TypeMP[k]='mono'
	for k in range(0,len(samples)):
		if TypeMP[k] != 'mono':
			Type[k] = ''
	#print len(Type)
	#print len(TypeMP)
	for k in range(0,len(samples)):
			try:
				if TypeMP[k] != 'mono' and Type[k] != 'Fib' and Type[k] != 'hESC':
					TypeMP[k]='poly'
			except(IndexError):
				TypeMP[k]='poly'
	days=getdays(samples)
	stats=getST(samples)
	return ({'ID':Type,'MvP':TypeMP,'day':days,'SEA4':stats['SEA4'],'TRA160':stats['TRA160'],'hESC':stats['hesc'],'Fib':stats['fib']})
					
#Type=getTypes(samples)
#print Type
def GetTypes2(n):
	Type=['']*len(n)
	for k in range(0,len(n)):
		if(n[k][len(n[0])-1]=='Fib'):
			Type[k]='Fib'
		if(n[k][len(n[0])-1]=='CDH1'):
			Type[k]='CDH1'
		if(n[k][len(n[0])-1]=='GFP'):
			Type[k]='GFP'
		if(n[k][len(n[0])-1]=='SSE4'):
			Type[k]='+-'
		if(n[k][len(n[0])-1]=='SSEA4 and Tra1-60'):
			Type='++'
		if(n[k][len(n[0])-1]=='hESC'):
			Type[k]='hESC'
		if(n[k][len(n[0])-1]=='SSEA4 and Tra1-60 Low GFP'):
			Type='++lgfp'
		if(n[k][len(n[0])-1]=='SSEA4 Low GFP'):
			Type[k]='+-lgfp'
	return(Type)
		


def EUC (BData,Type):
	#print Type
	#print BData
	fibAvg   =[0]*len(BData[0])
#	for k in BData:
#		if (len(BData[k])!= len(BData[0])):
#			print 'the fuck mang' 
#		else: print k
	fibCount =0
	hESCAvg  =[0]*len(BData[0])
	hESCCount=0
	#BData[1][len(BData[1])-1]
	#ke = BData.keys()
#	print ke
#	print len (BData)
#	print len (BData[ke[0]])
	for k in range(0,len(BData),1):						#Establish FIB and hESC average
		if(Type[k]=="Fib"):
			#print 'fib'
			fibCount+=1
			for k2 in range(0,len(BData[1])-1,1):
				if(is_number(BData[k][k2])):			
					fibAvg[k2]+=float(BData[k][k2])		
		if(Type[k]=="hESC"):
			#print 'hesc'
			hESCCount+=1
			for k2 in range(0,len(BData[k])-1,1):
				if(is_number(BData[k][k2])):
					hESCAvg[k2]+=float(BData[k][k2])
	#print fibCount
	#print hESCCount
	for k in range(0,len(fibAvg),1):					#fibAvg & hESCAvg have the same length
		fibAvg[k]=fibAvg[k]/fibCount
		hESCAvg[k]=hESCAvg[k]/hESCCount


	eucDistFib=[0]*len(BData)
	eucDisthESC=[0]*len(BData)
	fd=0
	hd=0

	for k in range(0,len(BData),1):
		fd=0
		hd=0
		for k2 in range(0,len(BData[k]),1):
			if(is_number(BData[k][k2])):
				fd+=(BData[k][k2]-fibAvg[k2])*(BData[k][k2]-fibAvg[k2])
				hd+=(BData[k][k2]-hESCAvg[k2])*(BData[k][k2]-hESCAvg[k2])
		eucDistFib[k]=numpy.sqrt(fd)
		eucDisthESC[k]=numpy.sqrt(hd)
	return(eucDistFib,eucDisthESC)

def PCA(data):
	try:
		NBD=numpy.zeros((len(data)-1,len(data[0])-1))
		for k in range(1,len(data),1):
			row=[]
			for k1 in range(1,len(data[0]),1):
				if(is_number(data[k][k1])):
					row.append(float(data[k][k1]))
			NBD[k-1]=row
	except(ValueError):
		NBD=numpy.zeros((len(data)-1,len(data[0])-1))
		n=0
		for k in data:
			row=[]
			for k1 in data[k]:
				if(is_number(k1)):
					row.append(float(k1))
			if len(row) == 98:
				NBD[n]=row
			else:
				print k
			n+=1
	pca=mdp.pca(NBD,svd=True)
	l1=numpy.zeros(len(pca))
	l2=numpy.zeros(len(pca))
	for k in range(0,len(pca),1):
		l1[k]=pca[k][0]
		l2[k]=pca[k][1]
	return l1,l2

def PCAfull(data):
	NBD=numpy.zeros((len(data)-1,len(data[0])-1))
	for k in range(1,len(data),1):
		row=[]
		for k1 in range(1,len(data[0]),1):
			if(is_number(data[k][k1])):
				row.append(float(data[k][k1]))
		NBD[k-1]=row
	pca=mdp.pca(NBD,svd=True)
	return pca

#print eucDistFib
#print eucDisthESC
#[eucDistFib,eucDisthESC] = EUC(BData,Type)

def GetTimes(eucDistFib,eucDisthESC,BData):
	FibFib=0								#Fibroblas center on Fibdist
	FibEsc=0								#Fibroblas center on HescDist
	EscFib=0								#Hesc center on Fib
	EscEsc=0								#Hesc senter on Hesc

	fibCount =0
	hESCCount=0
	Type = getTypes(getRowNames(BData))

	for k in range (0,len(BData),1):
		if Type[k]=="Fib" :
			fibCount+=1
			FibFib+=eucDistFib[k]
			FibEsc+=eucDisthESC[k]	
		if Type[k]=="hESC":
			hESCCount+=1
			EscFib+=eucDistFib[k]
			EscEsc+=eucDisthESC[k]

	EscFib/=hESCCount
	EscEsc/=hESCCount
	FibFib/=fibCount
	FibEsc/=fibCount
	#print EscFib
	#print EscEsc
	#print FibFib
	#print FibEsc
	FibD=EscFib-FibFib							#component on fibroblast axis of vector T in Euc space
	EscD=EscEsc-FibEsc							#component on Hesc axis of vector T in Euc space
	magT=numpy.power((numpy.power((FibD),2)+numpy.power((EscD),2)),0.5)

	times=[]
	for k in range(0,len(eucDistFib)):
		if is_number(eucDistFib[k]):
			times.append((eucDistFib[k]*FibD+eucDisthESC[k]*EscD)/magT)
		else:
			print "mmk"
	return(times)

def GetTimesNormalized(eucDistFib,eucDisthESC,BData):
	FibFib=0								#Fibroblas center on Fibdist
	FibEsc=0								#Fibroblas center on HescDist
	EscFib=0								#Hesc center on Fib
	EscEsc=0								#Hesc senter on Hesc

	fibCount =0
	hESCCount=0
	Type = getTypes(getRowNames(BData))

	for k in range (0,len(BData),1):
		if Type[k]=="Fib" :
			fibCount+=1
			FibFib+=eucDistFib[k]
			FibEsc+=eucDisthESC[k]	
		if Type[k]=="hESC":
			hESCCount+=1
			EscFib+=eucDistFib[k]
			EscEsc+=eucDisthESC[k]

	FibFib=FibFib/fibCount
	FibEsc=FibEsc/fibCount
	EscFib=EscFib/hESCCount
	EscEsc=EscEsc/hESCCount

	#print EscFib
	#print EscEsc
	#print FibFib
	#print FibEsc
	FibD=EscFib-FibFib							#component on fibroblast axis of vector T in Euc space
	#print FibD
	EscD=FibEsc-EscEsc							#component on Hesc axis of vector T in Euc space
	#print EscD
	magT=numpy.power((numpy.power(FibD,2)+numpy.power(EscD,2)),0.5)
	#print FibFib
	#print EscEsc
	#print magT

	times=[]
	for k in range(0,len(eucDistFib)):
		if is_number(eucDistFib[k]):
			times.append(((eucDistFib[k]-FibFib)*FibD+(FibEsc-eucDisthESC[k])*EscD)/(magT*magT))
		else:
			print "mmk"
	return(times)

#def euc():
#	pass

def GetTimesNormalizedkMeansMaximumVarienceBetweenCentroids(eucDistFib,eucDisthESC,BData,nmeans):
	FibFib=0								#Fibroblas center on Fibdist
	FibEsc=0								#Fibroblas center on HescDist
	EscFib=0								#Hesc center on Fib
	EscEsc=0								#Hesc senter on Hesc
	fibCount =0
	hESCCount=0								#Recomended numbe rof Means 3, creating 4 groups as fibs are also calculated for reference.
	Type = getTypes(getRowNames(BData))
	means={}
	for k in range(0,nmeans):
		temp=[]
		for k in range(0,96):
			if(numpy.random.ranf()>0.5):temp.append(1)
			else:temp.append(0)
		means[k]=temp
	Cdata=npClip(BData)
	def kmeans0(Cdata,means,itter):
		fmeans=means
		dist=[[0.]*nmeans]*len(Cdata)
		for k in range(0,len(Cdata)):		
			for k2 in range(0,nmeans):
				for k3 in range(0,len(Cdata[k])):
					if(is_number(Cdata[k][k2])):
						dists[k2][k]+=numpy.power(Cdata[k][k3]-means[k2][k3],2)
				dists[k2][k]=numpy.sqrt(dists[k2][k])
		#deesDemDoss z og
		clust=[]
		for k in range(0,len(Cdata)):
			temp=9000
			n=-1
			for k2 in range(0,nmeans):
				if(dists[k2][k]<temp):
					temp=dists[k2][k]
					n=k2
			clust.append[n]
		for k in range(0,nmeans):
			temp=[0.]*len(Cdata[0])
			n=0
			for k2 in range(0,len(Cdata)):
				if(clust[k2]==k):
					for k3 in range(0,len(Cdata[k2])):
						temp[k3]+=Cdata[k2][k3]
					n+=1
			if(n!=0):
				for k3 in range(0,len(Cdata[k2])):
					temp[k3]/=n
			means[k]=temp
		itter-=1
		repeat=False
		for k in range(0,len(fmeans)):
			for k2 in len(fmeans[k]):
				if(fmeans[k][k2]!=means):
					repeat=True
					break
		if(itter<0 or repeat==False):return(means)
		else:return(kmeans0(Cdata,means,itter))
	means=kmeans0(Cdata,means,10)
	dists='INcomplete'
			 

	for k in range (0,len(BData),1):
		if Type[k]=="fib" :
			fibCount+=1
			FibFib+=eucDistFib[k]
			FibEsc+=eucDisthESC[k]	
		if Type[k]=="end":
			hESCCount+=1
			EscFib+=eucDistFib[k]
			EscEsc+=eucDisthESC[k]

	FibFib=FibFib/fibCount
	FibEsc=FibEsc/fibCount
	EscFib=EscFib/hESCCount
	EscEsc=EscEsc/hESCCount

	#print EscFib
	#print EscEsc
	#print FibFib
	#print FibEsc
	FibD=EscFib-FibFib							#component on fibroblast axis of vector T in Euc space
	#print FibD
	EscD=FibEsc-EscEsc							#component on Hesc axis of vector T in Euc space
	#print EscD
	magT=numpy.power((numpy.power(FibD,2)+numpy.power(EscD,2)),0.5)
	#print FibFib
	#print EscEsc
	#print magT

	times=[]
	for k in range(0,len(eucDistFib)):
		if is_number(eucDistFib[k]):
			times.append(((eucDistFib[k]-FibFib)*FibD+(FibEsc-eucDisthESC[k])*EscD)/(magT*magT))
		else:
			print "mmk"
	return(times)
#times = GetTimes(eucDistFib,eucDisthESC,BData)


#freqESC=numpy.zeros((steps,len(BData[0])-3))
#freqFib=numpy.zeros((steps,len(BData[0])-3))
#freqtime=numpy.zeros((steps,len(BData[0])-3))

def Scat(eucDisthESC,eucDistFib):
	scat = plt.figure()
	p1 = scat.add_subplot(111)

	p1.scatter(eucDisthESC[1:],eucDistFib[1:],color='blue', s=7,edgecolor='none')
	p1.set_aspect(1./p1.get_data_ratio())
	p1.grid(True)

	plt.show()

def Scat2(eucDisthESC,eucDistFib,bdata,i,anote,invertX,invertY):
	plt.figure()
	x1=[]
	x2=[]
	y1=[]
	y2=[]
	try:
		for k in range(1,len(bdata),1):
			if (bdata[k][i]==1):
				y1.append(eucDisthESC[k])
				x1.append(eucDistFib[k])
			else:
				y2.append(eucDisthESC[k])
				x2.append(eucDistFib[k])
	except(TypeError):
		for k in range(0,len(bdata[0])):
			if(i==bdata[0][k]):
				i=k
				break				
		for k in range(1,len(bdata),1):
			if (bdata[k][i]==1):
				y1.append(eucDisthESC[k])
				x1.append(eucDistFib[k])
			else:
				y2.append(eucDisthESC[k])
				x2.append(eucDistFib[k])
	X1=numpy.zeros(len(x1))
	Y1=numpy.zeros(len(y1))
	X2=numpy.zeros(len(x2))
	Y2=numpy.zeros(len(y2))
	for k in range(0,len(x1),1):
		X1[k]=x1[k]
	for k in range(0,len(y1),1):
		Y1[k]=y1[k]
	for k in range(0,len(x2),1):
		X2[k]=x2[k]
	for k in range(0,len(y2),1):
		Y2[k]=y2[k]
	scat = plt.figure()
	p1 = scat.add_subplot(111)
	p1.scatter(X1,Y1,color='red', s=7,edgecolor='none')
	p1.scatter(X2,Y2,color='blue', s=7,edgecolor='none')
	plt.xlabel('EucDist from BJ')
	plt.ylabel('EucDist from hESC')
	if invertX:
		pl.gca().invert_xaxis()
	if invertY:
		pl.gca().invert_yaxis()
	p1.set_aspect(1./p1.get_data_ratio())
	p1.grid(True)
	plt.savefig(bdata[0][i]+anote+'.png')
	plt.clf()
	pass

def Scat3(eucDisthESC,eucDistFib,bdata,i,anote,invertX,invertY):
	plt.figure()
	x1=[]
	x2=[]
	y1=[]
	y2=[]
	try:
		for k in range(0,len(eucDisthESC),1):
			if (bdata[k+1][i]==1):
				y1.append(eucDisthESC[k])
				x1.append(eucDistFib[k])
			else:
				y2.append(eucDisthESC[k])
				x2.append(eucDistFib[k])
	except(TypeError):
		for k in range(0,len(bdata[0])):
			if(i==bdata[0][k]):
				i=k
				break
		for k in range(0,len(eucDisthESC),1):
			if (bdata[k+1][i]==1):
				y1.append(eucDisthESC[k])
				x1.append(eucDistFib[k])
			else:
				y2.append(eucDisthESC[k])
				x2.append(eucDistFib[k])	
	X1=numpy.zeros(len(x1))
	Y1=numpy.zeros(len(y1))
	X2=numpy.zeros(len(x2))
	Y2=numpy.zeros(len(y2))
	for k in range(0,len(x1),1):
		X1[k]=x1[k]
	for k in range(0,len(y1),1):
		Y1[k]=y1[k]
	for k in range(0,len(x2),1):
		X2[k]=x2[k]
	for k in range(0,len(y2),1):
		Y2[k]=y2[k]
	scat = plt.figure()
	p1 = scat.add_subplot(111)
	p1.scatter(X1,Y1,color='red', s=7,edgecolor='none')
	p1.scatter(X2,Y2,color='blue', s=7,edgecolor='none')
	plt.xlabel('PrinComp 1')
	plt.ylabel('PrinComp 2')
	if invertX:
		pl.gca().invert_xaxis()
	if invertY:
		pl.gca().invert_yaxis()
	p1.set_aspect(1./p1.get_data_ratio())
	p1.grid(True)
	plt.savefig(bdata[0][i]+anote+'.png')
	plt.clf()
	pass

def ScatA(eucDisthESC,eucDistFib,bdata,i,anote,invertX,invertY,title,n123):
	if n123==1 :
		scat = plt.figure()
		p1 = scat.add_subplot(111)

		p1.scatter(eucDisthESC[1:],eucDistFib[1:],color='blue', s=7,edgecolor='none')
		p1.set_aspect(1./p1.get_data_ratio())
		p1.grid(True)
		plt.suptitle(title)
		plt.show()
	elif n123 == 2 :
		plt.figure()
		x1=[]
		x2=[]
		y1=[]
		y2=[]
		try:
			for k in range(1,len(bdata),1):
				if (bdata[k][i]==1):
					y1.append(eucDisthESC[k])
					x1.append(eucDistFib[k])
				else:
					y2.append(eucDisthESC[k])
					x2.append(eucDistFib[k])
		except(TypeError):
			for k in range(0,len(bdata[0])):
				if(i==bdata[0][k]):
					i=k
					break				
			for k in range(1,len(bdata),1):
				if (bdata[k][i]==1):
					y1.append(eucDisthESC[k])
					x1.append(eucDistFib[k])
				else:
					y2.append(eucDisthESC[k])
					x2.append(eucDistFib[k])
		X1=numpy.zeros(len(x1))
		Y1=numpy.zeros(len(y1))
		X2=numpy.zeros(len(x2))
		Y2=numpy.zeros(len(y2))
		for k in range(0,len(x1),1):
			X1[k]=x1[k]
		for k in range(0,len(y1),1):
			Y1[k]=y1[k]
		for k in range(0,len(x2),1):
			X2[k]=x2[k]
		for k in range(0,len(y2),1):
			Y2[k]=y2[k]
		scat = plt.figure()
		p1 = scat.add_subplot(111)
		p1.scatter(X1,Y1,color='red', s=7,edgecolor='none')
		p1.scatter(X2,Y2,color='blue', s=7,edgecolor='none')
		plt.xlabel('EucDist from BJ')
		plt.ylabel('EucDist from hESC')
		if invertX:
			pl.gca().invert_xaxis()
		if invertY:
			pl.gca().invert_yaxis()
		p1.set_aspect(1./p1.get_data_ratio())
		p1.grid(True)
		plt.suptitle(title)
		plt.savefig(bdata[0][i]+anote+'.png')
		plt.clf()
	elif n123 == 3:
		plt.figure()
		x1=[]
		x2=[]
		y1=[]
		y2=[]
		try:
			for k in range(0,len(eucDisthESC),1):
				if (bdata[k+1][i]==1):
					y1.append(eucDisthESC[k])
					x1.append(eucDistFib[k])
				else:
					y2.append(eucDisthESC[k])
					x2.append(eucDistFib[k])
		except(TypeError):
			for k in range(0,len(bdata[0])):
				if(i==bdata[0][k]):
					i=k
					break
			for k in range(0,len(eucDisthESC),1):
				if (bdata[k+1][i]==1):
					y1.append(eucDisthESC[k])
					x1.append(eucDistFib[k])
				else:
					y2.append(eucDisthESC[k])
					x2.append(eucDistFib[k])	
		X1=numpy.zeros(len(x1))
		Y1=numpy.zeros(len(y1))
		X2=numpy.zeros(len(x2))
		Y2=numpy.zeros(len(y2))
		for k in range(0,len(x1),1):
			X1[k]=x1[k]
		for k in range(0,len(y1),1):
			Y1[k]=y1[k]
		for k in range(0,len(x2),1):
			X2[k]=x2[k]
		for k in range(0,len(y2),1):
			Y2[k]=y2[k]
		scat = plt.figure()
		p1 = scat.add_subplot(111)
		p1.scatter(X1,Y1,color='red', s=7,edgecolor='none')
		p1.scatter(X2,Y2,color='blue', s=7,edgecolor='none')
		plt.xlabel('PrinComp 1')
		plt.ylabel('PrinComp 2')
		if invertX:
			pl.gca().invert_xaxis()
		if invertY:
			pl.gca().invert_yaxis()
		p1.set_aspect(1./p1.get_data_ratio())
		p1.grid(True)
		plt.suptitle(title)
		plt.savefig(bdata[0][i]+anote+'.png')
		plt.clf()

steps = 200
window=numpy.sqrt(48.0)/20
step=numpy.sqrt(48.0)/steps

def windowed(dists,dat,window,steps,step):
	center = 0.0
	inWin=0.0
	k0=0
	freq={}
	#inWin=0.0
	while center < numpy.sqrt(48.0):
		#print inWin
		#print center
		freq[k0]=numpy.zeros(len(dat[0]))
		inWin=0.0
		for k in range(0,len(dists)):
			if((dists[k]>center-window) & (dists[k]<center+window)):
				inWin+=1.0
				for k1 in range(0,len(dat[0])):
					if is_number(dat[k][k1]):
						if(dat[k][k1]!=0):
							freq[k0][k1]+=1.0
		if inWin!=0:
			for k1 in range(0,len(dat[0])):
				freq[k0][k1]=freq[k0][k1]/inWin
		center+=step
		k0+=1
		#print k0
	return(freq)

def windowedT(dists,dat,window,steps,step):
	center = -3.0
	inWin=0.0
	k0=0
	freq={}
	#inWin=0.0
	while center < 3.0:
		#print inWin
		#print center
		freq[k0]=numpy.zeros(len(dat[0]))
		inWin=0.0
		for k in range(0,len(dists)):
			if((dists[k]>center-window) & (dists[k]<center+window)):
				inWin+=1.0
				for k1 in range(0,len(dat[0])):
					if is_number(dat[k][k1]):
						if(dat[k][k1]!=0):
							freq[k0][k1]+=1.0
		if inWin!=0:
			for k1 in range(1,len(dat[0])):
				freq[k0][k1]=freq[k0][k1]/inWin
		center+=step
		k0+=1
		#print k0
	return(freq)
#print "here we go"
def HeatLable(Names,ax):
	for k in range(0,len(Names),1):
		annotate(Names[k],xy=(k,0),xycoords='data',xytext=(k,k*-1),textcoords='offset points')

#testFib = windowed(eucDistFib,BData,window,steps,step)
#print testFib
#testESC = windowed(eucDisthESC,BData,window,steps,step)
#print testESC
#testTime= windowedT(times,BData,window,steps,step)
#print(testTime)

def findOrder(test):
	total=numpy.zeros(len(test[0]))
	for k in range(0,len(test[0])):
		tot=0.0
		for k2 in range(0,len(test)):
			tot+=test[k2][k]
		total[k]=tot
	sortedTotal=numpy.sort(total)
	place=[0]*len(total)	
	for k in range (0,len(total)):
		for k2 in range (0,len(total)):
			if (sortedTotal[k]==total[k2]):
				place[k]=k2
	return(place)

def byPlace(test,place):
	sbp=test
	for k in range(0,len(test[0])):
		for k2 in range(0,len(test)):
			sbp[k2][k]=test[k2][place[k]]
	return(sbp)

#def fit1(test,window,steps,step):
#	a=numpy.zeros(len(test))
#	b=numpy.zeros(len(test))
#	for k in range(0,len(test),1):
#		a[k]=k*step
#	fits={}
#	for k in range(0,len(test[0])):
#		for k2 in range(0,len(test)):
#			b[k2]=test[k2][k]
#		fits[k]=scipy.interpolate.interp1d(a,b,kind='linear')
#	return(fits)


#def fit2(test,window,steps,step):
#	a=numpy.zeros(len(test))
#	b=numpy.zeros(len(test))
#	for k in range(0,len(test),1):
#		a[k]=k*step
#	fits={}
#	for k in range(0,len(test[0])):
#		print k
#		for k2 in range(0,len(test)):
#			b[k2]=test[k2][k]
#		fits[k]=scipy.interpolate.interp1d(a,b,kind='quadratic')
#	return(fits)


#def fit3(test,window,steps,step):
#	a=numpy.zeros(len(test))
#	b=numpy.zeros(len(test))
#	for k in range(0,len(test),1):
#		a[k]=k*step
#	fits={}
#	for k in range(0,len(test[0])):
#		print k
#		for k2 in range(0,len(test)):
#			b[k2]=test[k2][k]
#		fits[k]=scipy.interpolate.interp1d(a,b,kind='cubic')
#	return(fits)

#print "mmmm"
#escL=fit1(testESC,window,steps,step)
#print "."
#escP=fit2(testESC,window,steps,step)
#print "."
#escC=fit3(testESC,window,steps,step)
#print ".."
#fibL=fit1(testFib,window,steps,step)
#print "."
#fibP=fit2(testFib,window,steps,step)
#print "."
#fibC=fit3(testFib,window,steps,step)
#print "..."
#timeL=fit1(testTime,window,steps,step)
#print "."
#timeP=fit2(testTime,window,steps,step)
#print "."
#timeC=fit3(testTime,window,steps,step)
#print "k"										#Takes way to long... Not Worth it


def sigmoid(p,x):
    x0,y0,c,k=p
    y = c / (1 + numpy.exp(-k*(x-x0))) + y0
    return y

def residuals(p,x,y):
    return y - sigmoid(p,x)

def sigfit(x,y):
	p_guess=(numpy.median(x),numpy.median(y),1.0,1.0)
	p, cov, infodict, mesg, ier = scipy.optimize.leastsq(residuals,p_guess,args=(x,y),full_output=1)
	x0,y0,c,k=p
	return(x0,y0,c,k)

def FindMinWin(test):
	minWin=0
	Top=False
	First=True
	for k in range(0,len(test),1):
		for k2 in range(1,len(test[0]),1):
			if test[k][k2]!=0:
				if First:
					minWin=k
					First=False
					break;
					break;
	return (minWin)

def FindMaxWin(test):
	maxWin=0
	Top=False
	First=True
	for b in range(1,len(test)+1,1):
		k=len(test)-b
		for k2 in range(0,len(test[0]),1):
			if test[k][k2]!=0:
				if First:
					maxWin=k
					First=False
					break;
					break;
	return(maxWin)

def GetFits(test,step):
	
	maxWin=0
	minWin=0
	Bottom=False
	Top=False
	First=True
	for k in range(0,len(test),1):
		for k2 in range(0,len(test[0]),1):
			if test[k][k2]!=0:
				if First:
					minWin=k
					First=False
	First=True
	for b in range(1,len(test),1):
		k=len(test)-b
		for k2 in range(0,len(test[0]),1):
			if test[k][k2]!=0:
				if First:
					maxWin=k
					First=False
	dat=numpy.zeros((maxWin-minWin)+1)
	base=numpy.zeros((maxWin-minWin)+1)
	n1=0	
	for k2 in range(minWin,maxWin+1,1):
		base[n1]=k2*step
		n1+=1
	n1=0
	sigFits={}
	for k in range(0,len(test[0]),1):
		n1+=1		
		n=0
		for k2 in range(minWin,maxWin+1,1):
			dat[n]=test[k2][k]
			n+=1
		sigFits[n1]=sigfit(base,dat)
	return(sigFits)

#fitsFib=GetFits(testFib,step)
#fitsESC=GetFits(testESC,step)
#fitsTime=GetFits(testTime,step)

def kSort(test,taken):
	minWin=FindMinWin(test)
	maxWin=FindMaxWin(test)
	so = numpy.zeros((maxWin-minWin,len(test[0])))
	for k in range(0,len(taken)-1,1):
		for k2 in range(minWin,maxWin,1):
			so[k2-minWin][k]=test[k2][taken[k]]
	return(so)
def nSort(names,taken):
	No=[""]*len(names)
	for k in range(0,len(names),1):
		No[k]=names[take[k]]

def SortFit(test,fits):
	minWin=FindMinWin(test)
	maxWin=FindMaxWin(test)
	ordered=numpy.zeros((maxWin-minWin,len(test[0])))
	m=32000.
	n=-1
	taken=[-1]*len(fits)
	so = numpy.zeros((maxWin-minWin,len(test[0])))
	for k0 in range(0,len(fits),1):	
		m=32000.
		for k in range(1,len(fits),1):
			if(taken.count(k)==0):
				if(fits[k][0]<m):
					m=fits[k][0]
					n=k
		taken[k0]=n
	so=kSort(test,taken)
	return(so)
							
		
	#for k in range(1,)

#print fitsFib
#print fitsESC
#print fitsTime

########################################################################################
def tkSortFit(test,fits):
	minWin=FindMinWin(test)
	maxWin=FindMaxWin(test)
	ordered=numpy.zeros((maxWin-minWin,len(test[0])))
	m=32000.
	n=-1
	taken=[-1]*len(fits)
	#so = numpy.zeros((maxWin-minWin,len(test[0])))
	for k0 in range(0,len(fits),1):	
		m=32000.
		for k in range(1,len(fits),1):
			if(taken.count(k)==0):
				if(fits[k][0]<m):
					m=fits[k][0]
					n=k
		taken[k0]=n
	#so=kSort(test,taken)
	return(taken)
							
		
	#for k in range(1,)

#print fitsFib
#print fitsESC
#print fitsTime

def tkSortDifs(test):
	minWin=FindMinWin(test)
	maxWin=FindMaxWin(test)
	delt=[0.]*len(test[0])
	for k in range(0,len(test[0])):
		for k2 in range(1,len(test)):
			if(test[k2][k]-test[k2-1][k]>0):
				delt[k]+=test[k2][k]-test[k2-1][k]
	
	ordered=numpy.zeros((maxWin-minWin,len(test[0])))
	m=32000.
	n=-1
	taken=[-1]*len(test[0])
	#so = numpy.zeros((maxWin-minWin,len(test[0])))
	for k0 in range(0,len(delt),1):	
		m=32000.
		for k in range(1,len(delt),1):
			if(taken.count(k)==0):
				if(delt[k]<m):
					m=delt[k]
					n=k
		taken[k0]=n
	#so=kSort(test,taken)
	return(taken)

def tkSortOne(test):
	minWin=FindMinWin(test)
	maxWin=FindMaxWin(test)
	delt=[-1]*len(test[0])
	for k in range(0,len(test[0])):
		for k2 in range(2,len(test)-1,1):
			if(test[k2][k]==1 & (bool(test[k2+1][k]!=numpy.float64(1))|bool(test[k2-1][k]!=numpy.float64(1)))):
				delt[k]=k2
				break
	ordered=numpy.zeros((maxWin-minWin,len(test[0])))
	m=32000.
	n=-1
	taken=[-1]*len(test[0])
	#so = numpy.zeros((maxWin-minWin,len(test[0])))
	for k0 in range(0,len(delt),1):	
		m=32000.
		for k in range(1,len(delt),1):
			if(taken.count(k)==0):
				if(delt[k]<m):
					m=delt[k]
					n=k
		taken[k0]=n
	#so=kSort(test,taken)	
	return(taken)

def tkSortOneB(test):
	minWin=FindMinWin(test)
	maxWin=FindMaxWin(test)
	delt=[-1]*len(test[0])
	for k in range(0,len(test[0])):
		for k2 in range(2,len(test)-1,1):
			k3=len(test)-k2
			if(test[k3][k]==1 & (bool(test[k3+1][k]!=numpy.float64(1))|bool(test[k3-1][k]!=numpy.float64(1)))):
				delt[k]=k3
				break
	ordered=numpy.zeros((maxWin-minWin,len(test[0])))
	m=32000.
	n=-1
	taken=[-1]*len(test[0])
	#so = numpy.zeros((maxWin-minWin,len(test[0])))
	for k0 in range(0,len(delt),1):	
		m=32000.
		for k in range(1,len(delt),1):
			if(taken.count(k)==0):
				if(delt[k]<m):
					m=delt[k]
					n=k
		taken[k0]=n
	#so=kSort(test,taken)
	return(taken)

def tkSortHalf(test):
	minWin=FindMinWin(test)
	maxWin=FindMaxWin(test)
	delt=[0.000]*len(test[0])
	deltB=[0.001]*len(test[0])
	deltT=[0.001]*len(test[0])
	taken=[-1]*len(test[0])
	#so = numpy.zeros((maxWin-minWin,len(test[0])))
	for k in range(0,len(test[0]),1):
		bot=0.
		top=0.
		for k2 in range(0,len(test)/2,1):
			bot+=test[k2][k]
		for k3 in range(len(test)/2,len(test),1):
			top+=test[k2][k]
		delt[k]=(0.001+bot)/(0.001+top)
	for k0 in range(0,len(delt),1):	
		m=32000.
		for k in range(1,len(delt),1):
			if(taken.count(k)==0):
				if(delt[k]<m):
					m=delt[k]
					n=k
		taken[k0]=n
	#so=kSort(test,taken)
	return(taken)
	
	

def HeatMap(test):
	heat=numpy.zeros((len(test),len(test[0])))
	for k in range(0,len(test),1):
		for k2 in range(0,len(test[k]),1):
			heat[k][k2]=test[k][k2]
	pl.pcolor(heat)
	pl.colorbar()
	pl.show()

def pgHeatMap(test,genes):
	minWin=FindMinWin(test)
	maxWin=FindMaxWin(test)
	h=500
 	w=20
	pygame.init()
	window = pygame.display.set_mode((w,h))
	for k in range(0,len(test[0]),1):
		window.fill((250,250,250))
		for k2 in range(0,len(test),1):
			cut = (len(test)-minWin)+k2
			t=pygame.Rect(0,int((h/len(test))*k2),20,int((h/len(test))))
			pygame.draw.rect(window,(255*test[k2][k],0,255*(1-test[k2][k])),t,0)
		font = pygame.font.SysFont('none',17,bold=False,italic=False)
		window.blit(pygame.transform.rotate(font.render(genes[k],0,(0,0,0)),90),(5,h-100))
		pygame.display.flip()
		pygame.image.save(window,genes[k]+".png")

#pgHeatMap(testTime,names,r'Time')
def htMap(sort,test,stest,ssort):
		tk=sort(test)
		htmlDoc='<!DOCTYPE html> \n <html> \n'
		f=stest
		for k in range(0,len(tk),1):
			htmlDoc+='<body><a href=G/'+genesA[tk[k]]+'/'+genesA[tk[k]]+'.html><img src="G/'+f+'/'+genesA[tk[k]]+'.png" width="20" height="500"></body>'
		htmlDoc+='</body></html>'
		output=open(f+ssort+'.html',"w")
		output.write(htmlDoc)
		output.close()
def htMap2(sort,test,stest,ssort,genesA):
		tk=sort(test)
		htmlDoc='<!DOCTYPE html> \n <html> \n'
		f=stest
		for k in range(0,len(tk),1):
			htmlDoc+='<body><a href=G/'+genesA[tk[k]]+'/'+genesA[tk[k]]+'.html><img src="G/'+f+'/'+genesA[tk[k]]+'.png" width="20" height="500"></body>'
		htmlDoc+='</body></html>'
		output=open(f+ssort+'.html',"w")
		output.write(htmlDoc)
		output.close()

def HTML(data,Path):
	
	steps = 200
	window=numpy.sqrt(48.0)/20
	step=numpy.sqrt(48.0)/steps
	
	genes=getColNames(data)
	samples=getRowNames(data)
	bData=Boolienize(data)

	Type = getTypes(samples)	
	
	[eucDistFib,eucDisthESC] = EUC(bData,Type)
	times = GetTimes(eucDistFib,eucDisthESC,bData)
	testFib = windowed(eucDistFib,bData,window,steps,step)
	testESC = windowed(eucDisthESC,bData,window,steps,step)
	testTime= windowedT(times,bData,window,steps,step)
	pca1,pca2=PCA(data)
	#print pca1
	#print pca2
	if not os.path.exists(Path): os.makedirs(Path)
	os.chdir(Path)
	if not os.path.exists(r'G'): 
		os.makedirs(r'G')	
	os.chdir(r'G')
	for k in range(1,len(genes),1):
		if not os.path.exists(genes[k]): 
			os.makedirs(genes[k])
	for k in range(1,len(genes),1):
		os.chdir(genes[k])
		Scat2(eucDisthESC,eucDistFib,bData,k,'_EUC',False,False)
		Scat3(pca2,pca1,bData,k,'_PCA',True,True)
		output=open(genes[k]+'.html',"w")
		htmlDoc='<!DOCTYPE html> \n <html> \n'
		htmlDoc+='<body><img src="'+genes[k]+'_EUC.png" width="800" height="600">'
		htmlDoc+='<img src="'+genes[k]+'_PCA.png" width="800" height="600"></body>\n'
		htmlDoc+='<body><img src="../../Fits/'+genes[k]+'_LinearRegression.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Fits/'+genes[k]+'_Normals-1.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Fits/'+genes[k]+'_Normals-2.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Fits/'+genes[k]+'_Normals-3.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Fits/'+genes[k]+'_Uniform.jpg" width="480" height="480"></body>'
		htmlDoc+='<body><img src="../../Violins/'+genes[k]+'.png" width="800" height="600"></body>\n'
		htmlDoc+='</html>'
		output.write(htmlDoc)
		output.close()
		os.chdir("..")
	if not os.path.exists(r'Fib'): os.makedirs(r'Fib')
	if not os.path.exists(r'hESC'): os.makedirs(r'hESC')
	if not os.path.exists(r'Time'): os.makedirs(r'Time')
	os.chdir(r'Fib')
	pgHeatMap(testFib,genes)
	os.chdir("..")
	os.chdir(r'hESC')
	pgHeatMap(testESC,genes)
	os.chdir("..")
	os.chdir(r'Time')
	pgHeatMap(testTime,genes)
	os.chdir("..")
	os.chdir("..")
#def htMap(sort,test,stest,ssort):
#	tk=sort(testFib)
#	htmlDoc='<!DOCTYPE html> \n <html> \n'
#	f=stest
#	for k in range(0,len(tk),1):
#		htmlDoc+='<body><a href=G/'+genes[tk[k]]+'/'+genes[tk[k]]+'.html><img src="G/'+f+'/'+genes[tk[k]]+'.png" width="20" height="500"></body>'
#	htmlDoc+='</body></html>'
#	output=open(f+ssort+'.html',"w")
#	output.write(htmlDoc)
#	output.close()
	htMap(tkSortHalf,testFib,'Fib','SortHalf')
	htMap(tkSortHalf,testESC,'hESC','SortHalf')
	htMap(tkSortHalf,testTime,'Time','SortHalf')

	htMap(tkSortOne,testFib,'Fib','SortOne')
	htMap(tkSortOne,testESC,'hESC','SortOne')
	htMap(tkSortOne,testTime,'Time','SortOne')

	htMap(tkSortOneB,testFib,'Fib','SortOneB')
	htMap(tkSortOneB,testESC,'hESC','SortOneB')
	htMap(tkSortOneB,testTime,'Time','SortOneB')

	htMap(tkSortDifs,testFib,'Fib','SortDifs')
	htMap(tkSortDifs,testESC,'hESC','SortDifs')
	htMap(tkSortDifs,testTime,'Time','SortDifs')
	pass

def HTML2(data1,data1T,data2,data2T,Path):
	
	Agg={}					# Combined the 2 inputs into one table
	AggT={}
	n=0
	for k in data1:
		Agg[n]=data1[k]
		n+=1
	first = True
	for k in data2:
		if first:
			first=False
		else:	
			Agg[n]=data2[k]
			n+=1
	for k in data1T:
		AggT['1'+str(k)]=data1T[k]
	first = True
	for k in data2T:
		AggT['2'+str(k)]=data2T[k]
		n+=1
	#print Agg
#	for k in Agg:
#		Agg[k]=dtl(Agg[k])

	steps = 200
	window=numpy.sqrt(48.0)/20
	step=numpy.sqrt(48.0)/steps
	
	genes1=getColNames(data1)
	genes2=getColNames(data2)
	genesA=getColNames(Agg)

	samples1=getRowNames(data1)
	samples2=getRowNames(data2)
	samplesA=getRowNames(Agg)

	bData1=Boolienize(data1)
	bData2=Boolienize(data2)
	bDataA=Boolienize(Agg)

	Type1 = getTypes(samples1)
	Type2 = getTypes(samples2)
	TypeA = getTypes(samplesA)	
	
	[eucDistFib1,eucDisthESC1] = EUC(bData1,Type1)
	[eucDistFib2,eucDisthESC2] = EUC(bData2,Type2)
	[eucDistFibA,eucDisthESCA] = EUC(bDataA,TypeA)
	
	times1 = GetTimes(eucDistFib1,eucDisthESC1,bData1)
	times2 = GetTimes(eucDistFib2,eucDisthESC2,bData2)
	timesA = GetTimes(eucDistFibA,eucDisthESCA,bDataA)

	testFib1 = windowed(eucDistFib1,bData1,window,steps,step)
	testFib2 = windowed(eucDistFib2,bData2,window,steps,step)
	testFibA = windowed(eucDistFibA,bDataA,window,steps,step)

	testESC1 = windowed(eucDisthESC1,bData1,window,steps,step)
	testESC2 = windowed(eucDisthESC2,bData2,window,steps,step)
	testESCA = windowed(eucDisthESCA,bDataA,window,steps,step)
	
	testTime1= windowedT(times1,bData1,window,steps,step)
	testTime2= windowedT(times2,bData2,window,steps,step)
	testTimeA= windowedT(timesA,bDataA,window,steps,step)

	#print data1
	#print data2
	#print Agg

	pca1_1,pca2_1=PCA(data1)
	pca1_2,pca2_2=PCA(data2)
	pca1_A,pca2_A=PCA(Agg)
	
	#print pca1
	#print pca2
	if not os.path.exists(Path): os.makedirs(Path)
	os.chdir(Path)
	if not os.path.exists(r'G'): 
		os.makedirs(r'G')	
	os.chdir(r'G')
	for k in range(1,len(genes1),1):
		if not os.path.exists(genes1[k]): 
			os.makedirs(genes1[k])
	for k in range(1,len(genes1),1):
		os.chdir(genes1[k])
		Scat2(eucDisthESC1,eucDistFib1,bData1,k,'_EUC_1',False,False)
		Scat2(eucDisthESC2,eucDistFib2,bData2,k,'_EUC_2',False,False)
		Scat2(eucDisthESCA,eucDistFibA,bDataA,k,'_EUC_A',False,False)
		Scat3(pca2_1,pca1_1,bData1,k,'_PCA_1',True,True)
		Scat3(pca2_2,pca1_2,bData2,k,'_PCA_2',True,True)
		Scat3(pca2_A,pca1_A,bDataA,k,'_PCA_A',True,True)

		output=open(genes1[k]+'.html',"w")
		
		htmlDoc='<<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n<html smlns="http://www.w3.org/1999/xhtml">\n'
		htmlDoc+='<head>\n<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n<title>'+genes1[k]+'</title>\n<style type="text/css">\n<!--\nbody {\nfont: 100% Verdana, Arial, Helvetica, sans-serif;\nbackground: #666666;\n	margin: 0; /* its good practice to zero the margin and padding of the body element to account for differing browser defaults */\npadding: 0;\ntext-align: center; /* this centers the container in IE 5* browsers. The text is then set to the left aligned default in the #container selector */\ncolor: #000000;\n}\n.oneColLiqCtr #container {\nwidth: 500px;  /* this will create a container 80% of the browser width */\nbackground: #FFFFFF;\nmargin: 0 auto; /* the auto margins (in conjunction with a width) center the page */\nborder: 1px solid #000000;\ntext-align: left; /* this overrides the text-align: center on the body element. */\n}\n.oneColLiqCtr #mainContent {\npadding: 12px 20px 0 20px; /* remember that padding is the space inside the div box and margin is the space outside the div box */\n}\n-->\n</style>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryCollapsiblePanel.js" type="text/javascript"></script>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/Scripts/swfobject_modified.js" type="text/javascript"></script>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryMenuBar.js" type="text/javascript"></script>\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/CSS/Level3_3.css" rel="stylesheet" type="text/css">\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryCollapsiblePanel.css" rel="stylesheet" type="text/css">\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryMenuBarHorizontal_edited.css" rel="stylesheet" type="text/css">\n<style type="text/css" media="screen">#FlashID {visibility:hidden}</style></head>\n'
		htmlDoc+='<body>'
		
		htmlDoc+='<div id="Eucledian Plots" class="CollapsiblePanel CollapsiblePanelClosed">\n'
		htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">Eucledian Plots</div>\n'
		htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
		htmlDoc+='<img src="'+genes1[k]+'_EUC_1.png" width="800" height="600">\n'
		htmlDoc+='<img src="'+genes1[k]+'_EUC_2.png" width="800" height="600">\n'
		htmlDoc+='<img src="'+genes1[k]+'_EUC_A.png" width="800" height="600">\n'
		htmlDoc+='</div>'
		htmlDoc+='</div>'
		
		htmlDoc+='<div id="PCA Plots" class="CollapsiblePanel CollapsiblePanelClosed">\n'
		htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">PCA Plots</div>\n'
		htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
		htmlDoc+='<img src="'+genes1[k]+'_PCA_1.png" width="800" height="600">\n'
		htmlDoc+='<img src="'+genes1[k]+'_PCA_2.png" width="800" height="600">\n'
		htmlDoc+='<img src="'+genes1[k]+'_PCA_A.png" width="800" height="600">\n'
		htmlDoc+='</div>'
		htmlDoc+='</div>'
		htmlDoc+='</body>'

		htmlDoc+='<div id="FitCurvesPoly" class="CollapsiblePanel CollapsiblePanelClosed">\n'
		htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">FitCurvesPoly</div>\n'
		htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
		htmlDoc+='<body><img src="../../Rstuf/dat1FitCurvesfigures/'+genes1[k]+'_LinearRegression.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/dat1FitCurvesfigures/'+genes1[k]+'_Normals-1.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/dat1FitCurvesfigures/'+genes1[k]+'_Normals-2.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/dat1FitCurvesfigures/'+genes1[k]+'_Normals-3.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/dat1FitCurvesfigures/'+genes1[k]+'_Uniform.jpg" width="480" height="480"></body>'
		htmlDoc+='</div>'
		htmlDoc+='</div>'
		htmlDoc+='<div id="FitCurvesMono" class="CollapsiblePanel CollapsiblePanelClosed">\n'
		htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">FitCurvesMono</div>\n'
		htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
		htmlDoc+='<body><img src="../../Rstuf/dat2FitCurvesfigures/'+genes1[k]+'_LinearRegression.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/dat2FitCurvesfigures/'+genes1[k]+'_Normals-1.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/dat2FitCurvesfigures/'+genes1[k]+'_Normals-2.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/dat2FitCurvesfigures/'+genes1[k]+'_Normals-3.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/dat2FitCurvesfigures/'+genes1[k]+'_Uniform.jpg" width="480" height="480"></body>'
		htmlDoc+='</div>'
		htmlDoc+='</div>'
		htmlDoc+='<div id="FitCurvesAll" class="CollapsiblePanel CollapsiblePanelClosed">\n'
		htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">FitCurvesAll</div>\n'
		htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
		htmlDoc+='<body><img src="../../Rstuf/datAFitCurvesfigures/'+genes1[k]+'_LinearRegression.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/datAFitCurvesfigures/'+genes1[k]+'_Normals-1.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/datAFitCurvesfigures/'+genes1[k]+'_Normals-2.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/datAFitCurvesfigures/'+genes1[k]+'_Normals-3.jpg" width="480" height="480">'
		htmlDoc+='<img src="../../Rstuf/datAFitCurvesfigures/'+genes1[k]+'_Uniform.jpg" width="480" height="480"></body>'
		htmlDoc+='</div>'
		htmlDoc+='</div>'
		#htmlDoc+='<body><img src="../../Violins/'+genes1[k]+'.png" width="800" height="600"></body>\n'
		htmlDoc+='<script type="text/javascript">\n<!--\nvar CollapsiblePanel1 = new Spry.Widget.CollapsiblePanel("Eucledian Plots", {contentIsOpen:false});\nvar CollapsiblePanel2 = new Spry.Widget.CollapsiblePanel("PCA Plots", {contentIsOpen:false});\nvar CollapsiblePanel3 = new Spry.Widget.CollapsiblePanel("FitCurvesPoly", {contentIsOpen:false});\nvar CollapsiblePanel3 = new Spry.Widget.CollapsiblePanel("FitCurvesMono", {contentIsOpen:false});\nvar CollapsiblePanel3 = new Spry.Widget.CollapsiblePanel("FitCurvesAll", {contentIsOpen:false});\nswfobject.registerObject("FlashID");\nvar MenuBar1 = new Spry.Widget.MenuBar("MenuBar1", {imgDown:"SpryAssets/SpryMenuBarDownHover.gif", imgRight:"SpryAssets/SpryMenuBarRightHover.gif"});\n//-->\n</script>\n'
		htmlDoc+='<p>'+genes1[k]+'<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene='+genes1[k]+'>'+genes1[k]+'</a></p>\n'
		htmlDoc+='</body>'
		htmlDoc+='</html>'
		output.write(htmlDoc)
		output.close()
		os.chdir("..")
	if not os.path.exists(r'Fib'): os.makedirs(r'Fib')
	if not os.path.exists(r'hESC'): os.makedirs(r'hESC')
	if not os.path.exists(r'Time'): os.makedirs(r'Time')
	os.chdir(r'Fib')
	pgHeatMap(testFibA,genesA)
	os.chdir("..")
	os.chdir(r'hESC')
	pgHeatMap(testESCA,genesA)
	os.chdir("..")
	os.chdir(r'Time')
	pgHeatMap(testTimeA,genesA)
	os.chdir("..")
	os.chdir("..")

#	output=open(f+ssort+'.html',"w")
#	output.write(htmlDoc)
#	output.close()
	htMap(tkSortHalf,testFibA,'Fib','SortHalf')
	htMap(tkSortHalf,testESCA,'hESC','SortHalf')
	htMap(tkSortHalf,testTimeA,'Time','SortHalf')

	htMap(tkSortOne,testFibA,'Fib','SortOne')
	htMap(tkSortOne,testESCA,'hESC','SortOne')
	htMap(tkSortOne,testTimeA,'Time','SortOne')

	htMap(tkSortOneB,testFibA,'Fib','SortOneB')
	htMap(tkSortOneB,testESCA,'hESC','SortOneB')
	htMap(tkSortOneB,testTimeA,'Time','SortOneB')

	htMap(tkSortDifs,testFibA,'Fib','SortDifs')
	htMap(tkSortDifs,testESCA,'hESC','SortDifs')
	htMap(tkSortDifs,testTimeA,'Time','SortDifs')
	pass

def pgHeatMap2(test,genes,condition):
	minWin=FindMinWin(test)
	maxWin=FindMaxWin(test)
	h=500
 	w=20
	pygame.init()
	window = pygame.display.set_mode((w,h))
	for k in range(0,len(test[0]),1):
		window.fill((250,250,250))
		for k2 in range(0,len(test),1):
			cut = (len(test)-minWin)+k2
			t=pygame.Rect(0,int((h/len(test))*k2),20,int((h/len(test))))
			pygame.draw.rect(window,(255*test[k2][k],0,255*(1-test[k2][k])),t,0)
		font = pygame.font.SysFont('none',17,bold=False,italic=False)
		window.blit(pygame.transform.rotate(font.render(genes[k],0,(0,0,0)),90),(5,h-100))
		pygame.display.flip()
		pygame.image.save(window,condition+genes[k]+".png")
	pass

def htMap3(sort,ssort,data1,data2,dataA,genes,condition):
	tk={1:sort(data1),2:sort(data2),3:sort(dataA)}
	#htmlDoc='<!DOCTYPE html> \n <html> \n'
	htmlDoc='<body>'
	f=condition
	types={1:'Poly',2:'Mono',3:'Agg'}
	if condition == 'Time':
		con2 = 'Progression'
	else:
		con2 = condition
	for k in range(0,len(tk[3]),1):
		for k2 in types:
			print genes[k2][tk[k2][k]]
			print types[k2]
			htmlDoc+='<a href=./G/'+genes[k2][tk[1][k]]+'/'+genes[k2][tk[1][k]]+'.html><img src=./'+condition+'/'+types[k2]+con2+genes[k2][tk[1][k]]+'.png width="20" height="500"></a>'
	htmlDoc+='</body>'
	#htmlDoc+='</body></html>'
	return(htmlDoc)

def HTML3(data1,data1T,data2,data2T,Path):
	Agg={}					# Combined the 2 inputs into one table
	AggT={}
	n=0
	for k in data1:
		Agg[n]=data1[k]
		n+=1
	first = True
	for k in data2:
		if first:
			first=False
		else:	
			Agg[n]=data2[k]
			n+=1
	for k in data1T:
		AggT[str(k)+' 1']=data1T[k]
	first = True
	for k in data2T:
		AggT[str(k)+' 2']=data2T[k]
		n+=1
	#print Agg
#	for k in Agg:
#		Agg[k]=dtl(Agg[k])
	prDict(Agg)
	
	steps = 200
	print 'steps:'+str(steps)
	window=numpy.sqrt(48.0)/20
	print 'window'+str(window)
	step=numpy.sqrt(48.0)/steps
	print 'step'+str(step)
	
	genes1=getColNames(data1)
	prDict(genes1)
	genes2=getColNames(data2)
	prDict(genes2)
	genesA=getColNames(Agg)
	prDict(genesA)
	if genes2!=genes1 or genes1!=genesA:
		print 'Assay misallignmnet'

	samples1=getRowNames(data1)
	print 'samples1\n'+prDictStr(samples1)
	samples2=getRowNames(data2)
	print 'samples2\n'+prDictStr(samples2)
	samplesA=getRowNames(Agg)
	print 'samplesA\n'+prDictStr(samplesA)

	bData1=Boolienize(data1)
	print 'bData1'+prDictStr(bData1)
	bData2=Boolienize(data2)
	print 'bData2'+prDictStr(bData2)
	bDataA=Boolienize(Agg)
	print 'bDataA'+prDictStr(bDataA)
	
	Type1 = getTypes(samples1)
	print 'Types1'+prDictStr(Type1)
	Type2 = getTypes(samples2)
	print 'Types1'+prDictStr(Type1)
	TypeA = getTypes(samplesA)
	print 'Types1'+prDictStr(Type1)
	
	[eucDistFib1,eucDisthESC1] = EUC(bData1,Type1)
	[eucDistFib2,eucDisthESC2] = EUC(bData2,Type2)
	[eucDistFibA,eucDisthESCA] = EUC(bDataA,TypeA)
	
	times1 = GetTimes(eucDistFib1,eucDisthESC1,bData1)
	times2 = GetTimes(eucDistFib2,eucDisthESC2,bData2)
	timesA = GetTimes(eucDistFibA,eucDisthESCA,bDataA)
	print '\n\n'
	prDict({'eucDistfib1':eucDistFib1,'eucDisthESC1':eucDisthESC1,'times1':times1})
	print '\n\n'
	prDict({'eucDistfib2':eucDistFib1,'eucDisthESC2':eucDisthESC1,'times2':times1})
	print '\n\n'
	prDict({'eucDistfibA':eucDistFib1,'eucDisthESCA':eucDisthESC1,'timesA':times1})

	if not os.path.exists('Fib'): os.makedirs('Fib')
	os.chdir('Fib')
	testFib1 = windowed(eucDistFib1,bData1,window,steps,step)
	pgHeatMap2(testFib1,genes1,'PolyFib')
	testFib2 = windowed(eucDistFib2,bData2,window,steps,step)
	pgHeatMap2(testFib2,genes2,'MonoFib')
	testFibA = windowed(eucDistFibA,bDataA,window,steps,step)
	pgHeatMap2(testFibA,genesA,'AggFib')
	os.chdir("..")
	if not os.path.exists('hESC'): os.makedirs('hESC')
	os.chdir('hESC')
	testESC1 = windowed(eucDisthESC1,bData1,window,steps,step)
	pgHeatMap2(testESC1,genes1,'PolyhESC')
	testESC2 = windowed(eucDisthESC2,bData2,window,steps,step)
	pgHeatMap2(testESC2,genes2,'MonohESC')
	testESCA = windowed(eucDisthESCA,bDataA,window,steps,step)
	pgHeatMap2(testESCA,genesA,'AgghESC')
	os.chdir("..")
	if not os.path.exists('Time'): os.makedirs('Time')
	os.chdir('Time')
	testTime1= windowedT(times1,bData1,window,steps,step)
	pgHeatMap2(testTime1,genes1,'PolyProgression')
	testTime2= windowedT(times2,bData2,window,steps,step)
	pgHeatMap2(testTime2,genes2,'MonoProgression')
	testTimeA= windowedT(timesA,bDataA,window,steps,step)
	pgHeatMap2(testTimeA,genesA,'AggProgression')
	os.chdir("..")
	htmlDoc='<<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n<html smlns="http://www.w3.org/1999/xhtml">\n'
	htmlDoc+='<head>\n<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n<title>MonoVPoly</title>\n<style type="text/css">\n<!--\nbody {\nfont: 100% Verdana, Arial, Helvetica, sans-serif;\nbackground: #666666;\n	margin: 0; /* its good practice to zero the margin and padding of the body element to account for differing browser defaults */\npadding: 0;\ntext-align: center; /* this centers the container in IE 5* browsers. The text is then set to the left aligned default in the #container selector */\ncolor: #000000;\n}\n.oneColLiqCtr #container {\nwidth: 500px;  /* this will create a container 80% of the browser width */\nbackground: #FFFFFF;\nmargin: 0 auto; /* the auto margins (in conjunction with a width) center the page */\nborder: 1px solid #000000;\ntext-align: left; /* this overrides the text-align: center on the body element. */\n}\n.oneColLiqCtr #mainContent {\npadding: 12px 20px 0 20px; /* remember that padding is the space inside the div box and margin is the space outside the div box */\n}\n-->\n</style>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryCollapsiblePanel.js" type="text/javascript"></script>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/Scripts/swfobject_modified.js" type="text/javascript"></script>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryMenuBar.js" type="text/javascript"></script>\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/CSS/Level3_3.css" rel="stylesheet" type="text/css">\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryCollapsiblePanel.css" rel="stylesheet" type="text/css">\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryMenuBarHorizontal_edited.css" rel="stylesheet" type="text/css">\n<style type="text/css" media="screen">#FlashID {visibility:hidden}</style></head>\n'
	genes={1:genes1,2:genes2,3:genesA}
	htmlDoc+='<body>'

	htmlDoc+='<div id="Time" class="CollapsiblePanel CollapsiblePanelClosed">\n'
	htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">Time</div>\n'
	htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
	htmlDoc+=htMap3(tkSortOneB,'tkSortOneB',testTime1,testTime2,testTimeA,genes,'Time')
	htmlDoc+='</div>'
	htmlDoc+='</div>'

	htmlDoc+='<div id="hESCDist" class="CollapsiblePanel CollapsiblePanelClosed">\n'
	htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">hESCDist</div>\n'
	htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
	htmlDoc+=htMap3(tkSortOneB,'tkSortOneB',testESC1,testESC2,testESCA,genes,'hESC')
	htmlDoc+='</div>'
	htmlDoc+='</div>'
	
	htmlDoc+='<div id="FibDist" class="CollapsiblePanel CollapsiblePanelClosed">\n'
	htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">FibDist</div>\n'
	htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
	htmlDoc+=htMap3(tkSortOneB,'tkSortOneB',testFib1,testFib2,testFibA,genes,'Fib')
	htmlDoc+='</div>'
	htmlDoc+='</div>'

	htmlDoc+='<script type="text/javascript">\n<!--\nvar CollapsiblePanel1 = new Spry.Widget.CollapsiblePanel("FibDist", {contentIsOpen:false});\nvar CollapsiblePanel2 = new Spry.Widget.CollapsiblePanel("hESCDist", {contentIsOpen:false});\nvar CollapsiblePanel3 = new Spry.Widget.CollapsiblePanel("Time", {contentIsOpen:false});\nswfobject.registerObject("FlashID");\nvar MenuBar1 = new Spry.Widget.MenuBar("MenuBar1", {imgDown:"SpryAssets/SpryMenuBarDownHover.gif", imgRight:"SpryAssets/SpryMenuBarRightHover.gif"});\n//-->\n</script>\n'
	htmlDoc+='</body>'
	htmlDoc+='</html>'

	output=open('Comp.html',"w")
	output.write(htmlDoc)
	output.close()
	pass
#	htMap(tkSortFit(testFib,fitsFib))
#	htMap(tkSortFit(testESC,fitsESC))
#	htMap(tkSortFit(testTime,fitsTime))

#HTML(data,r'./helloAlbert')
#ak = sys.argv[2]
#print data
#HTML(data,ak) 							#<r'filename'>#enter the name of the folder of pipeline site


#HeatMap(tkSortHalf(testFib))
#HeatMap(tkSortHalf(testESC))
#HeatMap(tkSortHalf(testTime))

#HeatMap(tkSortOne(testFib))
#HeatMap(tkSortOne(testESC))
#HeatMap(tkSortOne(testTime))

#HeatMap(tkSortOneB(testFib))
#HeatMap(tkSortOneB(testESC))
#HeatMap(tkSortOneB(testTime))


#HeatMap(tkSortDifs(testFib))
#HeatMap(tkSortDifs(testESC))
#HeatMap(tkSortDifs(testTime))

#HeatMap(tkSortFit(testFib,fitsFib))
#HeatMap(tkSortFit(testESC,fitsESC))
#HeatMap(tkSortFit(testTime,fitsTime))

#HeatMap(testFib)
#HeatMap(testESC)
#HeatMap(testTime)

#HeatMap(byPlace(testFib,findOrder(testFib)))
#HeatMap(byPlace(testESC,findOrder(testESC)))
#HeatMap(byPlace(testTime,findOrder(testTime)))

#NBD=numpy.zeros((len(BData),len(BData[0])-3))
##BData
##print NBD
#for k in range(1,len(BData),1):
#	row=[]
#	k1=1
#	while(k1<97):
#		if(is_number(BData[k][k1])):
#			row.append(float(BData[k][k1]))
#		k1+=1
#	NBD[k]=row

##print NBD

#pca=mdp.pca(NBD)

#pca1=numpy.zeros(len(NBD))
#pca2=numpy.zeros(len(NBD))
#pca3=numpy.zeros(len(NBD))

#for k in range(0,len(NBD),1):						#compute the projection ofeach cell onto Principal comps
#	mpca1=0					#Magnitude PCA
#	mpca2=0	
#	mpca3=0
#	dpd1=0
#	dpd2=0 
#	dpd3=0					#Dot product NBD and PCA
#	for k1 in range(0,len(NBD[1]),1):
#		dpd1+=(NBD[k][k1]*pca[0][k1])
#		mpca1+=pca[0][k1]*pca[0][k1]
#		dpd2+=(NBD[k][k1]*pca[1][k1])
#		mpca2+=pca[1][k1]*pca[1][k1]
#		dpd3+=(NBD[k][k1]*pca[2][k1])
#		mpca3+=pca[2][k1]*pca[2][k1]
#	pca1[k]=dpd1/(numpy.sqrt(mpca1))
#	pca2[k]=dpd2/(numpy.sqrt(mpca2))

#def PCAScat(pca1,pca2):
#	scat = plt.figure()
#	p1 = scat.add_subplot(111)

#	p1.scatter(pca1,pca2,color='blue', s=7,edgecolor='none')
#	p1.set_aspect(1./p1.get_data_ratio())
#	p1.grid(True)

#	plt.show()

#PCAScat(pca1,pca2)
#Scat(eucDisthESC,eucDistFib)
##for row in pca:
##	print row
##print pca.shape

#if(is_bool(sys.argv[3])):				###TEST MODE
#	s,g = ok.BCM(npClip(data))
#	S,G = ok.BCM(npClip(Boolienize(data)))
def getHescFib(File1,FileWithFib_Hesc):
	known = getTypesMP(getRowNames(File1))
	hasHESC=False
	hasfib=False
	#print known
	for k in known['hESC']:
		if known['hESC'][k] == '+':
			hasHESC=True
	for k in known['Fib']:
		if known['Fib'][k] == '+':
			hasfib=True
	#print 'hasHESC:',hasHESC
	#print 'hasfib:',hasfib
	n = 0
	#print '!!!!!!!!'
	if(hasHESC == False):
		known2 = getTypesMP(getRowNames(FileWithFib_Hesc))
		k2 = 0
		for k in FileWithFib_Hesc:
			try:
				if known2['ID'][k] == 'hESC' or known2['hESC'][k] == '+':
					File1[len(File1)]=FileWithFib_Hesc[k]
			except(TypeError):
				if known2['ID'][k2] == 'hESC' or known2['hESC'][k2] == '+':
					File1.append(k)
			k2+=1
	if(hasfib == False):
		k2 = 0
		known2 = getTypesMP(getRowNames(FileWithFib_Hesc))
		for k in FileWithFib_Hesc:
			try:
				if known2['ID'][k] == 'Fib' or known2['Fib'][k] == '+':
					File1[len(File1)]=FileWithFib_Hesc[k]
			except(TypeError):
				if known2['ID'][k2] == 'Fib' or known2['Fib'][k2] == '+':
					File1.append(k)
			k2+=1
	return(File1)
							##############################################################
def prDict2(Dists):					# This is usful for printing csv files.#######################
	tail=""						#Here we will us it to construct a predetermined input for R##
	#n=0						#With a columns 'Time',gene1,gene2,gene3...###################
	for k in Dists:					#################val###val###val###val########################
		for k2 in Dists[k]:			#################val...val...val...val########################
			tail+=str(k2)+','		##############################################################
			#n2+=1
		tail=tail[0:len(tail)-1]
		tail+="\n"
	tail=tail[0:len(tail)-1]
	#Head+="\n"
	print tail
	pass

def prDict2wd(Dists,delim):				
	tail=""					
	#n=0						
	for k in Dists:					
		for k2 in Dists[k]:			
			tail+=str(k2)+delim		#DEF WHITTLE
			#n2+=1
		tail=tail[0:len(tail)-1]
		tail+="\n"
	tail=tail[0:len(tail)-1]
	#Head+="\n"
	#print tail
	return(tail)

def whittle(File1,FileWithFib_Hesc,T):
	#print '!!!!!!!!'				#################################################################
	known = getTypesMP(getRowNames(File1))		#To construct the 'time' determined by cells expression profile #
	hasHESC=False					#we must define the begining and end point in our data set	#
	hasfib=False					#################################################################
	#print known
	for k in known['hESC']:
		if known['hESC'][k] == '+':
			hasHESC=True
	for k in known['Fib']:
		if known['Fib'][k] == '+':
			hasfib=True
	#print 'hasHESC:',hasHESC
	#print 'hasfib:',hasfib
	n = 0
	#print '!!!!!!!!'
	if(hasHESC == False):
		known2 = getTypesMP(mk.getRowNames(FileWithFib_Hesc))
		for k in FileWithFib_Hesc:
			if known2['ID'][k] == 'hESC' or known2['hESC'][k] == '+':
				File1[len(File1)]=FileWithFib_Hesc[k]
	if(hasfib == False):
		known2 = getTypesMP(getRowNames(FileWithFib_Hesc))
		for k in FileWithFib_Hesc:
			if known2['ID'][k] == 'Fib' or known2['Fib'][k] == '+':
				File1[len(File1)]=FileWithFib_Hesc[k]
	#print '!!!!!!!!'
	#print File1
	rFile1={}
	for k in File1:
		rFile1[k]=[]
		rFile1[k].append(File1[k][0])
		#print k,File1[k][0]
	m=k
	for k in range(0,len(File1[m])):
		c=0
		for k2 in File1:
			if(is_number(File1[k2][k])):
				if File1[k2][k] != 0:
					c+=1
		if c > T:
			n=0
			for k3 in rFile1:
#				try:
				rFile1[k3].append(File1[n][k])
				n+=1
#				except(KeyError):
#					rFile1[k3]=[]
#					rFile1[k3].append=File1[k3][k]
	first = True
	for k in rFile1:
		if first:
			l = len(rFile1[k])
		else:
			if(len(rFile1[k])!= l):
				print 'wtf!!!' 
				print k
	#PR.prDict2(rFile1)
	#print mk.getTypes(mk.getRowNames(rFile1))
	[a,b] = EUC(rFile1,getTypes(getRowNames(rFile1)))
	times = GetTimesNormalized(a,b,rFile1)
	n=0
	first=True
	for k in rFile1:
		if first:
			rFile1[k][0]='Time'
			first=False
			n+=1
		else:
			rFile1[k][0]=times[n]
			n+=1
	#PR.prDict2(rFile1)
	return(rFile1)

def HTML4(Path):
	Data={}
	files = [f for f in os.listdir('.') if os.path.isfile(f)]
	print '1.1 - ReadingFiles:'
	for f in files:
		if f[0]=='_':
			Data[f]=boolienize(clip(readFileC(f)))		# Mark files for input with an | as the first charactor in their name
			print f
			if f[1]=='_':
				Twep=f					# Mark the TableWithEndPoints with ||

	Agg={}				# Combined the 2 inputs into one table
	AggT={}
	n=0
	first = True
	First = True
	for k in Data:
		if first:
			first=False				#The first row or heading is only required in the first file
			for k2 in Data[k]:
				Agg[n]=Data[k][k2]
				n+=1
		else:
			First=True
			for k2 in Data[k]:
				if First:
					First=False
				else:	#Make an agragate table - It is important to note that this assumes Assays are alligned across input tables.
					Agg[n]=Data[k][k2]
					n+=1
	AggT = Agg # IDK what I'm even doing at this point
	print '1.2 - Combined Tables Made'
	for k in Data:
		if k != Twep:
			Data[k]=getHescFib(Data[k],Data[Twep])
	Data['Agg']=Agg

	steps = 200
	window=numpy.sqrt(48.0)/20
	step=numpy.sqrt(48.0)/steps
	ds={}				#Ds will contain a dictionary for each column associated to a row of data.

	ds['genes']={}
	ds['samples']={}
	ds['type']={}
	ds['euc']={}
	ds['time']={}
	ds['testFib']={}
	ds['testESC']={}
	ds['testTime']={}
	ds['Rinput']={}
	ds['name']={}
	ds['pca']={}
	print '1.3 - Makeing R input tables:'
	for k in Data:
		ds['genes'][k]=getColNames(Data[k])
		ds['samples'][k]=getRowNames(Data[k])
		ds['type'][k]=getTypes(ds['samples'][k])
		ds['euc'][k]={}
		ds['euc'][k]=EUC(Data[k],ds['type'][k])
		ds['time'][k]=GetTimes(ds['euc'][k][0],ds['euc'][k][1],Data[k])
		ds['testFib'][k]=windowed(ds['euc'][k][0],Data[k],window,steps,step)
		ds['testESC'][k]=windowed(ds['euc'][k][1],Data[k],window,steps,step)
		ds['testTime'][k]=windowedT(ds['time'][k],Data[k],window,steps,step)
		#ds['pca'][k]={}
		ds['pca'][k]=PCA(Data[k])
		ds['Rinput'][k]=whittle(Data[k],Data[Twep],10)
		ds['name'][k]=k[0:len(k)-4]
		output=open('R'+ds['name'][k]+'.csv',"w")
		output.write(prDict2wd(ds['Rinput'][k],','))
		output.close()
		print '\t'+k

	print '1.4 - Writing next.sh'
	output = open('next.sh',"w")
	for k in Data:
		output.write('Rscript fit_testing_ajay.R \'R'+ds['name'][k]+'.csv\' Rstuf/"'+ds['name'][k]+'FitCurves" "log'+ds['name'][k]+'" &\n')
	output.close()

	if not os.path.exists(Path): os.makedirs(Path)
	os.chdir(Path)
	if not os.path.exists(r'G'):
		os.makedirs(r'G')
	os.chdir(r'G')
	print '1.5 - Writing HTML for gene pages'
	for k3 in range(1,len(ds['genes'][k]),1):
		if not os.path.exists(ds['genes'][k][k3]): 
			os.makedirs(ds['genes'][k][k3])
	k3 = k
	for k in range(1,len(ds['genes'][k3]),1):
		print ds['genes'][k3][k]
		os.chdir(str(ds['genes'][k3][k]))
		output=open(ds['genes'][k3][k]+'.html',"w")
		htmlDoc='<<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n<html smlns="http://www.w3.org/1999/xhtml">\n'
		htmlDoc+='<head>\n<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n<title>'+ds['genes'][k3][k]+'</title>\n<style type="text/css">\n<!--\nbody {\nfont: 100% Verdana, Arial, Helvetica, sans-serif;\nbackground: #666666;\n	margin: 0; /* its good practice to zero the margin and padding of the body element to account for differing browser defaults */\npadding: 0;\ntext-align: center; /* this centers the container in IE 5* browsers. The text is then set to the left aligned default in the #container selector */\ncolor: #000000;\n}\n.oneColLiqCtr #container {\nwidth: 500px;  /* this will create a container 80% of the browser width */\nbackground: #FFFFFF;\nmargin: 0 auto; /* the auto margins (in conjunction with a width) center the page */\nborder: 1px solid #000000;\ntext-align: left; /* this overrides the text-align: center on the body element. */\n}\n.oneColLiqCtr #mainContent {\npadding: 12px 20px 0 20px; /* remember that padding is the space inside the div box and margin is the space outside the div box */\n}\n-->\n</style>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryCollapsiblePanel.js" type="text/javascript"></script>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/Scripts/swfobject_modified.js" type="text/javascript"></script>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryMenuBar.js" type="text/javascript"></script>\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/CSS/Level3_3.css" rel="stylesheet" type="text/css">\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryCollapsiblePanel.css" rel="stylesheet" type="text/css">\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryMenuBarHorizontal_edited.css" rel="stylesheet" type="text/css">\n<style type="text/css" media="screen">#FlashID {visibility:hidden}</style></head>\n'
		htmlDoc+='<body>'
		Eucplots=''
		PCAplots=''
		fitCurves=''
		fitCurvesScriptBit=''
		n=2;
		for k2 in Data:
			ScatA(ds['euc'][k2][1],ds['euc'][k2][0],Data[k2],ds['genes'][k2][k],'_EUC_'+k2[0:len(k2)-4],False,False,ds['name'][k2],2)
			ScatA(ds['pca'][k2][1],ds['pca'][k2][0],Data[k2],ds['genes'][k2][k],'_PCA_'+k2[0:len(k2)-4],True,True,ds['name'][k2],3)
			Eucplots+='<img src="'+ds['genes'][k2][k]+'_EUC_'+k2[0:len(k2)-4]+'.png" width="800" height="600">\n'
			PCAplots+='<img src="'+ds['genes'][k2][k]+'_PCA_'+k2[0:len(k2)-4]+'.png" width="800" height="600">\n'
			
			fitCurves+='<div id="FitCurves'+k2+'" class="CollapsiblePanel CollapsiblePanelClosed">\n'
			fitCurves+='<div class="CollapsiblePanelTab" tabindex="0">FitCurves'+k2+'</div>\n'
			fitCurves+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">\n'
			fitCurves+='<img src="../../Rstuf/'+ds['name'][k2]+'FitCurvesfigures/'+ds['genes'][k3][k]+'_LinearRegression.jpg" width="480" height="480">\n'
			fitCurves+='<img src="../../Rstuf/'+ds['name'][k2]+'FitCurvesfigures/'+ds['genes'][k3][k]+'_Normals-1.jpg" width="480" height="480">\n'
			fitCurves+='<img src="../../Rstuf/'+ds['name'][k2]+'FitCurvesfigures/'+ds['genes'][k3][k]+'_Normals-2.jpg" width="480" height="480">\n'
			fitCurves+='<img src="../../Rstuf/'+ds['name'][k2]+'FitCurvesfigures/'+ds['genes'][k3][k]+'_Normals-3.jpg" width="480" height="480">\n'
			fitCurves+='<img src="../../Rstuf/'+ds['name'][k2]+'FitCurvesfigures/'+ds['genes'][k3][k]+'_Uniform.jpg" width="480" height="480">\n'
			fitCurves+='</div>\n'
			fitCurves+='</div>\n'
			n+=1
			fitCurvesScriptBit+='\nvar CollapsiblePanel'+str(n)+' = new Spry.Widget.CollapsiblePanel("FitCurves'+k2+'", {contentIsOpen:false});'

		htmlDoc+='<div id="Eucledian Plots" class="CollapsiblePanel CollapsiblePanelClosed">\n'
		htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">Eucledian Plots</div>\n'
		htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
		htmlDoc+=Eucplots
		htmlDoc+='</div>'
		htmlDoc+='</div>'
		htmlDoc+='<div id="PCA Plots" class="CollapsiblePanel CollapsiblePanelClosed">\n'
		htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">PCA Plots</div>\n'
		htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
		htmlDoc+=PCAplots
		htmlDoc+='</div>'
		htmlDoc+='</div>'
		htmlDoc+=fitCurves
		htmlDoc+='</body>'
		#htmlDoc+='<body><img src="../../Violins/'+genes1[k]+'.png" width="800" height="600"></body>\n'
		htmlDoc+='<script type="text/javascript">\n<!--\nvar CollapsiblePanel1 = new Spry.Widget.CollapsiblePanel("Eucledian Plots", {contentIsOpen:false});\nvar CollapsiblePanel2 = new Spry.Widget.CollapsiblePanel("PCA Plots", {contentIsOpen:false});\n'+fitCurvesScriptBit+'\nswfobject.registerObject("FlashID");\nvar MenuBar1 = new Spry.Widget.MenuBar("MenuBar1", {imgDown:"SpryAssets/SpryMenuBarDownHover.gif", imgRight:"SpryAssets/SpryMenuBarRightHover.gif"});\n//-->\n</script>\n'
		htmlDoc+='<p>'+ds['genes'][k3][k]+'<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene='+ds['genes'][k3][k]+'>'+ds['genes'][k3][k]+'</a></p>\n'
		htmlDoc+='</body>'
		htmlDoc+='</html>'
		output.write(htmlDoc)
		output.close()
		os.chdir('..')
	print '1.6 - Writing top level HTML'
	if not os.path.exists(r'Fib'): os.makedirs(r'Fib')
	if not os.path.exists(r'hESC'): os.makedirs(r'hESC')
	if not os.path.exists(r'Time'): os.makedirs(r'Time')
	os.chdir(r'Fib')
	pgHeatMap(ds['testFib']['Agg'],ds['genes']['Agg'])
	os.chdir("..")
	os.chdir(r'hESC')
	pgHeatMap(ds['testESC']['Agg'],ds['genes']['Agg'])
	os.chdir("..")
	os.chdir(r'Time')
	pgHeatMap(ds['testTime']['Agg'],ds['genes']['Agg'])
	os.chdir("..")
	os.chdir("..")

#	output=open(f+ssort+'.html',"w")
#	output.write(htmlDoc)
#	output.close()
	htMap2(tkSortHalf,ds['testFib']['Agg'],'Fib','SortHalf',ds['genes'][k3])
	htMap2(tkSortHalf,ds['testESC']['Agg'],'hESC','SortHalf',ds['genes'][k3])
	htMap2(tkSortHalf,ds['testTime']['Agg'],'Time','SortHalf',ds['genes'][k3])

	htMap2(tkSortOne,ds['testFib']['Agg'],'Fib','SortOne',ds['genes'][k3])
	htMap2(tkSortOne,ds['testESC']['Agg'],'hESC','SortOne',ds['genes'][k3])
	htMap2(tkSortOne,ds['testTime']['Agg'],'Time','SortOne',ds['genes'][k3])

	htMap2(tkSortOneB,ds['testFib']['Agg'],'Fib','SortOneB',ds['genes'][k3])
	htMap2(tkSortOneB,ds['testESC']['Agg'],'hESC','SortOneB',ds['genes'][k3])
	htMap2(tkSortOneB,ds['testTime']['Agg'],'Time','SortOneB',ds['genes'][k3])

	htMap2(tkSortDifs,ds['testFib']['Agg'],'Fib','SortDifs',ds['genes'][k3])
	htMap2(tkSortDifs,ds['testESC']['Agg'],'hESC','SortDifs',ds['genes'][k3])
	htMap2(tkSortDifs,ds['testTime']['Agg'],'Time','SortDifs',ds['genes'][k3])
	pass

def HTML5(Path):
	Data={}
	MarK={}
	files = [f for f in os.listdir('.') if os.path.isfile(f)]
	print '1.1 - ReadingFiles:'
	for f in files:
		if f[0]=='_':
			Data[f]=boolienize(clip(readFileC(f)))		# Mark files for input with an | as the first charactor in their name
			print f
			if f[1]=='_':
				Twep=f					# Mark the TableWithEndPoints with ||
			if f[1]=='~':mark[f]=True
			else:mark[f]=False
				

	Agg={}				# Combined the 2 inputs into one table
	AggT={}
	n=0
	first = True
	First = True
	for k in Data:
		if first:
			first=False				#The first row or heading is only required in the first file
			for k2 in Data[k]:
				Agg[n]=Data[k][k2]
				n+=1
		else:
			First=True
			for k2 in Data[k]:
				if First:
					First=False
				else:	#Make an agragate table - It is important to note that this assumes Assays are alligned across input tables.
					Agg[n]=Data[k][k2]
					n+=1
	AggT = Agg # IDK what I'm even doing at this point
	print '1.2 - Combined Tables Made'
	for k in Data:
		if MarK[k]:
			continue
		elif k != Twep:
			Data[k]=getHescFib(Data[k],Data[Twep])
	Data['Agg']=Agg

	steps = 200
	window=numpy.sqrt(48.0)/20
	step=numpy.sqrt(48.0)/steps
	ds={}				#Ds will contain a dictionary for each column associated to a row of data.

	ds['genes']={}
	ds['samples']={}
	ds['type']={}
	ds['euc']={}
	ds['time']={}
	ds['testFib']={}
	ds['testESC']={}
	ds['testTime']={}
	ds['Rinput']={}
	ds['name']={}
	ds['pca']={}
	print '1.3 - Makeing R input tables:'
	for k in Data:
		ds['genes'][k]=getColNames(Data[k])
		ds['samples'][k]=getRowNames(Data[k])
		ds['type'][k]=getTypes(ds['samples'][k])
		ds['euc'][k]={}
		ds['euc'][k]=EUC(Data[k],ds['type'][k])
		ds['time'][k]=GetTimes(ds['euc'][k][0],ds['euc'][k][1],Data[k])
		ds['testFib'][k]=windowed(ds['euc'][k][0],Data[k],window,steps,step)
		ds['testESC'][k]=windowed(ds['euc'][k][1],Data[k],window,steps,step)
		ds['testTime'][k]=windowedT(ds['time'][k],Data[k],window,steps,step)
		#ds['pca'][k]={}
		ds['pca'][k]=PCA(Data[k])
		ds['Rinput'][k]=whittle(Data[k],Data[Twep],10)
		ds['name'][k]=k[0:len(k)-4]
		output=open('R'+ds['name'][k]+'.csv',"w")
		output.write(prDict2wd(ds['Rinput'][k],','))
		output.close()
		print '\t'+k

	print '1.4 - Writing next.sh'
	output = open('next.sh',"w")
	for k in Data:
		output.write('Rscript fit_testing_ajay.R \'R'+ds['name'][k]+'.csv\' Rstuf/"'+ds['name'][k]+'FitCurves" "log'+ds['name'][k]+'" &\n')
	output.close()

	if not os.path.exists(Path): os.makedirs(Path)
	os.chdir(Path)
	if not os.path.exists(r'G'):
		os.makedirs(r'G')
	os.chdir(r'G')
	print '1.5 - Writing HTML for gene pages'
	for k3 in range(1,len(ds['genes'][k]),1):
		if not os.path.exists(ds['genes'][k][k3]): 
			os.makedirs(ds['genes'][k][k3])
	k3 = k
	for k in range(1,len(ds['genes'][k3]),1):
		print ds['genes'][k3][k]
		os.chdir(str(ds['genes'][k3][k]))
		output=open(ds['genes'][k3][k]+'.html',"w")
		htmlDoc='<<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n<html smlns="http://www.w3.org/1999/xhtml">\n'
		htmlDoc+='<head>\n<meta http-equiv="Content-Type" content="text/html; charset=utf-8">\n<title>'+ds['genes'][k3][k]+'</title>\n<style type="text/css">\n<!--\nbody {\nfont: 100% Verdana, Arial, Helvetica, sans-serif;\nbackground: #666666;\n	margin: 0; /* its good practice to zero the margin and padding of the body element to account for differing browser defaults */\npadding: 0;\ntext-align: center; /* this centers the container in IE 5* browsers. The text is then set to the left aligned default in the #container selector */\ncolor: #000000;\n}\n.oneColLiqCtr #container {\nwidth: 500px;  /* this will create a container 80% of the browser width */\nbackground: #FFFFFF;\nmargin: 0 auto; /* the auto margins (in conjunction with a width) center the page */\nborder: 1px solid #000000;\ntext-align: left; /* this overrides the text-align: center on the body element. */\n}\n.oneColLiqCtr #mainContent {\npadding: 12px 20px 0 20px; /* remember that padding is the space inside the div box and margin is the space outside the div box */\n}\n-->\n</style>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryCollapsiblePanel.js" type="text/javascript"></script>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/Scripts/swfobject_modified.js" type="text/javascript"></script>\n<script src="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryMenuBar.js" type="text/javascript"></script>\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/CSS/Level3_3.css" rel="stylesheet" type="text/css">\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryCollapsiblePanel.css" rel="stylesheet" type="text/css">\n<link href="http://lrrpublic.cli.det.nsw.edu.au/lrrSecure/Sites/LRRView/10378/applets/spry/SpryAssets/SpryMenuBarHorizontal_edited.css" rel="stylesheet" type="text/css">\n<style type="text/css" media="screen">#FlashID {visibility:hidden}</style></head>\n'
		htmlDoc+='<body>'
		Eucplots=''
		PCAplots=''
		fitCurves=''
		fitCurvesScriptBit=''
		n=2
		for k2 in Data:
			ScatA(ds['euc'][k2][1],ds['euc'][k2][0],Data[k2],ds['genes'][k2][k],'_EUC_'+k2[0:len(k2)-4],False,False,ds['name'][k2],2)
			ScatA(ds['pca'][k2][1],ds['pca'][k2][0],Data[k2],ds['genes'][k2][k],'_PCA_'+k2[0:len(k2)-4],True,True,ds['name'][k2],3)
			Eucplots+='<img src="'+ds['genes'][k2][k]+'_EUC_'+k2[0:len(k2)-4]+'.png" width="800" height="600">\n'
			PCAplots+='<img src="'+ds['genes'][k2][k]+'_PCA_'+k2[0:len(k2)-4]+'.png" width="800" height="600">\n'
			
			fitCurves+='<div id="FitCurves'+k2+'" class="CollapsiblePanel CollapsiblePanelClosed">\n'
			fitCurves+='<div class="CollapsiblePanelTab" tabindex="0">FitCurves'+k2+'</div>\n'
			fitCurves+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">\n'
			fitCurves+='<img src="../../Rstuf/'+ds['name'][k2]+'FitCurvesfigures/'+ds['genes'][k3][k]+'_LinearRegression.jpg" width="480" height="480">\n'
			fitCurves+='<img src="../../Rstuf/'+ds['name'][k2]+'FitCurvesfigures/'+ds['genes'][k3][k]+'_Normals-1.jpg" width="480" height="480">\n'
			fitCurves+='<img src="../../Rstuf/'+ds['name'][k2]+'FitCurvesfigures/'+ds['genes'][k3][k]+'_Normals-2.jpg" width="480" height="480">\n'
			fitCurves+='<img src="../../Rstuf/'+ds['name'][k2]+'FitCurvesfigures/'+ds['genes'][k3][k]+'_Normals-3.jpg" width="480" height="480">\n'
			fitCurves+='<img src="../../Rstuf/'+ds['name'][k2]+'FitCurvesfigures/'+ds['genes'][k3][k]+'_Uniform.jpg" width="480" height="480">\n'
			fitCurves+='</div>\n'
			fitCurves+='</div>\n'
			n+=1
			fitCurvesScriptBit+='\nvar CollapsiblePanel'+str(n)+' = new Spry.Widget.CollapsiblePanel("'+k2+'", {contentIsOpen:false});'

		htmlDoc+='<div id="Eucledian Plots" class="CollapsiblePanel CollapsiblePanelClosed">\n'
		htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">Eucledian Plots</div>\n'
		htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
		htmlDoc+=Eucplots
		htmlDoc+='</div>'
		htmlDoc+='</div>'
		htmlDoc+='<div id="PCA Plots" class="CollapsiblePanel CollapsiblePanelClosed">\n'
		htmlDoc+='<div class="CollapsiblePanelTab" tabindex="0">PCA Plots</div>\n'
		htmlDoc+='<div class="CollapsiblePanelContent" style="display: none; visibility: visible; height: 0px;">'
		htmlDoc+=PCAplots
		htmlDoc+='</div>'
		htmlDoc+='</div>'
		htmlDoc+=fitCurves
		htmlDoc+='</body>'
		#htmlDoc+='<body><img src="../../Violins/'+genes1[k]+'.png" width="800" height="600"></body>\n'
		htmlDoc+='<script type="text/javascript">\n<!--\nvar CollapsiblePanel1 = new Spry.Widget.CollapsiblePanel("Eucledian Plots", {contentIsOpen:false});\nvar CollapsiblePanel2 = new Spry.Widget.CollapsiblePanel("PCA Plots", {contentIsOpen:false});'+fitCurvesScriptBit+'\nswfobject.registerObject("FlashID");\nvar MenuBar1 = new Spry.Widget.MenuBar("MenuBar1", {imgDown:"SpryAssets/SpryMenuBarDownHover.gif", imgRight:"SpryAssets/SpryMenuBarRightHover.gif"});\n//-->\n</script>\n'
		htmlDoc+='<p>'+ds['genes'][k3][k]+'<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene='+ds['genes'][k3][k]+'>'+ds['genes'][k3][k]+'</a></p>\n'
		htmlDoc+='</body>'
		htmlDoc+='</html>'
		output.write(htmlDoc)
		output.close()
		os.chdir('..')
	print '1.6 - Writing top level HTML'
	if not os.path.exists(r'Fib'): os.makedirs(r'Fib')
	if not os.path.exists(r'hESC'): os.makedirs(r'hESC')
	if not os.path.exists(r'Time'): os.makedirs(r'Time')
	os.chdir(r'Fib')
	pgHeatMap(ds['testFib']['Agg'],ds['genes']['Agg'])
	os.chdir("..")
	os.chdir(r'hESC')
	pgHeatMap(ds['testESC']['Agg'],ds['genes']['Agg'])
	os.chdir("..")
	os.chdir(r'Time')
	pgHeatMap(ds['testTime']['Agg'],ds['genes']['Agg'])
	os.chdir("..")
	os.chdir("..")

#	output=open(f+ssort+'.html',"w")
#	output.write(htmlDoc)
#	output.close()
	htMap2(tkSortHalf,ds['testFib']['Agg'],'Fib','SortHalf',ds['genes'][k3])
	htMap2(tkSortHalf,ds['testESC']['Agg'],'hESC','SortHalf',ds['genes'][k3])
	htMap2(tkSortHalf,ds['testTime']['Agg'],'Time','SortHalf',ds['genes'][k3])

	htMap2(tkSortOne,ds['testFib']['Agg'],'Fib','SortOne',ds['genes'][k3])
	htMap2(tkSortOne,ds['testESC']['Agg'],'hESC','SortOne',ds['genes'][k3])
	htMap2(tkSortOne,ds['testTime']['Agg'],'Time','SortOne',ds['genes'][k3])

	htMap2(tkSortOneB,ds['testFib']['Agg'],'Fib','SortOneB',ds['genes'][k3])
	htMap2(tkSortOneB,ds['testESC']['Agg'],'hESC','SortOneB',ds['genes'][k3])
	htMap2(tkSortOneB,ds['testTime']['Agg'],'Time','SortOneB',ds['genes'][k3])

	htMap2(tkSortDifs,ds['testFib']['Agg'],'Fib','SortDifs',ds['genes'][k3])
	htMap2(tkSortDifs,ds['testESC']['Agg'],'hESC','SortDifs',ds['genes'][k3])
	htMap2(tkSortDifs,ds['testTime']['Agg'],'Time','SortDifs',ds['genes'][k3])
	pass

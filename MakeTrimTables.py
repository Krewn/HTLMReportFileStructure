import main as mk
import os

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

Agg={}				# Combined the 2 inputs into one table
AggT={}
n=0
first = True
First = True
for k in Data:
	temp=[]
	if first:
		first=False				#The first row or heading is only required in the first file
		for k2 in Data[k]:
			temp.append(Data[k][k2])
			n+=1
	else:
		First=True
		for k2 in Data[k]:
			if First:
				First=False
			else:	#Make an agragate table - It is important to note that this assumes Assays are alligned across input tables.
				temp.append(Data[k][k2])
				n+=1
	Data[k]=temp
AggT = Agg # IDK what I'm even doing at this point
Data['Agg']=Agg
#print Data.keys()
euc = {}
PCA = {}
types = {}
Days = {}
for k in Data:
	for k2 in Data[k]:
		if(len(k2)!=len(Data[k][0])):
			print '\n'
			print k2
			print len(Data[k][k2])
			print len(Data[k][0])
			print '\n'
	print k
	print Twep
	Data[k]=mk.getHescFib(Data[k],Data[Twep])
	types[k]=mk.getTypes(mk.getRowNames(Data[k]))
	#print Data[k].keys()
	#print len(types[k])
	#print len(Data[k])
	euc[k]=mk.EUC(Data[k],types[k])
	PCA[k]=mk.PCA(Data[k])
	Days[k]=mk.getdays(mk.getRowNames(Data[k]))
	first = True
	k3 = 0
	for k2 in Data[k]:
		if first:
			k2.append('Type')
			k2.append('EucDistFib')
			k2.append('EucDistHesc')
			k2.append('Pca1')
			k2.append('Pca2')
			k2.append('Day')
			first=False
		else:
			k2.append(types[k][k3])
			k2.append(euc[k][0][k3-1])
			k2.append(euc[k][1][k3-1])
			k2.append(PCA[k][0][k3-1])
			k2.append(PCA[k][1][k3-1])
			k2.append(Days[k][k3])
		k3+=1
for k in Data:
	output=open('ajay'+k[0:len(k)-4],"w")
	output.write(mk.prDict2wd(Data[k],','))
	output.close()

























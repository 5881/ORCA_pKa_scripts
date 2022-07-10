#!/bin/python
#Скрипт рассчёта pKa из данных о энергии молекул рассчитанных 
#квантовохимическими методами.
#На вход принимает файл csv имеющий следующий формат:
#Имя соединения, dG рассчитанное из DFT, dG-E(el) из DFT, E(el) DLPNO-CCSD(T)
#В файле по строкам чередуются данные для анионов и нейтральных молекул
#Первые две строки занимает DMSO (или другой стандарт) как референс

import numpy as np
import sys
import csv
HEADERS=1
T=298
R=8.314
PKAREF=35 #pKa стандарта
KAREF=pow(10,-PKAREF)
def read_inputfile(filename):
	with open(filename, 'r') as datafile:
		headers=datafile.readline().split(',')
	tabl_format="U10"+",f4"*(len(headers)-1)
	#print(tabl_format)
	#dt=np.dtype("U10,f4,f4,f4,f4,f4")
	dt=np.dtype(tabl_format)
	#print(dt)
	data=np.genfromtxt(filename, delimiter=',',dtype=dt, skip_header=HEADERS)
	#print(len(data))
	with open(filename, 'r') as datafile:
		headers=datafile.readline().split(',')
	return data,headers

def calc_pka_dft_method(data,method,headers):
	#DMSO- + acid = DMSO + anion
	#G=[G-Eel](DFT)+Eel(DLPNO-CCSD(T))
	#dG(A>A-)=G(A-)-G(A)
	#dGpKa=dG(X>X-)-dG(ref>ref-)
	#ka=KAREF*np.exp(-ddg*2.6e6/(R*T))
	#pka=-np.log10(ka)
	#print(">>>>>>",method,headers[method+4])
	G=[]
	for i in data:
		G+=[i[3]+i[method+4]]
	#print("G",G)	
	dg=[]
	i=0
	while(i<len(data)):
		dg+=[G[i]-G[i+1]]
		i+=2
	#print("dg",dg)
	ddg=[]
	i=0
	while(i<len(dg)):
		ddg+=[dg[i]-dg[0]]
		i+=1
	#print("ddg",ddg)
	ddg=np.array(ddg)
	ddg=np.float128(ddg)
	#print(ddg)
	ka=KAREF*np.exp(-ddg*2.6e6/(R*T))
	#print(ka)
	pka=-np.log10(ka)
	#print(pka)
	#*********REPORT*************
	print("\npKa from DFT/",headers[method+4],sep='')
	print("compound","Ka","pKa")
	i=0
	while(i<len(pka)):
		print(data[2*i+1][0],ka[i],pka[i])
		i+=1
	#****************************
inputfile=sys.argv[1]
data, headers=read_inputfile(inputfile)
print(headers)
for i in range(len(headers)-4):
	calc_pka_dft_method(data,i,headers)


#calc_pka_dft(data)
#calc_pka_dft_ccsdt(data)
#print(data)

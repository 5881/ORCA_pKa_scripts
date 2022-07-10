#!/bin/python
#Скрипт для оптимизации структуры в два захода 
import os,sys
ORCA=os.path.expanduser('~')+"/orca/orca"
REPORTFILE="opt_report_file.csv"
SOLVENT="DMSO"
MAXCORE="11000"
NPROC="22"

def output_singlepoint_analyse(outputfile):
	E=0
	with open(outputfile, 'r') as out_file:
		for s in out_file:
			if 'FINAL SINGLE POINT ENERGY' in s: E=float(s.split()[-1])
	print(outputfile, "E=", E)
	#print("eee",E+G_el)
	return E

def output_freq_analyse(outputfile):
	gibbs_done=0 
	vib_freq_done=0
	geom_ok=1
	G=0
	G_el=0
	E=0
	with open(outputfile, 'r') as out_file:
		for s in out_file:
			if 'FINAL SINGLE POINT ENERGY' in s: E=float(s.split()[-1])
			if 'VIBRATIONAL FREQUENCIES' in s: vib_freq_done=1
			if vib_freq_done and '***imaginary mode***' in s:
				print("Мнимые частоты, геометрия возможно неверна")
				geom_ok=0
			if 'GIBBS FREE ENERGY' in s : 
				gibbs_done=1
				vib_freq_done=0
			if gibbs_done and 'Final Gibbs free energy' in s:
				G=float(s.split()[-2])
			if gibbs_done and 'G-E(el)' in s:
				G_el=float(s.split()[-4])
	print(outputfile, "E=", E, "G=", G,"G-el=",G_el, "OK=",geom_ok)
	#print("eee",E+G_el)
	return E, G, G_el, geom_ok

def getxyz_from_file(filename):
	xyz=[]
	with open(filename, 'r', encoding="utf-8", errors="ignore") as xyz_file:
		for line in xyz_file.readlines():
			if (len(line.split()) == 4 and line.split()[-1][-1].isdigit()):
				xyz.append(line[:-1])
	return xyz

def copyfile(src, dst):
	cmd="cp "+src+" "+dst
	os.system(cmd)
def movefile(src, dst):
	cmd="mv "+src+" "+dst
	os.system(cmd)


def gfnxtb_opt(xyzfile):
	job="_gfnxtb_opt"
	basename=xyzfile.replace('.xyz', '').rsplit('_')[0]
	basename=basename+job
	dirname=basename
	inpfile=basename+".inp"
	outputfile=basename+".out"
	xyzoutputfile=basename+".xyz"
	charge=0
	xyz=getxyz_from_file(xyzfile)
	if(basename.find("-")!=-1): charge=-1
	if(basename.find("+")!=-1): charge=1
	if(os.path.isdir(dirname) is False):
		os.mkdir(dirname)
	os.chdir(dirname)
	if(os.path.isfile(outputfile) is False):
		inputfilecontent="""
%pal nprocs NPROC end
%maxcore MAXCORE
!gfn-xtb opt xyzfile
%geom MaxIter 500 end
"""
		inputfilecontent=inputfilecontent.replace('NPROC',NPROC)
		inputfilecontent=inputfilecontent.replace('MAXCORE',MAXCORE)
		with open(inpfile, 'w') as inp_file:
			print(inputfilecontent,file=inp_file)
			print("* xyz", charge,"1",file=inp_file)
			for s in xyz:
				print(s,file=inp_file)
			print("*",file=inp_file)
		orcacmd=ORCA+" "+inpfile+"|tee "+outputfile
		os.system(orcacmd)
		os.chdir("..")
		outputfile=dirname+"/"+outputfile
	else:
		os.chdir("..")
	xyzfile=dirname+"/"+xyzoutputfile
	return xyzfile

def pbeh3c_opt(xyzfile):
	job="_pbeh3c_opt"
	basename=xyzfile.replace('.xyz', '').rsplit('_')[0]
	basename=basename+job
	dirname=basename
	inpfile=basename+".inp"
	outputfile=basename+".out"
	xyzoutputfile=basename+".xyz"
	charge=0
	xyz=getxyz_from_file(xyzfile)
	if(basename.find("-")!=-1): charge=-1
	if(basename.find("+")!=-1): charge=1
	if(os.path.isdir(dirname) is False):
		os.mkdir(dirname)
	os.chdir(dirname)
	if(os.path.isfile(outputfile) is False):
		inputfilecontent="""
%pal nprocs NPROC end
%maxcore MAXCORE
!pbeh-3c opt cpcm(SOLVENT) xyzfile
%geom MaxIter 500 end
"""
		inputfilecontent=inputfilecontent.replace('NPROC',NPROC)
		inputfilecontent=inputfilecontent.replace('MAXCORE',MAXCORE)
		inputfilecontent=inputfilecontent.replace('SOLVENT',SOLVENT)
		with open(inpfile, 'w') as inp_file:
			print(inputfilecontent,file=inp_file)
			print("* xyz", charge,"1",file=inp_file)
			for s in xyz:
				print(s,file=inp_file)
			print("*",file=inp_file)
		orcacmd=ORCA+" "+inpfile+"|tee "+outputfile
		os.system(orcacmd)
		os.chdir("..")
		#outputfile=dirname+"/"+outputfile
	else:
		os.chdir("..")
	xyzoutputfile=dirname+"/"+xyzoutputfile
	outputfile=dirname+"/"+outputfile
	return xyzoutputfile, outputfile


#создаём файл отчёта если его нет, вносим данные для DMSO чтобы сэкономить время
if(os.path.isfile(REPORTFILE) is False):
	with open(REPORTFILE, 'w') as reportfile:
		print("xyzfile, E",file=reportfile)


for xyzinputfile in sys.argv[1:]:
	xyzoutputfile=xyzinputfile
	xyzinputfile=gfnxtb_opt(xyzinputfile)
	optxyzfile, outputfile=pbeh3c_opt(xyzinputfile)
	E=output_singlepoint_analyse(outputfile)
	with open(REPORTFILE, 'a') as reportfile:
		print(optxyzfile.rsplit('/')[-1],E,sep=', ',file=reportfile)
	copyfile(optxyzfile, xyzoutputfile)

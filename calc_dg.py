#!/bin/python
#Скрипт для рассчёта dG образования веществ в растворе
#Скрипт принимает на вход файлы xyz, проводит предварительную оптимизацию
#PBEH-3c затем B3LYP(D3BJ)/def2-tzvp, считает гессиан
#затем считается single point для полученой структуры методами
#RI-MP2/ma-def2-tzvp и dlpno-ccsd(t)/ma-def2-tzvp
#Полученые данные заносятся в файл
#перед использованием следует задать правильные значения переменных
#ORCA,SOLVENT,MAXCORE,NPROC,MAXCORE_DLPNOCCSDT,NPROC_DLPNOCCSDT
#Скрипт проверен с ORCA 5.0.3
import os,sys
ORCA=os.path.expanduser('~')+"/orca/orca"
REPORTFILE="energy_report_file.csv"
SOLVENT="DMSO"
MAXCORE="11000"
NPROC="22"
CONTINUE=True
#Для DLPNO выставляется отдельно тк метод капризный
MAXCORE_DLPNOCCSDT="11000"
NPROC_DLPNOCCSDT="22"

def output_terminate_status(outputfile):
	with open(outputfile, 'r') as out_file:
		for s in out_file:
			if '****ORCA TERMINATED NORMALLY****' in s:
				return True
	return False

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
	
def pbeh3c_opt(xyzfile):
	RUN=True
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
	
	if(os.path.isfile(outputfile) is True):
		if(output_terminate_status(outputfile) is False):
			if(CONTINUE is True): xyz=getxyz_from_file(xyzoutputfile)
			#cmd='find . -not -name "'+xyzoutputfile+'" -delete'
			os.system("rm *")
		else: RUN=False
	
	if(RUN is True):
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
	outputfile=dirname+"/"+outputfile
	xyzfile=dirname+"/"+xyzoutputfile
	return outputfile,xyzfile

def b3lyp_d3bj_finalopt_freq_solv(xyzfile):
	RUN=True
	job="_b3lyp_d3bj_tzvp_optfreq"
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
	if(os.path.isfile(outputfile) is True):
		if(output_terminate_status(outputfile) is False):
			if(CONTINUE is True): xyz=getxyz_from_file(xyzoutputfile)
			#cmd='find . -not -name "'+xyzoutputfile+'" -delete'
			os.system("rm *")
		else: RUN=False
	
	if(RUN is True):
		inputfilecontent="""
%pal nprocs NPROC end
%maxcore MAXCORE
!b3lyp d3bj def2-tzvp def2-tzvp/c def2/j rijcosx opt freq tightscf tightopt cpcm(SOLVENT) xyzfile
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
		#os.chdir("..")
	
	
	os.chdir("..")
	xyzfile=dirname+"/"+xyzoutputfile
	outputfile=dirname+"/"+outputfile
	return outputfile,xyzfile

def ri_mp2_singlepoint_solv(xyzfile):
	RUN=True
	job="_rimp2_matzvp_solv"
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
	
	if(os.path.isfile(outputfile) is True):
		if(output_terminate_status(outputfile) is False): 
			os.system("rm *")
		else: RUN=False
	
	if(RUN is True):
		inputfilecontent="""
%pal nprocs NPROC end
%maxcore MAXCORE
!RI-MP2 ma-def2-tzvp autoaux tightscf nofrozencore cpcm(SOLVENT)
%basis
PCDTrimBas Overlap # Trim the orbital basis in the overlap metric
PCDTrimAuxJ Coulomb # Trim the AuxJ basis in the Coulomb metric
PCDTrimAuxJK Coulomb # Trim the AuxJK basis in the Coulomb metric
PCDTrimAuxC Coulomb # Trim the AuxC basis in the Coulomb metric
PCDThresh -1 # Threshold for the PCD: chosen automatically if <0
end
%scf
Thresh 1e-12
Sthresh 1e-7
end
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
	#xyzfile=dirname+"/"+xyzoutputfile
	outputfile=dirname+"/"+outputfile
	return outputfile

def dlpno_ccsdt_singlepoint_solv(xyzfile):
	RUN=True
	job="_dlpnoccsdt_matzvp_solv"
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
	if(os.path.isfile(outputfile) is True):
		if(output_terminate_status(outputfile) is False): 
			os.system("rm *")
		else: RUN=False
	
	if(RUN is True):
		inputfilecontent="""
%pal nprocs NPROC end
%maxcore MAXCORE
!dlpno-ccsd(t) ma-def2-tzvp autoaux tightscf cpcm(SOLVENT)
%basis
PCDTrimBas Overlap # Trim the orbital basis in the overlap metric
PCDTrimAuxJ Coulomb # Trim the AuxJ basis in the Coulomb metric
PCDTrimAuxJK Coulomb # Trim the AuxJK basis in the Coulomb metric
PCDTrimAuxC Coulomb # Trim the AuxC basis in the Coulomb metric
PCDThresh -1 # Threshold for the PCD: chosen automatically if <0
end
%scf
Thresh 1e-12
Sthresh 1e-7
end
"""
		inputfilecontent=inputfilecontent.replace('NPROC',NPROC_DLPNOCCSDT)
		inputfilecontent=inputfilecontent.replace('MAXCORE',MAXCORE_DLPNOCCSDT)
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
	#xyzfile=dirname+"/"+xyzoutputfile
	outputfile=dirname+"/"+outputfile
	return outputfile


#создаём файл отчёта если его нет, вносим данные для DMSO чтобы сэкономить время
if(os.path.isfile(REPORTFILE) is False):
	with open(REPORTFILE, 'w') as reportfile:
				print(
"""xyzfile,geom OK, G, G-Eel, E_RI-MP2, E_DLPNO-CCSD(T)
dmso-_b3lyp_d3bj_tzvp_optfreq.xyz, 1, -552.60891359, 0.03723523, -552.067138692847, -551.97374571165
dmso_b3lyp_d3bj_tzvp_optfreq.xyz, 1, -553.12829837, 0.0508621, -552.59076969229, -552.502786965165"""
,file=reportfile)


for xyzinputfile in sys.argv[1:]:
	stepoutputfile,stepxyzfile=pbeh3c_opt(xyzinputfile)
	if(output_terminate_status(stepoutputfile) is False): continue
	stepoutputfile,stepxyzfile=b3lyp_d3bj_finalopt_freq_solv(stepxyzfile)
	if(output_terminate_status(stepoutputfile) is False): continue
	E, G, G_el, geom_ok=output_freq_analyse(stepoutputfile)
	stepoutputfile=ri_mp2_singlepoint_solv(stepxyzfile)
	EMP2=output_singlepoint_analyse(stepoutputfile)
	stepoutputfile=dlpno_ccsdt_singlepoint_solv(stepxyzfile)
	EDLPNO_CCSDT=output_singlepoint_analyse(stepoutputfile)
	
	with open(REPORTFILE, 'a') as reportfile:
				print(stepxyzfile.rsplit('/')[-1],geom_ok,G,G_el,EMP2,EDLPNO_CCSDT,sep=', ',file=reportfile)

#!/usr/bin/env python3

import sys,re,gzip

#function to find abundance of contig
def findAbundance(contigName):
	abundance=0
	with gzip.open(sys.argv[1],'rt') as infile:
		for abundanceProfile in infile:
			if contigName in abundanceProfile:
				abundance=float(re.search(r'\s+([-\w.]+)\n', abundanceProfile).group(1))
				return(abundance)


#Defining each defense system type as a variable
(Abi2,AbiEii,AbiH,Aditi,AVAST,Azaca,
Borvo,BREX,BstA,Bunzi,
Cas,Cas_Class1,Cas_Class2,CBASS,
DarTG,Dazbog,dCTPdeaminase,dGTPase,DISARM,Dnd,Dodola,DRT,Druantia,Dsr,Dynamins,
Gabija,Gao_Ape,Gao_Her,Gao_Hhe,Gao_Iet,Gao_Mza,Gao_Ppl,Gao_Qat,Gao_RL,Gao_TerY,Gao_Tmn,Gao_Upx,
Hachiman,Kiwa,Lamassu,LamassuFam,Lit,Menshen,Mokosh,
Nhi,NixI,Olokun,Pif,PrrC,PsyrTA,
RADAR,Retron,RexAB,RM,RosmerTA,
Rst_2TM_1TM_TIR,Rst_3HP,Rst_DprAPPRT,Rst_DUF4238,Rst_gop_beta_cll,Rst_HelicaseDUF2290,
Rst_Hydrolase3Tm,Rst_Old_Tin,Rst_PARIS,Rst_RTnitrilaseTm,Rst_TIRNLR,
SEFIR,Septu,Shango,Shedu,ShosTA,SoFIC,SspBCDE,Stk2,
Thoeris,Tiamat,Uzume,Viperin,Wadjet,Zorya)=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)


#Function for looking for each type
for line in sys.stdin:
  #A - 6 stk
	if re.search(r'Abi2\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Abi2+=foundAbundance
	elif re.search(r'AbiEii\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		AbiEii+=foundAbundance
	elif re.search(r'AbiH\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		AbiH+=foundAbundance
	elif re.search(r'Aditi\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Aditi+=foundAbundance
	elif re.search(r'AVAST\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		AVAST+=foundAbundance
	elif re.search(r'Azaca\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Azaca+=foundAbundance
  #B - 4 stk
	elif re.search(r'Borvo\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Borvo+=foundAbundance
	elif re.search(r'BREX\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		BREX+=foundAbundance
	elif re.search(r'BstA\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		BstA+=foundAbundance
	elif re.search(r'Bunzi\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Bunzi+=foundAbundance
  #C - 4 stk
	elif re.search(r'Cas\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Cas+=foundAbundance
	elif re.search(r'Cas_Class1\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Cas_Class1+=foundAbundance	
	elif re.search(r'Cas_Class2\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Cas_Class2+=foundAbundance
	elif re.search(r'CBASS\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		CBASS+=foundAbundance   
  #D - 11 stk
	elif re.search(r'DarTG\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		DarTG+=foundAbundance 
	elif re.search(r'Dazbog\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Dazbog+=foundAbundance  
	elif re.search(r'dCTPdeaminase\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		dCTPdeaminase+=foundAbundance
	elif re.search(r'dGTPase\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		dGTPase+=foundAbundance    
	elif re.search(r'DISARM\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		DISARM+=foundAbundance
	elif re.search(r'Dnd\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Dnd+=foundAbundance
	elif re.search(r'Dodola\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Dodola+=foundAbundance    
	elif re.search(r'DRT\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		DRT+=foundAbundance    
	elif re.search(r'Druantia\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Druantia+=foundAbundance    
	elif re.search(r'Dsr\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Dsr+=foundAbundance      
	elif re.search(r'Dynamins\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Dynamins+=foundAbundance 
  #G - 12 stk     
	elif re.search(r'Gabija\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gabija+=foundAbundance  
	elif re.search(r'Gao_Ape\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_Ape+=foundAbundance
	elif re.search(r'Gao_Her\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_Her+=foundAbundance    
	elif re.search(r'Gao_Hhe\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_Hhe+=foundAbundance    
	elif re.search(r'Gao_Iet\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_Iet+=foundAbundance    
	elif re.search(r'Gao_Mza\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_Mza+=foundAbundance    
	elif re.search(r'Gao_Ppl\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_Ppl+=foundAbundance
	elif re.search(r'Gao_Qat\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_Qat+=foundAbundance
	elif re.search(r'Gao_RL\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_RL+=foundAbundance
	elif re.search(r'Gao_TerY\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_TerY+=foundAbundance
	elif re.search(r'Gao_Tmn\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_Tmn+=foundAbundance
	elif re.search(r'Gao_Upx\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Gao_Upx+=foundAbundance
	#H,K,L,M  - 7 stk 
	elif re.search(r'Hachiman\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Hachiman+=foundAbundance  
	elif re.search(r'Kiwa\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Kiwa+=foundAbundance  
	elif re.search(r'Lamassu\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Lamassu+=foundAbundance  
	elif re.search(r'Lamassu\-Fam\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		LamassuFam+=foundAbundance  
	elif re.search(r'Lit\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Lit+=foundAbundance
	elif re.search(r'Menshen\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Menshen+=foundAbundance
	elif re.search(r'Mokosh\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Mokosh+=foundAbundance
	#N,O,P - 6 stk             
	elif re.search(r'Nhi\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Nhi+=foundAbundance
	elif re.search(r'NixI\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		NixI+=foundAbundance
	elif re.search(r'Olokun\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Olokun+=foundAbundance    
	elif re.search(r'Pif\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Pif+=foundAbundance
	elif re.search(r'PrrC\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		PrrC+=foundAbundance
	elif re.search(r'PsyrTA\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		PsyrTA+=foundAbundance    
	#R 1 - 5 stk    
	elif re.search(r'RADAR\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		RADAR+=foundAbundance
	elif re.search(r'Retron\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Retron+=foundAbundance
	elif re.search(r'RexAB\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		RexAB+=foundAbundance    
	elif re.search(r'RM\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		RM+=foundAbundance
	elif re.search(r'RosmerTA\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		RosmerTA+=foundAbundance    
	#R 2 - 6 stk           
	elif re.search(r'Rst_2TM_1TM_TIR\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_2TM_1TM_TIR+=foundAbundance
	elif re.search(r'Rst_3HP\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_3HP+=foundAbundance
	elif re.search(r'Rst_DprA',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_DprAPPRT+=foundAbundance    
	elif re.search(r'Rst_DUF4238\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_DUF4238+=foundAbundance
	elif re.search(r'Rst_gop_beta_cll\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_gop_beta_cll+=foundAbundance
	elif re.search(r'Rst_HelicaseDUF2290\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_HelicaseDUF2290+=foundAbundance     
	#R3  - 5 stk 
	elif re.search(r'Rst_Hydrolase',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_Hydrolase3Tm+=foundAbundance  
	elif re.search(r'Rst_Old_Tin\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_Old_Tin+=foundAbundance  
	elif re.search(r'Rst_PARIS\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_PARIS+=foundAbundance  
	elif re.search(r'Rst_RT',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_RTnitrilaseTm+=foundAbundance  
	elif re.search(r'Rst_TIR',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Rst_TIRNLR+=foundAbundance    
	#S  - 8 stk
	elif re.search(r'SEFIR\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		SEFIR+=foundAbundance  
	elif re.search(r'Septu\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Septu+=foundAbundance  
	elif re.search(r'Shango\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Shango+=foundAbundance  
	elif re.search(r'Shedu\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Shedu+=foundAbundance  
	elif re.search(r'ShosTA\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		ShosTA+=foundAbundance
	elif re.search(r'SoFIC\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		SoFIC+=foundAbundance
	elif re.search(r'SspBCDE\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		SspBCDE+=foundAbundance
	elif re.search(r'Stk2\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Stk2+=foundAbundance    
	#T,U,V,Z - 6 stk        
	elif re.search(r'Thoeris\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Thoeris+=foundAbundance
	elif re.search(r'Tiamat\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Tiamat+=foundAbundance
	elif re.search(r'Uzume',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Uzume+=foundAbundance    
	elif re.search(r'Viperin\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Viperin+=foundAbundance
	elif re.search(r'Wadjet\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Wadjet+=foundAbundance
	elif re.search(r'Zorya\s',line):
		contig=re.search(r'\s([-\w.]+)\s', line).group(1)
		foundAbundance=findAbundance(contig)
		Zorya+=foundAbundance

                                           
#Setting name of sample
SampleName=0
SampleName=re.search(r'results/(\S+)',sys.argv[2]).group(1)

#Printing results:
print("S"+SampleName+"\t",Abi2,"\t",AbiEii,"\t",AbiH,"\t",Aditi,"\t",AVAST,"\t",Azaca,"\t",
Borvo,"\t",BREX,"\t",BstA,"\t",Bunzi,"\t",
Cas,"\t",Cas_Class1,"\t",Cas_Class2,"\t",CBASS,"\t",
DarTG,"\t",Dazbog,"\t",dCTPdeaminase,"\t",dGTPase,"\t",DISARM,"\t",Dnd,"\t",Dodola,"\t",DRT,"\t",Druantia,"\t",Dsr,"\t",Dynamins,"\t",
Gabija,"\t",Gao_Ape,"\t",Gao_Her,"\t",Gao_Hhe,"\t",Gao_Iet,"\t",Gao_Mza,"\t",Gao_Ppl,"\t",Gao_Qat,"\t",Gao_RL,"\t",Gao_TerY,"\t",Gao_Tmn,"\t",Gao_Upx,"\t",
Hachiman,"\t",Kiwa,"\t",Lamassu,"\t",LamassuFam,"\t",Lit,"\t",Menshen,"\t",Mokosh,"\t",
Nhi,"\t",NixI,"\t",Olokun,"\t",Pif,"\t",PrrC,"\t",PsyrTA,"\t",
RADAR,"\t",Retron,"\t",RexAB,"\t",RM,"\t",RosmerTA,"\t",
Rst_2TM_1TM_TIR,"\t",Rst_3HP,"\t",Rst_DprAPPRT,"\t",Rst_DUF4238,"\t",Rst_gop_beta_cll,"\t",Rst_HelicaseDUF2290,"\t",
Rst_Hydrolase3Tm,"\t",Rst_Old_Tin,"\t",Rst_PARIS,"\t",Rst_RTnitrilaseTm,"\t",Rst_TIRNLR,"\t",
SEFIR,"\t",Septu,"\t",Shango,"\t",Shedu,"\t",ShosTA,"\t",SoFIC,"\t",SspBCDE,"\t",Stk2,"\t",
Thoeris,"\t",Tiamat,"\t",Uzume,"\t",Viperin,"\t",Wadjet,"\t",Zorya)


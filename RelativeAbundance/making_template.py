#!/usr/bin/env python3
import sys,re

#Only initial run: for printing header:
outfile=open('defensesystems_Abundance.tsv','w')


#Printing results:
	#Only initial run: printing header

print("sample_name","\t","Abi2","\t","AbiEii","\t","AbiH","\t","Aditi","\t","AVAST","\t","Azaca","\t",
"Borvo","\t","BREX","\t","BstA","\t","Bunzi","\t",
"Cas","\t","Cas_Class1","\t","Cas_Class2","\t","CBASS","\t",
"DarTG","\t","Dazbog","\t","dCTPdeaminase","\t","dGTPase","\t","DISARM","\t","Dnd","\t","Dodola","\t","DRT","\t","Druantia","\t","Dsr","\t","Dynamins","\t",
"Gabija","\t","Gao_Ape","\t","Gao_Her","\t","Gao_Hhe","\t","Gao_Iet","\t","Gao_Mza","\t","Gao_Ppl","\t","Gao_Qat","\t","Gao_RL","\t","Gao_TerY","\t","Gao_Tmn","\t","Gao_Upx","\t",
"Hachiman","\t","Kiwa","\t","Lamassu","\t","Lamassu-Fam","\t","Lit","\t","Menshen","\t","Mokosh","\t",
"Nhi","\t","NixI","\t","Olokun","\t","Pif","\t","PrrC","\t","PsyrTA","\t",
"RADAR","\t","Retron","\t","RexAB","\t","RM","\t","RosmerTA","\t",
"Rst_2TM_1TM_TIR","\t","Rst_3HP","\t","Rst_DprA-PPRT","\t","Rst_DUF4238","\t","Rst_gop_beta_cll","\t","Rst_HelicaseDUF2290","\t",
"Rst_Hydrolase-3Tm","\t","Rst_Old_Tin","\t","Rst_PARIS","\t","Rst_RT-nitrilase-Tm","\t","Rst_TIR-NLR","\t",
"SEFIR","\t","Septu","\t","Shango","\t","Shedu","\t","ShosTA","\t","SoFIC","\t","SspBCDE","\t","Stk2","\t",
"Thoeris","\t","Tiamat","\t","Uzume","\t","Viperin","\t","Wadjet","\t","Zorya", file=outfile)




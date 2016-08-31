def extract_pairs(variable, putArray, minmax):
    arr_num = 1 if (minmax == "min" or minmax == "Min") else 2
    cutname = variable[0].replace(minmax, "")
    inserted = False

    for i, name in enumerate(putArray):
        if name[0] == cutname:
            putArray[i][arr_num] = variable[1]
            inserted = True
            break
    if not inserted:
        if arr_num == 1: putArray.append([cutname, variable[1], ""])
        else: putArray.append([cutname,"", variable[1]])




import os.path
import sys


directory= "PartDet"
filename = sys.argv[1]

if not os.path.isfile(filename):
    print "file not found"
    exit(1)

f = open(filename, 'r')

Gen = []
Single = [[[], [], [], [], [], []], [[],[]], [[],[]], [[],[]]]
Smear = [[["MuonMatchingDeltaR","0.3"], ["ElectronMatchingDeltaR","0.3"], ["TauMatchingDeltaR","10"]], [], [], []]
Cuts = []
Dipart = [[[],[],[],[]], [[],[],[],[]], [[],[],[],[]]]
VBF = []
Run = []

file_org = [["Muon1Tau1","Muon1Tau2","Muon2Tau1","Muon2Tau2"], 
            ["Electron1Tau1","Electron1Tau2","Electron2Tau1","Electron2Tau2"],
            ["DiElectron","DiMuon","DiTau","DiJet"]]  

sing_org = [["Jet1", "Jet2", "CentralJet", "FirstLeadingJet", "SecondLeadingJet", "BJet"],
            ["Tau1", "Tau2"], ["Muon1", "Muon2"], ["Electron1", "Electron2"]]

smear_org = ["Jet", "Tau", "Muon", "Electron"]

run_stat = ["CalculatePUSystematics","DataHistos", "MCHistos", "isData", "Trigger1FirstRequirement", 
            "Trigger1SecondRequirement", "Trigger2FirstRequirement", "Trigger2SecondRequirement", "TreatMuonsAsNeutrinos",
            "JetPtForMhtAndHt", "JetEtaForMhtAndHt", "ApplyJetLooseIDforMhtAndHt"]
partID ="TauID 15.\nTauStatus 2.\n\nTopID 6.\nTopStatus 2.\n\nElectronID 11.\nElectronStatus 1.\n\nMuonID 13.\nMuonStatus 1.\n\nZID 23.\nZStatus 2.\n\nWID 24.\nWStatus 2.\n\nHiggsID 25.\nHiggsStatus 2.\n"


for line in f:
    variable = line.split()
    if len(variable) == 0:
        continue
         
##### Multiplicity #####  
    if "Nm" in variable[0]:
        variable[0] = variable[0].replace("N", "")
        if "min" in variable[0]: extract_pairs(variable, Cuts, "min")
        else: extract_pairs(variable, Cuts, "max")
        continue

    elif "Fill" in variable[0]: continue

##### Gen #####
    elif "GenTau" in variable[0]:
        variable[0] = variable[0].replace("Gen","")
        if "Min" in variable[0]: extract_pairs(variable, Gen, "Min")
        elif "Max" in variable[0]: extract_pairs(variable, Gen, "Max")
        else: Gen.append([variable[0], variable[1], ""])
        continue
           
    used = False
##### DiParticle #####
    for i, group in enumerate(file_org):
        for j, item in enumerate(group):
            if item in variable[0] and not ("Lead" in variable[0]):
                variable[0] = variable[0].replace(item,"").replace("Do","")
                if "Min" in variable[0]: extract_pairs(variable, Dipart[i][j], "Min")
                elif "Max" in variable[0]: extract_pairs(variable, Dipart[i][j], "Max")
                else: Dipart[i][j].append([variable[0], variable[1], ""])
                used = True
                break
        if used: break
    if used: continue


##### Single Particle #####
    for i, group in enumerate(sing_org):
        for j, item in enumerate(group):
            if item in variable[0]:
                variable[0] = variable[0].replace(item,"").replace("Reco","").replace("DoDiscrByIsZllCut","DiscrIfIsZdecay").replace("Selects","SelectTaus")
                variable[0] += "Cut" if variable[0] == "Pt" else ""
                if "Met" in variable[0]: variable[0] = variable[0].replace("Do","")
                if "MaxIsolation" in variable[0] or "MinIsolation" in variable[0]:
                    Single[i][j].append([variable[0], variable[1], ""])
                elif "Min" in variable[0]: extract_pairs(variable, Single[i][j], "Min")
                elif "Max" in variable[0]: extract_pairs(variable, Single[i][j], "Max")
                else: Single[i][j].append([variable[0], variable[1], ""])
                used = True
                break
        if used: break
    if used: continue

##### Scraps #####
    if "JetB" in variable[0] or variable[0] == "MatchBToGen":
        Single[0][5].append([variable[0], variable[1], ""])
        continue

##### Smear #####                
    if ("Smear" in variable[0] or "Offset" in variable[0] or "Match" in variable[0] or "Mother" in variable[0]) and not ("MatchBToGen" == variable[0]):
        for i, part in enumerate(smear_org):
            if part in variable[0]: 
                if part == "Jet" : part = ""
                
                Smear[i].append([variable[0].replace(part,""), variable[1]])
        continue
        
##### RUN #####
    for value in run_stat:
        if value in variable[0]:
            Run.append([variable[0], variable[1]])
            used = True
            break
    if used: continue

##### VBF #####
    variable[0] = variable[0].replace("DoSUSY","").replace("LeadDiJet","").replace("Do","")
    if "Min" in variable[0]: extract_pairs(variable, VBF, "Min")
    elif "Max" in variable[0]: extract_pairs(variable, VBF, "Max")
    elif "LooseID" in variable[0]:
        Single[0][3].append(["ApplyLooseID", variable[1], ""])
        Single[0][4].append(["ApplyLooseID", variable[1], ""])

    else: VBF.append([variable[0], variable[1], ""])




###########################################
##############  WRITE OUT  ################
###########################################

for i, nfile in enumerate(["Jet", "Tau", "Muon", "Electron"]):
    filer = open(directory + "/" + nfile + "_info.in", 'w')
    
    filer.write("Smear\n\n")
    
    for item in Smear[i]:
        item[0] += "Particle" if item[0] == "SmearThe" else ""
        if "ToGenMatchingDeltaR" in item[0]: item[0] = item[0].replace("To","")
        if "MotherId" in item[0]: 
            item[0] = item[0].replace("Id", "ID")
            if "MotherID" == item[0] and not "." in item[1]: item[1] += "."
        filer.write('{:27}  {:>6}\n'.format(item[0], item[1]))
    

    for j, group in enumerate(sing_org[i]):
        if nfile == "Electron": group = group.replace("tron","")
        filer.write("\n" + group + "\n\n")
        for item in Single[i][j]:
            item[0] += "ThisJet" if item[0] == "DoDiscrBy" else ""
            filer.write('{:27}  {:>6}    {}\n'.format(item[0], item[1], item[2]))

    filer.close()
    
    


for i, nfile in enumerate(["MuonTau", "ElectronTau", "DiParticle"]):
    filer = open(directory + "/" + nfile + "_info.in", 'w')
    for j, subs in enumerate(file_org[i]):
        filer.write(subs + "\n\n")
        foundmasscalc = 0
        for entry in Dipart[i][j]:
            if "Use" in entry[0]:
                foundmasscalc += 1
                cutname = entry[0].replace("Use","").replace("MassReco","")
                if entry[1] == '1': 
                    filer.write('{:25}  {}\n'.format("HowCalculateMassReco", cutname))
                    foundmasscalc = -1
                elif foundmasscalc == 2:
                    filer.write('{:25}  {}\n'.format("HowCalculateMassReco", "None"))
            else: filer.write('{:25}  {:>6}    {}\n'.format(entry[0], entry[1], entry[2]))
        filer.write("\n")
    filer.close()
        


cutf = open(directory + "/Cuts.in", 'w')

for entry in Cuts:
    maxval = entry[2] if (len(entry[2]) != 0 and len(entry[2]) < 3)  else -1
    cutf.write('N{:25}  {:>2}  {:2}\n'.format(entry[0], entry[1], maxval))
cutf.close()

runf = open(directory + "/Run_info.in", 'w')
runf.write("Run\n\n")
for entry in Run:
    runf.write('{:27}  {:>2}\n'.format(entry[0], entry[1]))
runf.close()


vbff = open(directory + "/VBFCuts_info.in", 'w')
vbff.write("VBFSUSY\n\n")
for entry in VBF:
    vbff.write('{:25}  {:>6}  {:>6}\n'.format(entry[0], entry[1], entry[2]))
vbff.close()

genf = open(directory + "/Gen_info.in", 'w')
genf.write("Gen\n\n")
for entry in Gen:
    genf.write('{:8}  {:>6}  {:>6}\n'.format(entry[0], entry[1], entry[2]))
genf.write(partID)
genf.close()


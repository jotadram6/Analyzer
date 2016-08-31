import os

f=open("PartDet/VBFCuts_info.in", 'r')
ap=open("PartDet/Run_info.in", 'a')

string=""

for line in f:
    if "Met" in line or "HT" in line or "Mht" in line or "Ht" in line:
        ap.write(line)
        continue
    string += line

f.close()
os.remove("PartDet/VBFCuts_info.in")

f2 = open("PartDet/VBFCuts_info.in", 'w')
f2.write(string)

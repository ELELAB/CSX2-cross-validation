#!/usr/bin/env python
# By Mads Nygaard
# Script for calculating chi-square differences between the output of PPM_one in the bb_details.dat and proton_details.dat file and experimental values in NMRstar format.
# Updated with 1/(n-1) for average of species (2016-06-13)

# Dependencies
import subprocess
import os.path
import argparse

# Global variables
species = ["H", "CA", "CB", "CG", "CD", "CD1", "CD2", "CE", "CG1", "CG2", "C", "N", "HA", "HA2", "HA3", "HB2", "HB3", "HB", "HD1", "HD2", "HD3", "HE", "HG1", "HG2", "HG3", "HD", "HE", "HZ"]
residues = ["Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr"]

ppmfilename = "results_ppm"
oldppmfilename = "results_oldppm"
CH3filename = "results_CH3"
Arfilename = "results_Ar"


#Not used
#factors = {"N":2.91,
#       "CA":1.06,
#       "CB":1.23,
#       "C":1.32,
#       "H":0.53,
#       "HB":1.0,
#       "HD1":1.0,
#       "HD2":1.0,
#       "HE":0.01,
#       "HG1":1.0,
#       "HG2":1.0}

## Functions

# Create dictionary with defined species
def createdict():
    newdict = {}
    for n in species:
        newdict[n] = {}
    return newdict


# Float checker
def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


# Runs newPPM_One in bash from frame begin to frame stop of the pdb file
def ppmbashrun(begin, stop, pdbfile):
    f = open("ppmlog.txt", "w")
    subprocess.call(["./ppm_linux", "-pdb", str(pdbfile), "-begin", str(begin), "-stop", str(stop), " > ppmlog.log"], stdout=f)
    f.close()


# Edits the command.cmd file nessecary to run CH3Shift and ArShift "number" of times
def editCommFile(number):
    f = open("command.cmd", "r")
    fileline = f.readlines()
    f.close()
    for idx, l in enumerate(fileline):
        if "examine" in l:
            line = l.split()
            s = list(line[1])
            s[1:-1] = str(number)          
            s.append("\n")
            line[1] = "".join(s)
            newline = "      ".join(line)
            fileline[idx] = newline
    f = open("command.cmd", "w")
    f.writelines(fileline)
    f.close()


# Counts the number of structures in a PDB file (counts occurance of ENDMDL)
def countModelnum(pdbfile):
    f = open(pdbfile, "r")
    filelines = f.readlines()
    f.close()
    num = 0
    for l in filelines:
        if "ENDMDL" in l:
            num += 1
    return num


# Runs ArShift.R or CH3Shift.R scripts
def runR(rscript):
    subprocess.call("R CMD BATCH " + rscript, shell=True)


# Adds output from ArShift or CH3Shift to dict
def sideChShift_pred(filename, ppmdict):
    f = open(filename, "r")
    lines = f.readlines()
    for l in lines:
        attachSideChAvg(l, ppmdict)
    f.close()
    return ppmdict


# Reading Line after line for attachlineavg
def attachfiledetail(filename, ppmdict):
    f = open(filename, "r")
    lines = f.readlines()
    for l in lines:
        attachlineavg(l, ppmdict)
    f.close()
    return ppmdict


# Read data line after line from bmrb file to dictionary            
def attachfile(filename, ppmdict):
    if os.path.splitext(filename)[1] != ".dat":
        filename = os.path.splitext(filename)[0]+".dat" 
    f = open(filename, "r")
    lines = f.readlines()
    for l in lines:
        attachline(l, ppmdict)
    f.close()
    return ppmdict


# Reading Line after line for attachlineavg
def attachfilenewPPM(filename, ppmdict):
    f = open(filename, "r")
    lines = f.readlines()
    for l in lines:
        attachlinenewPPM(l, ppmdict)
    f.close()
    return ppmdict


# Adds values from line to dict newPPM, works in place
def attachlinenewPPM(fileline, dictionary):
    count = 0
    ppm = 1
    linespecies = [""]
    for splitline in fileline.split():
        if splitline.isdigit() and count < 2:
            count += 1
            resnum = int(splitline)
        #elif (splitline not in residues) and splitline.isdigit() == False:
        elif splitline in species:
            linespecies.append(splitline)
        elif ('.' in splitline) and count == 2:
            ppm = float(splitline)
            break
    if max(linespecies, key=len) in species:
        dictionary[max(linespecies, key=len)].setdefault(resnum, list()).append(ppm)
    else:
        return None


# Adds values from line of ArShift or CH3Shift output to dict
def attachSideChAvg(fileline, dictionary):
    linesplit = fileline.split()
    for l in linesplit:
        if "NOTE:" in l:
            break
        elif "RESULT:" in l and isFloat(linesplit[12]):
            if int(linesplit[4]) in dictionary[linesplit[10]]:
                avg = (float(linesplit[12])+sum(dictionary[linesplit[10]][int(linesplit[4])]))/(len(dictionary[linesplit[10]][int(linesplit[4])])+1)
                dictionary[linesplit[10]][int(linesplit[4])].append(avg)
            else:
                dictionary[linesplit[10]][int(linesplit[4])] = [float(linesplit[12])]


# Runs given Rscript (ArShift.R or CH3Shift.R) N times, dependent on model N in PDB file, outputs in dict
def runSideChOnPdbInR(Rscript):
    nummod = countModelnum("example.pdb")
    numrange = range(1, nummod+1)
    print "Running %s times of %s" % (nummod, Rscript)
    predictedValues = createdict()
    while nummod > 0:
        print "Running %s of %s times" % (numrange[-nummod], numrange[-1])
        editCommFile(numrange[-nummod])
        runR(Rscript)
        sideChShift_pred("out.txt", predictedValues)
        nummod -= 1
    return predictedValues


def runNewPpmOneBash():
    nummod = countModelnum("example.pdb")
    numrange = range(1, nummod+1)
    print "Running %s times of PPM_One predictor" % (nummod)
    predictedValues = createdict()
    while nummod > 0:
        print "Running %s of %s times" % (numrange[-nummod], numrange[-1])
        ppmbashrun(0, numrange[-nummod], "example.pdb")
        attachfilenewPPM("bmrb_pre.dat", predictedValues)
        nummod -= 1
    return predictedValues


# Loading line from bb_details.dat cummulated averaging
def attachlineavg(fileline, dictionary):
    #count = 0
    totalppm = 0
    avglst = []
    splitline = fileline.split()
    residnr = int(splitline.pop(0))  # ppm_details file is shifted one residue??
    residspec = splitline.pop(0)
    splitline.pop(0)
    for idx, num in enumerate(splitline):
        num = float(num)
        totalppm += num
        avglst.append(totalppm/(idx+1))
    if args.v is True:
        print ("Pred. AVG of %s %s %s" % (residnr, residspec, str(avglst[-1])))
    dictionary[residspec][residnr] = avglst


# Read line from file insert to dictionary (used for loadig nmrSTAR file)
def attachline(fileline, dictionary):
    count = 0
    ppm = 1
    linespecies = [""]
    for splitline in fileline.split():
        if splitline.isdigit() and count < 2:
            count += 1
            resnum = int(splitline)+args.shift
        #elif (splitline not in residues) and splitline.isdigit() == False:
        elif splitline in species:
            linespecies.append(splitline)
        elif splitline.lower() in (res.lower() for res in residues):
            resid = splitline
        elif ('.' in splitline) and count == 2:
            ppm = float(splitline)
            break
    if max(linespecies, key=len) in species:
        dictionary[max(linespecies, key=len)].setdefault(resnum, list()).append(ppm)
        dictionary[max(linespecies, key=len)][resnum].append(resid)
    else:
        return None


# Calculate Chisquared for every value read
def calcchi(predicdict, meassdict):
    resultdict = createdict()
    for specieskey in predicdict:
        for ppmlistkey in predicdict[specieskey]:
            for idx, ppm in enumerate(predicdict[specieskey][ppmlistkey]):
                try:
                    chisquare = ((ppm - meassdict[specieskey][ppmlistkey][0]) ** 2)/abs(meassdict[specieskey][ppmlistkey][0])
                    resultdict[specieskey].setdefault(ppmlistkey, list()).append(chisquare)
                except KeyError:
                    pass
                try:
                    if args.v is True:
                        print ("Avg, last X^2 val of %s no %s  %s  %s" % (str(specieskey), str(ppmlistkey), 
                               str((sum(resultdict[specieskey][ppmlistkey])/len(resultdict[specieskey][ppmlistkey]))), 
                               str(resultdict[specieskey][ppmlistkey][-1])))  # prints avg Chisq of the different species in column 1 and Chisq value of the last data point in column 2
                except KeyError:
                    pass
    empty_keys = [k for k, v in resultdict.iteritems() if not v]  # Cleaning up empty keys
    for k in empty_keys:
            del resultdict[k]
    return resultdict


def writeBfactor(results, measdict, runtype, chains):
    with open("ChiSq_detail"+str(runtype)+".dat", "w") as f:
        f.write("#CHAIN\t<RESI>\t<RESN>\t[name]\t<data>\n")
        for cidx, chain in enumerate(chains):
            for specieskey in results:
                for ppmlistkey in results[specieskey]:
                    try:
                        f.write("%s\t%s\t%s\t%s\t%.12f\n" % (
                               (chain,
                                str(ppmlistkey),
                                str(measdict[specieskey][ppmlistkey][1]),
                                str(specieskey),
                                results[specieskey][ppmlistkey][-(len(chains)-cidx)])))
                    except KeyError:
                        pass


# Average Chisquared values for species
def chiavgspec(chidict):
    #print chidict
    chisquaredict = {}
    for key in chidict:
        chisquaredict[key] = []
    #Count no of prediction values (seems weird)
    try:
        prednum = len(next(iter(chidict[next(iter(chidict.keys()))].values()))) 
    except StopIteration:
        print "Check that your chemical shifts file correspond to the sequence numbers in the structure file"
        print chidict
    for idx in range(prednum):
        for specieskey in chidict:
            #speciesavg = 0.0
            specielist = []
            for resid in chidict[specieskey]:
                specielist.append(chidict[specieskey][resid][idx])
            if specielist != []:
                if len(specielist) == 1:
                    chisquaredict[specieskey].append(sum(specielist)/(len(specielist)))
                else:
                    chisquaredict[specieskey].append(sum(specielist)/(len(specielist)-1))
                print specielist, len(specielist), specieskey
            else:
                chisquaredict.pop(specieskey, 0)
    return chisquaredict


# Print to file results.xvg
def printtofile(dictionary, filename, chains):
    nmol = len(chains)
    for n, chain in enumerate(chains):
        f = open(filename + "_" + str(chain) + ".xvg", "w")
        f.write("""@    title "ChiSquare"\n@    xaxis  label "Timestep"\n@    yaxis  label "Sum of Chisquare"\n@ TYPE xy\n""")
        num = 0 
        for idx2, keys in enumerate(sorted(dictionary)):
            if len(dictionary[keys]) > 0:
                f.write("""@ s%s legend "%s"\n""" % (num, keys))
                num += 1
        alistlength = len(next(iter(dictionary.values())))/nmol
        for idx in range(alistlength):
            f.write(str(idx)+"\t")
            for spec in sorted(dictionary):
                try:
                    f.write("%.14f" % dictionary[spec][n::nmol][idx] + "\t")  # Take every nmol(2nd, 3rd...) chemical shift value and write
                except IndexError:  # Should not be nessecary any more
                    print "Something strange happened"
                    pass
            f.write("\n")
        f.close()
        print ("Output is located in %s_%s.xvg and can be plotted using xmgrace -nxy.") % (filename, n)


# Cleaning up empty keys in a dict
def keyCleanup(dictionary):
    empty_keys = [k for k, v in dictionary.iteritems() if not v]  # Cleaning up empty keys
    for k in empty_keys:
        del dictionary[k]
    return dictionary


def runppm_new(bmrbfile, chains):
    measdict = attachfile(bmrbfile, createdict())
    ppmpredict = runNewPpmOneBash()
    results = calcchi(ppmpredict, measdict)
    writeBfactor(results, measdict, "ppm", chains) 
    speciesresults = chiavgspec(results)
    return speciesresults


def runppm_onedetail(bmrbfile):
    measdict = attachfile(bmrbfile, createdict())
    print ("."),
    ppmprefromfile = attachfiledetail("bb_details.dat", createdict())
    print ("."),
    ppmprefromfile = attachfiledetail("proton_details.dat", ppmprefromfile)
    print ("."),
    results = calcchi(ppmprefromfile, measdict)
    print ("."),
    speciesresults = chiavgspec(results)
    print (".")
    return speciesresults


def runSideCh(bmrbfile, CHorAR):
    measdict = attachfile(bmrbfile, createdict())
    preddict = runSideChOnPdbInR(CHorAR)
    keyCleanup(preddict)    
    results = calcchi(preddict, measdict)
    speciesresults = chiavgspec(results)
    print (".")
    return speciesresults


if __name__ == "__main__":
    #Argument Parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, 
        description=('''\
This script calculates X^2 values based on experimental values and predicted values.
Chisquare is calculated as:

     ((CS_pred - CS_exp)**2)/(CS_exp) 

..where CS_pred is the cummulative average predicted chemical shifts and CS_exp is the experimentally obtained chemical shifts.
Predictors can be PPM (old version, needs pre-run of ppm_linux generating proton_details.dat and bb_details.dat files (recommended for first run)). 
PPM_One new version (ppm_linux needed in folder and example.pdb structure file).

'''))
    parser.add_argument('-f', metavar='file.dat', default='ChemicalShifts.dat',
                        help='Experimental chemical shifts file in NMRstar format (Default: ChemicalShifts.dat)', type=str)
    parser.add_argument('-shift', metavar='N', default=0, type=int,
                        help='Shift the sequence number of the experimental chemical shifts by N residues')
    parser.add_argument('-chains', nargs="+", metavar='A B', default="A", type=str, help='Chains in the pdb (only supported for PPM_new)')
    parser.add_argument('--ppm', action='store_true', help='Run new PPM_One prediction')
    parser.add_argument('--ppmold', action='store_true', help='Run old PPM prediction (needs bb_details.dat and proton_details.dat)')
    parser.add_argument('--v', action='store_true', help='Run the script in verbose')
    args = parser.parse_args()
    parser.print_help()

    #Run ppm_one script
    if args.ppm is True:
        if len(args.chains) == 1:
            print "\n NOTE:Remember to specify -chain if number of molecules in .pdb is more than 1\n"
        chisqppm = runppm_new(args.f, args.chains)
        printtofile(chisqppm, ppmfilename, args.chains)

    #Run ppm old script
    if args.ppmold is True:
        chisqold = runppm_onedetail(args.f)
        printtofile(chisqold, oldppmfilename, args.chains)

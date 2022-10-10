#!/usr/local/bin/python3
 
import os
import re
import shutil
from collections import defaultdict
import statistics
import matplotlib.pyplot as plt
import pandas as pd

def pathmaker(token):
    mainpath = os.path.dirname(os.path.dirname(__file__))
    print(mainpath)
    global scriptspath
    scriptspath = os.path.join(mainpath,"Scripts")
    champ = input("What is the path to the results directory you want to analyse ?")
    global resultspath
    resultspath = os.path.join(os.path.join(champ,"Python"),token)
    global samplespath
    samplespath = os.path.join(os.path.join(champ,"ChAMP"),token)
    return scriptspath,resultspath,samplespath

def check():
    if os.path.exists(resultspath) :
        rmcheck = input("Results directory already exists : do you want to remove it ? (Yes/No)")
        while rmcheck != "Yes" or "No":
            if rmcheck == "Yes" :
                shutil.rmtree(resultspath)
                break
            elif rmcheck == "No" :
                print("Exciting script to allow backup")
                exit()
            else:
                print("Error : wrong input!")
                rmcheck = input("Results directory already exists : do you want to remove it ? (Yes/No)")
    os.makedirs(resultspath)

def concat():
    for files in os.walk(samplespath):
        samplesfiles = files[2]
    global header
    header = ["sample"]+[i.strip("\n").strip("\"") for i in open(os.path.join(samplespath,samplesfiles[0]),"r").readline().split(",")]+["methylation\n"]
    header[18] = "betaAv-CT"
    header[19] = "betaAv-Patient"
    header = ",".join(header)
    global concatfile
    concatfile = os.path.join(resultspath,"DMR-All.csv")
    open(concatfile,"w").writelines(header)
    open(concatfile,"a").writelines([",".join([filein,",".join([i.strip("\n").strip("\"") for i in line.split(",")]),"hypomethylation\n"]) for filein in samplesfiles for line in open(os.path.join(samplespath,filein),"r") if not line.startswith("\"seqnames\"") and float((",".join([i.strip("\n").strip("\"") for i in line.split("\t")])).split(",")[19]) < 0])
    open(concatfile,"a").writelines([",".join([filein,",".join([i.strip("\n").strip("\"") for i in line.split(",")]),"hypermethylation\n"]) for filein in samplesfiles for line in open(os.path.join(samplespath,filein),"r") if not line.startswith("\"seqnames\"") and float((",".join([i.strip("\n").strip("\"") for i in line.split("\t")])).split(",")[19]) > 0])
    global patientlist
    patientlist = sorted(set([line.split(",")[0] for line in open(concatfile,"r").readlines()[1:]]))
    rmnogenes = input("Do you want to remove DMR without genes ? (Yes/No)")
    while rmnogenes != "No":
        if rmnogenes == "Yes" :
            open(os.path.join(resultspath,"DMR-All-NoEmptyGene.csv"), "w").writelines([line for line in open(concatfile, "r") if not line.split(",")[17] == ""])
            concatfile = os.path.join(resultspath,"DMR-All-NoEmptyGene.csv")
            break
        else :
            print("Error : wrong input!")
            rmnogenes = input("Do you want to remove DMR without genes ? (Yes/No)")
    for i in patientlist:
        open(os.path.join(resultspath,"%s.csv" %i),"w").writelines(header)
        linelist = [line for line in open(concatfile,"r").readlines()[1:] if line.split(",")[0] == i]
        linelist.sort(key=lambda x: x.split(",")[17])
        open(os.path.join(resultspath,"%s.csv" %i),"a").writelines(line for line in linelist)

def globalmeth():
    globalpath = os.path.join(resultspath,"Global")
    if not os.path.exists(globalpath):
        os.makedirs(globalpath)
    genevalues = []
    genetablevalues = []
    dmrvalues = []
    dmrtablevalues =[]
    for i in patientlist:
        meth = [(line.split(",")[21].strip("\n"),line.split(",")[17]) for line in open(os.path.join(resultspath,"%s.csv" %i),"r").readlines()[1:]]
        dmrnum = len(meth)
        hypodmr = sum(1 for j in meth if j[0] == "hypomethylation")
        hyperdmr = sum(1 for j in meth if j[0] == "hypermethylation")        
        genesnum = {"hypomethylation":[],"hypermethylation":[]}
        for k in set(meth):
            genesnum[k[0]].extend(k[1].split(";"))
        bothgenes = set([j for j in genesnum.get("hypomethylation") if j in genesnum.get("hypermethylation")])
        lenhypogenes = len(set(genesnum.get("hypomethylation")))
        lenhypergenes = len(set(genesnum.get("hypermethylation")))
        lenbothgenes = len(bothgenes)
        lengenesnum = len(set([j for key,values in genesnum.items() for j in values]))
        genevalues.append([lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes])
        genetablevalues.append([str(lenhypogenes-lenbothgenes)+" ("+str(round(100*(lenhypogenes-lenbothgenes)/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",\
                                str(lenhypergenes-lenbothgenes)+" ("+str(round(100*(lenhypergenes-lenbothgenes)/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",\
                                str(lenbothgenes)+" ("+str(round(100*lenbothgenes/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",\
                                str(lenhypogenes+lenhypergenes-lenbothgenes)])
        dmrvalues.append([hypodmr,hyperdmr])
        dmrtablevalues.append([str(hypodmr)+" ("+str(round(100*hypodmr/(hypodmr+hyperdmr),2))+"%)",\
                               str(hyperdmr)+" ("+str(round(100*hyperdmr/(hypodmr+hyperdmr),2))+"%)",\
                               str(hypodmr+hyperdmr)])
        if lenbothgenes != 0:
            print("Beware : %s genes found both hypomethylated and hypermethylated for %s" %(lenbothgenes,i))
            open(os.path.join(globalpath,"%s-Global-Statistics-HypoHyper.csv" %i),"w").write(header)
            open(os.path.join(globalpath,"%s-Global-Statistics-HypoHyper.csv" %i),"a").writelines([line for line in open(os.path.join(resultspath,"%s.csv" %i),"r").readlines()[1:] if any(j in line.split(",")[17].split(";") for j in bothgenes)])
        df = pd.DataFrame([[lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes],[hypodmr,hyperdmr]],columns=["Hypomethylated","Hypermethylated","Both"],index=["Genes","DMR"])
        ax = df.plot.bar(stacked=True,title="Global Statistics of %s" %i,color=["#428bca","#d9534f","#5cb85c"])
        plt.ylabel("Number of genes/DMR")
        plt.grid()
        plt.xticks([],horizontalalignment="center",rotation="horizontal")
        plt.table(cellText=[[str(lenhypogenes-lenbothgenes)+" ("+str(round(100*(lenhypogenes-lenbothgenes)/lengenesnum,2))+"%)",str(hypodmr)+" ("+str(round(100*hypodmr/dmrnum,2))+"%)"],\
                            [str(lenhypergenes-lenbothgenes)+" ("+str(round(100*(lenhypergenes-lenbothgenes)/lengenesnum,2))+"%)",str(hyperdmr)+" ("+str(round(100*hyperdmr/dmrnum,2))+"%)"],\
                            [str(lenbothgenes)+" ("+str(round(100*lenbothgenes/lengenesnum,2))+"%)","X"],\
                            [str(lengenesnum),dmrnum]],\
                  colLabels=["Genes","DMR"],rowLabels=["Hypomethylated","Hypermethylated","Both","Total"],cellLoc = "center", rowLoc = "center")
        # for c in ax.containers:
        #     ax.bar_label(c,label_type="center",fontsize="small")
        plt.savefig(os.path.join(globalpath,"%s-Global-Statistics.png" %i),bbox_inches="tight",dpi=800)
        plt.close()
    genetablevaluesbis = [[j[i] for j in genetablevalues] for i in range(4)]
    df = pd.DataFrame(genevalues,columns=["Hypomethylated","Hypermethylated","Both"],index=list(patientlist))
    ax = df.plot.bar(stacked=True,title="Global Statistics for Genes",color=["#428bca","#d9534f","#5cb85c"])
    plt.ylabel("Number of genes")
    plt.grid()
    plt.xticks([],horizontalalignment="center",fontsize="small")
    plt.table(cellText=genetablevaluesbis,colLabels=list(patientlist),rowLabels=["Hypomethylated","Hypermethylated","Both","Total"],cellLoc="center",rowLoc="center")
    plt.savefig(os.path.join(globalpath,"Global-Statistics-Genes.png"),bbox_inches="tight",dpi=800)
    plt.close()
    dmrtablevaluesbis =[[j[i] for j in dmrtablevalues] for i in range(3)]    
    df = pd.DataFrame(dmrvalues,columns=["Hypomethylated","Hypermethylated"],index=list(patientlist))
    ax = df.plot.bar(stacked=True,title="Global Statistics for DMR",color=["#428bca","#d9534f"])
    plt.ylabel("Number of DMR")
    plt.grid()
    plt.xticks([],horizontalalignment="center",fontsize="small")
    plt.table(cellText=dmrtablevaluesbis,colLabels=list(patientlist),rowLabels=["Hypomethylated","Hypermethylated","Total"],cellLoc="center",rowLoc="center")
    plt.savefig(os.path.join(globalpath,"Global-Statistics-DMR.png"),bbox_inches="tight",dpi=800)
    plt.close()

def common():
    comparelist = input("Which patient do you want to compare ? (Choose in this list %s and separate with comma)" %patientlist).split(",")
    commonpath = os.path.join(resultspath,"Common")
    if not os.path.exists(os.path.join(commonpath,"%s" %comparelist)):
        os.makedirs(os.path.join(commonpath,"%s" %comparelist))
    patientgenes = list(set(tuple(k) for k in [[line.split(",")[i].strip("\n") for i in [0,17,21]] for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:]]))
    genes = set([i for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:] for i in line.split(",")[17].split(";")])
    patientdict = defaultdict(list)
    for k in patientgenes:
        patientdict[k[0]].extend(k[1].split(";"))
    commongenes = []
    for i in genes:
        temp = [key for key,values in patientdict.items() if i in values]
        if len(set(temp)) == len(set(comparelist)):
            commongenes.append(i)
    genevalues = []
    genetablevalues = []
    dmrvalues = []
    dmrtablevalues =[]
    for i in comparelist:
        meth = [(line.split(",")[21].strip("\n"),line.split(",")[17]) for line in open(os.path.join(resultspath,"%s.csv" %i),"r").readlines()[1:]]
        dmrnum = len(meth)
        hypodmr = sum(1 for j in meth if j[0] == "hypomethylation" and any(k in j[1].split(";") for k in commongenes))
        hyperdmr = sum(1 for j in meth  if j[0] == "hypermethylation" and any(k in j[1].split(";") for k in commongenes))
        genesnum = {"hypomethylation":[],"hypermethylation":[]}
        for k in set(meth):
            genesnum[k[0]].extend(k[1].split(";"))
        if os.path.exists(os.path.join(os.path.join(resultspath,"Global"),"%s-Global-Statistics-HypoHyper.csv" %i)):
            bothgenes = set([j for line in open(os.path.join(os.path.join(resultspath,"Global"),"%s-Global-Statistics-HypoHyper.csv" %i),"r").readlines()[1:] for j in line.split(",")[17].split(";") if j in commongenes])
        else:
            bothgenes = set([j for j in genesnum.get("hypomethylation") if j in genesnum.get("hypermethylation") and j in commongenes])
        lenhypogenes = len(set([j for j in genesnum.get("hypomethylation") if j in commongenes]))
        lenhypergenes = len(set([j for j in genesnum.get("hypermethylation") if j in commongenes]))
        lenbothgenes = len(bothgenes)
        lengenesnum = len(set([j for key,values in genesnum.items() for j in values]))
        if lenbothgenes != 0:
            print("Beware : %s genes found both hypomethylated and hypermethylated for %s" %(lenbothgenes,i))
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"%s-Common-Statistics-HypoHyper.csv" %i),"w").write(header)
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"%s-Common-Statistics-HypoHyper.csv" %i),"a").writelines([line for line in open(os.path.join(resultspath,"%s.csv" %i),"r").readlines()[1:] if any(j in line.split(",")[17].split(";") for j in bothgenes)])
        if lenhypogenes or lenhypergenes != 0:
            genevalues.append([lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes])
            genetablevalues.append([str(lenhypogenes-lenbothgenes)+" ("+str(round(100*(lenhypogenes-lenbothgenes)/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",\
                                    str(lenhypergenes-lenbothgenes)+" ("+str(round(100*(lenhypergenes-lenbothgenes)/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",\
                                    str(lenbothgenes)+" ("+str(round(100*lenbothgenes/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",\
                                    str(lenhypogenes+lenhypergenes-lenbothgenes)+"/"+str(lengenesnum)+" ("+str(round(100*(lenhypogenes+lenhypergenes-lenbothgenes)/lengenesnum,2))+"%)"])
            dmrvalues.append([hypodmr,hyperdmr])
            dmrtablevalues.append([str(hypodmr)+" ("+str(round(100*hypodmr/(hypodmr+hyperdmr),2))+"%)",\
                                str(hyperdmr)+" ("+str(round(100*hyperdmr/(hypodmr+hyperdmr),2))+"%)",\
                                str(hypodmr+hyperdmr)+"/"+str(dmrnum)+" ("+str(round(100*(hypodmr+hyperdmr)/dmrnum,2))+"%)"])
            df = pd.DataFrame([[lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes],[hypodmr,hyperdmr]],columns=["Hypomethylated","Hypermethylated","Both"],index=["Total Genes","Total DMR"])
            ax = df.plot.bar(stacked=True,title="Common Statistics of %s" %i,color=["#428bca","#d9534f","#5cb85c"])
            plt.ylabel("Number of genes/DMR")
            plt.grid()
            plt.xticks([],horizontalalignment="center",rotation="horizontal")
            plt.table(cellText=[[str(lenhypogenes-lenbothgenes)+" ("+str(round(100*(lenhypogenes-lenbothgenes)/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",str(hypodmr)+" ("+str(round(100*hypodmr/(hypodmr+hyperdmr),2))+"%)"],\
                                [str(lenhypergenes-lenbothgenes)+" ("+str(round(100*(lenhypergenes-lenbothgenes)/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",str(hyperdmr)+" ("+str(round(100*hyperdmr/(hypodmr+hyperdmr),2))+"%)"],\
                                [str(lenbothgenes)+" ("+str(round(100*lenbothgenes/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)","X"],\
                                [str(lenhypogenes+lenhypergenes-lenbothgenes)+"/"+str(lengenesnum)+" ("+str(round(100*(lenhypogenes+lenhypergenes-lenbothgenes)/lengenesnum,2))+"%)",str(hypodmr+hyperdmr)+"/"+str(dmrnum)+" ("+str(round(100*(hypodmr+hyperdmr)/dmrnum,2))+"%)"]],\
                    colLabels=["Common Genes","Common DMR"],rowLabels=["Hypomethylated","Hypermethylated","Both","Total"],cellLoc="center",rowLoc="center")
            plt.savefig(os.path.join(os.path.join(commonpath,"%s" %comparelist),"%s-Common-Statistics.png" %i),bbox_inches="tight",dpi=800)
            plt.close()
        else :
            genevalues.append([lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes])
            genetablevalues.append([str(lenhypogenes-lenbothgenes)+" ("+str(round(100*0,2))+"%)",\
                                    str(lenhypergenes-lenbothgenes)+" ("+str(round(100*0,2))+"%)",\
                                    str(lenbothgenes)+" ("+str(round(100*0,2))+"%)",\
                                    str(lenhypogenes+lenhypergenes-lenbothgenes)+"/"+str(lengenesnum)+" ("+str(round(100*(lenhypogenes+lenhypergenes-lenbothgenes)/lengenesnum,2))+"%)"])
            dmrvalues.append([hypodmr,hyperdmr])
            dmrtablevalues.append([str(hypodmr)+" ("+str(round(100*0,2))+"%)",\
                                str(hyperdmr)+" ("+str(round(100*0,2))+"%)",\
                                str(hypodmr+hyperdmr)+"/"+str(dmrnum)+" ("+str(round(100*(hypodmr+hyperdmr)/dmrnum,2))+"%)"])
            df = pd.DataFrame([[lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes],[hypodmr,hyperdmr]],columns=["Hypomethylated","Hypermethylated","Both"],index=["Total Genes","Total DMR"])
            ax = df.plot.bar(stacked=True,title="Common Statistics of %s" %i,color=["#428bca","#d9534f","#5cb85c"])
            plt.ylabel("Number of genes/DMR")
            plt.grid()
            plt.xticks([],horizontalalignment="center",rotation="horizontal")
            plt.table(cellText=[[str(lenhypogenes-lenbothgenes)+" ("+str(round(100*0,2))+"%)",str(hypodmr)+" ("+str(round(100*0,2))+"%)"],\
                                [str(lenhypergenes-lenbothgenes)+" ("+str(round(100*0,2))+"%)",str(hyperdmr)+" ("+str(round(100*0,2))+"%)"],\
                                [str(lenbothgenes)+" ("+str(round(100*0,2))+"%)","X"],\
                                [str(lenhypogenes+lenhypergenes-lenbothgenes)+"/"+str(lengenesnum)+" ("+str(round(100*(lenhypogenes+lenhypergenes-lenbothgenes)/lengenesnum,2))+"%)",str(hypodmr+hyperdmr)+"/"+str(dmrnum)+" ("+str(round(100*(hypodmr+hyperdmr)/dmrnum,2))+"%)"]],\
                    colLabels=["Common Genes","Common DMR"],rowLabels=["Hypomethylated","Hypermethylated","Both","Total"],cellLoc="center",rowLoc="center")
            plt.savefig(os.path.join(os.path.join(commonpath,"%s" %comparelist),"%s-Common-Statistics.png" %i),bbox_inches="tight",dpi=800)
            plt.close()
    genetablevaluesbis = [[j[i] for j in genetablevalues] for i in range(4)]
    df = pd.DataFrame(genevalues,columns=["Hypomethylated","Hypermethylated","Both"],index=comparelist)
    ax = df.plot.bar(stacked=True,title="Common Statistics for Genes",color=["#428bca","#d9534f","#5cb85c"])
    plt.ylabel("Number of genes")
    plt.grid()
    plt.xticks([],horizontalalignment="center",fontsize="small")
    plt.table(cellText=genetablevaluesbis,colLabels=comparelist,rowLabels=["Hypomethylated","Hypermethylated","Both","Total"],cellLoc="center",rowLoc="center")
    plt.savefig(os.path.join(os.path.join(commonpath,"%s" %comparelist),"Common-Statistics-Genes.png"),bbox_inches="tight",dpi=800)
    plt.close()
    dmrtablevaluesbis =[[j[i] for j in dmrtablevalues] for i in range(3)]
    df = pd.DataFrame(dmrvalues,columns=["Hypomethylated","Hypermethylated"],index=comparelist)
    ax = df.plot.bar(stacked=True,title="Common Statistics for DMR",color=["#428bca","#d9534f"])
    plt.ylabel("Number of DMR")
    plt.grid()
    plt.xticks([],horizontalalignment="center",fontsize="small")
    plt.table(cellText=dmrtablevaluesbis,colLabels=comparelist,rowLabels=["Hypomethylated","Hypermethylated","Total"],cellLoc="center",rowLoc="center")
    plt.savefig(os.path.join(os.path.join(commonpath,"%s" %comparelist),"Common-Statistics-DMR.png"),bbox_inches="tight",dpi=800)
    plt.close()
    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-%s.csv " %comparelist),"w").writelines(header)
    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-%s.csv " %comparelist),"a").writelines(line for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:] if any(k in line.split(",")[17].split(";") for k in commongenes))
    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-List-%s.csv " %comparelist),"w").writelines("geneSymbol\n")
    # open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-List-%s.csv " %comparelist),"a").writelines([",".join([k[0],l,k[2],"\n"]) for k in patientgenes for l in k[1].split(";") if l in commongenes])
    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-List-%s.csv " %comparelist),"a").writelines(set([",".join([l,"\n"]) for k in patientgenes for l in k[1].split(";") if l in commongenes]))
    samemeth = input("Do you want to trim for gene with only the same methylation profile ? (Yes/No)")
    while samemeth != "No":
        if samemeth == "Yes":
            same = []
            both = []
            diff = []
            for i in set(commongenes):
                temp = set(j[2] for j in patientgenes if i in j[1].split(";"))
                tempbis = set(tuple(k) for k in [[j[0],j[2]] for j in patientgenes if i in j[1].split(";")])
                if len(temp) == 1:
                    same.append([i,list(temp)])
                elif len(temp) != 1:
                    if len(tempbis) != len(comparelist)*2:
                        diff.append(i)
                    elif len(tempbis) == len(comparelist)*2:
                        both.append(i)
            print("Total number of common genes for %s: %s" %(comparelist,len(set(commongenes))))
            print("Total number of common genes with same methylation profile for %s: %s" %(comparelist,len(same)+len(set(both))))
            print("Total number of common genes with different methylation profile for %s: %s" %(comparelist,len(set(diff))))
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Same-Meth-%s.csv " %comparelist),"w").writelines(header)
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Same-Meth-%s.csv " %comparelist),"a").writelines(line for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:] if any(k[0] in line.split(",")[17].split(";") for k in same))
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Same-Meth-List-%s.csv " %comparelist),"w").writelines("geneSymbol,methylation\n")
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Same-Meth-List-%s.csv " %comparelist),"a").writelines([",".join([k[0],k[1][0]])+"\n" for k in same])
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Both-Meth-%s.csv " %comparelist),"w").writelines(header)
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Both-Meth-%s.csv " %comparelist),"a").writelines(line for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:] if any(k in line.split(",")[17].split(";") for k in both))
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Both-Meth-List-%s.csv " %comparelist),"w").writelines("sample,geneSymbol,methylation\n")
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Both-Meth-List-%s.csv " %comparelist),"a").writelines([",".join(k)+"\n" for k in patientgenes if k[1] in both])
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Diff-Meth-%s.csv " %comparelist),"w").writelines(header)
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Diff-Meth-%s.csv " %comparelist),"a").writelines(line for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:] if any(k in line.split(",")[17].split(";") for k in diff))
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Diff-Meth-List-%s.csv " %comparelist),"w").writelines("sample,geneSymbol,methylation\n")
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Diff-Meth-List-%s.csv " %comparelist),"a").writelines([",".join(k)+"\n" for k in patientgenes if k[1] in diff])
            if len(comparelist) == 2:
                if input("Do you want to look for gene methylation more affected in a situation ? (Yes/No)") == "Yes":
                    print("Beware : analysis not available for genes hypomethylated and hypermethylated in both samples")
                    meth = [[line.split(",")[0],line.split(",")[17],line.split(",")[20],line.split(",")[21].strip("\n")] for line in open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Same-Meth-%s.csv " %comparelist),"r").readlines()[1:]]
                    a = comparelist[0]
                    b = comparelist[1]
                    asupb = 0
                    bsupa = 0
                    asupbgene = []
                    bsupagene = []
                    for i in same:
                        temp = []
                        for j in comparelist:
                            tempbis = []
                            for k in meth:
                                if i[0] in k[1].split(";") and k[0] == j:
                                    tempbis.append(float(k[2]))
                            temp.append([j,abs(statistics.mean(tempbis))])
                        if temp[0][0] == a:
                            if temp[0][1]-temp[1][1] > 0.2:
                                asupb += 1
                                asupbgene.append(i[0])
                            elif temp[1][1]-temp[0][1] > 0.2:
                                bsupa += 1
                                bsupagene.append(i[0])
                        elif temp[0][0] == b:
                            if temp[0][1]-temp[1][1] > 0.2:
                                bsupa += 1
                                bsupagene.append(i[0])
                            elif temp[1][1]-temp[0][1] > 0.2:
                                asupb += 1
                                asupbgene.append(i[0])
                    # asupbdmr = sum(1 for j in meth if any(k in j[1].split(";") for k in asupbgene))
                    # bsupadmr = sum(1 for j in meth if any(k in j[1].split(";") for k in bsupagene))
                    print("Total number of common genes with more impact for %s than %s: %s" %(a,b,asupb))
                    # print("Total number of common DMR with more impact for %s than %s: %s" %(a,b,asupbdmr))
                    print("Total number of common genes with more impact for %s than %s: %s" %(b,a,bsupa))
                    # print("Total number of common DMR with more impact for %s than %s: %s" %(b,a,bsupadmr))
                    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Prog-%s-Sup-%s.csv " %(a,b)),"w").writelines(header)
                    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Prog-%s-Sup-%s.csv " %(a,b)),"a").writelines(line for line in open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Same-Meth-%s.csv " %comparelist),"r").readlines()[1:] if any(k in line.split(",")[17].split(";") for k in asupbgene))
                    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Prog-%s-Sup-%s.csv " %(b,a)),"w").writelines(header)
                    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Prog-%s-Sup-%s.csv " %(b,a)),"a").writelines(line for line in open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMR-Common-Genes-Same-Meth-%s.csv " %comparelist),"r").readlines()[1:] if any(k in line.split(",")[17].split(";") for k in bsupagene))
            elif len(comparelist) != 2:
                print("Cannot check difference in methylation modification : you need to reduce the number of sample to compare to 2")
            break
        elif samemeth != "Yes" or "No":
            print("Error : wrong input!")
            samemeth = input("Do you want to trim for gene with only the same methylation profile ? (Yes/No)")

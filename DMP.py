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
    header[7] = "CT-AVG"
    header[8] = "Patient-AVG"
    header = ",".join(header)
    global concatfile
    concatfile = os.path.join(resultspath,"DMP-All.csv")   
    open(concatfile,"w").writelines(header)
    open(concatfile,"a").writelines([",".join([filein,",".join([i.strip("\n").strip("\"") for i in line.split(",")]),"hypomethylation\n"]) for filein in samplesfiles for line in open(os.path.join(samplespath,filein),"r") if not line.startswith("\"logFC\"") and float((",".join([i.strip("\n").strip("\"") for i in line.split("\t")])).split(",")[8]) < 0])
    open(concatfile,"a").writelines([",".join([filein,",".join([i.strip("\n").strip("\"") for i in line.split(",")]),"hypermethylation\n"]) for filein in samplesfiles for line in open(os.path.join(samplespath,filein),"r") if not line.startswith("\"logFC\"") and float((",".join([i.strip("\n").strip("\"") for i in line.split("\t")])).split(",")[8]) > 0])
    global patientlist
    patientlist = sorted(set([line.split(",")[0] for line in open(concatfile,"r").readlines()[1:]]))
    rmnogenes = input("Do you want to remove DMP without genes ? (Yes/No)")
    while rmnogenes != "No":
            if rmnogenes == "Yes" :
                open(os.path.join(resultspath,"DMP-All-NoEmptyGene.csv"), "w").writelines([line for line in open(concatfile, "r") if not line.split(",")[14] == ""])
                concatfile = os.path.join(resultspath,"DMP-All-NoEmptyGene.csv")
                break
            else :
                print("Error : wrong input!")
                rmnogenes = input("Do you want to remove DMP without genes ? (Yes/No)")
    for i in patientlist:
        open(os.path.join(resultspath,"%s.csv" %i),"w").writelines(header)
        linelist = [line for line in open(concatfile,"r").readlines()[1:] if line.split(",")[0] == i]
        linelist.sort(key=lambda x: x.split(",")[14])
        open(os.path.join(resultspath,"%s.csv" %i),"a").writelines(line for line in linelist)

def globalmeth():
    globalpath = os.path.join(resultspath,"Global")
    if not os.path.exists(globalpath):
        os.makedirs(globalpath)
    genevalues = []
    genetablevalues = []
    dmpvalues = []
    dmptablevalues =[]
    for i in patientlist:
        meth = [(line.split(",")[22].strip("\n"),line.split(",")[14]) for line in open(os.path.join(resultspath,"%s.csv" %i),"r").readlines()[1:]]
        dmpnum = len(meth)
        hypodmp = sum(1 for j in meth if j[0] == "hypomethylation")
        hyperdmp = sum(1 for j in meth if j[0] == "hypermethylation")        
        genesnum = {"hypomethylation":[],"hypermethylation":[]}
        for k in set(meth):
            genesnum[k[0]].append(k[1])
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
        dmpvalues.append([hypodmp,hyperdmp])
        dmptablevalues.append([str(hypodmp)+" ("+str(round(100*hypodmp/(hypodmp+hyperdmp),2))+"%)",\
                               str(hyperdmp)+" ("+str(round(100*hyperdmp/(hypodmp+hyperdmp),2))+"%)",\
                               str(hypodmp+hyperdmp)])
        # if lenbothgenes != 0:
        #     print("Beware : %s genes found both hypomethylated and hypermethylated for %s" %(lenbothgenes,i))
        #     open(os.path.join(globalpath,"%s-Global-Statistics-HypoHyper.csv" %i),"w").write(header)
        #     open(os.path.join(globalpath,"%s-Global-Statistics-HypoHyper.csv" %i),"a").writelines([line for line in open(os.path.join(resultspath,"%s.csv" %i),"r").readlines()[1:] if line.split(",")[14] in bothgenes])
        df = pd.DataFrame([[lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes],[hypodmp,hyperdmp]],columns=["Hypomethylated","Hypermethylated","Both"],index=["Genes","DMP"])
        ax = df.plot.bar(stacked=True,title="Global Statistics of %s" %i,color=["#428bca","#d9534f","#5cb85c"])
        plt.ylabel("Number of genes/DMP")
        plt.grid()
        plt.xticks([],horizontalalignment="center",rotation="horizontal")
        plt.table(cellText=[[str(lenhypogenes-lenbothgenes)+" ("+str(round(100*(lenhypogenes-lenbothgenes)/lengenesnum,2))+"%)",str(hypodmp)+" ("+str(round(100*hypodmp/dmpnum,2))+"%)"],\
                            [str(lenhypergenes-lenbothgenes)+" ("+str(round(100*(lenhypergenes-lenbothgenes)/lengenesnum,2))+"%)",str(hyperdmp)+" ("+str(round(100*hyperdmp/dmpnum,2))+"%)"],\
                            [str(lenbothgenes)+" ("+str(round(100*lenbothgenes/lengenesnum,2))+"%)","X"],\
                            [str(lengenesnum),dmpnum]],\
                  colLabels=["Genes","DMP"],rowLabels=["Hypomethylated","Hypermethylated","Both","Total"],cellLoc = "center", rowLoc = "center")
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
    dmptablevaluesbis =[[j[i] for j in dmptablevalues] for i in range(3)]    
    df = pd.DataFrame(dmpvalues,columns=["Hypomethylated","Hypermethylated"],index=list(patientlist))
    ax = df.plot.bar(stacked=True,title="Global Statistics for DMP",color=["#428bca","#d9534f"])
    plt.ylabel("Number of DMP")
    plt.grid()
    plt.xticks([],horizontalalignment="center",fontsize="small")
    plt.table(cellText=dmptablevaluesbis,colLabels=list(patientlist),rowLabels=["Hypomethylated","Hypermethylated","Total"],cellLoc="center",rowLoc="center")
    plt.savefig(os.path.join(globalpath,"Global-Statistics-DMP.png"),bbox_inches="tight",dpi=800)
    plt.close()

def common():
    comparelist = input("Which patient do you want to compare ? (Choose in this list %s and separate with comma)" %patientlist).split(",")
    commonpath = os.path.join(resultspath,"Common")
    if not os.path.exists(os.path.join(commonpath,"%s" %comparelist)):
        os.makedirs(os.path.join(commonpath,"%s" %comparelist))
    patientgenes = list(set(tuple(k) for k in [[line.split(",")[i].strip("\n") for i in [0,14,22]] for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:]]))
    genes = set([line.split(",")[14] for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:]])
    patientdict = defaultdict(list)
    for k in patientgenes:
        patientdict[k[0]].append(k[1])
    commongenes = []
    for i in genes:
        temp = [key for key,values in patientdict.items() if i in values]
        if len(set(temp)) == len(set(comparelist)):
            commongenes.append(i)
    genevalues = []
    genetablevalues = []
    dmpvalues = []
    dmptablevalues =[]
    for i in comparelist:
        meth = [(line.split(",")[22].strip("\n"),line.split(",")[14]) for line in open(os.path.join(resultspath,"%s.csv" %i),"r").readlines()[1:]]
        dmpnum = len(meth)
        hypodmp = sum(1 for j in meth if j[0] == "hypomethylation" and j[1] in commongenes)
        hyperdmp = sum(1 for j in meth  if j[0] == "hypermethylation" and j[1] in commongenes) 
        genesnum = {"hypomethylation":[],"hypermethylation":[]}
        for k in set(meth):
            genesnum[k[0]].append(k[1])
        bothgenes = set([j for j in genesnum.get("hypomethylation") if j in genesnum.get("hypermethylation") and j in commongenes])
        lenhypogenes = len(set([j for j in genesnum.get("hypomethylation") if j in commongenes]))
        lenhypergenes = len(set([j for j in genesnum.get("hypermethylation") if j in commongenes]))
        lenbothgenes = len(bothgenes)
        lengenesnum = len(set([j for key,values in genesnum.items() for j in values]))
        if lenhypogenes or lenhypergenes !=0:
            genevalues.append([lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes])
            genetablevalues.append([str(lenhypogenes-lenbothgenes)+" ("+str(round(100*(lenhypogenes-lenbothgenes)/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",\
                                    str(lenhypergenes-lenbothgenes)+" ("+str(round(100*(lenhypergenes-lenbothgenes)/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",\
                                    str(lenbothgenes)+" ("+str(round(100*lenbothgenes/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",\
                                    str(lenhypogenes+lenhypergenes-lenbothgenes)+"/"+str(lengenesnum)+" ("+str(round(100*(lenhypogenes+lenhypergenes-lenbothgenes)/lengenesnum,2))+"%)"])
            dmpvalues.append([hypodmp,hyperdmp])
            dmptablevalues.append([str(hypodmp)+" ("+str(round(100*hypodmp/(hypodmp+hyperdmp),2))+"%)",\
                                str(hyperdmp)+" ("+str(round(100*hyperdmp/(hypodmp+hyperdmp),2))+"%)",\
                                str(hypodmp+hyperdmp)+"/"+str(dmpnum)+" ("+str(round(100*(hypodmp+hyperdmp)/dmpnum,2))+"%)"])
            df = pd.DataFrame([[lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes],[hypodmp,hyperdmp]],columns=["Hypomethylated","Hypermethylated","Both"],index=["Total Genes","Total DMP"])
            ax = df.plot.bar(stacked=True,title="Common Statistics of %s" %i,color=["#428bca","#d9534f","#5cb85c"])
            plt.ylabel("Number of genes/DMP")
            plt.grid()
            plt.xticks([],horizontalalignment="center",rotation="horizontal")
            plt.table(cellText=[[str(lenhypogenes-lenbothgenes)+" ("+str(round(100*(lenhypogenes-lenbothgenes)/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",str(hypodmp)+" ("+str(round(100*hypodmp/(hypodmp+hyperdmp),2))+"%)"],\
                                [str(lenhypergenes-lenbothgenes)+" ("+str(round(100*(lenhypergenes-lenbothgenes)/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",str(hyperdmp)+" ("+str(round(100*hyperdmp/(hypodmp+hyperdmp),2))+"%)"],\
                                [str(lenbothgenes)+" ("+str(round(100*lenbothgenes/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)","X"],\
                                [str(lenhypogenes+lenhypergenes-lenbothgenes)+"/"+str(lengenesnum)+" ("+str(round(100*(lenhypogenes+lenhypergenes-lenbothgenes)/lengenesnum,2))+"%)",str(hypodmp+hyperdmp)+"/"+str(dmpnum)+" ("+str(round(100*(hypodmp+hyperdmp)/dmpnum,2))+"%)"]],\
                    colLabels=["Common Genes","Common DMP"],rowLabels=["Hypomethylated","Hypermethylated","Both","Total"],cellLoc="center",rowLoc="center")
            plt.savefig(os.path.join(os.path.join(commonpath,"%s" %comparelist),"%s-Common-Statistics.png" %i),bbox_inches="tight",dpi=800)
            plt.close()
        else:
            genevalues.append([lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes])
            genetablevalues.append([str(lenhypogenes-lenbothgenes)+" ("+str(round(100*0,2))+"%)",\
                                    str(lenhypergenes-lenbothgenes)+" ("+str(round(100*0,2))+"%)",\
                                    str(lenbothgenes)+" ("+str(round(100*lenbothgenes/(lenhypogenes+lenhypergenes-lenbothgenes),2))+"%)",\
                                    str(lenhypogenes+lenhypergenes-lenbothgenes)+"/"+str(lengenesnum)+" ("+str(round(100*(lenhypogenes+lenhypergenes-lenbothgenes)/lengenesnum,2))+"%)"])
            dmpvalues.append([hypodmp,hyperdmp])
            dmptablevalues.append([str(hypodmp)+" ("+str(round(100*0,2))+"%)",\
                                str(hyperdmp)+" ("+str(round(100*0,2))+"%)",\
                                str(hypodmp+hyperdmp)+"/"+str(dmpnum)+" ("+str(round(100*(hypodmp+hyperdmp)/dmpnum,2))+"%)"])
            df = pd.DataFrame([[lenhypogenes-lenbothgenes,lenhypergenes-lenbothgenes,lenbothgenes],[hypodmp,hyperdmp]],columns=["Hypomethylated","Hypermethylated","Both"],index=["Total Genes","Total DMP"])
            ax = df.plot.bar(stacked=True,title="Common Statistics of %s" %i,color=["#428bca","#d9534f","#5cb85c"])
            plt.ylabel("Number of genes/DMP")
            plt.grid()
            plt.xticks([],horizontalalignment="center",rotation="horizontal")
            plt.table(cellText=[[str(lenhypogenes-lenbothgenes)+" ("+str(round(100*0,2))+"%)",str(hypodmp)+" ("+str(round(100*0,2))+"%)"],\
                                [str(lenhypergenes-lenbothgenes)+" ("+str(round(100*0,2))+"%)",str(hyperdmp)+" ("+str(round(100*0,2))+"%)"],\
                                [str(lenbothgenes)+" ("+str(round(100*0,2))+"%)","X"],\
                                [str(lenhypogenes+lenhypergenes-lenbothgenes)+"/"+str(lengenesnum)+" ("+str(round(100*(lenhypogenes+lenhypergenes-lenbothgenes)/lengenesnum,2))+"%)",str(hypodmp+hyperdmp)+"/"+str(dmpnum)+" ("+str(round(100*(hypodmp+hyperdmp)/dmpnum,2))+"%)"]],\
                    colLabels=["Common Genes","Common DMP"],rowLabels=["Hypomethylated","Hypermethylated","Both","Total"],cellLoc="center",rowLoc="center")
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
    dmptablevaluesbis =[[j[i] for j in dmptablevalues] for i in range(3)]
    df = pd.DataFrame(dmpvalues,columns=["Hypomethylated","Hypermethylated"],index=comparelist)
    ax = df.plot.bar(stacked=True,title="Common Statistics for DMP",color=["#428bca","#d9534f"])
    plt.ylabel("Number of DMP")
    plt.grid()
    plt.xticks([],horizontalalignment="center",fontsize="small")
    plt.table(cellText=dmptablevaluesbis,colLabels=comparelist,rowLabels=["Hypomethylated","Hypermethylated","Total"],cellLoc="center",rowLoc="center")
    plt.savefig(os.path.join(os.path.join(commonpath,"%s" %comparelist),"Common-Statistics-DMP.png"),bbox_inches="tight",dpi=800)
    plt.close()
    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-%s.csv " %comparelist),"w").writelines(header)
    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-%s.csv " %comparelist),"a").writelines(line for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:] if line.split(",")[14] in commongenes)
    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-List-%s.csv " %comparelist),"w").writelines("geneSymbol\n")
    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-List-%s.csv " %comparelist),"a").writelines(k+"\n" for k in commongenes)
    print("Total number of common genes for %s: %s" %(comparelist,len(commongenes)))
    samemeth = input("Do you want to trim for gene with only the same methylation profile ? (Yes/No)")
    while samemeth != "No":
        if samemeth == "Yes":
            samegenes = []
            for i in comparelist:
                temp = [[j[1],j[2]] for j in patientgenes if j[0] == i and j[1] in commongenes]
                for j in commongenes:
                    tempbis = [k[1] for k in temp if k[0] == j]
                    if len(set(tempbis)) == 1:
                        samegenes.append(j)
            same = []
            diff = []
            for i in set(samegenes):
                temp = set([j[2] for j in patientgenes if j[1] == i])
                if len(temp) == 1:
                    same.append([i,list(temp)])
                if len(temp) != 1:
                    diff.append(i)
            print("Total number of common genes with same methylation profile for %s: %s" %(comparelist,len(same)))
            print("Total number of common genes with different methylation profile for %s: %s" %(comparelist,len(diff)))
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Same-Meth-%s.csv " %comparelist),"w").writelines(header)
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Same-Meth-%s.csv " %comparelist),"a").writelines(line for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:] for k in same if k[0] == line.split(",")[14])
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Same-Meth-List-%s.csv " %comparelist),"w").writelines("geneSymbol,methylation\n")
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Same-Meth-List-%s.csv " %comparelist),"a").writelines([",".join([k[0],k[1][0]])+"\n" for k in same])
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Diff-Meth-%s.csv " %comparelist),"w").writelines(header)
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Diff-Meth-%s.csv " %comparelist),"a").writelines(line for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:] if line.split(",")[14] in diff)
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Diff-Meth-List-%s.csv " %comparelist),"w").writelines("sample,geneSymbol,methylation\n")
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Diff-Meth-List-%s.csv " %comparelist),"a").writelines([",".join(k)+"\n" for k in patientgenes if k[1] in diff])
            if len(comparelist) == 2:
                if input("Do you want to look for gene methylation more affected in a situation ? (Yes/No)") == "Yes":
                    print("Beware : analysis only available for genes with the same methylation profile in both samples")
                    meth = [[line.split(",")[0],line.split(",")[14],line.split(",")[9],line.split(",")[22].strip("\n")] for line in open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Same-Meth-%s.csv " %comparelist),"r").readlines()[1:]]
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
                                if i[0] == k[1] and k[0] == j:
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
                    print("Total number of common genes with more impact for %s than %s: %s" %(a,b,asupb))
                    print("Total number of common genes with more impact for %s than %s: %s" %(b,a,bsupa))
                    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Prog-%s-Sup-%s.csv " %(a,b)),"w").writelines(header)
                    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Prog-%s-Sup-%s.csv " %(a,b)),"a").writelines(line for line in open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Same-Meth-%s.csv " %comparelist),"r").readlines()[1:] if line.split(",")[14] in asupbgene)
                    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Prog-%s-Sup-%s.csv " %(b,a)),"w").writelines(header)
                    open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Prog-%s-Sup-%s.csv " %(b,a)),"a").writelines(line for line in open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Same-Meth-%s.csv " %comparelist),"r").readlines()[1:] if line.split(",")[14] in bsupagene)
            elif len(comparelist) != 2:
                print("Cannot check difference in methylation modification : you need to reduce the number of sample to compare to 2")
            break
        elif samemeth != "Yes" or "No":
            print("Error : wrong input!")
            samemeth = input("Do you want to trim for gene with only the same methylation profile ? (Yes/No)")
    episign = input("Do you want to trim for gene with only the same probes shared between patients ? (Yes/No)")
    while episign != "No":
        if episign == "Yes":
            meth = [[line.split(",")[0],line.split(",")[14],line.split(",")[21],line.split(",")[22].strip("\n")] for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:] if line.split(",")[14] in commongenes]
            commonprobes = []
            commonprobesgenes = []
            for i in commongenes:
                temp = [j[2] for j in meth if j[1] == i]
                for k in set(temp):
                    tempbis = [j[0] for j in meth if j[2] == k]
                    if len(set(tempbis)) == len(comparelist):
                        commonprobes.append(k)
                        commonprobesgenes.append(i)
            commonprobes = set(commonprobes)
            commonprobesgenes = set(commonprobesgenes)
            print("Total number of common probes for %s: %s" %(comparelist,len(commonprobes)))
            print("Total number of genes with at least one probe in common for %s: %s" %(comparelist,len(commonprobesgenes)))
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Probes-%s.csv" %comparelist),"w").writelines(header)
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Probes-%s.csv" %comparelist),"a").writelines(line for j in comparelist for line in open(os.path.join(resultspath,"%s.csv" %j),"r").readlines()[1:] if line.split(",")[14] in commonprobesgenes)
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Probes-List-%s.csv" %comparelist),"w").writelines("geneSymbol\n")
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Probes-List-%s.csv" %comparelist),"a").writelines([k+"\n" for k in commonprobesgenes])
            samemethgenes = []
            for i in commonprobesgenes:
                temp = [j[3] for j in meth if j[1] == i]
                if len(set(temp)) == 1:
                    samemethgenes.append(i)
            print("Total number of common genes with same methylation profile for %s: %s" %(comparelist,len(set(samemethgenes))))
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Probes-Same-Meth-List-%s.csv" %comparelist),"w").writelines("geneSymbol\n")
            open(os.path.join(os.path.join(commonpath,"%s" %comparelist),"DMP-Common-Genes-Probes-Same-Meth-List-%s.csv" %comparelist),"a").writelines([k+"\n" for k in set(samemethgenes)])
            break
        elif episign != "Yes" or "No":
            print("Error : wrong input!")
            episign = input("Do you want to trim for gene with only exactly the same probes shared between patients ? (Yes/No)")

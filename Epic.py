#!/usr/local/bin/python3

import os

samples = input("Do you want to generate data for DMR or DMP ? (DMR/DMP/Exit)")
while samples != "Exit":
    if samples == "DMR":
        import DMR
        logFC = input("What logFC did you used ?")
        DMR.pathmaker("DMR/DMR_"+logFC)
        DMR.check()
        DMR.concat()
        analysis = input("What analysis do you want to do ? (Global/Common/Back)")
        while analysis != "Back":
            if analysis == "Global":
                DMR.globalmeth()
                analysis = input("What analysis do you want to do ? (Global/Common/Back)")
            elif analysis == "Common":
                DMR.common()
                analysis = input("What analysis do you want to do ? (Global/Common/Back)")
            else:
                print("Error : wrong input!")
                analysis = input("What analysis do you want to do ? (Global/Common/Back)")
        samples = input("Do you want to generate data for DMR or DMP ? (DMR/DMP/Exit)")
    elif samples == "DMP" :
        import DMP
        logFC = input("What logFC did you used ?")
        DMP.pathmaker("DMP/DMP_"+logFC)
        DMP.check()
        DMP.concat()
        analysis = input("What analysis do you want to do ? (Global/Common/Back)")
        while analysis != "Back":
            if analysis == "Global":
                DMP.globalmeth()
                analysis = input("What analysis do you want to do ? (Global/Common/Back)")
            elif analysis == "Common":
                DMP.common()
                analysis = input("What analysis do you want to do ? (Global/Common/Back)")
            else:
                print("Error : wrong input!")
                analysis = input("What analysis do you want to do ? (Global/Common/Back)")
        samples = input("Do you want to generate data for DMR or DMP ? (DMR/DMP/Exit)")
    else :
        print("Error : wrong input!")
        samples = input("Do you want to generate data for DMR or DMP ? (DMR/DMP/Exit)")

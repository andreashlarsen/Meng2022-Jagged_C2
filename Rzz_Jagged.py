#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Rzz against distance
Andreas Haahr Larsen
"""

# import libraries
import numpy as np
import pylab
import matplotlib as mpl
from matplotlib import pyplot as plt

#fontsize
fontsize = 12
plt.rcParams.update({'font.size': fontsize   })
    
# folder (2 = jagged1, 3 = jagged2)
Folder = 2

# plot a single or selected trajectories
SELECTED_TRAJS = 0
#trajs = [0,1,2,3,4,5,7,8,9,10,11,13,14,15,16,17]
trajs = [14,21]

# PROT: Jagged2, LIP: PCPS
#trajs = [0,3,4,5,6,7,12,16,17,20,23,24]: mode 1, dominating(!) mode
#trajs = [1,8,9,11,13,22]: mode 3->1, mode 3 is most frequent intermediate mode
#trajs = [2,19]: mode 2->1, mode 2 is less frequent intermediate mode
#trajs = [10,15,18]: mode 3->2->1
#trajs = [14]: mode 4->1, mode 4 very(!) intermediate but physical mode
#trajs = [21]: mode 4->2->1

# PROT: Jagged1, LIP: PCPS
#trajs = [5,8,11,22] # mode 1,clash with membrane 
#trajs = [6,17,18,19] # mode 2 *rep17: first encounter is contructive mode - very intermediate!, mode 2: clashes
#trajs = [2] # mode 3
#trajs = [3,7,9,15] # mode 5
#trajs = [0,4,12] # mode 1 and 3
#trajs = [10,14] # mode 1 and 2
#trajs = [13,24] # mode 3 and 4, mode 4: clashes (N-term in membrane)
#trajs = [16,23] # mode 4 and 5
#trajs = [1,20] # mode 1 and 4 *rep1: ends in mode 4
#trajs = [21] # mode 1 and 2 and end in mode 4 *rep21: first encounter is contructive mode, intermediate, ca frames 50-250

# PROT: Jagged1, LIP: PCGM1GM2
#trajs = [0,1,13,15,17] # mode 1, 
#trajs = [9,10] # mode 2, clash (glycos and DSL-EGFs)
#trajs = [3] # mode 3, clash (glycos and DSL-EGFs)
#trajs = [8,14] # mode 4, Bin with highest density: Rzz = -0.72, Dist = 3.39, clash (glycos and DSLEGFs)
#trajs = [5,7,16] # mode 5, Bin with highest density: Rzz = -0.15, Dist = 3.83, clash (DSL-EGF)
#trajs = [4] # mode 6, upright mode, Bin with highest density: Rzz = -0.29, Dist = 4.29
#trajs = [11] # mode 7, (semi-)upright mode, Bin with highest density: Rzz = -0.52, Dist = 4.07, clash (glycos)
# 15 repeats in total, 9 clashes, 5 in mode 1 (lying), 1 in mode 4 (upright)

# bins
bins = 100

# find fram corresponding to point on Dist vs. Rzz plot
aTol = 20 # increase tolerance area by a factor 
RzzInt = 0.00
DistInt = 5.00
RzzTol = 0.002*aTol
DistTol = 0.002*aTol

# include only last frames (out of 5734 = 2 us)
LAST = 0
last_frames = 500

# include only range of frames
RANGE = 0
from_frame = 2000
to_frame = 3000

# include only specific frame
SPECIFIC = 0
specific_frame = 5729

# choose to save figures or not
SAVE = 0
# save high resolution (and without title and grid)
HIGH_RES = 0

# plot data from all repeats
PLOT_ALL = 0

Reps = [0,25,0]

if Folder == 2:
    Protein = "Jagged1"
elif Folder == 3:
    Protein = "Jagged2apo"

Elements = [\
            ["PCGM1GM3","red",Reps[0]],\
            ["PCPS","blue",Reps[1]],\
            ["PCPSGM1","green",Reps[2]]\
            ] 

for Element in Elements:
    Lip = Element[0]
    Color = Element[1]
    Repeats = Element[2]
    if Repeats > 0:
        RzzAll = []
        DistAll = []
        print('\n##################################################\nLipid: %s' % Lip)
        if SELECTED_TRAJS == 0:
            trajs = np.arange(Repeats)
        for i in trajs:
            # print('i = %d' % i)
            # Rzz
            File = str(Folder) + "/" + Lip + "_" + str(i) + "/Rzz.xvg"
            print('File: %s' % File)
            Time,Rxx,Rxy,Rxz,Ryx,Ryy,Ryz,Rzx,Rzy,Rzz = np.genfromtxt(File,skip_header=32,usecols=[0,1,2,3,4,5,6,7,8,9],unpack=True)
            Time /= 1000 # convert from ps to ns
            
            if LAST:
                Time = Time[-last_frames:]
                Rzz = Rzz[-last_frames:]
            if RANGE:
                Time = Time[from_frame:to_frame]
                Rzz = Rzz[from_frame:to_frame]
            if SPECIFIC:
                Time = Time[specific_frame]
                Rzz = Rzz[specific_frame]
                print('Time = %s, Rzz = %s' % (Time, Rzz))  
                bins = 20
                
            if PLOT_ALL:
                Figure1=pylab.figure(1)
                plot1=pylab.subplot(111)
                plot1.plot(Time,Rzz,label="Rzz",color="r")
                
                plot1.legend(loc=1,fontsize=9)
                plot1.set_xlabel("Time [ns]")
                plot1.set_ylabel("Rzz [nm]")
                pylab.suptitle("Rzz")
                pylab.show()
            
            RzzAll = np.append(RzzAll,Rzz)
        
            # Dist
            File = str(Folder) + "/" + Lip + "_" + str(i) + "/dist_com_cent.xvg"
            Time,Dist = np.genfromtxt(File,skip_header=17,usecols=[0,1],unpack=True)
            
            if LAST:
                Time = Time[-last_frames:]
                Dist = Dist[-last_frames:]
            if RANGE:
                Time = Time[from_frame:to_frame]
                Dist = Dist[from_frame:to_frame]    
            if SPECIFIC:
                Time = Time[specific_frame]
                Dist = Dist[specific_frame]
                print('Time = %s, Rzz = %s' % (Time, Dist))

            # find structure of interest
            Rzz_filter_index = np.where((Rzz > RzzInt-RzzTol) & (Rzz < RzzInt + RzzTol))
            Dist_filter_index = np.where((Dist > DistInt-DistTol) & (Dist < DistInt + DistTol))
            Common_index = np.intersect1d(Rzz_filter_index,Dist_filter_index)
            l_Common = len(Common_index)
            
            if l_Common:
                for j in range(l_Common):
                    print('Structure of interest: lip %s, repeat %d, frame %d' % (Lip,i,Common_index[j]))
            # for debubbing:
            #else:    
                #len_Rzz = len(Rzz_filter_index[0])
                #len_Dist = len(Dist_filter_index[0])
                #print('left Rzz values = %d, left Dist values = %d, left common values = %d' % (len_Rzz,len_Dist,l_Common) )
            
            if PLOT_ALL:  
                Figure1=pylab.figure(1)
                plot1=pylab.subplot(111)
                plot1.plot(Time,Dist,label="Dist",color="g")
                
                plot1.legend(loc=1,fontsize=9)
                plot1.set_xlabel("Time [ns]")
                plot1.set_ylabel("Dist [nm]")
                pylab.suptitle("Dist")
                pylab.show()
                
            DistAll = np.append(DistAll,Dist)
            
            # Together
            if PLOT_ALL:
                Figure1=pylab.figure(1)
                plot1=pylab.subplot(111)
                plot1.plot(Rzz,Dist,label="Dist vs Rzz",color="k",marker=".",linestyle="none")
                
                plot1.legend(loc=1,fontsize=9)
                plot1.set_xlabel("Rzz")
                plot1.set_ylabel("Dist")
                pylab.suptitle("Dist vs Rzz")
                pylab.show()
            
        # Histogram Rzz vs Dist
        hist=plt.hist2d(RzzAll,DistAll,bins=bins,norm=mpl.colors.LogNorm(),cmap='viridis')
        plt.xlabel('$R_{zz}$')
        plt.ylabel('Distance [nm]')
        plt.colorbar()
        plt.ylim((2.5,7.5))
        plt.xlim((-1.0,1.0))
        if HIGH_RES == 0:
            #plt.title("C2 from " + Protein + " and " + Lip + " membrane")
            plt.minorticks_on()
            plt.grid()
            plt.grid(which='minor')
        if SAVE:
            if LAST:
                plt.savefig('Rzz_%s_%s_last.png' % (Protein,Lip))
            elif SPECIFIC:
                pass
            elif SELECTED_TRAJS:
                plt.savefig('Rzz_%s_%s_rep%s.png' % (Protein,Lip,trajs[0]))
            else:
                if HIGH_RES:
                    plt.savefig('Rzz_%s_%s.eps' % (Protein,Lip), format='eps')
                else:
                    plt.savefig('Rzz_%s_%s.png' % (Protein,Lip))
        plt.show()

        # find bin with highest occupancy
        maxbin=np.amax(hist[0])
        indc = np.where(maxbin==hist[0])
        if len(indc) > 1:
            indxRzz = indc[0][0]
            indxDist = indc[1][0]
        else:
            indxRzz = indc[0]
            indxDist = indc[1]
        Rzz_indx_min = hist[1][indxRzz]
        Rzz_indx_max = hist[1][indxRzz+1]
        Dist_indx_min = hist[2][indxDist]
        Dist_indx_max = hist[2][indxDist+1]
        Rzz_indx_mean = (Rzz_indx_min + Rzz_indx_max)/2
        Dist_indx_mean = (Dist_indx_min + Dist_indx_max)/2
        
        print('Bin with highest density: Rzz = %1.2f, Dist = %1.2f' % (Rzz_indx_mean,Dist_indx_mean))
        

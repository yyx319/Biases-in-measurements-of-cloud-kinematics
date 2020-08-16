# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 09:40:39 2019

@author: YYX
"""
########
#set up#
########
import numpy as np
# scale parameter den = den*sca 
true= 0      
sca = 10
thin2 = 0
# simulation parameter
beta_a = ['0.2','2','20']
time_a = ['10','30','60']

dat_dir = '/avatar/yuxuan/tracer_ls_rel'
directory = '/home/yuxuan/tracer_ls_rel/quant'
#setting LOS
prj_ori_a = ['x','y','z'] #'x','y','z'

sp1 = ['co','co_4_3','c18o','pnh3','hcn','c18o_4_3','n2h']
sp2 = ['CO','CO_4_3','C18O','pNH3','HCN','C18O_4_3','N2H']

if true == 0:
    second_m = np.zeros((3,3,3,7))
    for i, beta in enumerate(beta_a):
        for j,time in enumerate(time_a):
            for k,prj_ori in enumerate(prj_ori_a):
                for s in range(len(sp1)):
                    if sca!=1:     
                        m0 = np.loadtxt("%s/rescale_moment_%.1f/beta%s_%s_%s/m0_%.1f_%s.txt"%(dat_dir,sca,beta,time,prj_ori,sca,sp1[s]))
                        m2 = np.loadtxt("%s/rescale_moment_%.1f/beta%s_%s_%s/m2_%.1f_%s.txt"%(dat_dir,sca,beta,time,prj_ori,sca,sp1[s]))
                    if sca==1:
                        m0 = np.loadtxt("%s/moment/beta%s_%s_%s/m0_%s.txt"%(dat_dir,beta,time,prj_ori,sp1[s]))
                        m2 = np.loadtxt("%s/moment/beta%s_%s_%s/m2_%s.txt"%(dat_dir,beta,time,prj_ori,sp1[s]))
                    lw = np.average(m2, weights=m0)
                    print(r'\beta=%s, t=%s t_{ff} orientation=%s, sp=%s: lw=%.2f'%(beta, time, prj_ori, sp2[s],lw))
                    second_m[i,j,k,s]=lw  
    if sca!=1:
        np.save("%s/linewidth_rescale_%.1f"%(directory,sca), second_m)   
    if sca==1:
        np.save("%s/linewidth"%(directory), second_m) 
        
        
if true == 1:
    second_m = np.zeros((3,3,3))
    for i, beta in enumerate(beta_a):
        for j,time in enumerate(time_a):
            for k,prj_ori in enumerate(prj_ori_a):
                m0 = np.loadtxt("%s/moment/beta%s_%s_%s/m0_true.txt"%(dat_dir,beta,time,prj_ori))
                m2 = np.loadtxt("%s/moment/beta%s_%s_%s/m2_true.txt"%(dat_dir,beta,time,prj_ori))
                lw = np.average(m2, weights=m0)
                print(r'\beta=%s, t=%s t_{ff} orientation %s: lw=%.2f'%(beta, time, prj_ori, lw))
                second_m[i,j,k]=lw   
    np.save("%s/linewidth_true"%(directory), second_m) 

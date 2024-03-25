#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 10:23:58 2023

@author: school
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sympy import symbols, solve
import glob
plt.ion()

class Noise:
    
    '''This class contains functions to perform Noise analysis of data at one IF Frequency
    ****Currently there is only code for Shot Noise Model.
    
    parameters: 
        
       filename of the hot and cold data at the same frequency i.e
       
       Noise(hot_data,cold_data)
    
    functions: 
        
        shot() : To get the L, T_m from shot noise model, T_R, 
                 and fitting paramters of the IV and power curves '''
      
    def __init__(self,filename_hot,filename_cold):
        
        self.hot = pd.read_csv(filename_hot) # reading data
        
        self.cold = pd.read_csv(filename_cold)
        
        self.maxout=np.where(self.hot['Vs']==self.hot['Vs'].max())[0]
        
        if self.hot['Vs'].max()>30:  # TO remove values with sudden hike in Voltage
            self.hot=self.hot.drop(self.maxout)
            self.cold=self.cold.drop(self.maxout)
        
        self.T_h = self.hot.T3.mean()
        
        self.T_c = self.cold.T7.mean()
        
        self.Y = self.hot.IFPower/self.cold.IFPower #Y factor
        
        self.T_r = ((self.T_h-self.Y*self.T_c)/(self.Y-1))
        
        idxs=np.where(self.hot['Vs'].between(7.5,9))[0]
        #positive_tem = self.T_r[self.T_r>0]
        
        self.cad = self.Y[idxs].idxmax()  # finding the max Y value
        
        self.noise = {'T_r' : self.T_r[self.cad] }
        
        #restric tot 15.5 to 18.8 mV
        
    def shot(self):
        
        '''This Function returns the noise output a
        ccording to the shot noise model: 
            
            return: 
                list: [noise,a,b]
                
                noise: dictonary with the following keys
                
                {L= Noise characteris
                T_r= reciever Noise
                T_m= Mixer Noise}
                
                a = fitting parameter of the linear part of IV curve
                   [slope,intercept]
                b = fitting parameter of the linear part of Power curve
                   [slope,intercept]
                
        '''
        
        
        # Power and temperature values at the Volate where |T_r| is minimum
        P_h = self.hot.IFPower[self.cad] 
        
        P_c = self.cold.IFPower[self.cad]
        
        self.V_Tr = self.hot.Vs[self.cad] #Volate where |T_r| is minimum
        
        Tc = self.cold.T7[self.cad]
        
        Th = self.hot.T3[self.cad] 
        
        ids = np.where(self.hot['Vs'].between(16,26))[0]  #Linear range
        idsp = np.where(self.hot['Vs'].between(16,24))[0]
        # Linear part of the VI anf PV curves
        
        Vs = (self.hot['Vs'][ids[0]:ids[-1]]+self.cold['Vs'][ids[0]:ids[-1]])/2
        
        Is = (self.hot['Is'][ids[0]:ids[-1]]+self.cold['Is'][ids[0]:ids[-1]])/2
        
        IFPower = (self.hot['IFPower'][idsp[0]:idsp[-1]]+self.cold['IFPower'][idsp[0]:idsp[-1]])/2
        
        VsP = (self.hot['Vs'][idsp[0]:idsp[-1]]+self.cold['Vs'][idsp[0]:idsp[-1]])/2
        
        #performing the fit
        
        a = np.polyfit(Vs,Is,deg=1) #IV curve
        bo = np.polyfit(VsP,IFPower,deg=1) # Power curve
        
        self.VI = -a[1]/a[0] #(x intercept= -yint/slope)
        
        self.alpha = bo[0]/5.8  #(slope/5.8)
        
        self.T_IF = 5.8*(bo[1]/bo[0])+5.8*self.VI #(5.8(Vi-xintercept))

        ## Solving L and T_m 
        
        L, T_m = symbols('x y')
        
        eq1 = (P_h/self.alpha-self.T_IF-L*(T_m+Th)) 
        eq2 = (P_c/self.alpha-self.T_IF-L*(T_m+Tc))
        
        noisey=solve((eq1,eq2), (L, T_m))
        
        # Noise Output of the data according to Shot Noise Model
        
        noise={'T_r':self.T_r[self.cad], 
               'L': noisey[0][0], 
               'T_m': noisey[0][1],
               'T_IF': (self.T_r[self.cad]-noisey[0][1])/noisey[0][0]}
        
        return (noise,a,bo,self.cad) #[P_h,P_c,Th,Tc,V_Tr,alpha]
   
'''def LO(freq,sis):
    
    Finds the T'm and P_IF using LO drive
    
    paramters:
        freq: frequency to check at
        hot: hot run
        cold:cold run
        
    return:
        T'm= mixer noise
        PIF= IFPower
        
    hotd=sorted(glob.glob(f'sis{sis}*{freq}*IF6_hot_hires.txt'))
        
    coldd=sorted(glob.glob(f'sis{sis}*{freq}*IF6_cold_hires.txt'))
    
    cord=[]
    for i in range(len(hotd)):
    
        fit=Noise(hotd[i],coldd[i])
        
        idz=fit.cad
    
        Ph,th=fit.hot["IFPower"][idz],fit.hot.T3[idz]
    
        Pc,tc=fit.cold["IFPower"][idz],fit.cold.T7[idz]
        
        cord.append((Ph,th,Pc,tc))
    
    # Getting the curve and intersection
    x=np.linspace(-35,th+5,100)
    
    l=(cord[1][0]-cord[1][2])/(cord[1][1]-cord[1][3])
    
    m=(cord[0][0]-cord[0][2])/(cord[0][1]-cord[0][3])
    
    y3=cord[1][0]
    
    x3=cord[1][1]
    
    x1= cord[0][1]
    
    y1=cord[0][0]
    
    inx= ((m*x1-l*x3)-(y1-y3))/(m-l)
    
    iny=l*(inx-x3)+y3
    
    yb=l*(x-x3)+y3
    
    ya=m*(x-x1)+y1
    
    
    plt.plot(x,yb/10**-7,label='Upper')
    plt.plot(x,ya/10**-7,label='Lower')
    plt.axvline(cord[0][1])
    plt.axvline(cord[0][3])
    plt.ylabel('Power')
    plt.xlabel('Temperature')
    plt.title(f'{freq} GHz')
    plt.legend()
    plt.show()
    
    #point_of_intersection = inx, iny

    print(cord)
    print(inx,iny,'upper loss/lower loss:',m/l )
    
    return(-inx,iny)
    
        
'''        
    
    
        
        
def plot(Files_hot, Files_cold, param ='T_r'):
    
    ''' Function to plot noise characteristics
        
        paramters: 
            
            Files_hot: either one file or group of files for hot data
            
            Files_cold: either one file or group of files for corresponding cold data
        ** 
        If one file is specified:
            The IV curve, Power curve, Y factor and T_r are plotted w.r.t Vs
            
            Note:- Power is scaled to a factor of 10^7
            
        If a group of files is chosen:
            
            We plot the relationship between the IF Frequency in GHz 
            
            and one of the following paramters:
                
            T_r, T_m or L 
            
            which is specified via the param argument 
            ** default value is T_r 
            
            
            '''
        
            
    
    if type(Files_hot) == str and type(Files_cold) == str:
        
        file = Noise(Files_hot,Files_cold)
        plt.figure()
        
        #Power curves
        plt.plot(file.hot['Vs'],file.hot['IFPower']/10**-8,'r.',label ='hot Power')
    
        plt.plot(file.cold['Vs'],file.cold['IFPower']/10**-8,'b.',label ='cold Power')
        
        # IV curves
        
        plt.plot(file.hot['Vs'],file.hot['Is'],'m.',label='hot Current')
        
        plt.plot(file.cold['Vs'],file.cold['Is'],'c.',label='cold Current')
        
        #Yfactor and TR
        plt.plot(file.hot['Vs'],file.T_r,'g*',label='Reciever noise')
        
        plt.plot(file.hot['Vs'],file.Y,'y*',label='Yfactor')
        
        #linear fits
        
        a = file.shot()[1]
        b = file.shot()[2]
        
        #IV
        x = np.linspace(-a[1]/a[0],28,100)
        y = a[0]*x+a[1]
        plt.plot(x,y,'--',label='fit IV')
        
        #power
        x1 = np.linspace(-b[1]/b[0],28,100)
        y1 = (b[0]*x1+b[1])*10**8
        
        plt.plot(x1,y1,label='fit power')
        
        #Formatting
        plt.ylim(0,900)
        plt.xlim(-5,28)
        plt.legend(fontsize='6',loc='upper left')
        plt.grid()
        print(file.hot['Vs'][file.cad])
    
    
    else:
        
        # extracting the set of noise output for corresponding frequency
        
        freq=[]
        
        T_r=[]
        
        T_m=[]
        
        L=[]
        
        T_IF=[]
        
        alpha =[]
        
        Vopt=[]
        
        VI=[]
        
        for i in range(len(Files_hot)):
            
            freq.append(Files_hot[i].split('_')[5][:-3])
            
            nz=Noise(Files_hot[i],Files_cold[i])
            
            k=nz.shot()[0]
            
            T_r.append(k['T_r'])
            
            T_m.append(k['T_m'])
            
            L.append(k['L'])
            
            T_IF.append(k['T_IF'])
            
            VI.append(nz.VI)
            
            Vopt.append(nz.V_Tr)
            
            alpha.append(nz.alpha)
            
        d={'freq':freq,'T_r':T_r,'T_m':T_m,'L':L,'Vopt':Vopt,
           'alpha':alpha,'VI':VI,'T_IF':T_IF}
        
        t={'freq (GHz)':freq,
           'T_r (K)':np.around(np.array(T_r).astype(float),2),
           'T_m (K)':np.around(np.array(T_m).astype(float),2),
           'L (Db)':np.around(10*np.log10(np.array(L).astype(float)),2),
           'Vopt':np.around(np.array(Vopt).astype(float),2),
           'alpha e-8':np.around(np.array(alpha).astype(float)*10**8,2),
           'VI ':np.around(np.array(VI).astype(float),2)}
        
        table= pd.DataFrame(t)
        data = pd.DataFrame(d)
        
        data = data.sort_values('freq') # sorting values accoring to frequency
        
        data['freq'] = pd.to_numeric(data['freq']) #Turning Frequency Numeric 
        
        table.to_excel("sis12.xlsx")  
        
        plt.figure()
        
        if param =='L':
            
            plt.plot(data['freq'],10*np.log10(data['L'].to_numpy().astype(float)),
                     'o--',
                     label='L vs IF')
            plt.ylabel('L (db)')
            plt.xlabel('LO frequency(GHz)')
            plt.title('Conversion loss SIS 1')
            plt.legend()
            
            
        elif param =='T_m':
            
            plt.plot(data['freq'],data['T_m'],'*-', label='mixer noise vs IF')
            plt.ylabel('T_m')
            plt.xlabel('LO frequency(GHz)')
            plt.title('Mixer Noise SIS 1')
            plt.legend()
            
        
        elif param == 'T_IF':
             
             TIF=data['T_IF'].to_numpy() 
             TIF=np.array(TIF, dtype=float)
             plt.hist(TIF,bins=5)
             #plt.scatter(TIF,np.arange(len(TIF)))
             plt.axvline(nz.T_IF,color='r',label='expected value')
             plt.ylabel('freq')
             plt.xlabel('T_IF')
             plt.title(' T_IF variation SIS1')
             plt.legend()
             print(nz.T_IF)
        else:
            
            plt.plot(data['freq'],data['T_r'],'o-', label='Reciever noise vs IF')
            plt.ylabel('T_r (K)')
            plt.xlabel('IF frequency Hz')
            plt.title("Reciever Noise SIS 1")
            plt.legend()
            
        #plt.ylim(0,900)
        plt.grid()
        return data
    
            
    
        
'''        
    
hotf='sis1_sep12_2023_282GHz_16mW_IF6_hot_hires.txt'
coldf='sis1_sep12_2023_282GHz_16mW_IF6_cold_hires.txt'
plt.figure()
k=Plots(hotf,coldf)
g=pd.read_csv(hotf)
'''
file_h=sorted(glob.glob('sis2*hot_hires.txt'))
file_c=sorted(glob.glob('sis2*cold_hires.txt'))
a=11
plot(file_h[0:9],file_c[0:9],'T_r')


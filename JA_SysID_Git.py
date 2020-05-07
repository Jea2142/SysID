# -*- coding: utf-8 -*-
"""
Spyder Editor

This file currently contains frequency sweep flight inputs

@title: ASDFSDAFSDAFSDAFSDAFSDAFASD 
@author: Johnee Angarita
@year: 2020-Feb-21

@description: This is the integration file from John's dissertation. 
""" 
 
import numpy as np
#import matplotlib.pyplot as plt
import math
from scipy.signal import chirp
#from matplotlib.pyplot import figure

class JA_SysID():

    def __init__(self, flightFrequency):
        # Test Parameters
        self.totalTestT= 1000;
        self.commandRate = flightFrequency;

        #####################################################################
        # FREQUENCY SWEEP parameters (scaleChirp to adjust)
        # (i.e.) Frequencies, length and scale of movements
        self.delta_chirp = 45
        self.delta_trim  = 5
        self.scaleChirp  = 5
        self.f_i = 0.1 #Hz
        self.f_f = 3 #Hz


        # Define ending/final time of movements
        self.t_chirp1 = self.delta_chirp
        self.t_trim1  = self.t_chirp1+self.delta_trim


        #Create time arrays of movments to submit into chirp
        self.tc1 = np.linspace(0, self.t_chirp1, self.delta_chirp*self.commandRate)
        self.tt1 = np.linspace(self.t_chirp1, self.t_trim1, self.delta_trim*self.commandRate)


        #Remove time series overlap
        self.tt1 = self.tt1[1:]



        # chirp takes in frequencies as Hertz
        self.w1 = self.scaleChirp*chirp(self.tc1, f0=self.f_i, f1=self.f_f, t1=self.t_chirp1, method='linear', phi=-90)
       

        # Setting zero value array for length of prescribed vector
        self.z1 = np.zeros(len(self.tt1))


        self.t_sweep = np.concatenate((self.tc1, self.tt1), axis=0)
        self.f_sweep = np.concatenate((self.w1, self.z1), axis=0)
        self.noise = np.random.normal(0,0.5,len(self.f_sweep))
        self.f_sweep = self.f_sweep+self.noise

        # Above doesn't quite work because must import as velocity commands 
        # so needs adjusting. For now just use chirp code in next section 

        #######################################################################
        # %% New section from Hebert (HMPIM), scale factor to adjust G

        # Create time segments
        self.delta_HMPIM = 10 #Total chirp time in sec
        self.t_Hebert = np.linspace(0, self.delta_HMPIM, self.delta_HMPIM*self.commandRate)
        self.delta_htrim = 5
        self.t_mt = np.linspace(self.delta_HMPIM, self.delta_HMPIM+self.delta_htrim, self.delta_htrim*self.commandRate)

        # Trim zeros for velocity commands
        self.u_Htrim = np.zeros(len(self.t_mt))

        # Defining parameters
        self.a = 2*math.pi*(self.f_f-self.f_i)/self.delta_HMPIM
        self.b = 2*math.pi*self.f_i
        self.n = 1.0 #extinction factor
        self.G = 15 #HMPIM Gain (i.e the scale factor)

        self.Sigma1 = 0.5*self.a*self.t_Hebert**2
        self.Sigma2 = self.b+0.5*self.a
        self.Sigma2ext = self.b+0.5*self.a*self.t_Hebert**self.n


        # input series values
        self.X_chirp = np.sin(self.Sigma1)/self.Sigma2
        self.X_chirp_dot = (self.b + self.a*self.t_Hebert) * np.cos(self.Sigma1)/self.Sigma2
        
        self.X_HMPIM = self.G*np.sin(self.Sigma1)/self.Sigma2ext
        self.X_HMPIM_dot = self.G*((self.b+self.a*self.t_Hebert)*np.cos(self.Sigma1)/self.Sigma2ext - self.a*self.n*self.t_Hebert**(self.n-1)*np.sin(self.Sigma1)/2/(self.Sigma2ext**2))

        self.u_Chirp_dot = np.append(self.X_chirp_dot, self.u_Htrim, axis=0)
        self.u_HMPIM_dot = np.append(self.X_HMPIM_dot, self.u_Htrim, axis=0)


        # %% Multisine
        # Omega values in Hertz but equation takes that into acount
        # Phase values must be changed into rad/s before substituting into eqn
        # First take values from Alabasi
        self.delta_multi = 10
        self.t_m = np.linspace(0, self.delta_multi, self.delta_multi*self.commandRate)
        self.delta_mtrim = 5
        self.t_mt = np.linspace(self.delta_multi, self.delta_multi+self.delta_mtrim, self.delta_mtrim*self.commandRate)

        # Trim zeros for velocity commands
        self.du_trim = np.zeros(len(self.t_mt))


        self.Omega1 = [0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.6]
        self.Omega2 = [x+0.1 for x in self.Omega1]
        self.Omega3 = [x+0.1 for x in self.Omega2]
        self.Omega4 = [x+0.1 for x in self.Omega3]
        self.Phase1 = np.array([90.0, 261.9, 69.5, 60.0, 224.6, 229.7, 200.3])
        self.Phase2 = np.array([92.5, 268.5, 70.2, 61.1, 227.2, 221.7, 204.6])
        self.Phase3 = np.array([111.9, 301.6, 118.3, 141.1, 324.9, 328.1, 331.2])
        self.Phase4 = np.array([187.7, 76.8, 318.3, 28.8, 271.1, 340.4, 42.7])
        self.Phase1 = math.pi*self.Phase1/180
        self.u1  = 0
        self.du1 = 0
        self.A = 0.50

        for x in range(0,len(self.Omega1)):
            self.u1  = self.u1+self.A*np.sin(self.Omega1[x]*2*math.pi*self.t_m + self.Phase1[x])
            self.du1 = self.du1+(self.A*2*math.pi*self.Omega1[x])*np.cos(self.Omega1[x]*2*math.pi*self.t_m + self.Phase1[x])
            
        # Now add trim value
        self.du2=np.append(self.du1, self.du_trim, axis=0)

        # Additional code for reading MATLAB phase outputs
#        PhaseImport = np.loadtxt(fname = "C:/Users/johna/Documents/MATLAB/multisineOptimalPhase.txt")
#        print(PhaseImport)
#        print(len(PhaseImport))


        # %% Define which case of input excitaion to use
         # Case 0: Chirp
         # Case 1: HMPIM
         # Case 2: Multisine
         # Case 3: Combine the three systems singular test
         
            # Then number of repetitions

        self.caseNum = 0
        self.repNum = 0
 
        if self.caseNum == 0 :
            self.finalInput = self.u_Chirp_dot
            
            
            for x in range(0,self.repNum):
                self.finalInput = np.append(self.finalInput, self.u_Chirp_dot, axis=0)
          
        elif self.caseNum == 1:
            self.finalInput = self.u_HMPIM_dot


            for x in range(0,self.repNum):
                 self.finalInput = np.append(self.finalInput, self.u_HMPIM_dot, axis=0)
          
        elif self.caseNum == 2:
            self.finalInput = self.du2
         
            for x in range(0,self.repNum):
                self.finalInput = np.append(self.finalInput, self.du2, axis=0)
        
        elif self.caseNum == 3:
            self.intInput = np.append( self.u_HMPIM_dot, self.du2, axis=0)
            self.finalInput = np.append(self.f_sweep, self.intInput, axis=0)
            
            for x in range(0,self.repNum):
                self.finalInput = np.append(self.f_sweep, self.intInput, axis=0)
                
        else: print("Improper case")
         


        # Put in form of [velx, vely, velz, ang]
        # Set up zeros and oop to dictate which axis to excite
        # Set i for what specific to excite
        self.axis = 3
        self.velCommands = np.zeros((len(self.finalInput),4))

        for x in range(0,len(self.finalInput)):
            self.velCommands[x][self.axis] = self.finalInput[x]
        
    def getAmountOfSteps(self):
        return len(self.velCommands)
        
    def getVelVector4(self, timeStep):
        return self.velCommands[timeStep]
    
    
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 10:36:38 2021

@author: zheng
"""

import numpy as np
# from sklearn.neighbors import DistanceMetric
from sklearn.metrics import pairwise_distances, mean_squared_error
import copy
import matplotlib.pyplot as plt
# import time

class MultiScaleKernelRidge():
    
    def __init__(self, gamma=1, 
                 densityFactor=1, KNN=10, 
                 singleScale=0, D=None, startScale=0,
                 usePINV=0, exactPoints=[]):
        """ 
        Multi-scale kernel ridge regression
        
        Reference: Adaptation of multiscale function extension to inexact matching. N. Duchateau et al. Proc. GSI, 2013
        
        Author: Benoit Freiche (CREATIS Lyon, France), modified by Fei Zheng (CREATIS Lyon, France)
        
        Parameters:
            gamma: smaller -> smoother
            KNN: number of nearest neighbors used to estimate the density of samples
            densityFactor: factor multiplies density which decide the accomplishment of 
                multi-scale fitting and stop the iteration
            singleScale: 1 if use single scare, 0 multi-scale
            D: bandwidth**2 = D**2/2**s, bandwidth the smaller the 
            startScale: the s value in bandwidth computation,
                if singleScale=1, s=startScale
                else, s = startScale + [0, 1, 2, etc]
            exactPoints:
        """
        self.gamma = gamma
        self.densityFactor = densityFactor 
        self.KNN = KNN
        self.singleScale = singleScale
        self.exactPoints = exactPoints # list of samples where to have exact matching
        self.X  = None
        self.NumFeatureY = None
        self.D  = D
        self.startScale = startScale
        self.C = None
        self.mse = None # mean squared error
        self.usePINV = usePINV
            
    def density(X, KNN):
        """ estimate density """
        Dg = pairwise_distances(X, metric = 'nan_euclidean')  # distance matrix, distance between any element and nan = nan
        DgSorted = np.sort(Dg + np.diag([np.inf]*X.shape[0])) # diagnols set to inf before sorting
        DgKNN = DgSorted[:, : min(KNN, (~np.isnan(DgSorted[0])).sum()-1)] # keep only KNN columns
        density = np.nanmean(DgKNN)
        return density, Dg
    
    def fit(self, X, Y, Xtest=None, iterMax=20, mse=0, NYSTROM=0):
        self.X = X; self.NumFeatureY = Y.shape[1]
        density, Dg = MultiScaleKernelRidge.density(X, self.KNN)  # Estimate data density
        if self.D is None: 
            self.D = np.nanmax(Dg) # if D is not given
        if self.singleScale: 
            iterMax = 1 # singleScale has only 1 loop
        self.scales = np.arange(self.startScale, self.startScale+iterMax)
        F_s = np.zeros(Y.shape)
        if Xtest is not None:
            Ypredict = np.zeros((Xtest.shape[0], self.NumFeatureY))
        else:
            self.C = []
        for s in self.scales:
            e_s = self.D**2/(2**s)
            if e_s < (self.densityFactor*density)**2:
                break
            else:
                print('\nscale: {:.3f}'.format(s))
                Ge_s = np.exp( -Dg**2 / e_s ) # kernel matrix
                pointsInexact = np.array(range(len(Ge_s)))
                pointsInexact[self.exactPoints] = [];
                c, f_s = self.fit_Bermanis(Ge_s, Y - F_s, self.X, e_s, 
                                           self.gamma, pointsInexact, usePINV=self.usePINV)
                if Xtest is not None:
                    Ypredict += self.predict_Bermanis(Xtest, self.X, e_s, c, NYSTROM=NYSTROM)
                else:
                    self.C.append(c) 
                F_s += f_s
        if mse:
            self.mse = mean_squared_error(Y, F_s)
        if Xtest is not None:
            return Ypredict
        else:
            return
    
    def predict(self, Xtest, NYSTROM=0):
        Ypredict = np.zeros((Xtest.shape[0], self.NumFeatureY));
        for i in range(len(self.C)):
            s = self.scales[i]
            e_s = self.D**2/(2**s)
            f_star_s = self.predict_Bermanis(Xtest, self.X, e_s, self.C[i],
                                             NYSTROM=NYSTROM)
            Ypredict += f_star_s
        return Ypredict
    
    @staticmethod
    def fit_Bermanis(B_s, f, x_s, e_s, gamma, i_list, usePINV=1):
        [n,l] = np.shape(B_s);
        if( (n == l) & (l == len(i_list)) ):
            EYE = np.eye(n)
        else :
            EYE = np.zeros((n,l))
            for j in range(len(i_list)):
                EYE[i_list[j],j] = 1
        K = ( B_s + (1/gamma) * EYE )
        #step 3
        if usePINV:
            # slower version for small number of attributes: 
            if n == l:
                invK = np.linalg.inv(K)
            else :
                invK = np.linalg.pinv(K)
            c = MultiScaleKernelRidge.nandot(f.T, invK).T 
        else:
            # t_start = time.time()
            invK = np.linalg.inv(K)
            # print("inv takes: {:.2f}".format(time.time()-t_start))
            # t_start = time.time()
            c = np.dot(f.T,invK).T
            # print("dot takes: {:.2f}".format(time.time()-t_start))
        # print(K.shape)
        f_s = MultiScaleKernelRidge.nandot(c.T , B_s).T
        return c, f_s
    
    @staticmethod
    def predict_Bermanis(x_star, x_s, e_s, c, NYSTROM=0):
        #step 4
        # if usePINV:
        #     c = MultiScaleKernelRidge.nandot(f.T, invK).T
        # else:
        #     c = np.dot(f.T,invK).T
        tmp2 = pairwise_distances(x_star, x_s, metric = 'nan_euclidean');
        if NYSTROM:
            tmp2 = np.exp( -tmp2**2 / (2*e_s**2) ) 
        else:
            tmp2 = np.exp( -tmp2**2 / e_s )  
        f_star_s = MultiScaleKernelRidge.nandot(c.T , tmp2.T).T
        return f_star_s
    
    @staticmethod
    def nandot(A, B, scale=1):
        A_bis = copy.copy(A)
        # classic dot product with 0 remplacing nans
        A_bis[np.isnan(A_bis)] = 0
        res = np.dot(A_bis,B)
        if scale: 
            scale_mx = np.dot(np.ones(np.shape(A)),B)/np.dot(~np.isnan(A),B)
            res = res*scale_mx
        mask = np.sum(np.isnan(A),axis = 1).astype('float')
        mask[mask == len(A.T)] = np.nan
        mask[~np.isnan(mask)] = 0
        res += np.tile(mask,(len(B.T),1)).T
        return res
    

# # ============== Example ================
if __name__ == "__main__" :
    plt.close('all')
    
    # train
    N = 200
    x = np.linspace(0, 2*np.pi, N).reshape(-1,1) 
    
    ### Create a multi-frequency sinus to simulate different densities between samples
    y1 = np.multiply( np.sin(x) , x<np.pi)
    y2 = np.multiply( np.sin(5*x) , x>=np.pi )
    y = y1 + y2
    
    y_noise = y + 0.05 * np.random.randn(N,1)
    
    # test
    N_test = 200
    x_test = np.linspace(0, 2*np.pi, N_test).reshape(-1,1) 
    
    gamma = 1
    # ###### Singlescale
    msrkSing = MultiScaleKernelRidge(gamma=gamma, 
                                      densityFactor=1, KNN=10, 
                                      singleScale=1, D=1, startScale=0, usePINV=0,
                                      exactPoints=[])
    ###### Single scale
    msrkSing.fit(x, y_noise)
    y_predict = msrkSing.predict(x_test)
    plt.figure(figsize=(10,10))
    plt.plot(x,y,color = 'r',linestyle = '--',label="ground truth")
    plt.plot(x,y_noise,'.',label="training")
    plt.plot(x_test,y_predict,color='green',label="prediction")
    plt.legend()
    plt.title('SingleScale \n Scale = 0')
    plt.show()
    
    ###### Multi scale
    msrkMulti = MultiScaleKernelRidge(gamma=gamma, 
                                      densityFactor=1, KNN=4, 
                                      singleScale=0, startScale=0, usePINV=0,
                                      exactPoints=[])
    msrkMulti.fit(x, y_noise)
    y_predict = msrkMulti.predict(x_test)
    plt.figure(figsize=(10,10))
    plt.plot(x,y,color = 'r',linestyle = '--',label="ground truth")
    plt.plot(x,y_noise,'.',label="training")
    plt.plot(x_test,y_predict,color='green',label="prediction")
    plt.legend()
    plt.title('MultiScale \n startScale = 0')
    plt.show()






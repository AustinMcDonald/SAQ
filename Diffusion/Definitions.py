import numpy as np


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def moving_average(a, n) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


######################################## #################### #################### #################### 
####  Magboltz array indexing for the iput file 
###       [0    , 1    ,2    , 3   , 4     , 5     , 6   , 7    , 8   , 9    , 10  ,     11  , 12    ,  13 , 14 ]
###Vals = [XePer, ArPer, Temp, Pres, Efield, Zdrift, Zerr, Tdiff, Terr, Ldiff, Lerr, LdiffTPC,LerrTPC, Mele, Merr]
###                  11    , 12    ,    13   ,    14  ,  15 , 16 ]
###              LdiffTPC  ,LerrTPC, TdiffTCP, TerrTPC, Mele, Merr]
############################ #################### #################### #################### ####################
#################### 
####################  Below is the Magboltz importing function for the velocity and both diffusion in both set of units
####################





def MB_V(data,x):
    Vz = data[x][:,5]
    P  = data[x][:,3]/760
    E  = data[x][:,4]
    xe = str(data[x][0][0])
    ar = str(data[x][0][1])
    lab = xe+'%Xe '+ar+'%He'
    X = E/P
    Y = Vz
    Yer =Vz*data[x][:,6]/100
    drop = np.where(Y == 0)[0]
    X = np.delete(X,drop)
    Y = np.delete(Y,drop)
    Yer=np.delete(Yer,drop)
    SORT = X.argsort()
    X = X[SORT]
    Y = Y[SORT]
    Yer = Yer[SORT]
    
    Xnew = X#np.linspace(0, 300, 10000)
    Ynew = Y#np.interp(Xnew, X, Y)
    return Xnew,Ynew,Yer,lab

def MB_mu(data,x):
    Vz = data[x][:,5]*1e5
    P  = data[x][:,3]/760
    E  = data[x][:,4]
    xe = str(data[x][0][0])
    ar = str(data[x][0][1])
    lab = xe+'%Xe '+ar+'%He'
    X = E/P
    Y = (Vz/E)*P
    Yer =Y*data[x][:,6]/100
    drop = np.where(Y == 0)[0]
    X = np.delete(X,drop)
    Y = np.delete(Y,drop)
    Yer=np.delete(Yer,drop)
    SORT = X.argsort()
    X = X[SORT]
    Y = Y[SORT]
    Yer = Yer[SORT]
    
    Xnew = np.linspace(0, 300, 10000)
    Ynew = np.interp(Xnew, X, Y)
    return Xnew,Ynew,Yer,lab

def MB_DLtpc(data,x):
    Vz = data[x][:,11]
    P  = data[x][:,3]/760
    E  = data[x][:,4]
    xe = str(data[x][0][0])
    ar = str(data[x][0][1])
    lab = xe+'%Xe '+ar+'%He'
    X = E/P
    Y = Vz*np.sqrt(P)
    Yer =Vz*data[x][:,12]/100
    drop = np.where(Y == 0)[0]
    X = np.delete(X,drop)
    Y = np.delete(Y,drop)
    Yer=np.delete(Yer,drop)
    SORT = X.argsort()
    X = X[SORT]
    Y = Y[SORT]
    Yer = Yer[SORT]

    Xnew = X#np.linspace(0, 300, 1000)
    Ynew = Y#np.interp(Xnew, X, Y)
    return Xnew,Ynew,Yer,lab

def MB_DTtpc(data,x):
    Vz = data[x][:,13]
    P  = data[x][:,3]/760
    E  = data[x][:,4]
    xe = str(data[x][0][0])
    ar = str(data[x][0][1])
    lab = xe+'%Xe '+ar+'%He'
    X = E/P
    Y = Vz*np.sqrt(P)
    Yer =Vz*data[x][:,14]/100
    drop = np.where(Y == 0)[0]
    X = np.delete(X,drop)
    Y = np.delete(Y,drop)
    Yer=np.delete(Yer,drop)
    SORT = X.argsort()
    X = X[SORT]
    Y = Y[SORT]
    Yer = Yer[SORT]

    Xnew = X#np.linspace(0, 300, 1000)
    Ynew = Y#np.interp(Xnew, X, Y)
    return Xnew,Ynew,Yer,lab


def MB_DL(data,x):
    Vz = data[x][:,9]
    P  = data[x][:,3]/760
    E  = data[x][:,4]
    xe = str(data[x][0][0])
    ar = str(data[x][0][1])
    lab = xe+'%Xe '+ar+'%He'
    X = E/P
    Y = Vz*P
    Yer =Vz*data[x][:,10]/100
    drop = np.where(Y == 0)[0]
    X = np.delete(X,drop)
    Y = np.delete(Y,drop)
    Yer=np.delete(Yer,drop)
    SORT = X.argsort()
    X = X[SORT]
    Y = Y[SORT]
    Yer = Yer[SORT]

    Xnew = X#np.linspace(0, 300, 1000)
    Ynew = Y#np.interp(Xnew, X, Y)
    return Xnew,Ynew,Yer,lab

def MB_DT(data,x):
    Vz = data[x][:,7]
    P  = data[x][:,3]/760
    E  = data[x][:,4]
    xe = str(data[x][0][0])
    ar = str(data[x][0][1])
    lab = xe+'%Xe '+ar+'%He'
    X = E/P
    Y = Vz*P
    Yer =Vz*data[x][:,8]/100
    drop = np.where(Y == 0)[0]
    X = np.delete(X,drop)
    Y = np.delete(Y,drop)
    Yer=np.delete(Yer,drop)
    SORT = X.argsort()
    X = X[SORT]
    Y = Y[SORT]
    Yer = Yer[SORT]

    Xnew = X#np.linspace(0, 300, 1000)
    Ynew = Y#np.interp(Xnew, X, Y)
    return Xnew,Ynew,Yer,lab


def MB_eff_EleL(data,x):
    Vz = data[x][:,5]
    Dl = data[x][:,9]
    P  = data[x][:,3]/760
    E  = data[x][:,4]
    xe = str(data[x][0][0])
    ar = str(data[x][0][1])
    lab = xe+'%Xe '+ar+'%He'
    X = E/P
    mu = (Vz*1e5)/E
    Y = Dl/mu
    Yer =Vz*data[x][:,6]/100
    drop = np.where(Y == 0)[0]
    X = np.delete(X,drop)
    Y = np.delete(Y,drop)
    Yer=np.delete(Yer,drop)
    SORT = X.argsort()
    X = X[SORT]
    Y = Y[SORT]
    Yer = Yer[SORT]

    Xnew = X#np.linspace(0, 300, 1000)
    Ynew = Y#np.interp(Xnew, X, Y)
    return Xnew,Ynew,Yer,lab
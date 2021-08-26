#autoslit precession routine

import numpy as np
from astropy.table import Table,Column
from astropy.io import ascii

#Taken from Autoslit and adapted for python
def precession(ra,dec,EPOCH,NEW_EPOCH):
    PI     = 3.14159265358979324
    RPD = PI/180.
    DPR = 1/RPD
    HALFPI = PI/2.
    TWOPI = 2.*PI
    DDEC = 0.0
    DRA = 0.0

    RA = ra*RPD
    DEC = dec*RPD
    DYEAR = NEW_EPOCH - EPOCH
    DELTADEC = DYEAR * DDEC
    DELTADEC = DELTADEC * RPD / 3600.
    DEC = DEC + DELTADEC

  #     Proper motion: RA is in radians; DRA is in sec of time/year.
    DELTARA = DYEAR * DRA
    DELTARA = DELTARA * 15. * RPD / 3600.
    RA = RA + DELTARA
  #     Calculate time parameters
    DT  = (EPOCH - NEW_EPOCH) * .001
    TAU = (EPOCH) * 0.001 - 1.900

    Z0 = (0.1117133 + 6.77527E-4 * TAU) * DT +     (1.4656E-4 - 1.31E-6    * TAU) * DT*DT + 8.726E-5 * DT**3

    Z = Z0 + (3.843118E-4 + 3.1998E-6 * TAU)*DT*DT + 1.5514E-6*DT**3

    ZJ = (0.097189874 - 4.13692E-4*TAU - 1.79381E-6*TAU*TAU) * DT -   (2.0687E-4 + 1.79381E-6*TAU) * DT*DT - 2.02652E-4 * DT**3

    Z = HALFPI + Z
    Z0 = Z0 - HALFPI
    ZJ = -ZJ
    XOLD = np.zeros(3)
    X1 = np.zeros(3)
    X2 = np.zeros(3)
  #     Precess coordinates using ROTATE().
    XOLD[0] = np.cos(RA) * np.cos(DEC)
    XOLD[1] = np.sin(RA) * np.cos(DEC)
    XOLD[2] = np.sin(DEC)
    X2 = ROTATE(2,Z,XOLD)
    X1 = ROTATE(0,ZJ,X2)
    X2 = ROTATE(2,Z0,X1)

    NEW_DEC = np.arcsin(X2[2])
    NEW_RA  = np.arccos(X2[0]/np.cos(NEW_DEC))

    if (X2[1] < 0.0):
        NEW_RA = TWOPI - NEW_RA
    RA = NEW_RA*DPR
    DEC = NEW_DEC*DPR

    return([RA,DEC])
#=================================================================================================================================
def ROTATE(M,ALPHA,XOLD):
    S = np.zeros(shape=(3,3))

    for i in range(0,3):
         for j in range(0,3):
            if ((i != M) and (j != M)):
               if (i != j):
                  if (i > j):
                      S[i,j] = (-1.)*np.sin(ALPHA)
                  if (i < j):
                      S[i,j] = np.sin(ALPHA)
               else:
                  S[i,j] = np.cos(ALPHA)
            else:
               S[i,j] = 0.0
               if (i == j):
                   S[i,j]=1.0

#     Rotate vector.
    XNEW = TRNSFM(S,XOLD)
    return(XNEW)

def TRNSFM(S,XOLD):
  XNEW = np.zeros(3)
  for i in range(0,3):
     XNEW[i]=0.0


#     Multiply matrix by vector.
  for i in range(0,3):
     for j in range(0,3):
        XNEW[i] = XNEW[i] + S[i,j]*XOLD[j]
  return(XNEW)

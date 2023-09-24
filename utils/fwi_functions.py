#This is a collection of functions including functions for the Canadian Fire Weather Index
#as well as 5 equations from (Lawson et al 1997) necessary to convert the Duff Moisture Code to
#an actual moisture percent value.
#
#FWI Functions:
#
#FFMC - takes temperature, relative humidity, wind speed, rain, and a previous FFMC value to produce the current FFMC value
#      - FFMC(17,42,25,0,85) = 87.692980092774448
#
#DMC - takes temperature, relative humidity, rainfall, previous DMC value, latitude, and current month to produce the current DMC value
#    - DMC(17,42,0,6,45.98,4) = 8.5450511359999997
#
#DC - takes temperature, rain, the previous DC value, latititude, and current month to produce the current DC value
#   - DC(17,0,15,45.98,4) = 19.013999999999999
#
#ISI - takes the wind speed and current FFMC value to produce the current ISI value
#    - ISI(25,87.692980092774448) = 10.853661073655068
#
#BUI - takes the current DMC and DC values to produce the current BUI value
#    - BUI(8.5450511359999997,19.013999999999999) = 8.4904265358371838
#
#FWI - takes the current ISI and BUI values to produce the current BUI value
#    - FWI(10.853661073655068,8.4904265358371838) = 10.096371392382368
#
#calcFWI - this function returns the current FWI value given all of the input values: 
#              month, temperature, relative humidity, wind speed, rain, previous FFMC, DMC, DC, and latitude
#        - calcFWI(4,17,42,25,0,85,6,15,45.98) = 10.096371392382368
#
# 
#
#Lawson equations:
#
#All of these equations take the current DMC and DC values and return moisture content as a % value
#
#LawsonEq1 - DMC National Standard and Coastal B.C. CWH (2.5-4 cm)^2
#          - LawsonEq1(8.5450511359999997)  = 250.7553985454235
#
#LawsonEq2 - Southern Interior B.C.3 (2-4 cm)^2
#          - LawsonEq2(8.5450511359999997)  = 194.93023948344205
#
#LawsonEq3 - Southern Yukon - Pine/White Spruce
#                             Feather moss, Sphagnum and Undifferentiated duff (2-4 cm)^2
#          - LawsonEq3(8.5450511359999997)  = 442.82109267231488
#
#LawsonEq4 - Southern Yukon - Pine/White Spruce
#                             Reindeer lichen (2-4 cm)^2
#          - LawsonEq4(8.5450511359999997)  = 746.02210402093272
#
#LawsonEq5 - Southern Yukon - White Spruce
#                             White spruce/feather moss (2-4 cm)^2
#          - LawsonEq5(8.5450511359999997)  = 853.2397847094652


import math
import numpy as np
#import xarray as xr
from numba import jit

class InvalidLatitude(Exception):
    """Exception to handle variables not covered by DataDict"""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value) + " is not a valid Latitude."

def FFMC(TEMP,RH,WIND,RAIN,FFMCPrev):
    '''Calculates today's Fine Fuel Moisture Code
    Parameters:
    TEMP is the 12:00 LST temperature in degrees celsius
    RH is the 12:00 LST relative humidity in %
    WIND is the 12:00 LST wind speed in kph
    RAIN is the 24-hour accumulated rainfall in mm, calculated at 12:00 LST
    FFMCPrev is the previous day's FFMC

    FFMC(17,42,25,0,85) = 87.692980092774448'''
    
    Nx = TEMP.shape[0]
    Ny = TEMP.shape[1]
    
    rf = np.zeros((Nx,Ny))*np.nan
    mr = np.zeros((Nx,Ny))*np.nan
    m = np.zeros((Nx,Ny))*np.nan
    ko = np.zeros((Nx,Ny))*np.nan
    kd = np.zeros((Nx,Ny))*np.nan
    k1 = np.zeros((Nx,Ny))*np.nan
    ew = np.zeros((Nx,Ny))*np.nan
    kw = np.zeros((Nx,Ny))*np.nan

    RH[np.where(RH>100.0)] = 100.0
    mo = 147.2 * (101.0 - FFMCPrev) / (59.5 + FFMCPrev)
    c_1 = np.where(RAIN>.5)
    rf[c_1] = RAIN[c_1]-.5
    c_2 = np.where(np.logical_and(RAIN>.5,mo<=150.))
    mr[c_2] = mo[c_2]+42.5 * rf[c_2] * np.exp(-100.0 / (251.0 - mo[c_2])) * (1.0-np.exp(-6.93 / rf[c_2]))
    c_3 = np.where(np.logical_and(RAIN>.5,mo>150.))
    mr[c_3] = mo[c_3]+42.5 * rf[c_3] * np.exp(-100.0 / (251.0 - mo[c_3])) * (1.0-np.exp(-6.93 / rf[c_3]))
    c_4 = np.where(np.logical_and(RAIN>.5,mr>250.))
    mr[c_4] = 250.0
    mo[c_1] = mr[c_1]
    ed = 0.942 * RH**0.679 + 11.0 * np.exp((RH - 100.0) / 10.0) + 0.18 * (21.1 - TEMP) * (1.0 - np.exp(-0.115 * RH))
    c_5 = np.where(mo>ed)
    ko[c_5] = 0.424 * (1.0 - (RH[c_5]/100.0)**1.7) + 0.0694 * WIND[c_5]**.5 * (1.0 - (RH[c_5] / 100.0)**8)
    kd[c_5] = ko[c_5] * 0.581 * np.exp(0.0365 * TEMP[c_5])
    m[c_5] = ed[c_5] + (mo[c_5] - ed[c_5]) * 10.0**(-kd[c_5])
    c_6 = np.where(mo<=ed)
    ew[c_6] = 0.618 * RH[c_6]**0.753 + 10.0 * np.exp((RH[c_6] - 100.0) / 10.0) + 0.18 * (21.1 - TEMP[c_6]) * (1.0 - np.exp(-0.115 * RH[c_6]))
    c_7 = np.where(np.logical_and(mo<=ed, mo < ew))
    k1[c_7] = 0.424 * (1.0 - ((100.0 - RH[c_7]) / 100.0)** 1.7) + 0.0694 * WIND[c_7]**.5 * (1.0 - (100.0 - RH[c_7]) / 100.0)** 8
    kw[c_7] = k1[c_7] * 0.581 * np.exp(0.0365 * TEMP[c_7])
    m[c_7] = ew[c_7] - (ew[c_7] - mo[c_7]) * 10.0**(-kw[c_7])
    c_8 = np.where(np.logical_and(mo<=ed, mo >= ew))
    m[c_8] = mo[c_8]
    return 59.5 * (250.0 - m) / (147.2 + m)


def DMC(TEMP,RH,RAIN,DMCPrev,LAT,MONTH):
    '''Calculates today's Duff Moisture Code
Parameters:
    TEMP is the 12:00 LST temperature in degrees celsius
    RH is the 12:00 LST relative humidity in %
    RAIN is the 24-hour accumulated rainfall in mm, calculated at 12:00 LST
    DMCPrev is the prevvious day's DMC
    Lat is the latitude in decimal degrees of the location for which calculations are being made
    Month is the month of Year (1..12) for the current day's calculations.

    DMC(17,42,0,6,45.98,4) = 8.5450511359999997'''
    import numpy as np
    Nx = TEMP.shape[0]
    Ny = TEMP.shape[1]
    
    re = np.zeros((Nx,Ny))*np.nan
    mo = np.zeros((Nx,Ny))*np.nan
    b = np.zeros((Nx,Ny))*np.nan
    k = np.zeros((Nx,Ny))*np.nan
    d1 = np.zeros((Nx,Ny))*np.nan
    mr = np.zeros((Nx,Ny))*np.nan
    pr = np.zeros((Nx,Ny))*np.nan
    
    c_0 = np.where(RH>100.)
    c_1 = np.where(RAIN>1.5)
    c_2 = np.where((RAIN>1.5) & (DMCPrev<=33.))
    c_3 = np.where((RAIN>1.5) & (DMCPrev>33.) & (DMCPrev<=65.))
    c_4 = np.where((RAIN>1.5) & (DMCPrev>33.) & (DMCPrev>65.))
    c_5 = np.where((RAIN>1.5) & (pr>0.))
    c_6 = np.where((RAIN>1.5) & (pr<=0.))
    c_7 = np.where((pr>0.))
    c_8 = np.where((pr<=0.))
    
    RH[c_0] = 100.
    re[c_1] = 0.92 * RAIN[c_1] - 1.27
    mo[c_1] = 20.0 + np.exp(5.6348 - DMCPrev[c_1] / 43.43)
    b[c_2] = 100.0 / (0.5 + 0.3 * DMCPrev[c_2])
    b[c_3] = 14.0 - 1.3 * np.log(DMCPrev[c_3])
    b[c_4] = 6.2 * np.log(DMCPrev[c_4]) - 17.2
    mr[c_1] = mo[c_1] + 1000.0 * re[c_1] / (48.77 + b[c_1] * re[c_1])
    pr[c_1] = 244.72 - 43.43 * np.log(mr[c_1] - 20.0)
    DMCPrev[c_7] = pr[c_7]
    DMCPrev[c_8] = 0
    if type(LAT)==np.ndarray:
        for i in range(Nx):
            for j in range(Ny):
                if TEMP[i,j] > -1.1:
                    d1[i,j] = DayLength(LAT[i,j],MONTH)
                    k[i,j] = 1.894 * (TEMP[i,j] + 1.1) * (100.0 - RH[i,j]) * d1[i,j] * 0.000001
                else:
                    k[i,j] = 0.0
    else:
        for i in range(Nx):
            for j in range(Ny):
                if TEMP[i,j] > -1.1:
                    d1[i,j] = DayLength(LAT,MONTH)
                    k[i,j] = 1.894 * (TEMP[i,j] + 1.1) * (100.0 - RH[i,j]) * d1[i,j] * 0.000001
                else:
                    k[i,j] = 0.0    
    return DMCPrev + 100.0 * k


def DC(TEMP,RAIN,DCPrev,LAT,MONTH):
    '''
    Calculates today's Drought Code
    Parameters:
    TEMP is the 12:00 LST temperature in degrees celsius
    RAIN is the 24-hour accumulated rainfall in mm, calculated at 12:00 LST
    DCPrev is the previous day's DC
    LAT is the latitude in decimal degrees of the location for which calculations are being made
    MONTH is the month of Year (1..12) for the current day's calculations.

    DC(17,0,15,45.98,4) = 19.013999999999999
    '''
    Nx = TEMP.shape[0]
    Ny = TEMP.shape[1]
    
    rd = np.zeros((Nx,Ny))*np.nan
    Qo = np.zeros((Nx,Ny))*np.nan
    Qr = np.zeros((Nx,Ny))*np.nan
    Dr = np.zeros((Nx,Ny))*np.nan
    Lf = np.zeros((Nx,Ny))*np.nan
    V = np.zeros((Nx,Ny))*np.nan    
    c_1 = np.where(RAIN>2.8)
    c_2 = np.where((RAIN>2.8) & (Dr>0.))
    c_3 = np.where((RAIN>2.8) & (Dr<=0.))

    rd[c_1] = 0.83 * RAIN[c_1] - 1.27
    Qo[c_1] = 800.0 * np.exp(-DCPrev[c_1] / 400.0)
    Qr[c_1] = Qo[c_1] + 3.937 * rd[c_1]
    Dr[c_1] = 400.0 * np.log(800.0 / Qr[c_1])
    DCPrev[c_2] = Dr[c_2]
    DCPrev[c_3] = 0.0
    
    if type(LAT)==np.ndarray:
        for i in range(Nx):
            for j in range(Ny):
                Lf[i,j] = DryingFactor(LAT[i,j], MONTH-1)
                if TEMP[i,j] > -2.8:
                    V[i,j] = 0.36 * (TEMP[i,j]+2.8) + Lf[i,j]
                else:
                    V[i,j] = Lf[i,j]
    else:
        for i in range(Nx):
            for j in range(Ny):
                Lf[i,j] = DryingFactor(LAT, MONTH-1)
                if TEMP[i,j] > -2.8:
                    V[i,j] = 0.36 * (TEMP[i,j]+2.8) + Lf[i,j]
                else:
                    V[i,j] = Lf[i,j]
                    
    c_4 = np.where(V < 0.)
    
    V[c_4] = 0.0

    D = DCPrev + 0.5 * V

    return D



def ISI(WIND,FFMC):
    '''Calculates today's Initial Spread Index
Parameters:
    WIND is the 12:00 LST wind speed in kph
    FFMC is the current day's FFMC

    ISI(25,87.692980092774448) = 10.853661073655068'''

    fWIND = np.exp(0.05039 * WIND)
    m = 147.2 * (101.0 - FFMC) / (59.5 + FFMC)
    fF = 91.9 * np.exp(-0.1386 * m) * (1.0 + m**5.31 / 49300000.0)
    return 0.208 * fWIND * fF



def BUI(DMC,DC):
    '''Calculates today's Buildup Index
Parameters:
    DMC is the current day's Duff Moisture Code
    DC is the current day's Drought Code

    BUI(8.5450511359999997,19.013999999999999) = 8.4904265358371838'''
    
    Nx = DMC.shape[0]
    Ny = DMC.shape[1]
    
    U = np.zeros((Nx,Ny))*np.nan
    
    c_1 = np.where(DMC <= 0.4 * DC)
    c_2 = np.where(DMC > 0.4 * DC)
    
    U[c_1] = 0.8 * DMC[c_1] * DC[c_1] / (DMC[c_1] + 0.4 * DC[c_1])
    U[c_2] = DMC[c_2] - (1.0 - 0.8 * DC[c_2] / (DMC[c_2] + 0.4 * DC[c_2])) * \
            (0.92 + (0.0114 * DMC[c_2])** 1.7)
                    
    U[np.where(U<0.0)] = 0.0
    return U



def FWI(ISI, BUI):
    '''Calculates today's Fire Weather Index
Paramteres:
    ISI is the current day's ISI
    BUI is the current day's BUI

    FWI(10.853661073655068,8.4904265358371838) = 10.096371392382368'''
    
    Nx = ISI.shape[0]
    Ny = BUI.shape[1]
    
    c_1 = np.where(BUI <= 80.0)
    c_2 = np.where(BUI > 80.0)

    fD = np.zeros((Nx,Ny))*np.nan
    S = np.zeros((Nx,Ny))*np.nan

    fD[c_1] = 0.626 * BUI[c_1]**0.809 + 2.0
    fD[c_2] = 1000.0 / (25.0 + 108.64 * np.exp(-0.023 * BUI[c_2]))

    B = 0.1 * ISI * fD
    
    c_3 = np.where(B > 1.0)
    c_4 = np.where(B <= 1.0)
    
    S[c_3] = np.exp(2.72 * (0.434 * np.log(B[c_3]))** 0.647)
    S[c_4] = B[c_4]

    return S



def DryingFactor(Latitude, Month):

    LfN = [-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5.0, 2.4, 0.4, -1.6, -1.6]
    LfS = [6.4, 5.0, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8]

    if Latitude > 0:
        retVal = LfN[Month]
    elif Latitude <= 0.0:
        retVal = LfS[Month]

    return retVal


def DayLength(Latitude, MONTH):
    '''Approximates the length of the day given month and latitude'''

    DayLength46N = [ 6.5,  7.5,  9.0, 12.8, 13.9, 13.9, 12.4, 10.9,  9.4,  8.0,  7.0,  6.0]
    DayLength20N = [ 7.9,  8.4,  8.9,  9.5,  9.9, 10.2, 10.1,  9.7,  9.1,  8.6,  8.1,  7.8]
    DayLength20S = [10.1,  9.6,  9.1,  8.5,  8.1,  7.8,  7.9,  8.3,  8.9,  9.4,  9.9, 10.2]
    DayLength40S = [11.5, 10.5,  9.2,  7.9,  6.8,  6.2,  6.5,  7.4,  8.7, 10.0, 11.2, 11.8]

    retVal = None

    if Latitude<= 90 and Latitude > 33:
        retVal = DayLength46N[MONTH-1]
    elif Latitude <= 33 and Latitude > 0.0:
        retVal = DayLength20N[MONTH-1]
    elif Latitude <= 0.0 and Latitude > -30.0:
        retVal = DayLength20S[MONTH-1]
    elif Latitude <= -30.0 and Latitude >= -90.0:
        retVal = DayLength40S[MONTH-1]

    if retVal==None:
        raise InvalidLatitude(Latitude)
    return retVal

def calcFWI(MONTH,TEMP,RH,WIND,RAIN,FFMCPrev,DMCPrev,DCPrev,LAT):
    '''Caclulates today's FWI
Parameters:
    TEMP is the 12:00 LST temperature in degrees celsius
    RH is the 12:00 LST relative humidity in %
    WIND is the 12:00 LST wind speed in kph
    RAIN is the 24-hour accumulated rainfall in mm, calculated at 12:00 LST
    FFMCPrev is the previous day's FFMC
    DMCPrev is the previous day's DCM
    DCPrev is the previous day's DC
    LAT is the latitude in decimal degrees of the location for which calculations are being made

    calcFWI(4,17,42,25,0,85,6,15,45.98) = 10.096371392382368'''
    
    ffmc = FFMC(TEMP,RH,WIND,RAIN,FFMCPrev)
    dmc = DMC(TEMP,RH,RAIN,DMCPrev,LAT,MONTH)
    dc = DC(TEMP,RAIN,DCPrev,LAT,MONTH)
    isi = ISI(WIND, ffmc)
    bui = BUI(dmc,dc)
    fwi = FWI(isi, bui)
    
    return fwi

def LawsonEq1(DMC):
    '''National Standard and Best-fit Non-linear Regression Equations
Linking DMC to Forest Floor Moisture Content in
Coastal B.C., Southern Interior B.C. and Southern Yukon

DMC National Standard and Coastal B.C. CWH (2.5-4 cm)^2

LawsonEq1(8.5450511359999997)  = 250.7553985454235'''

    return math.exp((DMC-244.7)/-43.4)+20.0

def LawsonEq2(DMC):
    '''National Standard and Best-fit Non-linear Regression Equations
Linking DMC to Forest Floor Moisture Content in
Coastal B.C., Southern Interior B.C. and Southern Yukon

Southern Interior B.C.3 (2-4 cm)^2

LawsonEq2(8.5450511359999997)  = 194.93023948344205'''
    return math.exp((DMC-223.9)/-41.7)+20.0

def LawsonEq3(DMC):
    '''National Standard and Best-fit Non-linear Regression Equations
Linking DMC to Forest Floor Moisture Content in
Coastal B.C., Southern Interior B.C. and Southern Yukon

Southern Yukon - Pine/White Spruce
Feather moss, Sphagnum and Undifferentiated duff (2-4 cm)^2

LawsonEq3(8.5450511359999997)  = 442.82109267231488'''
    return math.exp((DMC-157.3)/-24.6)+20

def LawsonEq4(DMC):
    '''National Standard and Best-fit Non-linear Regression Equations
Linking DMC to Forest Floor Moisture Content in
Coastal B.C., Southern Interior B.C. and Southern Yukon

Southern Yukon - Pine/White Spruce
Reindeer lichen (2-4 cm)^2

LawsonEq4(8.5450511359999997)  = 746.02210402093272'''
    return math.exp((DMC-106.7)/-14.9)+20.0

def LawsonEq5(DMC):
    '''National Standard and Best-fit Non-linear Regression Equations
Linking DMC to Forest Floor Moisture Content in
Coastal B.C., Southern Interior B.C. and Southern Yukon

Southern Yukon - White Spruce
White spruce/feather moss (2-4 cm)^2

LawsonEq5(8.5450511359999997)  = 853.2397847094652'''

    return math.exp((DMC-149.6)/-20.9)


def calc_rh(p,q,t):
    
    # constants
    svp1=611.2
    svp2=17.67
    svp3=29.65
    svpt0=273.15
    eps = 0.622
    
    # RH formula
    rh = 100 * (p*q/(q*(1.-eps) + eps))/(svp1*np.exp(svp2*(t-svpt0)/(t-svp3)))
    
    return rh





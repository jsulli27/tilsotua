#apply the refraction correction to the mask results

import numpy as np
from astropy.table import Table,Column
def refraction_calc(ra_center,dec_center):
    rph = np.pi*15./180 #radians per hour
    rpd = np.pi/180. #radians per degree
    dpr = 1/rpd #radians per degree
    rpam = rpd/60. #radiams per arcmin
    H = 0 * rph  #chosen hour angle
#=================================================================================================================================
    #generate grid of RA and Dec based on mask center position for the whole field of view of LRIS
    #field of view of LRIS is 6×7.8 arcmin
    delta =  4 * rpam

    ra_range = np.linspace(ra_center-delta,ra_center+delta, num = 100)
    dec_range = np.linspace(dec_center-delta,dec_center+delta, num = 100)
    ra_values = []
    dec_values = []
    ra_offsets = []
    dec_offsets = []
    offsets = []
#=================================================================================================================================
    #create the grid of ra and dec with the ranges there
    ra_grid = np.zeros(shape=(len(ra_range),len(dec_range)))
    dec_grid = np.zeros(shape=(len(ra_range),len(dec_range)))
    for i in range(len(ra_grid)):
        for j in range(len(dec_grid)):
            ra_grid[i:j]=ra_range[i]
            dec_grid[i:j]=dec_range[j]
#=================================================================================================================================
    #first calculate hour angle for the non center points
    HA = H + ra_center-ra_range

    #call function for the non center points
    #call function for the center points

    #set up the function that calculates the refraction

    def rad_refraction(ra,dec,HA):

        rph = np.pi*15./180 #radians per hour
        rpd = np.pi/180. #radians per degree
        dpr = 1/rpd #degrees per radian
        rpam = rpd/60. #radians per arcmin
        rpas = rpd/3600. #radians per arcsec
        pres = 486.  #atmos pressure
        wave = 6000./10000  #wavelength
        lat = 19.828 * rpd #keck latitude
        temp = 0

    #calculate alt and zenith distance
        sina = np.sin(lat)*np.sin(dec)+np.cos(lat)*np.cos(dec)*np.cos(HA)
        cosz = sina
        sinz = np.sqrt(1.-sina**2)
        #print(sinz)
        tanz = sinz/cosz


    #calculate parallactic angle
        if sinz != 0:
            sinq = np.cos(lat)*np.sin(HA)/sinz
            cosq = (np.cos(np.pi/2-lat) - cosz *np.cos(np.pi/2-dec))/(sinz*np.sin(np.pi/2-dec))

        else:
            sinq = 0
            cosq = 0

    #calculate refractive index
        wavers = 1/(wave**2)
        N = 64.328 + (29498.1/(146.-wavers)) + (255.4/(41.-wavers))
        N = (1*10**-6)*N


    #calculate temp and pressure correction
        tcorr = 1+0.003661*temp
        num = 720.88*tcorr
        N = N * ((pres * (1.+((1.049-0.0157*temp) * 10**-6 * pres))) / num)#N*((pres*(1.+(1.049-((0.0157*temp)*(1*10**-6)*pres))))/num)

    #calculate refraction
        R = N * 206265 * tanz

    #calculate correction to RA and Dec
        DA = R * sinq * rpas/np.cos(dec)
        DD = R * cosq * rpas
        return(ra+DA,dec+DD)

    aparent_ra = np.zeros(shape=(len(ra_range),len(dec_range)))
    aparent_dec = np.zeros(shape=(len(ra_range),len(dec_range)))

    for i in range(len(ra_range)):
        for j in range(len(dec_range)):
            temp = rad_refraction(ra_range[i],dec_range[j],HA[i])
            aparent_ra = temp[0]*dpr
            ra_values.append(aparent_ra)
            aparent_dec = temp[1]*dpr
            dec_values.append(aparent_dec)
            ra_offset = (aparent_ra-ra_range[i]*dpr)
            ra_offsets.append(ra_offset)
            dec_offset = (aparent_dec-dec_range[j]*dpr)
            dec_offsets.append(dec_offset)
            offset = np.sqrt(ra_offset**2+dec_offset**2)
            offsets.append(offset)
            
    #create an astropy table that I can then use as reference in the main code
    refraction_results= Table()
    ra_aparent = Column(ra_values, name='Aparent RA')
    dec_aparent = Column(dec_values, name='Aparent Dec')
    raoffsets = Column(ra_offsets, name='RA Offset')
    decoffsets = Column(dec_offsets, name='Dec Offset')

    refraction_results.add_column(ra_aparent)
    refraction_results.add_column(dec_aparent)
    refraction_results.add_column(raoffsets)
    refraction_results.add_column(decoffsets)

    return(refraction_results)

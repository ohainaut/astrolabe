#!/usr/bin/env python3
# coding: utf-8

# # Astrolabe
# 
# A python replacement to the good old astrolabe.for
# draw ``astrolabe'' grids (coordinates and an airmass overlay)
#
# usage:
#   Astrolabe.py
#      debug version, runs with hardcoded parametres
#   Astrolabe.py CGIobs code [projection] [view]
#      for a given observatory by code
#   Astrolabe.py CGIcoo long lat [projection] [view]
#      for a given observatory by coordinates
#   Astrolabe.py CGIband long [projection] [view]
#      for a band
#   code: 3digits (incl leading zero) code in observatory.dat
#     072 Brussels; 272 Munich; 144 VLT; 132 La Silla...
#   long, lat: coordinates, W>0, N>0
#   projection: stereographic; anything else is polar
#   view: outside; anything else is inside
#
#
# version:	sometime in early 92? early version for La Silla
# version:      Apr94, cleaned code
# version:	Jun96. kjm made code more "user friendly"
# version:      Feb20. oh ported to python
# input:
#     stars: in hr.dat ; each line:
#
#                ii, ij, ra, ram, ras, dec, dm, ds, vmag
#            if   ram=ras=0, the program considers that ra and dec are
#            given in degrees
#            if ij>0, the star is a solid circle, open otherwise
#	     ignore the ii - it is not useful.
#     observatory
#             |ctio  -30.0957800   70.4853600 2235.0000000 CTIO, Chile
#----------------------------------------------------------------------

import numpy as np
import matplotlib
matplotlib.use('Agg')  # no Display; must be before pyplot
import matplotlib.pyplot as plt
import sys
import os
from astropy.time import Time
from astropy.io import ascii

debug = False

doRete = True
doPlate = True
doRule = True
#======================================================================
if 0:
    astpath = '/home/ohainaut/public_html/astrolabe/'
    outpath = '/home/ohainaut/public_html/outsideWorld/'
    pubweb = 'https://www.eso.org/~ohainaut/outsideWorld/'
    astweb = 'https://www.eso.org/~ohainaut/astrolabe/'
    logfile = '/home/ohainaut/public_html/log/astrolabe.log'
    obsFile =  astpath+ 'src/observatory.dat'
    obsCleanFile =  astpath+ 'src/obsClean.dat'
    starFile = astpath+ 'src/hr.dat'
    constFile = astpath + 'src/const.dat'
    landFile = astpath+ 'src/landold.dat'
else:
    astpath = os.path.dirname(__file__)
    outpath = './'
    pubweb = './'
    astweb = './'
    obsFile =  astpath+'/observatory.dat'
    obsCleanFile =  astpath+'/obsClean.dat'
    logfile = 'astrolabe.log'
    starFile = astpath+'/hr.dat'
    constFile = astpath+'/const.dat'
    landFile = astpath+'/landold.dat'



#===========================================================================
#===========================================================================
#===========================================================================
def d2r(decd):
    # subroutine to convert Declination into plot radius
    # input: decd declination in degress
    # output: radius in plot units
    if stereo:
        r = np.tan( np.radians( (90.+decd)/2. )) * d2rnorm
        #+dec shoudl be -dec, but we move everything to S ->*-1
        # max: 173
    else:
        r = ( 90.+decd) *d2rnorm #decd in [-90,northlim]
    return r
#===========================================================================
def altaz2hadec(zr,azr): # zr and azr in rad
    # spherical trigo: converts ZenithDist, Azimuth
    # into HA (radians), Dec (dec)
    cosz = np.cos(zr)    
    sinDec = np.sin(np.radians(xlat))*cosz \
        + np.cos(np.radians(xlat))*np.sin(zr)*np.cos(azr)
    Decr = np.arcsin(sinDec)
    Dec = np.degrees(Decr)
    sinHA = -np.sin(azr)*np.sin(zr)/np.cos(Decr)
    cosHA = (cosz - sinDec*np.sin(np.radians(xlat))) \
        / np.cos(Decr)/np.cos(np.radians(xlat))
    
    HAr = np.arctan2(sinHA, cosHA)
    return HAr, Dec
    
#======================================================================
def suncoord(Lr):
    #Lr: longitude of the Sun in Rad
    #returns RAr, Dec of the sun
    
    X = np.cos(Lr) # position of the sun
    Y = np.cos(epsr)*np.sin(Lr)
    Z = np.sin(epsr)*np.sin(Lr)
    R = np.sqrt(1.0-Z*Z)
    deltaSun = np.degrees(np.arctan(Z/R)) #// in degrees
    RASun = (24./180.)*np.degrees(np.arctan2(Y,(X+R))) #// in hours
    RASunr = np.radians(15.*RASun)
    return RASunr, deltaSun
#======================================================================
def sunlong(JD):
    #Returns Sun longitude from Julian Day

    # number of Julian centuries since Jan 1, 2000, 12 UT
    T  = (JD-2451545.0) / 36525. 
    # mean anomaly, degree
    M = 357.52910 + 35999.05030*T - 0.0001559*T*T - 0.00000048*T*T*T
    # mean longitude, degree
    L0 = 280.46645 + 36000.76983*T + 0.0003032*T*T 
    # Sun's equation of center
    DL = (1.914600 - 0.004817*T - 0.000014*T*T)*np.sin(np.radians(M))+\
         (0.019993 - 0.000101*T)*np.sin(np.radians(2*M)) + 0.000290*np.sin(np.radians(3*M))
    # true longitude, degree
    L = L0 + DL
    return np.radians(L)
#===========================================================================
def plotAlmucantar(cosz, color):
    # plot a full almucantar corresponding to cosz
    zr = np.arccos(cosz)
    wHAr, wDec = altaz2hadec(zr, rad)
    HAr = wHAr[ wDec*isouth < northlim*isouth]
    Dec = wDec[ wDec*isouth < northlim*isouth]
    ax.plot(HAr, d2r(Dec*isouth) ,c=color, linewidth=linew)
    return HAr, Dec

#===========================================================================
def plotHalfmucantar(cosz, lw):
    # plot a partial almucantar corresponding to cosz
    zr = np.arccos(cosz)
    HAr, Dec = altaz2hadec(zr, rad)

    if isouth > 0:
        w1x = HAr[ deg > 30 ]
        w1y = Dec[ deg > 30 ]
        w1d = deg[ deg > 30 ] 
        w2x = HAr[ deg < -30 ]
        w2y = Dec[ deg < -30 ]
        w2d = deg[ deg < -30 ]
        wHAr = np.concatenate((w1x, w2x), axis=None)
        wDec = np.concatenate((w1y, w2y), axis=None)
        wd =  np.concatenate((w1d, w2d), axis=None)
    else:
        w1x = HAr[ deg > -150 ]
        w1y = Dec[ deg > -150 ]
        w1d = deg[ deg > -150 ] 

        wHAr = w1x[ w1d < 150 ]
        wDec = w1y[ w1d < 150 ]
        wd   = w1d[ w1d < 150 ]

    HAr = wHAr[ wDec*isouth < northlim*isouth]
    Dec = wDec[ wDec*isouth < northlim*isouth]
        
    ax.plot(HAr, d2r(Dec*isouth) ,c=airmasscol, linewidth=lw)
    return wHAr, wDec
#==============================================================================
def xy2ad(x,y):
    # finds the RA, radius for a position in X,Y
    a = np.arctan2(x,y)
    r =  np.sqrt(x**2+y**2)
    return a,r
#------------------------------------------------------------------------------
def write_legend_line(x,y,label,side):
    a,r = xy2ad(x,y)
    plt.text(a,r,label,
             color=legendc,
             fontsize=legends,
             horizontalalignment=side,
             verticalalignment='top')
#----------------------------------------------------------------------    
def write_legend(iwhat):
    dy = 8

    side = 'left'
    # writes the labels in the corner
    x = -200. * tradOrient*isouth*iwhat
    y = 200.
    write_legend_line(x,y,'www.eso.org/~ohainaut/bin/astrolabe.cgi',side)

    y -= dy
    a,r = xy2ad(x,y)
    write_legend_line(x,y,'ohainaut@eso.org',side)


    side = 'right'
    x = 200. * tradOrient*isouth*iwhat
    y = 200.
    a,r = xy2ad(x,y)
    write_legend_line(x,y,'Long, Lat: '+longlatlabel,side)
   
    y -= dy
    if iwhat < 0:
        write_legend_line(x,y,'Name: {} ({})'.format(obsname, obs),side)

    else:    
        write_legend_line(x,y,'Name: '+rete,side)

    y -= dy
    if stereo:
        wlab = 'stereographic'
    else:
        wlab = 'polar'
    write_legend_line(x,y,'Projection: '+wlab,side)

    y -= dy
    if tradOrient > 0:
        wlab = 'outside'
    else:
        wlab = 'inside'
    write_legend_line(x,y,'Orientation: '+wlab,side)
            
    y -= dy
    if isouth > 0:
        wlab = 'S'
    else:
        wlab = 'N'
    write_legend_line(x,y,'Hemisphere: '+wlab,side)

    y -= dy
    write_legend_line(x,y,'Limit dec={}deg'.format(northlim),side)

    #crosses at the corners
    side = 'center'
    xx = [-200,200]
    yy = [-200,200]
    for x in xx:
        for y in yy:
            write_legend_line(x,y,'+',side)

#==============================================================================
def plot_milkyway():
    #- galactic disc

    import astropy.units as u
    from astropy.coordinates import SkyCoord

    glong = np.arange(-180,181,1)

    # Bulge
    linew = linew_thin
    linea = linea_transp
    for i in np.arange(-10,11,1):
        glat = glong*0.+ 2*i* np.exp ( -(glong/30.)**2 )  #<- bulge flare
        glx = SkyCoord(l=glong*u.degree, b=glat*u.degree, frame='galactic')
        gra = glx.transform_to('icrs').ra.degree
        gde = glx.transform_to('icrs').dec.degree 
        linew = linew_thin
        linea = linea_normal
        ax.plot(np.radians(gra),d2r(gde*isouth),c=glxcol, linewidth=linew, alpha=linea)

    # +10 and -10 
    for i in np.arange(-10,11,10):
        if i == 0:
            linew = linew_medium
        else:
            linew = linew_thin

        glat = glong*0.+ i
        glx = SkyCoord(l=glong*u.degree, b=glat*u.degree, frame='galactic')
        gra = glx.transform_to('icrs').ra.degree
        gde = glx.transform_to('icrs').dec.degree 
        ax.plot(np.radians(gra),d2r(gde*isouth),c=glxcol, linewidth=linew, alpha=linea)


    # galactic long tickmarks
    glat =  np.arange(3,-1,-3)
    linew = linew_thin
    
    for i in np.arange(0,360,15):
        glong = glat*0. + i
        glx = SkyCoord(l=glong*u.degree, b=glat*u.degree, frame='galactic')
        gra = glx.transform_to('icrs').ra.degree
        gde = glx.transform_to('icrs').dec.degree 
        ax.plot(np.radians(gra),d2r(gde*isouth),c=glxcol, linewidth=linew, alpha=linea)

        if d2r(gde[0]) +2. < outring: # only if in view, label tick
            ax.text(np.radians(gra[0]),
                    d2r(gde[0]*isouth) +2.,
                    '{:3d}'.format(int(glong[0])),
                    rotation= -tradOrient*gra[0]*  isouth  +180.,
                    horizontalalignment='center',
                    verticalalignment='center',
                    color=glxcol,
                    fontsize=5
            )

#======================================================================
def plot_ecliptic():
    #====  ECLIPTIC
    # old simple ecliptic
    if False:
        xinc = 23.439
        delta =  xinc*np.sin(rad)  *isouth ## 
        linew = linew_thin
        ax.plot(rad,d2r(delta),c='g', linewidth=linew)
        
        #months on ecliptic
        myd = np.arange(0,12)*30.+9  # offset so that eqx happens on ~21
        myr = xinc * np.sin(np.radians(myd))  *isouth ## 
        ax.scatter(np.radians(myd),d2r(myr),s=10,c='g')
        
        for i in np.arange(0,len(Month)):
            ax.text(np.radians(myd[i]),
                    d2r(myr[i])+7,  ## myr already has been multiplied by isouth
                    Month[i],color='g',
                    rotation=-tradOrient*myd[i]*isouth + 180.,
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=7)

        # 15 day ticks on ecliptic
        myd = np.arange(0,360,5)+9
        myr = xinc * np.sin(np.radians(myd))  *isouth ## 
        ax.scatter(np.radians(myd),d2r(myr),s=.5,c='g')

    
    # FULL ecliptic

    # ecliptic circle
    linew = linew_thin
    linea = linea_normal
    Lr = np.radians(np.arange(0.,360.1,1.)) # longitude of the sun
    RASunr, deltaSun = suncoord(Lr)
    deltaSun = deltaSun * isouth

    ax.plot(RASunr,d2r(deltaSun),
            c=solcol, linewidth=linew, alpha=linea)
    ax.plot(RASunr,d2r(deltaSun)+0.25,
            c=solcol, linewidth=linew, alpha=linea)

    #ecliptic longitude ticks
    Li = np.arange(0,360)
    L = Li*1.
    Lr = np.radians(L) # longitude of the sun
    RASunr, deltaSun = suncoord(Lr)    
    deltaSun = deltaSun * isouth

    for i in Li: # 360 deg
        if i % 30 == 0:
            Sl=8
            ax.text(RASunr[i], d2r(deltaSun[i])-Sl-2,
                    '{:03d}'.format(i),
                    color=solcol, fontsize=5,
                    horizontalalignment='center',
                    verticalalignment='center',
                    rotation=-tradOrient*i*isouth+180.)
        elif i % 10 == 0:
            Sl=5
        elif i % 5 == 0:
            Sl=3
        else:
            Sl = 2
        
        ax.plot([RASunr[i],RASunr[i]],
            [d2r(deltaSun[i]),d2r(deltaSun[i])-Sl],
            c=solcol, linewidth=linew, alpha=linea)


    # Ecliptic Calendar
    for i in np.arange(0,len(SunJD)):
        SLab = SunDay[i]
        SJD =  SunJD[i]
        Lr = sunlong(SJD)
        RASunr, deltaSun = suncoord(Lr)    
        deltaSun = deltaSun * isouth

        # ecliptic calendar tick marks every 5d
        ax.plot([RASunr,RASunr],
            [d2r(deltaSun)+0,d2r(deltaSun)+5.],
            c=solcol, linewidth=linew, alpha=linea)

        # label the tickmarks
        ax.text(RASunr, d2r(deltaSun)+8,
            SLab,
            color=solcol, fontsize=8,
            horizontalalignment='center',
            verticalalignment='center',
            rotation= -tradOrient*np.degrees(RASunr)*  isouth+180.)

        # little tick marks between calendar dates
        if i == 0:
            RASunr0, deltaSun0 = suncoord(sunlong(SunJD[-1]))
            deltaSun0 = deltaSun0 * isouth
            #print(SunJD[-1], SunJD[len(SunJD)-1], RASunr0, RASunr)
            
        # ecliptic calendar
        if RASunr-RASunr0 < -2.: # catch the 23:59 - 00:00 crossing
            RASunr0 -= 2.*np.pi
        if RASunr-RASunr0 > 2.: # and in the other direction
            RASunr0 += 2.*np.pi
        for j in np.arange(1,int(SunStep)): # interpollate little ticks
            RASunrj = RASunr0 + j*(RASunr-RASunr0)/SunStep
            deltaSunj = deltaSun0 + j*(deltaSun-deltaSun0)/SunStep
            ax.plot([RASunrj,RASunrj],
                    [d2r(deltaSunj),d2r(deltaSunj)+2],
                    c=solcol, linewidth=linew, alpha=linea)
               

        RASunr0 = RASunr*1. # keep current value for next i, for interpollation
        deltaSun0 = deltaSun*1.


#======================================================================
def plot_stars():
    #====  STARS

    #- constellation

    # read the constellation files
    constData = ascii.read(constFile, delimiter=',')
    print( "Nr of constellation points:",len(constData))
    #Const,HR1,RA1,Dec1,Mag1,HR2,RA2,Dec2,Mag2

    racav = 0. # init, will be the average coords of the stars in constellation
    decav = 0.
    iicav = 0
    lacav = "XX"
    for i in np.arange(0, len(constData)):
        if constData['Dec1'][i] *isouth < northlim*isouth and \
           constData['Dec2'][i] *isouth < northlim*isouth: # is the line visible

            if constData['Const'][i] != lacav: # new constellation name

                if iicav > 0: # we finished a constellation that had lines
                              # let's label it
                    racav = racav/2./iicav
                    decav = decav/2./iicav*isouth
                    ax.text(np.radians(racav),d2r(decav), lacav,
                            color=colg,
                            rotation = -tradOrient*racav*isouth + 180.,
                            horizontalalignment='center',
                            verticalalignment='center')

                    iicav = 0 # and re-init the averages
                    racav = 0.
                    decav = 0.

            racav += constData['RA1'][i]  + constData['RA2'][i]
            decav += constData['Dec1'][i] + constData['Dec2'][i]
            iicav += 1
            lacav = constData['Const'][i] 
            
            #-draw the constellation line
            ax.plot([np.radians(constData['RA1'][i]), np.radians(constData['RA2'][i])],
                [d2r(constData['Dec1'][i] *isouth), d2r(constData['Dec2'][i] *isouth)],
                alpha=0.5, linewidth=.5,color=colk)

            
    #- stars
    starData = np.genfromtxt(starFile)  # read star data file
    # ii, ik, ra, junk, junk, dec, junk, junk, vmag
    starRA = starData[:,2]
    starDec = starData[:,5]
    starMag = starData[:,8]
    starSize = 15 - 3*starMag
    
    wrad =  d2r(starDec *isouth )
    selRAr = np.radians(starRA[wrad < outring])
    selrad = wrad[wrad < outring]
    selsize = starSize[wrad < outring]
    
    ax.scatter(selRAr, selrad,
               c=colk, s=selsize, 
               linewidth=0,
               alpha=1.)

    #- alternate source of stars
    if False:
        ax.scatter(np.radians(constData['RA1']), d2r(constData['Dec1'] *isouth), 
           c=colk, s=starSize, 
           linewidth=0,
           alpha=0.7)


    
#======================================================================
def plot_landold():
    #- landold standard fields
    landData = np.loadtxt(landFile,
            dtype={'names': ('j1','j2','rh','rm','rs','dd','dm','ds','j3','name'),
               'formats':(np.float,np.float,np.float,np.float,
                          np.float,np.float,np.float,np.float,np.float,'|S15')})
    landData = ascii.read(landFile, delimiter=' ')
    # ii, ik, RAH, RAM, RAS, DDE, DM, DS, j, ID
    landRA = 15.*(landData['col3'] + landData['col4']/60. + landData['col5']/3600.)
    landDec = landData['col6']
    landID = landData['col10']
    
    ax.scatter(np.radians(landRA), d2r(landDec *isouth),
               c=colr, s=10, alpha=0.5)
    for i in np.arange(0,len(landID)):
        ax.text(np.radians(landRA[i]),
                d2r( landDec[i]*isouth) +15.,landID[i],color=colr,
                rotation=landRA[i]+90,
                fontsize=5,
                horizontalalignment='center',
                verticalalignment='center')

#======================================================================
def plot_rete():
    # MAIN GRID PLOT - RETE
    if True :
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)

    #!!! polar expects y_limits to be so that center < edge
    #    limits(-90,30) will work, limits(90,-30) will not
    #    ===> all declinations are multiplied by isouth
    #         to reverse all Northern plots
    ax.set_ylim(0, outring +20.)
    ax.set_theta_direction(-isouth*tradOrient)
    ax.set_theta_zero_location("N")
    
    write_legend(1) # write stuff in the corners - 1=rete



    # Galactic disk - first, so that everything else ovewrites it.
    plot_milkyway()

    
    # RA tick
    linew= linew_thin
    linea = linea_normal
    calrad = outring + 10 # radius of the calendar
    
    for i in np.arange(0,360, 2.5):
        myr = np.array([calrad,calrad -5])
        myd = [i,i]
        ax.plot(np.radians(myd),(myr),c=colk, linewidth=linew, alpha=linea)
            
        myr = np.array([calrad,calrad -2])
        myd = [i+1.25,i+1.25] #5min
        ax.plot(np.radians(myd),(myr),c=colk, linewidth=linew, alpha=linea)


    if debug:
        print( 'dbg')
        plt.savefig(outpath+'debug.pdf')
        exit(1)

    # radial lines
    linew = linew_thin

    for h in hours:
        # RA lines
        dm=-50*isouth
        if h % 2 == 0:
            dm = -80*isouth
        else:
            dm = -90*isouth
            linew = linew_medium
        myd = np.array([h, h])*15.
        myr = np.array([dm, northlim]) *isouth ## 
        ax.plot(np.radians(myd),d2r(myr),c=colg, linewidth=linew, alpha=linea)

        # RA labels

        if h % 6 == 0:
            ralabel = 'RA= {:02d}'.format(h)
            ax.text(np.radians(h*15),
                    calrad -13.,
                    'UT for calendar',
                    color=calcol,
                    rotation=-tradOrient*h*15*isouth +180,
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=5)
        else:
            ralabel = '{:02d}'.format(h)

        ax.text(np.radians(h*15),
                    calrad -6,
                    ralabel,
                    color=colk,
                    rotation=-tradOrient*h*15*isouth +180,
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=10)

            
    # Declination circles
    for r in np.arange(-80, northlim*isouth +1, 10):
        if r == 0:
            linew = linew_thick # thicker equator
        else:
            linew = linew_thin
        myr = deg*0+ r  
        ax.plot(rad,d2r(myr), c=colg, linewidth=linew, alpha=linea)

    #label
    for r in np.arange(-60, northlim*isouth, 30):
        ax.text(np.radians(180.),d2r(r),
                '$\delta$={:+d}$^o$'.format(int(r*isouth)),
                color=colk,
                rotation=0.,
                horizontalalignment='center',
                verticalalignment='top',
                fontsize=7)

    
    #Outer ring
    ax.plot(rad,rad*0+ outring, c=colk, linewidth=linew, alpha=linea)
    ax.plot(rad,rad*0+ calrad , c=colk, linewidth=linew, alpha=linea)



    # the star catalogue and the constellations
    plot_stars()
    
    # the Ecliptic and its calendar
    plot_ecliptic()

#======================================================================
def plot_plate():
    # OVERLAY PLOT - MATER & PLATE


    #!!! polar expects y_limits to be so that center < edge
    #    limits(-90,30) will work, limits(90,-30) will not
    #    ===> all declinations are multiplied by isouth
    #         to reverse all Northern plots
   
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.set_ylim(0., outring +20.)
    ax.set_theta_direction(isouth*tradOrient)
    ax.set_theta_zero_location("N")
    
    
    write_legend(-1) # write stuff in the corners, -1 for mater
    linea = linea_normal


    #---------------------------------------------------------------
    #- Azimuth grid
    
    linea = linea_normal
    linec = colb
    
    for az in np.arange(0,360,5):
        if az % 90 == 0:
            linew = linew_thick
        else:
            linew = linew_thin
        
        if az % 90 == 0:
            zr = np.radians(np.arange(90.,-.01,-3.))
        elif az % 15 == 0:
            zr = np.radians(np.arange(90.,14.,-3.))
        else:
            zr = np.radians(np.arange(90.,86,-3.))

        azr = np.radians(az)
        wHAr, wDec = altaz2hadec(zr, azr)
        HAr = wHAr[ wDec*isouth < northlim*isouth]
        Dec = wDec[ wDec*isouth < northlim*isouth]
        ax.plot(HAr, d2r(Dec*isouth) , c=linec, linewidth=linew, alpha=linea)

        #labels
        if az % 15 == 0:
            zr = np.radians(np.arange(90.,-.01,-3.))

            x0 = np.cos(HAr[0])*d2r(Dec[0]*isouth)
            x1 = np.cos(HAr[1])*d2r(Dec[1]*isouth)
            
            y0 = np.sin(HAr[0])*d2r(Dec[0]*isouth)
            y1 = np.sin(HAr[1])*d2r(Dec[1]*isouth)
            labang = -tradOrient*np.degrees(np.arctan2( (y1-y0), isouth*(x1-x0)))

            if az == 90. or az == 270.:
                labang = -labang
                if isouth < 0:
                    labang += 180
                xl = x0 - 3.*(x1-x0)
                yl = y0 - 3.*(y1-y0)
                rl = np.sqrt(xl*xl+yl*yl)
                azl = np.arctan2(yl,xl)
                if az == 90:
                    letter = "E"
                else:
                    letter = "W"
                ax.text(azl,rl,
                    letter,
                    rotation=labang,
                    fontsize = 20,
                    color=colb,
                    horizontalalignment='center',
                    verticalalignment='center')

            if az < 180:
                labang = -labang -90
            else:
                labang = -labang +90
                    
            labrad = d2r(Dec[2]*isouth)
            if labrad < outring:
                ax.text(HAr[2], labrad,
                        '{:03.0f}'.format(az),
                        rotation=labang,
                        fontsize = 6,
                        horizontalalignment='center',
                        verticalalignment='center')


    
    #-  ZD grid--------------------------------------------------

    # airmass 2
    z = 60. 
    linew = linew_thick
    _ = plotAlmucantar(np.cos(np.radians(z)), zdcol)
    
    # altitude almucanars
    linew = linew_thin
    for z in np.arange(15,90.,15.):
        wh, wdelta = plotAlmucantar(np.cos(np.radians(z)), zdcol)
        if isouth < 0:
            myi = 0
        else:
            myi = int(len(wh)/2)    
        ax.text(wh[myi], d2r(wdelta[myi]*isouth+0.),
                '  h={:2d}$^o$ '.format(int(90-z)),
                color=linec,
                fontsize=4,
                rotation=180.,
                horizontalalignment='right',
                verticalalignment='top')


    #- airmass grid
    linew = linew_thin
    linea = linea_normal
    linec = collb
    for secz in np.array([1.1,1.2,1.3,1.6,1.8,2.5,3.,5., 7.]):
        cosz = 1./secz
        wh, wdelta = plotAlmucantar(cosz, airmasscol)
        if isouth < 0:
            myi = 0
        else:
            myi = int(len(wh)/2)
        ax.text(wh[myi],d2r(wdelta[myi]*isouth),
                    ' z={:2.1f}  '.format(secz),
                    color=linec,
                    fontsize=4,
                    rotation=180.,
                    horizontalalignment='left',
                    verticalalignment='top')

    #- Horizon
    linea = linea_normal
    z = 90. 
    linew = 3.
    _ = plotAlmucantar(np.cos(np.radians(z)), airmasscol)

    #- Twilights
    myhlabel = ['Horizon', 'Civil Twilight', 'Nautical Twilight', 'Astronomical Twilight']
    z = [90., 96., 102., 108.]
    for i in np.arange(0,len(z)):
        wh, wdelta = plotHalfmucantar( np.cos(np.radians(z[i])), linew )
        myi = int(len(wh)/2)
        ax.text(wh[myi],
            d2r(wdelta[myi]*isouth)+4.,
            myhlabel[i],color=linec,
            fontsize=4,
            horizontalalignment='center',
            verticalalignment='center')
        linew -= .7


    #- white circle to mask the tails of twilights
    myr = rad*0 + outring
    ax.plot(rad,myr, color='white', linewidth=12)

    
    #- Solar Time Circle
    # (formarly UT)

    linew = linew_thin
    linea = linea_normal
    colut = solcol

    #- Solar time ticks
    rin = outring -5
    for m in np.arange(0, 1440): #minutes
        myd = np.radians(np.array([m,m])/4.) # degrees, then rad
        if m % 60 == 0:
            myr = np.array([rin+7, rin]) ##
        elif m % 10 == 0:
            myr = np.array([rin+5, rin]) ##
        elif m % 5 == 0:
            myr = np.array([rin+3, rin]) ##
        else:
            myr = np.array([rin+2, rin]) ##

        ax.plot(myd,myr,c=colut, linewidth=linew, alpha=linea)


    #- circles
    linew = linew_medium
    myr = rad*0 + rin -5.
    ax.plot(rad,myr, color=colut, linewidth=0.3)
    myr = rad*0 + rin
    ax.plot(rad,myr, color=colut, linewidth=0.3)

    #- Hours labels
    for h in hours:
        if h == 12: # special marker for sidereal time
            ax.plot([np.radians(h*15. + 180.)], 
                    [rin +8.],
                    marker='^',
                    color=colk)
                    
        else:
            ax.text(np.radians(h*15. + 180.), 
                    rin +8.,
                    '{:02d}'.format(h), 
                    rotation=tradOrient*isouth*h*15. ,
                    color=colut,
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=8)


    #- "Solar Time" label
    for i in [90,180,270]:
        lang = i - tradOrient*isouth * 3.
        ax.text(np.radians(lang), 
            rin+7.,
            'Solar', 
            rotation= tradOrient*isouth*lang + 180,
            color=colut,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=5)
        lang = i + tradOrient* isouth * 3.
        ax.text(np.radians(lang), 
            rin+7.,
            'Time',
            rotation= tradOrient*isouth*lang + 180,
            color=colut,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=5 )

        
    #- "Sid Time" label
    lang =  - tradOrient*isouth * 4.
    ax.text(np.radians(lang), 
            rin+7.,
            'Sidereal', 
            rotation= tradOrient*isouth*lang + 180,
            color=colk,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=5)
    lang =  + tradOrient* isouth * 3.
    ax.text(np.radians(lang), 
            rin+7.,
            'Time',
            rotation= tradOrient*isouth*lang + 180,
            color=colk,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=5)



    #- Longitude circle

    collong = colg
    offlong = xlong
    rin = outring -10.
    for d in np.arange(-179,181): 
        myd = -np.radians(d+180.  -offlong) # degrees, then rad; 0 at bottom; # change orientation
        if d % 30 == 0:
            myr = np.array([rin-2, rin+5.]) ##
            if d != 180:
                ax.text(myd, rin-5,
                    '{:3d}$^o$'.format(d), 
                    rotation= - tradOrient*isouth*(d -offlong),
                    color=collong,
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=8)

        elif d % 10 == 0:
            myr = np.array([rin, rin+5.]) ##
        elif d % 5 == 0:
            myr = np.array([rin +2., rin+5.]) ##
        else:
            myr = np.array([rin +3, rin+5.]) ##
        ax.plot([myd,myd],myr,c=colk, linewidth=linew, alpha=linea)


    lang = 180 +offlong  - 6*isouth*tradOrient
    ax.text(np.radians(lang), 
            rin-2.,
            'Observatory', 
            rotation= tradOrient*isouth*lang + 180,
            color=collong,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=5)
    lang = 180 +offlong  +6 * isouth*tradOrient
    ax.text(np.radians(lang), 
            rin-2.,
            'longitude', 
            rotation= tradOrient*isouth*lang + 180,
            color=collong,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=5 )
    
    lang = 155 + offlong
    ax.text(np.radians(lang), 
            rin-2.,
            'West', 
            rotation= tradOrient*isouth*lang + 180,
            color=collong,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=5 )
    lang = 205 + offlong
    ax.text(np.radians(lang), 
            rin-2.,
            'East', 
            rotation= tradOrient*isouth*lang + 180,
            color=collong,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=5 )

    #- circle
    myr = rad*0 + rin +5.
    ax.plot(rad,myr, color=collong, linewidth=0.3)
    myr = rad*0 + rin
    ax.plot(rad,myr, color=collong, linewidth=0.3)


    #- Observatory name
    ax.text(np.radians(180.), 
            outring*0.55 , 
            obsname, 
            rotation=180,
            color=collb,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=20  )
    
    #- Observatory coords
    ax.text(np.radians(180.), 
            outring*0.65, 
            longlatlabel, 
            rotation=180,
            color=colb,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=10  )

    if band: #scan and plot relevant observatories
        with open(obsCleanFile) as in_file:
            #- Observatory name
            for line in in_file:
                lobs,w = line.split("%")
                lobsname = w.rstrip()
                lobsCode, lobslat, lobslong, lobsalt = lobs.split()
                if lobsCode[0] != '#':
                    wlat = float(lobslat)
                    wlong = -(float(lobslong) -offlong)

                    if abs(wlat-xlat) <= xlatband:
                        longr = np.radians(wlong +180.)

                        if isouth > 0: # south
                            if wlong >= 0:
                                rotation = tradOrient*isouth*wlong -90*tradOrient
                            else:
                                rotation = tradOrient*isouth*wlong +90*tradOrient
                        else: # North
                            if wlong >= 0:
                                rotation = tradOrient*isouth*wlong +90*tradOrient 
                            else:
                                rotation = tradOrient*isouth*wlong -90*tradOrient 

                        horalig = 'center'
                        vertalig = 'center'
                        ax.text(longr,
                                rin -10 -len(lobsname) , 
                                lobsname, 
                                rotation=rotation,
                                color=colg,
                                horizontalalignment=horalig,
                                verticalalignment=vertalig,
                                fontsize=5)
                        #print( lobsname, wlat, xlat, wlong)

                        ax.plot([longr,longr], [rin-10.,rin],
                                c=colg, linewidth=linew, alpha=linea)
                            
    

    #- N S E W
    if isouth == 1: # southern hmsph
        #ax.text(np.radians(0.), 
        #        outring-10., 
        #        "N", 
        #        rotation=180,
        #        color=colb,
        #        horizontalalignment='center',
        #        verticalalignment='center',
        #        fontsize=20 )
        ax.text(np.radians(180.), 
                outring-40 , 
                "S", 
                rotation=0,
                color=colb,
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=20 )

    else: # northern hmsph
        #ax.text(np.radians(0.), 
        #        outring -10. , 
        #        "S", 
        #        rotation=180,
        #        color=colb,
        #        horizontalalignment='center',
        #        verticalalignment='center',
        #        fontsize=20)
        ax.text(np.radians(180.), 
                outring-40, 
                "N", 
                rotation=0,
                color=colb,
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=20  )


    # N S E W lines
    linew = linew_thick
    for h in np.arange(0,24,6):
        myr = np.array([-89., -80])
        myd = myr*0. + h*15.
        ax.plot(np.radians(myd),d2r(myr),c=colr, linewidth=linew, alpha=linea)

    #little circle at the centre
    ax.plot(rad,d2r(rad*0.-89.),c=colr, linewidth=linew, alpha=linea)


    #- Calendar

    calrad = outring + 10 # radius of the calendar
    linew = linew_thin

    for i in np.arange(0,len(SunJD)):
        SLab = SunDay[i]
        SJD =  SunJD[i]
        Lr = sunlong(SJD)
        RASunr, deltaSun = suncoord(Lr)    
        deltaSun = deltaSun * isouth
        
        # outside ring calendar tick marks - shifted for Observatory longitude UT correction.
        # RAposition of calendar date (at 0UT) on the RA grid = local ST of UTdate0UT +12h
        # (12h shift because overlay grid is centred on 12h).
        # GMST = ST of Greenwhich Jan.1 at 0UT
        # long: longitude of the observatory
        # xst = GMST - long/15: ST of observatory on Jan1 at 0UT
        # For Jan1:  RA = xst +12h
        # For day:  RA = xst +12h + day/365.25 *24h
        
        # corrected for longitude - calendar is for UT@observatory
        RAr = (i*SunStep*2.*np.pi/365. + np.radians(xst*15.))
        
        ax.plot(np.array([RAr,RAr]),
                [calrad, calrad+6.],
                c=calcol, linewidth=linew, alpha=linea)

        # ... and their label
        ax.text(RAr, 
                calrad +8., 
                SLab, 
                rotation= tradOrient*np.degrees(RAr)*  isouth  +180.,
                horizontalalignment='center',
                verticalalignment='center',
                color=calcol,
                fontsize=5  )
        #  little tickmarks - by simple interpollation
        if i < len(SunJD)-2: # we skip the last one
            # outside ring
            for j in np.arange(1,7): # the deltaDate is 6 days
                RArj = RAr + j*2.*np.pi/365.
                ax.plot(np.array([RArj,RArj]) ,
                    [calrad, calrad+3],
                    c=calcol, linewidth=linew, alpha=linea)
        else:
            for j in np.arange(1,8): # the deltaDate is 6 days
                RArj = RAr + j*2.*np.pi/365.
                ax.plot(np.array([RArj,RArj]) ,
                    [calrad, calrad+3],
                    c=calcol, linewidth=linew, alpha=linea)

    #- calrad circle    
    ax.plot(rad,rad*0+calrad,
            c=calcol, linewidth=linew, alpha=linea)

    #- UT calendar label
    lang = 180 - tradOrient*isouth * 5.
    ax.text(np.radians(lang), 
            calrad-3,
            'UT calendar', 
            rotation= tradOrient*isouth*lang + 180,
            color=calcol,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=5 )
    lang = 180 + tradOrient* isouth * 5.
    ax.text(np.radians(lang), 
            calrad-3,
            'for Observatory',
            rotation= tradOrient*isouth*lang + 180,
            color=calcol,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=5    )

#======================================================================
def plot_rule():

    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)

    #!!! polar expects y_limits to be so that center < edge
    #    limits(-90,30) will work, limits(90,-30) will not
    #    ===> all declinations are multiplied by isouth
    #         to reverse all Northern plots
    ax.set_ylim(0, outring +20.)
    ax.set_theta_direction(-1) # does not depend on tradOrient
    ax.set_theta_zero_location("N")
    
    write_legend(1) # write stuff in the corners - 1=rete

    #little circle at the centre
    ax.plot(rad,d2r(rad*0.-89.),c=colr, linewidth=linew, alpha=linea)

    # dec scale
    for dec in np.arange(-90, northlim*isouth+0.1, 1.):
        if dec % 10 == 0:
            dy = 5
        elif dec % 5 == 0:
            dy = 4
        else:
            dy = 2.

        w = d2r(dec)
        x = np.array([w,w])
        y = np.array([0, -dy])
        a,r = xy2ad(x,y)
        ax.plot(a,r, linewidth=linew_thin, c=colk)

        a,r = xy2ad(w,-10)
        if dec % 10 == 0:
            ax.text(a,r,
                    '{:3d}$^o$'.format(int(isouth*dec)),
                    color=colk,
                    horizontalalignment='center',
                    verticalalignment='center',
                    fontsize=6,
                    rotation=90)

    # dec scale label
    w = d2r(0.)
    a,r = xy2ad(w,5.)
    ax.text(a,r,
            'Dec.$\delta$'.format(int(dec)),
            color=colk,
            fontsize=6,
            horizontalalignment='center',
            rotation=0)

    # dec scale line
    a = np.radians(np.array([90.,90.]))
    r = np.array([0., outring])
    ax.plot(a,r, linewidth=linew_thin, c=colk)

    # solar time line
    a = np.radians(np.array([0.,270.]))
    r = np.array([0.,outring-18])
    ax.plot(a,r, linewidth=linew_thin, c=colk)

    # solar pointer
    a = np.radians(np.array([270.]))
    r = np.array([outring-8])
    ax.plot(a,r,'<',color=solcol)

    # solar time label
    ax.text(a, outring-13, 'Solar Time',
            fontsize=6,
            color=solcol,
            rotation=270.,
            horizontalalignment='center',
            verticalalignment='center')

        
#===========================================================================
#===========================================================================
#===========================================================================
#================= MAIN ====================================================
#===========================================================================
#===========================================================================

#--Arguments----------------------------------------------------------------
stereo = False
band = True
xlatband = 10.
tradOrient = -1 # +1 = traditional (outside view), or -1 = direct (inside view)
obsCode = 0
obslat = 0.
obslong = 0

if len(sys.argv) == 1: # [0] is the name of the command...
    xlong = 0.
    xlat = 40.
    obs = 'TEST'
    obsname = 'DEBUG'
    stereo = False # stereo projection, or basic projection

    print( "Usage:")
    print( "   Astrolabe.py CGIcoo long lat")
    print( "   Astrolabe.py CGIobs obscode")
    print( "observatories in ", obsFile)
    exit(0)
    
elif sys.argv[1] == 'CGIband':
    xlong = 0.
    xlat =  np.float(sys.argv[2])
    xlatband = 5.
    band = True
    obs = ' '
    obsname = ' '
    if len(sys.argv) >= 4: #
        stereo = (sys.argv[3] == 'stereo')
        if len(sys.argv) == 5: #
            if sys.argv[4] == 'outside':
                tradOrient = 1

elif sys.argv[1] == 'CGIcoo':
    xlong = np.float(sys.argv[2])
    xlat =  np.float(sys.argv[3])
    obs = ' '
    obsname = ' '
    print( "<p>(AstrolabePY)<font size=-3 color=Gray><pre>")
    print( 'len GIcoo', len(sys.argv) )
    if len(sys.argv) >= 5: #
        stereo = (sys.argv[4] == 'stereo')
        if len(sys.argv) == 6: #
            if sys.argv[5] == 'outside':
                tradOrient = 1

        
elif sys.argv[1] == 'CGIobs':
    obs = sys.argv[2]
    with open(obsFile) as in_file:
        for line in in_file:
            lobs,w = line.split("%")
            obsname = w.rstrip()
            obsCode, obslat, obslong, obsalt = lobs.split()
            if obsCode == obs:
                break
    xlat = float(obslat)
    xlong = float(obslong)
    if len(sys.argv) >= 4: #
        stereo = (sys.argv[3] == 'stereo')
        if len(sys.argv) == 5: #
            if sys.argv[4] == 'outside':
                tradOrient = 1
else:
    exit('Expecting\n CGIcoo long lat\nor\nCGIobs code')

print( "Obsname: {}\ncode: {}\nlat: {}\nlong {}".format(
    obsname, obsCode, obslat, obslong))
print( "xlat: {}\nxlong: {}".format(
    xlat, xlong))

if tradOrient > 0:
    olab = 'out'
else:
    olab = 'in'


if stereo:
    print( 'proj: Stereo')
else:
    print( 'proj: Polar')

print( 'Orientation: ', tradOrient, olab)



# plot configuration ------- -------------------------------------------------

plotsize = 7.5 * 1.05 # in inch; will be reduced by 5% bc of margins

linew_thin = 0.3
linew_medium = 0.6
linew_thick = 1.2
linew = linew_medium *1.

linea_normal = 1.
linea_transp = 0.3
linea = linea_transp *1.

if True:
    colk = 'k'
    colr = 'crimson'
    colg = 'Grey'
    colc = 'paleturquoise' #'palevioletred'
    colb = 'mediumblue'
    collb = 'cornflowerblue'
else:
    colk = 'k'
    colr = 'k'
    colg = 'k'
    colc = 'k'
    colb = 'k'
    collb = 'k'


calcol = colb # colour of the calendar
solcol = colr # colour of the ecliptic and solar
glxcol = colc

airmasscol = collb
zdcol = colb
azcol = colb


linec = 'y'
legends = 6
legendc = colg

    
# static def----------------------------------------------------------------

epsr = np.radians(23.437) # obliquity of ecliptic; hardcoded for 2020

##Month = np.array(["IV", "V", "VI", "VII", "VIII", "IX", 
##                  "X", "XI", "XII", "I", "II", "III"])
##Day = ['Jan', '6', '11', '16', '21', '26', '31', 'Feb', '10', '15',
##    '20', '25', 'March', '7', '12', '17', '22', '27', 'April', '6', '11',
##    '16', '21', '26', 'May', '6', '11', '16', '21', '26', '31', 'June',
##    '10', '15', '20', '25', '30', 'July', '10', '15', '20', '25', '30',
##    'Aug', '9', '14', '19', '24', '29', 'Sept', '8', '13', '18', '23',
##    '28', 'Oct', '8', '13', '18', '23', '28', 'Nov', '7', '12', '17',
##    '22', '27', 'Dec', '7', '12', '17', '22', '27']

#- sun calendar

SunDay = np.array(['1Jan','08','15','22','29','5Feb','12','19','26','5Mar','12','19','26','2Apr','09','16','23','30','7May','14','21','28','4Jun','11','18','25','2Jul','09','16','23','30','6Aug','13','20','27','3Sep','10','17','24','1Oct','08','15','22','29','5Nov','12','19','26','3Dec','10','17','24'])

SunJD = np.array([ 2459215.5, 2459222.5, 2459229.5, 2459236.5, 2459243.5, 2459250.5, 2459257.5, 2459264.5, 2459271.5, 2459278.5, 2459285.5, 2459292.5, 2459299.5, 2459306.5, 2459313.5, 2459320.5, 2459327.5, 2459334.5, 2459341.5, 2459348.5, 2459355.5, 2459362.5, 2459369.5, 2459376.5, 2459383.5, 2459390.5, 2459397.5, 2459404.5, 2459411.5, 2459418.5, 2459425.5, 2459432.5, 2459439.5, 2459446.5, 2459453.5, 2459460.5, 2459467.5, 2459474.5, 2459481.5, 2459488.5, 2459495.5, 2459502.5, 2459509.5, 2459516.5, 2459523.5, 2459530.5, 2459537.5, 2459544.5, 2459551.5, 2459558.5, 2459565.5, 2459572.5])

SunStep = SunJD[2]-SunJD[1]
print( 'sunstep', SunStep)

hours = np.arange( 0, 24)
deg = np.linspace( -180, 180, 361, endpoint=True)
rad = np.radians(deg)



#=== SETUP

#--- hemisphere
if xlat > 0:
    isouth = -1
else:
    isouth = 1  # originally Astlb was for southern hemisphere

#--- northlim of the sky
# limit in Dec on the other side of equator
# (originally for southern hemisphere, generalized)

if stereo:
    northlim = 25. * isouth # we limit to tropic in case of stereo projectoin
    if isouth > 0:
        rete = 'stereoS'
    else:
        rete = 'stereoN'
else:
    northlim = np.rint((95.*isouth + xlat)/10.)*10 # rounded horizon to 10deg
    #northlim = 90.*isouth + xlat #horizon
    #northlim = 60.*isouth + xlat #airmass =2

    if np.abs(xlat) < 30:
        northlim = 75.*isouth
        if isouth > 0:
            rete = 'polarS'
        else:
            rete = 'polarN'
    else:
        northlim = 50.*isouth
        if isouth > 0:
            rete = 'polarSS'
        else:
            rete = 'polarNN'

    
if  northlim >= 90.:
    northlim =89.
if  northlim <= -90.:
    northlim =-89.


# compute the scale normalization 
# so that outring is always = 180.

d2rnorm = 1.
w = d2r( northlim*isouth) # compute normalizatoin
d2rnorm = 180./w #<- normalization factor set.
#print( 'd2rnorm', d2rnorm)
outring = d2r( northlim*isouth) # appliescd wo normalization

print( 'rete:' , rete)
print( 'northlim:', northlim)
print( 'outring:', outring)

#------------------------------------------------------------------------------
# longitude latitude label
if xlong < 0:
    slong = 'E'
    sxlong = -xlong
else:
    slong = 'W'
    sxlong = xlong
if xlat < 0:
    slat = 'S'
    sxlat = -xlat
else:
    slat = 'N'
    sxlat = xlat

longlatlabel = '{:5.2f}{} {:5.2f}{}'.format(sxlong, slong,sxlat,slat)
pdfuniquelabel = '{:03d}{}{:02d}{}'.format(int(sxlong), slong,int(sxlat),slat)


## compute LST at 0UT 
# GMST calc computed ref to epoch 2000.0
# from B6 of the 2000 Astronomical Almanac
# GST_zero = Time('2000-09-20T00:03:14',format='isot')
#    LST Greenwhich = 00:00

JAN01 = Time('2000-01-01T00:00:00',format='isot')
tu = (JAN01.jd - 2451545.0)/36525.0 # days*100

# gms: greenwich ST on JAN01
gms = ( 24110.54841 + 8640184.812866 * tu 
    + 0.093104       * tu**2 
    - 6.2e-6         * tu**3 ) # seconds
gmst = (gms % 86400)/3600. # hours

if (gmst < 0. ):
    gmst += 24.

# xst: LST on observatory on UT_ref=jan1 - in h
xst = gmst - xlong/15.
print( "Sidereal time:")
print( "  GMST: {:.3f}".format(gmst))
print( "  long {:.2f}, longH: {:.2f}".format( xlong, xlong/15.))
print( "  xST: LST on observatory on Jan.1.0: {:.3f}".format(xst))


    
#== do the work
retename = 'astrolabe_'+rete+olab+'_rete.pdf'
#retename = 'astrolabe_'+rete+olab+'_rete.svg'
if doRete:
    print( "Rete", retename)
    fig = plt.figure(figsize=(plotsize,plotsize))
    plt.subplots_adjust(left=0.025, right=.975, top=.975, bottom=0.025)
    ax = fig.add_subplot(111, projection='polar')
    plot_rete()
    plt.savefig(outpath+retename)

    
platename = 'astrolabe_'+pdfuniquelabel+rete+olab+'_plate.pdf'
#platename = 'astrolabe_'+pdfuniquelabel+rete+olab+'_plate.svg'
if doPlate:
    print( "Plate", platename)
    fig = plt.figure(figsize=(plotsize,plotsize))
    plt.subplots_adjust(left=0.025, right=.975, top=.975, bottom=0.025)
    ax = fig.add_subplot(111, projection='polar')
    plot_plate()
    plt.savefig(outpath+platename)

rulename = 'astrolabe_'+rete+'_rule.pdf'
if doRule:
    print( "Rule", rulename)
    fig = plt.figure(figsize=(plotsize,plotsize))
    plt.subplots_adjust(left=0.025, right=.975, top=.975, bottom=0.025)
    ax = fig.add_subplot(111, projection='polar')
    plot_rule()
    plt.savefig(outpath+rulename)


#== create the webpage snippet
print( '</font></pre><table>')
print('<tr><td><a href="'+pubweb+retename+
       '"><img src="'+astweb+
       'astrolabe_grid_s.gif"></a>'+
       '<td><a href="'+pubweb+retename+
       '">Star overlay (rete)</a>')

print ('<tr><td><a href="'+pubweb+platename+
       '"><img src="'+astweb+
       'astrolabe_over_s.gif"></a>'+
       '<td><a href="'+pubweb+platename+
       '">Background plate (and mater)</a>')


print ('<tr><td><a href="'+pubweb+rulename+
       '"><img width=160 src="'+astweb+
       'astrolabe_rule_s.gif"></a>'+
       '<td><a href="'+pubweb+rulename+
       '">Optional ruler</a>')

print( '</table>')

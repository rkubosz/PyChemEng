#!/usr/bin/env python
# distutils: language = c++
# cython: profile=True

import matplotlib.pyplot as plt
import numpy
def PonchonSavaritGraph(vlefunc, hliquid, hvapour, xbottom, ytop, xfeed, hfeed):
    xvals = numpy.linspace(0.00001, 0.999999, 50)
    #plot tie lines
    for x in xvals:
        y = vlefunc(x)
        plt.plot([x,y],[hliquid(x), hvapour(y)], 'b--')
    #plot the h lines
    plt.plot(xvals,hliquid(xvals), 'b')
    plt.plot(vlefunc(xvals),hvapour(vlefunc(xvals)), 'r')

    #Plot the key concentrations
    plt.plot([xfeed, xbottom, ytop],[hfeed, hliquid(xbottom), hvapour(ytop)], '.k')

    #annotate the feed point
    plt.annotate("$x_F$", xytext=(30, 0), bbox=dict(boxstyle="round", fc=(1.0, 1.0, 1.0), ec=(1., .5, .5)), arrowprops=dict(arrowstyle="->"), textcoords='offset points', xy=(xfeed, hfeed))
    #annotate the bottoms point
    plt.annotate("$x_W$", xytext=(0, -30), bbox=dict(boxstyle="round", fc=(1.0, 1.0, 1.0), ec=(1., .5, .5)), arrowprops=dict(arrowstyle="->"), textcoords='offset points', xy=(xbottom, hliquid(xbottom)))
    #annotate the top point
    plt.annotate("$x_D$", xytext=(0, 30), bbox=dict(boxstyle="round", fc=(1.0, 1.0, 1.0), ec=(1., .5, .5)), arrowprops=dict(arrowstyle="->"), textcoords='offset points', xy=(ytop, hvapour(ytop)))


def PonchonSavaritMinTrays(vlefunc, hliquid, hvapour, xbottom, ytop, xfeed, hfeed):
    PonchonSavaritGraph(vlefunc, hliquid, hvapour, xbottom, ytop, xfeed, hfeed)
    #Start from the bottom
    x=xbottom
    trays=0
    for iteration in range(100):
        y=vlefunc(x)
        trays+=1
        plt.plot([x,y],[hliquid(x),hvapour(y)], 'b')
        plt.plot([y,y],[hliquid(y),hvapour(y)], color='#367588')
        plt.annotate(str(trays), xy=(y,hvapour(y)))
        if y > ytop:
            plt.title("Minimum trays="+str(trays))
            return trays
        x=y
    raise Exception("More than 100 stages needed!")

def PonchonSavaritFindX(vlefunc, hliquid, hvapour,  xtarget, htarget):
    from scipy.optimize import brentq
    #search for the liquid x where its VLE tie line passes through the target point
    func = lambda x : numpy.poly1d(numpy.polyfit([x,vlefunc(x)], [hliquid(x), hvapour(vlefunc(x))], 1))(xtarget) - htarget
    return brentq(func, 0.00001, 0.99999)

def PonchonSavaritMinReflux(vlefunc, hliquid, hvapour, xbottom, ytop, xfeed, hfeed, n=20):
    import numpy
    PonchonSavaritGraph(vlefunc, hliquid, hvapour, xbottom, ytop, xfeed, hfeed)
    xtop = PonchonSavaritFindX(vlefunc, hliquid, hvapour, ytop, hvapour(ytop))
    xfeed_liquid =PonchonSavaritFindX(vlefunc, hliquid, hvapour, xfeed, hfeed)

    #Process the bottom points
    Prh=hliquid(xbottom)
    for x in numpy.linspace(xbottom, xfeed_liquid, n):
        y=vlefunc(x)
        #Grab the tie-line
        tieline = numpy.poly1d(numpy.polyfit([x,y], [hliquid(x), hvapour(y)], 1))
        #Extrapolate it to the PR
        Prh=min(Prh, tieline(xbottom))
        plt.plot([xbottom, y], [tieline(xbottom), tieline(y)], color='#367588')

    plt.annotate("$P_R^{(min)}$", xytext=(0, 30), bbox=dict(boxstyle="round", fc=(1.0, 1.0, 1.0), ec=(1., .5, .5)), arrowprops=dict(arrowstyle="->"), textcoords='offset points', xy=(xbottom, Prh))
    plt.plot([xbottom, xbottom], [hliquid(xbottom), Prh], color='#367588')

    #Process the top points
    Pch=hvapour(ytop)
    for x in numpy.linspace(xfeed_liquid, xtop, n):
        y=vlefunc(x)
        #Grab the tie-line
        tieline = numpy.poly1d(numpy.polyfit([x,y], [hliquid(x), hvapour(y)], 1))
        #Extrapolate it to the PR
        Pch=max(Pch, tieline(ytop))
        plt.plot([ytop, y], [tieline(ytop), tieline(y)], color='#367588')

    plt.annotate("$P_C^{(min)}$", xytext=(0, 30), bbox=dict(boxstyle="round", fc=(1.0, 1.0, 1.0), ec=(1., .5, .5)), arrowprops=dict(arrowstyle="->"), textcoords='offset points', xy=(ytop, Pch))
    plt.plot([ytop, ytop], [hvapour(ytop), Pch], color='#367588')

    #Now calculate if Pch is higher due to Prh.
    Pch_bottom = numpy.poly1d(numpy.polyfit([xbottom,xfeed], [Prh, hfeed], 1))(ytop)
    plt.plot([xbottom, ytop], [Prh, Pch_bottom], color='#367588')
    
    Pch = max(Pch, Pch_bottom)

    #Calculate the min R from the top section
    Rmin = (Pch-hliquid(ytop)) / (hvapour(ytop)-hliquid(ytop)) - 1
    plt.title("Minimum reflux ratio = "+str(Rmin))
    return Rmin

def PonchonSavarit(vlefunc, hliquid, hvapour, xbottom, ytop, xfeed, hfeed, R, nsteps_max=50):
    import numpy
    PonchonSavaritGraph(vlefunc, hliquid, hvapour, xbottom, ytop, xfeed, hfeed)
    Pch = (R+1) * (hvapour(ytop) - hliquid(ytop)) + hliquid(ytop)
    #Calculate Prh by taking a line through Pch and xfeed/hfeed
    Prh = numpy.poly1d(numpy.polyfit([xfeed,ytop], [hfeed,Pch], 1))(xbottom)

    plt.annotate("$P_R$", xytext=(0, 30), bbox=dict(boxstyle="round", fc=(1.0, 1.0, 1.0), ec=(1., .5, .5)), arrowprops=dict(arrowstyle="->"), textcoords='offset points', xy=(xbottom, Prh))
    plt.annotate("$P_C$", xytext=(0, 30), bbox=dict(boxstyle="round", fc=(1.0, 1.0, 1.0), ec=(1., .5, .5)), arrowprops=dict(arrowstyle="->"), textcoords='offset points', xy=(ytop, Pch))
    plt.plot([xbottom, xbottom], [hliquid(xbottom), Prh], color='#367588')
    plt.plot([ytop, ytop], [hvapour(ytop), Pch], color='#367588')
    plt.plot([xbottom, ytop], [Prh, Pch], color='#367588')

    from scipy.optimize import brentq
    ymaxfeedline = brentq(lambda x : numpy.poly1d(numpy.polyfit([xbottom,ytop], [Prh,Pch], 1))(x) - hvapour(x), xbottom, ytop)    
    
    #Start by calculating the lower tray count
    trays=[]
    nsteps=0
    x = xbottom
    while nsteps<nsteps_max:
        #Calculate the produced vapour
        y = vlefunc(x)
        trays.append((x,y))
        plt.plot([x, y], [hliquid(x), hvapour(y)], color='b')
        plt.annotate(str(len(trays)), xy=(y,hvapour(y)))
        #Check if we've crossed into the enrichment section:
        if y > ymaxfeedline:
            break
        #Fit an operating line between (xbottom, Prh) and (y, hvapour(y))
        op_line=numpy.poly1d(numpy.polyfit([xbottom,y], [Prh,hvapour(y)], 1))
        plt.plot([xbottom, y], [Prh, hvapour(y)], color='#367588')
        #Now solve for the x of the next tray by looking for the intercept of the hliquid and op line
        x=brentq(lambda x : op_line(x) - hliquid(x), xbottom, ytop)
    if nsteps >= nsteps_max:
        raise Exception("Max number of trays exceeded")

    y = vlefunc(x)
    while nsteps<nsteps_max:
        #Calculate the produced vapour of the last tray
        #Check if this is the end?
        if y > ytop:
            break            
        #Fit an operating line between (ytop, Pch) and (y, hvapour(y))
        op_line=numpy.poly1d(numpy.polyfit([ytop,y], [Pch,hvapour(y)], 1))
        #Calculate the new x
        x=brentq(lambda x : op_line(x) - hliquid(x), xbottom, ytop)
        plt.plot([ytop, x], [Pch, hliquid(x)], color='#367588')
        #Now solve for the x of the next tray by looking for the intercept of the hliquid and op line
        y=vlefunc(x)
        plt.plot([x, y], [hliquid(x), hvapour(y)], color='b')
        trays.append((x,y))
        plt.annotate(str(len(trays)), xy=(y,hvapour(y)))

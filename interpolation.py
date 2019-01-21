#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 17:05:29 2018
@author: Joe
"""

from tkinter import *
from math import *
import numpy
from scipy import stats

root = Tk()

################################################################################

def grabColumn(fileName,listName,index):
    file = open(fileName, 'r')
    text = file.readlines()
    file.close()
    for line in text:
        if line[0] == '2' and line.split()[1] != '03Labyrinth2':
            value = float(line.split()[index])
            listName.append(value)
        else:
            continue  
    #print(listName)
    return(listName)

################################################################################

gridXvals = []
gridYvals = []
gridXvals2 = []
gridYvals2 = []
def makeGrid(xVals,yVals,resolution):
    interval = (max(xVals)-min(xVals)) / resolution
    xVal = min(xVals)
    yVal = min(yVals)
    while xVal < max(xVals):
        gridXvals.append(xVal)
        xVal += interval*4
    while yVal < max(yVals):
        gridYvals.append(yVal)
        yVal += interval
    #print(gridXvals,gridYvals)
    for gridXval in gridXvals:
        for gridYval in gridYvals:
            gridYvals2.append(gridYval)
    for gridYval in gridYvals:
        for gridXval in gridXvals:
            gridXvals2.append(gridXval)
            
zValsList = []
def createZvalsKern(gridXvals2,gridYvals2,dataXvals,dataYvals,searchRadius):
    for gridX,gridY in zip(gridXvals2,gridYvals2):
        score = 0.0
        for dataX,dataY in zip(dataXvals,dataYvals):
            pi = 3.1415926
            distance = acos( (sin(gridY*pi/180)*sin(dataY*pi/180)) \
                            + (cos(gridY*pi/180)*cos(dataY*pi/180) \
                               *cos((((gridX-dataX)*pi/180)**2)**0.5)) )
            if distance < searchRadius:
                kern = (1.0/(searchRadius*len(dataXvals)))*((distance/(searchRadius)))
                score += kern
            else:
                continue
        zValsList.append(score)
    #print(zValsList)
    
# IDW algorithm from: https://www.e-education.psu.edu/geog486/node/1877
def createZvalsIDW(gridXvals2,gridYvals2,dataXvals,dataYvals,dataZvals,neighborhood):
    for gridX,gridY in zip(gridXvals2,gridYvals2):
        score = 0.0
        numer = 0.0
        denom = 0.0
        for dataX,dataY,dataZ in zip(dataXvals,dataYvals,dataZvals):
            pi = 3.1415926
            distance = acos( (sin(gridY*pi/180)*sin(dataY*pi/180)) \
                            + (cos(gridY*pi/180)*cos(dataY*pi/180) \
                               *cos((((gridX-dataX)*pi/180)**2)**0.5)) )
            if distance < neighborhood:
                weight = ( 1 / ( distance**2 ) )
                numer += ( weight * dataZ )
                denom += weight
            else:
                continue
        score = numer / denom
        zValsList.append(score)
    return(zValsList)
    
colors = []
def createColors(zValsList):
    ramp = ['#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58']
    print('max = ',max(zValsList))
    print('min = ',min(zValsList))
    print('range = ',max(zValsList)-min(zValsList))
    transform = (max(zValsList)-min(zValsList))/8
    print('transform = ',transform)
    for zVal in zValsList:
        colorIndex = int( round( ( zVal - min(zValsList) ) / transform ) ) 
        #print(colorIndex)
        color = ramp[colorIndex]
        colors.append(color)
    textValue = min(zValsList)
    rampLoc = canheight*0.8
    #for color in ramp:
    #    can1.create_text(canwidth*0.95,rampLoc, font="Arial 11",text=str(round(textValue,2)))
    #    box = [canwidth*0.97,rampLoc,canwidth*0.99,rampLoc,canwidth*0.99,rampLoc+1,canwidth*0.97,rampLoc+1]
    #    can1.create_line(box,fill=color)
    #    rampLoc -= (canheight*0.6/len(ramp))
    #    textValue += transform
    return(colors)
        
################################################################################

def pearsons(first,second):
    firstMean =numpy.mean(first)
    secondMean = numpy.mean(second)
    i = 0
    s = 0.0
    t = 0.0
    u = 0.0
    while i < len(first):
        s += ((first[i] - firstMean)*(second[i] - secondMean))
        t += ((first[i] - firstMean)**2)
        u += ((second[i] - secondMean)**2)
        i += 1
    r = (s)/((t**0.5)*(u**0.5))
    print('The Pearsons correlation coefficient is ',str(r))
    return(r)

def scatterPlotWithRegres(indepVarSet,depVarSet,indepName,depName,title):
    
    #xaxislabel = can1.create_text((canwidth/2),(yAxisStart+25),\
    #                     font="Arial 10",text=indepName)
    #yaxislabel = can1.create_text((25),(canheight/2),\
    #                     font="Arial 10",text=depName)
    
    indepMax = max(indepVarSet)
    indepMin = min(indepVarSet)     
    indepRange = indepMax - indepMin
    indepMean = numpy.mean(indepVarSet)
    
    depMax = max(depVarSet)
    depMin = min(depVarSet)     
    depRange = depMax - depMin
    depMean = numpy.mean(depVarSet)    
        
    u = 0.0
    v = 0.0
    w = 0.0
    x = 0.0
    y = 0.0
    z = 0.0
    n = len(indepVarSet)

    for ind,dep in zip(indepVarSet,depVarSet):
        u += dep
        v += ( ind ** 2 )
        w += ind
        x += (ind * dep )
        y += ( ind ** 2 )
        z += ind
    a = ( ( u * v ) - ( w * x ) ) / ( ( n * y ) - z ** 2 )
    b = ( ( n * x ) - (w * u ) ) / ( ( n * y ) - z ** 2 )
    slope, intercept, r_value, p_value, std_err = stats.linregress(indepVarSet,depVarSet)
    equationstring = 'y = {0}x + {1}\nr squared = {2}'.format(round(b,3),round(a,3),round(r_value**2,3))
    print(equationstring)
    regStart = ( ( indepMin * b ) + a )
    regEnd = ( ( indepMax * b ) + a )
    plotRegXStart =  xAxisStart + ( ((indepMin-indepMin)*(xAxisLength))/indepRange )
    plotRegYStart =  yAxisStart - ( ((regStart-depMin)*(yAxisLength))/depRange )  
    plotRegXEnd =  xAxisStart + ( ((indepMax-indepMin)*(xAxisLength))/indepRange )
    plotRegYEnd =  yAxisStart - ( ((regEnd-depMin)*(yAxisLength))/depRange ) 
    plotRegresLine = [plotRegXStart,plotRegYStart,plotRegXEnd,plotRegYEnd]
    #can1.create_line(plotRegresLine,fill='firebrick4',dash=(3, 4))
    #can1.create_text(xAxisEnd-(xAxisLength*1/10),(canheight/2),\
    #                 font="Arial 11", fill = 'firebrick4',text=equationstring)
    
    for ind,dep in zip(indepVarSet,depVarSet):    
        plotX1 =  xAxisStart + ( ((ind-indepMin)*(xAxisLength))/indepRange )
        plotY1 =  yAxisStart - ( ((dep-depMin)*(yAxisLength))/depRange )
        plotX = 400 + ( 0.866 * plotX1 ) - ( 0.866 * plotY1 )
        plotY = 900 - ( 0.500 * plotX1 ) - ( 0.500 * plotY1 ) 
        point = [plotX,plotY-1000,plotX,plotY+2000,plotX,plotY-1000]
        can1.create_line(point,fill='gray90')
        point = [plotX-2,plotY-2,plotX-2,plotY+2,plotX+2,plotY+2,plotX+2,plotY-2,plotX-2,plotY-2]
        can1.create_line(point,fill='gray75')
        
    j = xAxisStart
    jInterval = (indepRange) / xAxisLength 
    k = indepMin
    i = 0
    #while j < xAxisEnd:
        #if i % 150 == 0:
            #can1.create_text(j+7,(yAxisStart+8), font="Arial 10",text=str(round(k,2)))
            #can1.create_line( j, (yAxisStart-3), j, (yAxisStart+3) )
        #else:
        #    k = k
        #j += 1
        #k += jInterval
        #i += 1
        
    j = yAxisStart
    jInterval = yAxisLength / 10
    k = depMin
    kInterval = depRange / 10
    #while j >= yAxisEnd:
        #can1.create_line( (xAxisStart-3), j, (xAxisStart+3), j )
        #can1.create_text((xAxisStart-20),j, font="Arial 10",text=str(round(k,2)))
        #j -= jInterval
        #k += kInterval

    title = can1.create_text((canwidth/2),(yAxisEnd*0.4),\
                         font="Arial 15 bold",\
                         text=title )
    
    #can1.create_line( xAxisStart, yAxisStart, xAxisEnd, yAxisStart )
    #can1.create_line( xAxisStart, yAxisStart, xAxisStart, yAxisEnd )
    
def gridPlot(indepVarSet,depVarSet,Zs,colors,indepName,depName,title):
    
    #xaxislabel = can1.create_text((canwidth/2),(yAxisStart+25),\
    #                     font="Arial 12",text=indepName)
    #yaxislabel = can1.create_text((10),(canheight/2),\
    #                     font="Arial 12",text=depName)
    
    indepMax = max(indepVarSet)
    indepMin = min(indepVarSet)     
    indepRange = indepMax - indepMin
    indepMean = numpy.mean(indepVarSet)
    
    depMax = max(depVarSet)
    depMin = min(depVarSet)     
    depRange = depMax - depMin
    depMean = numpy.mean(depVarSet)    
    
    for ind,dep,color,Z in zip(indepVarSet,depVarSet,colors,Zs):    
        plotX1 =  xAxisStart + ( ((ind-indepMin)*(xAxisLength))/indepRange )
        plotY1 =  yAxisStart - ( ((dep-depMin)*(yAxisLength))/depRange )
        plotX = -300 + ( 1.35 * plotX1 ) + ( 0.3 * plotY1 )  #0.866
        plotY = 1800 - ( -0 * plotX1 ) + ( 0.5 * plotY1 ) + ( -220 * (Z**(0.9)) ) #0.5  ############################################
        size = 0.5+4*((Z-10)**2)
        point = [plotX,plotY-size,plotX,plotY,plotX+size,plotY,plotX+size,plotY-size,plotX,plotY-size]
        #point = [plotX-1,plotY-1,plotX+1,plotY-1,plotX+1,plotY+1,\
        #         plotX-1,plotY+1,plotX-1,plotY-1]
        can1.create_line(point,fill=color)
        
    #j = xAxisStart
    #jInterval = (indepRange) / xAxisLength 
    #k = indepMin
    #i = 0
    #while j < xAxisEnd:
    #    if i % 150 == 0:
    #        can1.create_text(j+7,(yAxisStart+8), font="Arial 12",text=str(round(k,2)))
    #        can1.create_line( j, (yAxisStart-3), j, (yAxisStart+3) )
    #    else:
    #        k = k
    #    j += 1
    #    k += jInterval
    #    i += 1
        
    j = yAxisStart
    jInterval = yAxisLength / 10
    k = depMin
    kInterval = depRange / 10
    #while j >= yAxisEnd:
        #can1.create_line( (xAxisStart-3), j, (xAxisStart+3), j )
        #can1.create_text((xAxisStart-20),j, font="Arial 12",text=str(round(k,2)))
        #j -= jInterval
        #k += kInterval

    #title = can1.create_text((canwidth/2),(yAxisEnd*0.4),\
                         #font="Arial 16 bold",\
                         #text=title )
    
    #can1.create_line( xAxisStart, yAxisStart, xAxisEnd, yAxisStart )
    #can1.create_line( xAxisStart, yAxisStart, xAxisStart, yAxisEnd )

################################################################################    

lons = []
lats = []
soilMoisture = []
chlorophyllA = []
organicCarbon = []
soilCond = []
soilSalinity = []
soilpH = []
elevation = []
grabColumn('DryValleyData.txt',lons,2)
grabColumn('DryValleyData.txt',lats,3)
grabColumn('DryValleyData.txt',soilMoisture,4)
grabColumn('DryValleyData.txt',chlorophyllA,5)
grabColumn('DryValleyData.txt',organicCarbon,6)
grabColumn('DryValleyData.txt',soilCond,7)
grabColumn('DryValleyData.txt',soilSalinity,8)
grabColumn('DryValleyData.txt',soilpH,9)
grabColumn('DryValleyData.txt',elevation,10)

print(elevation)

species = [8,4,3,3,5,3,3,4,7,7,9,11,7,10]

habitabilitys = []
for moisture,chloro,carbon,conductivity,salinity,ph,elev \
in zip(soilMoisture,chlorophyllA,organicCarbon,soilCond,soilSalinity,soilpH,elevation):
    habitability = (moisture*chloro*carbon*ph)/(salinity*conductivity*elev+0.0000001)
    habitabilitys.append(habitability)

weightHab = []
for moisture,chloro,carbon,conductivity,salinity,ph,elev \
in zip(soilMoisture,chlorophyllA,organicCarbon,soilCond,soilSalinity,soilpH,elevation):
    habitability = 1.274*0.565*moisture/max(soilMoisture)\
                    +0.062*0.123*chloro/max(chlorophyllA)\
                    +2.175*0.079*carbon/max(organicCarbon)\
                    +2.82*0.368*ph/max(soilpH)\
                    -1.437*0.193*salinity/max(soilSalinity)\
                    -0.832*0.198*conductivity/max(soilCond)\
                    -0.013*0.487*elev/max(elevation)
    weightHab.append(habitability)
print(weightHab)

weightHab2 = []
for moisture,chloro,carbon,conductivity,salinity,ph,elev \
in zip(soilMoisture,chlorophyllA,organicCarbon,soilCond,soilSalinity,soilpH,elevation):
    habitability = pearsons(species,soilMoisture)*moisture/max(soilMoisture)\
                    +pearsons(species,chlorophyllA)*chloro/max(chlorophyllA)\
                    +pearsons(species,organicCarbon)*carbon/max(organicCarbon)\
                    +pearsons(species,soilpH)*ph/max(soilpH)\
                    +pearsons(species,soilSalinity)*salinity/max(soilSalinity)\
                    +pearsons(species,soilCond)*conductivity/max(soilCond)\
                    +pearsons(species,elevation)*elev/max(elevation)
    weightHab2.append(habitability)
print(weightHab2)

canheight = 800
ratio = (((max(lons)-min(lons)) / (max(lats)-min(lats))**2)**0.5)*0.2
canwidth = canheight * ratio
xAxisFactor = 0.8
yAxisFactor = 0.8

can1 = Canvas(root, width = canwidth, height = canheight )
xAxisStart = canwidth - ( canwidth * xAxisFactor )
xAxisEnd = canwidth * xAxisFactor 
xAxisLength = xAxisEnd - xAxisStart
yAxisStart = canheight * yAxisFactor 
yAxisEnd =  canheight - ( canheight * yAxisFactor )
yAxisLength =  yAxisStart - yAxisEnd

#print(len(elevation))
#print(len(species))


makeGrid(lons,lats,600)
#print(elevation)
#createZvalsKern(gridXvals2,gridYvals2,lons,lats,0.0005)
#moistures = createZvalsIDW(gridXvals2,gridYvals2,lons,lats,soilMoisture,0.005)
#chloros = createZvalsIDW(gridXvals2,gridYvals2,lons,lats,chlorophyllA,0.005)
#carbons = createZvalsIDW(gridXvals2,gridYvals2,lons,lats,organicCarbon,0.005)
#conductivitys = createZvalsIDW(gridXvals2,gridYvals2,lons,lats,soilCond,0.005)
#salinitys = createZvalsIDW(gridXvals2,gridYvals2,lons,lats,soilSalinity,0.005)
phs = createZvalsIDW(gridXvals2,gridYvals2,lons,lats,soilpH,0.005)
#elevations = createZvalsIDW(gridXvals2,gridYvals2,lons,lats,elevation,0.02)
#biodiv = createZvalsIDW(gridXvals2,gridYvals2,lons,lats,species,0.005)
#createZvalsIDW(gridXvals2,gridYvals2,lons,lats,weightHab,0.005)
print(zValsList)
createColors(zValsList)
#print(colors)
gridPlot(gridXvals2,gridYvals2,zValsList,colors,'Lon','Lat','Antarctic Dry Valley Expected Habitability v2')

#scatterPlotWithRegres(lons,lats,\
#' ',' ',\
#' ')

can1.pack()
root.mainloop()

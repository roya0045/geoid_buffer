# -*- coding: utf-8 -*-

"""
/***************************************************************************
 EllipsoidBuffer
                                 A QGIS plugin
 Generate buffer in ellipsoid
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2023-05-26
        copyright            : (C) 2023 by Alex RL
        email                : contact on github
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'Alex RL'
__date__ = '2023-05-26'
__copyright__ = '(C) 2023 by Alex RL'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,QgsCoordinateTransform,QgsCoordinate,QgsCoordinateReferenceSystems,QgsGeometry,QgsGeometryCollection,
                       QgsFeatureSink,QgsProcessingParameterDistance,QgsProcessingParameterEnum,QgsProcessingParameterBoolean,
                       QgsProcessingException,QgsProcessingParameterCrs,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterFeatureSink)
from qgis import processing
from pyproj.geod import Geod #https://pyproj4.github.io/pyproj/stable/api/geod.html#pyproj.Geod.fwd
from geographiclib.geodesic import Geodesic #if fail run install script?


#//https://sourceforge.net/p/saga-gis/code/ci/master/tree/saga-gis/src/tools/shapes/shapes_tools/shapes_buffer.cpp

try:
	polyt= Qgis.WkbType.PolygonGeomettry
Except:
	polyt=QgsWkbTypes.Polygon

def buffer (geometry,distancem,srcCrs,destCrs,dissolve=True,flatcap=False):
	fwdtrsctx = QgsCoordinateTransform(srcCrs,destCrs )
	revtrsctx = QgsCoordinateTransform(destCrs,srcCrs )
	result=geometry.transform(fwdtrsctx)
    if wgs84:
        geodes = Geod()#Geodesic.WGS84
    else: #https://proj.org/en/9.2/usage/ellipsoids.html
        equatrad = 
        flattening = 
        geodes = Geod()#Geodesic(equatrad,flattening) #Geodesic(6378388, 1/297.0)
	if not(result):
		raise QgsException("Failed to transform")
	buffered=_buffer(geometry.asGeometryCollection(),distancem,geodes,flatcap)
	
	if ( dissolve):
		buffered= QgsGeometry.unaryUnion(buffered)
	result= buffered.transform(revtrsctx)
	if not(result):
		raise QgsExcepiton("could not transform back resulting geometry")
	return(buffered)

def _buffer (geometry,distancem:float,geoid:Geod,flat:Bool):
	if len(geometry>1):
		buffered_coll=[]
		for geom in  geometry:
			buffered=_buffer( geom,distancem,flat)
			buffered_coll.append(buffered)
		return(buffered_coll)
	geometry=geometry[0] //QgsGeometry
	if geometry.isMultipart():
		buffered_coll=[]
		for part in  geometry.parts():
			buffered=_buffer( part,distancem,flat)
			buffered_coll.append(buffered)
		return(buffered_coll)
	previousVertex=None
	previousAz = None
	buffered=list()
	for ix,vertex in enumerate(geometry.vertices()):

		if (previousVertex is None):
			continue
		newbuff=buff_line(peviousVertex,vertex,distancem,geoid,flatstart = flat and ix ==1,flatend=flat)
		buffered.append(newbuff)
		previousVertex=vertex

	if (geometry.wkbType() == polyt):
		if ( previousVertex != geoemtry.vertexAt(0) ):
			buffered.append(buffLine(previousVertex, geometry.vertexAt(0),distancem,geoid))
		buffered = QgsGeometry.unaryUnion(buffered)
		#//check outside and inside range? process as line and merge with polygon, if buffer is negative?
		#// buffer as line
		if distance < 0.0:
			buffered=geometry.difference(buffered)
		else:
			buffered= QgsGeometry.unaryUnion([buffered,geometry])
	elif previousVertex is None: // point
		buffered=make_arcs(srcPnt,distance)
		#make points at given interval/precision
	else:
		buffered =  QgsGeometry.unaryUnion(buffered)

	return(buffered) //need to reproject

def buff_line(p1,p2,distance,geoid:Geod,flatstart:Bool=False,flatend:Bool=False):
	az=geoid.inv(lat1, lon1, lat2, lon2,caps=512)[0]
	lim=abs(az)+90
	if flatend:
		precision = 180
	startarc= make_arc(p1,distance,geoid,lim,lim+180,180 if flatstart else precision)
	endarc= make_arc(p2,distance,geoid,lim,lim-180,180 if flatend else precision)
	#//join arcs
	return(polygon(startarc+endarc+startarc[0] ))
	

 #//https://geographiclib.sourceforge.io/Python/doc/code.html#geographiclib.geodesic.Geodesic.Direct
def make_arc(srcPnt,distance,geoid:Geod,start:float=0.o,end:float=360.0,precision:float=1.0):

	arc=[]
	for az in range(start,end+precision,precision):
		destPnt=geoid.fwd(srcPnt.lat,srcPnt.long,az,dist)
		arc.append(QgsPoint(destPnt.lat2,destPnt.long2))

	return(arc)


class EllipsoidBufferAlgorithm(QgsProcessingAlgorithm):
    """
    This is an example algorithm that takes a vector layer and
    creates a new identical one.

    It is meant to be used as an example of how to create your own
    algorithms and explain methods and variables used to do it. An
    algorithm like this will be available in all elements, and there
    is not need for additional work.

    All Processing algorithms should extend the QgsProcessingAlgorithm
    class.
    """

    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    OUTPUT = 'OUTPUT'
    DISTM = 'DISTM'
    ELLIPSOID = 'ELLIPSOID'
    INPUT = 'INPUT'
    DISSB = 'DISSB'
    ENDSTYLE = 'ENDSTYLE'
    Capstyle = ['Round','Flat']

    def initAlgorithm(self, config):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        # We add the input vector features source. It can have any kind of
        # geometry.
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr('Input layer'),
                [QgsProcessing.TypeVectorAnyGeometry]
            )
        )

        self.addParameter(QgsProcessingParameterDistance(self.DISTM,self.tr('Distance Meter'))
        )

        self.addParameter(QgsProcessingParameterEnum(self.ENDSTYLE,selt.tr('End style'),self.Capstyle,allowMultiple=False,
                                                      defaultValue=0))

        self.addParameter(QgsProcessingBoolean(self.DISSB,self.tr("Dissolve features"),False))

        self.addParameter(QgsProcessingParameterCrs(self.ELLIPSOID,self.tr("Ellipsoid to use"),defaultValue= ?
                        ,optional = True))

        # We add a feature sink in which to store our processed features (this
        # usually takes the form of a newly created vector layer when the
        # algorithm is run in QGIS).
        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                self.tr('Output layer')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """

        # Retrieve the feature source and sink. The 'dest_id' variable is used
        # to uniquely identify the feature sink, and must be included in the
        # dictionary returned by the processAlgorithm function.
        source = self.parameterAsSource(parameters, self.INPUT, context)
        (sink, dest_id) = self.parameterAsSink(parameters, self.OUTPUT,
                context, source.fields(), source.wkbType(), source.sourceCrs())

        distancem = self.parameterAsDouble(parameters, self.DISTM, context)
        interimCrs = self.parameterAsCrs(parameters, self.ELLIPSOID, context)
        capstyle = self.parameterAsEnum(parameters,self.ENDSTYLE,context)
        if capstyle == 1 :
            flatEnd =True
        else:
            flatEnd = False
        dissolveB = self.parameterAsBoolean(parameter,self.DISSB,context)
        # Compute the number of steps to display within the progress bar and
        # get features from source
        total = 100.0 / source.featureCount() if source.featureCount() else 0
        features = source.getFeatures()

        currentCrs =  source.sourceCrs()

        if interimCrs:
            Ellipsoid = interimCrs.ellipsoidAcronym()
        else:
            Ellipsoid = currentCrs.ellipsoidAcronym()

        for current, feature in enumerate(features):
            # Stop the algorithm if cancel button has been clicked
            if feedback.isCanceled():
                break
            oldgeometry=feature.geometry()
            buffered= buffer(oldgeometry,distancem,currentCrs, Ellipsoid,dissolveB,flatEnd)

            feature.setGeometry(buffered)

            # Add a feature in the sink
            sink.addFeature(feature, QgsFeatureSink.FastInsert)

            # Update the progress bar
            feedback.setProgress(int(current * total))

        # Return the results of the algorithm. In this case our only result is
        # the feature sink which contains the processed features, but some
        # algorithms may return multiple feature sinks, calculated numeric
        # statistics, etc. These should all be included in the returned
        # dictionary, with keys matching the feature corresponding parameter
        # or output names.
        return {self.OUTPUT: dest_id}

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Ellipsoid buffer'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr(self.name())

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr(self.groupId())

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Vector Geometry'

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return EllipsoidBufferAlgorithm()

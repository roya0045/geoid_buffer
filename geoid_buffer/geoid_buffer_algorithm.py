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

__author__ = "Alex RL"
__date__ = "2023-05-26"
__copyright__ = "(C) 2023 by Alex RL"

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = "$Format:%H$"

from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import (
    QgsProcessing,
    QgsCoordinateTransform,
    QgsCoordinateReferenceSystem,
    QgsGeometry,
    QgsGeometryCollection,
    QgsWkbTypes,
    QgsFeatureSink,
    QgsProcessingParameterNumber,
    QgsProcessingParameterEnum,
    QgsProcessingParameterBoolean,
    QgsProcessingException,
    QgsProcessingParameterCrs,
    QgsPropertyDefinition,
    QgsProject,
    QgsProcessingException,
    QgsMultiPolygon,
    QgsProcessingAlgorithm,
    QgsPointXY,
    QgsPolygon,
    QgsLineString,
    QgsPoint,
    QgsProcessingParameters,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSink,
)
from qgis import processing

# from geographiclib.geodesic import Geodesic
try:
    from pyproj import (
        Proj,
        Geod,
    )  # https://pyproj4.github.io/pyproj/stable/api/geod.html#pyproj.Geod.fwd
except:
    raise QgsProcessingException(
        "pyproj is not accessible, install the pyproj library via osgeo if available"
    )


# //https://sourceforge.net/p/saga-gis/code/ci/master/tree/saga-gis/src/tools/shapes/shapes_tools/shapes_buffer.cpp

parlist = ["a", "b", "k"]
geoparams = {"a": "a", "b": "b", "k": "f"}

try:
    polyt = Qgis.WkbType.PolygonGeomettry
    mpolyt = Qgis.WktType.MultiPolygon
except:
    polyt = QgsWkbTypes.Polygon
    mpolyt = QgsWkbTypes.MultiPolygon


def unigeom(geoms, feedback=None):
    return QgsGeometry.unaryUnion(geoms)
    # return(QgsGeometry.collectGeometry(geoms))


def pts2qgeom(ptslist, feedback, debug: bool = False):
    if debug:
        feedback.pushInfo(str(ptslist))
    multipoly = QgsMultiPolygon()
    if isinstance(ptslist[0], list):  # nested lists
        for ptys in ptslist:
            multipoly.addGeometry(QgsPolygon(QgsLineString(ptys)))
    else:
        multipoly.addGeometry(QgsPolygon(QgsLineString(ptslist)))
    return QgsGeometry.fromWkt(multipoly.asWkt())


def buffer(
    geometry,
    distancem,
    srcCrs,
    destCrs,
    feedback,
    dissolve=True,
    flatcap=False,
    precision=1.0,
):
    if distancem == 0.0:
        return geometry
    geoid = destCrs.toGeographicCrs()  # ellipsoidAcronym())
    if not (geoid.isValid):
        raise QgsProcessingException("invalid crs")
    fwdtrsctx = QgsCoordinateTransform(srcCrs, geoid, QgsProject.instance())
    revtrsctx = QgsCoordinateTransform(geoid, srcCrs, QgsProject.instance())
    result = geometry.transform(fwdtrsctx)
    if not (result or geometry.isGeosValid()):
        raise QgsProcessingException("Failed to transform")
    crsvars = dict(
        {tuple(i[1:].split("=")) for i in destCrs.toProj().split(" ") if "=" in i}
    )
    elkwg = dict(
        (geoparams.get(b), float(crsvars.get(b))) for b in parlist if b in crsvars
    )
    if "ellps" in crsvars:
        geodes = Geod(crsvars["ellps"])  # Geodesic.WGS84
    elif len(elkwg) == 0:
        geodes = Geod(ellps="WGS84")
    else:  # https://proj.org/en/9.2/usage/ellipsoids.html
        # feedback.pushInfo(str(elkwg))
        geodes = Geod(
            **elkwg
        )  # Geodesic(equatrad,flattening) #Geodesic(6378388, 1/297.0)

    buffered = _buffer(
        geometry.asGeometryCollection(), distancem, geodes, flatcap, feedback, precision
    )
    if not (isinstance(buffered, list)):
        # feedback.pushInfo(buffered.asWkt())
        buffered.transform(revtrsctx)
        if not (buffered.isGeosValid()):
            feedback.pushInfo("geom is invalid 93")
            buffered = buffered.makeValid()
        return buffered
    if dissolve:
        buffered = unigeom(buffered)
    geoms = list()
    # feedback.pushInfo(str(len(buffered)))
    for buff in buffered:
        # feedback.pushInfo(buff.asWkt())
        buff.transform(revtrsctx)
        if buff.isMultipart():
            ok = buff.convertToSingleType()
            # feedback.pushInfo(str(ok))
        # feedback.pushInfo(buff.asWkt())
        geoms.append(buff)

    retgeom = unigeom(geoms)
    # feedback.pushInfo(retgeom.asWkt())
    if not (retgeom.isGeosValid()):
        # feedback.pushInfo(retgeom.asWkt())
        feedback.pushInfo("geom is invalid 112")

        # feedback.pushInfo(retgeom.asWkt())
        retgeom = retgeom.makeValid()
        # raise QgsProcessingException("could not transform back resulting geometry")
    return retgeom


def _buffer(
    geometry,
    distancem: float,
    geoid: Geod,
    flat: bool,
    feedback,
    precision,
    debug: bool = False,
):
    if isinstance(geometry, list):
        if len(geometry) > 1:
            buffered_coll = []
            for geom in geometry:
                buffered = _buffer(geom, distancem, geoid, flat, feedback, precision)
                buffered_coll.append(buffered)
            return buffered_coll

        geometry = geometry[0]

    if geometry.isMultipart():
        buffered_coll = []
        for part in geometry.parts():
            buffered = _buffer(part, distancem, geoid, flat, feedback, precision)
            buffered_coll.append(buffered)
        return buffered_coll

    previousVertex = None
    previousAz = None
    buffered = list()
    for ix, vertex in enumerate(geometry.vertices()):
        # feedback.pushInfo(f" vtrx: {ix}")
        if ix == 0:
            previousVertex = vertex
            continue
        newbuff = buff_line(
            previousVertex,
            vertex,
            distancem,
            geoid,
            feedback,
            flatstart=flat and ix == 1,
            flatend=flat,
            precision=precision,
        )

        buffered.append(newbuff)
        previousVertex = vertex

    v0 = geometry.vertexAt(0)
    if geometry.wkbType() == polyt:
        if previousVertex != v0:
            buffered.append(
                buff_line(
                    previousVertex, v0, distancem, geoid, feedback, precision=precision
                )
            )
        buffered = unigeom(buffered)
        if distancem < 0.0:
            buffered = geometry.difference(buffered)
        else:
            buffered = unigeom([buffered, geometry])
    elif ix == 0:  # point
        if debug:
            feedback.pushInfo(str(v0.x()))
        buffered = make_arc(v0, distancem, geoid, feedback, precision=precision)
        buffered = pts2qgeom(buffered, feedback)
    else:
        buffered = unigeom(buffered)
    return buffered


def buff_line(
    p1,
    p2,
    distance,
    geoid: Geod,
    feedback,
    flatstart: bool = False,
    flatend: bool = False,
    precision=1.0,
    debug: bool = 0,
):
    az = geoid.inv(p1.x(), p1.y(), p2.x(), p2.y())[0]  # ,caps=512

    lim = abs((az + 90.0) % 360)

    startarc = make_arc(
        p1,
        distance,
        geoid,
        feedback,
        lim,
        lim + 180.0,
        180.0 if flatstart else precision,
    )
    endarc = make_arc(
        p2, distance, geoid, feedback, lim - 180.0, lim, 180.0 if flatend else precision
    )
    # endarc.reverse()
    if debug:
        feedback.pushInfo(str(lim))
        feedback.pushInfo(str(startarc))
        feedback.pushInfo(str(endarc))
    return pts2qgeom(
        [
            startarc
            + endarc
            + [
                startarc[0],
            ]
        ],
        feedback,
    )


# //https://geographiclib.sourceforge.io/Python/doc/code.html#geographiclib.geodesic.Geodesic.Direct
def make_arc(
    srcPnt,
    distance,
    geoid: Geod,
    feedback,
    start: float = 0.0,
    end: float = 360.0,
    precision: float = 1.0,
):
    arc = []
    steps = round((end - start) / precision)

    for az in range(abs(steps)):
        angle = (start + (az * precision)) % 360
        # feedback.pushInfo(str(angle))
        rlong, rlat, raz = geoid.fwd(
            srcPnt.x(), srcPnt.y(), angle, distance
        )  # ,return_back_azimuth =False)
        arc.append(QgsPointXY(rlong, rlat))
    rlong, rlat, raz = geoid.fwd(
        srcPnt.x(), srcPnt.y(), end % 360.0, distance
    )  # ,return_back_azimuth =False)
    arc.append(QgsPointXY(rlong, rlat))
    return arc


class GeoidBufferAlgorithm(QgsProcessingAlgorithm):
    OUTPUT = "OUTPUT"
    DISTM = "DISTM"
    GEOID = "GEOID"
    INPUT = "INPUT"
    DISSB = "DISSB"
    ENDSTYLE = "ENDSTYLE"
    PRECISION = "PRECISION"
    Capstyle = ["Round", "Flat"]

    def initAlgorithm(self, config=None):
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                self.tr("Input layer"),
                [QgsProcessing.TypeVectorAnyGeometry],
            )
        )

        dp = QgsProcessingParameterNumber(
            self.DISTM,
            self.tr("Distance in Meter"),
            QgsProcessingParameterNumber.Double,
            defaultValue=1000.0,
        )
        dp.setIsDynamic(True)
        dp.setDynamicPropertyDefinition(
            QgsPropertyDefinition(
                "Distance", self.tr("Distance in meters"), QgsPropertyDefinition.Integer
            )
        )
        dp.setDynamicLayerParameterName("INPUT")
        self.addParameter(dp)

        self.addParameter(
            QgsProcessingParameterEnum(
                self.ENDSTYLE,
                self.tr("End style"),
                self.Capstyle,
                allowMultiple=False,
                defaultValue=0,
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.DISSB, self.tr("Dissolve features"), False
            )
        )

        self.addParameter(
            QgsProcessingParameterCrs(
                self.GEOID,
                self.tr("Geoid to use"),
                defaultValue=None,
                optional=True,
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(self.OUTPUT, self.tr("Output layer"))
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.PRECISION,
                self.tr("Precision (in degree)"),
                QgsProcessingParameterNumber.Double,
                defaultValue=1.0,
            )
        )

    def prepareAlgorithm(self, parameters, context, feedback):
        self.distanceV = self.parameterAsDouble(parameters, "DISTM", context)

        self.dynamicDist = QgsProcessingParameters.isDynamic(parameters, "DISTM")
        if self.dynamicDist:
            self.distanceExp = parameters["DISTM"]
        else:
            self.distanceExp = ""
        return True

    def processAlgorithm(self, parameters, context, feedback):
        source = self.parameterAsSource(parameters, self.INPUT, context)
        (sink, dest_id) = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            source.fields(),
            mpolyt,
            source.sourceCrs(),
        )

        distancem = self.parameterAsDouble(parameters, self.DISTM, context)
        interimCrs = self.parameterAsCrs(parameters, self.GEOID, context)
        capstyle = self.parameterAsEnum(parameters, self.ENDSTYLE, context)
        dissolveB = self.parameterAsBoolean(parameters, self.DISSB, context)
        precision = self.parameterAsDouble(parameters, self.PRECISION, context)

        fc = source.featureCount() if source.featureCount() else 0
        total = 100.0 / fc
        features = source.getFeatures()

        currentCrs = source.sourceCrs()

        if interimCrs.isValid():
            Ellipsoid = interimCrs
        else:
            Ellipsoid = currentCrs
        expcont = self.createExpressionContext(parameters, context, source)

        for current, feature in enumerate(features):
            if feedback.isCanceled():
                break
            if not (feature.hasGeometry()):
                continue
            oldgeometry = feature.geometry()
            distancem = self.distanceV
            if self.dynamicDist:
                expcont.setFeature(feature)
                distancem, expeval = self.distanceExp.valueAsDouble(expcont, distancem)
                if not (expeval):
                    feedback.pushInfo("error evaluating expression")
                # feedback.pushInfo(str(distancem))

            buffered = buffer(
                oldgeometry,
                distancem,
                currentCrs,
                Ellipsoid,
                feedback,
                dissolveB,
                capstyle == 1,
                precision,
            )

            feature.setGeometry(buffered)

            sink.addFeature(feature, QgsFeatureSink.FastInsert)

            feedback.setProgress(int(current * total))

        return {self.OUTPUT: dest_id}

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return "Geoid buffer"

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
        return "Vector Geometry"

    def tr(self, string):
        return QCoreApplication.translate("Processing", string)

    def createInstance(self):
        return GeoidBufferAlgorithm()

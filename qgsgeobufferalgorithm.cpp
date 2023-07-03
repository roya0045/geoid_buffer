/***************************************************************************
                         qgsalgorithmbuffer.cpp
                         ---------------------
    begin                : April 2017
    copyright            : (C) 2017 by Nyall Dawson
    email                : nyall dot dawson at gmail dot com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "qgsalgorithmbuffer.h"
#include "qgsvectorlayer.h"


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



QgsGeometry ::pts2qgeom(QVector< QVector< QgsPointXY > > pointLists )
{
    QgsMultiPolygon multipoly = QgsMultiPolygon()
    QVector<QVector< QgsPointXY > >::const_iterator pointList = pointLists.constBegin();
    for ( ; pointList != pointLists.constEnd(); ++pointList)
    {
        multipoly.addGeometry( QgsPolygon( QgsLineString( *pointList ) ) )
    }
    return QgsGeometry.fromWkt( multipoly.asWkt() );
};

QgsGeometry ::pts2qgeom(QVector< QgsPointXY > ptslist )
{   QgsMultiPolygon multipoly = QgsMultiPolygon();
    multipoly.addGeometry(QgsPolygon( QgsLineString(ptslist) ) );
    return QgsGeometry.fromWkt( multipoly.asWkt() );
};

QgsGeometry ::buffer(
    QgsGeometry geometry,
    double distancem,
    QgsCoordinateReferenceSystem srcCrs,
    QgsCoordinateReferenceSystem destCrs,
    QgsProcessingFeedback *feedback,
    Bool dissolve=True,
    Bool flatcap=False,
    double precision=1.0,
)
{
    if distancem == 0.0:
        return geometry;
    
    QgsCoordinateReferenceSystem geoidCrs = destCrs.toGeographicCrs();
    if not (geoid.isValid):
        raise QgsProcessingException("invalid crs")

    QgsCoordinateTransform geoidTransformer = QgsCoordinateTransform(srcCrs, geoidCrs, QgsProject.instance());
    QgsCoordinateTransform reverseTransformer = QgsCoordinateTransform(geoidCrs, srcCrs, QgsProject.instance());

    result = geometry.transform(geoidTransformer)
    if !(result || geometry.isGeosValid()):
        raise QgsProcessingException("Failed to transform")

    geod_geodesic geoid = QgsEllipsoidUtils::getGeoid( geoidCrs );

    QVector<QgsGeometry> results = _buffer(
        geometry.asGeometryCollection(), distancem, geoid, flatcap, feedback, precision
    )
    if ( results.size() == 1 )
    {
        QgsGeometry resultsGeometry = results.at(0);
        bufferedGeometry.transform(reverseTransformer);
        if not (bufferedGeometry.isGeosValid()):
            return( bufferedGeometry.makeValid());
        return bufferedGeometry;
    }
    if dissolve:
        buffered = QgsGeometry.unaryUnion(results)
    QVector<QgsGeometry> reprojected;

    QVector<QgsGeometry>::const_iterator buff = results.constBegin();
  for ( ; buff != results.constEnd(); ++buff)
  {
        QgsGeometry 
        buff.transform(revtrsctx)
        if buff.isMultipart():
            ok = buff.convertToSingleType()

        reprojected << buff;
  }

    retgeom = QgsGeometry.unaryUnion(geoms)

    if not (retgeom.isGeosValid()):

        retgeom = retgeom.makeValid()
        # raise QgsProcessingException("could not transform back resulting geometry")
    return retgeom;
}

QVector<QgsGeometry> ::_buffer(
    QVector <QgsGeometry> geometries,
    double distancem,
    geod_geodesic Geod,
    flat: bool,
    QgsProcessingFeedback *feedback,
    double precision,
    debug: bool = False,
){

    QVector<QgsGeometry> results;
    QVector<QgsGeometry>::const_iterator geometry = geometries.constBegin();
    for ( ; geomery != geometries.constEnd(); ++geometry)
    {
        results << _buffer(*geometry, distancem, geoid, flat, feedback, precision);
    }
    return results;

};

QVector<QgsGeometry> ::_buffer(
    QgsGeometry geometry,
    double distancem,
    geod_geodesic Geod,
    flat: bool,
    QgsProcessingFeedback *feedback,
    double precision,
    debug: bool = False,
){

    if (geometry.isMultipart() )
{        QVector<QgsGeometry> buffered;
  QgsGeometryPartIterator parts = geometry->parts();
  while ( parts.hasNext() )
  {
    QgsGeometry *geompart = parts.next();
    buffered << _buffer( *geompart, distancem, geoid, flat, feedback, precision);
  }
  return buffered;
  }


    QgsPoint previousVertex; = None;
    previousAz = None;
    QVector<QgsGeometry> results;
    QgsGeometry buffered;
    int ix;

      QgsVertexIterator vertexIterator = geometry.vertices();
      while ( vertexIterator.hasNext() )
      {
        const QgsPoint &pt = vertexIterator.next();
        if (ix == 0)
        {
            previousVertex = vertex;buffered
            continue;
        }

        results << buff_line(
            previousVertex,
            vertex,
            distancem,
            geoid,
            precision=precision,
            flatstart=flat && ix == 1,
            flatend=flat
        );
        previousVertex = point;
        ix ++;
}

    QgsPoint v0 = geometry.vertexAt(0);
    if geometry.wkbType() == polyt:
        if !(previousVertex == v0)
            results << buff_line(previousVertex, v0, distancem, geoid, precision);
        buffered = QgsGeometry.unaryUnion(results);
        if distancem < 0.0:
            buffered = geometry.difference(buffered);
        else:
            buffered = QgsGeometry.unaryUnion([buffered, geometry])
    elif ix == 0:  # point
        return (pts2qgeom(make_arc(v0, distancem, geoid, precision=precision) ) );
    else:
        buffered = QgsGeometry.unaryUnion(results);
    return buffered
}


QgsGeometry ::buff_line(
    QgsPointXY p1, QgsPointXY p2, double distance,
    geod_geodesic Geod,double precision=1.0,bool flatstart= False,bool flatend= False
)
{
    az = geoid.inv(p1.x(), p1.y(), p2.x(), p2.y())[0]  # ,caps=512
    double lim = abs((az + 90.0) % 360)

    QVector< QgsPointXY >contour = make_arc(p1,distance,geoid, lim,lim + 180.0,( flatstart ) ? 180.0 : precision);
    contour << make_arc( p2, distance, geoid, lim - 180.0, lim, ( flatend ) ? 180.0 : precision);
    contour << contour.first();

    return pts2qgeom( contour )
}

# //https://geographiclib.sourceforge.io/Python/doc/code.html#geographiclib.geodesic.Geodesic.Direct
QVector< QgsPointXY > ::make_arc(
    srcPnt,
    distance,
    geod_geodesic Geod,
    double start = 0.0,
    double  end = 360.0,
    double precision = 1.0,
)
{
    double rlong, rlat, raz;
    double angle;
    QVector< QgsPointXY > points;
    int steps = std:abs(std::round((end - start) / precision));
    int step = 0;
    while ( step < steps ){
            angle = (start + (step * precision)) % 360
            rlong, rlat, raz = geoid.fwd( srcPnt.x(), srcPnt.y(), angle, distance)
            points << QgsPointXY(rlong, rlat);
            rlong, rlat, raz = geoid.fwd(srcPnt.x(), srcPnt.y(), end % 360.0, distance )
            points << QgsPointXY(rlong, rlat);
            step ++;
    }
    return points;
}


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

    def processAlgorithm(self, parameters, context, QgsProcessingFeedback *feedback):
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


///@cond PRIVATE

QString QgsBufferAlgorithm::name() const
{
  return QStringLiteral( "geoidbuffer" );
}

QString QgsBufferAlgorithm::displayName() const
{
  return QObject::tr( "Geoid Buffer" );
}

QStringList QgsBufferAlgorithm::tags() const
{
  return QObject::tr( "buffer,grow,fixed,variable,distance" ).split( ',' );
}

QString QgsBufferAlgorithm::group() const
{
  return QObject::tr( "Vector geometry" );
}

QString QgsBufferAlgorithm::groupId() const
{
  return QStringLiteral( "vectorgeometry" );
}

void QgsBufferAlgorithm::initAlgorithm( const QVariantMap & )
{
  addParameter( new QgsProcessingParameterFeatureSource( QStringLiteral( "INPUT" ), QObject::tr( "Input layer" ) ) );

  auto bufferParam = std::make_unique < QgsProcessingParameterDistance >( QStringLiteral( "DISTANCE" ), QObject::tr( "Distance" ), 10, QStringLiteral( "INPUT" ) );
  bufferParam->setIsDynamic( true );
  bufferParam->setDynamicPropertyDefinition( QgsPropertyDefinition( QStringLiteral( "Distance" ), QObject::tr( "Buffer distance" ), QgsPropertyDefinition::Double ) );
  bufferParam->setDynamicLayerParameterName( QStringLiteral( "INPUT" ) );
  addParameter( bufferParam.release() );
  auto segmentParam = std::make_unique < QgsProcessingParameterNumber >( QStringLiteral( "SEGMENTS" ), QObject::tr( "Segments" ), QgsProcessingParameterNumber::Integer, 5, false, 1 );
  segmentParam->setHelp( QObject::tr( "The segments parameter controls the number of line segments to use to approximate a quarter circle when creating rounded offsets." ) );
  addParameter( segmentParam.release() );
  addParameter( new QgsProcessingParameterEnum( QStringLiteral( "END_CAP_STYLE" ), QObject::tr( "End cap style" ), QStringList() << QObject::tr( "Round" ) << QObject::tr( "Flat" ) << QObject::tr( "Square" ), false, 0 ) );
  addParameter( new QgsProcessingParameterEnum( QStringLiteral( "JOIN_STYLE" ), QObject::tr( "Join style" ), QStringList() << QObject::tr( "Round" ) << QObject::tr( "Miter" ) << QObject::tr( "Bevel" ), false, 0 ) );
  addParameter( new QgsProcessingParameterNumber( QStringLiteral( "MITER_LIMIT" ), QObject::tr( "Miter limit" ), QgsProcessingParameterNumber::Double, 2, false, 1 ) );

  addParameter( new QgsProcessingParameterBoolean( QStringLiteral( "DISSOLVE" ), QObject::tr( "Dissolve result" ), false ) );

  auto keepDisjointParam = std::make_unique < QgsProcessingParameterBoolean >( QStringLiteral( "SEPARATE_DISJOINT" ), QObject::tr( "Keep disjoint results separate" ), false );
  keepDisjointParam->setFlags( keepDisjointParam->flags() | QgsProcessingParameterDefinition::FlagAdvanced );
  keepDisjointParam->setHelp( QObject::tr( "If checked, then any disjoint parts in the buffer results will be output as separate single-part features." ) );
  addParameter( keepDisjointParam.release() );

  addParameter( new QgsProcessingParameterFeatureSink( QStringLiteral( "OUTPUT" ), QObject::tr( "Buffered" ), QgsProcessing::TypeVectorPolygon, QVariant(), false, true, true ) );
}

QString QgsBufferAlgorithm::shortHelpString() const
{
  return QObject::tr( "This algorithm computes a buffer area for all the features in an input layer, using a fixed or dynamic distance.\n\n"
                      "The segments parameter controls the number of line segments to use to approximate a quarter circle when creating rounded offsets.\n\n"
                      "The end cap style parameter controls how line endings are handled in the buffer.\n\n"
                      "The join style parameter specifies whether round, miter or beveled joins should be used when offsetting corners in a line.\n\n"
                      "The miter limit parameter is only applicable for miter join styles, and controls the maximum distance from the offset curve to use when creating a mitered join." );
}

QgsBufferAlgorithm *QgsBufferAlgorithm::createInstance() const
{
  return new QgsBufferAlgorithm();
}

QVariantMap QgsBufferAlgorithm::processAlgorithm( const QVariantMap &parameters, QgsProcessingContext &context, QgsProcessingFeedback *feedback )
{
  std::unique_ptr< QgsProcessingFeatureSource > source( parameterAsSource( parameters, QStringLiteral( "INPUT" ), context ) );
  if ( !source )
    throw QgsProcessingException( invalidSourceError( parameters, QStringLiteral( "INPUT" ) ) );

  QString dest;
  std::unique_ptr< QgsFeatureSink > sink( parameterAsSink( parameters, QStringLiteral( "OUTPUT" ), context, dest, source->fields(), Qgis::WkbType::MultiPolygon, source->sourceCrs() ) );
  if ( !sink )
    throw QgsProcessingException( invalidSinkError( parameters, QStringLiteral( "OUTPUT" ) ) );

  // fixed parameters
  const bool dissolve = parameterAsBoolean( parameters, QStringLiteral( "DISSOLVE" ), context );
  const bool keepDisjointSeparate = parameterAsBoolean( parameters, QStringLiteral( "SEPARATE_DISJOINT" ), context );
  const int segments = parameterAsInt( parameters, QStringLiteral( "SEGMENTS" ), context );
  const Qgis::EndCapStyle endCapStyle = static_cast< Qgis::EndCapStyle >( 1 + parameterAsInt( parameters, QStringLiteral( "END_CAP_STYLE" ), context ) );
  const Qgis::JoinStyle joinStyle = static_cast< Qgis::JoinStyle>( 1 + parameterAsInt( parameters, QStringLiteral( "JOIN_STYLE" ), context ) );
  const double miterLimit = parameterAsDouble( parameters, QStringLiteral( "MITER_LIMIT" ), context );
  const double bufferDistance = parameterAsDouble( parameters, QStringLiteral( "DISTANCE" ), context );
  const bool dynamicBuffer = QgsProcessingParameters::isDynamic( parameters, QStringLiteral( "DISTANCE" ) );
  QgsExpressionContext expressionContext = createExpressionContext( parameters, context, source.get() );
  QgsProperty bufferProperty;
  if ( dynamicBuffer )
  {
    bufferProperty = parameters.value( QStringLiteral( "DISTANCE" ) ).value< QgsProperty >();
  }

  const long count = source->featureCount();

  QgsFeature f;
  // buffer doesn't care about invalid features, and buffering can be used to repair geometries
  QgsFeatureIterator it = source->getFeatures( QgsFeatureRequest(), QgsProcessingFeatureSource::FlagSkipGeometryValidityChecks );

  const double step = count > 0 ? 100.0 / count : 1;
  int current = 0;

  QVector< QgsGeometry > bufferedGeometriesForDissolve;
  QgsAttributes dissolveAttrs;

  while ( it.nextFeature( f ) )
  {
    if ( feedback->isCanceled() )
    {
      break;
    }
    if ( dissolveAttrs.isEmpty() )
      dissolveAttrs = f.attributes();

    QgsFeature out = f;
    if ( out.hasGeometry() )
    {
      double distance =  bufferDistance;
      if ( dynamicBuffer )
      {
        expressionContext.setFeature( f );
        distance = bufferProperty.valueAsDouble( expressionContext, bufferDistance );
      }

      QgsGeometry outputGeometry = f.geometry().buffer( distance, segments, endCapStyle, joinStyle, miterLimit );
      if ( outputGeometry.isNull() )
      {
        QgsMessageLog::logMessage( QObject::tr( "Error calculating buffer for feature %1" ).arg( f.id() ), QObject::tr( "Processing" ), Qgis::MessageLevel::Warning );
      }
      if ( dissolve )
      {
        bufferedGeometriesForDissolve << outputGeometry;
      }
      else
      {
        outputGeometry.convertToMultiType();

        if ( !keepDisjointSeparate )
        {
          out.setGeometry( outputGeometry );

          if ( !sink->addFeature( out, QgsFeatureSink::FastInsert ) )
            throw QgsProcessingException( writeFeatureError( sink.get(), parameters, QStringLiteral( "OUTPUT" ) ) );
        }
        else
        {
          for ( auto partIt = outputGeometry.const_parts_begin(); partIt != outputGeometry.const_parts_end(); ++partIt )
          {
            if ( const QgsAbstractGeometry *part = *partIt )
            {
              out.setGeometry( QgsGeometry( part->clone() ) );
              if ( !sink->addFeature( out, QgsFeatureSink::FastInsert ) )
                throw QgsProcessingException( writeFeatureError( sink.get(), parameters, QStringLiteral( "OUTPUT" ) ) );
            }
          }
        }
      }
    }
    else if ( !dissolve )
    {
      if ( !sink->addFeature( out, QgsFeatureSink::FastInsert ) )
        throw QgsProcessingException( writeFeatureError( sink.get(), parameters, QStringLiteral( "OUTPUT" ) ) );
    }

    feedback->setProgress( current * step );
    current++;
  }

  if ( dissolve && !bufferedGeometriesForDissolve.isEmpty() )
  {
    QgsGeometry finalGeometry = QgsGeometry::unaryUnion( bufferedGeometriesForDissolve );
    finalGeometry.convertToMultiType();
    QgsFeature f;
    f.setAttributes( dissolveAttrs );

    if ( !keepDisjointSeparate )
    {
      f.setGeometry( finalGeometry );
      if ( !sink->addFeature( f, QgsFeatureSink::FastInsert ) )
        throw QgsProcessingException( writeFeatureError( sink.get(), parameters, QStringLiteral( "OUTPUT" ) ) );
    }
    else
    {
      for ( auto partIt = finalGeometry.const_parts_begin(); partIt != finalGeometry.const_parts_end(); ++partIt )
      {
        if ( const QgsAbstractGeometry *part = *partIt )
        {
          f.setGeometry( QgsGeometry( part->clone() ) );
          if ( !sink->addFeature( f, QgsFeatureSink::FastInsert ) )
            throw QgsProcessingException( writeFeatureError( sink.get(), parameters, QStringLiteral( "OUTPUT" ) ) );
        }
      }
    }
  }

  QVariantMap outputs;
  outputs.insert( QStringLiteral( "OUTPUT" ), dest );
  return outputs;
}

QgsProcessingAlgorithm::Flags QgsBufferAlgorithm::flags() const
{
  Flags f = QgsProcessingAlgorithm::flags();
  f |= QgsProcessingAlgorithm::FlagSupportsInPlaceEdits;
  return f;
}

QgsProcessingAlgorithm::VectorProperties QgsBufferAlgorithm::sinkProperties( const QString &sink, const QVariantMap &parameters, QgsProcessingContext &context, const QMap<QString, QgsProcessingAlgorithm::VectorProperties> &sourceProperties ) const
{
  const bool keepDisjointSeparate = parameterAsBoolean( parameters, QStringLiteral( "SEPARATE_DISJOINT" ), context );

  QgsProcessingAlgorithm::VectorProperties result;
  if ( sink == QLatin1String( "OUTPUT" ) )
  {
    if ( sourceProperties.value( QStringLiteral( "INPUT" ) ).availability == QgsProcessingAlgorithm::Available )
    {
      const VectorProperties inputProps = sourceProperties.value( QStringLiteral( "INPUT" ) );
      result.fields = inputProps.fields;
      result.crs = inputProps.crs;
      result.wkbType = keepDisjointSeparate ? Qgis::WkbType::Polygon : Qgis::WkbType::MultiPolygon;
      result.availability = Available;
      return result;
    }
    else
    {
      std::unique_ptr< QgsProcessingFeatureSource > source( parameterAsSource( parameters, QStringLiteral( "INPUT" ), context ) );
      if ( source )
      {
        result.fields = source->fields();
        result.crs = source->sourceCrs();
        result.wkbType = keepDisjointSeparate ? Qgis::WkbType::Polygon : Qgis::WkbType::MultiPolygon;
        result.availability = Available;
        return result;
      }
    }
  }
  return result;
}

bool QgsBufferAlgorithm::supportInPlaceEdit( const QgsMapLayer *layer ) const
{
  const QgsVectorLayer *vlayer = qobject_cast< const QgsVectorLayer * >( layer );
  if ( !vlayer )
    return false;
  //Only Polygons
  return vlayer->wkbType() == Qgis::WkbType::Polygon || vlayer->wkbType() == Qgis::WkbType::MultiPolygon;
}

///@endcond
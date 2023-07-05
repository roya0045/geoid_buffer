/***************************************************************************
                         qgsalgorithmbuffer.cpp
                         ---------------------
    begin                : AJuly 2023
    copyright            : (C) 2023 by Alexis RL
    email                : roya0045 at github dot com
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


///@cond PRIVATE

QgsGeometry QgsBufferAlgorithm::pts2qgeom(QVector<QVector<QgsPointXY>> pointLists)
{
  QgsMultiPolygon multipoly = QgsMultiPolygon()
      QVector<QVector<QgsPointXY>>::const_iterator pointList = pointLists.constBegin();
  for (; pointList != pointLists.constEnd(); ++pointList)
  {
    multipoly.addGeometry(QgsPolygon(QgsLineString(*pointList)))
  }
  return QgsGeometry.fromWkt(multipoly.asWkt());
};

QgsGeometry QgsBufferAlgorithm::pts2qgeom(QVector<QgsPointXY> ptslist)
{
  QgsMultiPolygon multipoly = QgsMultiPolygon();
  multipoly.addGeometry(QgsPolygon(QgsLineString(ptslist)));
  return QgsGeometry.fromWkt(multipoly.asWkt());
};

QgsGeometry QgsBufferAlgorithm::buffer(
    QgsGeometry geometry,
    double distancem,
    QgsCoordinateReferenceSystem srcCrs,
    QgsCoordinateReferenceSystem destCrs,
    QgsProcessingFeedback *feedback,
    Bool dissolve = True,
    Bool flatcap = False,
    double precision = 1.0, )
{
  if distancem
    == 0.0 : return geometry;

  QgsCoordinateReferenceSystem geoidCrs = destCrs.toGeographicCrs();
  if not(geoid.isValid):
        raise QgsProcessingException("invalid crs")

  QgsCoordinateTransform geoidTransformer = QgsCoordinateTransform(srcCrs, geoidCrs, QgsProject.instance());
  QgsCoordinateTransform reverseTransformer = QgsCoordinateTransform(geoidCrs, srcCrs, QgsProject.instance());

  result = geometry.transform(geoidTransformer);
  if !(result || geometry.isGeosValid()) : 
   raise QgsProcessingException("Failed to transform");

  geod_geodesic geoid = QgsEllipsoidUtils::getGeoid(geoidCrs);

  QVector<QgsGeometry> results = _buffer(geometry.asGeometryCollection(), distancem, geoid, flatcap, feedback, precision);
  if (results.size() == 1)
  {
    QgsGeometry resultsGeometry = results.at(0);
    bufferedGeometry.transform(reverseTransformer);
    if not(bufferedGeometry.isGeosValid()):
            return( bufferedGeometry.makeValid());
    return bufferedGeometry;
  }
  if (dissolve)
    buffered = QgsGeometry.unaryUnion(results);
  QVector<QgsGeometry> reprojected;

  QVector<QgsGeometry>::const_iterator buff = results.constBegin();
  for (; buff != results.constEnd(); ++buff)
  {
    QgsGeometry earthBuffer = *buff;
    earthBuffer.transform(reverseTransformer);
    if (buff.isMultipart())
      earthBuffer.convertToSingleType();

    reprojected << earthBuffer;
  }

  retgeom = QgsGeometry.unaryUnion(geoms);

  if not(retgeom.isGeosValid())
    retgeom = retgeom.makeValid();

  return retgeom;
}

QVector<QgsGeometry> QgsBufferAlgorithm::_buffer(
    QVector<QgsGeometry> geometries,
    double distancem,
    geod_geodesic Geod,
    bool flat,
    QgsProcessingFeedback *feedback,
    double precision)
{

  QVector<QgsGeometry> results;
  QVector<QgsGeometry>::const_iterator geometry = geometries.constBegin();
  for (; geomery != geometries.constEnd(); ++geometry)
  {
    results << _buffer(*geometry, distancem, geoid, flat, feedback, precision);
  }
  return results;
};

QVector<QgsGeometry> QgsBufferAlgorithm::_buffer(
    QgsGeometry geometry,
    double distancem,
    geod_geodesic Geod,
    bool flat,
    QgsProcessingFeedback *feedback,
    double precision)
{

  if (geometry.isMultipart())
  {
    QVector<QgsGeometry> buffered;
    QgsGeometryPartIterator parts = geometry->parts();
    while (parts.hasNext())
    {
      QgsGeometry *geompart = parts.next();
      buffered << _buffer(*geompart, distancem, geoid, flat, feedback, precision);
    }
    return buffered;
  }

  QgsPoint previousVertex;
  QVector<QgsGeometry> results;
  QgsGeometry buffered;
  int ix;

  QgsVertexIterator vertexIterator = geometry.vertices();
  while (vertexIterator.hasNext())
  {
    const QgsPoint &pt = vertexIterator.next();
    if (ix == 0)
    {
      previousVertex = vertex;
      buffered continue;
    }

    results << buff_line(previousVertex, vertex, distancem, geoid, precision = precision, flatstart = flat && ix == 1, flat);
    previousVertex = point;
    ix++;
  }

  QgsPoint v0 = geometry.vertexAt(0);
  if (geometry.wkbType() == Qgis.WkbType.PolygonGeomettry)
  {
    if !(previousVertex == v0)
      results << buff_line(previousVertex, v0, distancem, geoid, precision);
    buffered = QgsGeometry.unaryUnion(results);
    if (distancem < 0.0)
      buffered = geometry.difference(buffered);
    else:
      buffered = QgsGeometry.unaryUnion([buffered, geometry]);
  }
  else if (ix == 0)
    return pts2qgeom(make_arc(v0, distancem, geoid, precision = precision));
  else
    buffered = QgsGeometry.unaryUnion(results);
  return buffered;
}

QgsGeometry QgsBufferAlgorithm::buff_line(
    QgsPointXY p1, QgsPointXY p2, double distance,
    geod_geodesic Geod, double precision = 1.0, bool flatstart = False, bool flatend = False)
{
  double az;
  QgsEllipsoidUtils::geoidInverseTransform(Geod, p1, p2, &az);
  double lim = std::abs((az + 90.0) % 360.0);

  QVector<QgsPointXY> contour = make_arc(p1, distance, Geod, lim, lim + 180.0, (flatstart) ? 180.0 : precision);
  contour << make_arc(p2, distance, Geod, lim - 180.0, lim, (flatend) ? 180.0 : precision);
  contour << contour.first();

  return pts2qgeom(contour)
}

#//https://geographiclib.sourceforge.io/Python/doc/code.html#geographiclib.geodesic.Geodesic.Direct
QVector<QgsPointXY> QgsBufferAlgorithm::make_arc(
    srcPnt,
    distance,
    geod_geodesic Geod,
    double start = 0.0,
    double end = 360.0,
    double precision = 1.0, )
{
  double rlong, rlat, raz;
  double angle;
  QVector<QgsPointXY> points;
  int steps = std : abs(std::round((end - start) / precision));
  int step = 0;
  while (step < steps)
  {
    angle = (start + (step * precision)) % 360;
    rlong, rlat, raz = QgsEllipsoidUtils::geoidDirectTransform(Geod, srcPnt, angle, distance);
    points << QgsPointXY(rlong, rlat);
    step++;
  }
  rlong, rlat, raz = QgsEllipsoidUtils::geoidDirectTransform(Geod, srcPnt, end % 360.0, distance);
  points << QgsPointXY(rlong, rlat);
  return points;
}

QString QgsBufferAlgorithm::name() const
{
  return QStringLiteral("geoidbuffer");
}

QString QgsBufferAlgorithm::displayName() const
{
  return QObject::tr("Geoid Buffer");
}

QStringList QgsBufferAlgorithm::tags() const
{
  return QObject::tr("buffer,grow,fixed,variable,distance").split(',');
}

QString QgsBufferAlgorithm::group() const
{
  return QObject::tr("Vector geometry");
}

QString QgsBufferAlgorithm::groupId() const
{
  return QStringLiteral("vectorgeometry");
}

void QgsBufferAlgorithm::initAlgorithm(const QVariantMap &)
{
  addParameter(new QgsProcessingParameterFeatureSource(QStringLiteral("INPUT"), QObject::tr("Input layer")));

  auto bufferParam = std::make_unique<QgsProcessingParameterNumber>(QStringLiteral("DISTANCE"), QObject::tr("Distance in meter"),QgsProcessingParameterNumber::Double, 1000);
  bufferParam->setIsDynamic(true);
  bufferParam->setDynamicPropertyDefinition(QgsPropertyDefinition(QStringLiteral("Distance"), QObject::tr("Buffer distance"), QgsPropertyDefinition::Double));
  bufferParam->setDynamicLayerParameterName(QStringLiteral("INPUT"));
  addParameter(bufferParam.release());
  addParameter(new QgsProcessingParameterEnum(QStringLiteral("END_CAP_STYLE"), QObject::tr("End cap style"), QStringList() << QObject::tr("Round") << QObject::tr("Flat") << QObject::tr("Square"), false, 0));
  addParameter(new QgsProcessingParameterNumber(QStringLiteral("PRECISION"), QObject::tr( "Precision (in degree)" ), QgsProcessingParameterNumber::Double, 1));
  addParameter(new QgsProcessingParameterBoolean(QStringLiteral("DISSOLVE"), QObject::tr("Dissolve result"), false));

  auto keepDisjointParam = std::make_unique<QgsProcessingParameterBoolean>(QStringLiteral("SEPARATE_DISJOINT"), QObject::tr("Keep disjoint results separate"), false);
  keepDisjointParam->setFlags(keepDisjointParam->flags() | QgsProcessingParameterDefinition::FlagAdvanced);
  keepDisjointParam->setHelp(QObject::tr("If checked, then any disjoint parts in the buffer results will be output as separate single-part features."));
  addParameter(keepDisjointParam.release());

  addParameter(new QgsProcessingParameterFeatureSink(QStringLiteral("OUTPUT"), QObject::tr("Buffered"), QgsProcessing::TypeVectorPolygon, QVariant(), false, true, true));
}

QString QgsBufferAlgorithm::shortHelpString() const
{
  return QObject::tr("This algorithm computes a buffer area for all the features in an input layer, using a fixed or dynamic distance.\n\n"
                     "The segments parameter controls the number of line segments to use to approximate a quarter circle when creating rounded offsets.\n\n"
                     "The end cap style parameter controls how line endings are handled in the buffer.\n\n"
                     "The join style parameter specifies whether round, miter or beveled joins should be used when offsetting corners in a line.\n\n"
                     "The miter limit parameter is only applicable for miter join styles, and controls the maximum distance from the offset curve to use when creating a mitered join.");
}

QgsBufferAlgorithm *QgsBufferAlgorithm::createInstance() const
{
  return new QgsBufferAlgorithm();
}

QVariantMap QgsBufferAlgorithm::processAlgorithm(const QVariantMap &parameters, QgsProcessingContext &context, QgsProcessingFeedback *feedback)
{
  std::unique_ptr<QgsProcessingFeatureSource> source(parameterAsSource(parameters, QStringLiteral("INPUT"), context));
  if (!source)
    throw QgsProcessingException(invalidSourceError(parameters, QStringLiteral("INPUT")));

  QString dest;
  std::unique_ptr<QgsFeatureSink> sink(parameterAsSink(parameters, QStringLiteral("OUTPUT"), context, dest, source->fields(), Qgis::WkbType::MultiPolygon, source->sourceCrs()));
  if (!sink)
    throw QgsProcessingException(invalidSinkError(parameters, QStringLiteral("OUTPUT")));

  // fixed parameters
  const bool dissolve = parameterAsBoolean(parameters, QStringLiteral("DISSOLVE"), context);
  const bool keepDisjointSeparate = parameterAsBoolean(parameters, QStringLiteral("SEPARATE_DISJOINT"), context);
  const Qgis::EndCapStyle endCapStyle = static_cast<Qgis::EndCapStyle>(1 + parameterAsInt(parameters, QStringLiteral("END_CAP_STYLE"), context));
  const double precision = parameterAsDouble(parameters, QStringLiteral("PRECISION"), context);
  const double bufferDistance = parameterAsDouble(parameters, QStringLiteral("DISTANCE"), context);
  const bool dynamicBuffer = QgsProcessingParameters::isDynamic(parameters, QStringLiteral("DISTANCE"));
  QgsExpressionContext expressionContext = createExpressionContext(parameters, context, source.get());
  QgsProperty bufferProperty;
  if (dynamicBuffer)
  {
    bufferProperty = parameters.value(QStringLiteral("DISTANCE")).value<QgsProperty>();
  }

  const long count = source->featureCount();

  QgsFeature f;
  // buffer doesn't care about invalid features, and buffering can be used to repair geometries
  QgsFeatureIterator it = source->getFeatures(QgsFeatureRequest(), QgsProcessingFeatureSource::FlagSkipGeometryValidityChecks);

  const double step = count > 0 ? 100.0 / count : 1;
  int current = 0;
  QgsCoordinateReferenceSystem currentCrs = source.sourceCrs();
  QVector<QgsGeometry> bufferedGeometriesForDissolve;
  QgsAttributes dissolveAttrs;

  while (it.nextFeature(f))
  {
    if (feedback->isCanceled())
    {
      break;
    }
    if (dissolveAttrs.isEmpty())
      dissolveAttrs = f.attributes();

    QgsFeature outputGeometry;
    if (f.hasGeometry())
    {
      double distance = bufferDistance;
      if (dynamicBuffer)
      {
        expressionContext.setFeature(f);
        distance = bufferProperty.valueAsDouble(expressionContext, bufferDistance);
      }

      outputGeometry = buffer(
                f,
                distancem,
                currentCrs,
                feedback,
                dissolve,
                capstyle == 1,
                precision,
            )
      if (outputGeometry.isNull())
      {
        QgsMessageLog::logMessage(QObject::tr("Error calculating buffer for feature %1").arg(f.id()), QObject::tr("Processing"), Qgis::MessageLevel::Warning);
      }
      if (dissolve)
      {
        bufferedGeometriesForDissolve << outputGeometry;
      }
      else
      {
        outputGeometry.convertToMultiType();

        if (!keepDisjointSeparate)
        {
          out.setGeometry(outputGeometry);

          if (!sink->addFeature(out, QgsFeatureSink::FastInsert))
            throw QgsProcessingException(writeFeatureError(sink.get(), parameters, QStringLiteral("OUTPUT")));
        }
        else
        {
          for (auto partIt = outputGeometry.const_parts_begin(); partIt != outputGeometry.const_parts_end(); ++partIt)
          {
            if (const QgsAbstractGeometry *part = *partIt)
            {
              out.setGeometry(QgsGeometry(part->clone()));
              if (!sink->addFeature(out, QgsFeatureSink::FastInsert))
                throw QgsProcessingException(writeFeatureError(sink.get(), parameters, QStringLiteral("OUTPUT")));
            }
          }
        }
      }
    }
    else if (!dissolve)
    {
      if (!sink->addFeature(out, QgsFeatureSink::FastInsert))
        throw QgsProcessingException(writeFeatureError(sink.get(), parameters, QStringLiteral("OUTPUT")));
    }

    feedback->setProgress(current * step);
    current++;
  }

  if (dissolve && !bufferedGeometriesForDissolve.isEmpty())
  {
    QgsGeometry finalGeometry = QgsGeometry::unaryUnion(bufferedGeometriesForDissolve);
    finalGeometry.convertToMultiType();
    QgsFeature f;
    f.setAttributes(dissolveAttrs);

    if (!keepDisjointSeparate)
    {
      f.setGeometry(finalGeometry);
      if (!sink->addFeature(f, QgsFeatureSink::FastInsert))
        throw QgsProcessingException(writeFeatureError(sink.get(), parameters, QStringLiteral("OUTPUT")));
    }
    else
    {
      for (auto partIt = finalGeometry.const_parts_begin(); partIt != finalGeometry.const_parts_end(); ++partIt)
      {
        if (const QgsAbstractGeometry *part = *partIt)
        {
          f.setGeometry(QgsGeometry(part->clone()));
          if (!sink->addFeature(f, QgsFeatureSink::FastInsert))
            throw QgsProcessingException(writeFeatureError(sink.get(), parameters, QStringLiteral("OUTPUT")));
        }
      }
    }
  }

  QVariantMap outputs;
  outputs.insert(QStringLiteral("OUTPUT"), dest);
  return outputs;
}

QgsProcessingAlgorithm::Flags QgsBufferAlgorithm::flags() const
{
  Flags f = QgsProcessingAlgorithm::flags();
  f |= QgsProcessingAlgorithm::FlagSupportsInPlaceEdits;
  return f;
}

QgsProcessingAlgorithm::VectorProperties QgsBufferAlgorithm::sinkProperties(const QString &sink, const QVariantMap &parameters, QgsProcessingContext &context, const QMap<QString, QgsProcessingAlgorithm::VectorProperties> &sourceProperties) const
{
  const bool keepDisjointSeparate = parameterAsBoolean(parameters, QStringLiteral("SEPARATE_DISJOINT"), context);

  QgsProcessingAlgorithm::VectorProperties result;
  if (sink == QLatin1String("OUTPUT"))
  {
    if (sourceProperties.value(QStringLiteral("INPUT")).availability == QgsProcessingAlgorithm::Available)
    {
      const VectorProperties inputProps = sourceProperties.value(QStringLiteral("INPUT"));
      result.fields = inputProps.fields;
      result.crs = inputProps.crs;
      result.wkbType = keepDisjointSeparate ? Qgis::WkbType::Polygon : Qgis::WkbType::MultiPolygon;
      result.availability = Available;
      return result;
    }
    else
    {
      std::unique_ptr<QgsProcessingFeatureSource> source(parameterAsSource(parameters, QStringLiteral("INPUT"), context));
      if (source)
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

bool QgsBufferAlgorithm::supportInPlaceEdit(const QgsMapLayer *layer) const
{
  const QgsVectorLayer *vlayer = qobject_cast<const QgsVectorLayer *>(layer);
  if (!vlayer)
    return false;
  // Only Polygons
  return vlayer->wkbType() == Qgis::WkbType::Polygon || vlayer->wkbType() == Qgis::WkbType::MultiPolygon;
}

///@endcond
package com.kygis.algorithmofdam;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.geotools.data.DefaultTransaction;
import org.geotools.data.Transaction;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import com.kygis.kuibamodel.po.SectionExactParam;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineSegment;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.LinearRing;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.operation.polygonize.Polygonizer;

final class ShpAnalyseTool {
	/**
	 * 将LineString按length进行分割。
	 * 
	 * @param ls
	 * @param length
	 * @return
	 */
	public static List<LineString> splitLineStringIntoParts(LineString ls, int length) {
		// list for linesegments from input linestring
		List<LineSegment> lineSegmentList = new ArrayList<LineSegment>();
		// create LineSegment objects from input linestring and add them to list
		for (int i = 1; i < ls.getCoordinates().length; i++) {
			lineSegmentList.add(new LineSegment(ls.getCoordinates()[i - 1], ls.getCoordinates()[i]));
		}

		return splitLineStringIntoParts(lineSegmentList, length);
	}

	/**
	 * @param lineSegmentList
	 *            连续的线段。每个线段的结束点是下一个线段的起始点。
	 * @param length
	 * @return
	 */
	public static List<LineString> splitLineStringIntoParts(List<LineSegment> lineSegmentList, int length) {
		List<LineString> resultList = new ArrayList<LineString>();
		LineString currentLineString = null;
		int neededLength = length;
		for (LineSegment s : lineSegmentList) {
			// int distance = (int)JTS.orthodromicDistance(s.p0, s.p1, crs);//经纬度下，两点间距离
			while (s.getLength() > 0) {
				int distance = (int) s.getLength();
				// case: current segment is small enough to be added to the linestring
				if (distance <= neededLength) {
					// create linestring if it does not exist
					if (currentLineString == null) {
						currentLineString = new GeometryFactory()
								.createLineString(new Coordinate[] { new Coordinate(s.p0), new Coordinate(s.p1) });
						// just add the new endpoint otherwise
					} else {
						Coordinate[] coords = new Coordinate[currentLineString.getCoordinates().length + 1];
						// copy old coordinates
						System.arraycopy(currentLineString.getCoordinates(), 0, coords, 0,
								currentLineString.getCoordinates().length);
						// add new coordinate at the end
						coords[coords.length - 1] = new Coordinate(s.p1);
						// create new linestring
						currentLineString = new GeometryFactory().createLineString(coords);
					}
					neededLength -= distance;
					s.setCoordinates(s.p1, s.p1);
					// add linestring to result list if needed length is 0
					if (neededLength == 0) {
						resultList.add(currentLineString);
						currentLineString = null;
						neededLength = length;
					}
					// current segment needs to be cut and added to the linestring
				} else {
					// get coordinate at desired distance (endpoint of linestring)
					Coordinate endPoint = s.pointAlong((double) neededLength / (double) distance);
					// create linestring if it does not exist
					if (currentLineString == null) {
						currentLineString = new GeometryFactory()
								.createLineString(new Coordinate[] { new Coordinate(s.p0), endPoint });
						// just add the new endpoint otherwise
					} else {
						// add new coordinate to linestring
						Coordinate[] coords = new Coordinate[currentLineString.getCoordinates().length + 1];
						// copy old coordinates
						System.arraycopy(currentLineString.getCoordinates(), 0, coords, 0,
								currentLineString.getCoordinates().length);
						// add new coordinate at the end
						coords[coords.length - 1] = endPoint;
						currentLineString = new GeometryFactory().createLineString(coords);
					}
					// add linestring to result list
					resultList.add(currentLineString);
					// reset needed length
					neededLength = length;
					// reset current linestring
					currentLineString = null;
					// adjust segment (calculated endpoint is the new startpoint)
					s.setCoordinates(endPoint, s.p1);
				}
			}
		}
		// add last linestring if there is a rest
		if (neededLength < length) {
			resultList.add(currentLineString);
		}

		// for(int i=0;i<resultList.size();i++)
		// {
		// System.out.println(resultList.get(i).toString());
		// }

		return resultList;
	}

	/**
	 * 将LineString按length进行分割，一直到lastPoint。
	 * 
	 * @param ls
	 * @param length
	 *            平面距离
	 * @param lastPoint
	 *            最后一个点
	 * @return
	 */
	public static List<LineString> splitLineStringIntoParts(LineString ls, int length, double[] lastPoint) {
		// list for linesegments from input linestring
		List<LineSegment> lineSegmentList = new ArrayList<LineSegment>();
		// create LineSegment objects from input linestring and add them to list

		Coordinate lastP = new Coordinate(lastPoint[0], lastPoint[1]);// 出口断面坐标
		Coordinate[] coords = ls.getCoordinates();
		for (int i = 1; i < coords.length; i++) {
			LineSegment tmpSeg = new LineSegment(coords[i - 1], coords[i]);
			double x = tmpSeg.distancePerpendicular(lastP);
			if (x == 0) {
				if (lastP.x != coords[i - 1].x || lastP.y != coords[i - 1].y) {
					lineSegmentList.add(new LineSegment(coords[i - 1], lastP));
				}
				break;
			} else {
				lineSegmentList.add(tmpSeg);
			}
		}
		return splitLineStringIntoParts(lineSegmentList, length);
	}

	public static List<LineString> splitLineStringIntoParts(LineString ls, int length, double[] startPoint,
			double[] lastPoint) {
		float precision = 50f;// TODO 判断起始点和结束点时的精度。因为坐标不一定为严格落在河道上。
		// list for linesegments from input linestring
		List<LineSegment> lineSegmentList = new ArrayList<LineSegment>();
		// create LineSegment objects from input linestring and add them to list
		Coordinate startP = new Coordinate(startPoint[0], startPoint[1]);
		Coordinate lastP = new Coordinate(lastPoint[0], lastPoint[1]);// 出口断面坐标
		Coordinate[] coords = ls.getCoordinates();
		boolean isStart = false;
		for (int i = 1; i < coords.length; i++) {
			LineSegment aimSeg = null;
			LineSegment tmpSeg = new LineSegment(coords[i - 1], coords[i]);
			double x = tmpSeg.distance(startP);
			double y = tmpSeg.distance(lastP);
			if (Double.compare(x, precision) < 0 && !isStart) {
				if (Double.compare(y, precision) < 0) {
					if ((lastP.x != coords[i - 1].x || lastP.y != coords[i - 1].y)
							&& (startP.x != coords[i].x || startP.y != coords[i].y)) {
						lineSegmentList.add(new LineSegment(startP, lastP));
					}
					break;
				} else {
					if (startP.x == coords[i].x && startP.y == coords[i].y) {
						continue;
					} else {
						aimSeg = new LineSegment(startP, coords[i]);
						isStart = true;
					}
				}
			} else if (Double.compare(y, precision) < 0) {
				if (lastP.x != coords[i - 1].x || lastP.y != coords[i - 1].y) {
					lineSegmentList.add(new LineSegment(coords[i - 1], lastP));
				}
				break;
			}

			if (null == aimSeg) {
				aimSeg = tmpSeg;
			}

			if (isStart && null != aimSeg) {
				lineSegmentList.add(aimSeg);
			}
		}
		// for(int i=0;i<lineSegmentList.size();i++)
		// {
		// System.out.println(lineSegmentList.get(i).toString());
		// }
		return splitLineStringIntoParts(lineSegmentList, length);
	}
	
	public static LineString[] splitLineStringIntoParts(LineString ls, double[] startPoint,
			double[] lastPoint) {
		float precision = 50f;// 判断起始点和结束点时的精度。因为坐标不一定为严格落在河道上。
		// list for linesegments from input linestring
		List<LineSegment> lineSegmentList = new ArrayList<LineSegment>();
		// create LineSegment objects from input linestring and add them to list
		Coordinate startP = new Coordinate(startPoint[0], startPoint[1]);
		Coordinate lastP = new Coordinate(lastPoint[0], lastPoint[1]);// 出口断面坐标
		Coordinate[] coords = ls.getCoordinates();
		
		List<Coordinate> pre = new ArrayList<Coordinate>();
		List<Coordinate> mid = new ArrayList<Coordinate>();
		List<Coordinate> aft = new ArrayList<Coordinate>();
		
		boolean isStart = false,isStartNew=false;		
		boolean isEnd=false,isEndNew=false;
		for (int i = 1; i < coords.length; i++) {
			LineSegment aimSeg = null;
			LineSegment tmpSeg = new LineSegment(coords[i - 1], coords[i]);
			double x = tmpSeg.distance(startP);
			double y = tmpSeg.distance(lastP);
			if (Double.compare(x, precision) < 0 && !isStart) {
				if (Double.compare(y, precision) < 0) {
					if ((lastP.x != coords[i - 1].x || lastP.y != coords[i - 1].y)
							&& (startP.x != coords[i].x || startP.y != coords[i].y)) {
						aimSeg=new LineSegment(startP, lastP);
						isEnd=true;
					}
					break;
				} else {
					if (startP.x == coords[i].x && startP.y == coords[i].y) {
						continue;
					} else {				
						Coordinate starPr = tmpSeg.project(startP);
						aimSeg = new LineSegment(starPr, coords[i]);
						isStart = true;
						if(starPr.x != coords[i-1].x || starPr.y != coords[i-1].y)
						{
							isStartNew=true;
						}
					}
				}
			} else if (Double.compare(y, precision) < 0) {
				/*if (lastP.x != coords[i - 1].x || lastP.y != coords[i - 1].y) {
					aimSeg=new LineSegment(coords[i - 1], lastP);
					isEnd=true;
				}*/
				Coordinate lastPr = tmpSeg.project(lastP);
				aimSeg=new LineSegment(lastPr,coords[i] );
				isEnd=true;			
				if(lastPr.x != coords[i].x || lastPr.y != coords[i].y)
				{
					isEndNew=true;
				}
			}

			if (null == aimSeg) {
				aimSeg = tmpSeg;
			}
			
			if(!isStart&&!isEnd)
			{
				pre.add(aimSeg.p0);
			}
			else if (isStart &&!isEnd) 
			{			
					mid.add(aimSeg.p0);
					if(isStartNew)
					{
						pre.add(coords[i - 1]);
						pre.add(aimSeg.p0);
						isStartNew=false;
					}
			}
			else if(isEnd)
			{
				aft.add(aimSeg.p0);
				if(isEndNew)
				{
					mid.add(coords[i-1]);
					mid.add(aimSeg.p0);
					isEndNew=false;
				}
			}
		}
		aft.add(coords[coords.length-1]);
		GeometryFactory geomFactory = new GeometryFactory();
		LineString[] threeLine = new LineString[3];
		threeLine[0]=geomFactory.createLineString(pre.toArray(new Coordinate[0]));
		threeLine[1]=geomFactory.createLineString(mid.toArray(new Coordinate[0]));
		threeLine[2]=geomFactory.createLineString(aft.toArray(new Coordinate[0]));
		
		return threeLine;
	}

	public static List<LineString> splitLineStringIntoParts(LineString ls, SectionExactParam sectionExactParam) {
		List<LineString> geoObjes;
		double[] endPoint = sectionExactParam.getEndPoint();
		double[] startPoint = sectionExactParam.getStartPoint();
		if (null == endPoint) {
			geoObjes = ShpAnalyseTool.splitLineStringIntoParts(ls, sectionExactParam.getStep());
		}
		geoObjes = ShpAnalyseTool.splitLineStringIntoParts(ls, sectionExactParam.getStep(), startPoint, endPoint);// 指定了出口断面位置。
		return geoObjes;
	}

	/**
	 * 线段末端点的正交线段。以线段末端为交点，做正交线段。线段长度与原线段长度相同。
	 * 
	 * @param c0
	 *            线段起点
	 * @param c1
	 *            线段终点
	 * @return
	 */
	public static LineSegment orthogonalLineOfLinesegment(Coordinate c0, Coordinate c1) {
		double x, y;

		double dx = c1.x - c0.x;
		double dy = c1.y - c0.y;

		x = c1.x - dy;
		y = c1.y + dx;

		return new LineSegment(new Coordinate(c1.x, c1.y), new Coordinate(x, y));
	}
	
	public static LineSegment orthogonalLineOfLinesegmentA(Coordinate c0, Coordinate c1) {
		double x, y;

		double dx = c0.x - c1.x;
		double dy = c0.y - c1.y;

		x = c0.x - dy;
		y = c0.y + dx;

		return new LineSegment(new Coordinate(c0.x, c0.y), new Coordinate(x, y));
	}

	/**
	 * 做LineString最后一个线段的正交线段。
	 * 
	 * @param line
	 * @return
	 */
	public static LineSegment orthogonalLineOfLinesegment(LineString line) {
		Coordinate[] coords = line.getCoordinates();

		Coordinate c1 = coords[coords.length - 1];
		Coordinate c0 = coords[coords.length - 2];

		return orthogonalLineOfLinesegment(c0, c1);
	}
	
	
	/**
	 * 做LineString第一个线段的正交线段，经过第0点。
	 * 
	 * @param line
	 * @return
	 */
	public static LineSegment orthogonalLineOfLinesegmentA(LineString line) {
		Coordinate[] coords = line.getCoordinates();

		Coordinate c1 = coords[1];
		Coordinate c0 = coords[0];

		return orthogonalLineOfLinesegment(c0, c1);
	}

	/**
	 * 在某线段延长线上按步长以线段首点为中心分别在两边取点。以河道为中心，取断面两边的点。线段首点，实质是河道上的点。
	 * 
	 * @param c1
	 *            线段的首端，也是河道上的点
	 * @param c2
	 *            线段的末端
	 * @param stepLength
	 *            步长
	 * @param stepNum
	 *            步数
	 * @return
	 */
	public static Coordinate[] coordsOfLinesegment(Coordinate c1, Coordinate c2, float stepLength, int stepNum) {
		double dx = c2.x - c1.x;
		double dy = c2.y - c1.y;
		double dc = Math.sqrt(dx * dx + dy * dy);

		Coordinate[] coords = new Coordinate[stepNum * 2 + 1];
		coords[stepNum] = c1;
		double sacle = stepLength / dc;
		for (int i = 1; i <= stepNum; i++) {
			double dyy = i * sacle * dy;
			double dxx = i * sacle * dx;

			coords[stepNum + i] = new Coordinate(c1.x + dxx, c1.y + dyy);
			coords[stepNum - i] = new Coordinate(c1.x - dxx, c1.y - dyy);
		}
		return coords;
	}

	/**
	 * 
	 * @param c1
	 *            开始点
	 * @param c2
	 *            结束点
	 * @param scale
	 *            缩放比例
	 * @return
	 */
	public static Coordinate coordOfLinesegment(Coordinate c1, Coordinate c2, float scale) {
		double dx = c2.x - c1.x;
		double dy = c2.y - c1.y;

		double dyy = scale * dy;
		double dxx = scale * dx;

		return new Coordinate(c1.x + dxx, c1.y + dyy);
	}

	public static void writeShpFile(SimpleFeatureCollection collection, String shpFile, CoordinateReferenceSystem crs) {
		try {
			File newFile = new File(shpFile);
			ShapefileDataStoreFactory dataStoreFactory = new ShapefileDataStoreFactory();
			Map<String, Serializable> params = new HashMap<String, Serializable>();
			params.put("url", newFile.toURI().toURL());
			params.put("create spatial index", Boolean.TRUE);

			ShapefileDataStore newDataStore = (ShapefileDataStore) dataStoreFactory.createNewDataStore(params);

			/*
			 * //定义属性 final SimpleFeatureType TYPE = DataUtilities.createType("Location",
			 * "location:Point," + // <- the geometry attribute: Point type "POIID:String,"
			 * + // <- a String attribute "MESHID:String," + // a number attribute
			 * "OWNER:String" );
			 * 
			 * 
			 * newDataStore.createSchema(TYPE);
			 */
			newDataStore.forceSchemaCRS(crs);

			Transaction transaction = new DefaultTransaction("create");
			String typeName = newDataStore.getTypeNames()[0];
			SimpleFeatureSource featureSource = newDataStore.getFeatureSource(typeName);

			if (featureSource instanceof SimpleFeatureStore) {
				SimpleFeatureStore featureStore = (SimpleFeatureStore) featureSource;
				featureStore.setTransaction(transaction);
				try {
					featureStore.addFeatures(collection);
					transaction.commit();
				} catch (Exception problem) {
					problem.printStackTrace();
					transaction.rollback();
				} finally {
					transaction.close();
				}
			} else {
				System.out.println(typeName + " does not support read/write access");
			}
		} catch (Exception e) {

		}
	}

	/**
	 * Get / create a valid version of the geometry given. If the geometry is a
	 * polygon or multi polygon, self intersections / inconsistencies are fixed.
	 * Otherwise the geometry is returned.
	 * 
	 * @param geom
	 * @return a geometry
	 */
	public static Geometry validate(Geometry geom) {
		if (geom instanceof Polygon) {
			if (geom.isValid()) {
				geom.normalize(); // validate does not pick up rings in the wrong order - this will fix that
				return geom; // If the polygon is valid just return it
			}
			Polygonizer polygonizer = new Polygonizer();
			addPolygon((Polygon) geom, polygonizer);
			return toPolygonGeometry(polygonizer.getPolygons(), geom.getFactory());
		} else if (geom instanceof MultiPolygon) {
			if (geom.isValid()) {
				geom.normalize(); // validate does not pick up rings in the wrong order - this will fix that
				return geom; // If the multipolygon is valid just return it
			}
			Polygonizer polygonizer = new Polygonizer();
			for (int n = geom.getNumGeometries(); n-- > 0;) {
				addPolygon((Polygon) geom.getGeometryN(n), polygonizer);
			}
			return toPolygonGeometry(polygonizer.getPolygons(), geom.getFactory());
		} else {
			return geom; // In my case, I only care about polygon / multipolygon geometries
		}
	}

	/**
	 * Add all line strings from the polygon given to the polygonizer given
	 * 
	 * @param polygon
	 *            polygon from which to extract line strings
	 * @param polygonizer
	 *            polygonizer
	 */
	static void addPolygon(Polygon polygon, Polygonizer polygonizer) {
		addLineString(polygon.getExteriorRing(), polygonizer);
		for (int n = polygon.getNumInteriorRing(); n-- > 0;) {
			addLineString(polygon.getInteriorRingN(n), polygonizer);
		}
	}

	/**
	 * Add the linestring given to the polygonizer
	 * 
	 * @param linestring
	 *            line string
	 * @param polygonizer
	 *            polygonizer
	 */
	static void addLineString(LineString lineString, Polygonizer polygonizer) {

		if (lineString instanceof LinearRing) { // LinearRings are treated differently to line strings : we need a
												// LineString NOT a LinearRing
			lineString = lineString.getFactory().createLineString(lineString.getCoordinateSequence());
		}

		// unioning the linestring with the point makes any self intersections explicit.
		Point point = lineString.getFactory().createPoint(lineString.getCoordinateN(0));
		Geometry toAdd = lineString.union(point);

		// Add result to polygonizer
		polygonizer.add(toAdd);
	}

	/**
	 * Get a geometry from a collection of polygons.
	 * 
	 * @param polygons
	 *            collection
	 * @param factory
	 *            factory to generate MultiPolygon if required
	 * @return null if there were no polygons, the polygon if there was only one, or
	 *         a MultiPolygon containing all polygons otherwise
	 */
	static Geometry toPolygonGeometry(Collection<Polygon> polygons, GeometryFactory factory) {
		switch (polygons.size()) {
		case 0:
			return null; // No valid polygons!
		case 1:
			return polygons.iterator().next(); // single polygon - no need to wrap
		default:
			return factory.createMultiPolygon(polygons.toArray(new Polygon[polygons.size()])); // multiple polygons -
																								// wrap them
		}
	}
	
	/**
	 * 将线段lineSeg延长至maxMeter米
	 * @param lineSeg
	 * @return
	 */
	public static LineSegment prolongLineSegment(LineSegment lineSeg,int maxMeter)
	{
		double dx = lineSeg.p1.x - lineSeg.p0.x;
		double dy = lineSeg.p1.y - lineSeg.p0.y;
		double dc = Math.sqrt(dx * dx + dy * dy);
		double sacle = 1 / dc;			
		double deltX = sacle * dx,deltY =  sacle * dy;
        double xx = maxMeter*deltX;
        double yy = maxMeter*deltY;
		return new LineSegment(lineSeg.p0,new Coordinate(xx,yy));
	}
}

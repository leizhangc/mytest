package com.kygis.algorithmofdam;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.List;
import java.util.stream.IntStream;

import javax.imageio.ImageIO;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridEnvelope2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.coverage.processing.CoverageProcessor;
import org.geotools.data.FeatureSource;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.gce.geotiff.GeoTiffReader;
import org.geotools.gce.geotiff.GeoTiffWriter;
import org.geotools.geometry.DirectPosition2D;
import org.opengis.coverage.PointOutsideCoverageException;
import org.opengis.coverage.grid.GridCoverage;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.geometry.DirectPosition;
import org.opengis.parameter.ParameterValueGroup;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.TransformException;

import com.alibaba.fastjson.JSONArray;
import com.alibaba.fastjson.JSONObject;
import com.kygis.kuibamodel.po.Point;
import com.kygis.kuibamodel.po.Section;
import com.kygis.kuibamodel.po.SectionExactParam;
import com.kygis.kuibamodel.po.WaterDepthMetaInfo;
import com.kygis.kuibamodel.po.WaterSurface;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.Geometry;
import com.vividsolutions.jts.geom.GeometryCollection;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineSegment;
import com.vividsolutions.jts.geom.LineString;
import com.vividsolutions.jts.geom.MultiLineString;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.io.ParseException;

/**
 * @author wheeler
 *
 */
final class DamAnalyseServiceImpl implements DamAnalyseService {
	private String demRaster;// dem文件路径
	private String reverwayShp;// 河道shp文件路径
	private String boundShp;// 河道边界shp文件路径
	private GridCoverage2D coverage;
	private CoordinateReferenceSystem crs;
	private LineString geoObj; // 河道
	private LineString boundObj; // 河道边界
	private Polygon boundObjPolygon;
	private List<Section> sectionList = new ArrayList<Section>();

	private void loadBaseData() throws Exception {
		Logger.info("loading dem file " + demRaster);
		// 加载dem
		File demFile = new File(demRaster);
		GeoTiffReader tifReader = new GeoTiffReader(demFile);
		coverage = tifReader.read(null);
		crs = coverage.getCoordinateReferenceSystem2D();

		Logger.info("loading river file " + reverwayShp);
		// 加载河道
		FeatureIterator<SimpleFeature> itertor = readShp(reverwayShp);
		MultiLineString tmpGeoObj = (MultiLineString) itertor.next().getDefaultGeometry();
		geoObj = (LineString) tmpGeoObj.getGeometryN(0);

		Logger.info("loading bound file " + boundShp);
		// 加载边界
		FeatureIterator<SimpleFeature> itertor1 = readShp(boundShp);
		Geometry tmpBoundObj = (Geometry) itertor1.next().getDefaultGeometry();
		if (tmpBoundObj instanceof MultiLineString) {
			MultiLineString tmpGeoObj1 = (MultiLineString) tmpBoundObj;
			boundObj = (LineString) tmpGeoObj1.getGeometryN(0);
			boundObjPolygon = new GeometryFactory().createPolygon(boundObj.getCoordinates());
		} else if (tmpBoundObj instanceof MultiPolygon) {
			MultiPolygon tmpGeoObj1 = (MultiPolygon) tmpBoundObj;
			boundObjPolygon = (Polygon) tmpGeoObj1.getGeometryN(0);
			boundObj = (LineString) boundObjPolygon.getBoundary();
		}

	}

	private FeatureIterator<SimpleFeature> readShp(String shpFile) throws IOException {
		File riverWayFile = new File(shpFile);
		ShapefileDataStore shpDataStore = null;
		shpDataStore = new ShapefileDataStore(riverWayFile.toURI().toURL());
		String typeName = shpDataStore.getTypeNames()[0];
		FeatureSource<SimpleFeatureType, SimpleFeature> featureSource = null;
		featureSource = (FeatureSource<SimpleFeatureType, SimpleFeature>) shpDataStore.getFeatureSource(typeName);
		FeatureCollection<SimpleFeatureType, SimpleFeature> result = featureSource.getFeatures();

		FeatureIterator<SimpleFeature> itertor = result.features();
		return itertor;
	}

	@Override
	public List<Section> sectionAnalyse(SectionExactParam sectionExactParam) {
		int step = sectionExactParam.getStep();// 断面间距
		int distance = sectionExactParam.getDistance();// 断面上点的间距
		float topDem = sectionExactParam.getTopDem() + 0.1f;// 断面最大高程
		// 等间距取断面点
		// List<LineString> geoObjes =
		// ShpAnalyseTool.splitLineStringIntoParts(geoObj,step);//全河道分析
		// List<LineString> geoObjes = ShpAnalyseTool.splitLineStringIntoParts(geoObj,
		// step,
		// sectionExactParam.getEndPoint());// 指定了出口断面位置。
		List<LineString> geoObjes = ShpAnalyseTool.splitLineStringIntoParts(geoObj, sectionExactParam);

		// List<Section> sectionList = new ArrayList<Section>();
		Logger.info(String.format("在河道上每隔%s米取一个断面点，共取断面个数%s，各断面取的样点坐标及高程如下：", step, geoObjes.size()));
		GeometryFactory geometryFactory = new GeometryFactory();
		int numOfEverySide = sectionExactParam.getPointCount();
		for (int i = 0; i < geoObjes.size(); i++) {
			List<Point> pointsOfSection = new ArrayList<Point>();
			Logger.info(String.format("----第%s个断面，在河道两边，每隔%s米取1个点，共取%s点，连上河道上的点，共%s个点，各点高程如下：", i, distance,
					numOfEverySide, numOfEverySide * 2 + 1));
			LineSegment tmp = ShpAnalyseTool.orthogonalLineOfLinesegment(geoObjes.get(i));

			double dx = tmp.p1.x - tmp.p0.x;
			double dy = tmp.p1.y - tmp.p0.y;
			double dc = Math.sqrt(dx * dx + dy * dy);

			Point[] pointAray = new Point[numOfEverySide * 2 + 1];
			float dem = readDemOfDirectPosition(tmp.p0);
			if (dem > topDem) {
				Logger.error("河道上的点竟然超过了最大高程，属于数据异常，程序退出。");
				return null;
			}
			// Logger.debug(
			// String.format("--------第%04d个点,高程:%f,经度：%s,维度:%s,", numOfEverySide, dem,
			// tmp.p0.x, tmp.p0.y));
			pointAray[numOfEverySide] = new Point(tmp.p0.x, tmp.p0.y, dem);
			double sacle = distance / dc;
			boolean left = true, right = true;
			for (int j = 1; j <= numOfEverySide; j++) {
				double dyy = j * sacle * dy;
				double dxx = j * sacle * dx;

				if (right) {
					Coordinate tmpPoint = new Coordinate(tmp.p0.x + dxx, tmp.p0.y + dyy);
					try {
						float tmpDem = readDemOfDirectPosition(tmpPoint);
						if (tmpDem < topDem) {
							pointAray[numOfEverySide + j] = new Point(tmpPoint.x, tmpPoint.y, tmpDem);
							// Logger.debug(String.format("--------第%04d个点,高程:%f,经度：%s,维度:%s",
							// (numOfEverySide + j),
							// tmpDem, tmpPoint.x, tmpPoint.y));
						} else {
							right = false;
						}
					} catch (PointOutsideCoverageException e) {
						left = false;
					}
				}
				if (left) {
					Coordinate tmpPoint = new Coordinate(tmp.p0.x - dxx, tmp.p0.y - dyy);
					float tmpDem;
					try {
						tmpDem = readDemOfDirectPosition(tmpPoint);

						if (tmpDem < topDem) {
							pointAray[numOfEverySide - j] = new Point(tmpPoint.x, tmpPoint.y, tmpDem);
							// Logger.debug(String.format("--------第%04d个点,高程:%f,经度：%s,维度:%s",
							// (numOfEverySide - j),
							// tmpDem, tmpPoint.x, tmpPoint.y));
						} else {
							left = false;
						}
					} catch (PointOutsideCoverageException e) {
						left = false;
					}
				}
				if ((!left) && (!right)) {
					break;
				}
			}

			for (int m = 0; m < pointAray.length; m++) {
				if (pointAray[m] != null) {
					if (m == numOfEverySide) {
						pointAray[m].setPosiInSection(pointsOfSection.size());
					}
					pointsOfSection.add(pointAray[m]);
				}
			}

			Coordinate[] corrdis = new Coordinate[pointsOfSection.size()];
			for (int m = 0; m < pointsOfSection.size(); m++) {
				Point tmpPoint = pointsOfSection.get(m);
				corrdis[m] = new Coordinate(tmpPoint.getLon(), tmpPoint.getLat());
			}

			Section section = new Section(pointsOfSection);
			section.setNid(String.valueOf(System.nanoTime()));
			section.setSectionParamID(sectionExactParam.getNid());
			LineString geoLine = geometryFactory.createLineString(corrdis);
			section.setLineString(geoLine);
			section.setDistanceOffStart(sectionExactParam.getStep() * (i + 1));
			sectionList.add(section);
			Point riverPoint = new Point(tmp.p0.x, tmp.p0.y,
					readDemOfDirectPosition(new Coordinate(tmp.p0.x, tmp.p0.y)));
			riverPoint.setPosiInSection(pointAray[numOfEverySide].getPosiInSection());
			section.setRiverPoint(riverPoint);
		}
		return sectionList;
	}

	/**
	 * 设定最多的点数maxNum，据河道点riverPoint做垂线，每隔distance取点，作为断面的点，返回点数组。使用时注意数组顺序，是以河道点为开始的。
	 * 
	 * @param riverPoint
	 *            河道点
	 * @param deltX
	 *            代尔塔x
	 * @param deltY
	 *            代尔塔y
	 * @param maxNum
	 *            最大点数，比如distance为1m时，maxNum为12000，就是12公里。
	 * @param distance
	 *            断面上点的平面距离
	 * @param isLeft
	 *            向左还是向右。
	 * @return
	 */
	private Point[] findPointAlongSection(Coordinate riverPoint, double deltX, double deltY, double maxNum,
			double distance, boolean isLeft) {
		int intLeft = isLeft ? -1 : 1;

		GeometryFactory geometryFactory = new GeometryFactory();
		Coordinate coordinate = intersectionWithBoundShp(riverPoint,
				new Coordinate(riverPoint.x + intLeft * maxNum * deltX, riverPoint.y + intLeft * maxNum * deltY));
		if (coordinate == null) {
			return new Point[0];
		}
		com.vividsolutions.jts.geom.Point interPoint = geometryFactory.createPoint(coordinate);
		double allMeter = riverPoint.distance(interPoint.getCoordinate());
		int leftNum = (int) Math.ceil(allMeter / distance);
		Point[] pointArray = new Point[leftNum];
		int length = pointArray.length;
		for (int i = 0; i < length - 1; i++) {
			double dyyy = (i + 1) * deltY;
			double dxxx = (i + 1) * deltX;

			Coordinate tmpPoint = new Coordinate(riverPoint.x + intLeft * dxxx, riverPoint.y + intLeft * dyyy);
			float tmpDem = readDemOfDirectPosition(tmpPoint);
			pointArray[i] = new Point(tmpPoint.x, tmpPoint.y, tmpDem);
		}
		pointArray[leftNum - 1] = new Point(interPoint.getX(), interPoint.getY(),
				readDemOfDirectPosition(interPoint.getCoordinate()));
		return pointArray;
	}

	private Point[] findPointAlongSection(Coordinate stPoint, Coordinate edPoint, double distance) {

		double[] deltXY = computeDeltXY(stPoint, edPoint);
		double allMeter = stPoint.distance(edPoint);
		int leftNum = (int) Math.ceil(allMeter / distance);
		Point[] pointArray = new Point[leftNum];// TODO 少算一个点，应该不影响结果。
		int length = pointArray.length;
		for (int i = 0; i < length; i++) {
			double dyyy = (i) * distance * deltXY[1];
			double dxxx = (i) * distance * deltXY[0];

			Coordinate tmpPoint = new Coordinate(stPoint.x + dxxx, stPoint.y + dyyy);
			float tmpDem = readDemOfDirectPosition(tmpPoint);
			pointArray[i] = new Point(tmpPoint.x, tmpPoint.y, tmpDem);
		}

		return pointArray;
	}

	/**
	 * 通过边界截取断面
	 * 
	 * @param sectionExactParam
	 * @return
	 */
	@Override
	public List<Section> sectionAnalyseWidthBound(SectionExactParam sectionExactParam) {
		int distance = sectionExactParam.getDistance();// 断面上点的间距
		List<LineString> geoObjes = ShpAnalyseTool.splitLineStringIntoParts(geoObj, sectionExactParam);// 河道分段后的折线集合
		GeometryFactory geometryFactory = new GeometryFactory();
		for (int i = 0; i < geoObjes.size(); i++) {
			List<Point> pointsOfSection = new ArrayList<Point>();
			LineSegment tmp = ShpAnalyseTool.orthogonalLineOfLinesegment(geoObjes.get(i));// 折线的正交线段

			double dx = tmp.p1.x - tmp.p0.x;
			double dy = tmp.p1.y - tmp.p0.y;
			double dc = Math.sqrt(dx * dx + dy * dy);
			double sacle = distance / dc;
			double deltX = sacle * dx, deltY = sacle * dy;
			int maxNum = 50000000;

			// 河道点
			float dem = readDemOfDirectPosition(tmp.p0);
			Point riverPoint = new Point(tmp.p0.x, tmp.p0.y, dem);

			// 取河道左侧点，并装入
			Point[] leftPoints = findPointAlongSection(tmp.p0, deltX, deltY, maxNum, distance, true);
			for (int m = leftPoints.length - 1; m > -1; m--) {
				pointsOfSection.add(leftPoints[m]);
			}

			// 装入河道点
			riverPoint.setPosiInSection(leftPoints.length);
			pointsOfSection.add(riverPoint);

			// 取河道右侧点，并装入
			Point[] rightPoints = findPointAlongSection(tmp.p0, deltX, deltY, maxNum, distance, false);
			for (int m = 0; m < rightPoints.length; m++) {
				pointsOfSection.add(rightPoints[m]);
			}

			// 形成坐标集合
			Coordinate[] corrdis = new Coordinate[pointsOfSection.size()];
			for (int m = 0; m < pointsOfSection.size(); m++) {
				Point tmpPoint = pointsOfSection.get(m);
				corrdis[m] = new Coordinate(tmpPoint.getLon(), tmpPoint.getLat());
			}

			Section section = new Section(pointsOfSection);
			section.setNid(String.valueOf(System.nanoTime()));
			section.setSectionParamID(sectionExactParam.getNid());
			LineString geoLine = geometryFactory.createLineString(corrdis);
			section.setLineString(geoLine);
			section.setDistanceOffStart(sectionExactParam.getStep() * (i + 1));
			section.setRiverPoint(riverPoint);
			sectionList.add(section);

		}
		return sectionList;
	}

	/**
	 * 通过边界截取平行断面，以河道为参考点
	 * 
	 * @param sectionExactParam
	 * @return
	 */
	public List<Section> sectionParallelAnalyseWidthBound(SectionExactParam sectionExactParam) {
		int distance = sectionExactParam.getDistance();// 断面上点的间距
		int step = sectionExactParam.getStep();
		LineString aimLine = ShpAnalyseTool.splitLineStringIntoParts(geoObj, sectionExactParam.getStartPoint(),
				sectionExactParam.getEndPoint())[1];
		Coordinate[] allcoods = aimLine.getCoordinates();
		LineSegment aimLineSeg = new LineSegment(allcoods[0], allcoods[1]);// 河道第一条线段。

		double dx = aimLineSeg.p1.x - aimLineSeg.p0.x;
		double dy = aimLineSeg.p1.y - aimLineSeg.p0.y;
		double dc = Math.sqrt(dx * dx + dy * dy);
		double sacle = 1 / dc;
		double deltX = sacle * dx, deltY = sacle * dy;

		GeometryFactory geometryFactory = new GeometryFactory();
		boolean isFinished = false;
		int i = 1;
		while (!isFinished) {
			List<Point> pointsOfSection = new ArrayList<Point>();
			Coordinate p = new Coordinate(aimLineSeg.p0.x + i * step * deltX, aimLineSeg.p0.y + i * step * deltY);
			LineSegment aimLineSegX = new LineSegment(aimLineSeg.p0, p);
			LineSegment tmp1 = ShpAnalyseTool.orthogonalLineOfLinesegment(aimLineSegX.p0, aimLineSegX.p1);// 第一条断面的初始线段。后面把线段延长寻找与河道交点。

			double dx1 = tmp1.p1.x - tmp1.p0.x;
			double dy1 = tmp1.p1.y - tmp1.p0.y;
			double dc1 = Math.sqrt(dx1 * dx1 + dy1 * dy1);
			double sacle1 = 1 / dc1;
			double deltX1 = sacle1 * dx1, deltY1 = sacle1 * dy1;
			int maxNum = 50000000;

			/*
			 * Coordinate riverCoord =
			 * interWithRiverWayShp(geometryFactory.createLineString(new Coordinate[] {
			 * tmp1.p0, new Coordinate(tmp1.p0.x + maxNum * deltX1, tmp1.p0.y + maxNum *
			 * deltY1) })); if (null == riverCoord) { riverCoord =
			 * interWithRiverWayShp(geometryFactory.createLineString(new Coordinate[] {
			 * tmp1.p0, new Coordinate(tmp1.p0.x - maxNum * deltX1, tmp1.p0.y - maxNum *
			 * deltY1) })); }
			 */
			// 在起止点之内
			Coordinate riverCoord = interWithEachOther(aimLine, geometryFactory.createLineString(new Coordinate[] {
					tmp1.p0, new Coordinate(tmp1.p0.x + maxNum * deltX1, tmp1.p0.y + maxNum * deltY1) }));
			if (null == riverCoord) {
				riverCoord = interWithEachOther(aimLine, geometryFactory.createLineString(new Coordinate[] { tmp1.p0,
						new Coordinate(tmp1.p0.x - maxNum * deltX1, tmp1.p0.y - maxNum * deltY1) }));
			}

			if (null == riverCoord) {
				isFinished = true;
				break;
			}
			LineSegment tmp = new LineSegment(riverCoord, tmp1.p1);// 与河道相交的点。

			// 河道点
			float dem = readDemOfDirectPosition(tmp.p0);
			Point riverPoint = new Point(tmp.p0.x, tmp.p0.y, dem);

			// 取河道左侧点，并装入
			Point[] leftPoints = findPointAlongSection(tmp.p0, deltX1, deltY1, maxNum, distance, true);
			for (int m = leftPoints.length - 1; m > -1; m--) {
				pointsOfSection.add(leftPoints[m]);
			}

			// 装入河道点
			riverPoint.setPosiInSection(leftPoints.length);
			pointsOfSection.add(riverPoint);

			// 取河道右侧点，并装入
			Point[] rightPoints = findPointAlongSection(tmp.p0, deltX1, deltY1, maxNum, distance, false);
			for (int m = 0; m < rightPoints.length; m++) {
				pointsOfSection.add(rightPoints[m]);
			}

			// 形成坐标集合
			Coordinate[] corrdis = new Coordinate[pointsOfSection.size()];
			for (int m = 0; m < pointsOfSection.size(); m++) {
				Point tmpPoint = pointsOfSection.get(m);
				corrdis[m] = new Coordinate(tmpPoint.getLon(), tmpPoint.getLat());
			}

			Section section = new Section(pointsOfSection);
			section.setNid(String.valueOf(System.nanoTime()));
			section.setSectionParamID(sectionExactParam.getNid());
			LineString geoLine = geometryFactory.createLineString(corrdis);
			section.setLineString(geoLine);
			section.setDistanceOffStart(sectionExactParam.getStep() * (i));
			section.setRiverPoint(riverPoint);
			sectionList.add(section);
			// System.out.println(section.getLineString());
			i++;
			Logger.debug(section.getLineString().toText());
		}

		return sectionList;
	}

	/**
	 * 通过边界截取平行断面，不考虑河道点
	 * 
	 * @param sectionExactParam
	 * @return
	 */
	public List<Section> sectionParallelAnalyseWidthBoundX(SectionExactParam sectionExactParam) {
		int distance = sectionExactParam.getDistance();// 断面上点的间距
		int step = sectionExactParam.getStep();
		LineString aimLine = ShpAnalyseTool.splitLineStringIntoParts(geoObj, sectionExactParam.getStartPoint(),
				sectionExactParam.getEndPoint())[1];
		Coordinate[] allcoods = aimLine.getCoordinates();
		LineSegment aimLineSeg = new LineSegment(allcoods[0], allcoods[1]);// 河道第一条线段。

		double[] deltXY = computeDeltXY(aimLineSeg.p0, aimLineSeg.p1);

		GeometryFactory geometryFactory = new GeometryFactory();
		boolean isFinished = false;
		int i = 0;
		while (!isFinished) {
			List<Point> pointsOfSection = new ArrayList<Point>();
			Coordinate p = new Coordinate(aimLineSeg.p0.x + (i + 1) * step * deltXY[0],
					aimLineSeg.p0.y + (i + 1) * step * deltXY[1]);
			LineSegment aimLineSegX = new LineSegment(aimLineSeg.p0, p);

			LineSegment tmp1 = ShpAnalyseTool.orthogonalLineOfLinesegment(aimLineSegX.p0, aimLineSegX.p1);// 每条断面的初始线段。后面把线段延长寻找与边界交点。
			double[] deltXY1 = computeDeltXY(tmp1.p0, tmp1.p1);
			int maxNum = 5000000;
			// 向左延长寻找交点
			Coordinate[] leftInterPoint = this.intersectionWithBoundShpX(tmp1.p0,
					new Coordinate(tmp1.p0.x - maxNum * deltXY1[0], tmp1.p0.y - maxNum * deltXY1[1]));

			// 向右延长寻找交点
			Coordinate[] rightInterPoint = this.intersectionWithBoundShpX(tmp1.p0,
					new Coordinate(tmp1.p0.x + maxNum * deltXY1[0], tmp1.p0.y + maxNum * deltXY1[1]));

			if (null == leftInterPoint && null == rightInterPoint) {
				isFinished = true;
				break;
			}
			Coordinate[] leftAndRight = null;
			if (deltXY1[0] > 0) {
				leftAndRight = combineLeftAndRight(leftInterPoint, rightInterPoint);
			} else if (deltXY1[0] < 0) {
				leftAndRight = combineLeftAndRight(rightInterPoint, leftInterPoint);
			}
			// System.out.println(i+"...1..."+aimLineSegX);
			// System.out.println(i+"...2..."+geometryFactory.createLineString(leftAndRight));
			// System.out.println(i+"...3..."+tmp1);

			for (int m = 1; m < leftAndRight.length; m++) {
				Coordinate innerSt = leftAndRight[m - 1];
				Coordinate innerEd = leftAndRight[m];
				LineSegment tmpSeg = new LineSegment(innerSt, innerEd);
				if (boundObjPolygon.contains(geometryFactory.createPoint(tmpSeg.midPoint()))) {
					Point[] points = findPointAlongSection(innerSt, innerEd, distance);
					for (int n = 0; n < points.length; n++) {
						pointsOfSection.add(points[n]);
					}
				}
			}

			// 形成坐标集合
			Coordinate[] corrdis = new Coordinate[pointsOfSection.size()];
			for (int m = 0; m < pointsOfSection.size(); m++) {
				Point tmpPoint = pointsOfSection.get(m);
				corrdis[m] = new Coordinate(tmpPoint.getLon(), tmpPoint.getLat());
			}

			Section section = new Section(pointsOfSection);
			section.setNid(String.valueOf(System.nanoTime()));
			section.setSectionParamID(sectionExactParam.getNid());
			LineString geoLine = geometryFactory.createLineString(corrdis);
			section.setLineString(geoLine);
			section.setDistanceOffStart(sectionExactParam.getStep() * (i + 1));
			sectionList.add(section);
			i++;
		}

		return sectionList;
	}

	public List<Section> sectionParallelAnalyseWidthBoundXX(SectionExactParam sectionExactParam) {
		int distance = sectionExactParam.getDistance();// 断面上点的间距
		int step = sectionExactParam.getStep();

		GeometryFactory geometryFactory = new GeometryFactory();
		Coordinate stPoint = new Coordinate(sectionExactParam.getStartPoint()[0], sectionExactParam.getStartPoint()[1]);
		Coordinate edPoint = new Coordinate(sectionExactParam.getEndPoint()[0], sectionExactParam.getEndPoint()[1]);// 出口断面坐标
		double[] deltXY = computeDeltXY(stPoint, edPoint);
		double allMeter = stPoint.distance(edPoint);
		int leftNum = (int) Math.floor(allMeter / step);
		Point[] pointArray = new Point[leftNum];
		int length = pointArray.length;
		for (int i = 0; i < length; i++) {
			List<Point> pointsOfSection = new ArrayList<Point>();

			Coordinate tmpPoint = new Coordinate(stPoint.x + (i + 1) * step * deltXY[0],
					stPoint.y + (i + 1) * step * deltXY[1]);

			LineSegment aimLineSegX = new LineSegment(stPoint, tmpPoint);

			LineSegment tmp1 = ShpAnalyseTool.orthogonalLineOfLinesegment(aimLineSegX.p0, aimLineSegX.p1);// 每条断面的初始线段。后面把线段延长寻找与边界交点。
			double[] deltXY1 = computeDeltXY(tmp1.p0, tmp1.p1);
			int maxNum = 5000000;
			// 向左延长寻找交点
			Coordinate[] leftInterPoint = this.intersectionWithBoundShpX(tmp1.p0,
					new Coordinate(tmp1.p0.x - maxNum * deltXY1[0], tmp1.p0.y - maxNum * deltXY1[1]));

			// 向右延长寻找交点
			Coordinate[] rightInterPoint = this.intersectionWithBoundShpX(tmp1.p0,
					new Coordinate(tmp1.p0.x + maxNum * deltXY1[0], tmp1.p0.y + maxNum * deltXY1[1]));

			if (null == leftInterPoint && null == rightInterPoint) {
				break;
			}
			Coordinate[] leftAndRight = null;
			if (deltXY1[0] > 0) {
				leftAndRight = combineLeftAndRight(leftInterPoint, rightInterPoint);
			} else if (deltXY1[0] < 0) {
				leftAndRight = combineLeftAndRight(rightInterPoint, leftInterPoint);
			}
			// System.out.println(i+"...1..."+aimLineSegX);
			// System.out.println(i+"...2..."+geometryFactory.createLineString(leftAndRight));
			// System.out.println(i+"...3..."+tmp1);

			for (int m = 1; m < leftAndRight.length; m++) {
				Coordinate innerSt = leftAndRight[m - 1];
				Coordinate innerEd = leftAndRight[m];
				LineSegment tmpSeg = new LineSegment(innerSt, innerEd);
				if (boundObjPolygon.contains(geometryFactory.createPoint(tmpSeg.midPoint()))) {
					Point[] points = findPointAlongSection(innerSt, innerEd, distance);
					for (int n = 0; n < points.length; n++) {
						pointsOfSection.add(points[n]);
					}
				}
			}

			// 形成坐标集合
			Coordinate[] corrdis = new Coordinate[pointsOfSection.size()];
			Point riverPoint = null;
			List<Point> riverPoints = new ArrayList<>();
			for (int m = 0; m < pointsOfSection.size(); m++) {
				Point tmpPoint1 = pointsOfSection.get(m);
				corrdis[m] = new Coordinate(tmpPoint1.getLon(), tmpPoint1.getLat());
				Geometry p = geometryFactory.createPoint(corrdis[m]);
				// if (riverPoint == null && geoObj.distance(p) < distance) {
				if (geoObj.distance(p) < distance) {
					Point tPoint = new Point(tmpPoint1.getLon(), tmpPoint1.getLat(), tmpPoint1.getDem());
					tPoint.setPosiInSection(m);
					tPoint.tag = geoObj.distance(p);
					riverPoints.add(tPoint);
				}
			}
			int altCount = riverPoints.size();
			if (altCount == 1) {
				riverPoint = riverPoints.get(0);
			} else if (altCount > 1) {
				riverPoints.sort(new Comparator<Point>() {

					@Override
					public int compare(Point o1, Point o2) {
						return Double.compare(Double.parseDouble(o1.tag + ""), Double.parseDouble(o2.tag + ""));
					}
				});
				riverPoint = riverPoints.get(0);
			}

			Section section = new Section(pointsOfSection);
			section.setRiverPoint(riverPoint);
			section.setNid(String.valueOf(System.nanoTime()));
			section.setSectionParamID(sectionExactParam.getNid());
			LineString geoLine = geometryFactory.createLineString(corrdis);
			section.setLineString(geoLine);
			section.setDistanceOffStart(sectionExactParam.getStep() * (i + 1));
			sectionList.add(section);
		}
		return sectionList;
	}

	/**
	 * 计算deltx 和 delty
	 * 
	 * @param start
	 * @param end
	 * @param distance
	 * @return
	 */
	private double[] computeDeltXY(Coordinate start, Coordinate end) {
		double dx = end.x - start.x;
		double dy = end.y - start.y;
		double dc = Math.sqrt(dx * dx + dy * dy);
		double sacle = 1 / dc;
		double deltX = sacle * dx, deltY = sacle * dy;
		return new double[] { deltX, deltY };
	}

	private Coordinate[] combineLeftAndRight(Coordinate[] leftInterPoint, Coordinate[] rightInterPoint) {
		if (null == leftInterPoint && null == rightInterPoint)
			return null;

		if (null == leftInterPoint)
			return rightInterPoint;

		if (null == rightInterPoint)
			return leftInterPoint;

		Coordinate[] leftAndRight = null;
		if (leftInterPoint[leftInterPoint.length - 1].equals(rightInterPoint[0])) {
			leftAndRight = new Coordinate[leftInterPoint.length + rightInterPoint.length - 1];
			System.arraycopy(leftInterPoint, 0, leftAndRight, 0, leftInterPoint.length - 1);
			System.arraycopy(rightInterPoint, 0, leftAndRight, leftInterPoint.length - 1, rightInterPoint.length);
		} else {
			leftAndRight = new Coordinate[leftInterPoint.length + rightInterPoint.length];
			System.arraycopy(leftInterPoint, 0, leftAndRight, 0, leftInterPoint.length);
			System.arraycopy(rightInterPoint, 0, leftAndRight, leftInterPoint.length, rightInterPoint.length);
		}
		return leftAndRight;
	}

	/**
	 * 判断与边界的交点 返回一个点
	 * 
	 * @param p0
	 *            起始点
	 * @param p1
	 *            终点
	 * @return
	 */
	private Coordinate intersectionWithBoundShp(Coordinate p0, Coordinate p1) {
		LineString tmpLineStr = new GeometryFactory().createLineString(new Coordinate[] { p0, p1 });

		Geometry tmpPo = boundObj.intersection(tmpLineStr);
		if (tmpPo instanceof com.vividsolutions.jts.geom.Point) {
			return ((com.vividsolutions.jts.geom.Point) tmpPo).getCoordinate();
		} else if (tmpPo instanceof com.vividsolutions.jts.geom.MultiPoint) {
			return ((com.vividsolutions.jts.geom.MultiPoint) tmpPo).getCoordinate();
		}
		return null;
	}

	/**
	 * 判断与边界的交点 返回全部的点
	 * 
	 * @param p0
	 * @param p1
	 * @return
	 */
	private Coordinate[] intersectionWithBoundShpX(Coordinate p0, Coordinate p1) {
		LineString tmpLineStr = new GeometryFactory().createLineString(new Coordinate[] { p0, p1 });

		Geometry tmpPo = boundObj.intersection(tmpLineStr);
		if (tmpPo instanceof com.vividsolutions.jts.geom.Point) {
			return new Coordinate[] { ((com.vividsolutions.jts.geom.Point) tmpPo).getCoordinate() };
		} else if (tmpPo instanceof com.vividsolutions.jts.geom.MultiPoint) {
			return ((com.vividsolutions.jts.geom.MultiPoint) tmpPo).getCoordinates();
		}
		return null;
	}

	private Coordinate interWithRiverWayShp(LineString tmpLineStr) {
		Geometry tmpPo = this.geoObj.intersection(tmpLineStr);
		if (tmpPo instanceof com.vividsolutions.jts.geom.Point) {
			return ((com.vividsolutions.jts.geom.Point) tmpPo).getCoordinate();
		} else if (tmpPo instanceof com.vividsolutions.jts.geom.MultiPoint) {
			return ((com.vividsolutions.jts.geom.MultiPoint) tmpPo).getCoordinate();
		}
		return null;
	}

	private Coordinate interWithEachOther(LineString subRiverLineStr, LineString tmpLineStr) {
		Geometry tmpPo = subRiverLineStr.intersection(tmpLineStr);
		if (tmpPo instanceof com.vividsolutions.jts.geom.Point) {
			return ((com.vividsolutions.jts.geom.Point) tmpPo).getCoordinate();
		} else if (tmpPo instanceof com.vividsolutions.jts.geom.MultiPoint) {
			return ((com.vividsolutions.jts.geom.MultiPoint) tmpPo).getCoordinate();
		}
		return null;
	}

	private float readDemOfDirectPosition(Coordinate p) {
		return readDemOfDirectPosition(p.x, p.y);
	}

	private float readDemOfDirectPosition(double lon, double lat) {
		DirectPosition position = new DirectPosition2D(crs, lon, lat);

		float[] results = (float[]) coverage.evaluate(position); // assume float
		// resample with the same array
		results = coverage.evaluate(position, results);
		return results[0];
	}

	@Override
	public WaterSurface waterSurfaceAnalyse(double[] waterDemArray, List<Section> allSection) {
		try {
			Coordinate[] coords = DestoriedDamAlgorithm.waterSurfaceCoordinates1(waterDemArray, allSection);
			GeometryFactory geometryFactory = new GeometryFactory();
			Polygon waterSurfaceWkt = geometryFactory.createPolygon(coords);
			Logger.debug(waterSurfaceWkt.toText());

			Polygon theBiggest = (Polygon) waterSurfaceWkt.convexHull();

			Polygon minPoly = null;
			Geometry mutilGeos = ShpAnalyseTool.validate(waterSurfaceWkt);
			if (mutilGeos instanceof MultiPolygon) {
				MultiPolygon multiObj = ((MultiPolygon) mutilGeos);
				int size = multiObj.getNumGeometries();
				minPoly = (Polygon) multiObj.getGeometryN(0);
				for (int i = 1; i < size; i++) {
					Polygon the = (Polygon) multiObj.getGeometryN(i);
					if (the.getArea() > minPoly.getArea()) {
						minPoly = the;
					}
				}
			} else {
				minPoly = (Polygon) mutilGeos;
			}

			WaterSurface waterSurface = new WaterSurface();
			waterSurface.setNid(LocalDateTime.now().format(DateTimeFormatter.ofPattern("yyyyMMddHHmmss")));
			Section tmp = allSection.get(0);
			String damParamID = tmp.getDamParamID();
			String sectionParamID = tmp.getSectionParamID();
			waterSurface.setDamParamID(damParamID);
			waterSurface.setSectionParamID(sectionParamID);
			waterSurface.setWkt(minPoly);
			waterSurface.setWktConVex(theBiggest);
			waterSurface.setWktOriginal(waterSurfaceWkt);
			return waterSurface;
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}

	public Polygon waterSurfaceBoundaryAnalyse(double[] waterDemArray, List<Section> allSection) {
		try {
			Coordinate[] coords = DestoriedDamAlgorithm.waterSurfaceCoordinates1(waterDemArray, allSection);
			GeometryFactory geometryFactory = new GeometryFactory();
			Polygon waterSurfaceWkt = (Polygon) geometryFactory.createPolygon(coords).convexHull();
			return waterSurfaceWkt;
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;

	}

	public WaterDepthMetaInfo waterDepthAnalyse(WaterSurface waterSurface, float bufferDistance,
			String storedAsFileName) {
		CoverageProcessor processor = new CoverageProcessor();
		ParameterValueGroup params = processor.getOperation("CoverageCrop").getParameters();
		params.parameter("Source").setValue(coverage);
		// params.parameter("ENVELOPE").setValue(bounds);
		Polygon bigBoundary = waterSurface.getWktConVex();
		params.parameter("ROI").setValue(bigBoundary);
		params.parameter("ForceMosaic").setValue(true);

		GridCoverage2D res = (GridCoverage2D) processor.doOperation(params);
		GridGeometry2D geometry = res.getGridGeometry();

		RenderedImage img = res.getRenderedImage();
		Raster raster = img.getData();
		int width = raster.getWidth();
		int height = raster.getHeight();

		BufferedImage resultImg = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

		ColorModel cm = img.getColorModel();
		WritableRaster writableRaster = cm.createCompatibleWritableRaster(width, height);
		Hashtable<String, Object> properties = new Hashtable<String, Object>();
		String[] keys = img.getPropertyNames();
		if (keys != null) {
			for (int i = 0; i < keys.length; i++) {
				properties.put(keys[i], img.getProperty(keys[i]));
			}
		}
		BufferedImage resultImg_gray = new BufferedImage(cm, writableRaster, true, null);

		int minX = raster.getMinX();
		int minY = raster.getMinY();
		float maxDepth = 0, minDepth = 0;
		try {
			GeometryFactory geometryFactory = new GeometryFactory();
			for (int i = 0; i < width; i++) {
				final int ii = i;
				IntStream.rangeClosed(0, height - 1).parallel().forEach((j) -> {
					org.geotools.geometry.Envelope2D pixelEnvelop = null;
					try {
						pixelEnvelop = geometry.gridToWorld(new GridEnvelope2D(ii + minX, j + minY, 1, 1));
					} catch (TransformException e) {
						e.printStackTrace();
					}
					double lat = pixelEnvelop.getCenterY();
					double lon = pixelEnvelop.getCenterX();
					float x = ((float[]) raster.getDataElements(ii + minX, j + minY, (float[]) null))[0];
					com.vividsolutions.jts.geom.Point p = geometryFactory.createPoint(new Coordinate(lon, lat));
					if (x == 0 || !bigBoundary.contains(p)) {
						Color color = new Color(0, 0, 0, 0);
						resultImg.setRGB(ii, j, color.getRGB());
						resultImg_gray.setRGB(ii, j, color.getRGB());
					} else {
						double maxZ = 0;
						try {
							maxZ = findMaxZOfSection2(new Coordinate(lon, lat), bufferDistance);
						} catch (ParseException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						float watDepth = (float) (maxZ - x);

						/*
						 * if (watDepth < minDepth) { minDepth = watDepth; } else if (watDepth >
						 * maxDepth) { maxDepth = watDepth; }
						 */
						double[] level = new double[] { 2, 4, 7, 10, 16 };
						// double[] level = new double[] {0.5,1,1.5,2.5,5.0};
						Color color = null;
						if (watDepth <= 0) {
							color = new Color(0, 0, 0, 0);
							watDepth = 0;

						} else if (watDepth <= level[0]) {
							color = new Color(230, 255, 255);
						} else if (watDepth > level[0] && watDepth <= level[1]) {
							color = new Color(110, 255, 255);
						} else if (watDepth > level[1] && watDepth <= level[2]) {
							color = new Color(0, 230, 230);
						} else if (watDepth > level[2] && watDepth <= level[3]) {
							color = new Color(0, 191, 191);
						} else if (watDepth > level[3] && watDepth <= level[4]) {
							color = new Color(0, 147, 217);
						} else if (watDepth > level[4]) {
							color = new Color(0, 80, 170);
						}
						resultImg.setRGB(ii, j, color.getRGB());
						writableRaster.setPixel(ii, j, new float[] { watDepth });
					}
				}

				);
			}

			GridCoverage2D tmp = RasterAnalyseTool.bufferedImageToGridCoverage2D(resultImg_gray, res.getEnvelope());
			// 写tif
			File tst = new File(storedAsFileName + ".tif");
			tst.createNewFile();
			GeoTiffWriter writer;
			writer = new GeoTiffWriter(tst);
			writer.write((GridCoverage) tmp, null);
			writer.dispose();

			// 写png
			File tst1 = new File(storedAsFileName + ".png");
			tst1.createNewFile();
			ImageIO.write(resultImg, "png", tst1);

			WaterDepthMetaInfo waterDepthMetaInfo = new WaterDepthMetaInfo();
			waterDepthMetaInfo.setNid(String.valueOf(System.nanoTime()));
			waterDepthMetaInfo.setWaterSurfaceID(waterSurface.getNid());
			waterDepthMetaInfo.setDamParamID(waterSurface.getDamParamID());
			waterDepthMetaInfo.setSectionParamID(waterSurface.getSectionParamID());
			// waterDepthMetaInfo.setMaxDepth(maxDepth);
			// waterDepthMetaInfo.setMinDepth(minDepth);
			waterDepthMetaInfo.setPath(storedAsFileName);
			waterDepthMetaInfo.setArea(waterSurface.getWkt().getArea());
			Envelope rect = waterSurface.getWkt().getEnvelopeInternal();
			waterDepthMetaInfo.setMinX(rect.getMinX());
			waterDepthMetaInfo.setMinY(rect.getMinY());
			waterDepthMetaInfo.setMaxX(rect.getMaxX());
			waterDepthMetaInfo.setMaxY(rect.getMaxY());
			return waterDepthMetaInfo;
		} catch (IOException e) {
			e.printStackTrace();
		} // catch (TransformException e) {
			// // TODO Auto-generated catch block
			// e.printStackTrace();
			// } catch (ParseException e) {
			// // TODO Auto-generated catch block
			// e.printStackTrace();
			// }
		return null;
	}

	/**
	 * 两断面之间插值来计算水位，推算水深。
	 * 
	 * @param sectionExactParam
	 * 
	 * @param waterSurface
	 * @param storedAsFileName
	 * @return
	 */
	public WaterDepthMetaInfo waterDepthAnalyse(SectionExactParam sectionExactParam, WaterSurface waterSurface,
			String storedAsFileName, JSONArray colorLevel1) {

		if (null == colorLevel1)
			colorLevel1 = getColorLevel();
		GeometryFactory geometryFactory = new GeometryFactory();
		final JSONArray colorLevel = colorLevel1;
		CoverageProcessor processor = new CoverageProcessor();
		ParameterValueGroup params = processor.getOperation("CoverageCrop").getParameters();
		params.parameter("Source").setValue(coverage);
		// params.parameter("ENVELOPE").setValue(bounds);
		Polygon bigBoundary = waterSurface.getWktConVex();
		Geometry spGeometry = geometryFactory.createPoint(
				new Coordinate(sectionExactParam.getStartPoint()[0], sectionExactParam.getStartPoint()[1]));
		Geometry epGeometry = geometryFactory
				.createPoint(new Coordinate(sectionExactParam.getEndPoint()[0], sectionExactParam.getEndPoint()[1]));
		Geometry unionGeometry = bigBoundary.union(spGeometry);
		if (unionGeometry instanceof GeometryCollection) {
			unionGeometry = ((GeometryCollection) unionGeometry).convexHull();
		}
		unionGeometry = unionGeometry.union(epGeometry);
		if (unionGeometry instanceof GeometryCollection) {
			unionGeometry = ((GeometryCollection) unionGeometry).convexHull();
		}
		params.parameter("ROI").setValue(unionGeometry);
		// params.parameter("ROI").setValue(boundObjPolygon);
		params.parameter("ForceMosaic").setValue(true);
		Logger.trace("ConVex : " + bigBoundary.toText());
		Logger.trace("union : " + unionGeometry.toText());
		GridCoverage2D res = (GridCoverage2D) processor.doOperation(params);
		GridGeometry2D geometry = res.getGridGeometry();

		RenderedImage img = res.getRenderedImage();
		Raster raster = img.getData();
		int width = raster.getWidth();
		int height = raster.getHeight();

		BufferedImage resultImg = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

		ColorModel cm = img.getColorModel();
		WritableRaster writableRaster = cm.createCompatibleWritableRaster(width, height);
		Hashtable<String, Object> properties = new Hashtable<String, Object>();
		String[] keys = img.getPropertyNames();
		if (keys != null) {
			for (int i = 0; i < keys.length; i++) {
				properties.put(keys[i], img.getProperty(keys[i]));
			}
		}
		BufferedImage resultImg_gray = new BufferedImage(cm, writableRaster, true, null);

		int minX = raster.getMinX();
		int minY = raster.getMinY();
		// List<Float> values = new ArrayList<>();
		float maxDepth = Float.MIN_VALUE;
		float minDepth = Float.MAX_VALUE;
		try {
			for (int i = 0; i < width; i++) {
				final int ii = i;
				final Geometry bgBoundary = unionGeometry;
				IntStream.rangeClosed(0, height - 1).parallel().forEach((j) -> {
					org.geotools.geometry.Envelope2D pixelEnvelop = null;
					try {
						pixelEnvelop = geometry.gridToWorld(new GridEnvelope2D(ii + minX, j + minY, 1, 1));
					} catch (TransformException e) {
						e.printStackTrace();
					}
					double lat = pixelEnvelop.getCenterY();
					double lon = pixelEnvelop.getCenterX();
					float x = ((float[]) raster.getDataElements(ii + minX, j + minY, (float[]) null))[0];
					com.vividsolutions.jts.geom.Point p = geometryFactory.createPoint(new Coordinate(lon, lat));
					if (Float.compare(x, 0f) == 0 || Float.compare(x, Float.MIN_VALUE) == 0
							|| Float.compare(x, Float.MAX_VALUE) == 0 || Float.compare(x, -Float.MAX_VALUE) == 0
							|| !boundObjPolygon.contains(p) || !bgBoundary.contains(p)) {
						Color color = new Color(0, 0, 0, 0);
						resultImg.setRGB(ii, j, color.getRGB());
						resultImg_gray.setRGB(ii, j, color.getRGB());
					} else {
						double maxZ = 0;
						try {
							maxZ = findMaxZOfSection(new Coordinate(lon, lat), sectionExactParam);
						} catch (ParseException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						float watDepth = (float) (maxZ - x);
						double[] level = new double[] { 2, 4, 7, 10, 16 };
						// double[] level = new double[] {0.5,1,1.5,2.5,5.0};
						Color color = null;
						if (watDepth <= 0) {
							color = new Color(0, 0, 0, 0);
							watDepth = 0;

						} else if (watDepth <= level[0]) {
							color = new Color(230, 255, 255);
						} else if (watDepth > level[0] && watDepth <= level[1]) {
							color = new Color(110, 255, 255);
						} else if (watDepth > level[1] && watDepth <= level[2]) {
							color = new Color(0, 230, 230);
						} else if (watDepth > level[2] && watDepth <= level[3]) {
							color = new Color(0, 191, 191);
						} else if (watDepth > level[3] && watDepth <= level[4]) {
							color = new Color(0, 147, 217);
						} else if (watDepth > level[4]) {
							color = new Color(0, 80, 170);
						}
						// values.add(watDepth);
						/*
						 * Color color = null; if (watDepth <= 0) { color = new Color(0, 0, 0, 0);
						 * watDepth = 0; } else { int size = colorLevel.size(); float valst =
						 * ((JSONObject) colorLevel.get(0)).getFloat("value"); float valed =
						 * ((JSONObject) colorLevel.get(size - 1)).getFloat("value"); if (watDepth <=
						 * valst) { String[] rgb = ((JSONObject)
						 * colorLevel.get(0)).getString("color").split(","); color = changeColor(new
						 * String[] { "255", "255", "255" }, rgb, valst, watDepth); } else if (watDepth
						 * > valed) { String[] rgb = ((JSONObject) colorLevel.get(size -
						 * 1)).getString("color").split(","); int R = Integer.valueOf(rgb[0]), G =
						 * Integer.valueOf(rgb[1]), B = Integer.valueOf(rgb[2]); color = new Color(R, G,
						 * B); } else { for (int ci = 1; ci < size - 1; ci++) { float val0 =
						 * ((JSONObject) colorLevel.get(ci - 1)).getFloat("value"); float val1 =
						 * ((JSONObject) colorLevel.get(ci)).getFloat("value"); if (watDepth > val0 &&
						 * watDepth <= val1) { String[] rgb0 = ((JSONObject) colorLevel.get(ci -
						 * 1)).getString("color") .split(","); String[] rgb1 = ((JSONObject)
						 * colorLevel.get(ci)).getString("color").split(","); color = changeColor(rgb0,
						 * rgb1, val1, watDepth); } } } }
						 */
						resultImg.setRGB(ii, j, color.getRGB());
						writableRaster.setPixel(ii, j, new float[] { watDepth });
					}
				}

				);
			}
			// values.sort(new Comparator<Float>() {
			// @Override
			// public int compare(Float o1, Float o2) {
			// return o1.compareTo(o2);
			// }
			//
			// });
			// minDepth = values.get(0);
			// maxDepth = values.get(values.size() - 1);
			GridCoverage2D tmp = RasterAnalyseTool.bufferedImageToGridCoverage2D(resultImg_gray, res.getEnvelope());
			// 写tif
			File tst = new File(storedAsFileName + ".tif");
			tst.createNewFile();
			GeoTiffWriter writer;
			writer = new GeoTiffWriter(tst);
			// writer.setMetadataValue(Integer.toString(BaselineTIFFTagSet.TAG_MAX_SAMPLE_VALUE),
			// maxDepth + "");
			// writer.setMetadataValue(Integer.toString(BaselineTIFFTagSet.TAG_MIN_SAMPLE_VALUE),
			// minDepth + "");
			writer.write((GridCoverage) tmp, null);
			writer.dispose();

			// 写png
			File tst1 = new File(storedAsFileName + ".png");
			tst1.createNewFile();
			ImageIO.write(resultImg, "png", tst1);

			WaterDepthMetaInfo waterDepthMetaInfo = new WaterDepthMetaInfo();
			waterDepthMetaInfo.setNid(String.valueOf(System.nanoTime()));
			waterDepthMetaInfo.setWaterSurfaceID(waterSurface.getNid());
			waterDepthMetaInfo.setDamParamID(waterSurface.getDamParamID());
			waterDepthMetaInfo.setSectionParamID(waterSurface.getSectionParamID());
			// waterDepthMetaInfo.setMaxDepth(maxDepth);
			// waterDepthMetaInfo.setMinDepth(minDepth);
			waterDepthMetaInfo.setPath(storedAsFileName);
			waterDepthMetaInfo.setArea(waterSurface.getWkt().getArea());
			Envelope rect = waterSurface.getWkt().getEnvelopeInternal();
			waterDepthMetaInfo.setMinX(rect.getMinX());
			waterDepthMetaInfo.setMinY(rect.getMinY());
			waterDepthMetaInfo.setMaxX(rect.getMaxX());
			waterDepthMetaInfo.setMaxY(rect.getMaxY());
			return waterDepthMetaInfo;
		} catch (IOException e) {
			e.printStackTrace();
		} // catch (TransformException e) {
			// // TODO Auto-generated catch block
			// e.printStackTrace();
			// } catch (ParseException e) {
			// // TODO Auto-generated catch block
			// e.printStackTrace();
			// }
		return null;
	}

	/**
	 * @param waterSurface
	 * @param bufferDistance
	 * @param storedAsFileName
	 * @param colorLevel
	 *            [{value:color:}]
	 * @return
	 */
	public WaterDepthMetaInfo waterDepthAnalyse(WaterSurface waterSurface, float bufferDistance,
			String storedAsFileName, JSONArray colorLevel) {
		if (null == colorLevel)
			colorLevel = getColorLevel();
		final JSONArray colorLevels = colorLevel;
		CoverageProcessor processor = new CoverageProcessor();
		ParameterValueGroup params = processor.getOperation("CoverageCrop").getParameters();
		params.parameter("Source").setValue(coverage);
		// params.parameter("ENVELOPE").setValue(bounds);
		Polygon bigBoundary = waterSurface.getWktConVex();
		Geometry buffered = bigBoundary.buffer(100);
		params.parameter("ROI").setValue(buffered);
		params.parameter("ForceMosaic").setValue(true);

		GridCoverage2D res = (GridCoverage2D) processor.doOperation(params);
		GridGeometry2D geometry = res.getGridGeometry();

		RenderedImage img = res.getRenderedImage();
		Raster raster = img.getData();
		int width = raster.getWidth();
		int height = raster.getHeight();

		BufferedImage resultImg = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);

		ColorModel cm = img.getColorModel();
		WritableRaster writableRaster = cm.createCompatibleWritableRaster(width, height);
		Hashtable<String, Object> properties = new Hashtable<String, Object>();
		String[] keys = img.getPropertyNames();
		if (keys != null) {
			for (int i = 0; i < keys.length; i++) {
				properties.put(keys[i], img.getProperty(keys[i]));
			}
		}
		BufferedImage resultImg_gray = new BufferedImage(cm, writableRaster, true, null);

		int minX = raster.getMinX();
		int minY = raster.getMinY();
		float maxDepth = 0, minDepth = 0;
		try {
			GeometryFactory geometryFactory = new GeometryFactory();
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {

					org.geotools.geometry.Envelope2D pixelEnvelop = null;
					try {
						pixelEnvelop = geometry.gridToWorld(new GridEnvelope2D(i + minX, j + minY, 1, 1));
					} catch (TransformException e) {
						e.printStackTrace();
					}
					double lat = pixelEnvelop.getCenterY();
					double lon = pixelEnvelop.getCenterX();
					float x = ((float[]) raster.getDataElements(i + minX, j + minY, (float[]) null))[0];
					com.vividsolutions.jts.geom.Point p = geometryFactory.createPoint(new Coordinate(lon, lat));
					if (x == 0 || !bigBoundary.contains(p)) {
						Color color = new Color(0, 0, 0, 0);
						resultImg.setRGB(i, j, color.getRGB());
						resultImg_gray.setRGB(i, j, color.getRGB());
					} else {
						double maxZ = 0;
						try {
							maxZ = findMaxZOfSection2(new Coordinate(lon, lat), bufferDistance);
						} catch (ParseException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						float watDepth = (float) (maxZ - x);

						if (watDepth < minDepth) {
							minDepth = watDepth;
						} else if (watDepth > maxDepth) {
							maxDepth = watDepth;
						}

						Color color = null;
						if (watDepth <= 0) {
							color = new Color(0, 0, 0, 0);
							watDepth = 0;
						} else {
							int size = colorLevel.size();
							float valst = ((JSONObject) colorLevel.get(0)).getFloat("value");
							float valed = ((JSONObject) colorLevel.get(size - 1)).getFloat("value");
							if (watDepth <= valst) {
								String[] rgb = ((JSONObject) colorLevel.get(0)).getString("color").split(",");
								color = changeColor(new String[] { "255", "255", "255" }, rgb, valst, watDepth);
							} else if (watDepth > valed) {
								String[] rgb = ((JSONObject) colorLevel.get(size - 1)).getString("color").split(",");
								int R = Integer.valueOf(rgb[0]), G = Integer.valueOf(rgb[1]),
										B = Integer.valueOf(rgb[2]);
								color = new Color(R, G, B);
							} else {
								for (int ci = 1; ci < size - 1; ci++) {
									float val0 = ((JSONObject) colorLevel.get(ci - 1)).getFloat("value");
									float val1 = ((JSONObject) colorLevel.get(ci)).getFloat("value");
									if (watDepth > val0 && watDepth <= val1) {
										String[] rgb0 = ((JSONObject) colorLevel.get(ci - 1)).getString("color")
												.split(",");
										String[] rgb1 = ((JSONObject) colorLevel.get(ci)).getString("color").split(",");
										color = changeColor(rgb0, rgb1, val1, watDepth);
										/*
										 * int R = Integer.valueOf(rgb[0]),G =
										 * Integer.valueOf(rgb[1]),B=Integer.valueOf(rgb[2]); color = new Color(R, G,
										 * B); while(watDepth<=(val1=-0.5f)&&watDepth > val0) {
										 * if(R+colv>=255||G+colv>=255|| B+colv>=255) { break; } color = new
										 * Color(R+colv, G+colv, B+colv); }
										 */
									}
								}
							}
						}
						resultImg.setRGB(i, j, color.getRGB());
						writableRaster.setPixel(i, j, new float[] { watDepth });
					}
				}
			}

			GridCoverage2D tmp = RasterAnalyseTool.bufferedImageToGridCoverage2D(resultImg_gray, res.getEnvelope());
			// 写tif
			File tst = new File(storedAsFileName + ".tif");
			tst.createNewFile();
			GeoTiffWriter writer;
			writer = new GeoTiffWriter(tst);
			writer.write((GridCoverage) tmp, null);
			writer.dispose();

			// 写png
			File tst1 = new File(storedAsFileName + ".png");
			tst1.createNewFile();
			ImageIO.write(resultImg, "png", tst1);

			WaterDepthMetaInfo waterDepthMetaInfo = new WaterDepthMetaInfo();
			waterDepthMetaInfo.setNid(String.valueOf(System.nanoTime()));
			waterDepthMetaInfo.setWaterSurfaceID(waterSurface.getNid());
			waterDepthMetaInfo.setDamParamID(waterSurface.getDamParamID());
			waterDepthMetaInfo.setSectionParamID(waterSurface.getSectionParamID());
			// waterDepthMetaInfo.setMaxDepth(maxDepth);
			// waterDepthMetaInfo.setMinDepth(minDepth);
			waterDepthMetaInfo.setPath(storedAsFileName);
			waterDepthMetaInfo.setArea(waterSurface.getWkt().getArea());
			Envelope rect = waterSurface.getWkt().getEnvelopeInternal();
			waterDepthMetaInfo.setMinX(rect.getMinX());
			waterDepthMetaInfo.setMinY(rect.getMinY());
			waterDepthMetaInfo.setMaxX(rect.getMaxX());
			waterDepthMetaInfo.setMaxY(rect.getMaxY());
			return waterDepthMetaInfo;
		} catch (IOException e) {
			e.printStackTrace();
		} // catch (TransformException e) {
			// // TODO Auto-generated catch block
			// e.printStackTrace();
			// } catch (ParseException e) {
			// // TODO Auto-generated catch block
			// e.printStackTrace();
			// }
		return null;
	}

	private Color changeColor(String[] rgb0, String[] rgb1, float maxVal, float watDepth) {
		// String[] rgb = ((JSONObject)colorLevel.get(0)).getString("color").split(",");
		int colv = 1;
		int R0 = Integer.valueOf(rgb0[0]), G0 = Integer.valueOf(rgb0[1]), B0 = Integer.valueOf(rgb0[2]);
		int R1 = Integer.valueOf(rgb1[0]), G1 = Integer.valueOf(rgb1[1]), B1 = Integer.valueOf(rgb1[2]);

		Color color = new Color(R1, G1, B1);
		while (watDepth <= (maxVal -= 0.01f)) {
			color = new Color(R1 + colv < R0 ? R1 + colv : R0 - 1, G1 + colv < G0 ? G1 + colv : G0 - 1,
					B1 + colv <= B0 ? B1 + colv : B0 - 1);
			colv++;
		}
		return color;
	}

	private static JSONArray getColorLevel() {
		String color = "[{\"value\":5,\"color\":\"230,255,255\"},{\"value\":15,\"color\":\"110,255,255\"},{\"value\":30,\"color\":\"0,230,230\"},{\"value\":40,\"color\":\"0,191,191\"},{\"value\":50,\"color\":\"0,147,217\"},{\"value\":50,\"color\":\"0,80,170\"}]";
		JSONArray colorLevel = JSONArray.parseArray(color);
		return colorLevel;
		/*
		 * double[] level = new double[] { 5, 15, 30, 40, 50 }; // double[] level = new
		 * double[] {0.5,1,1.5,2.5,5.0}; Color color = null; if (watDepth <= 0) { color
		 * = new Color(0, 0, 0, 0); watDepth=0;
		 * 
		 * } else if (watDepth <= level[0]) { color = new Color(230, 255, 255); } else
		 * if (watDepth > level[0] && watDepth <= level[1]) { color = new Color(110,
		 * 255, 255); } else if (watDepth > level[1] && watDepth <= level[2]) { color =
		 * new Color(0, 230, 230); } else if (watDepth > level[2] && watDepth <=
		 * level[3]) { color = new Color(0, 191, 191); } else if (watDepth > level[3] &&
		 * watDepth <= level[4]) { color = new Color(0, 147, 217); } else if (watDepth >
		 * level[4]) { color = new Color(0, 80, 170); }
		 */
	}

	/**
	 * 在各断面的bufferDistance缓冲范围内，读取该坐标对应的水位。
	 * 
	 * @param cood
	 * @param bufferDistance
	 * @return
	 * @throws ParseException
	 */
	private double findMaxZOfSection2(Coordinate cood, float bufferDistance) throws ParseException {
		GeometryFactory geometryFactory = new GeometryFactory();
		com.vividsolutions.jts.geom.Point p = geometryFactory.createPoint(cood);
		int size = sectionList.size();
		int num = 0, bufferNum = 5;
		Section aimSec = null;
		for (int i = 1; i < size; i++) {
			if (num == bufferNum)
				break;
			Section pre = sectionList.get(i - 1);
			Point p1 = pre.getSectionLine().get(pre.getStartIndex());
			Point p2 = pre.getSectionLine().get(pre.getEndIndex());
			aimSec = sectionList.get(i);
			Point p3 = aimSec.getSectionLine().get(aimSec.getStartIndex());
			Point p4 = aimSec.getSectionLine().get(aimSec.getEndIndex());
			Polygon polygon = null;
			/*
			 * LineString l1 = geometryFactory.createLineString(new Coordinate[] {
			 * p1.getCoordinate(), p2.getCoordinate()}); LineString l2 =
			 * geometryFactory.createLineString(new Coordinate[] { p3.getCoordinate(),
			 * p4.getCoordinate()}); if(l1.intersects(l2)) { polygon =
			 * geometryFactory.createPolygon(new Coordinate[] { p1.getCoordinate(),
			 * p3.getCoordinate(), p2.getCoordinate(), p4.getCoordinate(),
			 * p1.getCoordinate() }); } else { polygon = geometryFactory.createPolygon(new
			 * Coordinate[] { p1.getCoordinate(), p3.getCoordinate(), p4.getCoordinate(),
			 * p2.getCoordinate(), p1.getCoordinate() }); }
			 */

			polygon = geometryFactory.createPolygon(new Coordinate[] { p1.getCoordinate(), p3.getCoordinate(),
					p4.getCoordinate(), p2.getCoordinate(), p1.getCoordinate() });
			// polygon = geometryFactory.createPolygon((LinearRing) polygon.getBoundary());
			if (!polygon.contains(p)) {
				aimSec = null;
				continue;
			}
			int stIndex = aimSec.getStartIndex();
			int endIndex = aimSec.getEndIndex();
			for (int j = stIndex; j <= endIndex; j++) {
				Point point = aimSec.getPointOfSectionLine(j);
				Coordinate coor = new Coordinate(point.getLon(), point.getLat());
				double d = coor.distance(cood);
				if (d < bufferDistance) {
					num++;
					if (num >= bufferNum) {
						break;
					}
				}
			}
			if (num < bufferNum) {
				aimSec = pre;
				num = bufferNum;
			}
		}
		if (null == aimSec)
			return -1;// 说明该点在水面范围之外，所以应该不算。
		// aimSec = sectionList.get(0);
		return aimSec.getMaxZ();
	}

	/**
	 * 以各断面水位为基础，用插值法，计算各坐标点的水位。
	 * 
	 * @param cood
	 * @return
	 * @throws ParseException
	 */
	private double findMaxZOfSection(Coordinate cood, SectionExactParam param) throws ParseException {
		GeometryFactory geometryFactory = new GeometryFactory();
		com.vividsolutions.jts.geom.Point p = geometryFactory.createPoint(cood);
		int size = sectionList.size();
		double aimZ = -1;
		Section aimSec = null;
		for (int i = 1; i < size; i++) {
			Section pre = sectionList.get(i - 1);
			Point p1 = pre.getSectionLine().get(pre.getStartIndex());
			Point p2 = pre.getSectionLine().get(pre.getEndIndex());
			aimSec = sectionList.get(i);
			Point p3 = aimSec.getSectionLine().get(aimSec.getStartIndex());
			Point p4 = aimSec.getSectionLine().get(aimSec.getEndIndex());
			Polygon polygon = null;
			LineSegment l1 = new LineSegment(p1.getCoordinate(), p2.getCoordinate());
			LineSegment l2 = new LineSegment(p3.getCoordinate(), p4.getCoordinate());
			polygon = geometryFactory.createPolygon(new Coordinate[] { p1.getCoordinate(), p3.getCoordinate(),
					p4.getCoordinate(), p2.getCoordinate(), p1.getCoordinate() });
			// double buffer = param.getStep();
			// Geometry geometry = polygon.buffer(buffer);
			Geometry geometry = polygon;
			if (!geometry.contains(p)) {
				aimSec = null;
				continue;
			}
			double distance1 = l1.distance(cood);
			double distance2 = l2.distance(cood);
			double distance3 = l2.distance(l1);
			aimZ = Interpolation.linear(pre.getDistanceOffStart(), pre.getMaxZ(), aimSec.getDistanceOffStart(),
					aimSec.getMaxZ(), pre.getDistanceOffStart() + distance1);
			// double zx = pre.getMaxZ() - aimSec.getMaxZ();
			// aimZ = pre.getMaxZ() - zx * (distance1 / (distance1 + distance2));
			break;
		}
		if (aimSec == null) {
			double[] start = param.getStartPoint();
			double[] end = param.getEndPoint();
			double z = getZFromOutter(sectionList.get(0), p, start);
			if (Double.compare(z, -1) == 0) {
				z = getZFromOutter(sectionList.get(sectionList.size() - 1), p, end);
				if (Double.compare(z, -1) > 0) {
					aimZ = z;
				}
			} else {
				aimZ = z;
			}
		}
		return aimZ;
	}

	private double getZFromOutter(Section section, com.vividsolutions.jts.geom.Point p, double[] point) {
		double result = -1;
		GeometryFactory geometryFactory = new GeometryFactory();
		Point p1 = section.getSectionLine().get(section.getStartIndex());
		Point p2 = section.getSectionLine().get(section.getEndIndex());
		LineSegment l1 = new LineSegment(p1.getCoordinate(), p2.getCoordinate());
		Polygon union = (Polygon) geometryFactory.createPolygon(new Coordinate[] { p1.getCoordinate(),
				p2.getCoordinate(), new Coordinate(point[0], point[1]), p1.getCoordinate() });
		if (union.contains(p)) {
			double distance = l1.distance(p.getCoordinate());
			result = section.getMaxZ();
		}
		return result;
	}

	/**
	 * 将断面集合中，每相邻两个断面组成一个多边形，形成一个多边形数组。在计算每个点的水位时，会用到。
	 * 
	 * @param cood
	 * @param bufferDistance
	 * @return
	 * @throws ParseException
	 */
	private List<Polygon> mkSectionZoneList(Coordinate cood, float bufferDistance) throws ParseException {
		List<Polygon> lp = new ArrayList<Polygon>();
		GeometryFactory geometryFactory = new GeometryFactory();
		com.vividsolutions.jts.geom.Point p = geometryFactory.createPoint(cood);
		int size = sectionList.size();
		Section aimSec = null;
		for (int i = 0; i < size; i++) {
			Section pre = sectionList.get(i - 1);
			Point p1 = pre.getSectionLine().get(pre.getStartIndex());
			Point p2 = pre.getSectionLine().get(pre.getEndIndex());
			aimSec = sectionList.get(i);
			Point p3 = aimSec.getSectionLine().get(aimSec.getStartIndex());
			Point p4 = aimSec.getSectionLine().get(aimSec.getEndIndex());
			Polygon polygon = geometryFactory.createPolygon(new Coordinate[] { p1.getCoordinate(), p3.getCoordinate(),
					p4.getCoordinate(), p2.getCoordinate(), p1.getCoordinate() });
			lp.add(polygon);
		}
		return lp;
	}

	public String getDemRaster() {
		return demRaster;
	}

	public String getReverwayShp() {
		return reverwayShp;
	}

	public void setDemRaster(String demRaster) {
		this.demRaster = demRaster;
	}

	public void setReverwayShp(String reverwayShp) {
		this.reverwayShp = reverwayShp;
	}

	@Override
	public boolean configBaseData(String demRaster, String riverWayShp, String boundShp) {
		boolean isGood = false;
		this.demRaster = demRaster;
		this.reverwayShp = riverWayShp;
		this.boundShp = boundShp;
		try {
			loadBaseData();
			isGood = true;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return isGood;
	}

	@Override
	public String getUUID() {
		return "88dbce03-38b9-41d4-b266-15bc232be0ca";
	}

	public CoordinateReferenceSystem getCrs() {
		return crs;
	}
}
package com.kygis.algorithmofdam;

import java.awt.Color;
import java.awt.color.ColorSpace;
import java.awt.image.BufferedImage;
import java.awt.image.ColorConvertOp;
import java.awt.image.Raster;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridEnvelope2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.coverage.processing.CoverageProcessor;
import org.geotools.data.FeatureSource;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureIterator;
import org.geotools.feature.FeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.gce.geotiff.GeoTiffReader;
import org.geotools.gce.geotiff.GeoTiffWriter;
import org.opengis.coverage.Coverage;
import org.opengis.coverage.grid.GridCoverage;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.parameter.ParameterValueGroup;
import org.opengis.referencing.crs.CoordinateReferenceSystem;
import org.opengis.referencing.operation.TransformException;

import com.kunyuan.framework.component.CmptContext;
import com.kunyuan.framework.component.CmptDescription;
import com.kunyuan.framework.component.CmptResult;
import com.kunyuan.framework.component.DefaultCmptContext;
import com.kunyuan.framework.component.DefaultCmptResult;
import com.kunyuan.framework.component.KYComponent;
import com.kygis.kuibamodel.po.DestroiedDamParam;
import com.kygis.kuibamodel.po.Point;
import com.kygis.kuibamodel.po.Section;
import com.kygis.kuibamodel.po.SectionExactParam;
import com.kygis.kuibamodel.po.WaterDepthMetaInfo;
import com.kygis.kuibamodel.po.WaterSurface;
import com.kygis.kuibamodel.po.ZQ;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.MultiPolygon;
import com.vividsolutions.jts.geom.Polygon;
import com.vividsolutions.jts.io.ParseException;

/**
 * 仅仅是lz用的测试文件，其他人可不用关注它。
 * 
 * @author wheeler
 *
 */
public class MyMain implements KYComponent {

	public static void main(String[] args) throws TransformException, ParseException, IOException {

		test();

		// MyMain.waterDepthAnalyseTest();
		// MyMain.waterDepthAnalyse();
		DefaultCmptContext context = new DefaultCmptContext();
		context.set("width", 150d);
		context.set("height", 55d);
		context.set("z", 400f);
		context.set("w", 50000000d);
		context.set("n", 0.035d);
		context.set("v", 2.5d);

		context.set("distance", 1);
		context.set("step", 1000);
		context.set("max", 370f);
		context.set("pc", 12000);
		context.set("ep", new double[] { 14042141, 5280013 });
		// context.set("ep", new double[] { 14047261, 5279700 });
		context.set("sp", new double[] { 14015867.2, 5292649.7 });
		// context.set("cacc", "0.6");
		context.set("demfile", "E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/panhai_5wan_3857_geoserver.tif");
		// context.set("demfile", "E:\\kunyuan\\map\\panhai_5wan_3857_geoserver.tif");
		context.set("riverfile", "E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/2亮子河/亮子河_3857_line.shp");
		context.set("bound", "E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/边界/riverbound.shp");

		context.set("output", "E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/test1/");
		MyMain main = new MyMain();
		main.execute(context);
		System.exit(0);

		/*
		 * boolean isGood = damAnalyseService.configBaseData(
		 * "/home/wheeler/workspace/mygood/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/4免费30米数据/jl_dem_30_3857_float_geoserver.tif",
		 * "/home/wheeler/workspace/mygood/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/2亮子河/亮子河_3857_line.shp",
		 * "/home/wheeler/workspace/mygood/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/边界/盘海水库坝址下游最大淹没范围_3857.shp"
		 * );
		 */
		/*
		 * boolean isGood = damAnalyseService.configBaseData(
		 * "/home/wheeler/workspace/mygood/gisdata/fh/tiantang_dem_3857_geoserver.tif",
		 * "/home/wheeler/workspace/mygood/gisdata/fh/河道线.shp",
		 * "/home/wheeler/workspace/mygood/gisdata/fh/边界.shp");
		 */

		/*
		 * tiantang sectionExactParam.setStartPoint(new
		 * double[]{12870221.641732594,3643432.1892489414});
		 * sectionExactParam.setEndPoint(new double[]
		 * {12840392.033327293,3629921.9181774827 });
		 */
		// 盘海
		// sectionExactParam.setStartPoint(new
		// double[]{14010458.854106868,5299577.678112299});
		// sectionExactParam.setEndPoint(new double[]
		// {14047264.132526495,5279716.239901049});

	}

	@Override
	public CmptDescription getDesc() {
		return null;
	}

	@Override
	public CmptResult execute(CmptContext context) {
		DefaultCmptResult result = new DefaultCmptResult();
		DefaultCmptContext ctx = (DefaultCmptContext) context;
		DamAnalyseServiceImpl damAnalyseService = new DamAnalyseServiceImpl();
		// 1、输入参数：溃坝参数、断面提取参数
		DestroiedDamParam destroiedDamParam = new DestroiedDamParam();
		destroiedDamParam.setWidth(ctx.get("width")); // 溃口宽
		destroiedDamParam.setHeight(ctx.get("height"));// 溃口深
		destroiedDamParam.setDamTopDem(ctx.get("z"));// 起溃时水位
		destroiedDamParam.setWaterVolume(ctx.get("w"));// 起溃时下泻水量
		destroiedDamParam.setRoughRate(ctx.get("n"));
		destroiedDamParam.setNid(ctx.get("cid"));
		if (ctx.get("cacc") != null) {
			destroiedDamParam.setCacc(Double.parseDouble(ctx.get("cacc")));
		}

		SectionExactParam sectionExactParam = new SectionExactParam();
		sectionExactParam.setDistance(ctx.get("distance"));// 断面上的点之间距离1米
		sectionExactParam.setStep(ctx.get("step"));// 断面间距1公里
		sectionExactParam.setTopDem(ctx.get("max"));// 断面最大高程
		sectionExactParam.setEndPoint(ctx.get("ep"));
		sectionExactParam.setStartPoint(ctx.get("sp"));
		sectionExactParam.setPointCount((int) ctx.get("pc"));

		boolean isGood = damAnalyseService.configBaseData(ctx.get("demfile"), ctx.get("riverfile"), ctx.get("bound"));

		if (!isGood)
			return null;// 如果数据配置和加载失败，无法继续执行。

		// 2、计算溃口最大流量
		double waterThroughput = Formula.broadCrestedWeir(destroiedDamParam.getHeight(), destroiedDamParam.getWidth());
		destroiedDamParam.setWaterThroughput(waterThroughput);
		System.out.println("溃口最大流量 " + waterThroughput);
		result.set("maxQ", waterThroughput);
		// 3、提取切断面
		List<Section> sectionList = damAnalyseService.sectionAnalyse(sectionExactParam);// 用dao将sectionList入库
		// List<Section> sectionList =
		// damAnalyseService.sectionAnalyseWidthBound(sectionExactParam, false);
		for (int i = 0; i < sectionList.size(); i++) {
			System.out.println(sectionList.get(i).getLineString());
		}

		// 4、计算非出口断面的最大流量
		for (int i = 0; i < sectionList.size(); i++) {
			Section section = sectionList.get(i);
			double wt = Formula.liszt(waterThroughput, section.getDistanceOffStart(),
					destroiedDamParam.getWaterVolume(), ctx.get("v"), 1);
			section.setMaxQ(wt);
			section.setRoughRate(destroiedDamParam.getRoughRate());
			System.out.println(String.format("断面%s的最大流量为%s", i, wt));
			// System.out.println(sectionList.get(i).getLineString());
		}

		// 5、计算出口断面水位流量关系
		Section laSection = sectionList.get(sectionList.size() - 1);
		List<ZQ> zqList = DestoriedDamAlgorithm.calcZQListByManning(sectionList, sectionList.size() - 1,
				waterThroughput);
		/*
		 * for (ZQ zq : zqList) { System.out.format("水位%s,流量%s", zq.getZ(), zq.getQ());
		 * System.out.println(); }
		 */
		laSection.setWaterVolumeMap(zqList);
		laSection.setMaxZ(DestoriedDamAlgorithm.calcZByQ(zqList, laSection.getMaxQ()));
		Object[] sResult = DestoriedDamAlgorithm.getSectionWidthByZ(laSection, laSection.getMaxZ(),0);
		laSection.setStartIndex(((Point) sResult[1]).getPosiInSection());
		laSection.setEndIndex(((Point) sResult[2]).getPosiInSection());
		DestoriedDamAlgorithm.cfg.put("cacc", destroiedDamParam.getCacc());
		// 6、计算各断面水位
		long st1 = System.currentTimeMillis();
		double[] waterDemArray = DestoriedDamAlgorithm.calcZByWaterEnergyBalance(sectionList);
		System.out.println("计算各断面水位耗时ms： " + (System.currentTimeMillis() - st1));
		// 7、形成水面线
		WaterSurface waterSurface = damAnalyseService.waterSurfaceAnalyse(waterDemArray, sectionList);
		// Polygon boundary =
		// damAnalyseService.waterSurfaceBoundaryAnalyse(waterDemArray, sectionList);
		System.out.println(waterSurface.getWkt().toText());
		// 8、生成水深栅格图
		long st = System.currentTimeMillis();
		System.out.println(" begin render ");
		// String dataBaseDir = "/home/wheeler/workspace/mygood/gisdata/fh/";
		String dataBaseDir = ctx.get("output");
		// String fileName = waterSurface.getNid() + "_" + destroiedDamParam.getNid();
		String storedAsFileName = dataBaseDir + System.nanoTime();
		damAnalyseService.waterDepthAnalyse(waterSurface, sectionExactParam.getStep() / 2, storedAsFileName);
		// damAnalyseService.waterDepthAnalyse(waterSurface, storedAsFileName);
		System.out.println(" end render 耗时ms： " + (System.currentTimeMillis() - st));
		/*
		 * List<Map<String, Object>> sList = new ArrayList<>(); for (Section s :
		 * sectionList) { List<Point> points = s.getSectionLine(); Map<String, Object>
		 * map = new HashMap<>(); sList.add(map); map.put("distance",
		 * s.getDistanceOffStart()); map.put("x", s.getRiverPoint().getLon());
		 * map.put("y", s.getRiverPoint().getLat()); map.put("q", s.getMaxQ());
		 * map.put("z", s.getMaxZ()); map.put("n", s.getRoughRate()); map.put("si",
		 * s.getStartIndex()); map.put("ei", s.getEndIndex()); List<Map<String, Double>>
		 * ps = new ArrayList<>(points.size()); map.put("points", ps); for (Point p :
		 * points) { Map<String, Double> pm = new HashMap<>(); pm.put("x", p.getLon());
		 * pm.put("y", p.getLat()); pm.put("z", (double) p.getDem()); ps.add(pm); } }
		 * result.set("imgOut", fileName + ".png"); result.set("demOut", fileName +
		 * ".tif"); result.set("areaOut", fileName + ".shp"); result.set("area",
		 * metaInfo.getArea()); result.set("minx", metaInfo.getMinX());
		 * result.set("miny", metaInfo.getMinY()); result.set("maxx",
		 * metaInfo.getMaxX()); result.set("maxy", metaInfo.getMaxY());
		 * result.set("section", sList);
		 */
		return result;
	}

	public static Coverage waterDepthAnalyseTest() throws TransformException, ParseException, IOException {

		File demFile = new File(
				"E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/4免费30米数据/jl_dem_30_3857_float_geoserver.tif");
		GeoTiffReader tifReader = new GeoTiffReader(demFile);
		GridCoverage2D coverage = tifReader.read(null);
		CoordinateReferenceSystem crs = coverage.getCoordinateReferenceSystem2D();

		// 加载边界

		File riverWayFile = new File("E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/边界/bound.shp");
		ShapefileDataStore shpDataStore = null;
		shpDataStore = new ShapefileDataStore(riverWayFile.toURI().toURL());
		String typeName = shpDataStore.getTypeNames()[0];
		FeatureSource<SimpleFeatureType, SimpleFeature> featureSource = null;
		featureSource = (FeatureSource<SimpleFeatureType, SimpleFeature>) shpDataStore.getFeatureSource(typeName);
		FeatureCollection<SimpleFeatureType, SimpleFeature> result = featureSource.getFeatures();

		FeatureIterator<SimpleFeature> itertor = result.features();
		MultiPolygon tmpGeoObj1 = (MultiPolygon) itertor.next().getDefaultGeometry();
		Polygon boundObj = (Polygon) tmpGeoObj1.getGeometryN(0);

		CoverageProcessor processor = new CoverageProcessor();
		ParameterValueGroup params = processor.getOperation("CoverageCrop").getParameters();
		params.parameter("Source").setValue(coverage);
		// params.parameter("ENVELOPE").setValue(bounds);
		params.parameter("ROI").setValue(boundObj);
		params.parameter("ForceMosaic").setValue(true);

		GridCoverage2D res = (GridCoverage2D) processor.doOperation(params);

		write(res, "test1_" + System.currentTimeMillis());

		GridGeometry2D geometry = res.getGridGeometry();

		RenderedImage img = res.getRenderedImage();
		Raster raster = img.getData();
		int width = raster.getWidth();
		int height = raster.getHeight();

		BufferedImage resultImg = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		int minX = raster.getMinX();
		int minY = raster.getMinY();
		float maxDepth = 0, minDepth = 0;
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				float x = ((float[]) raster.getDataElements(i + minX, j + minY, (float[]) null))[0];
				if (x == 0) {
					Color color = new Color(0, 0, 0, 0);
					resultImg.setRGB(i, j, color.getRGB());
				} else {

					org.geotools.geometry.Envelope2D pixelEnvelop = geometry
							.gridToWorld(new GridEnvelope2D(i + minX, j + minY, 1, 1));

					double lat = pixelEnvelop.getCenterY();
					double lon = pixelEnvelop.getCenterX();

					double maxZ = 380; // 应该根据lonlat，从sectionList里找到最靠近section，得到其水位。

					float watDepth = (float) (maxZ - x);

					if (watDepth < minDepth) {
						minDepth = watDepth;
					} else if (watDepth > maxDepth) {
						maxDepth = watDepth;
					}

					Color color = null;
					if (watDepth <= 0) {
						color = new Color(0, 0, 0, 0);
					} else if (watDepth <= 0.5) {
						color = new Color(230, 255, 255);
					} else if (watDepth > 0.5 && watDepth <= 1.0) {
						color = new Color(110, 255, 255);
					} else if (watDepth > 1.0 && watDepth <= 1.5) {
						color = new Color(0, 230, 230);
					} else if (watDepth > 1.5 && watDepth <= 2.5) {
						color = new Color(0, 191, 191);
					} else if (watDepth > 2.5 && watDepth <= 5.0) {
						color = new Color(0, 147, 217);
					} else if (watDepth > 5.0) {
						color = new Color(0, 80, 170);
					}
					resultImg.setRGB(i, j, color.getRGB());
				}
			}
		}
		GridCoverage2D tmp = RasterAnalyseTool.bufferedImageToGridCoverage2D(resultImg, res.getEnvelope());
		GeoTiffWriter writer;
		try {
			File tst = new File(
					"E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/test/" + System.currentTimeMillis() + ".tif");
			tst.createNewFile();
			writer = new GeoTiffWriter(tst);
			writer.write((GridCoverage) tmp, null);
			writer.dispose();

			File tst1 = new File(
					"E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/test/" + System.currentTimeMillis() + ".png");
			tst1.createNewFile();
			ImageIO.write(resultImg, "png", tst1);

		} catch (IOException e) {
			e.printStackTrace();
		}
		return tmp;
	}

	public static void write(GridCoverage tmp, String name) {
		GeoTiffWriter writer;
		try {
			File tst = new File("E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/test/" + name + ".tif");
			tst.createNewFile();
			writer = new GeoTiffWriter(tst);
			writer.write(tmp, null);
			writer.dispose();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static void test() throws IOException {
		File demFile = new File("E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/test3/24558133591955_null.tif");
		GeoTiffReader tifReader = new GeoTiffReader(demFile);
		GridCoverage2D coverage = tifReader.read(null);
		GridGeometry2D geometry = coverage.getGridGeometry();

		RenderedImage img = coverage.getRenderedImage();
		Raster raster = img.getData();
		int width = raster.getWidth();
		int height = raster.getHeight();
		int minX = raster.getMinX();
		int minY = raster.getMinY();

		List<Coordinate> boundaryPoint = new ArrayList<Coordinate>();
		Coordinate preCoor;
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
				if (x <= 0) {
					preCoor = new Coordinate(lon, lat);
				} else {
					boundaryPoint.add(new Coordinate(lon, lat));
				}
			}
		}
	}

	public static Coverage waterDepthAnalyse() throws IOException {

		File demFile = new File(
				"E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/4免费30米数据/jl_dem_30_3857_float_geoserver.tif");
		GeoTiffReader tifReader = new GeoTiffReader(demFile);
		GridCoverage2D coverage = tifReader.read(null);

		// 加载边界
		File riverWayFile = new File("E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/1盘海水库位置/盘海水库位置_3857.shp");
		ShapefileDataStore shpDataStore = null;
		shpDataStore = new ShapefileDataStore(riverWayFile.toURI().toURL());
		String typeName = shpDataStore.getTypeNames()[0];
		FeatureSource<SimpleFeatureType, SimpleFeature> featureSource = null;
		featureSource = (FeatureSource<SimpleFeatureType, SimpleFeature>) shpDataStore.getFeatureSource(typeName);
		FeatureCollection<SimpleFeatureType, SimpleFeature> result = featureSource.getFeatures();

		FeatureIterator<SimpleFeature> itertor = result.features();
		MultiPolygon tmpGeoObj1 = (MultiPolygon) itertor.next().getDefaultGeometry();
		Polygon boundObj = (Polygon) tmpGeoObj1.getGeometryN(0);

		CoverageProcessor processor = new CoverageProcessor();
		ParameterValueGroup params = processor.getOperation("CoverageCrop").getParameters();
		params.parameter("Source").setValue(coverage);
		params.parameter("ROI").setValue(boundObj);
		params.parameter("ForceMosaic").setValue(true);
		GridCoverage2D res = (GridCoverage2D) processor.doOperation(params);

		RenderedImage redImg = res.getRenderedImage();
		BufferedImage img = RasterAnalyseTool.toBufferedImage(redImg);
		Raster raster = img.getData();
		for (int i = 0; i < img.getWidth(); i++) {
			for (int j = 0; j < img.getHeight(); j++) {
				// float x = raster.getSample(i, j, 1); raster.getPixel(i,j,(double[])null);
				float v = ((float[]) raster.getDataElements(3, 120, (float[]) null))[0];
				Color c = new Color(img.getRGB(i, j));
				if (c.equals(Color.BLACK)) {
					Color color = new Color(0, 0, 0, 0);
					img.setRGB(i, j, color.getRGB());
				}
			}
		}
		GridCoverage2D tmp = RasterAnalyseTool.bufferedImageToGridCoverage2D(img, res.getEnvelope());

		GeoTiffWriter writer;
		try {
			File tst = new File(
					"E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/test/" + System.currentTimeMillis() + ".tif");
			tst.createNewFile();
			writer = new GeoTiffWriter(tst);
			writer.write((GridCoverage) tmp, null);
			writer.dispose();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return tmp;
	}

}

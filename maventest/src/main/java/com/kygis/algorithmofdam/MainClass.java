package com.kygis.algorithmofdam;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.geotools.data.DataUtilities;
import org.geotools.feature.SchemaException;
import org.opengis.feature.simple.SimpleFeatureType;

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
import com.kygis.util.ShapeTools;

public class MainClass implements KYComponent {

	public static void main(String[] args) {
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
		context.set("sp", new double[] { 14015867.2, 5292649.7 });
		context.set("ep", new double[] { 14042141, 5280013 });
		context.set("demfile", "E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/panhai_5wan_3857_geoserver.tif");
		context.set("riverfile", "E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/2亮子河/亮子河_3857_line.shp");
		context.set("bound", "E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/边界/riverbound.shp");
		context.set("output", "E:/gisdata/data_3857_盘海水库_V20180527_FLOAT_GEOSERVE/test3/");
		MainClass main = new MainClass();
		main.execute(context);
		System.exit(0);
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
		Logger.info("溃口最大流量 " + waterThroughput);
		result.set("maxQ", waterThroughput);
		// 3、提取切断面
		// List<Section> sectionList
		// =damAnalyseService.sectionAnalyse(sectionExactParam);// 用dao将sectionList入库
		// List<Section> sectionList =
		// damAnalyseService.sectionAnalyseWidthBound(sectionExactParam);
		// List<Section> sectionList =
		// damAnalyseService.sectionParallelAnalyseWidthBound(sectionExactParam);//
		// 平行断面顺着河道版本
		// List<Section> sectionList =
		// damAnalyseService.sectionParallelAnalyseWidthBoundX(sectionExactParam); //
		// 平行断面全覆盖边界版本
		List<Section> sectionList = damAnalyseService.sectionParallelAnalyseWidthBoundXX(sectionExactParam);

		for (int i = 0; i < sectionList.size(); i++) {
			Logger.trace(sectionList.get(i).getLineString().toText());
		}

		// for (int i = 0; i < sectionList.size(); i++) {
		// List<Point> points = sectionList.get(i).getSectionLine();
		// for (int j = 0; j < points.size(); j++) {
		// Logger.debug(String.format("%s,%s,%s,%s,%s,%S", i,
		// sectionList.get(i).getMaxQ(), j,
		// points.get(j).getLon(), points.get(j).getLat(), points.get(j).getDem()));
		// }
		// }

		// 4、计算非出口断面的最大流量
		for (int i = 0; i < sectionList.size(); i++) {
			Section section = sectionList.get(i);
			double wt = Formula.liszt(waterThroughput, section.getDistanceOffStart(),
					destroiedDamParam.getWaterVolume(), ctx.get("v"), 1);
			section.setMaxQ(wt);
			section.setRoughRate(destroiedDamParam.getRoughRate());
			Logger.info(String.format("断面%s的最大流量为%s", i, wt));
		}
		// sectionList.get(0).setMaxQ(39030);
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

		Object[] sResult = DestoriedDamAlgorithm.getSectionWidthByZ(laSection, laSection.getMaxZ(), 0);

		laSection.setStartIndex(((Point) sResult[1]).getPosiInSection());
		laSection.setEndIndex(((Point) sResult[2]).getPosiInSection());
		DestoriedDamAlgorithm.cfg.put("cacc", destroiedDamParam.getCacc());
		// 6、计算各断面水位
		double[] waterDemArray = DestoriedDamAlgorithm.calcZByWaterEnergyBalance(sectionList);

		// int count = sectionList.size();
		// int size = sectionList.size() - 1;
		// for (int i = size; i >= 0; i--) {
		// Section sec = sectionList.get(i);
		// sec.restore();
		// Logger.trace(sec.getLineString().toText());
		// // System.out.println(String.format("%s,%s,1,%s,%s,0.035,0.035,0.035", count
		// -
		// // i, sec.getSectionLine().size(),
		// // sec.getSectionLine().size(), (count - i - 1) * 1000));
		// // for (int j = 0; j < sec.getSectionLine().size(); j++) {
		// // // Logger.debug(String.format("%s,%s", j,
		// // // sec.getSectionLine().get(j).getDem()));
		// // System.out.println(String.format("%s,%s", j,
		// // sec.getSectionLine().get(j).getDem()));
		// // }
		// }
		// 7、形成水面线
		WaterSurface waterSurface = damAnalyseService.waterSurfaceAnalyse(waterDemArray, sectionList);

		String dataBaseDir = ctx.get("output");
		String fileName = ctx.get("Ennmcd") + "_" + waterSurface.getNid() + "_" + destroiedDamParam.getNid();
		String storedAsFileName = dataBaseDir + fileName;
		new Thread(new Runnable() {
			@Override
			public void run() {
				String filepath = storedAsFileName;
				SimpleFeatureType type;
				try {

					type = DataUtilities.createType("watersurface", "the_geom:Polygon," + // <- the geometry attribute:
																							// Polygon type
					"NID:String," + // <- a String attribute
					"damParaID:String," + // a number attribute
					"sectionID:String," + "reserName:String");
					// 直接连点
					String f1 = filepath + "_original.shp";
					ShapeTools.write(type,
							new Object[] { waterSurface.getWktOriginal(), waterSurface.getNid(),
									waterSurface.getDamParamID(), waterSurface.getSectionParamID(), "panhaishuiku" },
							damAnalyseService.getCrs(), f1);
					// 取多边形
					String f2 = filepath + "_convex.shp";
					ShapeTools.write(type,
							new Object[] { waterSurface.getWktConVex(), waterSurface.getNid(),
									waterSurface.getDamParamID(), waterSurface.getSectionParamID(), "panhaishuiku" },
							damAnalyseService.getCrs(), f2);
					// 排除交叉点后的最小图形
					String f3 = filepath + "_min.shp";
					ShapeTools.write(
							type, new Object[] { waterSurface.getWkt(), waterSurface.getNid(),
									waterSurface.getDamParamID(), waterSurface.getSectionParamID(), "panhaishuiku" },
							damAnalyseService.getCrs(), f3);

				} catch (SchemaException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}).start();

		// 8、生成水深栅格图
		Logger.info(" begin render ");
		// String dataBaseDir = "/home/wheeler/workspace/mygood/gisdata/fh/";
		// String dataBaseDir = ctx.get("output");
		// String fileName = ctx.get("Ennmcd") + "_" + waterSurface.getNid() + "_" +
		// destroiedDamParam.getNid();

		new File(dataBaseDir).mkdirs();
		WaterDepthMetaInfo metaInfo = damAnalyseService.waterDepthAnalyse(sectionExactParam, waterSurface,
				storedAsFileName, null);
		// damAnalyseService.waterDepthAnalyse(waterSurface, storedAsFileName);
		Logger.info(" end render " + storedAsFileName);

		List<Map<String, Object>> sList = new ArrayList<>();
		for (Section s : sectionList) {
			List<Point> points = s.getSectionLine();
			Map<String, Object> map = new HashMap<>();
			sList.add(map);
			map.put("distance", s.getDistanceOffStart());
			map.put("x", s.getMinAltitudePoint().getLon());
			map.put("y", s.getMinAltitudePoint().getLat());
			map.put("q", s.getMaxQ());
			map.put("z", s.getMaxZ());
			map.put("n", s.getRoughRate());
			map.put("si", s.getStartIndex());
			map.put("ei", s.getEndIndex());
			List<Map<String, Double>> ps = new ArrayList<>(points.size());
			map.put("points", ps);
			for (Point p : points) {
				Map<String, Double> pm = new HashMap<>();
				pm.put("x", p.getLon());
				pm.put("y", p.getLat());
				pm.put("z", (double) p.getDem());
				ps.add(pm);
			}
		}
		result.set("imgOut", fileName + ".png");
		result.set("demOut", fileName + ".tif");
		result.set("areaOut", fileName + ".shp");
		result.set("area", metaInfo.getArea());
		result.set("minx", metaInfo.getMinX());
		result.set("miny", metaInfo.getMinY());
		result.set("maxx", metaInfo.getMaxX());
		result.set("maxy", metaInfo.getMaxY());
		result.set("section", sList);
		return result;
	}
}
package com.kygis.algorithmofdam;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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

public class LongBoMain implements KYComponent {

	public static void main(String[] args) {
		DefaultCmptContext context = new DefaultCmptContext();
		context.set("width", 150d);
		context.set("height", 20d);
		context.set("z", 550f);
		context.set("w", 40000000d);
		context.set("n", 0.035d);
		context.set("v", 2.5d);

		context.set("distance", 1);
		context.set("step", 1000);
		context.set("max", 550f);
		context.set("pc", 13000);
		context.set("sp", new double[] { 12428239, 4111275 });
		context.set("ep", new double[] { 12444042, 4089768 });
		context.set("demfile",
				"/home/wheeler/workspace/mygood/gisdata/longbo/ASTGTM2_N34E111_dem_3857_float_geoserver.tif");
		context.set("riverfile", "/home/wheeler/workspace/mygood/gisdata/longbo/river.shp");
		context.set("bound", "/home/wheeler/workspace/mygood/gisdata/longbo/lb_fw_line.shp");
		context.set("output", "/home/wheeler/workspace/mygood/gisdata/longbo/test/");
		LongBoMain main = new LongBoMain();
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
		List<Section> sectionList = damAnalyseService.sectionParallelAnalyseWidthBoundXX(sectionExactParam); // 平行断面全覆盖边界版本

		/*
		 * for (int i = 0; i < sectionList.size(); i++) {
		 * Logger.trace(sectionList.get(i).getLineString().toText()); }
		 */

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

		for (int i = 0; i < sectionList.size(); i++) {
			Section sec = sectionList.get(i);
			sec.restore();
			System.out.println(sec.getLineString());
		}

		// 7、形成水面线
		WaterSurface waterSurface = damAnalyseService.waterSurfaceAnalyse(waterDemArray, sectionList);

		// 8、生成水深栅格图
		Logger.info(" begin render ");
		// String dataBaseDir = "/home/wheeler/workspace/mygood/gisdata/fh/";
		String dataBaseDir = ctx.get("output");
		String fileName = waterSurface.getNid() + "_" + destroiedDamParam.getNid();
		String storedAsFileName = dataBaseDir + fileName;
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

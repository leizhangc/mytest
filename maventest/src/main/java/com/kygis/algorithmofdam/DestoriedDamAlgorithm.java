package com.kygis.algorithmofdam;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.kygis.kuibamodel.po.Point;
import com.kygis.kuibamodel.po.Section;
import com.kygis.kuibamodel.po.WaterSectionFeature;
import com.kygis.kuibamodel.po.ZQ;
import com.vividsolutions.jts.algorithm.distance.PointPairDistance;
import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.io.ParseException;

/**
 * 溃坝分析计算公式类
 * 
 * @author KAIFA01
 *
 */
final class DestoriedDamAlgorithm {

	public static Map<String, Object> cfg = new HashMap<>();
	private static int lastStartIndex;
	private static int lastEndIndex;

	/**
	 * 过水断面特征值计算
	 * 
	 * @param section
	 *            断面
	 * @param waterDem
	 *            水位
	 * @return
	 */
	public static WaterSectionFeature waterSectionAnalyse(Section section, double waterDem) {
		Coordinate[] corrdis = new Coordinate[2];
		double d = waterDem;
		List<Point> ps = section.getSectionLine();
		int size = ps.size();
		boolean left = true, right = true;
		for (int j = 1; j < size; j++) {
			if (left || right) {
				Point p0 = ps.get(j - 1);
				Point p1 = ps.get(j);
				if (d > p1.getDem() && left) {
					float dd = p0.getDem() - p1.getDem();
					corrdis[0] = ShpAnalyseTool.coordOfLinesegment(p1.getCoordinate(), p0.getCoordinate(), dd);
					left = false;
				}

				Point pp0 = ps.get(size - j);
				Point pp1 = ps.get(size - j - 1);
				if (d > pp1.getDem() && right) {
					float dd = pp0.getDem() - pp1.getDem();
					corrdis[1] = ShpAnalyseTool.coordOfLinesegment(pp1.getCoordinate(), pp0.getCoordinate(), dd);
					right = false;
				}
			} else {
				break;
			}
		}

		WaterSectionFeature ws = new WaterSectionFeature();
		ws.setWidth(calcDistance(corrdis[0], corrdis[1]));
		// TODO 有点重复与getSectionWidthByZ
		return null;
	}

	/**
	 * 过水断面流速
	 * 
	 * @param waterSectionArea
	 *            过水面积
	 * @param waterVolume
	 *            流量
	 * @return
	 */
	public static double waterSpeedAnalyse(double waterSectionArea, double waterVolume) {
		return 0;
	}

	/**
	 * 非出口断面比降计算
	 * 
	 * @param sectionList
	 *            全部断面
	 * @return 非出口断面比降数组
	 */
	public static float[] otherSectionGradientAnalyse(List<Section> sectionList) {
		return null;
	}

	/**
	 * 根据水位流量关系计算水位
	 * 
	 * @param list
	 *            水位流量列表
	 * @param q
	 *            流量
	 * @return 水位
	 */
	public static double calcZByQ(List<ZQ> list, double q) {
		double result = 0;
		double minZ, minQ, maxZ, maxQ;
		minZ = minQ = maxZ = maxQ = 0;

		for (int i = 0; i < list.size() - 1; i++) {
			ZQ current = list.get(i);
			ZQ next = list.get(i + 1);
			if (Double.compare(current.getQ(), q) < 0 && Double.compare(next.getQ(), q) > 0) {
				minZ = current.getZ();
				minQ = current.getQ();
				maxZ = next.getZ();
				maxQ = next.getQ();
				break;
			} else if (Double.compare(current.getQ(), q) == 0) {
				result = current.getZ();
				break;
			} else if (Double.compare(next.getQ(), q) == 0) {
				result = next.getZ();
				break;
			}
		}
		if (Double.compare(0d, result) == 0 && Double.compare(0d, minZ) < 0 && Double.compare(0d, maxZ) < 0) {
			result = Interpolation.linear(minQ, minZ, maxQ, maxZ, q);
		}
		if (Double.compare(0d, result) == 0) {
			Logger.warn("cannot get z by q {} ,min {} max {}", q, list.get(0).getQ(), list.get(list.size() - 1).getQ()); // TODO
			if (Double.compare(q, list.get(0).getQ()) < 0) { // 补救一下
				maxZ = list.get(0).getZ();
				maxQ = list.get(0).getQ();
				result = Interpolation.linear(0, 0, maxQ, maxZ, q);
			}
		}
		return result;
	}

	/**
	 * 根据曼宁公式计算某一断面水位流量关系
	 * 
	 * @param sections
	 * @param index
	 * @param qMax
	 * @param n
	 * @return
	 */
	public static List<ZQ> calcZQListByManning(List<Section> sections, int index, double qMax) {
		List<ZQ> result = new ArrayList<>();
		Section section = sections.get(index);
		double slope = calcAvgSlopeByZ(sections, index, section.getRoughRate());// 计算比降
		section.setGradient(slope);
		double z = 0.01d + section.getMinAltitudePoint().getDem();
		double q = 0;
		int testCount = 0;
		do {
			Object[] widthResult = getSectionWidthByZ(section, z, 0);
			double width = (double) widthResult[0];
			Point lPoint = (Point) widthResult[1];
			Point rPoint = (Point) widthResult[2];
			Object[] areaResult = calcSectionArea(section, lPoint, rPoint, z);
			double area = (double) areaResult[0];
			double h = area / width;
			q = Formula.manning(section.getRoughRate(), h, slope, width);
			result.add(new ZQ((double) widthResult[3], q, z));
			z += 0.05; // TODO
			testCount++;
			if (testCount > 6000) {
				Logger.warn(String.format("Cannot calc ZQ List By Manning, specify q %s ,current %s", qMax, q));
				break;
			}

		} while (q <= qMax);

		return result;
	}

	/**
	 * 根据河底高程计算平均比降
	 * 
	 * @param sections
	 *            断面序列
	 * @param index
	 *            需要计算比降的断面索引
	 * @return
	 */
	@SuppressWarnings("unused")
	private static double calcAvgSlopeByDEM(List<Section> sections, int index) {
		Section section = sections.get(index);
		float dem = section.getMinAltitudePoint().getDem();
		double slope = 0;
		int count = 6; // TODO
		int start = index - count;
		if (start < 0) {
			start = 0;
		}
		if (index == 0) {
			index = 1;
		}
		for (int i = start; i < index; i++) {
			Section uSection = sections.get(i);
			float uDem = uSection.getRiverPoint().getDem();
			slope += Formula.riverSlope(uDem, dem, section.getDistanceOffStart() - uSection.getDistanceOffStart());
		}
		slope /= count;
		if (slope <= 0 || Double.isNaN(slope) || Double.isInfinite(slope)) {
			Logger.warn("slope less than zero " + slope);
			slope = calcAvgSlopeByDEM(sections, index + 1); // TODO 需保证最后一个断面比降大于0
		}
		return slope;
	}

	/**
	 * 根据推求的水位计算比降
	 * 
	 * @param sections
	 * @param index
	 * @param n
	 * @return
	 */
	private static double calcAvgSlopeByZ(List<Section> sections, int index, double n) {
		double v = 2d; // TODO
		Section section = sections.get(index);
		// float dem = section.getRiverPoint().getDem();
		int count = 7; // TODO
		int start = index - count;
		if (start < 0) {
			start = 0;
		}
		if (index < count) {
			index = count;
		}
		double q = sections.get(index).getMaxQ();
		double area = q / v;
		double slope = 0d;
		double h = calcZByArea(section, area, n);
		// double z = h+dem;
		double z = h;
		for (int i = start; i < index; i++) {
			Section uSection = sections.get(i);
			// float uDem = uSection.getRiverPoint().getDem();
			// double uZ = calcZByArea(uSection, area, n) + uDem;
			double uZ = calcZByArea(uSection, area, n);
			double r = Formula.riverSlope(uZ, z, section.getDistanceOffStart() - uSection.getDistanceOffStart());
			if (r < 0 || Double.isNaN(r) || Double.isInfinite(r)) {
				count--;
			} else {
				slope += r;
			}
		}
		slope /= count;
		if (Double.compare(0d, slope) == 0 || Double.isNaN(slope)) {
			Logger.warn("slope is zero");
			slope = 0.00001;
		}
		return slope;
	}

	/**
	 * 根据面积求水深
	 * 
	 * @param section
	 * @param area
	 * @param n
	 * @return
	 */
	private static double calcZByArea(Section section, double area, double n) {
		double presion = 0.1d; // TODO
		double offset = 0.0001d; // TODO
		double result = 0d;
		double min = 0.01d;
		double max = 1000d;
		double h = min;
		double calcArea = 0d;
		int testCount = 0;
		do {
			if (Double.compare(area, calcArea) > 0) {
				min = h;
			} else {
				max = h;
			}
			h = min + (max - min) / 2;
			Object[] widthResult = getSectionWidthByZ(section, h, 0);
			Point leftPoint = (Point) widthResult[1];
			Point rightPoint = (Point) widthResult[2];

			Object[] areaResult = calcSectionArea(section, leftPoint, rightPoint, h);
			calcArea = (double) areaResult[0];
			if (Double.compare(max - min, offset) <= 0) {
				testCount++;
			}
			if (testCount > 20) { // min,max非常接近的情况下，尝试20次,不行就退出
				break;
			}
		} while (Double.compare(Math.abs(area - calcArea), presion) > 0);
		result = h;
		return result;
	}

	/**
	 * 根据水深计算断面宽度
	 * 
	 * @param section
	 *            断面
	 * @param z
	 *            水位
	 * @param type
	 *            截取类型
	 * @return [0]断面宽度<br>
	 *         [1]左岸点<br>
	 *         [2]右岸点
	 */
	public static Object[] getSectionWidthByZ(Section section, double z, int type) {
		if (type == 0) {
			return getSectionWidthByZFromSide(section, z);
		}
		return getSectionWidthByZFromCenter(section, z);
	}

	public static Object[] getSectionWidthByZFromCenter(Section section, double z) {
		/*
		 * double precision = 0.1; // TODO Object[] result = new Object[3];
		 * 
		 * List<Point> points = section.getSectionLine(); Point point =
		 * getMinAltitudePoint(points); double maxAltitude = point.getDem() + z;
		 * 
		 * int count = points.size() / 2; Point lPoint = null; Point rPoint = null; for
		 * (int i = 0; i < count; i++) { Point lp = points.get(i); Point rp =
		 * points.get(points.size()-1 - i); if (Double.compare(lp.getDem() - precision,
		 * maxAltitude) <= 0 && Double.compare(lp.getDem() + precision, maxAltitude) >=
		 * 0) { lPoint = lp; } if (Double.compare(rp.getDem() - precision, maxAltitude)
		 * <= 0 && Double.compare(rp.getDem() + precision, maxAltitude) >= 0) { rPoint =
		 * rp; } } result[0] = calcDistance(lPoint, rPoint); result[1] = lPoint;
		 * result[2] = rPoint; return result;
		 */
		// TODO 应做交线
		Coordinate[] corrdis = new Coordinate[2];
		Object[] result = new Object[4];

		List<Point> ps = section.getSectionLine();
		Point point = section.getRiverPoint();
		double precision = 1.0d;
		double d = z + precision;

		Point lPoint = null, rPoint = null;
		int size = ps.size();
		boolean left, right;
		left = right = false;
		int pos = point.getPosiInSection();
		for (int i = 1; i < size; i++) {
			Point lp = null;
			if (pos - i >= 0) {
				lp = ps.get(pos - i);
			} else {
				left = true;
			}
			Point rp = null;
			if (pos + i < ps.size()) {
				rp = ps.get(pos + i);
			} else {
				right = true;
			}
			if (!left && lp.getDem() > d) {
				Point p = ps.get(pos - i + 1);
				float dd = (float) ((d - p.getDem()) / (lp.getDem() - p.getDem()));
				if (Float.isInfinite(dd) || Float.isNaN(dd)) {
					dd = 0f;
				}
				left = true;
				corrdis[0] = ShpAnalyseTool.coordOfLinesegment(p.getCoordinate(), lp.getCoordinate(), dd);
				lPoint = new Point(corrdis[0].x, corrdis[0].y, (float) d);

				lPoint.setPosiInSection(pos - i);
			}
			if (!right && rp.getDem() > d) {
				Point p = ps.get(pos + i - 1);
				float dd = (float) ((d - p.getDem()) / (rp.getDem() - p.getDem()));
				if (Float.isInfinite(dd) || Float.isNaN(dd)) {
					dd = 0f;
				}
				right = true;
				corrdis[0] = ShpAnalyseTool.coordOfLinesegment(p.getCoordinate(), rp.getCoordinate(), dd);
				rPoint = new Point(corrdis[0].x, corrdis[0].y, (float) d);
				rPoint.setPosiInSection(pos + i);

			}
			if (left && right) {
				break;
			}
		}

		// for (int j = 1; j < size / 2; j++) {
		// Point p0 = ps.get(j - 1);
		// Point p1 = ps.get(j);
		// int a = 0;
		// if (d > p1.getDem()) {
		// a++;
		// }
		// if (d > p1.getDem() && d < p0.getDem()) {
		// float dd = (float) ((d - p1.getDem()) / (p0.getDem() - p1.getDem()));
		// if (Float.isInfinite(dd) || Float.isNaN(dd)) {
		// dd = 0f;
		// }
		// corrdis[0] = ShpAnalyseTool.coordOfLinesegment(p1.getCoordinate(),
		// p0.getCoordinate(), dd);
		// lPoint = new Point(corrdis[0].x, corrdis[0].y, (float) d);
		// lPoint.setPosiInSection(j - 1);
		// left = false;
		// }
		// Point pp0 = ps.get(size - j);
		// Point pp1 = ps.get(size - j - 1);
		// if (d > pp1.getDem() && d < pp0.getDem()) {
		// float dd = (float) (d - pp1.getDem()) / (pp0.getDem() - pp1.getDem());
		// if (Float.isInfinite(dd) || Float.isNaN(dd)) {
		// dd = 0f;
		// }
		// corrdis[1] = ShpAnalyseTool.coordOfLinesegment(pp1.getCoordinate(),
		// pp0.getCoordinate(), dd);
		// rPoint = new Point(corrdis[1].x, corrdis[1].y, (float) d);
		// rPoint.setPosiInSection(size - j - 1);
		// right = false;
		// }
		// }

		if (lPoint == null) {
			lPoint = ps.get(0);
			lPoint.setPosiInSection(0);
		}
		if (rPoint == null) {
			rPoint = ps.get(size - 1);
			rPoint.setPosiInSection(size - 1);
		}
		result[0] = calcDistance(lPoint, rPoint);
		result[1] = lPoint;
		result[2] = rPoint;
		result[3] = d;
		return result;
	}

	public static Object[] getSectionWidthByZFromSide(Section section, double z) {
		Coordinate[] corrdis = new Coordinate[2];
		Object[] result = new Object[4];

		List<Point> ps = section.getSectionLine();
		double d = z;

		Point lPoint = null, rPoint = null;
		int size = ps.size() - 1;
		int riverIndex = section.getRiverPoint().getPosiInSection();
		boolean left, right;
		left = right = false;
		for (int i = 0; i < size; i++) {
			Point lp = ps.get(i);
			Point lpNext = null;
			if (i + 1 > riverIndex) {
				// if (rPoint != null && i + 1 > rPoint.getPosiInSection()) {
				left = true;
			} else {
				lpNext = ps.get(i + 1);
			}
			Point rp = null;
			Point rpPre = null;
			if (size - i < ps.size() && size - i - 1 > riverIndex) {
				// if (size - i < ps.size() && (lPoint == null || size - i - 1 >
				// lPoint.getPosiInSection())) {
				rp = ps.get(size - i);
				rpPre = ps.get(size - i - 1);
			} else {
				right = true;
			}
			if (!left && Double.compare(lp.getDem(), d) >= 0 && Double.compare(lpNext.getDem(), d) <= 0) {
				float dd = (float) ((d - lpNext.getDem()) / (lp.getDem() - lpNext.getDem()));
				if (Float.isInfinite(dd) || Float.isNaN(dd)) {
					dd = 0f;
				}
				left = true;
				corrdis[0] = ShpAnalyseTool.coordOfLinesegment(lpNext.getCoordinate(), lp.getCoordinate(), dd);
				lPoint = new Point(corrdis[0].x, corrdis[0].y, (float) d);

				lPoint.setPosiInSection(i);
			}
			if (!right && Double.compare(rp.getDem(), d) >= 0 && Double.compare(rpPre.getDem(), d) <= 0) {
				float dd = (float) ((d - rpPre.getDem()) / (rp.getDem() - rpPre.getDem()));
				if (Float.isInfinite(dd) || Float.isNaN(dd)) {
					dd = 0f;
				}
				right = true;
				corrdis[0] = ShpAnalyseTool.coordOfLinesegment(rpPre.getCoordinate(), rp.getCoordinate(), dd);
				rPoint = new Point(corrdis[0].x, corrdis[0].y, (float) d);
				rPoint.setPosiInSection(size - i);

			}
			if (left && right) {
				break;
			}
		}
		if (lPoint == null) {
			lPoint = ps.get(0);
			lPoint.setPosiInSection(0);
		}
		if (rPoint == null) {
			rPoint = ps.get(size);
			rPoint.setPosiInSection(size);
		}
		double length = 0;
		int startIndex = lPoint.getPosiInSection();
		int endIndex = rPoint.getPosiInSection();

		length += calcDistance(lPoint, section.getPointOfSectionLine(startIndex + 1));
		length += calcDistance(rPoint, section.getPointOfSectionLine(endIndex - 1));

		for (int i = startIndex + 1; i < endIndex - 1; i++) {
			Point point = section.getPointOfSectionLine(i);
			Point nPoint = section.getPointOfSectionLine(i + 1);
			// 高出水面
			if (Double.compare(point.getDem(), d) > 0 && Double.compare(nPoint.getDem(), d) > 0) {
				continue;
			}
			if (Double.compare(point.getDem(), d) >= 0) {
				Point p = nPoint;
				float dd = (float) ((d - p.getDem()) / (point.getDem() - p.getDem()));
				if (Float.isInfinite(dd) || Float.isNaN(dd)) {
					dd = 0f;
				}
				Coordinate coor = ShpAnalyseTool.coordOfLinesegment(point.getCoordinate(), p.getCoordinate(), dd);
				length += calcDistance(coor, p.getCoordinate());

			} else if (Double.compare(nPoint.getDem(), d) >= 0) {
				Point p = point;
				float dd = (float) ((d - p.getDem()) / (nPoint.getDem() - p.getDem()));
				if (Float.isInfinite(dd) || Float.isNaN(dd)) {
					dd = 0f;
				}
				Coordinate coor = ShpAnalyseTool.coordOfLinesegment(p.getCoordinate(), nPoint.getCoordinate(), dd);
				length += calcDistance(coor, p.getCoordinate());
			} else {
				length += calcDistance(point, nPoint);
			}
		}

		result[0] = length;
		result[1] = lPoint;
		result[2] = rPoint;
		result[3] = d;
		return result;
	}

	/**
	 * 两点距离
	 * 
	 * @param lPoint
	 * @param rPoint
	 * @return
	 */
	private static double calcDistance(Point lPoint, Point rPoint) {
		PointPairDistance pd = new PointPairDistance();
		pd.initialize(new Coordinate(lPoint.getLon(), lPoint.getLat()),
				new Coordinate(rPoint.getLon(), rPoint.getLat()));
		double len = pd.getDistance();
		return len;
	}

	private static double calcDistance(Coordinate p0, Coordinate p1) {
		PointPairDistance pd = new PointPairDistance();
		pd.initialize(p0, p1);
		double len = pd.getDistance();
		return len;
	}

	/**
	 * 计算复式断面面积 <br>
	 * 计算面积时,将断面根据断面点切割成两个三角形和多个梯形
	 * 
	 * @param section
	 * @param z
	 * @return [0] 复式断面面积之和<br>
	 *         [1] 复式断面流量模数之和<br>
	 *         [2] 复式断面面积与流量模数比率之和
	 */
	public static Object[] calcSectionArea(Section section, Point leftPoint, Point rightPoint, double z) {
		/*
		 * xiaochen Object result[] = new Object[3]; double area = 0; double modulus =
		 * 0; double modulusRatio = 0; List<Point> points = section.getSectionLine();
		 * 
		 * Point firstPoint = points.get(0); Point lastPoint = points.get(points.size()
		 * - 1);
		 * 
		 * float avgDem = (leftPoint.getDem() + rightPoint.getDem()) / 2; // TODO
		 * 应求水位线高程
		 * 
		 * double lDistance = calcDistance(firstPoint, leftPoint); //TODO 感觉这两个三角形多余算。
		 * double rDistance = calcDistance(lastPoint, rightPoint);
		 * 
		 * double a = Formula.triangleArea(firstPoint.getDem()-avgDem, lDistance); area
		 * += a; double k = Formula.hydromodulus(1, a / lDistance, a); double kr = k * k
		 * * k / (a * a); modulus += k; modulusRatio += kr; a =
		 * Formula.triangleArea(lastPoint.getDem()-avgDem, rDistance); area += a; k =
		 * Formula.hydromodulus(1, a / rDistance, a); kr = k * k * k / (a * a); modulus
		 * += k; modulusRatio += kr; int leftIndex = leftPoint.getPosiInSection(); int
		 * rightIndex = rightPoint.getPosiInSection(); for (int i = leftIndex+1; i <
		 * rightIndex; i++)//TODO 少算了2个小的三角形，暂时先不算。 { Point lp = points.get(i); Point rp
		 * = points.get(i+1); double h = calcDistance(lp, rp); float lh = avgDem -
		 * lp.getDem(); float rh = avgDem - rp.getDem(); a = Formula.trapezoidArea(lh,
		 * rh, h); area += a; k = Formula.hydromodulus(1, a / h, a); kr = k * k * k / (a
		 * * a); modulus += k; modulusRatio += kr;
		 * 
		 * } result[0] = area; result[1] = modulus; result[2] = modulusRatio; return
		 * result;
		 */
		double n = section.getRoughRate();
		Object result[] = new Object[3];
		double area = 0; // 总面积
		double modulus = 0;
		double modulusRatio = 0;
		List<Point> points = section.getSectionLine();
		int leftIndex = leftPoint.getPosiInSection() + 1;
		int rightIndex = rightPoint.getPosiInSection() - 1;
		// float avgDem = (leftPoint.getDem() + rightPoint.getDem()) / 2;
		float avgDem = (float) z;

		double a, k, kr;

		Point firstPoint = points.get(leftIndex);
		Point lastPoint = points.get(rightIndex);

		double lDistance = calcDistance(firstPoint, leftPoint);
		double rDistance = calcDistance(lastPoint, rightPoint);

		a = Formula.triangleArea(avgDem - firstPoint.getDem(), lDistance);
		if (Double.compare(a, 0d) > 0)// 小于0说明中间有超出水位的点
		{
			area += a;
			k = Formula.hydromodulus(n, a / lDistance, a);
			kr = k * k * k / (a * a);
			modulus += k;
			modulusRatio += kr;
		}
		a = Formula.triangleArea(avgDem - lastPoint.getDem(), rDistance);
		if (Double.compare(a, 0d) > 0) {
			area += a;
			k = Formula.hydromodulus(n, a / rDistance, a);
			kr = k * k * k / (a * a);
			modulus += k;
			modulusRatio += kr;
		}

		for (int i = leftIndex; i < rightIndex; i++) {
			Point lp = points.get(i);
			Point rp = points.get(i + 1);
			float lh = avgDem - lp.getDem();
			float rh = avgDem - rp.getDem();
			if (lh < 0 || rh < 0) // 在河道中央出现超出水面的点，所以用三角形算。
			{

			} else // 梯形
			{
				double h = calcDistance(lp, rp);
				a = Formula.trapezoidArea(lh, rh, h);
				if (Double.compare(a, 0d) > 0) {
					area += a;
					k = Formula.hydromodulus(n, a / h, a);
					kr = k * k * k / (a * a);
					modulus += k;
					modulusRatio += kr;
				}
			}
		}
		result[0] = area;
		result[1] = modulus;
		result[2] = modulusRatio;
		return result;
	}

	/**
	 * 根据水动能平衡公式推求水位
	 * 
	 * @param sections
	 * @param n
	 * @return
	 */
	public static double[] calcZByWaterEnergyBalance(List<Section> sections) {
		double[] result = new double[sections.size()];
		result[sections.size() - 1] = sections.get(sections.size() - 1).getMaxZ(); // 出口断面水位
		Logger.info(String.format("出口断面的水位%s，最低点高程%s", result[sections.size() - 1],
				sections.get(sections.size() - 1).getMinAltitudePoint().getDem()));

		for (int i = sections.size() - 2; i >= 0; i--) { // 倒数第二个开始
			Logger.info(String.format("============断面 %s================", i));
			Section section = sections.get(i);
			double z = calcZByWaterEnergyBalance(sections, i, section.getMaxQ());
			sections.get(i).setMaxZ(z);
			sections.get(i).setStartIndex(lastStartIndex);
			sections.get(i).setEndIndex(lastEndIndex);
			result[i] = z;
			Logger.info(String.format("============断面 %s================", i));
		}

		for (int i = sections.size() - 2; i >= 0; i--) { // 倒数第二个开始
			Section section = sections.get(i);
			Logger.info(String.format("断面%s的水位%s，最低点高程%s", i, section.getMaxZ(),
					sections.get(i).getMinAltitudePoint().getDem()));
		}
		return result;
	}

	private static double calcZByWaterEnergyBalance(List<Section> sections, int index, double q) {
		double precision = (double) cfg.get("cacc");
		double result = 0;
		Section section = sections.get(index);
		if (index == sections.size() - 1) {
			return calcZByQ(section.getWaterVolumeMap(), q);
		}
		Section nextSection = sections.get(index + 1);
		Point nextMinAltitudePoint = nextSection.getMinAltitudePoint();
		Point minAltitudePoint = section.getMinAltitudePoint();
		// Point minAltitudePoint = nextSection.getRiverPoint();
		double nextSectionZ = calcZByWaterEnergyBalance(sections, index + 1, q);
		// if (nextSectionZ < 0) // 说明下游断面的水位还没有达到本断面的最低点高程。
		// {
		// return result;
		// }
		double z = nextSection.getMaxZ(); // TODO // 初始值
		z += sections.get(sections.size() - 1).getGradient() * 1000;
		if (z < 0 || z < section.getMinAltitudePoint().getDem()) {
			z = section.getMinAltitudePoint().getDem() + 0.1d;
		}
		Logger.debug(String.format("断面%s,初始水深%s,指定流量%s", index, z, q));
		// double z = 0.1;
		double e1 = 0, e2 = 0;
		double deltaZ = 0;

		Object[] nextWidthResult = getSectionWidthByZ(nextSection, nextSectionZ, 0); // TODO
																						// 需要优化，应该从map里去查询，不应该反复计算。
		double nextSectionWidth = (double) nextWidthResult[0];
		Object[] areaResult = calcSectionArea(nextSection, (Point) nextWidthResult[1], (Point) nextWidthResult[2],
				nextSectionZ);
		double nextSectionArea = (double) areaResult[0];
		double nextSectionQ = q;
		double nextSectionV = nextSectionQ / nextSectionArea;
		double nextSectionModulus = (double) areaResult[1];
		double nextFactor = Formula.kineticEnergyFactor(nextSectionArea, (double) areaResult[2], nextSectionModulus);
		double nextSectionH = nextSectionArea / nextSectionWidth;
		double j2 = Formula.hydraulicSlope(nextSection.getRoughRate(), nextSectionV, nextSectionH);

		Logger.debug(String.format(
				"next section width %s ,secation area %s,section v %s,factor %s, section H %s ,section z %s , bottom dem %s",
				nextSectionWidth, nextSectionArea, nextSectionV, nextFactor, nextSectionH, nextSectionZ,
				nextMinAltitudePoint.getDem()));

		// Logger.debug(
		// String.format("next section z %s , bottom dem %s , next section area %s
		// ,secation width %s,factor %s ",
		// nextSectionZ, nextMinAltitudePoint.getDem(), nextSectionArea,
		// nextSectionWidth, nextFactor));
		Object[] widthResult = null;
		int testCount = 0;
		int testCount2 = 0;
		int testCount3 = 0;
		boolean con = false;
		int calcCount = 0;
		do {
			if (testCount3 > 2000) { // 连续多次，都算不出来
				List<ZQ> zqList = section.getWaterVolumeMap();
				if (zqList == null || zqList.size() == 0) {
					zqList = calcZQListByManning(sections, index, sections.get(0).getMaxQ() * 2);
					section.setWaterVolumeMap(zqList);
				}
				z = calcZByQ(zqList, q);
				widthResult = getSectionWidthByZ(section, z, 0);
				Logger.warn("calc multi times get current z by manning " + index);
				break;
			}
			calcCount++;
			z -= (con ? 0 : deltaZ);
			if (z < 0) {
				z = Math.abs(deltaZ);
				testCount2++;
				if (testCount2 > 100) {
					Logger.warn("negative z ,calc z by Manning " + index);
					List<ZQ> zqList = section.getWaterVolumeMap();
					if (zqList == null || zqList.size() == 0) {
						zqList = calcZQListByManning(sections, index, sections.get(0).getMaxQ() * 2);
						section.setWaterVolumeMap(zqList);
					}
					z = calcZByQ(zqList, q);
					widthResult = getSectionWidthByZ(section, z, 0);
					break;
				}
			} else {
				testCount2 = 0;
			}
			if (con) {
				con = false;
			}
			widthResult = getSectionWidthByZ(section, z, 0);
			double sectionWidth = (double) widthResult[0];
			areaResult = calcSectionArea(section, (Point) widthResult[1], (Point) widthResult[2], z);
			double sectionArea = (double) areaResult[0];
			if (Double.compare(sectionArea, 0d) == 0) { // 假设水位低于断面高程
				z += 50;
				testCount3++;
				con = true;
				continue;
			}
			double sectionQ = q;
			double sectionV = sectionQ / sectionArea;
			double sectionModulus = (double) areaResult[1];
			double factor = Formula.kineticEnergyFactor(sectionArea, (double) areaResult[2], sectionModulus);
			double sectionH = sectionArea / sectionWidth;
			double j1 = Formula.hydraulicSlope(section.getRoughRate(), sectionV, sectionH);
			Logger.debug(String.format(
					"section width %s ,secation area %s,section v %s,factor %s, section H %s ,section Z %s , bottom dem %s",
					sectionWidth, sectionArea, sectionV, factor, sectionH, z, minAltitudePoint.getDem()));
			// Logger.debug(String.format("section z %s , bottom dem %s , section area %s
			// ,secation width %s,factor %s ",
			// z, minAltitudePoint.getDem(), sectionArea, sectionWidth, factor));
			double fr = Formula.froudeNumber(sectionH, sectionV);
			if (fr >= 1) {// 急流情况下，使用曼宁公式推求
				// List<ZQ> zqList = section.getWaterVolumeMap();
				// if (zqList == null || zqList.size() == 0) {
				// zqList = calcZQListByManning(sections, index, sections.get(0).getMaxQ() * 2);
				// section.setWaterVolumeMap(zqList);
				// }
				// testCount++;
				Logger.debug(String.format("fr>1 ,z %s , fr %s ", z, fr));
				if (testCount > 1000) {
					Logger.warn("multi fr>1 calc z by Manning " + index);
					List<ZQ> zqList = section.getWaterVolumeMap();
					if (zqList == null || zqList.size() == 0) {
						zqList = calcZQListByManning(sections, index, sections.get(0).getMaxQ() * 2);
						section.setWaterVolumeMap(zqList);
					}
					z = calcZByQ(zqList, q);
					widthResult = getSectionWidthByZ(section, z, 0);
					break;
				}
				z += 0.5;
				testCount++;
				con = true;
				continue;

			} else {
				// testCount = 0;
			}

			double j = (j1 + j2) / 2;

			Point nextPoint = nextSection.getMinAltitudePoint();
			Point point = section.getMinAltitudePoint();
			double l = calcDistance(nextPoint, point);
			double deltaX = nextSection.getDistanceOffStart() - section.getDistanceOffStart();
			double hf = Formula.frictionalHeadLoss(j, deltaX);
			double hj = 0; // TODO

			e1 = Formula.waterEnergyBalance(z, factor, sectionV, 0, 0);
			e2 = Formula.waterEnergyBalance(nextSectionZ, nextFactor, nextSectionV, hf, hj);

			double deltaE = e1 - e2;
			deltaZ = deltaE / (1 - factor * fr * fr + 3 * j * deltaX / (2 * sectionH));
			Logger.debug(String.format(
					"j %s, j1 %s, j2 %s , deltaE %s , deltaZ %s , z %s , e1 %s, e2 %s , hf %s , fr %s ,l %s", j, j1, j2,
					deltaE, deltaZ, z, e1, e2, hf, fr, l));
			testCount3++;

		} while (Double.compare(Math.abs(e2 - e1), precision) > 0 || Double.compare(e1, 0) == 0);
		// result = section.getMinAltitudePoint().getDem() + z;
		// result = section.getRiverPoint().getDem() + z;
		result = z;
		Logger.debug(String.format("本次试算%s次,试算水位%s", calcCount, z));

		lastStartIndex = ((Point) widthResult[1]).getPosiInSection();
		lastEndIndex = ((Point) widthResult[2]).getPosiInSection();
		return result;
	}

	public static Coordinate[] waterSurfaceCoordinates1(double[] waterDemArray, List<Section> sectionList)
			throws ParseException {
		Coordinate[] coor = new Coordinate[sectionList.size() * 2 + 1];
		for (int i = 0; i < sectionList.size(); i++) {
			Section s = sectionList.get(i);
			int st = s.getStartIndex();
			int ed = s.getEndIndex();
			Point sp = s.getPointOfSectionLine(st);
			Point ep = s.getPointOfSectionLine(ed);
			coor[i] = new Coordinate(sp.getLon(), sp.getLat());
			coor[coor.length - i - 2] = new Coordinate(ep.getLon(), ep.getLat());
			// GeometryFactory geometryFactory = new GeometryFactory();
			// Logger.debug(geometryFactory.createPoint(coor[i]).toText());
			// Logger.debug(geometryFactory.createPoint(coor[coor.length - i -
			// 2]).toText());
		}
		coor[coor.length - 1] = coor[0];
		return coor;
	}

	/**
	 * 根据各断面水位，在各断面上取点（是水面线的点序列）
	 * 
	 * @param waterDemArray
	 * @param sectionList
	 * @return
	 * @throws ParseException
	 */
	public static Coordinate[] waterSurfaceCoordinates(double[] waterDemArray, List<Section> sectionList)
			throws ParseException {
		Coordinate[] corrdis = new Coordinate[waterDemArray.length * 2 + 1];
		int length = corrdis.length - 1;
		for (int i = 0; i < waterDemArray.length; i++) {
			double d = waterDemArray[i];
			List<Point> ps = sectionList.get(i).getSectionLine();
			int size = ps.size();
			boolean left = true, right = true;
			for (int j = 1; j < size; j++) {
				if (left || right) {
					Point p0 = ps.get(j - 1);
					Point p1 = ps.get(j);
					if (d > p1.getDem() && left) {
						float dd = p0.getDem() - p1.getDem();
						corrdis[i] = ShpAnalyseTool.coordOfLinesegment(p1.getCoordinate(), p0.getCoordinate(), dd);
						left = false;
					}

					Point pp0 = ps.get(size - j);
					Point pp1 = ps.get(size - j - 1);
					if (d > pp1.getDem() && right) {
						float dd = pp0.getDem() - pp1.getDem();
						corrdis[length - 1 - i] = ShpAnalyseTool.coordOfLinesegment(pp1.getCoordinate(),
								pp0.getCoordinate(), dd);
						right = false;
					}
				} else {
					break;
				}
			}
		}
		corrdis[length] = corrdis[0];
		return corrdis;
	}
}
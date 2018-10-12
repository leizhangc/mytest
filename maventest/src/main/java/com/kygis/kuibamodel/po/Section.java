package com.kygis.kuibamodel.po;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.LineString;

/**
 * 断面类
 * 
 * @author wheeler
 */
public class Section {
	private String nid;
	private LineString lineString;// 断面平面线（不包含高程），存成LineString类型。
	private String damParamID;// 溃口记录ID
	private String sectionParamID;// 断面提取参数记录ID
	private List<Point> sectionLine;// 断面底线 。用GSON（google json）是否可以直接转，得实验一下。
	private double gradient;// 作为断面的比降
	private Map<Float, WaterSectionFeature> waterSectionFeatureMap = new HashMap<Float, WaterSectionFeature>();// 水位及对应的过水断面特征
	private List<ZQ> waterVolumeMap = new ArrayList<>();// 水位-流量关系,用于出口断面或急流断面
	private double distanceOffStart;// 离溃口的距离。=断面间距*断面序号。
	private double maxQ; // 溃口最大流量
	private double maxZ; // 断面最高水位
	private double roughRate; // 糙率
	private int startIndex = -1; // 最高水位左端点
	private int endIndex = -1; // 最高水位右端点
	private Point riverPoint; // 河道点

	public Section(List<Point> sectionLine) {
		this.sectionLine = sectionLine;
	}

	public Point getPointOfSectionLine(int i) {
		return sectionLine.get(i);
	}

	public int getIndexOfPoint(Point other) {
		int size = sectionLine.size();
		for (int i = 0; i < size; i++) {
			if (sectionLine.get(i).equalesOfLonLat(other)) {
				return i;
			}
		}
		return -1;
	}

	public List<Point> getSectionLine() {
		return sectionLine;
	}

	public String getNid() {
		return nid;
	}

	public String getDamParamID() {
		return damParamID;
	}

	public String getSectionParamID() {
		return sectionParamID;
	}

	public void setNid(String nid) {
		this.nid = nid;
	}

	public void setDamParamID(String damParamID) {
		this.damParamID = damParamID;
	}

	public void setSectionParamID(String sectionParamID) {
		this.sectionParamID = sectionParamID;
	}

	public void setSectionLine(List<Point> sectionLine) {
		this.sectionLine = sectionLine;
	}

	public LineString getLineString() {
		return lineString;
	}

	public void setLineString(LineString lineString) {
		this.lineString = lineString;
	}

	public double getGradient() {
		return gradient;
	}

	public void setGradient(double gradient) {
		this.gradient = gradient;
	}

	public Map<Float, WaterSectionFeature> getWaterSectionFeatureMap() {
		return waterSectionFeatureMap;
	}

	public WaterSectionFeature getWaterSectionFeatureMap(Float waterDem) {
		return waterSectionFeatureMap.get(waterDem);
	}

	public void addWaterSectionFeature(Float waterDem, WaterSectionFeature waterSectionFeature) {
		waterSectionFeatureMap.put(waterDem, waterSectionFeature);
	}

	public List<ZQ> getWaterVolumeMap() {
		return waterVolumeMap;
	}

	public void setWaterVolumeMap(List<ZQ> list) {
		waterVolumeMap = list;
	}

	public double getDistanceOffStart() {
		return distanceOffStart;
	}

	public void setDistanceOffStart(double distanceOffStart) {
		this.distanceOffStart = distanceOffStart;
	}

	public double getMaxQ() {
		return maxQ;
	}

	public void setMaxQ(double maxQ) {
		this.maxQ = maxQ;
	}

	public double getMaxZ() {
		return maxZ;
	}

	public void setMaxZ(double maxZ) {
		this.maxZ = maxZ;
	}

	/**
	 * 获取断面点序列中海拔最低点
	 * 
	 * @param points
	 * @return 最低点
	 */
	public Point getMinAltitudePoint() {
		if (startIndex > -1 && endIndex > -1) {
			return getMinAltitudePoint(startIndex, endIndex);
		} else {
			return getMinAltitudePoint(0, sectionLine.size() - 1);
		}
	}

	public Point getMinAltitudePoint(int start, int end) {
		Point result = sectionLine.get(start);
		for (int i = start + 1; i < end + 1; i++) {
			float h = sectionLine.get(i).getDem();
			if (Float.compare(h, result.getDem()) < 0) {
				result = sectionLine.get(i);
			}
		}
		return result;
	}

	public double getRoughRate() {
		return roughRate;
	}

	public void setRoughRate(double roughRate) {
		this.roughRate = roughRate;
	}

	public int getStartIndex() {
		return startIndex;
	}

	public void setStartIndex(int startIndex) {
		/*
		 * if(this.startIndex<=-1) { this.startIndex = startIndex; }
		 */
		this.startIndex = startIndex;
	}

	public int getEndIndex() {
		return endIndex;
	}

	public void setEndIndex(int endIndex) {
		/*
		 * if(this.endIndex<=-1) { this.endIndex = endIndex; }
		 */
		this.endIndex = endIndex;
	}

	public Point getRiverPoint() {
		return riverPoint;
	}

	public void setRiverPoint(Point riverPoint) {
		this.riverPoint = riverPoint;
	}

	public void printPoint(int start, int end) {
		for (int i = start; i < end + 1; i++) {
			Point p = sectionLine.get(i);
			System.out.println(String.format("%s,%s,%s,%s", i, p.getLon(), p.getLat(), p.getDem()));
		}
	}

	/**
	 * 
	 * @return
	 */
	public boolean restoreAsStartAndEnd(int start, int end) {
		List<Point> sectionLineX = new ArrayList<Point>();
		Coordinate[] coordinates = new Coordinate[end - start + 1];
		for (int i = start; i < end + 1; i++) {
			sectionLineX.add(sectionLine.get(i));
			coordinates[i - start] = new Coordinate(sectionLine.get(i).getLon(), sectionLine.get(i).getLat());
		}
		sectionLine = sectionLineX;

		// Coordinate[] coods = lineString.getCoordinates();
		// Coordinate[] coodsX = new Coordinate[end - start + 1];
		// System.arraycopy(coods, start, coodsX, 0, coodsX.length);
		lineString = new GeometryFactory().createLineString(coordinates);

		startIndex = 0;
		endIndex = end - start;
		return true;
	}

	/**
	 * 根据断面两端点序号进行裁剪。
	 * 
	 * @return
	 */
	public boolean restore() {
		if (startIndex > -1 && endIndex > -1) {
			restoreAsStartAndEnd(startIndex, endIndex);
		}
		return true;
	}

}
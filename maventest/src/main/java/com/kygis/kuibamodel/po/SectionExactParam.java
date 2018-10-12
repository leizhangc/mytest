package com.kygis.kuibamodel.po;

/**
 * 断面分析提取时参数配置类
 * 
 * @author wheeler
 */
public class SectionExactParam {
	private String nid;
	private int step;// 断面间距，单位公里
	private int distance;// 在平面上看，断面线上点的距离，单位m
	private double[] endPoint;// 在河道上的坐标位置，数组长度为2，出口断面位置
	private double[] startPoint;
	private float topDem;// 断面最高高程。//一般等于大坝水位
	private float[] onlyOneSideLowParam;// 一边低断面修正参数，一个是与最低点高程的差、一个是低一边继续延长的长度。
	private int pointCount;

	public int getStep() {
		return step;
	}

	public int getDistance() {
		return distance;
	}

	public double[] getEndPoint() {
		return endPoint;
	}

	public void setStep(int step) {
		this.step = step;
	}

	public void setDistance(int distance) {
		this.distance = distance;
	}

	public void setEndPoint(double[] endPoint) {
		this.endPoint = endPoint;
	}

	public String getNid() {
		return nid;
	}

	public float getTopDem() {
		return topDem;
	}

	public void setNid(String nid) {
		this.nid = nid;
	}

	public void setTopDem(float topDem) {
		this.topDem = topDem;
	}

	public double[] getStartPoint() {
		return startPoint;
	}

	public void setStartPoint(double[] startPoint) {
		this.startPoint = startPoint;
	}

	public float[] getOnlyOneSideLowParam() {
		return onlyOneSideLowParam;
	}

	public void setOnlyOneSideLowParam(float[] onlyOneSideLowParam) {
		this.onlyOneSideLowParam = onlyOneSideLowParam;
	}

	public int getPointCount() {
		return pointCount;
	}

	public void setPointCount(int pointCount) {
		this.pointCount = pointCount;
	}
}

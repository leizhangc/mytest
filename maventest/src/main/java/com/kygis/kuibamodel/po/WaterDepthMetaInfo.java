package com.kygis.kuibamodel.po;

public class WaterDepthMetaInfo {
	private String nid;// 记录ID
	private String damParamID;// 溃口记录ID
	private String sectionParamID;// 断面提取参数记录ID
	private String waterSurfaceID;// 水面线记录ID
	private String path;// 水深栅格文件存放路径
	private String url;// 水深栅格文件服务URL（以后会用到）
	private float maxDepth;// 最深值
	private float minDepth;// 最浅值
	private double area; // 淹没范围面积
	private double minX; // 淹没范围最小X
	private double minY; // 淹没范围最小Y
	private double maxX; // 淹没范围最大X
	private double maxY; // 淹没范围最大Y

	public String getNid() {
		return nid;
	}

	public void setNid(String nid) {
		this.nid = nid;
	}

	public String getDamParamID() {
		return damParamID;
	}

	public void setDamParamID(String damParamID) {
		this.damParamID = damParamID;
	}

	public String getSectionParamID() {
		return sectionParamID;
	}

	public void setSectionParamID(String sectionParamID) {
		this.sectionParamID = sectionParamID;
	}

	public String getWaterSurfaceID() {
		return waterSurfaceID;
	}

	public void setWaterSurfaceID(String waterSurfaceID) {
		this.waterSurfaceID = waterSurfaceID;
	}

	public String getPath() {
		return path;
	}

	public void setPath(String path) {
		this.path = path;
	}

	public String getUrl() {
		return url;
	}

	public void setUrl(String url) {
		this.url = url;
	}

	public float getMaxDepth() {
		return maxDepth;
	}

	public void setMaxDepth(float maxDepth) {
		this.maxDepth = maxDepth;
	}

	public float getMinDepth() {
		return minDepth;
	}

	public void setMinDepth(float minDepth) {
		this.minDepth = minDepth;
	}

	public double getArea() {
		return area;
	}

	public void setArea(double area) {
		this.area = area;
	}

	public double getMinX() {
		return minX;
	}

	public void setMinX(double minX) {
		this.minX = minX;
	}

	public double getMinY() {
		return minY;
	}

	public void setMinY(double minY) {
		this.minY = minY;
	}

	public double getMaxX() {
		return maxX;
	}

	public void setMaxX(double maxX) {
		this.maxX = maxX;
	}

	public double getMaxY() {
		return maxY;
	}

	public void setMaxY(double maxY) {
		this.maxY = maxY;
	}
}
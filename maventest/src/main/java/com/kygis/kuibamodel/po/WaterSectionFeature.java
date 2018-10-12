package com.kygis.kuibamodel.po;

/**
 * 过水断面特征值类
 * 
 * @author KAIFA01
 *
 */
public class WaterSectionFeature {
	private double width;// 水面宽
	private double averageDepth;// 平均水深
	private double z; // 水深   TODO  为什么需要这个水深呢？？？？ @晓晨
	private double area;// 过水面积

	public double getWidth() {
		return width;
	}

	public void setWidth(double width) {
		this.width = width;
	}

	public double getAverageDepth() {
		return averageDepth;
	}

	public void setAverageDepth(double averageDepth) {
		this.averageDepth = averageDepth;
	}

	public double getArea() {
		return area;
	}

	public void setArea(double area) {
		this.area = area;
	}

	public double getZ() {
		return z;
	}

	public void setZ(double z) {
		this.z = z;
	}

}

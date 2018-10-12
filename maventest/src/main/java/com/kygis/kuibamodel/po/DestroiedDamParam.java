package com.kygis.kuibamodel.po;

/**
 * 溃坝参数类
 * 
 * @author wheeler
 *
 */
public class DestroiedDamParam {
	private String nid;// 记录ID
	private double width;// 溃口宽度
	private double height;// 溃口高度
	private float damTopDem;// 坝顶高程
	private double waterVolume;// 溃坝时水量(是指李斯特公式里Vw为水库溃坝后的下泄水量体积)
	private double waterThroughput;// 溃口最大流量（用宽口堰公式计算得到）
	private double roughRate;// 糙率
	private double cacc=0.001; // 收敛精度

	public String getNid() {
		return nid;
	}

	public void setNid(String nid) {
		this.nid = nid;
	}

	public double getWidth() {
		return width;
	}

	public double getHeight() {
		return height;
	}

	public double getDamTopDem() {
		return damTopDem;
	}

	public double getWaterVolume() {
		return waterVolume;
	}

	public void setWidth(double width) {
		this.width = width;
	}

	public void setHeight(double height) {
		this.height = height;
	}

	public void setDamTopDem(float damTopDem) {
		this.damTopDem = damTopDem;
	}

	public void setWaterVolume(double waterVolume) {
		this.waterVolume = waterVolume;
	}

	public double getWaterThroughput() {
		return waterThroughput;
	}

	public void setWaterThroughput(double waterThroughput) {
		this.waterThroughput = waterThroughput;
	}

	public double getRoughRate() {
		return roughRate;
	}

	public void setRoughRate(double roughRate) {
		this.roughRate = roughRate;
	}

	public double getCacc() {
		return cacc;
	}

	public void setCacc(double cacc) {
		this.cacc = cacc;
	}
}

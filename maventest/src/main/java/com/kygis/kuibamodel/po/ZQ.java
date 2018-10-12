package com.kygis.kuibamodel.po;

/**
 * 水位流量对应关系
 * 
 * @author nico
 *
 */
public class ZQ {
	private double z;
	private double q;
	private double h;

	public ZQ(double z, double q, double h) {
		this.z = z;
		this.q = q;
		this.h = h;
	}

	public double getZ() {
		return z;
	}

	public void setZ(double z) {
		this.z = z;
	}

	public double getQ() {
		return q;
	}

	public void setQ(double q) {
		this.q = q;
	}

	public double getH() {
		return h;
	}

	public void setH(double h) {
		this.h = h;
	}
}

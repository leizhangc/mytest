package com.kygis.algorithmofdam;

/**
 * 插值
 * 
 * @author nico
 *
 */
final class Interpolation {

	/**
	 * 线性插值
	 * 
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 * @param x
	 * @return
	 */
	public static double linear(double x1, double y1, double x2, double y2, double x) {
		return (x - x2) / (x1 - x2) * y1 + (x - x1) / (x2 - x1) * y2;
	}
}

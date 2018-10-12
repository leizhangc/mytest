package com.kygis.algorithmofdam;

/**
 * 公式
 * 
 * @author nico
 *
 */
final class Formula {
	public static double g = 9.81d;

	/**
	 * 宽顶堰公式
	 * 
	 * @param h
	 *            溃口高度
	 * @param lb
	 *            溃口宽度
	 * @return 最大流量
	 */
	public static double broadCrestedWeir(double h, double lb) {
		return 8d / 27d * lb * h * Math.sqrt(Formula.g * h);
	}

	/**
	 * τ 经验系数
	 * 
	 * @param v
	 *            断面平均流速，山区采用3.0，山前区用2.0～3.0，平原用1.0～2.0
	 * @param k0
	 *            为调整系数，山前区采用1.0，平原用0.8～0.9
	 * @return 系数值
	 */
	public static double tau(double v, double k0) {
		return 1 / (v * k0);
	}

	/**
	 * 李斯特公式
	 * 
	 * @param qmax
	 *            溃口最大流量
	 * @param distance
	 *            距溃决点距离
	 * @param vw
	 *            水库溃坝后的下泄水量体积
	 * @param t
	 *            经验系数
	 * @return 距溃决点为distance处的溃决洪水最大流量
	 */
	public static double liszt(double qmax, double distance, double vw, double t) {
		return qmax * 1 / (1 + ((t * qmax) / vw) * distance);
	}

	/**
	 * 李斯特公式
	 * 
	 * @param qmax
	 *            溃口最大流量
	 * @param distance
	 *            距溃决点距离
	 * @param vw
	 *            水库溃坝后的下泄水量体积
	 * @param v
	 *            断面平均流速
	 * @param k
	 *            调整系数
	 * @return 距溃决点为distance处的溃决洪水最大流量
	 */
	public static double liszt(double qmax, double distance, double vw, double v, double k) {
		double t = tau(v, k);
		return liszt(qmax, distance, vw, t);
	}

	/**
	 * 曼宁公式
	 * 
	 * @param n
	 *            断面的糙率
	 * @param h
	 *            平均水深
	 * @param j
	 *            水面比降,计算时取计算河段的河床平均比降
	 * @param b
	 *            水面宽度
	 * @return 流量
	 */
	public static double manning(double n, double h, double j, double b) {
		return 1 / n * Math.pow(h, 5d / 3d) * Math.pow(j, 1d / 2d) * b;
	}

	/**
	 * 河道比降
	 * 
	 * @param z1
	 *            上游河道高程
	 * @param z2
	 *            下游河道高程
	 * @param distance
	 *            z1,z2间距离
	 * @return
	 */
	public static double riverSlope(double z1, double z2, double distance) {
		return (z1 - z2) / distance;
	}

	/**
	 * 水力坡度
	 * 
	 * @param n
	 *            糙率
	 * @param v
	 *            水流速度
	 * @param r
	 *            水力半径
	 * @return 水力坡度
	 */
	public static double hydraulicSlope(double n, double v, double r) {
		return n * n * v * v / Math.pow(r, 4d / 3d);
	}

	/**
	 * 求梯形面积
	 * 
	 * @param top
	 * @param bottom
	 * @param h
	 * @return
	 */
	public static double trapezoidArea(double top, double bottom, double h) {
		return (top + bottom) * h / 2;
	}

	/**
	 * 求三角形面积
	 * 
	 * @param bottom
	 * @param h
	 * @return
	 */
	public static double triangleArea(double bottom, double h) {
		return bottom * h / 2;
	}

	/**
	 * 动能平衡公式
	 * 
	 * @param z
	 *            水位
	 * @param a
	 *            动能修正系数
	 * @param v
	 *            水流速
	 * @param hf
	 *            沿程水头损失
	 * @param hj
	 *            局部水头损失
	 * @return 上游水位
	 */
	public static double waterEnergyBalance(double z, double a, double v, double hf, double hj) {
		double result = 0;
		result = z + a * v * v / (2 * Formula.g) + hf + hj;
		return result;
	}

	/**
	 * 动能修正系数α
	 * 
	 * @param A
	 *            断面面积
	 * @param KR
	 *            流量模数与面积之比和
	 * @param K
	 *            断面各条流量模数之和
	 * @return 动能修正系数
	 */
	public static double kineticEnergyFactor(double A, double KR, double K) {
		double result = 0;
		result = A * A * KR / (K * K * K);
		return result;
	}

	/**
	 * 流量模数
	 * 
	 * @param n
	 *            糙率
	 * @param r
	 *            水力半径
	 * @param a
	 *            面积
	 * @return 流量模数
	 */
	public static double hydromodulus(double n, double r, double a) {
		return (1 / n) * Math.pow(r, 2d / 3d) * a;
	}

	/**
	 * 沿程水头损失
	 * 
	 * @param j
	 *            平均水力坡度
	 * @param l
	 *            断面距离
	 * @return
	 */
	public static double frictionalHeadLoss(double j, double l) {
		return j * l;
	}

	/**
	 * 弗汝德数
	 * 
	 * @param r
	 *            水力半径
	 * @param v
	 *            水流速度
	 * @return 弗汝德数
	 */
	public static double froudeNumber(double r, double v) {
		return v / Math.sqrt(Formula.g * r);
	}
}

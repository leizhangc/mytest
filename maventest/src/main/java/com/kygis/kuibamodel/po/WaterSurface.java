package com.kygis.kuibamodel.po;

import com.vividsolutions.jts.geom.Polygon;

public class WaterSurface
{
	private String nid;//记录ID
	private String damParamID;//溃口记录ID
	private String sectionParamID;//断面提取参数记录ID
	private Polygon wkt;//删除小碎图斑之后
	private Polygon wktConVex;//获取凸多边形
	private Polygon wktOriginal;//原始图
	public String getNid()
	{
		return nid;
	}
	public String getDamParamID()
	{
		return damParamID;
	}
	public String getSectionParamID()
	{
		return sectionParamID;
	}
	public Polygon getWkt()
	{
		return wkt;
	}
	public void setNid(String nid)
	{
		this.nid = nid;
	}
	public void setDamParamID(String damParamID)
	{
		this.damParamID = damParamID;
	}
	public void setSectionParamID(String sectionParamID)
	{
		this.sectionParamID = sectionParamID;
	}
	public void setWkt(Polygon wkt)
	{
		this.wkt = wkt;
	}
	public Polygon getWktConVex()
	{
		return wktConVex;
	}
	public void setWktConVex(Polygon wktConVex)
	{
		this.wktConVex = wktConVex;
	}
	public Polygon getWktOriginal()
	{
		return wktOriginal;
	}
	public void setWktOriginal(Polygon wktOriginal)
	{
		this.wktOriginal = wktOriginal;
	}
	
}

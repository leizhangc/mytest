package com.kygis.kuibamodel.po;

import com.vividsolutions.jts.geom.Coordinate;

public class Point
{
	private double lon;
	private double lat;
	private float dem;
	private int offsetPosiInSection=-1;//它在所属断面的点集合中，是在那个点之后。根据不同水位算出点时，这个点在断面上某个点之后，某个点的位置（集合中的序号）就是它。
	public Object tag;

	public Point(double lon,double lat,float dem)
	{
		this.lon=lon;
		this.lat=lat;
		this.dem=dem;
	}
	
	
	

	public int getPosiInSection()
	{
		return offsetPosiInSection;
	}




	public void setPosiInSection(int posiInSection)
	{
		this.offsetPosiInSection = posiInSection;
	}




	public double getLon()
	{
		return lon;
	}

	public double getLat()
	{
		return lat;
	}

	public float getDem()
	{
		return dem;
	}
	
	public Coordinate getCoordinate()
	{
		return new Coordinate(lon,lat); 
	}
	

	public boolean equalesOfLonLat(Point other)
	{
		if(this.lon==other.getLon()&&this.lat==other.getLat())
		{
			return true;
		}
		return false;
	}
	
}
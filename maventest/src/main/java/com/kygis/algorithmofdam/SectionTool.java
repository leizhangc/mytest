package com.kygis.algorithmofdam;

import java.util.List;

import org.geotools.coverage.grid.GridCoverage2D;
import com.kygis.kuibamodel.po.Section;
/**
 * 断面工具类。主要是解决断面难题（比如如何在10多公里断面上找出相对合理的那段凹槽来），对断面进行优化。
 * 因为不可控因素太多，所以优化结果不一定合理。
 * @author KAIFA01
 */
class SectionTool
{
	public boolean restoreSections(List<Section> sectionList,GridCoverage2D coverage)
	{
		//1、去掉多余的点
		
		return true;
	}
}

package com.kygis.algorithmofdam;

import java.util.List;

import org.opengis.referencing.operation.TransformException;

import com.kunyuan.framework.service.KYService;
import com.kygis.kuibamodel.po.Section;
import com.kygis.kuibamodel.po.SectionExactParam;
import com.kygis.kuibamodel.po.WaterDepthMetaInfo;
import com.kygis.kuibamodel.po.WaterSurface;

/**
 * 溃坝分析服务接口
 * 
 * @author wheeler
 *
 */
interface DamAnalyseService extends KYService {
	public boolean configBaseData(String demRaster, String riverWayShp, String boundShp);

	/**
	 * 断面提取
	 * 
	 * @param sectionExactParam
	 * @return
	 */
	public List<Section> sectionAnalyse(SectionExactParam sectionExactParam);

	/**
	 * 水面线提取
	 * 
	 * @param destroiedDamParam
	 * @return
	 */
	public WaterSurface waterSurfaceAnalyse(double[] waterDemArray, List<Section> allSection);

	/**
	 * 水位栅格图生成
	 * 
	 * @param destroiedDamParam
	 * @return
	 * @throws TransformException
	 */
	public WaterDepthMetaInfo waterDepthAnalyse(WaterSurface waterSurface, float bufferDistance,
			String storedAsFileName);

	public List<Section> sectionAnalyseWidthBound(SectionExactParam sectionExactParam);

}

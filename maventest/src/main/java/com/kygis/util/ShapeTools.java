package com.kygis.util;

import java.io.File;
import java.io.Serializable;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.Map;

import org.geotools.data.DataUtilities;
import org.geotools.data.DefaultTransaction;
import org.geotools.data.FeatureWriter;
import org.geotools.data.Transaction;
import org.geotools.data.shapefile.ShapefileDataStore;
import org.geotools.data.shapefile.ShapefileDataStoreFactory;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.data.simple.SimpleFeatureSource;
import org.geotools.data.simple.SimpleFeatureStore;
import org.geotools.feature.DefaultFeatureCollection;
import org.geotools.feature.FeatureCollections;
import org.geotools.feature.simple.SimpleFeatureBuilder;
import org.geotools.feature.simple.SimpleFeatureTypeBuilder;
import org.geotools.referencing.crs.DefaultGeographicCRS;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.feature.simple.SimpleFeatureType;
import org.opengis.referencing.crs.CoordinateReferenceSystem;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.geom.Polygon;

public class ShapeTools
{
	/**
	 * @param filepath
	 */
	public static void write(SimpleFeatureType TYPE, Object[] data,CoordinateReferenceSystem crs,String filepath)
	{
		 try{    
		        //定义属性  
		/*        final SimpleFeatureType TYPE = DataUtilities.createType("watersurface",  
		            "the_geom:Polygon," + // <- the geometry attribute: Polygon type  
		            "NID:String," + // <- a String attribute  
		            "damParamID:String," + // a number attribute  
		            "sectionParamID:String,"+
		            "reservoirName:String"
		        );  */
		        
		        DefaultFeatureCollection collection = new DefaultFeatureCollection();  
		        
		        SimpleFeatureBuilder featureBuilder = new SimpleFeatureBuilder(TYPE);  
		        /*
		         * 
		        GeometryFactory geometryFactory = new GeometryFactory();  
		        double latitude = Double.parseDouble("116.123456789");  
		        double longitude = Double.parseDouble("39.120001");  
		        String POIID = "2050003092";  
		        String MESHID = "0";  
		        String OWNER = "340881";  
		  
		        Point point = geometryFactory.createPoint(new Coordinate(longitude, latitude)); */ 
		        
		        
		       // Object[] obj = {polygon, POIID, MESHID, OWNER};  
		        SimpleFeature feature = featureBuilder.buildFeature(null, data);  
		        collection.add(feature);  
		        feature = featureBuilder.buildFeature(null, data);  
		        collection.add(feature);  
		  
		        File newFile = new File(filepath);  
		        ShapefileDataStoreFactory dataStoreFactory = new ShapefileDataStoreFactory();  
		        Map<String, Serializable> params = new HashMap<String, Serializable>();  
		        params.put("url", newFile.toURI().toURL());  
		        params.put("create spatial index", Boolean.FALSE);  
		        ShapefileDataStore newDataStore = (ShapefileDataStore) dataStoreFactory.createNewDataStore(params);  
		        newDataStore.createSchema(TYPE);  
		        //newDataStore.forceSchemaCRS(DefaultGeographicCRS.WGS84);
		        newDataStore.forceSchemaCRS(crs);
		  
		        Transaction transaction = new DefaultTransaction("create");  
		        String typeName = newDataStore.getTypeNames()[0];  
		        SimpleFeatureSource featureSource = newDataStore.getFeatureSource(typeName);  
		  
		        if (featureSource instanceof SimpleFeatureStore) {  
		            SimpleFeatureStore featureStore = (SimpleFeatureStore) featureSource;  
		            featureStore.setTransaction(transaction);  
		            try {  
		                featureStore.addFeatures(collection);  
		                transaction.commit();  
		            } catch (Exception problem) {  
		                problem.printStackTrace();  
		            transaction.rollback();  
		            } finally {  
		                transaction.close();  
		            }  
		        } else {  
		            System.out.println(typeName + " does not support read/write access");  
		        }  
		    } catch (Exception e) {  
		        e.printStackTrace();  
		    }  
	}
}

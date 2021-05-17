var fc_3spc = ee.FeatureCollection("users/rriggs/na_sj_using_R_min_100")
var gauge_q = ee.FeatureCollection("users/rriggs/qSpatial")
var gauge_w = ee.FeatureCollection("users/rriggs/wSpatial")
var missSD = ee.FeatureCollection("users/rriggs/eSpatialCal")

//Read in RivWidthCloud functions. 
var x = require('users/eeProject/RivWidthCloudPaper:functions_Landsat578/functions_waterClassification_Jones2019.js')
var fls = require('users/eeProject/RivWidthCloudPaper:functions_Landsat578/functions_landsat.js');
var flsh = require('users/eeProject/RivWidthCloudPaper:rwc_landsat.js');
var fnsLandsat = require('users/eeProject/RivWidthCloudPaper:functions_Landsat578/functions_landsat.js');
var lsFun = require('users/eeProject/RivWidthCloudPaper:functions_Landsat578/functions_landsat.js');
var riverFun = require('users/eeProject/RivWidthCloudPaper:functions_river.js');
var grwl_cline = ee.FeatureCollection('users/eeProject/GRWL_summaryStats')

var gauges = gauge_q
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions. 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Function to convert line to point. 
var line2pt =  function(f){
        var f = f.set({'lineGeometry': f.geometry()})
        var l = f.geometry().coordinates()
        var g = ee.Geometry.MultiPoint(l, 'EPSG:4326')
        return(f.setGeometry(g))
}


//Function to convert point to line. 
var pt2line = function(f){
        //return(f.setGeometry(f.get('lineGeometry')).set(ee.Dictionary(['lineGeometry'])))
        ////return(ee.Feature(ee.Geometry.LineString(f.geometry().coordinates())))
        
        var ln = ee.Feature(ee.Geometry.LineString(f.geometry().coordinates()))
        var lnc = ln.copyProperties(f)
        return(lnc)
}

//Calculate distance of water pixels from GRWL centerline for filtering out non river pixels. 
riverFun.ExtractChannel = function(image) {
  // extract the channel water bodies from the water mask, based on connectivity to the reference centerline.
  var connectedToCl = image.not().cumulativeCost({
    source: ee.Image().toByte().paint(grwl_cline, 1).and(image), // only use the centerline that overlaps with the water mask
    maxDistance: 4000,
    geodeticDistance: false
  }).eq(0);
  var channel = image.updateMask(connectedToCl).unmask(0).updateMask(image.gte(0)).rename(['channelMask']);
  //Map.addLayer(channel, {}, 'river mask/channel')
  return (channel);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Set up GEE tool layout and text. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
var mapLayers = Map.layers(); 
var grwlLayer = ui.Map.Layer(gauges, {color: 'yellow'}, 'Locations', true, 1);
Map.setOptions('satellite');
Map.style().set({cursor: 'crosshair'});
var instructionTitle = ui.Label('INSTRUCTIONS', {fontWeight: 'bold'});
var textStyle = {margin: '1px 15px'};
var L1 = ui.Label('Click on a river gauge (yellow) to calculate river discharge time series;', textStyle);
var L12 = ui.Label('(Optional): Adjust the time period and cloud cover prior to clicking on the map;', textStyle);
//var L2 = ui.Label('Click on data point in the time series will show the corresponding satellite image;', textStyle);
var L3 = ui.Label('Click on the tabout buttom to save data as csv.', textStyle);
var L4 = ui.Label('Result will appear after 1–3 min.', textStyle);
var citeTitle = ui.Label('Remotely Observed Discharge from Effective width Occurrence (RODEO)', {fontWeight: 'bold'});
var L5 = ui.Label('Riggs, R., Allen, G.H., David, C.H., Lin, P., Pan, M., Yang, X., "RODEO: An algorithm and Google Earth Engine application for river discharge retrieval from Landsat."', {margin: '1px 15px 10px 15px'}, 'https://github.com/Ryan-Riggs/RODEO');
var instructionPanel = ui.Panel([instructionTitle, L1, L12, L3, L4], ui.Panel.Layout.flow('vertical'));
var citationPanel = ui.Panel([citeTitle, L5], ui.Panel.Layout.flow('vertical'));
mapLayers.set(0, grwlLayer);
var controlTitle = ui.Label('ADJUST DEFAULT PARAMETER', {fontWeight: 'bold'});
var startLabel = ui.Label('Start year:', textStyle);
var startYear = ui.Slider(1984, 2021, 2015, 1);
startYear.style().set({minWidth: '300px', margin: '1px 15px'});
var endLabel = ui.Label('End year:', textStyle);
var endYear = ui.Slider(1984, 2021, 2021, 1);
endYear.style().set({minWidth: '300px', margin: '1px 15px 10px 15px'});
var cloudLabel = ui.Label('Cloud cover less than (%):', textStyle)
var cloud = ui.Slider(1, 100,10, 1)
cloud.style().set({minWidth: '300px', margin: '1px 15px 10px 15px'});
var controlPanel = ui.Panel(
  [controlTitle, cloudLabel, cloud, startLabel, startYear, endLabel, endYear], 
  ui.Panel.Layout.flow('vertical'),
  {border: '1px dashed black'});
  
var rawTitle = ui.Label('Raw estimated errors:', {fontWeight:'bold'})
var correctedTitle = ui.Label('Corrected estimated errors:', {fontWeight:'bold'})
var panel = ui.Panel(
  [instructionPanel, citationPanel, controlPanel], 
  ui.Panel.Layout.flow('vertical'), 
  {
    position: 'bottom-right', 
    width: '500px'
  });
var widgetList = panel.widgets();
ui.root.add(panel);

//Click on a location and filter through rating curves by location. 
var test = Map.onClick(function(coords) {
  widgetList.set(6, ui.Label(''));
  var thisPoint = ee.Geometry.Point([coords.lon, coords.lat]);
  var aoi = thisPoint.buffer(10000).bounds();
  var aoiLayer = ui.Map.Layer(ee.Geometry.LinearRing(ee.List(aoi.coordinates().get(0))), {color: 'red'}, 'AOI', true, 1);
  var containFilter = ee.Filter.contains({
    leftField: 'geometry', 
    rightValue: aoi});
  var gaugesROI = gauges.filterBounds(aoi).map(function(f){
    var dist = f.geometry().distance(thisPoint)
    var out = f.set('distance', dist)
    return(out)
  })
  var gauge = gaugesROI.sort('distance').first().get('Sttn_Nm')
  var gauge_filter = gauges.filter(ee.Filter.eq('Sttn_Nm', gauge)) // Changed from 'SITE_NUM

//Find corresponding rating curve. 
  var i = gauge
  var qrc = gauge_q.filter(ee.Filter.eq('Sttn_Nm', i)).first()
  var wrc = gauge_w.filter(ee.Filter.eq('Sttn_Nm', i)).first()
  var stats = missSD.filter(ee.Filter.eq('Sttn_Nm', i)).first()

//Convert rating curve tables to only width/discharge columns. 
var qrc_names = qrc.propertyNames().filter(ee.Filter.stringStartsWith('item','Q'))
var wrc_names = wrc.propertyNames().filter(ee.Filter.stringStartsWith('item','W'))
var qrc_names = qrc_names.filter(ee.Filter.neq('item', 'Q_50'))
var qrcFilter = ee.Feature(qrc).select(qrc_names)
var wrcFilter = ee.Feature(wrc).select(wrc_names)
var qrcFilter = qrcFilter.toDictionary()
var wrcFilter = wrcFilter.toDictionary()
//Convert Width/discharge value into arrays.
var wArray = wrcFilter.toArray().sort()
var qArray = qrcFilter.toArray().sort()
var wList = wArray.toList()
var wList = wList.reverse()
var qArray = qArray.toList().reverse()
var qArray = ee.Array(qArray)
var wDistinct = wList.distinct()
var qList = wDistinct.map(function(f){
  var b = wList.indexOf(f)
  return(qArray.toList().getNumber(b))
})
var wArray = ee.Array(wDistinct)
var qArray = ee.Array(qList)
var bias = ee.Array(ee.Number(stats.get('Bias')).abs())
var stde = ee.Array(stats.get('STDE'))
var rmse = ee.Array(stats.get('RMSE'))
var data = ee.List([bias, stde, rmse]);
var labs = ee.List(['Absolute Bias', 'STDE', 'RMSE'])
var statsChart = ui.Chart.array.values(data, 0, labs)
  .setChartType('BarChart').setOptions({legend: {position: 'none'}, hAxis:{'title':'Discharge(cms)'}});

//Set values for Rating curve chart. 
var arrayChart = ui.Chart.array.values(qArray, 0, wArray).setChartType('LineChart')
  arrayChart.setSeriesNames(['Rating Curve (1984-2013)']).setOptions({title: 'Rating Curve', hAxis: {'title': 'Width (m)'}, 
  vAxis: {'title': 'Discharge (cms)'},
  lineWidth: 5,
  pointSize:0, series: [{color: 'black'}]
})

//Determine nearest grwl cross section to point location. 
  var distance_fun = function(rr){
  var l = rr.geometry().distance(gauge_filter.geometry())
  var f = rr.set('distance', l)
  return(f)
}

//Create a 2 km buffer around the point location. 
  var buff1 = gauge_filter.map(function(f) {
  return f.buffer(2000)
})

var fi1 = buff1.map(function(f) {
  return ee.FeatureCollection(f)
})

var filt1 = fc_3spc.filterBounds(fi1.flatten().geometry())
var fff = filt1.map(distance_fun)
var filt = fff.filter(ee.Filter.lte('distance',1000))

//Closest xsection to center point. 
var fc_id = filt.limit(1,'distance')

//Filter grwl centerline to region of interest to limit computation time. 
var grwl_filt = grwl_cline.filterBounds(filt)

//Create a GRWL buffer based on grwl mean width. 
var grwl_filt_geom = grwl_filt.geometry().buffer(ee.Number(fc_id.first().get('width_m')).multiply(1.5))

  var buff1 = gauge_filter.map(function(f) {
  return f.buffer(ee.Number(fc_id.first().get('width_m')).multiply(3))
})
//Intersecting geometry of the GRWL buffer with the 2 km buffer around the point location. 
var intersect = grwl_filt_geom.intersection(buff1)
var intersectBuffer = intersect.buffer(30).intersection(buff1)
var buff1 = intersectBuffer

//Add map of flagging buffer. 
Map.addLayer(buff1, {color:'blue'}, 'Flagging buffer')
var polygon = intersect

//Add map of reach location geometry. 
Map.centerObject(polygon);
Map.addLayer(polygon, {}, 'Reach location')

/*
//Determine the GRADES flowlines found within the intersecting geometry. 
var distinct_comid = filt1.filterBounds(polygon).sort('distance').distinct('COMID')
//Map.addLayer(distinct_comid, {}, 'xsections within polygon', false)
var distinct_comid = distinct_comid.reduceColumns(ee.Reducer.toList(), ['COMID'])

//Create a dictionary of unique COMID. 
var distinct_length = ee.List(distinct_comid.get('list'))
var numb_comid = distinct_length.size()
var seq = ee.List.sequence(1, numb_comid)
var comid_fun = function(f){
  var a = ee.String('COMID').cat(ee.String(f)).slice(0,6)
  return(a)
}
var comid_seq = seq.map(comid_fun)
var comid_dict = ee.Dictionary.fromLists(comid_seq, distinct_length)
*/

//Function to clip to the 2 km point location buffer. 
var connected = function(f){
  var con = f.clip(buff1)
  return(con)
}

//Function to clip to the intesecting geometry. 
var grwl_filter = function(f){
  var con = f.clip(polygon)
  return(con)
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Landsat filtering. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Match common names to the sensor-specific bands:
var LT5_BANDS = ['B1',   'B2',    'B3',  'B4',  'B5',    'B7',    'B6', 'pixel_qa'];
var LE7_BANDS = ['B1',   'B2',    'B3',  'B4',  'B5',    'B7',    'B6', 'pixel_qa'];
var LC8_BANDS = ['B2',   'B3',    'B4',  'B5',  'B6',    'B7',    'B10', 'pixel_qa'];
var STD_NAMES = ['Blue', 'Green', 'Red', 'Nir', 'Swir1', 'Swir2', 'Temp', 'BQA'];

// load Landsat 5,7,8 collections:
// TODO(GHA): combine 5, 7, and 8 collections:
var LT5 = ee.ImageCollection('LANDSAT/LT5_L1T_SR')
    .select(LT5_BANDS, STD_NAMES);
var LT5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR')
    .select(LT5_BANDS, STD_NAMES); 
var LE7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR')
.select(LE7_BANDS, STD_NAMES).filterDate('1999-01-01', '2003-05-30');
var LC8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR')
   .select(LC8_BANDS, STD_NAMES);
var collection = LC8.merge(LT5).merge(LE7)
var filtered = collection
    .filterDate(startYear.getValue() + '-01-01', endYear.getValue() + '-12-31')
    .sort('system:time_start', false)
    .filterBounds(filt);
var ls = filtered

//Clip all landat images by 2 km point location region. 
var filtered = filtered.map(connected)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Following RivWidthCloud steps to classify clouds, water and river pixels. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Add Fmask values to the images for cloud filtering. 
var Unpack = function(qualityBand, startingBit, bitWidth) {
  // unpacking bit information from the quality band given the bit position and width
  // see: https://groups.google.com/forum/#!starred/google-earth-engine-developers/iSV4LwzIW7A
  return(qualityBand
  .rightShift(startingBit)
  .bitwiseAnd(ee.Number(2).pow(bitWidth).subtract(1).int()));
};
var UnpackAllSR = function(bitBand) {
  // apply Unpack function for multiple pixel qualities
  var bitInfoSR = {
    'Cloud': [5, 1],
    'CloudShadow': [3, 1], 
    'SnowIce': [4, 1],
    'Water': [2, 1]
  };
  var unpackedImage = ee.Image();
  for (var key in bitInfoSR) {
    unpackedImage = ee.Image.cat(
      unpackedImage, Unpack(bitBand, bitInfoSR[key][0], bitInfoSR[key][1])
      .rename(key));
  }
  return(unpackedImage.select(Object.keys(bitInfoSR)));
};
var AddFmaskSR = function(image) {
  // add fmask as a separate band to the input image
  var temp = UnpackAllSR(image.select(['BQA']));

  // construct the fmask
  var fmask = (temp.select(['Water']).rename(['fmask'])
  .where(temp.select(['SnowIce']), ee.Image(3))
  .where(temp.select(['CloudShadow']), ee.Image(2))
  .where(temp.select(['Cloud']), ee.Image(4)))
  .mask(temp.select(['Cloud']).gte(0)); // mask the fmask so that it has the same footprint as the quality (BQA) band

  return(image.addBands(fmask));
};
var flagged = filtered.map(AddFmaskSR)

//Clipped Landsat scenes with Fmask values attached. 
var filtered = flagged

//Clip Grwl centerline to the exact length of the 2 km buffer. This allows for a consistent length in calculating the effective widths. 
var cline_updated = grwl_filt.geometry().intersection(polygon)

//Determine the Fmask value of landsat scenes along the grwl centerline. 
var cloudFunction = function(f){
  var cld = f.select('fmask').gt(2)
  var min = cld.mask(cld)
  var min = min.reduceRegion(ee.Reducer.sum().unweighted(), polygon) //polygon
  var area = polygon.area()
  var div= (ee.Number(min.get('fmask')).multiply(900)).divide(ee.Number(area))
  var add = f.set({'cloud': div.multiply(100)})
  return(add)
}

//Apply the above function and filter out any scenes with clouds or cloudshadows along the centerline. 
var filtered = filtered.map(cloudFunction)//.filterMetadata('cloud', 'less_than', 2)
var cldFilter = filtered.filter(ee.Filter.lte('cloud', cloud.getValue()))
var filtered = cldFilter

//Apply the DSWE water classification on all clipped Landsat scenes. 
var waterJones = function(f){
  var img = x.ClassifyWaterJones2019(f)
  var out = img.copyProperties(f)
  return(out.set({'system:time_start': f.get('system:time_start')}))
}
var waterMask = filtered.map(waterJones)//(x.ClassifyWaterJones2019)

//Determine connectivity of water pixels to grwl centerline and filter out any water pixels more than 4km away from centerline. 
riverFun.GetCenterline = function(clDataset, bound) {
  // filter the GRWL centerline based on area of interest
  var cl = clDataset.filterBounds(bound); 
  return(cl);
};
riverFun.GetCenterline(grwl_cline, filt)
var riverMask = waterMask.map(riverFun.ExtractChannel)

//Function to mask out non river pixels. 
var MaskFunction = function(f){
  var a = f.eq(1)
  var b = f.mask(a)
  return(b)
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Calculate river widths of inner GRWL buffer and flagging buffer. 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

var connected_mask = riverMask
//Clip the cirlce mask to the intersecting geometry. 
var connected_filter = connected_mask.map(grwl_filter)
//Map.addLayer(connected_filter.first(), {}, 'connected mask grwl filter', false)
var connected_test = connected_filter

//Function to determine river length. 
var effective_width_1=function (f){
  var lng = f.reduceRegions({
  collection: grwl_cline.filterBounds(filt), 
  reducer: ee.Reducer.count(),
  })
  return(lng) 
}
var eff_width1 = connected_test.map(effective_width_1)

//Function to count all water pixels. 
var effective_width_2=function(con){
  var sum = con.select('channelMask').eq(1).reduceRegions({
  collection: buff1, 
  reducer: ee.Reducer.sum().unweighted(),
  })
  return(sum.copyProperties(con).set({'system:time_start': con.get('system:time_start')}))
}

//Mask out non water pixels after length has been calculated. 
var connected_test = connected_test.map(MaskFunction)
var connected_mask = connected_mask.map(MaskFunction)

//Determine amount of river pixels in entire 2 km buffer. 
var eff_width2 = connected_test.map(effective_width_2)

//Determine amount of river pixels in intersecting geometry. 
var eff_width2_filt = connected_mask.map(effective_width_2)

//Function to convert number of river pixels into river area. 
var Area_fun = function (feature) {
  var d = ee.Number(feature.get('sum'));
  var id = feature.getString('system:index').slice(0,24)
  return feature.set({Area: ee.Number(d).multiply(900), id: id, 'system:time_start':feature.get('system:time_start'), 'length': cline_updated.length()});
}
Map.addLayer(cline_updated, {color:'yellow'}, 'GRWL centerline')
//Determine area of river in 2 km buffer. 
var eff_width2 = eff_width2.map(function(f){
  var ft = ee.FeatureCollection(f).first()
  return (ee.Feature(ft).copyProperties(f).set({'system:time_start': f.get('system:time_start')}))
})
var area_map = eff_width2.map(Area_fun)
//Join the 2 km buffer analysis and intersecting geometry analysis to determine if the intersecting geometry is missing river pixels. 
var toyFilter = ee.Filter.equals({
  leftField: 'system:time_start',
  rightField: 'system:time_start'
});

// Define the join.
var innerJoin = ee.Join.inner()

//Determine area of river in intersecting geometry. 
var eff_width2_filt = eff_width2_filt.map(function(f){
  var ft = ee.FeatureCollection(f).first()
  return (ee.Feature(ft).copyProperties(f).set({'system:time_start': f.get('system:time_start')}))
})
var circle = eff_width2_filt.map(Area_fun)

var fc_function = function(f){
  var p = ee.Feature(f.get('primary'))
  var s = ee.Feature(f.get('secondary'))
  var properties = p.copyProperties(s)
  return(properties)
}

var fc_test = area_map
var difference = innerJoin.apply(circle, fc_test, toyFilter)

//Function to determine the different amount of river pixels between the two geometries.  
var difference_function = function(f){
  var p = ee.Feature(f.get('primary'))
  var a = ee.Number(p.get('Area'))
  var s = ee.Feature(f.get('secondary'))
  var b = ee.Number(s.get('Area'))
  return s.set({Difference: a.subtract(b)})
}
var diff_t = difference.map(difference_function)
var fc_test = diff_t

//Function to calculate effective width. 
var effW_fun = function (feature) {
  var Area = ee.Number(feature.get('Area'))
  var Length = ee.Number(feature.get('length'))
  return feature.set({Effective_width: Area.divide(Length)});
}
var testing = fc_test.map(effW_fun)
var output = testing

//Add in the closest xsection ID so that all locations are unique. 
var fields = function(f){
  var a = f.set({ID: fc_id.first().get('ID_2')})
  return(a)
  }
var output1 = output.map(fields)
var selection = output1.select(['Difference', 'Effective_width', 'ID', 'id', 'COMID', 'width_m', 'change', 'system:time_start', 'cloud'])
var final = selection

var composeprop = function(f) {
  var x = f
  return(ee.Image(x))
}
var matches = final.map(composeprop)

var final_filtered = matches.filter(ee.Filter.eq('Difference', 0))
var final_filtered = final_filtered.filter(ee.Filter.gte('Effective_width', wArray.toList().sort().get(0)))
var final_filtered = final_filtered.randomColumn()
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create charts and add results to panels. 
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Create chart layout for effective width timeserires. 
var chart = ui.Chart.feature.byFeature({
    features: final_filtered, 
    xProperty: 'system:time_start', 
    yProperties: 'Effective_width'});
  
  chart.setOptions({
    title: '', 
    hAxis: {title: 'Date'},
    vAxis: {title: 'Width (m)'},
    series: [
    {color: 'black', visibleInLegend: false, pointSize: 3, lineWidth: 0.5}
  ]//,trendlines: {0: {color: 'red'}},
  });
  
var finalArray = final_filtered.aggregate_array('Effective_width')
var img = final_filtered.reduceToImage(['Effective_width'], ee.Reducer.allNonZero())

//Steps to interpolate from our rating curves and allow for width to discharge estimates to be made. 
var floatFun = function(f){
  var f = ee.Number.parse(f)
  var a = ee.Number(f).float()
  return (ee.Image(a))
}

//
var qFunction = function(f){
  var img = f.get('Effective_width')
  var int = floatFun(img).interpolate(wArray.toList().sort(), qArray.toList().sort(), 'mask')
  var qTrend = int.reduceRegion({reducer: ee.Reducer.max().forEachBand(int), geometry:aoi, scale: 30, bestEffort: true})
  var qImg = f.set({'Discharge': qTrend.get('constant')});
  return(qImg)
}
var finalQ = final_filtered.map(qFunction)

//Create discharge time series chart. 
var qChart = ui.Chart.feature.byFeature({
    features: finalQ, 
    xProperty: 'system:time_start', 
    yProperties: ['Discharge']});
  
  qChart.setOptions({
    title: '', 
    hAxis: {title: 'Date'},
    vAxis: {title: 'Discharge (cms)'},
    series: [
    {color: 'black', visibleInLegend: false, pointSize: 3, lineWidth: 0.5},
    {color:'red', visibleInLegend:true, pointSize:0, lineWidth:0.25}
  ],//trendlines: {0: {color: 'red'}},
  });
  
//Assign previously designed charts to panel. 
  var StatsPanel =ui.Panel([ui.Label('Uncertainty', {fontWeight: 'bold'}), statsChart], ui.Panel.Layout.flow('vertical'));

    widgetList
  .set(2, StatsPanel)

  
  var ChartPanel = ui.Panel([ui.Label('Width Time Series', {fontWeight: 'bold'}), chart], ui.Panel.Layout.flow('vertical'));

    widgetList
  .set(4, ChartPanel)



  var arrayPanel = ui.Panel([ui.Label('Rating Curve', {fontWeight: 'bold'}), arrayChart], ui.Panel.Layout.flow('vertical'));
  
  widgetList
  .set(3, arrayPanel)
  
  var ChartPanel = ui.Panel([ui.Label('Discharge Time Series', {fontWeight: 'bold'}), qChart], ui.Panel.Layout.flow('vertical'));

    widgetList
  .set(5, ChartPanel)  

  var parameterPanel = controlPanel

    widgetList
  .set(1, parameterPanel) 

})

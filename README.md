## Contouring

Contouring flat-file time-series Pressure Transducer data and conversion to a web uploadable format for visualization in Kepler.gl
 

## Uses
<b> Python </b>
- Pandas
- Geopandas
- Shapely
- Numpy
- Matplotlib
- Scipy


## Usage
```
from Contour import ContourData

cd = ContourData(
        filepath = 'WellElevationSubsetSmaller.xlsx',
        xcol = 'X',
        ycol = 'Y',
        zcol = "Water Elevation(ft NAVD88)",
        tcol = "DateTimeRounded",
        crs = 2926,
        method = 'linear',
        steps = 50
)

cd.subsetData()
cd.contourData(subdivisions=5)
cd.extractGeometry()
cd.mergeGeometry()
cd.reproject()
cd.toGeoJSON("test.json")
```
## Output Looks like

![Alt Text](contourmpl.png)

## In Kepler.gl

![Alt Text](Full-Year.gif)




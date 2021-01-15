
import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from shapely import geometry as geo
import numpy as np
from scipy.interpolate import griddata

class ContourData():

    """
    Contour Creation Data class

    Attributes
    ----------
    xcol : str
        x-coordinate column header
    ycol : str
        y-coordinate column header
    zcol : str
        z-coordinate column header
    tcol : str
        time stamp column header
    crs : int
        coordinate reference system of x,y coordinate data  - EPSG code
            ex. WGS84 - Lat/Long = 4326 
     method : str
        surface interpolation/creation method:
            - cubic (default)
            - linear
            - triangulation        
    steps : int
        surface interpolation grid resolution in crs units 
            - 2 (default)    

    """

    def __init__(self,filepath,xcol,ycol,zcol,tcol,crs,method='cubic',steps=2):
        
        self.filePath = filepath
        self.xcol = xcol
        self.ycol = ycol
        self.zcol = zcol
        self.tcol = tcol
        self.crs = crs
        self.data = pd.read_excel(filepath)
        self.method = method
        self.steps = steps

    def subsetData(self):
        '''
        Subsets the dataframe by each unique timestamp

            Returns:
                    self.dtTime (df): 
                    Pandas dataframe of the form
                    DateTime | DataFrame
                        DateTime (dt) : datetime of the contours
                        DataFrame (df) : DataFrame of orginal data at DateTime
        '''
        
        times = self.data[self.tcol].unique()
        dfsubs = []
        
        for t in times:
            dft = self.data[self.data[self.tcol] == t]
            dfsubs.append(dft)
        
        self.dfTimes = pd.DataFrame(
                        data = [times,dfsubs]
                        )
        self.dfTimes = self.dfTimes.T
        self.dfTimes.columns = ['DateTime','DataFrames']
    
    def interpSurface(self,steps=2,method='cubic'):
        '''
        Interpolates a surface from the sparse points

            Parameters:
                    type (str) : interpolation method, options include:
                                - linear
                                - cubic

            Returns:
                    self.dfSurfaces (dataframe):
                    Pandas dataframe of the form
                    DateTime | Surface
                        DateTime (dt) : datetime of the contours
                        Surface (dataframe) : dataframe of x,y,z coordinates at DateTime
        '''

        surfaces = []
        
        for index, row in self.dfTimes.iterrows():
            
            x = row['DataFrames'][self.xcol].values
            y = row['DataFrames'][self.ycol].values
            z = row['DataFrames'][self.zcol].values

            points = np.array([x,y]).T
            xmin = x.min()
            xmax = x.max()
            ymin = y.min()
            ymax = y.max()
            grid_x, grid_y = np.meshgrid(np.arange(xmin,xmax,steps), np.arange(ymin,ymax,steps))
            grid_z = griddata(points, z, (grid_x, grid_y), method=method)

            # remove non-finite values (not interpolated)
            mask = np.isfinite(grid_z.flatten())

            gfx = grid_x.flatten()
            gfy = grid_y.flatten()
            gfz = grid_z.flatten()

            grid_xn = gfx[mask]
            grid_yn = gfy[mask]
            grid_zn = gfz[mask]
            xyz = pd.DataFrame([grid_xn,grid_yn,grid_zn])
            xyz = xyz.T
            xyz.columns = [self.xcol,self.ycol,self.zcol]
            surfaces.append(xyz)

        self.dfTimes['Surfaces'] = surfaces
            
        
    def contourIntervals(self,zvalues,delta=0.25):
        '''
        Defines the intervals at which to create contours given the z-values and the interval delta 

            Parameters:
                    zvalues (array): A decimal integer
                    delta (float): contour interval delta 

            Returns:
                    vals (array): array of values betwee z_max and z_min spaced by delta units
        '''
        
        min_ = np.min(zvalues)
        max_ = np.max(zvalues)
        vals = np.arange(min_,max_,delta)
    
        return vals
        
    def contourData(self,interval=None,subdivisions=5):
        
        '''
        Contours each self.dfTimes dataframe using using matplotlib tricontour at a given interval and subdivisions

            Parameters:
                    interval (float): contour interval delta 
                    subdivisions (int): number of subdivided triangulations to compute before contouring

            Returns:
                    self.dfContours (dataframe): DataFrame of the form
                    DateTime | DataFrame | Contours
                    
                        DateTime (datetime) : datetime of the contours
                        DataFrame (dataframe) : DataFrame of orginal data at DateTime
                        Contours (cs obj) : matplotlib collection of the computer contours at DateTime
        '''

        if self.method in ['cubic','linear']:

            self.interpSurface(steps= self.steps, method=self.method)
            col = 'Surfaces'

        else:
            col = 'DataFrames'

        csc = []
        
        for index, row in self.dfTimes.iterrows():
            
            x = row[col][self.xcol].values
            y = row[col][self.ycol].values
            z = row[col][self.zcol].values
            
            if interval != None:
                vals = self.contourIntervals(z,interval)
            else:
                vals = None
                
            triang = tri.Triangulation(x, y)
            refiner = tri.UniformTriRefiner(triang)
            tri_refi, z_test_refi = refiner.refine_field(z, subdiv=subdivisions)
            
            cs = plt.tricontour(tri_refi, z_test_refi, levels=vals)
            csc.append(cs)
        
        self.dfContours = self.dfTimes.copy()
        self.dfContours['CS'] = csc
    
    def dateTimeForm(self,dt):
            
        '''
        Converts DateTime format to string and formats for use in GeoJSON

            Parameters:
                    dt (datetime): contour interval delta 

            Returns:
                    dtjson (str): str of the form "yyyy-mm-ddThh:mm:ss"

        '''
        
        dt = str(dt)
        dtjson = dt.replace(" ","T")
        
        return dtjson
        
    def extractGeometry(self):
    
        '''
        Converts the matplotlib contour objects into shapely LineStrings and returns a GeoDataFrame of 
        Contours at each time


            Returns:
                    self.dfContours (dataframe): DataFrame of the form
                    
                    DateTime | Contour | gdf
                    
                        DateTime (datetime) : datetime of the contours
                        Contours (cs obj) : matplotlib collection of the computer contours at DateTime
                        gdf (geodataframe) : GeoDataFrame of contours at DateTime of the form:
                            
                            DateTime | Contour | geometry
                            
                            DateTime (datetime) : datetime of the contours
                            Contour (float) : value of the contour
                            geometry : shapely LineString geometry of the object.
                            
            
        '''
        
        cntrdfs = []
        
        for index, row in self.dfContours.iterrows():
            
            cs = row['CS']
            linestrings = []
            times = []
            for clc in cs.collections:
                paths = clc.get_paths()
                
                if len(paths) > 0:
                    path = paths[0].vertices
                    if len(path) > 1:
                        ls = geo.LineString(path)
                        linestrings.append(ls)
              
            times = [row['DateTime'] for i in range(len(linestrings))]
            lvls = cs.levels[0:len(linestrings)]
            cntrdf = gpd.GeoDataFrame(
                        data = [lvls,times,linestrings],
                        )
            cntrdf = cntrdf.T
            cntrdf.columns = ['Contour','DateTime','geometry']
            cntrdfs.append(cntrdf)
            cntrdf['DateTime'] = pd.to_datetime(cntrdf['DateTime'])
            cntrdf['DateTime'] = [self.dateTimeForm(dt) for dt in cntrdf['DateTime']]
        self.dfContours['gdf'] = cntrdfs
    
    def mergeGeometry(self):
    
        '''
        Merges all self.dfContours into a single geodataframe

            Returns:
                    self.MergedGeo (GeoDataFrame): GeoDataFrame of the form"
                    
                         DateTime | Contour | geometry
                            
                            DateTime (datetime) : datetime of the contour
                            Contour (float) : value of the contour
                            geometry : shapely LineString geometry of the contour.

        '''
        
        gdfs = self.dfContours['gdf'].tolist()
        geos = pd.concat(gdfs)
        geo = gpd.GeoDataFrame(geos)
        geo.crs = geo.crs = self.crs
        self.MergedGeo = geo 
        
    def reproject(self,crsNew=4326):
         
       
        '''
        Reprojects the data to a new coordinate reference system

            Parameters:
                    crsNew (int) : EPSG code of the new CRS

            Returns:
                    self.MergedGeo (GeoDataFrame): GeoDataFrame of the form"

                         DateTime | Contour | geometry

                            DateTime (datetime) : datetime of the contour
                            Contour (float) : value of the contour
                            geometry : shapely LineString geometry of the contour.
        '''
        self.MergedGeo = self.MergedGeo.to_crs(crsNew)
            
        

    def toGeoJSON(self,filepath):
        

        '''
        Exports the file to GEOJSON

            Parameters:
                    filepath (str) : filepath of the output File

        '''
        self.MergedGeo.to_file(filepath)
    
  
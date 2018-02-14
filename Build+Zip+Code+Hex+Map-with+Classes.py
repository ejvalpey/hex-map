
# coding: utf-8

# In[1]:

### Import required packages

#system packages
import os   # operating system library
import odbc # for open database connectivity connections
import gc   # garbage collector

#data science packages
import pandas as pd  # best for working with dataframes
import numpy as np
from matplotlib import pyplot

import itertools
from scipy.spatial import distance_matrix
from operator import itemgetter

### Set Parameters

#change the folder we're working with
os.chdir('C:\Users\us44532\Documents\Python Scripts')

chunkSize = 50000
states = ['WA'] # to be updated with more states and to be looped through in _main_

### Connect to Netezza

cnxn = odbc.odbc('DSN=CUE') # depends on your named connections in your ODBC Manager



# FACTOR = 0.3


# In[2]:

# set Header Query
def setHDRSQL(state, cursor, chunkSize):
    sqlHDRString = """
    SELECT   CITY
           , ZIP_CODE
           , LONGITUDE
           , LATITUDE
      FROM IDS_LAB..ZIP_LAT_LONGS
     WHERE
         STATE_CD = '{}'
         AND ZIP_CODE > {}
     LIMIT {}""".format(state, str(cursor), chunkSize)
    return sqlHDRString


# In[3]:

def readStateZips(state):
    cursor = 0
    chunkID = 1
    continueQuery = True
    first = True
    while continueQuery:
        try:
            sqlHDRString = setHDRSQL(state, cursor, chunkSize)
            #print (sqlHDRString)
            tempDf = pd.read_sql(sqlHDRString,cnxn)
            continueQuery = tempDf.shape[0]>0
            cursor = tempDf.ZIP_CODE.max()
            if first:
                df = tempDf
                first = False
            else:
                df.append(tempDf)
                tempDf = None
            return df
        except:
            print ('error in read_sql')
            continueQuery = False
            return None


# In[60]:

def normalize(arr):
    arr = arr+(np.random.randn(len(arr))-0.5)/1000000
    return arr #(arr-arr.min())/(arr.max()-arr.min())


# In[61]:

try:
    # read directly from the flat file if it exists in our directory
    df = pd.read_csv('raw_zip_codes.csv')
    df.drop([df.columns[0]], axis = 1, inplace = True)
    print('Read from .csv')
except:
    # otherwise execute a SQL script and read from the database    
    state = states[0]
    df = readStateZips(state)
    print('Read from database')
df = pd.DataFrame(zip(df.CITY,df.ZIP_CODE,df.LONGITUDE, df.LATITUDE))
df.columns = (['CITY','ZIP_CODE','LONGITUDE','LATITUDE'])


# In[62]:

e = np.array(zip(df.CITY #cityEncode.transform(df.CITY)
                 , df.ZIP_CODE #zipEncode.transform(df.ZIP_CODE)
                 , normalize(df.LONGITUDE)
                 , normalize(df.LATITUDE)
                ))
# %matplotlib inline
e


# ## Define Functions

# In[63]:

class GeoHex(object):
    """ A GeoHex is a unit mapping spatial coordinates (primarily in longitude/latitude terms)
    to a Hexagonal grid while carrying a set of hierarchical properties 
    
    Attributes:
        hierarchy of categorical variables (such as country, state, county, city, zip_code)
        longitude, latitude: float coordinates
        distance from some point in long/lat space 
        x, y coordinates in hex-grid space 
        vacant: boolvacancy (a GeoHex might be initialized as empty)
        placement order
        """
    neighbors = set([(4,0),(2,3),(-2,3),(-4,0),(-2,-3),(2,-3)]) # shared class variable 
    
    def __init__(self, limbo = True, vacant_grid = True, name = ''):
        self.name = name
        self.hier = []
        self.position = (None, None) # initializes nowhere with limbo being True
        self.vacant_grid = vacant_grid
        self.grid_xy = (None, None) # initializes vacant
        self.limbo = limbo
        self.local_neighbors = set() # set of tuples
        self.cursor_distance = None
        """Return a GeoHex object """
        
    def __repr__(self):
        return str(self.__dict__)
        
    def add_hier(self, hierarchy):
        for h in hierarchy:
            self.hier.append(h)
            
    def set_position(self, position):
        self.position = position
        self.limbo = False
        if self.grid_xy != (None, None):
            vacant_grid = False
    
    def set_grid(self, grid_xy):
        self.grid_xy = grid_xy
        self.vacant_grid = False
        if self.position != (None, None):
            limbo = False
        self.local_neighbors = set({(self.grid_xy[0] + neighbor[0], 
                                     self.grid_xy[1] + neighbor[1]) 
                                    for neighbor in GeoHex.neighbors})


# In[64]:

class HexNeighborhood(object):
    """ A HexNeighborhood object is a collection of GeoHex objects organized in a hexagonal grid
    
    Attributes:
        name:    name of the neighborhood
        hex_count:    the number of placed hexagons in the neighborhood
        map_center:    the weighted average lat/long of the placed hexagons
        hex_center:    the weighted average of the placed hex coordinates
        placed:   set of placed hexagons (GeoHex objects)
        unplaced:   set of unplaced hexagons (GeoHex objects)
        available_neighbors:   set of available neighbors (tuples of x, y coordinates)
        stop_loop:   an internal status flag for debugging loops
        position_array:   a distance matrix used to measure density (how to know where to start placement)
    
    Global Attributes:
        factor: a multiplier used to tune the relationship between lat/long distance and hex grid distance
        epsilon: a positive non-zero potentially helpful in avoiding div/0 errors
        
    Methods:
        
     """
    factor = 0.9
    epsilon = 0.00000005
    
    def __init__(self, name):
        self.name = name
        self.hex_count = 0
        self.map_center = (0.0,0.0)
        self.hex_center = (0.0,0.0)
        self.placed = set() #set of GeoHex objects
        self.unplaced = set() #set of GeoHex objects
        self.available_neighbors = set() #set of tuples
        self.stop_loop = False
        self.position_array = np.empty(0)
    
#     def __repr__(self):
# #         print(self.__dict__)
#         pass
        
    def add(self, geohex):
        """adds GeoHex object to the set of unplaced items"""
#         if not geohex.limbo:
#             self.map_center = tuple([(self.hex_count * self.map_center[i] + geohex.position[i])
#                                      / (self.hex_count + 1) for i in range(2)])
#         if not geohex.vacant_grid:
#             self.hex_center = tuple([(self.hex_count * self.map_center[i] + geohex.grid_xy[i])
#                                      / (self.hex_count + 1) for i in range(2)])
#         self.hex_count += 1
        self.unplaced.add(geohex)
    
    def set_dist_matrix(self): # the distance matrix is used to measure density
        """establishes a distance matrix between all added nodes"""
        self.position_array = np.vstack({node.position for node in self.unplaced})
        self.dist_matrix = distance_matrix(self.position_array, self.position_array) # distance_matrix is imported from scipy.spatial
        return self.dist_matrix
    
    @staticmethod
    def inv_sum(row, epsilon = 0.005): 
        return np.mean([1/(x+epsilon) for x in row])
        
    def start_node(self, epsilon = 0.04):
        try:
            density = np.apply_along_axis(HexNeighborhood.inv_sum, 0, self.dist_matrix, epsilon)
            highest = self.position_array[np.argmax(density)]
            return [geohex for geohex in self.unplaced if geohex.position == (highest[0], highest[1])][0]
        except:
            if not self.stop_loop:
                self.set_dist_matrix()
                self.stop_loop = True
                return self.start_node()
            else:
                return None
            
    def update_available_neighbors(self, geohex):
        new_neighbors = geohex.local_neighbors
        for neighbor in new_neighbors:
            self.available_neighbors.add(neighbor)
        for preplaced in new_neighbors.intersection(set(node.grid_xy for node in self.placed)):
            self.available_neighbors.remove(preplaced)
        return self.available_neighbors
    
    def update_squared_dist(self):
        for node in self.unplaced:
            node.cursor_distance = HexNeighborhood.squared_distance(node.position, self.map_center)
        
    def set_start(self): # run once after nodes are loaded into unplaced set
        start_node = self.start_node()
        start_node.set_grid((0,0))
        start_node.cursor_distance = 0.0
        self.map_center = start_node.position
        self.placed.add(start_node)
        self.unplaced.remove(start_node)
        self.update_squared_dist()
        self.hex_count += 1
        return start_node
    
    def min_node(self, cursor_hier):
        x = 1e10
        outnode = None
        num_updates = 0
        for node in self.unplaced:
            if node.cursor_distance < x and np.all([node.hier[i]==cursor_hier[i] for i in range(len(cursor_hier))]):
#                 print(cursor_hier, [node.hier[i]==cursor_hier[i] for i in range(len(cursor_hier))])
                x = node.cursor_distance
                num_updates += 1
                outnode = node
#         print num_updates
        return outnode
    
    def squared_distances_to_geo_center(self, hier):
        counter = 0
        for node in self.unplaced: 
            if np.all([node.hier[i] == hier[i] for i in range(len(hier))]):
                node.cursor_distance = HexNeighborhood.squared_distance(self.map_center, node.position)
                counter += 1
        sorted(self.unplaced, key = self.dist)
        return counter
        
    def get_close_hex(self, geohex):
        hex_space = [HexNeighborhood.factor*4*(geohex.position[i] - self.map_center[i])
                     /(HexNeighborhood.epsilon + np.sqrt(geohex.cursor_distance/self.hex_count)) 
                     for i in range(2)]
#         print(geohex.cursor_distance, self.hex_count, hex_space)
        return min([(node, HexNeighborhood.squared_distance(node, hex_space))
                    for node in self.available_neighbors], key = itemgetter(1))#[1][1:]
    
    @staticmethod
    def squared_distance(node1, node2):
        return sum([(node1[i]-node2[i])**2 for i in range(len(node1))])
    
    def dist(self, x):
        yield [node.cursor_distance for node in self.unplaced]
        
    def place(self, geohex):
#         print self.available_neighbors
        geohex.set_grid(self.get_close_hex(geohex)[0])
        self.placed.add(geohex)
        self.unplaced.remove(geohex)
        self.update_available_neighbors(geohex) 
        self.hex_count += 1
        self.map_center = tuple([(self.hex_count * self.map_center[i] + geohex.position[i])
                                / (self.hex_count + 1.0) for i in range(2)])
        self.hex_center = tuple([(self.hex_count * self.hex_center[i] + np.float(geohex.grid_xy[i]))
                                / (self.hex_count + 1.0) for i in range(2)])
#         print self.available_neighbors

#     def plotUnplaced(self):
#         arr = np.array([[t.position[0], t.position[1]] for t in self.unplaced])
#         pyplot.scatter(arr[:100,0], arr[:100,1])
#         #Look, a map of Washington!
#     def plotPlaced(self):
#         arrXy = np.array([[t.position[0],t.position[1]] for t in self.placed])
#         pyplot.scatter(arrXy[:100,0],arrXy[:100,1])
#         #Look, a scatter map of placed in Washington!
#     def plotHexPlaced(self):
#         arr2 = np.array([[t.grid_xy[0],t.grid_xy[1]] for t in self.placed])
#         pyplot.scatter(arr2[:100,0],arr2[:100,0])
#         #Look, a scatter map of placed in Washington!


# In[65]:

test_neighborhood = HexNeighborhood('Washington Zip Codes')
for i in range(len(e)-1):
    test = GeoHex(name = e[i][1])
    test.set_position((np.float(e[i][2]),np.float(e[i][3])))
    test.add_hier(e[i][:1])
    test_neighborhood.add(test)
# test_neighborhood.set_dist_matrix()
start_node = test_neighborhood.set_start()
start_node.set_grid((0,0))
test_neighborhood.update_available_neighbors(start_node)
#, [node.__dict__ for node in test_neighborhood.unplaced]
cursor_hier = start_node.hier
plot = pyplot


# In[66]:

# %matplotlib
for i in range(len(e)):
    distance_count = test_neighborhood.squared_distances_to_geo_center(cursor_hier)
#     print("top", cursor_hier, distance_count)
    loop_counter = 0
    while distance_count == 0 and loop_counter < 5:
        cursor_hier = cursor_hier[:-1]
        distance_count = test_neighborhood.squared_distances_to_geo_center(cursor_hier)
        loop_counter += 1
#     print("bottom", cursor_hier, distance_count)
    next_place = test_neighborhood.min_node(cursor_hier)
#     print("out", next_place.hier)
    if len(test_neighborhood.unplaced) > 0:
        test_neighborhood.place(next_place)
        cursor_hier = next_place.hier
        test_neighborhood.update_squared_dist()
#     pyplot.subplot(1,2,1)
#     test_neighborhood.plotUnplaced()
        #Look, a map of Washington!
#     test_neighborhood.plotPlaced()
        #Look, a scatter map of placed in Washington!
#     pyplot.subplot(1,2,2)
#     test_neighborhood.plotHexPlaced()


# In[67]:

columns = []
for a in ['name', ['hier {}'.format(i+1) for i in range(len(start_node.hier))]
                    , 'longitude', 'latitude'
                    , 'path', 'x_grid', 'y_grid']:
    if type(a) == type(columns):
        for i in range(len(a)):
            columns.append(a[i])
    else:
        columns.append(a)
verticies = [(2,1),(0,2),(-2,1),(-2,-1),(0,-2),(2,-1)]
out_array = [[[c] for c in range(len(columns))] for r in range(6*len(test_neighborhood.placed))]
# out_array = pd.DataFrame(np_array, columns = columns)
row = 0
hier_size = len(start_node.hier)
out_set = test_neighborhood.placed


# In[68]:

row = 0
for node in out_set:
    for k in range(6):
        out_array[row][0] = node.name
        for i in range(hier_size):
            try:
                out_array[row][1 + i] = node.hier[i]
            except:
                pass
        out_array[row][1 + hier_size] = node.position[0]
        out_array[row][2 + hier_size] = node.position[1]
        out_array[row][3 + hier_size] = 1 + k
        out_array[row][4 + hier_size] = node.grid_xy[0] + verticies[k][0]
        out_array[row][5 + hier_size] = node.grid_xy[1] + verticies[k][1]
        row += 1
    if row%50 == 1:
        print(node.__dict__)
#     break
#     out_set.remove(node)
pd.DataFrame(out_array, columns = columns).to_csv('Hex Polygons.csv')


# In[54]:

len(test_neighborhood.unplaced)


# In[86]:

# updateDistances(center,unplaced)
pyplot.subplot(1,2,1)
plotUnplaced(unplaced)
plotPlaced(placed)
pyplot.subplot(1,2,2)
plotHexPlaced(placed)


# In[11]:

df[['ZIP_CODE','CITY', 'LONGITUDE', 'LATITUDE']].to_csv('raw_zip_codes.csv')


# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 10:19:30 2018

@author: mora
"""


# Before you get started, download file from databricks and rename it.
# The file can be accessed with the web browser with the URL looking something
# like this:
# https://community.cloud.databricks.com/files/tables/plot_data.parquet/part-00000-tid-5476989714034886627-381941b3-5284-40e3-82b2-1cde74ace252-20-c000.snappy.parquet?o=4523053627268386
# where the thing at the end is the session ID and needs to be replaced with
# your current session ID

import pandas as pd
from pyproj import Proj, transform
import datashader as ds
import datashader.transfer_functions as tf
from colorcet import fire
import cartopy.io.img_tiles as cimgt
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

# fix issue with pyproj
import os
os.environ['PROJ_LIB'] = ('C:\\Users\\mora\\AppData\\Local\\Continuum\\'
                          'anaconda2\\envs\\crm_env_36\\Lib\\site-packages'
                          '\\pyproj\\data')


# transform data into Mercator coordinate system to be later used by datashader
def toWebMercatorLon(xLon):
    mercator = transform(Proj(init='epsg:4326'), Proj(init='epsg:3857'),
                         xLon,
                         0)
    # longitude first, latitude second.
    return mercator[0]


def toWebMercatorLat(yLat):
    mercator = transform(Proj(init='epsg:4326'), Proj(init='epsg:3857'),
                         0,
                         yLat)
    return mercator[1]


def create_datashader_plot(df, x_col_name, y_col_name):
    '''
    Plot fancy graph with datashader
    '''

    # calculate web mercator coordinates
    df_merc = pd.DataFrame()
    df_merc['Long_clean'] = df['Long_clean'].map(toWebMercatorLon)
    df_merc['Lat_clean'] = df['Lat_clean'].map(toWebMercatorLat)

    # define Canvas
    x_center = (df_merc[y_col_name].max()+df_merc[y_col_name].min())/2
    y_center = (df_merc[x_col_name].max()+df_merc[x_col_name].min())/2

    # definition of limits
    x_half_range = (df_merc[y_col_name].max()-df_merc[y_col_name].min())/2
    y_half_range = (df_merc[x_col_name].max()-df_merc[x_col_name].min())/2
    x_range, y_range = ((x_center - x_half_range, x_center + x_half_range),
                        (y_center-y_half_range, y_center+y_half_range))
    plot_width = 1000
    plot_height = int(plot_width/(x_half_range/y_half_range))

    cvs = ds.Canvas(plot_width=plot_width,
                    plot_height=plot_height,
                    x_range=x_range,
                    y_range=y_range)

    # plotting
    coord_agg = cvs.points(df_merc, y_col_name, x_col_name)
    coord_img = tf.shade(coord_agg, cmap=fire, how='eq_hist')
    coord_img = tf.set_background(coord_img, 'black')
    return coord_img


def plot_on_map(df,
                x_col_name,
                y_col_name,
                sub_n=None,
                sub_m=None,
                sub_i=None):
    '''
    Plot data on map
    '''
    if sub_n is None or sub_m is None or sub_i is None:
        plt.figure(figsize=(10, 10))
        sub_n = 1
        sub_m = 1
        sub_i = 1

    minLat = min(df[y_col_name])
    minLon = min(df[x_col_name])
    maxLat = max(df[y_col_name])
    maxLon = max(df[x_col_name])

    # create a Stamen Terrain instance.
    stamen_terrain = cimgt.StamenTerrain()

    # create a GeoAxes in the tile's projection.
    ax = plt.subplot(sub_n, sub_m, sub_i, projection=stamen_terrain.crs)

    # limit the extent of the map to a small longitude/latitude range.
    ax.set_extent([minLon, maxLon, minLat, maxLat])

    quality = 14
    ax.add_image(stamen_terrain, quality)
    ax.plot(df[x_col_name],
            df[y_col_name],
            'r.',
            transform=ccrs.Geodetic(),
            markersize=6)
    return ax


def load_and_clean_data(name):
    '''
    Load data and clean up lat/lon
    '''
    df = pd.read_parquet(name)

    df.fillna(value=pd.np.nan, inplace=True)
    df['Long_clean'] = df['Long_clean'].map(float)
    df['Lat_clean'] = df['Lat_clean'].map(float)
    return df


###############################################################################
# all crimes

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# read and prepare data
file_name = os.path.join(parent_dir, 'data', 'plot_data.parquet')
coords = load_and_clean_data(file_name)

# plot fancy graph with datashader
img = create_datashader_plot(coords, 'Lat_clean', 'Long_clean')
img

ax = plot_on_map(coords, 'Long_clean', 'Lat_clean')
plt.title('All crimes, whole period')

x_lim = ax.get_xlim()
y_lim = ax.get_ylim()

###############################################################################
# shootings

# read and prepare data
file_name = os.path.join(parent_dir, 'data', 'plot_data_shooting.parquet')
coords = load_and_clean_data(file_name)

years = coords.Year.unique()
plt.figure(figsize=(10, 10))
i = 1
for year in years:
    plt.subplot(2, 2, i)
    ax = plot_on_map(coords[coords.Year == year],
                     'Long_clean',
                     'Lat_clean',
                     sub_n=2,
                     sub_m=2,
                     sub_i=i)
    plt.title('Shootings in the year ' + year)
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    i += 1

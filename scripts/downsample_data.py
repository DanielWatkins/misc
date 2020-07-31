"""This script applies the downsampling operations outlined in the notebook
'Analysis of inversion strength-LTS relationship' to the Arctic stations in the
IGRA derived dataset. Downsampled files are created for all stations with data 
extending to 2020 north of 65 latitude. Downsampling methods used are the following:
1. Interpolation to height grids 50m, 100m, 150m
2. Interpolation to 25 hPa pressure grid
3. Selection and (where needed) interpolation to mandatory levels
4. Iterative downsampling to significant levels.

Right now, it assumes that the data has already been downloaded and saved to CSV files, 
but that could be included as a step.
"""

##### PARAMETERS ######
# Data folder
saveloc = '/Users/watkdani/Documents/Research/IGRA2 Paper/Data/'
dataloc = '/Users/watkdani/Documents/Research/IGRA2 Paper/Data/IGRA_derived_soundings/'
method1_dir = 'height_grid/'
method2_dir = 'pressure_grid/'
method3_dir = 'mandatory_levels/'
method4_dir = 'significant_levels/'

# Grid spacing for method 1
zgrid_dz = [50, 100, 150] 
zmax = 3000 # Max Z above ground level

# Grid spacing for method 2
pgrid_dp = [25]
pmax = 500 # min p in hPa

# Grid for method 3. Pressure levels to retain.
mandatory_levels = [1000, 925, 850, 700, 500]

# Threshold for method 4. Maximum allowed error in interpolated temperature profile.
eps = [0.5, 1.0]
##### IMPORTS ##########
import pandas as pd
import numpy as np
from siphon.simplewebservice.igra2 import IGRAUpperAir



##### HELPER FUNCTIONS ############
def adjust_time(date):
    """If launch time > 22, assign to 00 the next day. If launch time is between 10 and 2, assign to 12."""
    
    if date.hour <= 2:
        return (date - pd.Timedelta(hours=date.hour)).round('1h')
    if date.hour >= 22:
        adj = np.abs(date.hour - 24)
        return (date + pd.Timedelta(hours=adj)).round('1h')
    if np.abs(date.hour - 12) <= 2:
        adj = date.hour - 12
        return (date - pd.Timedelta(hours=adj)).round('1h')
    else:
        return date
    
def interp_fun_z(df, min_z, dz):
    """Interpolates the lowest 3000 meters of data to the grid given in dz"""
    max_z = 3000
    z = np.arange(min_z, max_z + min_z, dz)
    z0 = df.height
    p0 = df.pressure
    t0 = df.temperature
    idx_notnan = ~np.isnan(z0 + p0 + t0)

    z0 = z0[idx_notnan]
    p0 = p0[idx_notnan]
    t0 = t0[idx_notnan]
    
    df = pd.DataFrame({'height': z, 'pressure': np.exp(np.interp(z, z0, np.log(p0))), 
                         'temperature': np.interp(z, z0, t0)})
    return df


def interp_fun_p(df, dp, p1=None):
    """Interpolates data to the grid given in dp"""
    if p1 is None:
        p1 = np.arange(500,1001, dp)
    
    
    #
    
    z0 = df.height.values[::-1]
    p0 = df.pressure.values[::-1]
    t0 = df.temperature.values[::-1]
    
    idx_notnan = ~np.isnan(z0 + p0 + t0)
    z0 = z0[idx_notnan]
    t0 = t0[idx_notnan]
    p0 = p0[idx_notnan]
    
    z1 = np.interp(np.log(p1), np.log(p0), z0)
    t1 = np.interp(np.log(p1), np.log(p0), t0)
    
    
    df = pd.DataFrame({'height': z1[::-1], 'pressure': p1[::-1], 
                         'temperature': t1[::-1]})

    return df


def significant_levels(df, t_thresh):
    """Retains only enough points in the profile to recreate the profile
    with an error of t_thresh. Intended use is to downsample high resolution
    radiosondes for more even comparison of profiles.
    Beginning with just the top and bottom values, the code repeatedly checks
    the value of the departure from the linear interpolation. The point with
    greatest error is added to the significant levels set. 
    """
    import warnings
    warnings.filterwarnings('ignore', category=FutureWarning)

    traw = df.temperature.values
    zraw = df.height.values
    praw = df.pressure.values
    keep_idx = []
    keep_idx.append(0)
    keep_idx.append(len(traw)-1)
    keep_idx.append(np.random.choice(np.arange(len(traw))))
    t0 = np.array(traw[keep_idx])
    z0 = np.array(zraw[keep_idx])
    p0 = np.array(praw[keep_idx])
    thresh = 0.5
    count = 0
    max_iter = 100
    t_error = np.abs(np.interp(zraw, z0, t0) - traw)
    while (np.max(t_error) > thresh) & (count < max_iter):
        t_error = np.abs(np.interp(zraw, z0, t0) - traw)
        max_idx = np.argmax(t_error)

        keep_idx.append(max_idx)
        keep_idx.sort()

        t0 = np.array(traw[keep_idx])
        z0 = np.array(zraw[keep_idx])
        p0 = np.array(praw[keep_idx])
        count += 1
        if count > max_iter:
            break
            
    return pd.DataFrame({'height':z0, 
                         'temperature':t0, 
                         'pressure':p0})


##### MAIN LOOP #######
station_list = pd.read_fwf('/Users/watkdani/Documents/Research/IGRA2 Paper/igra2-station-list.txt',
                          header=None, names=['station_id', 'lat', 'lon', 'elevation', 'name',
                                              'start_year', 'end_year', 'number'])
arctic_stations = station_list.loc[(station_list.lat >= 65) & (station_list.end_year == 2020)]

for station in arctic_stations.station_id.values[:13]:
    print(station)
    df = pd.read_csv(dataloc + station + '-igra2-derived.csv', index_col=False)
    df.drop('Unnamed: 0', axis=1, inplace=True)
    df.date = pd.to_datetime(df.date)
    df.date = [adjust_time(d) for d in df.date]
    df.rename({'calculated_height': 'height'}, axis=1, inplace=True)

    # identify station elevation
    min_z = np.min(df['height'])
    max_z = 3000
    try:
        # Height grid interpolation
        for dz in zgrid_dz:
            df_zgrid = df.groupby('date').apply(interp_fun_z, min_z, dz)
            df_zgrid.reset_index(inplace=True)
            df_zgrid.drop('level_1', axis=1, inplace=True)
            df_zgrid.to_csv(saveloc + method1_dir + station + '_zgrid_' + str(dz) + '-m.csv', index=False)
            del df_zgrid

    # Pressure grid interpolation
        for dp in pgrid_dp:
            df_pgrid = df.groupby('date').apply(interp_fun_p, dp)
            df_pgrid.reset_index(inplace=True)
            df_pgrid.drop('level_1', axis=1, inplace=True)
            df_pgrid.to_csv(saveloc + method2_dir + station + '_pgrid_' + str(dp) + '-hPa.csv', index=False)
            del df_pgrid
    
    # Mandatory levels
        df_mgrid = df.groupby('date').apply(interp_fun_p, dp, [1000, 925, 850, 700, 500])
        df_mgrid.reset_index(inplace=True)
        df_mgrid.drop('level_1', axis=1, inplace=True)
        df_mgrid.to_csv(saveloc + method3_dir + station + '_mlevels_' + '.csv', index=False)
        del df_mgrid
    
    # Significant levels
        for epsilon in eps:
            df_sl = df.groupby('date').apply(significant_levels, epsilon)
            df_sl.reset_index(inplace=True)
            df_sl.drop('level_1', axis=1, inplace=True)
            df_sl.to_csv(saveloc + method4_dir + station + '_slevels_' + str(epsilon) + '-K.csv', index=False)
            del df_sl
    except:
        print('Failed to finish ' + station)
    del df
    
    
    
    
    
    
    
    
    
    
    
    
    
    

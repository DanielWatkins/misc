"""Downloads all Arctic (Lat > 65) radiosonde data that extend at least from 1990-2020"""
import pandas as pd
import numpy as np
from siphon.simplewebservice.igra2 import IGRAUpperAir
from datetime import datetime

save_loc_soundings = '/Users/watkdani/Documents/Research/IGRA2 Paper/Data/IGRA_derived_soundings/'
save_loc_headers = '/Users/watkdani/Documents/Research/IGRA2 Paper/Data/IGRA_derived_headers/'
station_list_loc = '/Users/watkdani/Documents/Research/IGRA2 Paper/igra2-station-list.txt'



station_list = pd.read_fwf(station_list_loc,
                          header=None, names=['station_id', 'lat', 'lon', 'elevation', 'name', 'start_year', 'end_year', 'number'])
idx = (station_list.lat >= 65) & (station_list.end_year == 2020)
idx = idx & (station_list.begin_year <= 1990)
arctic_stations = station_list.loc[idx]
for station in arctic_stations.station_id:
    dates = [datetime(1990, 1, 1, 0), datetime(2020, 7, 1, 23)]
    df, header = IGRAUpperAir.request_data(dates, station, derived=True)
    df.to_csv(save_loc_soundings + station + '-igra2-derived.csv')
    header.to_csv(save_loc_headers + station + '-igra2-derived.csv')
    print(station + ' finished')
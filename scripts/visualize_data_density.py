"""Uses proplot to make a heatmap chart of available soundings as a percentage of expected soundings."""
import pandas as pd
import numpy as np
import proplot as plot
data_loc = 'Data/IGRA_derived_headers/'
image_loc = 'Images/'
station_list = pd.read_fwf('/Users/watkdani/Documents/Research/IGRA2 Paper/igra2-station-list.txt',
                          header=None, names=['station_id', 'lat', 'lon', 'elevation',
                                              'name', 'start_year', 'end_year', 'number'])
arctic_stations = station_list.loc[(station_list.lat >= 65) & (station_list.end_year == 2020)]

# Goal: heatmap plot 
full_list = []
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
for station in arctic_stations.station_id:
    df = pd.read_csv(data_loc + station + '-igra2-derived.csv')
    df['date'] = pd.to_datetime(df.date)
    df['date'] = [adjust_time(d) for d in df.date]
    df['year'] = df.date.dt.year
    df['month'] = df.date.dt.month
    df['hour'] = df.date.dt.hour
    df = df.loc[(df.hour == 0) | (df.hour == 12)].copy()
    counts = df.groupby(['year', 'month', 'hour']).count()['day'].reset_index()
    counts['station'] = station
    full_list.append(counts)
df = pd.concat(full_list)
df.rename({'inv_height': 'count'}, axis=1, inplace=True)
df.reset_index(drop=True, inplace=True)
df['date'] = [pd.to_datetime(x) for x in df['year'].map(str) + df['month'].map(lambda x: '-' + str(x).zfill(2)) + '-01']
df.rename({'day': 'count'}, axis=1, inplace=True)
df = df.loc[df.year < 2020]
df00 = df.loc[df.hour==0].pivot_table(index='station', values='count', columns='date')
df12 = df.loc[df.hour==12].pivot_table(index='station', values='count', columns='date')

df00.fillna(0, inplace=True)
df12.fillna(0, inplace=True)
df00 = df00/df00.columns.days_in_month*100
df00.astype(int)

df12 = df12/df12.columns.days_in_month*100
df12.astype(int)

df00[df00.values==0] = np.nan
df12[df12.values==0] = np.nan

fig, axs = plot.subplots(ncols=1, span=False, width=10, height=6)


data = df00

x = data.columns.values
y = data.index.values
y = np.arange(len(data.index))
pcm = axs[0].pcolor(data.values, cmap='PuBuGn', levels=10)
axs.colorbar(pcm, 'r', title='Percentage of daily soundings available')
axs.format(
    xlabel='Year', ylabel='Station code', 
    suptitle='00Z Sounding availability at Arctic stations, 1990-2019',
    ylocator=np.arange(len(data.index)),
    yticklabels=list(data.index.values),
    ytickminor=False,
    xlocator=np.arange(len(data.columns))[::12],
    xticklabels=[str(x) for x in np.arange(1990,2020)],
    xtickminor=False,
    xrotation=90)

fig.save(image_loc + 'available_soundings_00z.pdf')

fig, axs = plot.subplots(ncols=1, span=False, width=10, height=6)
data = df12
x = data.columns.values
y = data.index.values
y = np.arange(len(data.index))
pcm = axs[0].pcolor(data.values, cmap='PuBuGn', levels=10)
axs.colorbar(pcm, 'r', title='Percentage of daily soundings available')
axs.format(
    xlabel='Year', ylabel='Station code', 
    suptitle='12Z Sounding availability at Arctic stations, 1990-2019',
    ylocator=np.arange(len(data.index)),
    yticklabels=list(data.index.values),
    ytickminor=False,
    xlocator=np.arange(len(data.columns))[::12],
    xticklabels=[str(x) for x in np.arange(1990,2020)],
    xtickminor=False,
    xrotation=90)

fig.save(image_loc + 'available_soundings_12z.pdf')









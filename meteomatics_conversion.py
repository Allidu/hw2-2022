import pandas as pd
import xarray as xr
import scipy
df = pd.read_csv('meteomatics_output_example.csv')
df.set_index(["Group Start","Group End","Time Start","Time End","Temperature (deg C)", "Dew Point", "Relative Humidity", "Wind Speed (m per s)","Wind Direction","Pressure (hPa)","Cloud Cover","Total Radiation","Diffuse Radiation","Direct Radiation"])
xr = df.to_xarray()
nc=xr.to_netcdf('my_netcdf.nc')

'''
CHANGE IT SO THAT THE GROUPS OUTPUT THE LONGITUDE LATITUDE WE INPUT SO THE RESULTING NETCDF IS VALID
Command line arguments to run this script are: <# of days in the future you want prediction to go>
<title of location csv (excluding .csv)> <timezone offset> <output folder name>

Version 1
'''

import csv
import os
import shutil
import argparse
import pandas as pd
import meteomatics_conversion.api as api
from math import sqrt, atan2, degrees, log
from datetime import datetime, timedelta

def query_meteomatics(tz_offset,folder_name, num_days, location_file):
    # Credentials:
    username = 'umsolarcar_shaheen'
    password = 'OTN2nIbD0j7gs'
    output_file = 'Meteomatics_' + location_file + '_' + datetime.now().strftime("%B_%d_%y")

    # Timeseries Variables (must be in UTC):
    start_date = datetime.utcnow()
    if start_date.day < datetime.now().day:
        start_date += timedelta(days = 1)
    start_date = start_date.replace(hour=0, minute=30, second=0, microsecond=0) - timedelta(hours = tz_offset)
    end_date = start_date + timedelta(days=1 + num_days)
    interval = timedelta(hours=1)

    # Can only query 8 coordinates at a time, can query up to 10 time a day for a total of 80 points
    latlong = read_location_file(location_file)

    #for coords in range(0, len(latlong), 8):

    parameters = ['t_2m:C',
                  'dew_point_2m:C',
                  'relative_humidity_2m:p',
                  'relative_humidity_700hPa:p',
                  'relative_humidity_500hPa:p',
                  'relative_humidity_300hPa:p',
                  'wind_speed_u_10m:ms',
                  'wind_speed_v_10m:ms',
                  'sfc_pressure:hPa',
                  'air_density_2m:kgm3',
                  'low_cloud_cover:p',
                  'medium_cloud_cover:p',
                  'high_cloud_cover:p',
                  'effective_cloud_cover:p',
                  'diffuse_rad:W',
                  'direct_rad:W',
                  'global_rad:W',
                  'clear_sky_rad:W',
                  'precip_1h:mm',
                  'prob_precip_1h:p',
                  'fosberg_fire_weather_index:idx',
                  'sunrise:sql',
                  'sunset:sql',
                  'mixing_ratio_1000hPa:kgkg',
                  'mixing_ratio_700hPa:kgkg',
                  't_1000hPa:C',
                  't_700hPa:C',
                  'geopotential_height_1000hPa:m',
                  'geopotential_height_700hPa:m',
                  'cloud_base_agl:m',
                  'cape:Jkg',
                  'lifted_index:K']

    model = 'mix'
    ens_select = None
    cluster_select = None
    interp_select = 'gradient_interpolation'

    try:
        # resultDF = api.query_time_series(latlong, start_date, end_date, interval, parameters,
        #                               username, password, model, ens_select, interp_select,
        #                               cluster_select = cluster_select)

        convert_to_sim_input(latlong, output_file, folder_name, resultDF, start_date + timedelta(hours = tz_offset), end_date+ timedelta(hours = tz_offset), tz_offset)
    except Exception as e:
        print("Failed, the exception is {}".format(e))

def convert_to_sim_input(latlong, output_file, folder_name, tempDF, start_date, end_date, tz_offset):
    # Weather file data lists
    group = []
    lat = []
    long = []
    time_start = []
    time_end = []
    temp = []
    dew = []
    rel_hum = []
    rel_hum_700mb = []
    rel_hum_500mb = []
    rel_hum_300mb = []
    wind_speed = []
    wind_dir = []
    pressure = []
    density = []
    total_cloud = []
    low_cloud = []
    mid_cloud = []
    high_cloud = []
    ghi = []
    diffuse = []
    dni = []
    clear_sky = []
    wind = []
    precip = []
    prob_precip = []
    bushfire = []
    sunrise = []
    sunset = []
    elr = []
    dalr = []
    salr= []
    lcl = []
    cape = []
    lifted_index = []

    hours_passed = 0
    wx_group = 0
    for i in range(0, len(tempDF['t_2m:C'])) :
        # Restart date after each weather group is complete
        if end_date < (start_date + timedelta(hours=hours_passed)):
            hours_passed = 0
            wx_group += 1
        lat.append(latlong[wx_group][0])
        long.append(latlong[wx_group][1])

        time_start.append(start_date + timedelta(hours=hours_passed))
        hours_passed += 1
        time_end.append(start_date + timedelta(hours=hours_passed))
        temp.append(tempDF['t_2m:C'][i])
        dew.append(tempDF['dew_point_2m:C'][i])
        rel_hum.append(float(tempDF['relative_humidity_2m:p'][i]) / 100)
        rel_hum_700mb.append(float(tempDF['relative_humidity_700hPa:p'][i]) / 100)
        rel_hum_500mb.append(float(tempDF['relative_humidity_500hPa:p'][i]) / 100)
        rel_hum_300mb.append(float(tempDF['relative_humidity_300hPa:p'][i]) / 100)
        wind.append(convert_wind(float(tempDF['wind_speed_u_10m:ms'][i]), float(tempDF['wind_speed_v_10m:ms'][i])))
        wind_speed.append(wind[-1][0])
        wind_dir.append(wind[-1][1])
        pressure.append(tempDF['sfc_pressure:hPa'][i])
        density.append(tempDF['air_density_2m:kgm3'][i])
        low_cloud.append(float(tempDF['low_cloud_cover:p'][i]) / 100)
        mid_cloud.append(float(tempDF['medium_cloud_cover:p'][i]) / 100)
        high_cloud.append(float(tempDF['high_cloud_cover:p'][i]) / 100)
        total_cloud.append(float(tempDF['effective_cloud_cover:p'][i]) / 100)
        ghi.append(tempDF['global_rad:W'][i])
        diffuse.append(tempDF['diffuse_rad:W'][i])
        dni.append(tempDF['direct_rad:W'][i])
        clear_sky.append(tempDF['clear_sky_rad:W'][i])
        precip.append(tempDF['precip_1h:mm'][i])
        prob_precip.append(float(tempDF['prob_precip_1h:p'][i]) / 100)
        bushfire.append(float(tempDF['fosberg_fire_weather_index:idx'][i]) / 100)
        sunrise.append(tempDF['sunrise:sql'][i] + timedelta(hours=tz_offset))
        sunset.append(tempDF['sunset:sql'][i] + timedelta(hours=tz_offset))
        elr.append(calculate_ELR(float(tempDF['t_1000hPa:C'][i]), float(tempDF['t_700hPa:C'][i]), float(tempDF['geopotential_height_1000hPa:m'][i]), float(tempDF['geopotential_height_700hPa:m'][i])))
        dalr.append(-9.8)
        salr.append(calculate_SALR(float(tempDF['mixing_ratio_1000hPa:kgkg'][i]), float(tempDF['mixing_ratio_700hPa:kgkg'][i]), float(tempDF['t_1000hPa:C'][i]), float(tempDF['t_700hPa:C'][i])))
        lcl.append(float(tempDF['cloud_base_agl:m'][i]) / 1000)
        cape.append(tempDF['cape:Jkg'][i])
        lifted_index.append(tempDF['lifted_index:K'][i])

    # Simulator file
    concatenated_data = {'Latitude': lat, 'Latitude': long, 'Time Start': time_start,
             'Time End': time_end, 'Temperature (deg C)': temp, 'Dew Point': dew, 'Relative Humidity': rel_hum,
             'Wind Speed (m/s)': wind_speed, 'Wind Direction': wind_dir, 'Pressure (hPa)': pressure, 'Air Density': density,
             'Cloud Cover': total_cloud,'Total Radiation': ghi, 'Diffuse Radiation': diffuse, 'Direct Radiation': dni}
    weatherDF = pd.DataFrame(concatenated_data)
    weatherDF.to_csv(folder_name + '/' + output_file + '.csv', na_rep = "", index = False)

    # Literally everything else that could be useful file
    concatenated_data = {'Latitude': lat, 'Longitude': long, 'Time Start': time_start, 'Time End': time_end,
             'Projected Sunrise': sunrise, 'Projected Sunset': sunset,'Temperature (deg C)': temp,
             'Dew Point': dew,'Relative Humidity': rel_hum, 'Relative Humidity 700mb': rel_hum_700mb,
             'Relative Humidity 500mb': rel_hum_500mb, 'Relative Humidity 300mb': rel_hum_300mb,
             'Wind Speed (m/s)': wind_speed, 'Wind Direction': wind_dir, 'Pressure (hPa)': pressure,
             'Low Cloud Cover': low_cloud, 'Mid Cloud Cover': mid_cloud, 'High Cloud Cover': high_cloud, 'Effective Cloud Cover': total_cloud,
             'Total Radiation': ghi, 'Diffuse Radiation': diffuse, 'Direct Radiation': dni, 'Clear Sky Radiation': clear_sky,
             'Precipitation (mm/hr)' : precip, 'Precip Probability': prob_precip, 'Bushfire Probability': bushfire, 'CAPE (J/kg)': cape,
             'Lifted Index (K)': lifted_index, 'ELR (degC / km)': elr, 'DALR (degC / km)': dalr, 'SALR (degC / km)': salr, 'LCL (km)': lcl}

    weatherDF = pd.DataFrame(concatenated_data)
    weatherDF.to_csv(folder_name + '/' + output_file + '_everything_else.csv', na_rep = "", index = False)

# Environmental Lapse Rate
def calculate_ELR(temp_1000mb, temp_700mb, height_1000mb, height_700mb):
    return (temp_700mb - temp_1000mb) / (height_700mb / 1000 - height_1000mb / 1000)

# Saturated Adiabatic Lapse Rate
def calculate_SALR(r_1000mb, r_700mb, t_1000mb, t_700mb):
    h_vap = 2453000 #J/kg
    Cpd = 1003.5 #J/kg degC

    return(-9.8 / (1 + (h_vap / Cpd) * ((r_700mb - r_1000mb) / (t_700mb - t_1000mb))))

# Converts wind at 10m to wind speed and direction at 1m and returns it as a tuple
def convert_wind(wind_u, wind_v):
    # Find wind direction and wind speed from u and v vector components
    theta = degrees(atan2(wind_v, wind_u))
    phi = 90 - theta
    wind_direction = str((180 + phi) % 360)
    wind_speed = sqrt((wind_v)**2 + (wind_u)**2)

    # Meters above the ground at which zero wind speed is achieved as a
	# result of flow obstacles such as trees or building. It can be
	# approximated as 2/3 to 3/4 of the average height of obstacles
	# Note: The average height in the outback for plants with a slight
	# addition for the height of walls is around 0.6-0.8 meters
    d = 0.5
	# surface roughness in meters
	# 0.03: Open flat terrain; grass, few isolated obstacles
	# 0.10: low crops; occational large obstacles
	# 0.25: high crops; scattered obstacles
    z0 = 0.06
	# original height of the wind in meters
    z1 = 10
	# converted height of the wind in meters
    z2 = 1
    wind_speed = wind_speed * (log((z2 - d) / z0) / log((z1 - d) / z0))
    return (wind_speed, wind_direction)

# Takes in a index/lat/long csv and outputs a list of lat/longs
def read_location_file(location_file):
    location_file += '.csv'
    location_list = []

    with open(location_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            else:
                location_list.append([float(row[1]), float(row[2])])

    return location_list

if __name__ == "__main__":

	#Argument Parser
	parser = argparse.ArgumentParser()
	# Adds arguments to the parser
	parser.add_argument("-n", "--num_days", type=float, default=None,
		required=True, help="Number of days in the future you want predictions to go")
	parser.add_argument("-f", "--loc_file", type=str, default=None,
		required=True, help="Title of location csv (excluding .csv)")
	parser.add_argument("-t", "--tz_offset", type=float, default=None,
		required=True, help="Timezone offset. Integers and floats allowed.")
	parser.add_argument("-o", "--output_folder", type=str, default=None,
		required=False, help="Name of the folder to output csv to")

        # Parses command line arguments
	args = parser.parse_args()
	output_dir = ""

	if args.output_folder:
        # Make a folder to export data if it doesn't already exist
		output_dir = args.output_folder
		if not os.path.exists(output_dir):
			os.mkdir(output_dir)

	else:
		output_dir = os.path.dirname(os.path.realpath(__file__)) + os.sep + "AUS_WEATHER_DATA" + os.sep
		if (not os.path.exists(output_dir)):
			os.makedirs(output_dir)

		output_dir += "BOM_DOWNLOADS" + os.sep
		if (not os.path.exists(output_dir)):
			os.makedirs(output_dir)

		output_dir += "BOM_{}".format((datetime.utcnow() + timedelta(hours = 9.5)).strftime("%Y%m%d_%H%M")) + os.sep
		if os.path.exists(output_dir):
            # Removes the directory and all folders/files within it
            # MAKE SURE YOU ARE INPUTTING THE RIGHT FOLDER
			shutil.rmtree(output_dir, ignore_errors=True)
		os.mkdir(output_dir)

    # FUNCTION CALL
	query_meteomatics(args.tz_offset, output_dir, args.num_days, args.loc_file)

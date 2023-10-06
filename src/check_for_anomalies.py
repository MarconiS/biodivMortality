import ee
import pandas as pd
import matplotlib.pyplot as plt

ee.Initialize()

# Function to extract values
def extract_values(image_collection, geom):
    values = image_collection.getRegion(geom, 500).getInfo()
    return values


def get_anomalies(aoi, start_date, end_date):
    # Load temperature data
    temp_data = ee.ImageCollection('MODIS/006/MOD11A1') \
        .filterBounds(aoi) \
        .filterDate(ee.Date(start_date), ee.Date(end_date)) \
        .select('LST_Day_1km')

    # Load precipitation data
    precip_data = ee.ImageCollection('UCSB-CHG/CHIRPS/PENTAD') \
        .filterBounds(aoi) \
        .filterDate(ee.Date(start_date), ee.Date(end_date)) \
        .select('precipitation')
    
    # Extract temperature and precipitation values
    temp_values = extract_values(temp_data, aoi)
    precip_values = extract_values(precip_data, aoi)

    # Convert to DataFrame
    temp_df = pd.DataFrame(temp_values[1:], columns=temp_values[0])
    precip_df = pd.DataFrame(precip_values[1:], columns=precip_values[0])

    # Convert time to datetime format and set as index
    temp_df['time'] = pd.to_datetime(temp_df['time'], unit='ms')
    temp_df.set_index('time', inplace=True)

    precip_df['time'] = pd.to_datetime(precip_df['time'], unit='ms')
    precip_df.set_index('time', inplace=True)

    '''
    # Plotting the data
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.plot(temp_df.index, temp_df['LST_Day_1km'])
    plt.title('Temperature Over Time')
    plt.xlabel('Time')
    plt.ylabel('Temperature (K)')

    plt.subplot(1, 2, 2)                
    plt.plot(precip_df.index, precip_df['precipitation'])
    plt.title('Precipitation Over Time')
    plt.xlabel('Time')
    plt.ylabel('Precipitation (mm)')

    plt.tight_layout()
    plt.show()
    '''



    mean_temp = temp_df['LST_Day_1km'].mean()
    std_temp = temp_df['LST_Day_1km'].std()

    threshold_high = mean_temp + 3 * std_temp
    threshold_low = mean_temp - 3 * std_temp

    extreme_high_temp_events = temp_df[temp_df['LST_Day_1km'] > threshold_high]
    extreme_low_temp_events = temp_df[temp_df['LST_Day_1km'] < threshold_low]

    if extreme_high_temp_events.empty:
        print("No extreme high temperature events found.")
    else:
        print("Extreme high temperature events found:")
        print(extreme_high_temp_events)
    if extreme_low_temp_events.empty:
        print("No extreme low temperature events found.")
    else:
        print("Extreme low temperature events found:")
        print(extreme_low_temp_events)

    # same with precipitation
    mean_prec = precip_df['precipitation'].mean()
    std_prec = precip_df['precipitation'].std()
    #recalculate using only dates with precipitation
    precip_df2 = precip_df[precip_df['precipitation'] > 0]
    mean_prec = precip_df2['precipitation'].mean()
    std_prec = precip_df2['precipitation'].std()
    #recalculate using only dates between april and september
    precip_df3 = precip_df2[precip_df2.index.month > 3]
    precip_df3 = precip_df3[precip_df3.index.month < 10]
    mean_prec = precip_df3['precipitation'].mean()
    std_prec = precip_df3['precipitation'].std()

    threshold_high = mean_prec + 3 * std_prec
    threshold_low = mean_prec - 3 * std_prec

    extreme_low_prec_events = precip_df[precip_df['precipitation'] < threshold_low]

    if extreme_low_prec_events.empty:
        print("No extreme low prec events found.")
    else:
        print("Extreme low prec events found:")
        print(extreme_low_prec_events)

# Define the area of interest (AOI) using a point (longitude, latitude)
aoi = ee.Geometry.Point([-72.17265, 42.53691])  # Example coordinates
# Define the time range
start_date = '1990-01-01'
end_date = '2020-12-31'
get_anomalies(aoi, start_date, end_date)
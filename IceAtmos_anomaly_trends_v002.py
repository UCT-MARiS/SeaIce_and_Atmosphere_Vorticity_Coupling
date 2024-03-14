#%%
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import matplotlib.dates as mdates
import datetime as datetime
import numpy as np
from matplotlib.lines import Line2D
 
if __name__ == '__main__':
    
    sea             = 'weddell'
    version         = 'v001'

    #column_name1    = 'ice_pos_count'
    #column_name2    = 'ice_neg_count'
    
    years           = mdates.YearLocator()   # every year
    months          = mdates.MonthLocator()  # every month
    years_fmt       = mdates.DateFormatter('%Y')
    months_fmt      = mdates.DateFormatter('%M')

    #iterate through corresponding csv files, concat into single dataframe
    df_main       = pd.DataFrame()
    for year in range(1991,2021):
        
        fdir       = f'/home/waynedj/Projects/seaicedrift_correlation/intermediate_data/24hr/vorticity/iceatmos/{version}/{sea}/'
        fname      = f'VortFieldOverlap_{sea}_{year}_{version}.csv'
        df         = pd.read_csv(fdir+fname)
        df_main    = pd.concat([df_main, df])

    # Convert date column to datetime objects
    df_main['date'] = pd.to_datetime(df_main['date'])

    #df_main = df_main[df_main['date'].dt.month != 4]


    # Create year and month column
    df_main['year']  = df_main['date'].dt.year
    df_main['month'] = df_main['date'].dt.month


    df_main['atmos_zeros'] = df_main['atmos_total_count'] - (df_main['atmos_pos_count']+df_main['atmos_neg_count'])
    df_main['ice_zeros']   = df_main['ice_total_count'] - (df_main['ice_pos_count']+df_main['ice_neg_count'])

    column_name1    = 'atmos_zeros'
    column_name2    = 'ice_zeros'

    # Compute monthly average
    df_main['monthlymean1'] = df_main.groupby('month')[column_name1].transform('mean')
    df_main['monthlymean2'] = df_main.groupby('month')[column_name2].transform('mean')

    # Compute daily anomaly
    df_main['anomaly1'] = df_main[column_name1] - df_main['monthlymean1']
    df_main['anomaly2'] = df_main[column_name2] - df_main['monthlymean2']

    
    # Extract year and month
    
    #df_main['month'] = df_main['date'].dt.month

    # Calculate average measurement per month per year
    monthly_anom1 = df_main.groupby(['year', 'month'])['anomaly1'].mean()
    monthly_anom2 = df_main.groupby(['year', 'month'])['anomaly2'].mean()

    # generate x-axis values (1st is arbitary choice)
    date_list = []
    for year in range(1991,2021,1):
        for month in range(4,12,1):
            date_list.append(datetime.datetime(year,month,1,0,0,0))

    # Set the colors based on positive or negative values
    colors1 = ['red' if val >= 0 else 'blue' for val in monthly_anom1]
    colors2 = ['red' if val >= 0 else 'blue' for val in monthly_anom2]

    #Figure modification
    fig = plt.figure(figsize=(20,10), facecolor='w')
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    # Plot the bar graph
    ax1.bar(date_list, monthly_anom1, color=colors1, width=10)
    ax1.set_title(f'{sea}', fontsize=20)
    ax1.set_ylabel(column_name1, fontsize = 20)
    ax1.xaxis.set_major_locator(years)
    ax1.xaxis.set_major_formatter(years_fmt)
    ax1.xaxis.set_minor_locator(months)
    ax1.format_xdata = mdates.DateFormatter('%Y-%m-%d')
    #ax1.set_ylim([-0.0000005,0.0000005])
    ax1.set_xlim([datetime.datetime(1991,1,1), datetime.datetime(2021,1,1)])
    ax1.tick_params(axis='both', which='major')
    ax1.tick_params(axis='x', which='major', labelsize=16, labelrotation = 45)
    ax1.grid(True)
    #plt.xticks(rotation='vertical')


    # Plot the bar graph
    ax2.bar(date_list, monthly_anom2, color=colors2, width=10)
    #ax2.set_xlabel("Date", fontsize = 20)
    ax2.set_ylabel(column_name2, fontsize = 20)
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_major_formatter(years_fmt)
    ax2.xaxis.set_minor_locator(months)
    ax2.format_xdata = mdates.DateFormatter('%Y-%m-%d')
    #ax1.set_ylim([-0.00000015,0.00000015])
    ax2.set_xlim([datetime.datetime(1991,1,1), datetime.datetime(2021,1,1)])
    ax2.tick_params(axis='both', which='major')
    ax2.tick_params(axis='x', which='major', labelsize=16, labelrotation = 45)
    ax2.grid(True)
    #plt.xticks(rotation='vertical')
    plt.subplots_adjust(hspace=0.4)

    #df_main.to_csv(f'/home/waynedj/Projects/seaicedrift_correlation/images/iceatmos_monthlyanomaly/data_iceatmos_{sea}_{version}.csv', index=False)
    #plt.savefig(f'/home/waynedj/Projects/seaicedrift_correlation/images/iceatmos_monthlyanomaly/iceatmos_{sea}_col1-{column_name1}_col2-{column_name2}_{version}.png', bbox_inches = 'tight',dpi=400)
    

    #pearson          = np.corrcoef(monthly_anom1, monthly_anom2, rowvar=True)#[0][1]
    #print(pearson)

# %%

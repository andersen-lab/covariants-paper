import pandas as pd 
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns 
import matplotlib.dates as mdates
import matplotlib.pyplot as plt 
import datetime as dt
import numpy as np
import matplotlib.dates as mdates
import numpy as np

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


# import dash
import pandas as pd
import pickle
import json
import requests
from scipy import signal
from datetime import date,timedelta
import yaml
import copy



def non_uniform_savgol(x, y, window, polynom):
    """
    Applies a Savitzky-Golay filter to y with non-uniform spacing
    as defined in x

    ----------
    x : array_like
        List of floats representing the x values of the data (MUST BE ORDERED)
    y : array_like
        List of floats representing the y values. Must have same length
        as x
    window : int (odd)
        Window length of datapoints. Must be odd and smaller than x
    polynom : int
        The order of polynom used. Must be smaller than the window size

    Returns
    -------
    np.array of float
        The smoothed y values
    """
    if len(x) != len(y):
        raise ValueError('"x" and "y" must be of the same size')

    if len(x) < window:
        raise ValueError('The data size must be larger than the window size')

    if type(window) is not int:
        raise TypeError('"window" must be an integer')

    if window % 2 == 0:
        raise ValueError('The "window" must be an odd integer')

    if type(polynom) is not int:
        raise TypeError('"polynom" must be an integer')

    if polynom >= window:
        raise ValueError('"polynom" must be less than "window"')

    half_window = window // 2
    polynom += 1

    # Initialize variables
    A = np.empty((window, polynom))     # Matrix
    tA = np.empty((polynom, window))    # Transposed matrix
    t = np.empty(window)                # Local x variables
    y_smoothed = np.full(len(y), np.nan)

    # Start smoothing
    for i in range(half_window, len(x) - half_window, 1):
        # Center a window of x values on x[i]
        for j in range(0, window, 1):
            t[j] = x[i + j - half_window] - x[i]

        # Create the initial matrix A and its transposed form tA
        for j in range(0, window, 1):
            r = 1.0
            for k in range(0, polynom, 1):
                A[j, k] = r
                tA[k, j] = r
                r *= t[j]

        # Multiply the two matrices
        tAA = np.matmul(tA, A)

        # Invert the product of the matrices
        tAA = np.linalg.inv(tAA)

        # Calculate the pseudoinverse of the design matrix
        coeffs = np.matmul(tAA, tA)

        # Calculate c0 which is also the y value for y[i]
        y_smoothed[i] = 0
        for j in range(0, window, 1):
            y_smoothed[i] += coeffs[0, j] * y[i + j - half_window]

        # If at the end or beginning, store all coefficients for the polynom
        if i == half_window:
            first_coeffs = np.zeros(polynom)
            for j in range(0, window, 1):
                for k in range(polynom):
                    first_coeffs[k] += coeffs[k, j] * y[j]
        elif i == len(x) - half_window - 1:
            last_coeffs = np.zeros(polynom)
            for j in range(0, window, 1):
                for k in range(polynom):
                    last_coeffs[k] += coeffs[k, j] * y[len(y) - window + j]

    # Interpolate the result at the left border
    for i in range(0, half_window, 1):
        y_smoothed[i] = 0
        x_i = 1
        for j in range(0, polynom, 1):
            y_smoothed[i] += first_coeffs[j] * x_i
            x_i *= x[i] - x[half_window]

    # Interpolate the result at the right border
    for i in range(len(x) - half_window, len(x), 1):
        y_smoothed[i] = 0
        x_i = 1
        for j in range(0, polynom, 1):
            y_smoothed[i] += last_coeffs[j] * x_i
            x_i *= x[i] - x[-half_window - 1]

    return y_smoothed

ww_dict = {'Point Loma':['../data/lineage-prevalence/PointLoma_sewage_seqs.csv','../data/qPCR/PointLoma_sewage_qPCR.csv'],
           'Encina':['../data/lineage-prevalence/Encina_sewage_seqs.csv','../data/qPCR/Encina_sewage_qPCR.csv'],
           'South Bay':['../data/lineage-prevalence/SouthBay_sewage_seqs.csv','../data/qPCR/SouthBay_sewage_qPCR.csv']}
for site, files in zip(ww_dict.keys(),ww_dict.values()):
    df = pd.read_csv(f'{files[0]}')

    df['Date'] = pd.to_datetime(df['Date'])
    df.columns = [dfc.split(' (')[0] for dfc in df.columns]
    df =df.set_index('Date')
    df = df[df.index>='2022-12-31']
    df = df[df.index<='2025-01-01']
    df = df.dropna(axis = 0, how = 'all')
    df = df.fillna(0)
    df = df/100.
    # df['Other'] = 1. - df['Delta']- df['Omicron']- df['Alpha']
    # df = df.drop(columns=['Alpha','Delta','Omicron','Other'])
    df = df.drop(columns=['Other'])
    df = df[df.columns[df.sum(axis=0) > 0.01]]

    cdf = pd.read_csv(f'{files[1]}')
    cdf['Sample_Date'] = pd.to_datetime(cdf['Sample_Date'])
    cdf =cdf.set_index('Sample_Date')
    sharedInds = np.sort(list(set(cdf.index) & set(df.index)))
    cdf = cdf.loc[sharedInds]
    df = df.loc[sharedInds]
    df = df[~df.index.duplicated(keep='last')]
    scaleddf = df.mul(cdf['Mean viral gene copies/L'],axis=0)



    xlim = pd.date_range(df.index.min(), df.index.max())
    numberDates = [dvi.value/10**11 for dvi in df.index]
    for c in scaleddf.columns:
        scaleddf[c] = non_uniform_savgol(numberDates,scaleddf[c],7,1)
    scaleddf[scaleddf<0]=0.
    df = scaleddf.divide(np.sum(scaleddf,axis=1),axis=0)

    #------ build colormaps

    with open( "../data/plot_config.yml", "r" ) as f :
        plot_config = yaml.load( f, Loader=yaml.FullLoader )

    #--- borrowed from SEARCH wastewater surveillance dashboard, coordinated by Scripps Research.--#
    def convert_rbg_to_tuple( rgb ):
        rgb = rgb.lstrip( "#" )
        return tuple( int( rgb[i :i + 2], 16 ) for i in (0, 2, 4) )
    def convert_tuple_to_rgb( r, g, b ):
        return '#%02x%02x%02x' % (int(r), int(g), int(b))
    def lighten_field( value, alpha, gamma=2.2 ):
        return pow( pow(255, gamma) * (1 - alpha) + pow( value, gamma ) * alpha, 1 / gamma)
    def lighten_color( r, g, b, alpha, gamma=2.2 ):
        return lighten_field(r, alpha, gamma ), lighten_field( g, alpha, gamma ), lighten_field( b, alpha, gamma )

    children_dict = dict()
    delta = 0.15
    for key in reversed( list( plot_config.keys() ) ):
        for value in ["name", "members"]:
            assert value in plot_config[key], f"YAML entry {key} is not complete. Does not contain '{value}' entry."
        if "color" not in plot_config[key]:
            assert "parent" in plot_config[key], f"YAML entry {key} is incomplete. Must specify either a 'color' or 'parent' entry."
            if plot_config[key]["parent"] in children_dict:
                children_dict[plot_config[key]["parent"]] += 1
            else:
                children_dict[plot_config[key]["parent"]] = 1
            child_idx = children_dict[plot_config[key]["parent"]]
            parent_color = plot_config[plot_config[key]["parent"]]["color"]
            parent_color = convert_rbg_to_tuple( parent_color )
            plot_config[key]["color"] = convert_tuple_to_rgb( *lighten_color( *parent_color, alpha=1.0-(delta*child_idx) ) )

    #----#

    def inGroup(plot_config0,linName):
        for key in plot_config0.keys():
            if linName in plot_config0[key]['members']:
                return key
            else:
                for member in plot_config[key]['members']:
                    if col.startswith(member.split('.X')[0]):
                        return key
    alreadySeen =[]
    for col in df.columns:
        grp0 = inGroup(plot_config,col)
        if col!= grp0:
            try:
                
                df[grp0] = df[[grp0,col]].sum(axis=1)

            except:
                continue
            df = df.drop(columns=col)
        print(col,inGroup(plot_config,col))

    df = df[[key for key in plot_config.keys() if key in df.columns]]

    scaleddf = scaleddf[df.columns[::-1]]

    # newcolOrder = list(scaleddf.columns[0:scaleddf.shape[1]-2]) + [scaleddf.columns[-1]]+ [scaleddf.columns[-2]]
    # scaleddf = scaleddf[newcolOrder]

    colors = [plot_config[linName]['color'] for linName in scaleddf.columns]
    # quickPalette = quickPalette[::-1]
    # quickPalette = quickPalette[0:len(quickPalette)-2] + [quickPalette[-1]]+ [quickPalette[-2]]
    scaleddf = scaleddf[scaleddf.columns]
    fig, ax = plt.subplots(figsize=(13,7))
    ax.stackplot(scaleddf.index, scaleddf.to_numpy().T,
                labels=scaleddf.columns, colors=colors,edgecolor='k',linewidth=0.1)
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 4})
    ax.set_ylabel('Lineage copies/L')
    # ax.set_ylim([0, 100])

    # locator =mdates.WeekdayLocator()#bymonthday=1)
    locator = mdates.MonthLocator(bymonthday=1)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    ax.set_xlim([df.index.min(),df.index.max()])
    # plt.setp(ax.get_xticklabels(), rotation=90)
    ax.grid(axis='x',alpha=0.25)
    # ax.set_aspect(1e-5)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(0.3, -0.2),ncol=6)
    fig.tight_layout()
    plt.savefig(f'plots/scaled_lineage_deconv{site}.pdf')
    plt.close()

    df = df[scaleddf.columns]
    fig, ax = plt.subplots(figsize=(18,7))
    ax.stackplot(df.index,df.to_numpy().T,
                labels=df.columns, colors=colors,edgecolor='k',linewidth=0.1)
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 4})
    ax.set_ylabel('Variant Prevalence')
    ax.set_ylim([0, 1.])
    # plt.axvline(dt.datetime(2021, 12, 8),linestyle='--',color='black')
    # plt.axvline(dt.datetime(2021, 11, 27),linestyle='--',color='grey')

    # locator =mdates.WeekdayLocator()#bymonthday=1)
    locator = mdates.MonthLocator(bymonthday=1)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(locator))
    ax.set_xlim([df.index.min(),df.index.max()])
    ax.tick_params(axis='both', which='major', labelsize=16)

    # plt.setp(ax.get_xticklabels(), rotation=90)
    # ax.set_aspect(200)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1],loc='center left', bbox_to_anchor=(1.02, 0.5),ncol=1,fontsize=16)
    fig.tight_layout()
    plt.savefig(f'plots/lineage_deconv{site}.pdf',bbox_inches='tight',)
    plt.close()

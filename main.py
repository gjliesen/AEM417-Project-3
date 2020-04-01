# Gregory Liesen
# AEM 417 Project 3

# Imports
import math as m
import datetime
import time
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import signal
from scipy import linalg as la
from scipy.interpolate import interp1d
import georinex as gr
from astropy.time import Time

# CONSTANTS
U = 3.986005e14 # m^3/s^2 WGS84 value of Earth's Gravitational Constant
C = 2.99792458e8 # m/s GPS Value for speed of light
O_e = 7.2921151467e-5 # rad/s WGS84 Value of Earth's rotation rate

def read_rinex(file):
    brdc = gr.load(file)
    return brdc

def sv_position(br_df):
    global U
    br_df['A'] = br_df['sqrtA']**2
    br_df['N'] = np.sqrt(U * br_df['A']**-1) + br_df['DeltaN']
    return br_df

def get_sat_data(rx, sv):
    df = rx.sel(sv=sv).to_dataframe()
    df = df.fillna(method='ffill')
    return df


def read_icp(file):
    df = pd.read_csv(file , delimiter='\t', names = ['Time', '2', '3', '4', '5', '6', '7', 'Pseudorange'])
    df = df.drop(['2', '3', '4', '5', '6', '7'], axis=1)
    return df


def get_base_position(file):
    df = pd.read_csv(file, delimiter='\t', names = ['Time', 'GPS_Week', 'X_ECEF', 'Y_ECEF', 'Z_ECEF',
                                                    'VX', 'VY', 'VZ', '1' , '2', '3', '4', '5', '6',
                                                    '7', '8', '9', '10'])
    df = df.drop(['GPS_Week', 'VX', 'VY', 'VZ', '1' , '2', '3', '4', '5', '6', '7', '8', '9', '10'], axis=1)
    x_base = df['X_ECEF'].mean()
    y_base = df['Y_ECEF'].mean()
    z_base = df['Z_ECEF'].mean()
    return [x_base, y_base, z_base]


def gen_least_sqr():
    print('Generalized Least Squares Position')


def vertical_dop():
    print('VDOP')


def horizontal_dop():
    print('HDOP')


def position_dop():
    print('PDOP')


def mapping(longitude_min, longitude_max, latitude_min, latitude_max, longitude, latitude):
    image_box = ((longitude_min, longitude_max, latitude_min, latitude_max))
    fig, ax = plt.subplots(figsize = (8,7))
    ax.scatter(longitude, latitude, zorder=1, alpha=0.2, c ='b', s=10)
    ax.set_title('Plotting Spatial data on riyadh Map')
    ax.set_xlim(longitude_min, longitude_max)
    ax.set_ylim(latitude_min, latitude_max)


def main():
    brdc = read_rinex('brdc2930.11n')
    brG02_df= get_sat_data(brdc, 'G02')
    brG04_df= get_sat_data(brdc, 'G04')
    brG05_df= get_sat_data(brdc, 'G05')
    brG09_df= get_sat_data(brdc, 'G09')
    brG10_df= get_sat_data(brdc, 'G10')
    brG12_df= get_sat_data(brdc, 'G12')
    brG17_df= get_sat_data(brdc, 'G17')
    brG23_df= get_sat_data(brdc, 'G23')
    brG25_df= get_sat_data(brdc, 'G25')

    #Data Base
    b_sat2_df = read_icp('data_base/icp_sat2.txt')
    b_sat4_df = read_icp('data_base/icp_sat4.txt')
    b_sat5_df = read_icp('data_base/icp_sat5.txt')
    b_sat9_df = read_icp('data_base/icp_sat9.txt')
    b_sat10_df = read_icp('data_base/icp_sat10.txt')
    b_sat12_df = read_icp('data_base/icp_sat12.txt')
    b_sat17_df = read_icp('data_base/icp_sat17.txt')
    b_sat23_df = read_icp('data_base/icp_sat23.txt')
    b_sat25_df = read_icp('data_base/icp_sat25.txt')

    #Data Rover
    r_sat2_df = read_icp('data_rover/icp_sat2.txt')
    r_sat4_df = read_icp('data_rover/icp_sat4.txt')
    r_sat5_df = read_icp('data_rover/icp_sat5.txt')
    r_sat9_df = read_icp('data_rover/icp_sat9.txt')
    r_sat10_df = read_icp('data_rover/icp_sat10.txt')
    r_sat12_df = read_icp('data_rover/icp_sat12.txt')
    r_sat17_df = read_icp('data_rover/icp_sat17.txt')
    r_sat23_df = read_icp('data_rover/icp_sat23.txt')
    r_sat25_df = read_icp('data_rover/icp_sat25.txt')

    #Base Position
    Base_Vector = get_base_position('data_base/ecef_rx0.txt')
    print(Base_Vector)

    print('Done')
main()

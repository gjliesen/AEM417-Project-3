# Gregory Liesen
# AEM 417 Project 3

# Imports
import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import signal
from scipy import linalg as la
from scipy.interpolate import interp1d
import georinex as gr

# CONSTANTS
U = 3.986005e14 # m^3/s^2 WGS84 value of Earth's Gravitational Constant
C = 2.99792458e8 # m/s GPS Value for speed of light
O_e = 7.2921151467e-5 # rad/s WGS84 Value of Earth's rotation rate

def read_rinex(file):
    df = gr.load(file).to_dataframe()
    return df

def sv_position(br_df):
    global U
    br_df['A'] = br_df['sqrtA']**2
    br_df['corr_mean_motion'] = np.sqrt(U * br_df['A']**-1) + br_df['DeltaN']


def gen_least_sqr():
    print('Generalized Least Squares Position')


def vdop():
    print('VDOP')


def hdop():
    print('HDOP')


def pdop():
    print('PDOP')


def mapping(longitude_min, longitude_max, latitude_min, latitude_max, longitude, latitude):
    image_box = ((longitude_min, longitude_max, latitude_min, latitude_max))
    fig, ax = plt.subplots(figsize = (8,7))
    ax.scatter(longitude, latitude, zorder=1, alpha=0.2, c ='b', s=10)
    ax.set_title('Plotting Spatial data on riyadh Map')
    ax.set_xlim(longitude_min, longitude_max)
    ax.set_ylim(latitude_min, latitude_max)


def main():
    br_df = read_rinex('brdc2930.11n')
    igs_df = read_rinex('igs16584.sp3')
    sv_position(br_df)
    print('Done')
main()

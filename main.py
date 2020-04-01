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
OmegaDote = 7.2921151467e-5 # rad/s WGS84 Value of Earth's rotation rate
To = 417600.000
f = 1/298.257223563
e = np.sqrt(f*(2-f))
Ro = 6378137
Rp = Ro * (1-f)


def read_rinex(file):
    brdc = gr.load(file)
    times = gr.gettime(file)
    return brdc


def get_sat_data(rx, sv):
    df = rx.sel(sv=sv).to_dataframe()
    df = df.fillna(method='ffill')
    return df


def read_icp(file, name):
    df = pd.read_csv(file , delimiter='\t', names = ['Time', '2', '3', '4', '5', '6', '7', name])
    df = df.drop(['2', '3', '4', '5', '6', '7'], axis=1)
    df = df.set_index('Time')
    df = df.loc[417136:417984]
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


def sv_position(rover_icp_data, br_data):
    global U, To, OmegaDote
    temp = pd.DataFrame(index=rover_icp_data.index)
    temp['A'] = br_data.get('sqrtA')[0] ** 2
    temp['U/A^3'] = ((temp['A']**3) ** -1) * U
    temp['sqrt(U/A^3)'] = np.sqrt(temp['U/A^3'])
    temp['N'] = temp['sqrt(U/A^3)'] + br_data.get('DeltaN')[0]
    temp['te'] = temp.index
    temp['Tk'] = temp['te'] - To
    temp['Mk'] = br_data.get('M0')[0] + temp['N'] * temp['Tk']
    temp['E'] = temp['Mk']
    temp['Ratio'] = 1
    temp['RatioAbs'] = temp['Ratio']
    tolerance = 10e-8
    e = br_data.get('Eccentricity')[0]
    while temp['RatioAbs'].max() > tolerance:
        temp['CosE'] = np.cos(temp['E'])
        temp['SinE'] = np.sin(temp['E'])
        temp['E_error'] = temp['E'] - temp['SinE'].multiply(e) - temp['Mk']
        temp['E_derivative'] = temp['CosE'] * e * -1 + 1
        temp['Ratio'] = temp['E_error'].div(temp['E_derivative'])
        temp['E'] = temp['E'] - temp['Ratio']
        temp['RatioAbs'] = temp['Ratio'].abs()
    temp['E'] = temp['E'] + temp['Ratio']
    temp['SinE'] = np.sin(temp['E'])
    temp['CosE'] = np.cos(temp['E'])
    temp['nums'] = np.sqrt(1 - e**2) * temp['SinE']
    temp['dens'] = temp['CosE'] * e * -1 + 1
    temp['Sinv'] = temp['nums'] / temp['dens']
    temp['numc'] = temp['CosE'] - e
    temp['denc'] = temp['CosE'] * e * -1 + 1
    temp['Cosv'] = temp['numc'] / temp['denc']
    temp['vk'] = np.arctan2(temp['Sinv'],temp['Cosv'])
    temp['Phi_k'] = temp['vk'] + br_data.get('omega')[0]
    temp['Uk'] = temp['Phi_k'] + br_data.get('Cus')[0] * np.sin(temp['Phi_k'] * 2) + br_data.get('Cuc')[0] * \
                 np.cos(temp['Phi_k'] * 2)
    temp['rk'] = temp['A'] * (temp['CosE'] * e * -1 + 1) + br_data.get('Crs')[0]
    temp['ik'] = br_data.get('Io')[0] + br_data.get('IDOT')[0] * temp['Tk'] + br_data.get('Cis')[0] * \
                 np.sin(temp['Phi_k'] * 2) \
                 + br_data.get('Cic')[0] * np.cos(temp['Phi_k'] * 2)
    temp['Omega_k'] = temp['Tk'] * (br_data.get('OmegaDot')[0] - OmegaDote) + br_data.get('Omega0')[0] - OmegaDote * To
    temp['xk_prime'] = temp['rk'] * np.cos(temp['Uk'])
    temp['yk_prime'] = temp['rk'] * np.sin(temp['Uk'])

    #Sattelite Position calculations
    sat_pos = pd.DataFrame(index=rover_icp_data.index)
    sat_pos['x'] = temp['xk_prime'] * np.cos(temp['Omega_k']) - temp['yk_prime'] * np.cos(temp['ik']) * \
                   np.sin(temp['Omega_k'])
    sat_pos['y'] = temp['xk_prime'] * np.sin(temp['Omega_k']) + temp['yk_prime'] * np.cos(temp['ik']) * \
                   np.cos(temp['Omega_k'])
    sat_pos['z'] = temp['yk_prime'] * np.sin(temp['ik'])
    return sat_pos


def wgs_to_lla(Base_Vector):
    global e, f, Ro, Rp
    [x_base, y_base, z_base] = Base_Vector
    longitude = np.arctan2(y_base,x_base)
    p= np.sqrt(x_base**2 + y_base**2)
    E = np.sqrt(Ro**2 - Rp**2)
    f = 54 * (Rp * z_base)**2
    G = p**2 + (1 - e**2) * z_base**2 - (e * E)**2
    c_num = e**4 * f * p**2
    c_den = G**3
    c = c_num/c_den
    s = (1 + c + np.sqrt(c**2 + 2 * c)*(1/3))
    P_num = f * (3 * G**2)**-1
    P_den = (s + s**-1 + 1)**2
    P = P_num / P_den
    Q = np.sqrt(1 + 2 * e**4 * P)
    k1 = (-1 * P * e ** 2 * p) / (1 + Q)
    k2 = 0.5 * Ro ** 2 * (1 + 1 / Q)
    k3 = -1 * P * (1 - e ** 2) * (z_base ** 2 / (Q * (1 + Q)))
    k4 = -1 * 0.5 * P * p ** 2
    k5 = p - e ** 2 * (k1 + np.sqrt(k2 + k3 + k4))
    u_lla = np.sqrt(k5 ** 2 + z_base ** 2)
    V_lla = np.sqrt(k5 ** 2 + (1 - e ** 2) * z_base ** 2)
    h = u_lla * (1 - (Rp ** 2 / (Ro * V_lla)))
    zo = (Rp ** 2 * z_base) / (Ro * V_lla)
    ep = (Ro * e) / Rp
    phi = np.arctan((z_base + zo * ep ** 2) / p)
    return [phi, longitude]


def ecef_to_ned(x,y,z,R):
    print('placeholer')


def line_of_sight(icp_data, base_vector):
    print('Line of Sight')


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
    sat = ['2', '4', '5', '9', '10', '12', '17', '23', '25']
    sv = ['G02', 'G04', 'G05', 'G09', 'G10', 'G12', 'G17', 'G23', 'G25']

    # Broadcast Rinex File
    brdc = read_rinex('brdc2930.11n')
    brdc_data = []
    for i in range(9):
        brdc_data.append(get_sat_data(brdc, sv[i]).reset_index())

    # Data Base and Rover
    bdf_list = []
    for i in range(9):
        file = 'data_base/icp_sat' + sat[i] + '.txt'
        name = 'Base Pseudorange ' + sat[i]
        temp = read_icp(file,name)
        if i == 0:
            base_icp_data = temp
        else:
            base_icp_data = base_icp_data.join(temp, sort=True)

    # Data Rover
    rdf_list = []
    for i in range(9):
        file = 'data_rover/icp_sat' + sat[i] + '.txt'
        name = 'Rover Pseudorange ' + sat[i]
        temp = read_icp(file, name)
        if i == 0:
            rover_icp_data = temp
        else:
            rover_icp_data = rover_icp_data.join(temp, sort=True)


    index = [33,33,33,33,33,33,33,33,33]
    for i in range(9):
        name = 'Rover Pseudorange ' + sat[i]
        sat_pos = sv_position(rover_icp_data[name], brdc_data[i].iloc[[index[i]]].to_dict('list'))


    # Base Position
    Base_Vector = get_base_position('data_base/ecef_rx0.txt')
    [lat, long] = wgs_to_lla(Base_Vector)

    R = np.array([[-np.sin(lat) * np.cos(long), -np.sin(long), -np.cos(lat) * np.cos(long)],
                   [-np.cos(lat) * np.sin(long), np.cos(long), -np.cos(lat) * np.sin(long)],
                  [np.cos(lat), 0, -np.sin(lat)]])


    print('Done')


main()

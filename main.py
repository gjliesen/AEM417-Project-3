# Gregory Liesen
# AEM 417 Project 3

# Imports
import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import linalg as la
import georinex as gr
import navpy

# CONSTANTS
C = 2.99792458e8 # m/s GPS Value for speed of light


def read_rinex(file):
    brdc = gr.load(file)
    times = gr.gettime(file)
    return brdc


def get_sat_data(rx, sv):
    df = rx.sel(sv=sv).to_dataframe()
    df = df.fillna(method='ffill')
    return df


def get_brdc_data(file):
    sv = ['G02', 'G04', 'G05', 'G09', 'G10', 'G12', 'G17', 'G23', 'G25']
    brdc = read_rinex(file)
    brdc_data = []
    for i in range(9):
        brdc_data.append(get_sat_data(brdc, sv[i]).reset_index())
    return brdc_data


def read_icp(file, name):
    df = pd.read_csv(file , delimiter='\t', names = ['Time', '2', '3', '4', '5', '6', '7', name])
    df = df.drop(['2', '3', '4', '5', '6', '7'], axis=1)
    df = df.set_index('Time')
    df = df.loc[417136:417984]
    return df


def get_icp_data(type, sat):
    icp_data = pd.DataFrame()
    for i in range(9):
        file = 'data_' + type + '/icp_sat' + sat[i] + '.txt'
        name = type + ' pseudorange ' + sat[i]
        temp = read_icp(file, name)
        if i == 0:
            icp_data = temp
        else:
            icp_data = icp_data.join(temp, sort=True)
    return icp_data


def match_indices(base, rover):
    base.reset_index(inplace=True)
    base = base.drop(['Time'], axis=1)
    base.set_index(rover.index, inplace=True)
    return base


def get_ned_origin(file):
    df = pd.read_csv(file, delimiter='\t', names = ['Time', 'GPS_Week', 'X_ECEF', 'Y_ECEF', 'Z_ECEF',
                                                    'VX', 'VY', 'VZ', '1' , '2', '3', '4', '5', '6',
                                                    '7', '8', '9', '10'])
    df = df.drop(['GPS_Week', 'VX', 'VY', 'VZ', '1' , '2', '3', '4', '5', '6', '7', '8', '9', '10'], axis=1)
    x_base = df['X_ECEF'].mean()
    y_base = df['Y_ECEF'].mean()
    z_base = df['Z_ECEF'].mean()
    return [x_base, y_base, z_base]


def get_sv_position(brdc, rover, sat):
    index = [32, 32, 32, 32, 34, 32, 32, 32, 33]
    sat_pos = []
    for i in range(9):
        name = 'rover pseudorange ' + sat[i]
        sat_pos.append(calc_sv_position(rover[name], brdc[i].iloc[[index[i]]].to_dict('list')))
    return sat_pos


def calc_sv_position(rover_icp_data, br_data):
    # Constants
    U = 3.986005e14  # m^3/s^2 WGS84 value of Earth's Gravitational Constant
    OmegaDote = 7.2921151467e-5  # rad/s WGS84 Value of Earth's rotation rate
    To = 417600.000

    # DataFrame declaration
    temp = pd.DataFrame(index=rover_icp_data.index)
    # temp2 = pd.DataFrame(index=rover_icp_data.index, columns=['A'])
    # rover_icp_data = rover_icp_data.fillna(0)
    # for Time in temp.index:
    #     if rover_icp_data.loc[Time] != 0:
    #         temp2['A'].loc[Time] = br_data.get('sqrtA')[0]**2
    #     else:
    #         temp2['A'].loc[Time] = np.nan
    # Check for empty data
    rover_icp_data.fillna(0)

    # Calculated Columns with Logic
    temp['flag'] = np.where(rover_icp_data==0,False,True)
    temp['A'] = np.where(temp['flag'], br_data.get('sqrtA')[0]**2, np.nan)
    # temp['A'] = temp2['A']
    # Calculations
    temp['U/A^3'] = ((temp['A']**3)**-1) * U
    temp['sqrt(U/A^3)'] = np.sqrt(temp['U/A^3'].astype(float))
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
    sat_pos = sat_pos.fillna(0)
    return sat_pos


def convert_wgs_to_lla(base_vector):
    # Constants
    f = 1 / 298.257223563
    e = np.sqrt(f * (2 - f))
    Ro = 6378137
    Rp = Ro * (1 - f)

    # Calculations
    [x_base, y_base, z_base] = base_vector
    long = np.arctan2(y_base,x_base)
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
    lat = np.arctan((z_base + zo * ep ** 2) / p)

    # Rotation Matrix
    R = np.array([[-np.sin(lat) * np.cos(long), -np.sin(long), -np.cos(lat) * np.cos(long)],
                  [-np.cos(lat) * np.sin(long), np.cos(long), -np.cos(lat) * np.sin(long)],
                  [np.cos(lat), 0, -np.sin(lat)]])

    return [lat, long, h, R]


def convert_ecef_to_ned(x, y, z, R):
    x_ecef = np.array([[x], [y], [z]])
    x_temp = R.T @ (x_ecef)
    x_ned = [x_temp.item(0), x_temp.item(1), x_temp.item(2)]
    return x_ned


def calc_los_positions(sat_pos, ned_origin):
    [x_base, y_base, z_base] = ned_origin
    los_df = pd.DataFrame(index=sat_pos.index, columns=['x_los', 'y_los', 'z_los', 'elevation'])
    for i, series in sat_pos.iterrows():
        x = sat_pos['x'][i]
        y = sat_pos['y'][i]
        z = sat_pos['z'][i]
        if x != 0:
            magnitude = la.norm([x - x_base, y - y_base, z - z_base])
            x_los = (x - x_base) / magnitude
            y_los = (y - y_base) / magnitude
            z_los = (z - z_base) / magnitude
            los_df.loc[i, 'x_los'] = x_los
            los_df.loc[i, 'y_los'] = y_los
            los_df.loc[i, 'z_los'] = z_los
    return los_df


def calc_los_elevations(sat_pos, los_df, R):
    for i, series in sat_pos.iterrows():
        x_los = los_df['x_los'][i]
        y_los = los_df['y_los'][i]
        z_los = los_df['z_los'][i]
        if x_los != 0:
            x_ned = convert_ecef_to_ned(x_los, y_los, z_los, R)
            temp = m.asin(-1 * x_ned[2]) * (180 / m.pi)
            los_df.loc[i, 'elevation'] = temp


def get_los_positions(sat_pos, ned_origin):
    los_df = []
    for i in range(9):
        los_df.append(calc_los_positions(sat_pos[i], ned_origin))
    return los_df


def get_los_elevations(sat_pos, los_df, R):
    for i in range(9):
        calc_los_elevations(sat_pos[i], los_df[i], R)
    return los_df


def get_ref_sat(base_icp_data, los_df):
    high_sat = pd.DataFrame(index=los_df[0].index)
    high_sat = high_sat.join(los_df[0]['elevation'], sort=True, rsuffix='_0')
    high_sat = high_sat.join(los_df[1]['elevation'], sort=True, rsuffix='_1')
    ref = pd.DataFrame(index=los_df[0].index, columns=['Ref'])
    for Time in high_sat.index:
        if high_sat.loc[Time, 'elevation'] > high_sat.loc[Time, 'elevation_1']:
            ref.loc[Time, 'Ref'] = 0
        else:
            ref.loc[Time, 'Ref'] = 1
    base_icp_data = base_icp_data.join(ref, sort=True)
    return base_icp_data


def calc_dilution_of_precisions(pseudo_inv):
    global c
    sigma_x_sqr = pseudo_inv[0][0]
    sigma_y_sqr = pseudo_inv[1][1]
    sigma_z_sqr = pseudo_inv[2][2]

    VDOP = np.sqrt(sigma_z_sqr)
    HDOP = np.sqrt(sigma_x_sqr + sigma_y_sqr + sigma_z_sqr)
    DOP = [VDOP, HDOP]
    return DOP


def vector(x_ref):
    vector = []
    for i in range(3):
        vector.append(x_ref[i] / la.norm(x_ref))
    return vector


def calc_least_squares(data, H, rho, lat, long, h):
    rover_pos = pd.DataFrame(index=data.index, columns=['x','y','z','lat','long','alt'])
    pseudo_inv = []
    for i in range(848):
        p_inv = la.inv(H[i].T @ H[i])
        x_hat = p_inv @ H[i].T @ rho[i]
        NED = [float(x_hat[0]), float(x_hat[1]), float(x_hat[2])]
        rover_pos.iloc[i, 0] = float(x_hat[0])
        rover_pos.iloc[i, 1] = float(x_hat[1])
        rover_pos.iloc[i, 2] = float(x_hat[2])
        coord = navpy.ned2lla(NED, lat * 180 / m.pi, long * 180 / m.pi, h, latlon_unit='deg', alt_unit='m',
                                 model='wgs84')
        rover_pos.iloc[i, 3] = coord[0]
        rover_pos.iloc[i, 4] = coord[1]
        rover_pos.iloc[i, 5] = coord[2]
        pseudo_inv.append(p_inv)
    return [rover_pos, pseudo_inv]


def p_range_multi(base_icp_data, sat_pos, rover_icp_data, R, sat):
    H_list = []
    rho_list = []
    for Time in base_icp_data.index:
        rho = np.array([[]])
        H = np.array([[]])
        iter = 0
        flag = True
        x = sat_pos[base_icp_data.at[Time, 'Ref']].loc[Time, 'x']
        y = sat_pos[base_icp_data.at[Time, 'Ref']].loc[Time, 'y']
        z = sat_pos[base_icp_data.at[Time, 'Ref']].loc[Time, 'z']
        x_ref = convert_ecef_to_ned(x, y, z, R)
        vec = vector(x_ref)
        for i in sat:
            if sat_pos[i].at[Time, 'x'] != 0:
                x_temp = sat_pos[i].loc[Time, 'x']
                y_temp = sat_pos[i].loc[Time, 'y']
                z_temp = sat_pos[i].loc[Time, 'z']
                sat_temp = convert_ecef_to_ned(x_temp, y_temp, z_temp, R)
                vec_temp = vector(sat_temp)
                H_temp = np.array([vec[0] - vec_temp[0],
                              vec[1] - vec_temp[1],
                              vec[2] - vec_temp[2]])
                phi_rover_1 = rover_icp_data.iloc[iter, i]
                phi_base_1 = base_icp_data.iloc[iter, i]
                phi_rover_2 = rover_icp_data.iloc[iter, base_icp_data.at[Time, 'Ref']]
                phi_base_2 = rover_icp_data.iloc[iter, base_icp_data.at[Time, 'Ref']]
                rho_temp = (phi_rover_1 - phi_base_1) - (phi_rover_2 - phi_base_2)
                if flag:
                    H = H_temp
                    rho = np.array([[rho_temp]])
                    flag = False
                else:
                    H = np.vstack((H, H_temp))
                    rho = np.vstack((rho, rho_temp))
        iter+=1
        H_list.append(H)
        rho_list.append(rho)
    return [H_list, rho_list]


def get_least_squares(base_icp_data, sat_pos, rover_icp_data, R, lat, long, h):
    ind = [0, 2, 3, 4, 5, 6, 7, 8]
    [H, rho] = p_range_multi(base_icp_data, sat_pos, rover_icp_data, R, ind)
    return calc_least_squares(base_icp_data, H, rho, lat, long, h)



def mapping(longitude_min, longitude_max, latitude_min, latitude_max, longitude, latitude):
    image_box = ((longitude_min, longitude_max, latitude_min, latitude_max))
    fig, ax = plt.subplots(figsize = (8,7))
    ax.scatter(longitude, latitude, zorder=1, alpha=0.2, c ='b', s=10)
    ax.set_title('Plotting Spatial data on riyadh Map')
    ax.set_xlim(longitude_min, longitude_max)
    ax.set_ylim(latitude_min, latitude_max)


def main():
    sat = ['2', '4', '5', '9', '10', '12', '17', '23', '25']

    # Broadcast
    brdc_data = get_brdc_data('brdc2930.11n')

    # Data Base
    base_icp_data = get_icp_data('base', sat)

    # Data Rover
    rover_icp_data = get_icp_data('rover', sat)

    # Matching Indices
    base_icp_data = match_indices(base_icp_data, rover_icp_data)

    # Calculating Satellite Position
    sat_pos = get_sv_position(brdc_data, rover_icp_data, sat)

    # Use ecef to find ned origin
    ned_origin = get_ned_origin('data_base/ecef_rx0.txt')

    # Convert WGS to LLA
    [lat, long, h, R] = convert_wgs_to_lla(ned_origin)

    # Line of Sight
    los_df = get_los_positions(sat_pos, ned_origin)
    los_df = get_los_elevations(sat_pos, los_df, R)

    # Calculate Reference Satellite
    base_icp_data = get_ref_sat(base_icp_data, los_df)

    # Iterative Least Squares
    [rover_pos, pseudo_inv] = get_least_squares(base_icp_data, sat_pos, rover_icp_data, R, lat, long, h)

    # Dilutions of Precision
    DOP = calc_dilution_of_precisions(pseudo_inv)

    # Plotting
    plt.plot(rover_pos['x'], rover_pos['y'])
    plt.grid(b=None, which='major', axis='both')
    plt.xlabel('N')
    plt.ylabel('E')
    plt.title('Rover')
    plt.show()

    print('End')


main()

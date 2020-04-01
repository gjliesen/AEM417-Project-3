#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 11:55:11 2020

@author: jacobmontgomery

"""

import pandas as pd
import csv
import math as m
import numpy as np
from numpy import linalg as la
import navpy
from matplotlib import pyplot as plt

# %% Givens

mu = 3.986005 * 10 ** 14
c = 2.99792458 * 10 ** 8
omegadote = 7.2921151467 * 10 ** (-5)
pi = 3.1415926535898

# %% Data Parsing

i = 0
j = 0
current_sat = 0
sat_1_index = -1
sat_2_index = -1
sat_3_index = -1
sat_4_index = -1
sat_5_index = -1
sat_6_index = -1
sat_7_index = -1
sat_8_index = -1
sat_9_index = -1

columns = ['time1', 'axis1', 'i01', 'eccentricity1', 'omega1', 'OMEGA1', 'MO1',
           'IDOT1', 'OMEGADOT1', 'DEL_n1', 'Cuc1', 'Cus1', 'Crc1', 'Crs1',
           'Cic1', 'Cis1',
           'time2', 'axis2', 'i02', 'eccentricity2', 'omega2', 'OMEGA2',
           'MO2', 'IDOT2', 'OMEGADOT2', 'DEL_n2', 'Cuc2', 'Cus2', 'Crc2',
           'Crs2', 'Cic2', 'Cis2',
           'time3', 'axis3', 'i03', 'eccentricity3', 'omega3', 'OMEGA3', 'MO3',
           'IDOT3', 'OMEGADOT3', 'DEL_n3', 'Cuc3', 'Cus3', 'Crc3', 'Crs3',
           'Cic3', 'Cis3',
           'time4', 'axis4', 'i04', 'eccentricity4', 'omega4', 'OMEGA4', 'MO4',
           'IDOT4', 'OMEGADOT4', 'DEL_n4', 'Cuc4', 'Cus4', 'Crc4', 'Crs4',
           'Cic4', 'Cis4',
           'time5', 'axis5', 'i05', 'eccentricity5', 'omega5', 'OMEGA5', 'MO5',
           'IDOT5', 'OMEGADOT5', 'DEL_n5', 'Cuc5', 'Cus5', 'Crc5', 'Crs5',
           'Cic5', 'Cis5',
           'time6', 'axis6', 'i06', 'eccentricity6', 'omega6', 'OMEGA6', 'MO6',
           'IDOT6', 'OMEGADOT6', 'DEL_n6', 'Cuc6', 'Cus6', 'Crc6', 'Crs6',
           'Cic6', 'Cis6',
           'time7', 'axis7', 'i07', 'eccentricity7', 'omega7', 'OMEGA7', 'MO7',
           'IDOT7', 'OMEGADOT7', 'DEL_n7', 'Cuc7', 'Cus7', 'Crc7', 'Crs7',
           'Cic7', 'Cis7',
           'time8', 'axis8', 'i08', 'eccentricity8', 'omega8', 'OMEGA8', 'MO8',
           'IDOT8', 'OMEGADOT8', 'DEL_n8', 'Cuc8', 'Cus8', 'Crc8', 'Crs8',
           'Cic8', 'Cis8',
           'time9', 'axis9', 'i09', 'eccentricity9', 'omega9', 'OMEGA9', 'MO9',
           'IDOT9', 'OMEGADOT9', 'DEL_n9', 'Cuc9', 'Cus9', 'Crc9', 'Crs9',
           'Cic9', 'Cis9', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E8',
           'E9']

time_drifts_cols = ['bias1', 'drift1', 'drift_rate1',
                    'bias2', 'drift2', 'drift_rate2',
                    'bias3', 'drift3', 'drift_rate3',
                    'bias4', 'drift4', 'drift_rate4',
                    'bias5', 'drift5', 'drift_rate5',
                    'bias6', 'drift6', 'drift_rate6',
                    'bias7', 'drift7', 'drift_rate7',
                    'bias8', 'drift8', 'drift_rate8',
                    'bias9', 'drift9', 'drift_rate9']
index = range(14)

All_data = pd.DataFrame(index=index, columns=columns)
All_data = All_data.fillna(0)

with open('brdc2930.11n', newline='') as rinex:
    rinex_read = csv.reader(rinex, delimiter='\t')
    for lines in rinex_read:
        i = i + 1
        if i > 8:
            if j % 8 == 0:
                y = lines[0].split()
                if float(y[0].replace('D', 'e')) == 2:
                    sat_1_index = sat_1_index + 1
                    current_sat = 1
                    index = sat_1_index
                elif float(y[0].replace('D', 'e')) == 4:
                    sat_2_index = sat_2_index + 1
                    current_sat = 2
                    index = sat_2_index
                elif float(y[0].replace('D', 'e')) == 5:
                    sat_3_index = sat_3_index + 1
                    current_sat = 3
                    index = sat_3_index
                elif float(y[0].replace('D', 'e')) == 9:
                    sat_4_index = sat_4_index + 1
                    current_sat = 4
                    index = sat_4_index
                elif float(y[0].replace('D', 'e')) == 10:
                    sat_5_index = sat_5_index + 1
                    current_sat = 5
                    index = sat_5_index
                elif float(y[0].replace('D', 'e')) == 12:
                    sat_6_index = sat_6_index + 1
                    current_sat = 6
                    index = sat_6_index
                elif float(y[0].replace('D', 'e')) == 17:
                    sat_7_index = sat_7_index + 1
                    current_sat = 7
                    index = sat_7_index
                elif float(y[0].replace('D', 'e')) == 23:
                    sat_8_index = sat_8_index + 1
                    current_sat = 8
                    index = sat_8_index
                elif float(y[0].replace('D', 'e')) == 25:
                    sat_9_index = sat_9_index + 1
                    current_sat = 9
                    index = sat_9_index
                else:
                    current_sat == 0
                k = lines[0]

            if j % 8 == 1 and current_sat != 0:
                k = lines[0]
                y = [float(k[3:22].replace('D', 'e')),
                     float(k[22:41].replace('D', 'e')),
                     float(k[41:60].replace('D', 'e')),
                     float(k[60:79].replace('D', 'e'))]
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 13] = float(y[1])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 9] = float(y[2])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 6] = float(y[3])

            if j % 8 == 2 and current_sat != 0:
                k = lines[0]
                y = [float(k[3:22].replace('D', 'e')),
                     float(k[22:41].replace('D', 'e')),
                     float(k[41:60].replace('D', 'e')),
                     float(k[60:79].replace('D', 'e'))]
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 10] = float(y[0])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 3] = float(y[1])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 11] = float(y[2])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 1] = float(y[3] ** 2)

            if j % 8 == 3 and current_sat != 0:
                k = lines[0]
                y = [float(k[3:22].replace('D', 'e')),
                     float(k[22:41].replace('D', 'e')),
                     float(k[41:60].replace('D', 'e')),
                     float(k[60:79].replace('D', 'e'))]
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 0] = float(y[0])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 14] = float(y[1])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 5] = float(y[2])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 15] = float(y[3])

            if j % 8 == 4 and current_sat != 0:
                k = lines[0]
                y = [float(k[3:22].replace('D', 'e')),
                     float(k[22:41].replace('D', 'e')),
                     float(k[41:60].replace('D', 'e')),
                     float(k[60:79].replace('D', 'e'))]
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 2] = float(y[0])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 12] = float(y[1])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 4] = float(y[2])
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 8] = float(y[3])

            if j % 8 == 5 and current_sat != 0:
                k = lines[0]
                y = [float(k[3:22].replace('D', 'e')),
                     float(k[22:41].replace('D', 'e')),
                     float(k[41:60].replace('D', 'e')),
                     float(k[60:79].replace('D', 'e'))]
                All_data.iloc[index,
                              (current_sat - 1) * 16 + 7] = float(y[0])
                current_sat = 0

        j = j + 1


def icp_read_base(icp):
    with open('data_base/icp_sat' + str(icp) + '.txt', newline='') as sat:
        sat_read = csv.reader(sat, delimiter='\t')
        time = []
        psudoRange = []
        for sats in sat_read:
            time.append(round(float(sats[0])))
            psudoRange.append(float(sats[7]))
        sat_pd = pd.DataFrame(
            {'time': time, 'psudoRange': psudoRange})
        return sat_pd


def icp_read_rover(icp):
    with open('data_rover/icp_sat' + str(icp) + '.txt', newline='') as sat:
        sat_read = csv.reader(sat, delimiter='\t')
        time = []
        psudoRange = []
        for sats in sat_read:
            time.append(str(round(float(sats[0]))))
            psudoRange.append(float(sats[7]))
        sat_pd = pd.DataFrame({'time': time, 'psudoRange': psudoRange})
        return sat_pd


sat_1_base = icp_read_base(2)
sat_2_base = icp_read_base(4)
sat_3_base = icp_read_base(5)
sat_4_base = icp_read_base(9)
sat_5_base = icp_read_base(10)
sat_6_base = icp_read_base(12)
sat_7_base = icp_read_base(17)
sat_8_base = icp_read_base(23)
sat_9_base = icp_read_base(25)

sat_1_rover = icp_read_rover(2)
sat_2_rover = icp_read_rover(4)
sat_3_rover = icp_read_rover(5)
sat_4_rover = icp_read_rover(9)
sat_5_rover = icp_read_rover(10)
sat_6_rover = icp_read_rover(12)
sat_7_rover = icp_read_rover(17)
sat_8_rover = icp_read_rover(23)
sat_9_rover = icp_read_rover(25)

# %% Base Position
with open('data_base/ecef_rx0.txt', newline='') as sat:
    sat_read = csv.reader(sat, delimiter='\t')
    x_base = []
    y_base = []
    z_base = []
    for sats in sat_read:
        x_base.append(float(sats[2]))
        y_base.append(float(sats[3]))
        z_base.append(float(sats[4]))
    x_base = sum(x_base) / len(x_base)
    y_base = sum(y_base) / len(y_base)
    z_base = sum(z_base) / len(z_base)

# %% Putting Pseudoranges in all_data

cols3 = ['ref_sat', 'rover1', 'rover2', 'rover3', 'rover4', 'rover5', 'rover6', 'rover7', 'rover8', 'rover9',
         'base1', 'base2', 'base3', 'base4', 'base5', 'base6', 'base7', 'base8', 'base9']
index2 = range(417092, 417985)

Alll_data = pd.DataFrame(index=index2, columns=cols3)
Alll_data = Alll_data.fillna(0)


def psudo_put(pandas, sat):
    p = 0
    for times in pandas.time:
        Alll_data.loc[(int(times), 'rover' + str(sat))] = pandas.loc[(p, 'psudoRange')]
        p = p + 1


def psudo_put_base(pandas, sat):
    p = 0
    for times in pandas.time:
        Alll_data.loc[(int(times), 'base' + str(sat))] = pandas.loc[(p, 'psudoRange')]
        p = p + 1


psudo_put(sat_1_rover, 1)
psudo_put(sat_2_rover, 2)
psudo_put(sat_3_rover, 3)
psudo_put(sat_4_rover, 4)
psudo_put(sat_5_rover, 5)
psudo_put(sat_6_rover, 6)
psudo_put(sat_7_rover, 7)
psudo_put(sat_8_rover, 8)
psudo_put(sat_9_rover, 9)

psudo_put_base(sat_1_base, 1)
psudo_put_base(sat_2_base, 2)
psudo_put_base(sat_3_base, 3)
psudo_put_base(sat_4_base, 4)
psudo_put_base(sat_5_base, 5)
psudo_put_base(sat_6_base, 6)
psudo_put_base(sat_7_base, 7)
psudo_put_base(sat_8_base, 8)
psudo_put_base(sat_9_base, 9)

# %% Sat Positions

cols2 = ['x1', 'y1', 'z1', 'x2', 'y2', 'z2', 'x3', 'y3', 'z3',
         'x4', 'y4', 'z4', 'x5', 'y5', 'z5', 'x6', 'y6', 'z6',
         'x7', 'y7', 'z7', 'x8', 'y8', 'z8', 'x9', 'y9', 'z9']

sat_pos = pd.DataFrame(index=index2, columns=cols2)
sat_pos = sat_pos.fillna(0)

current_sat = 1
j = 0

# THIS IS THE ALGORITHM
def sat_post(i, current_sat):
    if current_sat in [1, 2, 3, 4, 6, 7, 8]:
        index = 11
    elif current_sat == 5:
        index = 12
    elif current_sat == 9:
        index = 9
    t0 = 417600
    a = All_data.loc[index, 'axis' + str(current_sat)]
    n = m.sqrt(mu / a ** 3) + All_data.loc[index, 'DEL_n' + str(current_sat)]
    tk = i - t0
    Mk = All_data.loc[index, 'MO' + str(current_sat)] + n * tk
    e = All_data.loc[index, 'eccentricity' + str(current_sat)]
    E = Mk
    tol = 10 ** (-8)
    ratio = .5
    while abs(ratio) > tol:
        E = E - ratio
        E_error = E - e * m.sin(E) - Mk
        E_derivative = 1 - e * m.cos(E)
        ratio = E_error / E_derivative
    All_data.loc[index, 'E' + str(current_sat)] = E
    num = (m.sqrt(1 - e ** 2) * m.sin(E)) / (1 - e * m.cos(E))
    den = (m.cos(E) - e) / (1 - e * m.cos(E))
    v = m.atan2(num, den)
    phi_k = v + All_data.loc[index, 'omega' + str(current_sat)]
    u = phi_k + All_data.loc[index, 'Cus' + str(current_sat)] * m.sin(2 * phi_k) \
        + All_data.loc[index, 'Cuc' + str(current_sat)] * m.cos(2 * phi_k)
    rk = a * (1 - e * m.cos(E)) + All_data.loc[index, 'Crs' + str(current_sat)] * \
         m.sin(2 * phi_k) + All_data.loc[index, 'Crc' + str(current_sat)] * \
         m.cos(2 * phi_k)
    i0 = All_data.loc[index, 'i0' + str(current_sat)]
    IDOT = All_data.loc[index, 'IDOT' + str(current_sat)]
    Cis = All_data.loc[index, 'Cis' + str(current_sat)]
    Cic = All_data.loc[index, 'Cic' + str(current_sat)]
    ik = i0 + IDOT * tk + Cis * m.sin(2 * phi_k) + Cic * m.cos(2 * phi_k)
    omega_dot = All_data.loc[index, 'OMEGADOT' + str(current_sat)]
    OMEGA = All_data.loc[index, 'OMEGA' + str(current_sat)]
    OMEGA_k = OMEGA + (omega_dot - omegadote) * tk - omegadote * t0
    x_p = rk * m.cos(u)
    y_p = rk * m.sin(u)

    x = x_p * m.cos(OMEGA_k) - y_p * m.cos(ik) * m.sin(OMEGA_k)
    y = x_p * m.sin(OMEGA_k) + y_p * m.cos(ik) * m.cos(OMEGA_k)
    z = y_p * m.sin(ik)

    sat_pos.loc[i, 'x' + str(current_sat)] = x
    sat_pos.loc[i, 'y' + str(current_sat)] = y
    sat_pos.loc[i, 'z' + str(current_sat)] = z


for times in sat_1_rover.time:
    if Alll_data.loc[(int(times), 'base1')] != 0:
        sat_post(int(times), 1)

for times in sat_2_rover.time:
    if Alll_data.loc[(int(times), 'base2')] != 0:
        sat_post(int(times), 2)

for times in sat_3_rover.time:
    if Alll_data.loc[(int(times), 'base3')] != 0:
        sat_post(int(times), 3)

for times in sat_4_rover.time:
    if Alll_data.loc[(int(times), 'base4')] != 0:
        sat_post(int(times), 4)

for times in sat_5_rover.time:
    if Alll_data.loc[(int(times), 'base5')] != 0:
        sat_post(int(times), 5)

for times in sat_6_rover.time:
    if Alll_data.loc[(int(times), 'base6')] != 0:
        sat_post(int(times), 6)

for times in sat_7_rover.time:
    if Alll_data.loc[(int(times), 'base7')] != 0:
        sat_post(int(times), 7)

for times in sat_8_rover.time:
    if Alll_data.loc[(int(times), 'base8')] != 0:
        sat_post(int(times), 8)

for times in sat_9_rover.time:
    if Alll_data.loc[(int(times), 'base9')] != 0:
        sat_post(int(times), 9)

# %% Convert to LLA
x = x_base
y = y_base
z = z_base
f = 1 / 298.257223563
e = m.sqrt(f * (2 - f))
R0 = 6378137
Rp = R0 * (1 - f)
long = m.atan2(y, x)
p = m.sqrt(x ** 2 + y ** 2)
E = m.sqrt(R0 ** 2 - Rp ** 2)
F = 54 * (Rp * z) ** 2
G = p ** 2 + (1 - e ** 2) * z ** 2 - (e * E) ** 2
c = (e ** 4 * F * p ** 2) / G ** 3
s = (1 + c + m.sqrt(c ** 2 + 2 * c)) ** (1 / 3)
P = (F / (3 * G ** 2)) / (s + 1 / s + 1) ** 2
Q = m.sqrt(1 + 2 * e ** 4 * P)
k1 = (-P * e ** 2 * p) / (1 + Q)
k2 = .5 * R0 ** 2 * (1 + 1 / Q)
k3 = -P * (1 - e ** 2) * (z ** 2 / (Q * (1 + Q)))
k4 = -.5 * P * p ** 2
k5 = p - e ** 2 * (k1 + m.sqrt(k2 + k3 + k4))
U = m.sqrt(k5 ** 2 + z ** 2)
V = m.sqrt(k5 ** 2 + (1 - e ** 2) * z ** 2)
h = U * (1 - (Rp ** 2 / (R0 * V)))
z0 = (Rp ** 2 * z) / (R0 * V)
ep = (R0 / Rp) * e
lat = m.atan((z + z0 * ep ** 2) / p)
R = np.array([[-m.sin(lat) * m.cos(long), -m.sin(long),
               -m.cos(lat) * m.cos(long)],
              [-m.sin(lat) * m.sin(long), m.cos(long),
               -m.cos(lat) * m.sin(long)],
              [m.cos(lat), 0, -m.sin(lat)]])


def ECEF_to_NED(x, y, z, R):
    vals = R.T @ np.array([[x], [y], [z]])
    valss = [vals.item(0), vals.item(1), vals.item(2)]
    return valss


# %% Line of Sight Vectors Base to Sats ECEF

Line_of_sight = pd.DataFrame(index=index2, columns=cols2)
Line_of_sight = Line_of_sight.fillna(0)
sats = [1, 2, 3, 4, 5, 6, 7, 8, 9]

for sat in sats:
    for times in sat_pos.index:
        x = sat_pos.loc[(times, 'x' + str(sat))]
        y = sat_pos.loc[(times, 'y' + str(sat))]
        z = sat_pos.loc[(times, 'z' + str(sat))]
        if x != 0:
            mag = la.norm([x - x_base, y - y_base, z - z_base])
            Line_of_sight.loc[(times, 'x' + str(sat))] = (x - x_base) / mag
            Line_of_sight.loc[(times, 'y' + str(sat))] = (y - y_base) / mag
            Line_of_sight.loc[(times, 'z' + str(sat))] = (z - z_base) / mag

# %% Elevation Angles
for times in sat_pos.index:
    elev = []
    sat_elev = []
    for sat in sats:
        x = Line_of_sight.loc[(times, 'x' + str(sat))]
        y = Line_of_sight.loc[(times, 'y' + str(sat))]
        z = Line_of_sight.loc[(times, 'z' + str(sat))]
        if x != 0:
            vals = ECEF_to_NED(x, y, z, R)
            elev.append((m.asin(-vals[2])) * 180 / m.pi)
            sat_elev.append(sat)
    if len(elev) > 3:
        Alll_data.loc[(times, 'ref_sat')] = sat_elev[elev.index(max(elev))]


# %% Constructing H Matrix
def unit_vector(vector):
    mag = la.norm(vector)
    vec = [0, 0, 0]
    vec[0] = vector[0] / mag
    vec[1] = vector[1] / mag
    vec[2] = vector[2] / mag
    return vec


Hs = []
psudos = []
ii = 0
j = -1
p = 0

for times in sat_pos.index:
    ii = ii + 1
    if Alll_data.loc[(times, 'ref_sat')] != 0:
        sat_options = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        ref_sat = Alll_data.loc[(times, 'ref_sat')]
        ref_sat = int(ref_sat)
        del sat_options[sat_options.index(ref_sat)]
        x = sat_pos.loc[(times, 'x' + str(ref_sat))]
        y = sat_pos.loc[(times, 'y' + str(ref_sat))]
        z = sat_pos.loc[(times, 'z' + str(ref_sat))]
        ref = ECEF_to_NED(x, y, z, R)
        ref_unit = unit_vector(ref)
        b = 0
        for sate in sat_options:
            if sat_pos.loc[(times, 'x' + str(sate))] != 0:
                x1 = sat_pos.loc[(times, 'x' + str(sate))]
                y1 = sat_pos.loc[(times, 'y' + str(sate))]
                z1 = sat_pos.loc[(times, 'z' + str(sate))]
                Cur_sat = ECEF_to_NED(x1, y1, z1, R)
                Cur_unit = unit_vector(Cur_sat)
                H = np.array([ref_unit[0] - Cur_unit[0],
                              ref_unit[1] - Cur_unit[1],
                              ref_unit[2] - Cur_unit[2]])
                psu_i_rov = Alll_data.loc[(times, 'rover' + str(sate))]
                psu_i_base = Alll_data.loc[(times, 'base' + str(sate))]
                psu_j_rov = Alll_data.loc[(times, 'rover' + str(ref_sat))]
                psu_j_base = Alll_data.loc[(times, 'base' + str(ref_sat))]
                psudo = (psu_i_rov - psu_i_base) - (psu_j_rov - psu_j_base)
                if b == 0:
                    Cur_H = H
                    psudo1 = np.array([[psudo]])
                    b = 1
                else:
                    Cur_H = np.vstack((Cur_H, H))
                    psudo1 = np.vstack((psudo1, psudo))
        Hs.append(Cur_H)
        psudos.append(psudo1)
        Cur_H = 0

# %% Least Squares
x = []
y = []
z = []
Lat = []
Long = []
alt = []
i = 0

while i < len(Hs):
    pos = la.inv(Hs[i].T @ Hs[i]) @ Hs[i].T @ psudos[i]
    x.append(float(pos[0]))
    y.append(float((pos[1])))
    z.append(float(pos[2]))
    NED = [float(pos[0]), float(pos[1]), float(pos[2])]
    LAT_LONG = navpy.ned2lla(NED, lat * 180 / m.pi, long * 180 / m.pi, h, latlon_unit='deg', alt_unit='m',
                             model='wgs84')
    Lat.append(LAT_LONG[0])
    Long.append(LAT_LONG[1])
    alt.append(LAT_LONG[2])
    i = i + 1

index5 = range(417098, 417985)

rover_pos = pd.DataFrame(index=index5, columns=['x', 'y'])
rover_pos = rover_pos.fillna(0)

i = 0

for times in rover_pos.index:
    rover_pos.loc[(times, 'x')] = x[i]
    rover_pos.loc[(times, 'y')] = y[i]
    i = i + 1

# %% For lat long plot

plt.plot(x, y)
plt.grid(b=None, which='major', axis='both')
plt.xlabel('N')
plt.ylabel('E')
plt.title('Rover')
# plt.ylim(min(Long), max(Long))
# plt.xlim(min(Lat), max(Lat))
plt.show()

# %% VDOP and HDOP

VDOP = []
HDOP = []

for H in Hs:
    mat = la.inv(H.T @ H)
    VDOP.append(m.sqrt(mat[2][2]))
    HDOP.append(m.sqrt(mat[1][1] + mat[0][0]))

plt.plot(VDOP)
plt.grid(b=None, which='major', axis='both')
plt.xlabel('Time')
plt.ylabel('VDOP')
plt.title('VDOP')
plt.show()

plt.plot(z)
plt.grid(b=None, which='major', axis='both')
plt.xlabel('Time')
plt.ylabel('Down')
plt.title('Down')
plt.show()

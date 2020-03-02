# Gregory Liesen
# AEM 417 Project 3

# Imports
import math as m
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import signal
from scipy import linalg as la
import georinex as gr

def read_rinex():
    df = gr.load('brdc2930.11n').to_dataframe()
    return df
def gen_least_sqr():
    print('Generalized Least Squares Position')


def vdop():
    print('VDOP')


def hdop():
    print('HDOP')

def pdop():
    print('PDOP')

def main():
    brdc_df = read_rinex()
    print('Done')
main()

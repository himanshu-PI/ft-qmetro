import stim
import numpy as np
import pymatching
import sinter
from typing import List
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import matplotlib
matplotlib.use("TkAgg")  # or "Qt5Agg"


import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['text.color'] = 'k'
plt.rcParams['axes.labelcolor'] = 'k'
plt.rcParams['text.usetex'] = False

pal = ["003049","d62828","f77f00","fcbf49","eae2b7"]
pal = ['#' + i for i in pal]

marks = [
    '.',  # point
    ',',  # pixel
    'o',  # circle
    'v',  # triangle_down
    '^',  # triangle_up
    '<',  # triangle_left
    '>',  # triangle_right
    '1',  # tri_down
    '2',  # tri_up
    '3',  # tri_left
    '4',  # tri_right
    '8',  # octagon
    's',  # square
    'p',  # pentagon
    'P',  # plus (filled)
    '*',  # star
    'h',  # hexagon1
    'H',  # hexagon2
    '+',  # plus
    'x',  # x
    'X',  # x (filled)
    'D',  # diamond
    'd',  # thin_diamond
    '|',  # vline
    '_',  # hline
]

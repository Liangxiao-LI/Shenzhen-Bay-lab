#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 20:52:44 2023

@author: ryan
"""
import pandas as pd

df = pd.read_excel('/Users/ryan/Documents/GitHub/Shenzhen Bay lab/spike_protein_3structures_markerd_with_variants_from_cncb/6vxx_variants.xls')

import pysal as py
import numpy as np

coordinates = df[['X', 'Y', 'Z']].values

w = pysal.weights.DistanceBand(coordinates, threshold=100, binary=False)

import libpysal
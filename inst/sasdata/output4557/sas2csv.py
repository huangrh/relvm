#!/usr/bin/env python3

import pandas as pd
from glob import glob
import os
from os.path import basename

path="./"
filenames = glob(path + "/*.sas7bdat")

for filename in filenames:
    csvbase = basename(filename)
    csvbase = os.path.splitext(csvbase)[0]
    csvfile = ".".join([csvbase,"csv"])
    # print(csvfile)
    dat = pd.read_sas(filename)
    dat.to_csv(csvfile,index=False)

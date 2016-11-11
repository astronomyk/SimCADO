"""
A code snippet to generate the Data Listings document
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import openpyxl
from astropy.io import ascii, fits

home_dir = "./"
plot_dir = "./source/images/"
site_plot_dir = "./images/"
data_dir = "../data/"
data_ext_dir = "../../data_ext/"

data = openpyxl.load_workbook(home_dir + "SimCADO_defaults.xlsx").worksheets[0]

text = """
# SimCADO configuration for MICADO

SimCADO's main purpose is to simulate the optical train comprising of the E-ELT and MICADO. As such the base configuration is for this instrument/telescope combination.

Although SimCADO allows the user to configure the optical train as they like, it is *highly* recommended to leave the default parameters alone if the goal of the simulations is comparability between scientific use cases.

This document lists the current configuration which best describes MICADO. Each section describes a physical effect included by SimCADO and **gives the Keyword/Value pair** contained in the default configuration file **for MICADO**.

**Data Source**

The data for this table is generated from this [Excel sheet](https://github.com/gastronomyk/SimCADO/blob/master/SimCADO/simcado/docs/SimCADO_defaults.xlsx)

"""

j = 1
while data.rows[j][0].value != "EOF":
    # Section
    if data.rows[j][0].value is not None:
        text += "## " + str(data.rows[j][0].value) + "\n\n"
    
    # Subsection
    if data.rows[j][1].value is not None:
        text += "------ \n\n"
        text += "### " + str(data.rows[j][1].value) + "\n\n"
    
    # Description of source
    if data.rows[j][7].value is not None:
        text += "**Description:** " + str(data.rows[j][7].value) + "\n\n "

    # Keyword-default value pair
    if data.rows[j][3].value is not None and data.rows[j][4].value is not None:
        text += "`" + str(data.rows[j][3].value) + \
                " = " + str(data.rows[j][4].value) + "`          "
        
    # Source and last update date of data, with comments
    if data.rows[j][5].value is not None and data.rows[j][6].value is not None:
        text += "Source: " + str(data.rows[j][5].value) + \
                ", updated on " + str(data.rows[j][6].value)
        if data.rows[j][8].value is not None:
            text += ", **" + str(data.rows[j][8].value) + "**"
    text += "\n\n"
    
    
    # if the default is a transmission curve, generate a plot
    if data.rows[j][9].value is not None:
        fname = os.path.split(data.rows[j][4].value)[1]
        
        if data.rows[j][9].value == "y":
            print(data.rows[j][4].value)
            tmp = ascii.read(os.path.join(data_dir, fname))
            plt.figure(figsize=(6,3))
            plt.plot(tmp[tmp.colnames[0]].data, tmp[tmp.colnames[1]].data)
            plt.xlabel("Wavelength [um]")
            plt.xlim(0.7,2.5)
            plt.title(fname + "\n" + str(data.rows[j][6].value))
            plt.tight_layout()

        elif data.rows[j][9].value in ("ec", "tc"):
            #print(data.rows[j][4].value)
            print(os.path.join(data_dir, fname))
            tmp = fits.getdata(os.path.join(data_dir, fname))
            plt.figure(figsize=(6,3))
            if data.rows[j][9].value == "ec":
                plt.plot(tmp["lam"].data, tmp["flux"].data)
                fname += "[emission]"
            elif data.rows[j][9].value == "tc":
                plt.plot(tmp["lam"].data, tmp["trans"].data)
                fname += "[transmission]"
            plt.xlabel("Wavelength [um]")
            plt.xlim(0.7,2.5)
            plt.semilogy()
            plt.ylim(ymin=0.1)
            plt.title(fname + "\n" + str(data.rows[j][6].value))
            plt.tight_layout()

        elif data.rows[j][9].value == "im":
            print(data.rows[j][4].value, fname)
            tmp = fits.getdata(os.path.join(data_ext_dir, fname))
            w, h = tmp.shape
            n=32
            fig, sp = plt.subplots(1, 2, figsize=(6,4))
            sp[0].imshow(tmp[w//2-n:w//2+n, h//2-n:h//2+n], norm=LogNorm(), cmap="hot")
            sp[1].imshow(tmp, norm=LogNorm(), cmap="hot", vmax=0.01*np.max(tmp))
            #sp[1].colorbar()
            fig.suptitle(fname + "\n" + str(data.rows[j][6].value))
            fig.tight_layout()

        hname = os.path.join(plot_dir, fname+".png")
        plt.savefig(format="png", filename=hname)

        gname = os.path.join(site_plot_dir, fname+".png")
        text += "![" + gname + "](" + gname + ") \n\n"
        print(hname, gname)
    j += 1

f = open("./source/SimCADO_defaults.md", "w")
f.write(text) 
f.close()
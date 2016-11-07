import os
import matplotlib.pyplot as plt
import openpyxl
from astropy.io import ascii

home_dir = "./"
plot_dir = "./source/images/"
site_plot_dir = "./images/"
data_dir = "../data/"

data = openpyxl.load_workbook(home_dir + "SimCADO_defaults.xlsx").worksheets[0]

text = """
# SimCADO configuration for MICADO

SimCADO's main purpose is to simulate the optical train comprising of the E-ELT and MICADO. As such the base configuration is for this instrument/telescope combination.

Although SimCADO allows the user to configure the optical train as they like, it is *highly* recommended to leave the default parameters alone if the goal of the simulations is comparability between scientific use cases.

This document lists the current configuration which best describes MICADO. Each section describes a physical effect included by SimCADO and **gives the Keyword/Value pair** contained in the default configuration file **for MICADO**.

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
        
    # Keyword-default value pair
    if data.rows[j][3].value is not None and data.rows[j][4].value is not None:
        text += "`" + str(data.rows[j][3].value) + " = " + str(data.rows[j][4].value) + "`          "
        
        # Source and last update date of data, with comments
        if data.rows[j][5].value is not None and data.rows[j][6].value is not None:
            text += "Source: " + str(data.rows[j][5].value) + ", updated on " + str(data.rows[j][6].value)
            if data.rows[j][8].value is not None:
                text += ", **" + str(data.rows[j][8].value) + "**"
        text += "\n\n"
        
        
        # if the default is a transmission curve, generate a plot
        if data.rows[j][9].value is not None and "y" in data.rows[j][9].value:
            fname = os.path.split(data.rows[j][4].value)[1]
            tmp = ascii.read(os.path.join(data_dir, fname))
            plt.figure(figsize=(6,3))
            plt.plot(tmp[tmp.colnames[0]].data, tmp[tmp.colnames[1]].data)
            plt.xlabel("Wavelength [um]")
            plt.title(fname + "\n" + str(data.rows[j][6].value))
            plt.tight_layout()
            
            hname = os.path.join(plot_dir, fname+".png")
            plt.savefig(format="png", filename=hname)

            gname = os.path.join(site_plot_dir, fname+".png")
            text += "![" + gname + "](" + gname + ") \n\n"
    
    # Description of source
    if data.rows[j][7].value is not None:
        text += "**Description:** " + str(data.rows[j][7].value) + "\n\n "

    j += 1

f = open("./source/SimCADO_defaults.md", "w")
f.write(text) 
f.close()
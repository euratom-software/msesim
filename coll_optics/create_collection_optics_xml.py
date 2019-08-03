import xml.etree.ElementTree
import numpy as np
from xml.dom import minidom

import matplotlib.pyplot as plt

# Open the collection optics xml file template:
tree = xml.etree.ElementTree.parse('xml_template.xml')
root = tree.getroot()

#width and height of whole fiber array in mm
fibre_height = 37.5
fibre_width =  7.8
n_channels = 25
n_fibers = 6 #number fibers per channel

x_length = np.linspace(-1*fibre_width/2, fibre_width/2, n_fibers)
y_length = np.linspace(-1*fibre_height/2,fibre_height/2, n_channels)

fname = "JET_MSE_collection_optics.xml"

grid_x, grid_y = np.meshgrid(x_length, y_length)

pixel_names = []
attrib_x = dict()
attrib_y = dict()
coords = []

for i in range(len(grid_x[:,0])):
    for j in range(len(grid_y[0,:])):
        coords.append(([grid_x[i,j],grid_y[i,j]]))


names = []

for i in range(len(coords)):
    names.append('p' + str([i][0]))

xml.etree.ElementTree.SubElement(tree.find('coll'), 'bundleID',
                                     attrib={'description': '', 'type': 'string', 'value':', '.join(names)})

for i in range(len(coords)):
    xml.etree.ElementTree.SubElement(tree.find('coll'), 'p'+str([i][0])+'l', attrib={'description':'', 'type':'str', 'value':str(coords[i][0])})
    xml.etree.ElementTree.SubElement(tree.find('coll'), 'p'+str([i][0])+'m', attrib={'description':'', 'type':'str', 'value':str(coords[i][1])})

xmlstr = minidom.parseString(xml.etree.ElementTree.tostring(root)).toprettyxml(indent="   ")
with open(fname, "w") as f:
    f.write(xmlstr)


# # #input the sensor details:
# nx = int(input('Number of x pixels in sensor array?'))
# ny = int(input('Number of y pixels in sensor array?'))
# pixel_size = float(input('Pixel size in mm?'))

# x_length = np.arange(-int(nx/2), int(nx/2)+0.5, 1)*pixel_size
# y_length = np.arange(-int(ny/2), int(ny/2)+0.5, 1)*pixel_size

#Simulate every nth pixel:
# x_length = x_length[0::50]
# y_length = y_length[0::50]
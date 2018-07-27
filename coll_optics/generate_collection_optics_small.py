import xml.etree.ElementTree
import numpy as np
from xml.dom import minidom

# Open the collection optics xml file template:
tree = xml.etree.ElementTree.parse('xml_template.xml')
root = tree.getroot()

x_length = np.linspace(-10.23,10.23,32)
y_length = np.linspace(-10.23,10.23,32)

fname = "imse_2d_32x32_centerpixel.xml"

grid_x, grid_y = np.meshgrid(x_length, y_length)

# grid_x = grid_x[::32]
# grid_y = grid_y[::32]

pixel_names = []
attrib_x = dict()
attrib_y = dict()
pixel_number = np.arange(0, len(x_length) ** 2 + 1, 1)
coords = []

for i in range(len(grid_x)):
    for j in range(len(grid_y)):
        coords.append(([grid_x[i,j], grid_y[i,j]]))

names = []

for i in range(len(coords)):
    names.append('p' + str([i][0]))

xml.etree.ElementTree.SubElement(tree.find('coll'), 'bundleID', attrib={'description': '', 'type': 'string', 'value': ', '.join(names)})

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
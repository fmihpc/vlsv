import ast
import hashlib
import numpy as np
import struct
import sys
import xml.etree.ElementTree as ET
import matplotlib
import matplotlib.pyplot as plt

### User parameters ###
scale_size = 1e6 # scaling of file sizes
cmap = matplotlib.cm.hsv # colormap to use
### End of User parameters ###

if (len(sys.argv) != 2):
   print "This script takes a vlsv file name as input!"
   sys.exit()

file_name = sys.argv[1]

max_xml_size = 1000000

fptr = open(file_name,"rb")
# Eight first bytes indicate whether the system is big_endianness or something else
endianness_offset = 8
fptr.seek(endianness_offset)
# Read 8 bytes as unsigned long long (uint64_t in this case) after endianness, this tells the offset of the XML file.
uint64_byte_amount = 8
(offset,) = struct.unpack("Q", fptr.read(uint64_byte_amount))
# Move to the xml offset
fptr.seek(offset)
# Read the xml data
xml_data = fptr.read(max_xml_size)
# Read the xml as string
(xml_string,) = struct.unpack("%ds" % len(xml_data), xml_data)
# Input the xml data into xml_root
xml_root = ET.fromstring(xml_string)

fptr.close()


positions = []
names = []
sizes = []

for child in xml_root:
   positions.append(ast.literal_eval(child.text) / scale_size)
   size = ast.literal_eval(child.attrib["vectorsize"]) * ast.literal_eval(child.attrib["arraysize"]) * ast.literal_eval(child.attrib["datasize"]) / scale_size
   sizes.append(size)
   names.append(child.tag)
   if ("name" in child.attrib):
      names[-1] = names[-1]+"_"+child.attrib["name"]
   names[-1] = names[-1]+"_"+"%.1e" % size

positions = np.array(positions)
sizes = np.array(sizes)
names = np.array(names)

positions_sorted = positions[positions.argsort()]
sizes_sorted = sizes[positions.argsort()]
names_sorted = names[positions.argsort()]

total_size = np.sum(sizes)

colors=[int(hashlib.md5(i).hexdigest(), 16) for i in names_sorted]
norm = matplotlib.colors.Normalize(vmin=min(colors), vmax=max(colors))
m = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
colors = m.to_rgba(colors)

fig = plt.figure()
plt.figtext(0.5, 0.95, file_name+" memory layout in MB (total size %.1e MB)" % total_size , fontsize=36, ha="center")

plt.subplot(121)
ax = fig.gca()

#left=0
bottom=0

for i,val in enumerate(positions_sorted):
#   patch = ax.barh(0, sizes_sorted[i], height=1.0, left=val, color=colors[i])
   patch = ax.bar(0, sizes_sorted[i], width=1.0, bottom=val, color=colors[i])
   
   coords = patch[0].get_xy()
   x = 1.05*patch[0].get_width() + coords[0]
   y = 0.5*patch[0].get_height() + coords[1]
   
   ax.text(x, y, names_sorted[i], rotation=0, ha="left", va="center")

#ax.set_aspect(0.05*np.sum(sizes))
ax.set_aspect(1.0/(0.05*total_size))


plt.subplot(122)
ax = fig.gca()

patches, texts = ax.pie(sizes_sorted, labels=names_sorted, explode=0.0*np.ones(len(sizes)), colors=colors)
for i,patch in enumerate(patches):
   texts[i].set_horizontalalignment("left")
   texts[i].set_verticalalignment("center")
   texts[i].set_rotation_mode("anchor")
   angle = 0.5*(patch.theta1 + patch.theta2)
   if (angle > 90.0 and angle < 270.0):
      angle = angle - 180.0
      texts[i].set_horizontalalignment("right")
   texts[i].set_rotation(angle)
ax.axis('equal')


plt.show()

# A script that populates a MarkDown file to contain all the images

from elements_to_simulate import elements
import openmc


# TODO replace the "<committers>" tag with the github committers found with https://github.com/PyGithub/PyGithub

# Read in the template file
with open('report_template.md', 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('<openmc_version>', openmc.__version__)
filedata = filedata.replace('<nuc_data>', 'ENDF/B 8.0') # TODO replace with user arg


string_to_add = ''
for element in elements:
  string_to_add += f'![{element} atoms](figs/atoms/element_{element}.png)\n'
  string_to_add += f'![{element} activity](figs/activity/element_{element}.png)\n'

filedata = filedata.replace('<elements>', string_to_add)

# Write the file out again
with open('report.md', 'w') as file:
  file.write(filedata)

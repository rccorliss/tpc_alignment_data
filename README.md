# tpc_alignment_data
raw files for TPC alignment

Going clockwise from the chalice hole on petal 54, the petals in use are:
north, going clockwise from the special 54-CH:  54-CH,58,69,55-CH,60,61,50-CH,62,63,51-CH,64,59,52-CH,65,66,53-CH,67,57
south, going clockwise from the special 54-CH:  78-CH,48,76,92-CH,75,88,80-CH,86,85,93-CH,45,68,79-CH,87,44,83-CH,74,46


OGP contains the data from the OGP scans of the central membrane petals.
.DAT files are detailed outlines, in global coordinates, of all stripes on the petal, as well as the four arcs that are the visible portions of the survey marks
.txt files are a separate survey of the centers and areas of each petal.  The first four places surveyed in each pattern are the survey marks themselves, which are silver circles with black outlines.  The circles are irregular on this scale.

survey contains the various survey data files.
survey/PETAL 1 COORDS-mm.txt   :   the top face of the CM as lifted into the air on the assembly table
survey/PETAL BACKSIDE COORDS-mm.txt   :  the bottom face of the CM as lifted into the air on the assembly table

tpc_magnet_axis_calculation.C -- finds the orientation of the tpc wrt the magnet, for purposes of rotating the magnetic field map into annular fieldsim local coordinates.
ogp_to_tree.C -- converts the edges.DAT files into .ogp_tree.root files, which contains a ntuple with petal#, stripe#, and tvectors for each stripe
opc_to_tree.C -- functions to determine area and geometric centers.


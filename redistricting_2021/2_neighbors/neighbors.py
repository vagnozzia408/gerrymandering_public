#!/usr/bin/python3

# To run in the terminal: ./scriptname.py input_clockwise_file.csv output_clockwise_cleaned.csv output_neighbors_file.csv

# Anna Marie Vagnozzi
# Last Updated: November 6, 2021
# Written in Python3

# Any items marked as "TO DO" are things to check and adjust as necessary before running script on a new data set.

###################################################################################################

# Import necessary libraries.
import io
import csv
import sys

# Open the CSV file containing the clockwise point data.
# This is created in ArcGIS Pro by Tabulating the Intersection of precinct lines and point features.

# TO DO: Set the filepath for the CSV file.
infile = io.open(sys.argv[1],newline='')
reader = csv.reader(infile)
header = next(reader)

# Read in SRC_PSN (line), NBR_PSN (line), Shared_Perim (line), PSN (point), CLCKWS_IDX (point).
# This step may take a while.
# TO DO: Adjust column numbers as needed.
clockwise_data = [[eval(row[1]),eval(row[2]),eval(row[3]),eval(row[4]),eval(row[5])] for row in reader]
infile.close()

# TO DO: Make sure the total number of records matches the files generated in ArcMap.
assert(len(clockwise_data)==2728766)

## TEST ##
print("TEST: Clockwise precinct data read correctly...")
print("SRC_PSN, NBR_PSN, Shape_Length, PSN, CLCKWS_IDX")
for i in range(5):
    print(clockwise_data[i])
print("")

# We only care about instances where SRC_PSN and PSN are the same.
# SRC_PSN is the source precinct with respect to the line segment.
# PSN is the source precinct with respect to the point.

# Clean the data frame to only keep instances where SRC_PSN == PSN.
clockwise_data_cleaned = list()

for line in clockwise_data:
    # If the source precincts match...
    if line[0] == line[3]:
        # Include SRC_PSN, NBR_PSN, Shared_Perim, and CLCKWS_IDX in one line of data.
        clockwise_data_cleaned.append([line[0],line[1],line[2],line[4]])
		
## TEST ##
print("TEST: Clockwise data cleaned properly...")
print("SRC_PSN, NBR_PSN, Shared_Perim, CLCKWS_IDX")
for i in range(5):
    print(clockwise_data_cleaned[i])
print("")

# Write the cleaned data to a CSV file for easy reference.
# TO DO: Set the appropriate filepath.
outfile1 = open(sys.argv[2],'w',newline='') #Python 3
#outfile1 = open(sys.argv[2],'wb')
writer1 = csv.writer(outfile1)
writer1.writerow(["SRC_PSN","NBR_PSN","Shared_Perim","CLCKWS_IDX"])

for row in clockwise_data_cleaned:
    writer1.writerow(row)
    
outfile1.close()

# Create a dictionary of clockwise points along each line segment.
# Keys: (SRC_PSN, NBR_PSN, Shared_Perim) - uniquely identifies a line segment.
# Values: [CLCKWS_IDX 1, CLCKWS_IDX 2, ...] - list of clockwise indices for the points along the line segment.

# The length of each line segment is also the length of the shared perimeter between SRC_PSN and NBR_PSN.
# Points are clockwise with respect to the source precinct.
# Note that due to the way ArcGIS generates the file, "node neighbors" (neighbors only at a point) are not included.

clockwise_dict = dict()

for line in clockwise_data_cleaned:
    SRC_PSN = line[0]
    NBR_PSN = line[1]
    Shared_Perim = line[2]
    CLCKWS_IDX = line[3]
    
    # Create the unique line segment key.
    line_seg = (SRC_PSN, NBR_PSN, Shared_Perim)
    
    # If the line segment is already stored in the dictionary...
    if line_seg in clockwise_dict:
        # Add the clockwise ID in this line to the key's value list.
        clockwise_dict[line_seg].append(CLCKWS_IDX)
    # If this is the first time we've seen this line segment...
    # i.e. if line_seg not in clockwise_dict:
    else:
        # Add the line segment key to the dictionary and add the clockwise ID to its value list.
        clockwise_dict[line_seg] = [CLCKWS_IDX]
        
# TO DO: Set this assertion statement to the number of line segments. (See number of records in ArcGIS lines feature.)
# If an error occurs, it's likely due to a mapping error in the shapefile boundaries, implying that there are gaps, slivers, or overlaps.
#print("Number of Line Segments?")
#print(len(clockwise_dict))
assert(len(clockwise_dict)==13532)

## TEST ##
# Note: clockwise_dict.items() returns a dict_items object with (key,value) pairs
print("TEST: Clockwise dictionary created properly...")
print("Line Segment: List of Clockwise Indices")
print("")
for i in range(5):
    print(str(list(clockwise_dict.items())[i][0]) + ": " + str(list(clockwise_dict.items())[i][1]))
    print("")

# Ensure that clockwise indices are sorted in increasing order for each neighbor.
for line_seg in clockwise_dict:
    clockwise_dict[line_seg] = sorted(clockwise_dict[line_seg])

# "Problem segments" in this algorithm are those with only two vertices which are NOT sequential.
# Identify any problem segments. If nothing happens, everything is fine and the algorithm will work.

problem_segments = list()

for line_seg in clockwise_dict:
    # If there are only two vertices AND they are not sequential...
    if len(clockwise_dict[line_seg])==2 and (clockwise_dict[line_seg][1]-clockwise_dict[line_seg][0]) != 1:
        # Keep a record of this problem segment.
        problem_segments.append(line_seg)

assert(len(problem_segments)==0)

# Check that EITHER of the following (but not both) is true:
# First and last indices are off by more than one
# Last two indices are off by more than one

# I haven't actually proved this, but I have reason to believe that these
# two clockwise indices are duplicates of endpoints and one can be deleted.

for line_seg in clockwise_dict:
    idx_list = clockwise_dict[line_seg]
    n = len(idx_list)
    if idx_list[1]-idx_list[0]!=1 and idx_list[n-1]-idx_list[n-2]!=1:
        print(str(line_seg) + ": " + str(idx_list))
		
# If nothing printed, then delete any duplicate endpoints.

for line_seg in clockwise_dict:
    # Store the list of clockwise indices for the given line segment.
    index_list = clockwise_dict[line_seg]
    # Store the number of clockwise indices for the line segment.
    n = len(index_list)
    
    # If the first two indices are not sequential...
    if index_list[1] - index_list[0] != 1:
        # The first one is a duplicate point and can be deleted.
        del index_list[0]
        assert(len(index_list) == n-1)
    # If the last two indices are not sequential...
    elif index_list[n-1] - index_list[n-2] != 1:
        # The last one is a duplicate point and can be deleted.
        del index_list[n-1]
        assert(len(index_list) == n-1)
    # Otherwise, everything is good to go...
    else:
        # So do nothing!
        pass
		
## TEST ##
# We can test to make sure that there are no weird endpoints.
#for item in sorted(clockwise_dict.items()):
#    print(str(item[0]) + ": " + str(item[1]))
#    print("")

# Now we need the endpoints for each line segment so we can match them.
# Keys: SRC_PSN - source precincts
# Values: (NBR_PSN, Shared_Perim, First_Endpoint, Last_Endpoint)

endpoint_dict = dict()

for line_seg in clockwise_dict:
    SRC_PSN = line_seg[0]
    NBR_PSN = line_seg[1]
    perim = line_seg[2]
    
    idx_list = clockwise_dict[line_seg]
    
    first_endpt = idx_list[0]
    last_endpt = idx_list[-1]
    
    assert(last_endpt>first_endpt)
    
    if SRC_PSN in endpoint_dict:
        endpoint_dict[SRC_PSN].append((NBR_PSN, perim, first_endpt, last_endpt))
    else:
        endpoint_dict[SRC_PSN] = [(NBR_PSN, perim, first_endpt, last_endpt)]

## TEST ##
print("TEST: Endpoint dictionary created properly...")
print("Source Precinct: List of Neighboring Segments")
print("")
for i in range(5):
    print(str(list(endpoint_dict.items())[i][0]) + ": " + str(list(endpoint_dict.items())[i][1]))
    print("")
	
## TEST ##
print("Number of line segments around the border of the state: " + str(len(endpoint_dict[-1])))
print("")

## TEST ##
# We can verify that the line segments appear in clockwise order around the state boundary.
# for each in sorted(endpoint_dict[-1], key = lambda tup: tup[2], reverse = True):
#     print(each)
	
# Generate the neighbor lists in clockwise order.
# Create the neighbor_lists data frame.

all_neighbors = list()
boundary_neighbors = list()

# For each source precinct:
for SRC_PSN in endpoint_dict:
    # For the outer neighbors, we need to reverse the order.
    if SRC_PSN == -1:
        # Sort the potential neighbors into a list, sorted in reverse order.
        border_precincts = sorted(endpoint_dict[SRC_PSN], key = lambda nbr: nbr[2], reverse = True)
        
        # Make sure that the endpoints match sequentially.
        for i in range(len(border_precincts)-1):
            assert(border_precincts[i][2] == border_precincts[i+1][3])
            
        all_boundary_precincts = list()
        
        for bdr_pct in border_precincts:
            pct = bdr_pct[0]
            perim = bdr_pct[1]
            
            # If this is the first boundary precinct:
            if len(all_boundary_precincts) == 0:
                all_boundary_precincts.append([pct,perim])
            # If it isn't the first neighbor we've seen...
            else:
                # If we have two of the same border precinct in a row...
                # (Look at the last item in the list of boundary precincts...)
                if all_boundary_precincts[-1][0] == pct:
                    assert(all_boundary_precincts[-1][1] != perim)
                    # Then we have two adjacent line segments along the state boundary for a single precinct.
                    # Don't change the neighbor, but add the shared perimeter.
                    all_boundary_precincts[-1][1] = all_boundary_precincts[-1][1] + perim
                # But if it's a new neighbor...
                else:
                    all_boundary_precincts.append([pct,perim])
        print("New number of line segments around the state border: " + str(len(all_boundary_precincts)))
        # Note: This includes multi-adjacencies, so this number is greater than the number of border precincts.
                    
    # For all other neighbors:
    else:
        # Sort the potential neighbors into a list, sorted by their starting endpoint.
        sorted_neighbors = sorted(endpoint_dict[SRC_PSN], key = lambda nbr: nbr[2])
        
        # Make sure the endpoints match sequentially.
        for i in range(len(sorted_neighbors)-1):
            # If this assertion error fails, it implies an issue with the vertices. Use the lines of code below to pinpoint the issue. 
            # An error may happen if two polygons intersect at a point, and the bordering polygon is missing a vertex at that point.
            assert(sorted_neighbors[i][3] == sorted_neighbors[i+1][2])
            #assert(sorted_neighbors[i][3] == sorted_neighbors[i+1][2]),("Errors: "+str(sorted_neighbors[i])+" and "+str(sorted_neighbors[i+1]))
            # "Errors: "+str(sorted_neighbors[i])+" and "+str(sorted_neighbors[i+1])
            #if sorted_neighbors[i][3] != sorted_neighbors[i+1][2]:
                #print(str(sorted_neighbors[i])+" and "+str(sorted_neighbors[i+1]))
            
        pending_neighbor_list = list()
        pending_perim_list = list()
        
        for neighbor in sorted_neighbors:
            nbr = neighbor[0]
            perim = neighbor[1]
            
            # If this is the first neighbor we come across...
            if len(pending_neighbor_list) == 0 and len(pending_perim_list) == 0:
                pending_neighbor_list.append(nbr)
                pending_perim_list.append(perim)
            # If it isn't the first neighbor we've seen...
            else:
                # If we have two of the same neighbor in a row...
                if pending_neighbor_list[-1] == nbr:
                    assert(pending_perim_list[-1] != perim)
                    # Then we have two adjacent line segments along the same neighbor boundary.
                    # Don't change the pending_neighbor_list.
                    # Add the shared perimeter to the previous one.
                    pending_perim_list[-1] = pending_perim_list[-1] + perim
                # But if it's a new neighbor...
                else:
                    pending_neighbor_list.append(nbr)
                    pending_perim_list.append(perim)
                    
        all_neighbors.append([SRC_PSN, pending_neighbor_list, pending_perim_list])
        
print("Number of precincts total: " + str(len(all_neighbors)))

## TEST ##
print("TEST: Boundary precincts properly calculated...")
for i in range(5):
    print(all_boundary_precincts[i])
print("")

## TEST ##
print("TEST: All precincts properly calculated...")
print("Source Precinct, List of Neighbors, List of Shared Perimeters")
print("")
for i in range(5):
    print(all_neighbors[i])
    print("")
	
# Clean the main neighbor list.
for line in all_neighbors:
    src_PSN = line[0]
    nbrs_list = line[1]
    sp_list = line[2]
    
    # If the first and last precincts in the neighbor list are the same...
    if nbrs_list[0] == nbrs_list[-1]:
        # Verify that the shared perimeters are different; i.e. that they are different line segments.
        assert(sp_list[0] != sp_list[-1])
        # Remove last instance of the neighbor.
        nbrs_list.pop()
        # Verify that the neighbors list contains no more duplicates.
        assert(nbrs_list[0] != nbrs_list[-1])
        # Overwrite the old neighbor list with the new one.
        line[1] = nbrs_list
        
        # Add the last shared perimeter length to the first instance of the neighbor's shared boundary.
        to_add = sp_list[-1]
        sp_list.pop()
        sp_list[0] = sp_list[0] + to_add
        line[2] = sp_list
		
## TEST ##
print("TEST: Neighbor lists contain no duplicates...")
print("")
for i in range(5):
    print(all_neighbors[i])
    print("")

###########################################################################
# NOTE: No reason in particular why this code is here. It can probably be moved up but I was in a hurry.
# If the first and last boundary precinct in the list are the same, the border precincts were counted starting in the middle of a precinct...
if all_boundary_precincts[0][0] == all_boundary_precincts[-1][0]:
    # Move the last item in the list to the front.
    all_boundary_precincts.insert(0, all_boundary_precincts.pop())
    # Make sure the first two boundary precincts in the list are the same.
    assert(all_boundary_precincts[0][0]==all_boundary_precincts[1][0])
    # Add the perimeters.
    all_boundary_precincts[1][1] += all_boundary_precincts[0][1]
    # Remove the first boundary precinct, because now it's a duplicate.
    all_boundary_precincts.pop(0)
    
# Print number of boundary precincts.
print("Number of boundary precincts: " + str(len(all_boundary_precincts)))
print("")
###########################################################################
	
# Give each border precinct a location serial number ranging from -n, -(n-1), ..., -1.
# In the dictionary...
# Keys: Border Precinct PSNs.
# Values: Border Precinct Location SNs.

n = len(all_boundary_precincts)

bdr_pct_dict = dict()

for line in all_boundary_precincts:
    pct = line[0]
    perim = line[1]
    # If the boundary precinct is already in the dictionary...
    # i.e. if we have a multi-adjacency
    if pct in bdr_pct_dict:
        # Add the perimeter and corresponding location SN to the border precinct's dictionary.
        bdr_pct_dict[pct][perim] = -n
    else: # if pct not in bdr_pct_dict:
        # Add the precinct as a key to the dictionary.
        # Give it a dictionary with perim as the key and the location SN as the value.
        bdr_pct_dict[pct] = {perim : -n}
        
    n = n-1
	
## TEST ##
print("TEST: Check border precincts...")
for i in range(5):
    print(all_boundary_precincts[i])
	
## TEST ##
print("TEST: Verify border precinct location serial numbers...")
for i in range(5):
    print(str(sorted(bdr_pct_dict.keys())[i]) + ": " + str(bdr_pct_dict[sorted(bdr_pct_dict.keys())[i]]))
print("")
	
## TEST ##
# We can test to check and see if some precincts have multi-adjacency.
# for item in bdr_pct_dict.items():
#     print(item)	
	
print("Total number of border precincts: " + str(len(bdr_pct_dict)))
print("")

## EXAMPLES ##
# print(all_neighbors[0])
# print("")
# print(all_neighbors[4])
# print("")
# print(all_neighbors[2127])
# print("")

#for line in bdr_pct_dict.items():
#    print(line)

# If a precinct in the neighbor list data frame contains -1 in its list of neighbors, it is a border precinct.
# Replace the -1's with the location SNs that give the precinct's relative position around the state border.
for line in all_neighbors:
    src_PSN = line[0]
    nbrs_list = line[1]
    sp_list = line[2]
    if src_PSN in bdr_pct_dict: # if not, do nothing
        # Find the -1 in the neighbor list.
        # Note that this may be either the only instance or merely the first instance.
        idx = nbrs_list.index(-1)
        # Find the corresponding boundary length.
        perim = sp_list[idx]
        # Note the location SN.
        loc_SN = bdr_pct_dict[src_PSN][perim]
        # Replace the instance of -1 in the neighbor list with the location SN.
        nbrs_list[idx] = loc_SN
        
        if loc_SN != -1:
            # Check to see if there is another instance of -1 in the neighbor list.
            if -1 in nbrs_list: # if not, do nothing
                idx2 = nbrs_list.index(-1)
                perim2 = sp_list[idx2]
                loc_SN2 = bdr_pct_dict[src_PSN][perim2]
                # Replace the other instance of -1 with the location SN.
                nbrs_list[idx2] = loc_SN2
                # Make sure this is the last instance of -1.
                assert(-1 not in nbrs_list)
        
# Now we not only know from all_neighbors which precincts border the edge of the state,
# but also their relative positions to the center of the state.
		
## TEST ##
print("TEST: Verify that outerstate neighbors have correct indices...")
print("")
for i in range(5):
    print(all_neighbors[i])
    print("")
	
## TEST ##
print("TEST: Verify that no neighbors appear first AND last in a neighbor list...")

for line in all_neighbors:
    nbrs_list = line[1]
    if nbrs_list[0] == nbrs_list[-1]:
        print(line)
		
# Create a new data frame.
# This one will contain each unique precinct's neighbors and shared perimeters
# as STRINGS of comma-separated values.

neighbor_list = list()

for line in all_neighbors:
    pct = line[0]
    nbrs = line[1]
    perims = line[2]
    
    nbr_list = ''
    sp_list = ''
    
    for nbr in nbrs:
        nbr_list = nbr_list + str(nbr) + ','
    for sp in perims:
        sp_list = sp_list + str(sp) + ','
        
    nbr_list = nbr_list[:-1]
    sp_list = sp_list[:-1]
    
    neighbor_list.append([pct,nbr_list,sp_list])
	
# Write the final data frame to a CSV.

outfile2 = open(sys.argv[3],'w')
writer = csv.writer(outfile2)
writer.writerow(["PSN","nb","sp"])

for row in neighbor_list:
    writer.writerow(row)
    
outfile2.close()

print("Success! Your neighbor list has been generated.")
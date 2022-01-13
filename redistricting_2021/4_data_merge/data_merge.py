#!/usr/bin/python

# TO DO: Edit the output files you want based on which districtings you're looking at (e.g. Congressional, Senate, House).
# To run in the terminal: ./data_merge.py map_input.csv neighbor_input.csv final_output1.csv final_output2.csv final_output3.csv

# Anna Marie Vagnozzi

###################################################################################################

# Import necessary libraries.
import io
import csv
import sys

# Load the map data CSV.
# TO DO: Set the filepath for the CSV file.
infile1 = io.open(sys.argv[1],newline='')
reader1 = csv.reader(infile1)
header1 = next(reader1)

# Read in PSN, unshared, area, pop, voteA, voteB, CongDrft2S, county
# TO DO: Edit column numbers appropriately based on map data.
map_data = [[eval(row[0]),eval(row[1]),eval(row[2]),eval(row[3]),eval(row[4]),eval(row[5]),eval(row[16]),row[21]] for row in reader1]
infile1.close()

# Load the neighbor data CSV.
# TO DO: Set the filepath for the CSV file.
infile2 = io.open(sys.argv[2],newline='')
reader2 = csv.reader(infile2)
header2 = next(reader2)

# Read in PSN, nb, sp.
neighbor_data = [[eval(row[0]),row[1],row[2]] for row in reader2]
infile2.close()

# TO DO: Adjust the assertion statement based on the map.
assert(len(map_data)==2260)
assert(len(neighbor_data)==2260)

InputSCCong = list()
#InputSCSen = list()
#InputSCHouse = list()

for i in range(len(neighbor_data)):
	# Read in a row of data.
	map_line =  map_data[i]
	nbr_line = neighbor_data[i]
	
	# Make sure the precincts match in the map/neighbor data.
	assert(map_line[0] == nbr_line[0])
	
	# To the final data frames, add the columns required for the chain program.
	# Congressional: PSN, nb, sp, unshared, area, pop, voteA, voteB, congD, county
	# Senate: PSN, nb, sp, unshared, area, pop, voteA, voteB, senD, county
	# House: PSN, nb, sp, unshared, area, pop, voteA, voteB, houseD, county
	InputSCCong.append([nbr_line[0],nbr_line[1],nbr_line[2],map_line[1],map_line[2],map_line[3],map_line[4],map_line[5],map_line[6],map_line[7]])
	#InputSCSen.append([nbr_line[0],nbr_line[1],nbr_line[2],map_line[1],map_line[2],map_line[3],map_line[4],map_line[5],map_line[6],map_line[7]])
	#InputSCHouse.append([nbr_line[0],nbr_line[1],nbr_line[2],map_line[1],map_line[2],map_line[3],map_line[4],map_line[5],map_line[6],map_line[7]])

# Write the final datas frame to CSV files.

# TO DO: Specify output files based on the map you are evaluating.

### SC Congressional Map ###

outfile1 = open(sys.argv[3],'wb')
writer = csv.writer(outfile1)
writer.writerow(["","nb","sp","unshared","area","pop","voteA","voteB","congD","county"])

for row in InputSCCong:
    writer.writerow(row)
    
outfile1.close()

### SC Senate Map ###

#outfile2 = open(sys.argv[3],'wb')
#writer = csv.writer(outfile2)
#writer.writerow(["","nb","sp","unshared","area","pop","voteA","voteB","senD","county"])

#for row in InputSCSen:
#    writer.writerow(row)
    
#outfile2.close()

### SC House Map ###

#outfile3 = open(sys.argv[3],'wb')
#writer = csv.writer(outfile3)
#writer.writerow(["","nb","sp","unshared","area","pop","voteA","voteB","houseD","county"])

#for row in InputSCHouse:
#    writer.writerow(row)
    
#outfile3.close()

print("Success! Your Markov Chain input files have been generated.")
# fastq_glue.py
# Made by L. Percival-Alwyn, Twitter: LP_Alwyn
# Joins fastq files together at sequences and quality scores, enabling
# triming from leftside sequence and qual score and takes the header
# from either the left or right side sequence.

# Useful where fastq file reads are derrived from next to one another
# but have been separated due to a custom recipe.
# E.g 14 dark cycles, 136 Read1(b) cycles, rehyb, 14 Read1(a) cycles....
# Read1(a) sequences could be joined to Read1(b) with header from read 1(a) file.

# python fastq_glue.py leftsidefilename.fastq rightsidefilename.fastq <-L/-R> <trim left by interger> joinedfilename.fastq

# E.g python fastq_glue.py leftsidefilename.fastq rightsidefilename.fastq -L 0 joinedfilename.fastq

# import modules
import sys
import os

# Get input files, header side, trim back number and outfile name
try:
    leftside = open(sys.argv[1])
    rightside = open(sys.argv[2])
    use_header_from = sys.argv[3]
    trim_left = int(sys.argv[4])
    outfile = sys.argv[5]

# If no input file print help
except:
    print("\nfastq_glue.py joins fastq files together at sequences and quality scores,\n"
          "enabling triming from leftside sequence and qual score and takes the header\n"
          "from either the left or right side sequence\n"
          "\nUse as follows: \npython fastq_glue.py leftsidefilename.fastq rightsidefilename.fastq <-L/-R> <trim left by interger> joinedfilename.fastq")
    print("\nE.g\npython fastq_glue.py leftsidefilename.fastq rightsidefilename.fastq -L 0 joinedfilename.fastq\n")
    sys.exit()

# Iterate through leftside fastq and add every seq and qual line to the seq_qual list
# Iterate through leftside fastq and add every header to the header list
list_leftside_seq_qual = []
list_leftside_header = []
counter = 0
for x in leftside:
    if (counter-1) % 2 == 0:
        if trim_left < 0:
            list_leftside_seq_qual.append(x.replace("\n","")[:trim_left])
        elif trim_left == 0:
            list_leftside_seq_qual.append(x.replace("\n",""))
        elif trim_left > 0:
            list_leftside_seq_qual.append(x.replace("\n","")[:trim_left * -1])
    if (counter) % 4 == 0:
        list_leftside_header.append(x)
    counter += 1
    
leftside.close()
print(sys.argv[1] + "headers and seq/quals saved to memory")

# Iterate through rightside fastq and add every seq and qual line to the seq_qual list
# Iterate through rightside fastq and add every header to the header list
list_rightside_seq_qual = []
list_rightside_header = []
counter = 0
for x in rightside:
    if (counter-1) % 2 == 0:
        list_rightside_seq_qual.append(x.replace("\n",""))
    if (counter) % 4 == 0:
        list_rightside_header.append(x)
    counter += 1

rightside.close()
print(sys.argv[2] + "headers and seq/quals saved to memory")

# Determine header choice
if use_header_from == "-L":
    header_choice = list_leftside_header

elif use_header_from == "-R":
    header_choice = list_rightside_header
    
# Iterate through the header list (left or right choice) and add header and joined seq/qual plus spacers to new fastq
joined_fastq = ""
counter = 0
        
for x in header_choice:
    try:
        joined_fastq += x
        joined_fastq += list_leftside_seq_qual[counter] + list_rightside_seq_qual[counter] + "\n+\n"
        joined_fastq += list_leftside_seq_qual[counter+1] + list_rightside_seq_qual[counter+1] + "\n"
        counter += 2
    except:
        print("Something went wrong.")

# Write joined fastq to joinedfilename.fastq   
fout = open(outfile, "w")
fout.write(joined_fastq)
fout.close()

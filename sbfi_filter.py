# sbfi_filter.py
# Removes fastqs that do not contain SbfI site.
# 

try:
    from itertools import izip # Python 2 compatible
except:
    izip = zip # Python 3 compatible
    
import sys

try:
    R1_in = open(sys.argv[1])
    R2_in = open(sys.argv[2])
    R1_out = open(sys.argv[3], "w")
    R2_out = open(sys.argv[4], "w")
except:
    print("\nsbfi_filter.py")
    print("Removes fastqs that do not contain SbfI site in the first 6 nucleotides.")
    print("\nUse as follows:")
    print("\npython sbfi_filter.py input_Read1.fastq input_Read2.fastq filtered_output_Read1.fastq filtered_output_Read2.fastq\n\n")
    sys.exit()
    
sbfi = "TGCAGG" # Sequence expected in SbfI cut reads (read 1)

counter = 0
fastqR1 = ""
fastqR2 = ""

# Iterate through each line of the fastq
for R1_line, R2_line in zip(R1_in, R2_in):
    if (counter -1) % 4 == 0:
        if R1_line[:6] == sbfi:
            fastqR1 += R1_line
            fastqR2 += R2_line
            add_fastq = True
        else:
            add_fastq = False
    elif (counter -3) % 4 == 0:
        fastqR1 += R1_line
        fastqR2 += R2_line
        if add_fastq == True:
            R1_out.write(fastqR1)
            R2_out.write(fastqR2)
            fastqR1 = ""
            fastqR2 = ""
        else:
            fastqR1 = ""
            fastqR2 = ""
    else:
        fastqR1 += R1_line
        fastqR2 += R2_line
    counter += 1

R1_in.close()
R2_in.close()
R1_out.close()
R2_out.close()

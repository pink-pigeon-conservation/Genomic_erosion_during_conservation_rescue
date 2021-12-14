# fastq_to_bam_wapper_LSF.py
# Made by L. Percival-Alwyn, Twitter: LP_Alwyn
# Aligns reads to a bwa (-0.7.7) indexed reference using LSF job scheduler.
# It can be used for both single or paired end reads.
# Use as follows:
#       place fastq_to_bam_wapper_LSF.py in the same folder as fastq files
#	python fastq_to_bam_wapper_LSF.py <PE/SE> reference.fa <fastq_list> <options(BSR1XD)> <LSF queue> <optional prefix>
# Where:
# Choosing PE or SE denotes paired or single end sequencing
# fastq_list is a list of fastq files:
#	SE - a simple list (ls *.fastq > fastq_list.txt)
#	PE - list of read 1 and 2 separated by a tab:
#		ls *R1*.fastq > R1fastq_list.txt;
#		ls *R2*.fastq > R2fastq_list.txt;
#		paste R1fastq_list.txt R2fastq_list.txt > fastq_list.txt
# Options are as follows:
#	B or b - make bam
#	BS or bs - make sorted bam
#	BSX or bsx - make sorted indexed bam
#	BSXD or bsxd - make sorted indexed bam and make depth_files (_depth.txt)
#	BSRX or bsrx - make sorted indexed bam with duplicates removed
#	BSRXD or bsrxd - make sorted indexed bam with duplicates removed and make depth_files (_depth.txt)
#	BSR1X or bsr1x - make sorted indexed bam with duplicates removed of only read 1 of PE
#	BSR1XD or bsr1xd - make sorted indexed bam with duplicates removed of only read 1 of PE and make depth_files (_depth.txt)

# Import modules
import sys
import os
import datetime

# Get current time/date and make time and date stamp string with only numbers (means each wrapper run is unique)
now = datetime.datetime.now()
this_moment = str(now)
t_d_stamp = this_moment.replace(" ","").replace("-","").replace(".","").replace(":","")

# Get input files and options
try:
    read_type = sys.argv[1]
    ref = sys.argv[2]
    fastq_files_file = open(sys.argv[3])
    options = sys.argv[4]
    queue = sys.argv[5]
    
    try:
        prefix = sys.argv[6]
    except:
        prefix = ""

# If no input files print help
except:
    print("\nfastq_to_bam_wapper_LSF.py aligns reads to a bwa (-0.7.7) indexed reference using LSF job scheduler.\n")
    print("It can be used for both single or paired end reads.\n")
    print("Use as follows:\n\tplace fastq_to_bam_wapper_LSF.py in the same folder as fastq files")
    print("\tpython fastq_to_bam_wapper_LSF.py <PE/SE> reference.fa <fastq_list> <options(BSR1XD)> <LSF queue> <optional prefix>")
    print("\nWhere:\nChoosing PE or SE denotes paired or single end sequencing")
    print("fastq_list is a list of fastq files:\n\tSE - a simple list (ls *.fastq > fastq_list.txt)")
    print("\tPE - list of read 1 and 2 separated by a tab:\n\t\tls *R1*.fastq > R1fastq_list.txt;")
    print("\t\tls *R2*.fastq > R2fastq_list.txt;")
    print("\t\tpaste R1fastq_list.txt R2fastq_list.txt > fastq_list.txt")
    print("\nOptions are as follows:")
    print("\tB or b - make bam")
    print("\tBS or bs - make sorted bam")
    print("\tBSX or bsx - make sorted indexed bam")
    print("\tBSXD or bsxd - make sorted indexed bam and make depth_files (_depth.txt)")
    print("\tBSRX or bsrx - make sorted indexed bam with duplicates removed")
    print("\tBSRXD or bsrxd - make sorted indexed bam with duplicates removed and make depth_files (_depth.txt)")
    print("\tBSR1X or bsr1x - make sorted indexed bam with duplicates removed of only read 1 of PE")
    print("\tBSR1XD or bsr1xd - make sorted indexed bam with duplicates removed of only read 1 of PE and make depth_files (_depth.txt)\n")
    sys.exit()

# Read fastq file names into list and close file
fastq_files_list = fastq_files_file.readlines()
fastq_files_file.close()

# Make new directory for output file
os.system('mkdir LSF_outputs' + t_d_stamp[:12] + prefix)

# Iterate through the fastq file names list
for lines in fastq_files_list:
    # If paired end reads add the tab separated file names into list make final name from combined names minus.fastq
    if read_type.upper() == "PE" and "B" in options.upper():
        pair = lines.split()
        fastq1  = pair[0].replace(".fastq","")
        fastq2  = pair[1].replace(".fastq","")
        final_name = fastq1 + fastq2
        rmdup_option = ""
        # Align and make PE .sam
        os.system('bsub -o LSF_outputs' + t_d_stamp[:12] + prefix + '/mem' + final_name + '.lsf.txt -R \"rusage[mem=5000] span[ptile=8]\" -J mem' +
                  t_d_stamp + final_name + ' -q ' + queue + ' "source bwa-0.7.7; bwa mem -t 8 ' + ref + ' ' + pair[0] + ' ' + pair[1] +
                  ' > ' + prefix + final_name + '.sam"')
    # If single end reads add file names into list make final name from combined names minus .fastq
    elif read_type.upper() == "SE" and "B" in options.upper():
        single = lines.split()
        fastq1  = single[0].replace(".fastq","")
        final_name = fastq1
        rmdup_option = "-s "
        # Align and make SE .sam
        os.system('bsub -o LSF_outputs' + t_d_stamp[:12] + prefix + '/mem' + final_name + '.lsf.txt -R \"rusage[mem=5000] span[ptile=8]\" -J mem' +
                  t_d_stamp + final_name + ' -q ' + queue + ' "source bwa-0.7.7; bwa mem -t 8 ' + ref + ' ' + single[0] +
                  ' > ' + prefix + final_name + '.sam"')

    # If no PE/SE print out how to get to help and exit
    else:
        print("\nFor help on how to use the wrapper, enter the following command:\n python fastq_to_bam_wapper_LSF.py\n")
        sys.exit()

    # If B in run options make .bam
    if "B" in options.upper():
        os.system('bsub -o LSF_outputs' + t_d_stamp[:12] + prefix + '/makebam' + final_name + '.lsf.txt -w \'ended(mem' + t_d_stamp +
                  final_name + ')\' -R \"rusage[mem=5000]\" -J makebam' + t_d_stamp + final_name + ' -q ' + queue + ' "source samtools-0.1.19; samtools view -q 41 -b -o ' +
                  prefix + final_name + '.bam -S -t ' + ref + ' ' + prefix + final_name + '.sam"')
        wait_for = "makebam" + t_d_stamp

    # If S in run options make .sorted.bam
    if "S" in options.upper():
        os.system('bsub -o LSF_outputs' + t_d_stamp[:12] + prefix + '/sortbam' + final_name + '.lsf.txt -w \'ended(' + wait_for + final_name +
                  ')\' -R \"rusage[mem=5000]\" -J sortbam' + t_d_stamp + final_name + ' -q ' + queue + ' "source samtools-0.1.19; samtools sort ' +
                  prefix + final_name + '.bam ' + prefix + final_name + '.sorted"')
        file_for_indexing = prefix + final_name + '.sorted.bam'
        wait_for = "sortbam" + t_d_stamp

    # If R in run options remove duplicates from .sorted.bam
    if "R" in options.upper():
        os.system('bsub -o LSF_outputs' + t_d_stamp[:12] + prefix + '/rmdupbam' + final_name + '.lsf.txt -w \'ended(' + wait_for + final_name +
                  ')\' -R \"rusage[mem=5000]\" -J rmdupbam' + t_d_stamp + final_name + ' -q ' + queue + ' "source samtools-0.1.19; samtools rmdup ' + rmdup_option +
                  prefix + final_name + '.sorted.bam ' + prefix + final_name + 'rmdup.sorted.bam"')
        file_for_indexing = prefix + final_name + 'rmdup.sorted.bam'
        wait_for = "rmdupbam" + t_d_stamp

    # If 1 in run options make bam of read 1 only
    if "1" in options.upper():
        os.system('bsub -o LSF_outputs' + t_d_stamp[:12] + prefix + '/read1bam' + final_name + '.lsf.txt -w \'ended(' + wait_for + final_name +
                  ')\' -R \"rusage[mem=5000]\" -J read1bam' + t_d_stamp + final_name + ' -q ' + queue + ' "source samtools-0.1.19; samtools view -bh -f 0x0040 ' +
                  prefix + final_name + 'rmdup.sorted.bam > ' + prefix + final_name + 'rmdupread1.sorted.bam"')
        file_for_indexing = prefix + final_name + 'rmdupread1.sorted.bam'
        wait_for = "read1bam" + t_d_stamp

    # If 1 in run options index sorted/rmdup/read1 .bam
    if "X" in options.upper():
        os.system('bsub -o LSF_outputs' + t_d_stamp[:12] + prefix + '/indexbam' + final_name + '.lsf.txt -w \'ended(' + wait_for + final_name +
                  ')\' -R \"rusage[mem=5000]\" -J indexbam' + t_d_stamp + final_name + ' -q ' + queue + ' "source samtools-0.1.19; samtools index ' +
                  file_for_indexing + '"')
        wait_for = "indexbam" + t_d_stamp

    # If D in run options make depth.txt files from indexed sorted/rmdup/read1 .bam
    if "D" in options.upper():
        os.system('bsub -o LSF_outputs' + t_d_stamp[:12] + prefix + '/depthbam' + final_name + '.lsf.txt -w \'ended(' + wait_for + final_name +
                  ')\' -R \"rusage[mem=5000]\" -J depthbam' + t_d_stamp + final_name + ' -q ' + queue + ' "source samtools-0.1.19; samtools depth ' +
                  file_for_indexing + ' > ' + prefix + fastq1 + '_depth.txt"')
        wait_for = "depthbam" + t_d_stamp
        
# When run is complete, search all lsf output files generated by the job, for key words fail and exit (not case sensitive)
# Then email search results and announce that job is complete
os.system('bsub -w \'ended(' + wait_for + '*)\' -q ' + queue + ' "echo All jobs started at ' +
          t_d_stamp + ' with the prefix ' + prefix + ' are finished!; echo The number of "fail"s appearing in output files is;'
          'grep -i fail LSF_outputs' + t_d_stamp[:12] + prefix + '/* | wc -l;'
          'echo The number of "exit"s appearing in output files is;grep -i exit LSF_outputs' + t_d_stamp[:12] + prefix + '/* | wc -l"')

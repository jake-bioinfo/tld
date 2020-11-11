#!/usr/bin/env python

import os, os.path, sys, getopt
from Bio import SeqIO

# Defining variables
# Usage function
def usage():
  print ("\nThis is the usage function:")
  print ('Usage: '+sys.argv[0]+' -h <help> -i <inputFile> -o <outputFile> -w <windowSize> -s <stepSize>\n')

# Parsing arguments
def main(argv):
  in_file = ''
  out_file = ''
  win_sz = ''
  step_sz = ''

  try:
    opts, args = getopt.getopt(argv,"hi:o:w:s:",["help","inputFile=","outputFile=","windowSize=","stepSize="])
    if not opts:
      print ('\nError:\nNo options supplied')
      usage()
      sys.exit()
  except getopt.GetoptError as err:
    print ("\nError:\n", str(err))
    usage()
    sys.exit(2)
  for opt, arg in opts:
    if opt in ("-h", "--help"):
      usage()
      sys.exit()
    elif opt in ("-i", "--inputFile"):
      in_file = arg
    elif opt in ("-o", "--outputFile"):
      out_file = arg
    elif opt in ("-w", "--windowSize"):
      win_sz = int(arg)
    elif opt in ("-s", "--stepSize"):
      step_sz = int(arg)

  print ('Variables list\n'  )
  print ('Input file is :', in_file)
  print ('Output file is :', out_file)
  print ('Window size is :', win_sz)
  print ('Step size is :', step_sz, '\n')

# Sliding window work
  with open(out_file,"w") as f:
    for seq_record in SeqIO.parse(in_file, "fasta"):
      i = 0
      l = len(seq_record.seq)
      while (i < l):
        seq_start = i

        if (i + int(win_sz)) < l:
          seq_end = i + int(win_sz)

        else:
          seq_end = i + (l - i)

        seq_start_nm = str(seq_start)
        seq_end_nm = str(seq_end)

        if seq_start_nm == "0":
          seq_start_nm = str("00000")
        if (len(seq_start_nm)) <= 1:
          seq_start_nm = str("0000" + seq_start_nm)
        if (len(seq_start_nm)) <= 2:
          seq_start_nm = str("000" + seq_start_nm)
        if (len(seq_start_nm)) <= 3:
          seq_start_nm = str("00" + seq_start_nm)
        if (len(seq_start_nm)) <= 4:
          seq_start_nm = str("0" + seq_start_nm)
        if seq_end_nm == "0":
          seq_end_nm = str("00000")
        if (len(seq_end_nm)) <= 1:
          seq_end_nm = str("0000" + seq_end_nm)
        if (len(seq_end_nm)) <=2:
          seq_end_nm = str("000" + seq_end_nm)
        if (len(seq_end_nm)) <= 3:
          seq_end_nm = str("00" + seq_end_nm)
        if (len(seq_end_nm)) <= 4:
          seq_end_nm = str("0" + seq_end_nm)

        base=os.path.basename(in_file)
        file_str=base.split('.')[0]        
 
        f.write(">" + file_str + "__" + seq_record.id + "/seqRng:" + seq_start_nm + "-" + seq_end_nm + "\n")
        f.write(str(seq_record.seq[seq_start:seq_end]) + "\n")
        
        if seq_end == l:
          i = l	
        
        else:
          i += int(step_sz)

if __name__ == "__main__":
  main(sys.argv[1:])

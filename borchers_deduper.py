#Sam file to be deduped MUST be sorted by chromosome
#import packages to be used later
import re
import argparse
import gzip

#Implement argparse so that users can specify the file to be deduped and the list of umis.
def get_args():
    parser = argparse.ArgumentParser(description="DemultiplexingProgram")
    parser.add_argument("-f", "--file", help="Type -f followed by the path to the SAM file to be deduped.", required=True)
    parser.add_argument("-u", "--umi", help="Type -u followed by the path to a file containing the UMIs. Each UMI should be on a new line.", required=False)
    parser.add_argument("-p", "--paired", help="paired-end reads are not supported in this program.")
    return parser.parse_args()

args = get_args()
sam_file = args.file
umi_file = args.umi
paired_end = args.paired

#Tell user that paired end functionality is not supported
if paired_end:
    print("Paired end reads are not supported.")
    quit()

#take provided filename and add deduped to it, this will be ouput file title
fname = sam_file
new_filename = fname + "_" + "deduped"

#make an empty dictionary for each umi and a lower digit code to be assigned later
UMI_dic = {}

#This is a function to parse a CIGAR string, provided as a variable called CIGARstring, and return the 5' mapping position for a reverse strand
def get_reverse_pos(CIGARstring):
    CIGARlist = re.findall(r'[0-9]+[A-Z]', CIGARstring)
    length_spanned = 0
    for index in range(len(CIGARlist)):
        if "M" in CIGARlist[index]:
            variable1 = str(CIGARlist[index])
            variable2 = variable1.strip("M")
            length_spanned += int(variable2)
        elif "N" in CIGARlist[index]:
            variable1 = str(CIGARlist[index])
            variable2 = variable1.strip("N")
            length_spanned += int(variable2)
        elif "D" in CIGARlist[index]:
            variable1 = str(CIGARlist[index])
            variable2 = variable1.strip("D")
            length_spanned += int(variable2)
        elif "S" in CIGARlist[index]:
            if index >> 1:
                variable1 = str(CIGARlist[index])
                variable2 = variable1.strip("S")
                length_spanned += int(variable2)
    mapping_pos = length_spanned + int(align_pos)
    return(mapping_pos)

#open the provided list of umis and assign each one a numerical identifier, to lighten memory requirements for the program
with open (umi_file, "r") as fh:
    count = 0
    for line in fh:
        count += 1
        umi = line.strip()
        UMI_dic[umi] = count

#Begin looping through sam file to be deduped
with open (sam_file, "r") as f1:
    new_sam = open(new_filename, "a+")
    align_dic = {}
    chr = 0
    #set each regex variable to None so that if any regex matches nothing, the variable still exists and the loop can continue on the next line
    for line in f1:
        a=None
        b=None
        c=None
        d=None
        e=None
        info = line
        #ignore the header lines
        header = re.search(r'^[@][\S]+', info)
        if header:
            continue
        #extract the chromosome
        a = re.search(r'^[\S]+\t[\S]+\t([0-9]+)', info)
        if a is None:
            continue
        chrm = a.group(1)
        if chrm == chr:
            pass
        #if the chromosome doesn't match the last seen chromosome, empty the dictionary of unduplicated read info. This improves memory demand.
        else:
            chr = chrm
            align_dic = {}
        #extract UMI in sam file
        b = re.search(r'^[\S]+[:]([A-Z]+)\t', info)
        if b is None:
            continue
        umi = b.group(1)
        #check if observed UMI is in the provided list of umis
        if umi in UMI_dic:
            #extract alignment position according to SAM file
            c = re.search(r'^[\S]+\t[\S]+\t[\S]+\t([0-9]+)\t', info)
            if c is None:
                continue
            align_pos = c.group(1)
            #extract flag from SAM file
            d = re.search(r'^[\S]+\t([0-9]+)\t', info)
            if d is None:
                continue
            flag = d.group(1)
            #check if the flag indicates forward or reverse strand
            if int(flag) == 16:
                strnd = "rev"
                #extract CIGAR string
                e = re.search(r'^[\S]+\t[\S]+\t[\S]+\t[\S]+\t[\S]+\t([0-9,A-Z]+)', info)
                if e is None:
                    continue
                CIGARstring = e.group(1)
                #parse CIGAR string using get_reverse_pos function to determine the adjusted mapping position
                get_reverse_pos(CIGARstring)
                store_umi = UMI_dic[umi]
                #If the umi, adjusted mapping position, and strand value observed in the sam file line are already in the dictionary,the line(read) is a PCR duplicate. If a duplicated, move to the next line and keep looping. If not, store in dictionary
                if store_umi in align_dic:
                    if align_dic[store_umi][0] == mapping_pos:
                        if align_dic[store_umi][1] == strnd:
                            continue
                else:
                    align_dic[store_umi] = mapping_pos, strnd
                new_sam.write(info)
            else:
                strnd = "fwd"
                #extract CIGAR string
                e = re.search(r'^[\S]+\t[\S]+\t[\S]+\t[\S]+\t[\S]+\t([0-9,A-Z]+)', info)
                CIGARstring = e.group(1)
                CIGARlist = re.findall(r'[0-9]+[A-Z]', CIGARstring)
                #Parse CIGAR string for front clipping only, return adusted mapping position
                if "S" in CIGARlist[0]:
                    clipping_string = str(CIGARlist[0])
                    clipping_int = clipping_string.strip("S")
                    mapping_pos = int(align_pos) - int(clipping_int)
                    align_pos = mapping_pos
                store_umi = UMI_dic[umi]
                #If the umi, adjusted mapping position, and strand value observed in the sam file line are already in the dictionary,the line(read) is a PCR duplicate. If a duplicated, move to the next line and keep looping. If not, store in dictionary
                if store_umi in align_dic:
                    if align_dic[store_umi][0] == align_pos:
                        if align_dic[store_umi][1] == strnd:
                            continue
                else:
                    align_dic[store_umi] = align_pos, strnd
                new_sam.write(info)
        else:
            continue

# Codon-Counting
#This python project identifies the number of sequences in a file with a corresponding name. It also generates a list of A, T, G, C codons (sets of 3) so that the frequency of each amino acid in each sequence can be tallied.
#				ASSIGNMENT WEEK 5
#
#19232985

#Priyanka Kooverjee

#Biochemistry 324

#Assignment 5

#Submitted: 18/03/2018

#Due: 31/03/2018

class Read_FastA:

    #==================================================================
    # This function return a tuple of 2 lists 'sequence_names' and
    #'sequences' containing the sequence names and the sequences of the
    # entries in the supplied fastA file.  There should be no spaces,
    # newlines or '>' characters in either 'sequence_names' or
    # 'sequences'. Sequence_names must be a list composed of string
    # items corresponding to each sequence name sequences must be a
    # list of string items, each item corresponding to one sequence
    #
    #==================================================================

    def Read_FastA_Names_And_Sequences(self, filepath):
        print("Reading fastA sequences...")
        self.sequence_names = []
        self.sequences = []
        f = open(filepath, 'r')
        self.counter = 0
        for i in f:
            if (i[0] == '>'):
                self.counter += 1
                self.sequence_names.append(i[1:].replace('\n', ''))
                self.sequences.append(str())
            else:
                self.sequences[self.counter - 1] = self.sequences[self.counter - 1] + i.replace('\n', '')
        f.close()
        return (self.sequence_names, self.sequences)


class Read_GFF:

    def Get_Gene_Positions(self, list_of_gff_features, filepath, feature):
        filepath = open(filepath, 'r')
        data = filepath.read().splitlines()
        hold = []
        count = 0
        for i in data:
            hold.append(i.split())
        for i in hold:
            if i.__contains__(feature):
                tupl = i[0], i[3], i[4], i[7]
                list_of_gff_features.append(tupl)

        return (list_of_gff_features)

    #==================================================================
    # This function should be passed a list 'list_of_gff_features' to
    # which you append a tuple (seqID, start, end, offset) that
    # contains the information from each line from the GFF file
    # corresponding to a coding sequences (CDS).  Thus, the information
    # from the line 
    # 'chrI	SGD	CDS	335	649	.	+	0
    #	Parent=YAL069W_mRNA;Name=YAL069W_CDS;orf_classificat...'
    # must be appended to the list_of_gff_features as
    # ('chrI','335','649','0').  Step through each line of the GFF 
    # file, selecting only the ones with column 2 == 'CDS', and append
    # the seqID, start, end, offset information of each as a tuple to
    # 'list_of_gff_features'. Filepath is the full path to the
    # 'saccharomyces_cerevisiae_2018.gff' file and feature == 'CDS'
    #==================================================================



class My_Codons:

    #==================================================================
    # The function generates a list of all possible combinations of 3
    #  of the 4 nucleotides G, A, T and C, and returns the list as
    # 'codons'
    #==================================================================

    def Make_List_Of_Codons(self):
        self.nucleotides = ['G','A','T','C']
        self.codons = []
        self.tempcodon = ''
        for a in self.nucleotides:
            for b in self.nucleotides:
                for c in self.nucleotides:
                    self.tempcodon = a+b+c
                    self.codons.append(self.tempcodon)
        return(self.codons)



    def Count_Codons(self, sequence, codons, number_of_occurrences, offset=0):
        seq = []
        count = 0
        for i in range(0, len(sequence)-2, 3):
            seq.append(sequence[i])
            seq[count] += sequence[i + 1]
            seq[count] += sequence[i +2]
            count = count +1
        counter = 0
        for i in codons:
            for j in seq:
                if i == j:
                    number_of_occurrences[counter]= number_of_occurrences[counter]+1
            counter = counter +1

        return (number_of_occurrences)

                #==================================================================
    # The string.count(substring,start,end) method may look like a
    # suitable function to use, but it finds the occurrence of a
    # substring in ANY FRAME.  Codons are arranged in non-overlapping
    # groups of three.  So we have to begin searching groups of three,
    # starting at the beginning of the sequence, taking care of any
    # offset defined in the gff file, and then jumping by three bases
    # after each comparison.  We write our own function to do exactly
    # this.
    #==================================================================



path_of_gff_file = 'saccharomyces_cerevisiae_2018.gff' # change the path string if yours is different
path_of_fasta_file = 'saccharomyces_cerevisiae_2018.fna' # change the path string if yours is different

# Get the positions and offset of all codings sequences (CDS) in the yeast genome

GFF_file_object = Read_GFF()
list_of_gff_features=[]
total_sequence_length = 0
list_of_gff_features = GFF_file_object.Get_Gene_Positions(list_of_gff_features, path_of_gff_file,'CDS')

# make a list of all 64 possible codons

codon_object = My_Codons()
codons = codon_object.Make_List_Of_Codons()

# Read the chromosome sequences

FASTA_file_object = Read_FastA()
sequence_name, sequences = FASTA_file_object.Read_FastA_Names_And_Sequences(path_of_fasta_file)

# Loop over list_of_gff_features, using one entry at a time

number_of_occurrences =[0]*64
print('Counting codons...')
for gff_line in list_of_gff_features:

# get chromosome and slice the gene sequence of the chromosome with the calculated index
    chromosome_sequence = sequences[sequence_name.index(gff_line[0])]
    number_of_occurrences = codon_object.Count_Codons(chromosome_sequence[int(gff_line[1])-1:int(gff_line[2])], codons, number_of_occurrences, int(gff_line[3]))

# Print out the total codons
total_codons = sum(number_of_occurrences)
for i in range(0,64):
    if(i == 0):
        print('Codon','Number','Frequency (/1000)')
    print(codons[i],number_of_occurrences[i],1000*number_of_occurrences[i]/total_codons)

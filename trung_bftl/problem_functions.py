"""
Rosalind Problem 1-40, write as a lib function
"""
# Problem 1 - Counting DNA Nucleotides
def count_dna_base(dna):
    """
    Counts the number of each nucleotide in a DNA sequence.
    Args:
        dna (str): The input DNA sequence.
    Returns:
        dict: A dictionary with the counts of each base type.
    """
    dna = dna.upper()
    count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for base in dna:
        count[base] += 1
    # Convert dict in to a series of number divided by space
    count_freq = ' '.join([str(count[base]) for base in count])
    return count, count_freq

# Problem 2 - Transcribing DNA into RNA
def transcribe_dna_to_rna(dna):
    """
    Converts a DNA sequence to an RNA sequence.
    Args:
        dna (str): The input DNA sequence.
    Returns
        str: The transcribed RNA sequence.
    """
    return dna.replace('T', 'U')

# Problem 3 - Complementing a Strand of DNA
def reverse_complement(dna):
    """
    Finds the reverse complement of a DNA sequence.
    Args:
        dna (str): The input DNA sequence.
    Returns:
        str: The reverse complement of the DNA sequence.
    """
    dna = dna.upper()
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in reversed(dna)])

# Problem 4 - Rabbits and Recurrence Relations
def fibonacci(n, k):
    """
    Computes the nth Fibonacci number given the initial number of pairs and the number of offspring produced.
    Args:
        n (int): The number of months to simulate.
        k (int): The number of offspring produced by each pair.
    Returns:
        int: The number of rabbit pairs after n months.
    """
    if n == 1:
        return 1
    if n == 2:
        return 1
    return fibonacci(n - 1, k) + k * fibonacci(n - 2, k)

# Problem 5 - Counting Point Mutations
def hamming_distance(dna1, dna2):
    """
    Computes the Hamming distance between two DNA sequences.
    Args:
        dna1 (str): The first DNA sequence.
        dna2 (str): The second DNA sequence.
    Returns:
        int: The Hamming distance between the two DNA sequences.
    """
    dna1 = dna1.upper()
    dna2 = dna2.upper()
    return sum([1 for i in range(len(dna1)) if dna1[i] != dna2[i]])

# Problem 6 - Mendel's First Law
def mendel_first_law(k, m, n):
    """
    Computes the probability that two randomly selected mating organisms will produce an individual possessing a dominant allele.
    Args:
        k (int): The number of homozygous dominant individuals.
        m (int): The number of heterozygous individuals.
        n (int): The number of homozygous recessive individuals.
    Returns:
        float: The probability of producing an individual with a dominant allele.
    """
    total = k + m + n

    # Probabilities of each combination
    prob_kk = (k / total) * ((k - 1) / (total - 1))
    prob_km = (k / total) * (m / (total - 1)) + (m / total) * (k / (total - 1))
    prob_kn = (k / total) * (n / (total - 1)) + (n / total) * (k / (total - 1))
    prob_mm = (m / total) * ((m - 1) / (total - 1)) * 0.75
    prob_mn = (m / total) * (n / (total - 1)) * 0.5 + (n / total) * (m / (total - 1)) * 0.5

    # Total probability of dominant allele
    prob_dominant = prob_kk + prob_km + prob_kn + prob_mm + prob_mn

    return prob_dominant.__round__(5)

# Problem 7 - Translating RNA into Protein
def rnaTranslate(string):
    """
    Translates an RNA sequence to a protein sequence.
    Args:
        string (str): The input RNA sequence.
    Returns:
        str: The translated protein sequence.
    """
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    protein_seq = ''
    for i in range(0, len(string), 3):
        codon = string[i:i+3]
        protein_seq += codon_table[codon]
        if protein_seq[-1] == '*':
            break
        # remove the * from string
    protein_seq = protein_seq.replace('*', '')
    return protein_seq

# Problem 8 - Finding a Motif in DNA
def find_motif(dna, motif):
    """
    Finds all occurrences of a motif in a DNA sequence.
    Args:
        dna (str): The input DNA sequence.
        motif (str): The motif to search for.
    Returns:
        list: The starting positions of the motif in the DNA sequence.
    """
    positions = []
    for i in range(len(dna) - len(motif) + 1):
        if dna[i:i+len(motif)] == motif:
            positions.append(i + 1)
    return positions.__str__().replace('[', '').replace(']', '').replace(',', '')

# Problem 9 - Consensus and Profile
def consensus_profile(dnas):
    """
    Computes the consensus sequence and profile matrix of a list of DNA sequences.
    Args:
        dnas (list): A list of DNA sequences.
    Returns:
        tuple: A tuple containing the consensus sequence and the profile matrix.
    """
    profile = {'A': [], 'C': [], 'G': [], 'T': []}
    # loop through all the dna sequences in the list
    for i in range(len(dnas[0])):
        count, _ = count_dna_base(''.join([dna[i] for dna in dnas]))
        for base in count:
            if base in profile:
                profile[base].append(count[base]) # append the count of each base to the profile
    consensus = ''
    for i in range(len(dnas[0])):
        max_count = 0
        max_base = ''
        for base in profile:
            if profile[base][i] > max_count:
                max_count = profile[base][i]
                max_base = base
        consensus += max_base
        # We need to format the dict output to paste to Rosalind
    profile_output = ''
    for base in profile:
        profile_output += f'{base}: {" ".join([str(count) for count in profile[base]])}\n'
    return consensus, profile, profile_output

# Problem 10 - Mortal Fibonacci Rabbits
def mortal_fibonacci(n, m):
    """
    Computes the number of rabbit pairs after n months given that each pair of rabbits dies after m months.
    Args:
        n (int): The number of months to simulate.
        m (int): The number of months each pair of rabbits lives.
    Returns:
        int: The number of rabbit pairs after n months.
    """
    rabbits = [1, 1]
    for i in range(2,n): # start from 2 because we already know the first 2 months
        if i < m: # if i is less than m, we add the last 2 months
            rabbits.append(rabbits[-1] + rabbits[-2])
        elif i == m or i == m+1: # if i is equal to m or m+1, we add the last 2 months and subtract 1
            rabbits.append(rabbits[-1] + rabbits[-2] - 1)
        else: # if i is greater than m, we add the last 2 months and subtract the m-th month
            rabbits.append(rabbits[-1] + rabbits[-2] - rabbits[-(m+1)])
    return rabbits[-1].__int__()

# Problem 11 - Computing GC Content
def gc_content(dna):
    """
    Computes the GC content of a DNA sequence.
    Args:
        dna (str): The input DNA sequence.
    Returns:
        float: The maximum percentage of GC content in the DNA sequence.
    """
    dna = dna.upper()
    count = {'G': 0, 'C': 0}
    for base in dna:
        if base in count:
            count[base] += 1
    percent_gc = (count['G'] + count['C']) / len(dna) * 100
    # Rosalind allows for a default error of 0.001 in all decimal answers unless otherwise stated
    formatted_percent = "{:.6f}".format(percent_gc)
    return formatted_percent

# Problem 12 - Calculating Expected Offspring
def expected_offspring(couples):
    """
    Computes the expected number of offspring displaying the dominant phenotype.
    Args:
        couples (list): A list sep by space of the number of couples with each genotype AA-AA, AA-Aa, AA-aa, Aa-Aa, Aa-aa, aa-aa.
    Returns:
        int: The expected number of offspring displaying the dominant phenotype.
    """
    expected_offspring = 0
    for i in range(len(couples)):
        expected_offspring += couples[i] * [1, 1, 1, 0.75, 0.5, 0][i]
    # we multiply by 2 because each couple is assumed to have 2 offspring
    final_offspring = expected_offspring * 2
    return final_offspring
# Problem 13 - Finding a Shared Motif
def find_shared_motif(dnas):
    """
    Finds the longest common motif shared by a list of DNA sequences.
    Args:
        dnas (list): A list of DNA sequences.
    Returns:
        str: The longest common motif shared by the DNA sequences.
    """
    dnas.sort(key=len)
    for i in range(len(dnas[0]), 0, -1):
        for j in range(len(dnas[0]) - i + 1):
            motif = dnas[0][j:j+i]
            if all([motif in dna for dna in dnas]):
                return motif

# Problem 14 - Independent Alleles
def independent_alleles(k, n):
    """
    Computes the probability that at least n AaBb organisms will belong to the k-th generation.
    Args:
        k (int): The generation number.
        n (int): The minimum number of AaBb organisms.
    Returns:
        float: The probability that at least n AaBb organisms will belong to the k-th generation.
    """
    from math import comb
    total = 2 ** k
    prob = 0
    for i in range(n, total + 1):
        prob += comb(total, i) * 0.25 ** i * 0.75 ** (total - i)
    return prob.__round__(3)

# Problem 15 - Finding a Protein Motif
def find_protein_motif(protein):
    """
    Finds the locations of the N-glycosylation motif in a protein sequence.
    Args:
        protein (str): The input protein sequence.
    Returns:
        list: The locations of the N-glycosylation motif in the protein sequence.
    """
    import re
    motif = re.compile(r'(?=(N[^P][ST][^P]))')
    return [m.start() + 1 for m in motif.finditer(protein)]
# Problem 16 - Inferring mRNA from Protein
def mrna_from_protein(protein):
    """
    Computes the number of RNA strings from which the protein could have been translated.
    Args:
        protein (str): The input protein sequence.
    Returns:
        int: The number of RNA strings that could have been translated to the protein.
    """
    codon_table = {
        'F': 2, 'L': 6, 'S': 6, 'Y': 2, '*': 3,
        'C': 2, 'W': 1, 'P': 4, 'H': 2, 'Q': 2,
        'R': 6, 'I': 3, 'M': 1, 'T': 4, 'N': 2,
        'K': 2, 'V': 4, 'A': 4, 'D': 2, 'E': 2,
        'G': 4
    }
    return codon_table['*'] * 3 % 1000000

# Problem 17 - ORF Open Reading Frames
def find_orfs(dna):
    """
    Finds all open reading frames in a DNA sequence.
    Args:
        dna (str): The input DNA sequence.
    Returns:
        list: The open reading frames in the DNA sequence
    """
    dna = dna.upper()
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    for i in range(len(dna) - 2):
        if dna[i:i+3] == start_codon:
            for j in range(i + 3, len(dna), 3):
                if dna[j:j+3] in stop_codons:
                    orfs.append(dna[i:j])
                    break
    return orfs

# Problem 18 - Enumerating Gene Orders
def enumerate_gene_orders(n):
    """
    Enumerates all possible permutations of a list of numbers.
    Args:
        n (int): The number of elements in the list.
    Returns:
        list: A list of all possible permutations of the list.
    """
    from itertools import permutations
    return list(permutations(range(1, n + 1)))

# Problem 19 - Calculating Protein Mass
def protein_mass(protein):
    """
    Computes the mass of a protein sequence.
    Args:
        protein (str): The input protein sequence.
    Returns:
        float: The mass of the protein sequence.
    """
    mass_table = {
        'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
        'F': 147.06841, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
        'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
        'P': 97.05276, 'Q': 128.05858, 'R': 156.10111, 'S': 87.03203,
        'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333
    }
    return sum([mass_table[aa] for aa in protein])

# Problem 20 - Locating Restriction Sites
def find_restriction_sites(dna):
    """
    Finds the locations of reverse palindromic restriction sites in a DNA sequence.
    Args:
        dna (str): The input DNA sequence.
    Returns:
        list: The locations of the reverse palindromic restriction sites in the DNA sequence.
    """
    dna = dna.upper()
    revcomp = reverse_complement(dna)
    sites = []
    for i in range(len(dna)):
        for j in range(4, 13):
            if i + j > len(dna):
                break
            if dna[i:i+j] == revcomp[-(i+j):-(i)]:
                sites.append((i + 1, j))
    return sites

# Problem 21 - RNA Splicing
def rna_splicing(dna, introns):
    """
    Splices the introns out of a DNA sequence and transcribes the exons to an RNA sequence.
    Args:
        dna (str): The input DNA sequence.
        introns (list): A list of intron sequences to be spliced out.
    Returns:
        str: The transcribed RNA sequence after splicing out the introns.
    """
    for intron in introns:
        dna = dna.replace(intron, '')
    return transcribe_dna_to_rna(dna)
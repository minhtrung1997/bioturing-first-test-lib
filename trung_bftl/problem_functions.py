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
    total_rna = 1
    # Attach the * to the protein sequence
    protein += '*'
    for aa in protein:
        total_rna *= codon_table[aa]
    return total_rna % 1000000

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
    dna_rc = reverse_complement(dna)
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    for i in range(len(dna) - 2):
        # Find orf on the sense strand
        if dna[i:i+3] == start_codon:
            for j in range(i + 3, len(dna), 3):
                if dna[j:j+3] in stop_codons:
                    orfs.append(dna[i:j])
                    break
        # Find orf on the antisense strand
        if dna_rc[i:i+3] == start_codon:
            for j in range(i + 3, len(dna_rc), 3):
                if dna_rc[j:j+3] in stop_codons:
                    orfs.append(dna_rc[i:j])
    # Deduplicate the list of orfs
    orfs = list(set(orfs))
    # Translate the DNA sequence to protein
    proteins = []
    for orf in orfs:
        proteins.append(rnaTranslate(transcribe_dna_to_rna(orf)))
    # Deduplicate the list of proteins
    proteins = list(set(proteins))
    return orfs, proteins

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
    mass = sum([mass_table[aa] for aa in protein]).__round__(3)
    return mass

# Problem 20 - Locating Restriction Sites
def find_restriction_sites(dna):
    """
    Finds the locations of reverse palindromic restriction sites in a DNA sequence.
    Args:
        dna (str): The input DNA sequence.
    Returns:
        list: The locations of the reverse palindromic restriction sites in the DNA sequence,
        having length between 4 and 12.
    """
    dna = dna.upper()
    revcomp = reverse_complement(dna)
    sites = []
    # Check specifically for position 1
    for j in range(4, 13):
        if dna[0:j] == dna[-1:-(j+1):-1]:
            sites.append((1, j))
    # Check for all other positions
    for i in range(0, len(dna), 1):
        for j in range(4, 13):
            if i + j > len(dna):
                break
            if dna[i:i+j] == revcomp[-(i+j):-(i)]:
                if dna[0:j] == revcomp[-(i+j):-(i)]: # Check if the site is a reverse palindrome of the start sequence position
                    sites.append((1, j))
                sites.append((i + 1, len(dna[i:i+j])))
    # deduplicate the list of sites
    sites = list(set(sites))
    # Sort according to the starting position
    sites.sort(key=lambda x: x[0])
    formatted_sites = '\n'.join([f'{site[0]} {site[1]}' for site in sites])

    return sites, formatted_sites   
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

# Problem 22 - Introduction to Random Strings
def random_string_probability(dna, gc_content):
    """
    Computes the probability of a random string having the given GC content.
    Args:
        dna (str): The input DNA sequence.
        gc_content (float): The GC content of the random string.
    Returns:
        float: The probability of a random string having the given GC content.
    """
    import math
    gc_prob = gc_content / 2
    at_prob = (1 - gc_content) / 2
    prob = 1
    for base in dna:
        if base in 'GC':
            prob *= gc_prob
        else:
            prob *= at_prob
    # Convert the probability to logaritmic form
    log_prob = math.log10(prob)
    # format to 3 decimal places and trailing zeros
    log_prob = "{:.3f}".format(log_prob)
    return log_prob

# Problem 23 - Overlap Graphs
def is_overlap_graphs(dna1, dna2, n):
    """
    Checks if two DNA sequences overlap with each other.
    Args:
        dna1 (str): The first DNA sequence.
        dna2 (str): The second DNA sequence.
        n (int): The length of the overlap.
    Returns:
        bool: True if the DNA sequences overlap, False otherwise.
    """
    return dna1[-n:] == dna2[:n]

# Problem 24 - Enumerating k-mers Lexicographically
def enumerate_kmers_lexicographically(alphabet, n):
    """
    Enumerates all possible k-mers of a given alphabet and length lexicographically.
    Args:
        alphabet (str): The alphabet to use for the k-mers.
        n (int): The length of the k-mers.
    Returns:
        list: A list of all possible k-mers of the given alphabet and length.
    """
    from itertools import product
    kmers = [''.join(kmer) for kmer in product(alphabet, repeat=n)]
    return kmers

# Problem 25 - Longest Increasing Subsequence
def longest_increasing_subsequence(sequence):
    """
    Computes the longest increasing subsequence of a sequence of numbers.
    Args:
        sequence (list): A list of numbers.
    Returns:
        list: The longest increasing subsequence of the sequence.
    """
    n = len(sequence)
    lis = [1] * n
    for i in range(1, n):
        for j in range(0, i):
            if sequence[i] > sequence[j]:
                lis[i] = max(lis[i], lis[j] + 1)
    length = max(lis)
    subsequence = []
    for i in range(n - 1, -1, -1):
        if lis[i] == length:
            subsequence.append(sequence[i])
            length -= 1
    subsequence.reverse()
    return subsequence
def longest_decreasing_subsequence(sequence):
    """
    Computes the longest decreasing subsequence of a sequence of numbers.
    Args:
        sequence (list): A list of numbers.
    Returns:
        list: The longest decreasing subsequence of the sequence.
    """
    n = len(sequence)
    lds = [1] * n
    for i in range(1, n):
        for j in range(0, i):
            if sequence[i] < sequence[j]:
                lds[i] = max(lds[i], lds[j] + 1)
    length = max(lds)
    subsequence = []
    for i in range(n - 1, -1, -1):
        if lds[i] == length:
            subsequence.append(sequence[i])
            length -= 1
    subsequence.reverse()
    return subsequence

# Problem 26 - Genome Assembly as Shortest Superstring
def shortest_superstring(dnas):
    """
    Computes the shortest superstring that contains all the given DNA sequences.
    Args:
        dnas (list): A list of DNA sequences.
    Returns:
        str: The shortest superstring that contains all the DNA sequences.
    Note: Only gluing pairs of strings that overlap by at least half of their length.
    """
    def overlap(a, b):
        """
        Computes the overlap between two strings.
        Args:
            a (str): The first string.
            b (str): The second string.
        Returns:
            int: The overlap between the two strings.
        """
        max_overlap = 0
        for i in range(1, len(a)):
            if b.startswith(a[i:]):
                max_overlap = len(a) - i
                break
        return max_overlap

    while len(dnas) > 1:
        max_overlap = -1
        best_pair = (0, 0)
        best_string = ""
        
        for i in range(len(dnas)):
            for j in range(len(dnas)):
                if i != j:
                    ov = overlap(dnas[i], dnas[j])
                    if ov > max_overlap:
                        max_overlap = ov
                        best_pair = (i, j)
                        best_string = dnas[i] + dnas[j][ov:]
        
        i, j = best_pair
        # Remove the string with higher index first
        # Remove the string with the higher index first
        if i > j:
            dnas.pop(i)
            dnas[j] = best_string
        else:
            dnas.pop(j)
            dnas[i] = best_string
        print(f'Merging strings {i if i > j else j} with {j if i > j else i} with overlap {max_overlap}')

    return dnas[0]

# Problem 27 - Perfect Matchings and RNA Secondary Structures
def perfect_matchings(rna):
    """
    Computes the number of perfect matchings in an RNA sequence.
    Args:
        rna (str): The input RNA sequence.
    Returns:
        int: The number of perfect matchings in the RNA sequence.
    """
    from math import factorial
    au = rna.count('A')
    gc = rna.count('G')
    return factorial(au) * factorial(gc)

# Problem 28 - Partial Permutations
def partial_permutations(n, k):
    """
    Computes the number of partial permutations of k elements from a set of n elements.
    Args:
        n (int): The total number of elements in the set.
        k (int): The number of elements to permute.
    Returns:
        int: The number of partial permutations of k elements from a set of n elements.
    """
    from math import factorial
    return (factorial(n) // factorial
            (n - k)) % 1000000

# Problem 29 - Enumerating Oriented Gene Orderings
def enumerate_oriented_gene_orderings(n):
    """
    Enumerates all possible oriented gene orderings of a list of numbers.
    Args:
        n (int): The number of elements in the list.
    Returns:
        list: A list of all possible oriented gene orderings of the list.
    """
    from itertools import permutations
    from itertools import product
    orientations = ['+', '-']
    permutations = list(permutations(range(1, n + 1)))
    oriented_permutations = []
    for perm in permutations:
        for orientation in product(orientations, repeat=n):
            # Append as a positive or negative number based on the orientation
            oriented_permutations.append(' '.join([f'{o}{p}' for o, p in zip(orientation , perm)]))
    return oriented_permutations

# Problem 30 - Finding a Spliced Motif
def find_spliced_motif(dna, motif):
    """
    Finds the locations of a motif in a DNA sequence after splicing out the introns.
    Args:
        dna (str): The input DNA sequence.
        motif (str): The motif to search for.
    Returns:
        list: The starting positions of the motif in the DNA sequence after splicing out the introns.
    """
    start = 0
    positions = []
    for base in motif:
        pos = dna.find(base, start)
        if pos == -1:
            return positions
        positions.append(pos + 1)
        start = pos + 1
    return positions

# Problem 31 - Transitions and Transversions
def transitions_transversions(dna1, dna2):
    """
    Computes the ratio of transitions to transversions between two DNA sequences.
    Args:
        dna1 (str): The first DNA sequence.
        dna2 (str): The second DNA sequence.
    Returns:
        float: The ratio of
    """
    transitions = 0
    transversions = 0
    for i in range(len(dna1)):
        if dna1[i] == dna2[i]:
            continue
        if dna1[i] in 'AG' and dna2[i] in 'AG':
            transitions += 1
        elif dna1[i] in 'CT' and dna2[i] in 'CT':
            transitions += 1
        else:
            transversions += 1
    return transitions / transversions

# Problem 32 - Completing a Tree
def completing_tree(n, edges):
    """
    Computes the number of edges required to complete a tree with n nodes.
    Args:
        n (int): The number of nodes in the tree.
        edges (list): A list of edges in the tree.
    Returns:
        int: The number of edges required to complete the tree.
    """
    return n - len(edges) - 1

# Problem 33 - Catalan Numbers and RNA Secondary Structures
def get_catalan_numbers(s, nodes, catalan_memo={}):
    n = int(nodes/2)
    if n <= 1:
        return 1
    if catalan_memo.get((s, nodes),0):
        return catalan_memo[(s, nodes)]
    Cn = 0
    for k in range(1, 2*n, 2):
        a, u, c, g = s[1:k].count("A"), s[1:k].count("U"), s[1:k].count("C"), s[1:k].count("G")
        if a==u and c==g and (s[0], s[k]) in [("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")]:
            Cn += get_catalan_numbers(s[1:k], k-1, catalan_memo) * get_catalan_numbers(s[k+1:], 2*n-k-1, catalan_memo)
    #  Memorize calculated Catalan Numbers values
    catalan_memo[(s, nodes)] = Cn
    return Cn
# Problem 35 - Counting Phylogenetic Ancestors
def count_phylogenetic_ancestors(n):
    """
    Computes the number of internal nodes in a phylogenetic tree with n leaves.
    Args:
        n (int): The number of leaves in the phylogenetic tree.
    Returns:
        int: The number of internal nodes in the phylogenetic tree.
    """
    return n - 2

# Problem 36 - k-Mer Composition
def kmer_composition(dna, k):
    """
    Computes the k-mer composition of a DNA sequence.
    Args:
        dna (str): The input DNA sequence.
        k (int): The length of the k-mers.
    Returns:
        dict: The k-mers in the DNA sequence and their frequencies.
    """
    # Prepare a dict of k-mers from A, T, G, C and their frequencies = 0 initially, sorted lexicographically    
    kmers = enumerate_kmers_lexicographically('ATGC', k)
    kmer_freq = {kmer: 0 for kmer in kmers}
    for i in range(len(dna) - k + 1):
        # Increment the frequency of the k-mer in the dict
        kmer_freq[dna[i:i + k]] += 1
    return kmer_freq

# Problem 37 - Speeding Up Motif Finding
def get_failure_array(s):
    failure_array = [0] * len(s)
    longest_motif_length = 0 # this will store the length of the longest motif
    for i in range(1, len(s)):
        for j in range(1, len(s)-i+1):
            if s[:i] == s[j:j+i]:
                failure_array[j+i-1] = len(s[:i])
                longest_motif_length = len(s[:i])

        # If the length of the longest motif is less than the length of the current motif, break
        if longest_motif_length < len(s[:i]):
            break
    return failure_array

# Problem 38 - Finding a Shared Spliced Motif
def find_shared_spliced_motif(dnas):
    """
    Finds the longest common motif shared by a list of DNA sequences after splicing out the introns.
    Args:
        dnas (list): A list of DNA sequences.
    Returns:
        str: A longest common subsequence of dna string (If more than one solution exists, you may return any one.)
    """
    s = dnas[0]
    t = dnas[1]
    lengths = [[0 for j in range(len(t) + 1)] for i in range(len(s) + 1)] # this will store the length of the longest common subsequence
    #creates array of len(s) containing arrays of len(t) filled with 0
    for i, x in enumerate(s): # enumerate s
        for j, y in enumerate(t): # enumerate t
            if x == y:
                lengths[i + 1][j + 1] = lengths[i][j] + 1 # if the characters are equal, add 1 to the length of the longest common subsequence
            else:
                lengths[i + 1][j + 1] = max(lengths[i + 1][j], lengths[i][j + 1]) # if the characters are not equal, take the maximum of the length of the longest common subsequence
    spliced_motif = ''
    x, y = len(s), len(t)
    while x * y != 0: # while length of s and t is not 0
        if lengths[x][y] == lengths[x - 1][y]:
            x -= 1 # if the length of s is not equal to the length of s - 1, decrease the length of s
        elif lengths[x][y] == lengths[x][y - 1]:
            y -= 1 # if the length of t is not equal to the length of t - 1, decrease the length of t
        else:
            spliced_motif = s[x - 1] + spliced_motif # if the length of s and t are equal, add the character to the spliced_motif
            x -= 1 # decrease the length of s
            y -= 1 # decrease the length of t
    return spliced_motif

# Problem 39 - Ordering Strings of Varying Length Lexicographically
def order_strings_lexicographically(alphabet, n):
    """
    Orders all possible strings of varying length lexicographically.
    Args:
        alphabet (str): The alphabet to use for the strings.
        n (int): The maximum length of the strings.
    Returns:
        list: A list of all possible strings of varying length lexicographically.
    """
    strings = []
    for i in range(1, n + 1):
        strings.extend(enumerate_kmers_lexicographically(alphabet, i))
    # Sort according to the order of the original alphabet
    strings.sort(key=lambda x: [alphabet.index(c) for c in x])
    return strings

# Problem 40 - Maximum Matchings and RNA Secondary Structures
def maximum_matchings(rna):
    """
    Computes the number of maximum matchings in an RNA sequence.
    Args:
        rna (str): The input RNA sequence.
    Returns:
        int: The number of maximum matchings in the RNA sequence.
    """
    from math import factorial
    au_mm = factorial(max(rna.count('A'), rna.count('U')))//factorial(abs(rna.count('A') - rna.count('U')))
    gc_mm = factorial(max(rna.count('G'), rna.count('C')))//factorial(abs(rna.count('G') - rna.count('C')))
    return au_mm * gc_mm

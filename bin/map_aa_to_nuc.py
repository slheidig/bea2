#!/usr/bin/env python3
"""Map an amino-acid alignment onto its unaligned CDS (back-translation).

Adapted from bioenvada v1 bin/mapAAtoNUC.py. Used instead of `cdskit backalign`,
which aborts the whole OG on two things that are normal in bacterial data:
alternative start codons (GTG/TTG/CTG, ~14% of genes here, annotated as M) and
a CDS whose length is not a multiple of three (partial gene at a contig edge).

Usage: map_aa_to_nuc.py <aa_alignment> <unaligned_cds> <og>   ->  <og>.codon.fa
"""
import sys

aa_file, nuc_file, og = sys.argv[1], sys.argv[2], sys.argv[3]

# NCBI translation table 11. The codon->AA map is identical to table 1;
# table 11 differs only in allowing more initiation codons (`starts` below).
gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }

# Table 11 initiation codons: these are translated as M at position 1 only.
starts = {'TTG', 'CTG', 'ATT', 'ATC', 'ATA', 'ATG', 'GTG'}


def read_fasta(path):
    """Parse FASTA, joining wrapped lines and keeping labels verbatim."""
    seqs, order, name = {}, [], None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                order.append(name)
                seqs[name] = []
            elif line and name:
                seqs[name].append(line)
    return order, {k: ''.join(v) for k, v in seqs.items()}


def tripletify(label, nuc_seq):
    nuc_seq = nuc_seq.replace('-', '').replace('.', '').upper()
    extra = len(nuc_seq) % 3
    if extra:
        print('WARNING:', label, 'is not divisable by 3, trimmed', extra, 'trailing nt')
        nuc_seq = nuc_seq[:len(nuc_seq) - extra]
    return [nuc_seq[i:i + 3] for i in range(0, len(nuc_seq), 3)]


def map_gaps(label, triplets, aseq):
    mapped = ''
    j = 0
    for i in range(len(aseq)):
        if aseq[i] in '-.':     # if aa is a gap add '---'
            mapped += '---'
            continue
        if j >= len(triplets):
            raise ValueError('The CDS of ' + label + ' ran out of codons at position ' + str(i + 1))
        codon = triplets[j]     # for every aa add the next 3 nuc
        aa = aseq[i].upper()
        translation = gencode.get(codon)
        if translation is None:
            print('Undefined codons detected!', codon)
        elif translation != aa and aa not in 'X?':
            if j == 0 and aa == 'M' and codon in starts:
                pass            # alternative initiation codon, annotated as M
            else:
                raise ValueError('Missmatch between nucleotide triplet and amino acid detected! '
                                 + label + ' position ' + str(i + 1) + ': ' + aa + ' vs ' + translation)
        mapped += codon
        j += 1
    return mapped               # trailing codons (e.g. the stop) are dropped


aa_order, aa_seqs = read_fasta(aa_file)
nuc_order, nuc_seqs = read_fasta(nuc_file)

with open(og + '.codon.fa', 'w') as out:
    for label in aa_order:
        triplets = tripletify(label, nuc_seqs[label])
        out.write('>' + label + '\n' + map_gaps(label, triplets, aa_seqs[label]) + '\n')

print(og + ':', len(aa_order), 'sequences back-translated to', 3 * len(aa_seqs[aa_order[0]]), 'nt')

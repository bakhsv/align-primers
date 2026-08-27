#!/usr/bin/python3
"""Search for primers in reference sequence and append them to alignment.

Takes in an alignment that has a reference sequence, i.e. the sequence
long enough to have the primer binding sites, and a file with primers.
You should point the script to the reference sequence by providing its
uniform sequence address. For example, if you have an alignment file
called `alignment.fasta` and the reference sequence inside it called `ref`,
it looks like this:

python3 align_primers.py -s alignment.fasta:ref -i coi.primers > aligment_with_primers.fasta

The contents of the `coi.primers` file in this example are:

LCO/HCO ggtcaacaaatcataaagatattgg taaacttcagggtgaccaaaaaatca
jgLCO/jgHCO TNTCNACNAAYCAYAARGAYATTGG TANACYTCNGGRTGNCCRAARAAYCA
dgLCO/dgHCO GGTCAACAAATCATAAAGAYATYGG TAAACTTCAGGGTGACCAAARAAYCA
UCOIR/UCOIF ACWAAYCAYAAAGAYATYGG TAWACTTCDGGRTGRCCRAAAAAYCA 
ZplankF1/ZplankR1 TCTASWAATCATAARGATATTGG TTCAGGRTGRCCRAARAATCA

"""

import sys
import subprocess
import re
import argparse
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import Alignment
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    '-s', '--seqall',
    help='Uniform Sequence Address of a sequence (see Emboss documentation)',
    required=True
)
parser.add_argument(
    '-i', '--infile',
    help='file with primers, see primersearch(1)',
    required=True,
)
parser.add_argument(
    '-o', '--outfile',
    help='output file, default stdout',
    default='-'
)
parser.add_argument(
    '--mismatchpercent',
    help='passed through to primersearch(1)',
    type=int, default=20,
)
parser.add_argument(
    '--print-primersearch',
    action='store_true',
    default=False,
    help='prints primersearch(1) output to stderr'
)


name_re = re.compile(r'Primer name (.+)')
amplimer_re = re.compile(r'Amplimer (\d+)')
hits_re = re.compile(
    r'\s*\S+ hits (forward|reverse) strand at \[?(\d+)\]? with \d+ mismatches'
)
usa_re = re.compile(r'(?P<format>\w+::)?(?P<file>[^:]+)(?P<entry>:\w+)?')


def parse_usa(usa: str) -> dict[str, str]:
    m = usa_re.match(usa)
    if not m:
        raise ValueError(f'{usa} is not in accepted subset of USA.')
    out = m.groupdict()
    if not out['format']:
        out['format'] = 'fasta'
    out['format'] = out['format'].strip(':')
    if out['entry']:
        out['entry'] = out['entry'].strip(':')
    return out


def parse_primersearch(primersearch_result: str) -> list[dict]:
    """Parses the output of `primersearch` program.

    Returns a list of dicts of primer hits with the following keys:
        - name: name of a primer pair.
        - id: one primer pair can have multiple matches. This is the id of
              the current one.
        - f: position at which the forward primer matches
        - r: position at which the reverse primer matches, starting from
             the end of the DNA sequence.
    """
    amplimers = []
    for line in primersearch_result.split('\n'):

        name_match = name_re.match(line)
        if name_match:
            name = name_match.group(1)

        amplimer_match = amplimer_re.match(line)
        if amplimer_match:
            amplimer_id = amplimer_match.group(1)

        hit_match = hits_re.match(line)
        if hit_match:
            if hit_match.group(1) == 'forward':
                curr_ampl = {'name': name,
                             'id': amplimer_id,
                             'f': int(hit_match.group(2))}
            else:  # hit_match.group(2) == 'reverse'
                curr_ampl['r'] = int(hit_match.group(2))
                amplimers.append(curr_ampl)
    return amplimers


comp = {
    'a': 't',
    'c': 'g',
    'g': 'c',
    't': 'a',
    'r': 'y',
    'y': 'r',
    's': 's',
    'w': 'w',
    'k': 'm',
    'm': 'k',
    'b': 'v',
    'v': 'b',
    'd': 'h',
    'h': 'd',
    'n': 'n'
}


def parse_primers(primers_path: str) -> dict[str, tuple[str, str]]:
    """Parses the input primers file."""
    primers = {}
    with open(primers_path) as f:
        for line in f:
            if line[0] != '#':
                if not line.strip():
                    continue
                line = line.split('#', maxsplit=1)[0]
                name, forward, reverse = line.strip().split()
                forward = forward.lower()
                reverse = reverse.lower()
                rrevcomp = ''.join(comp[a] for a in reversed(reverse))
                primers[name] = forward, rrevcomp
    return primers


def convert_primersearch_to_pos(
    f_hit: int, r_hit: int,
    f_primer: str, r_primer: str, n: int
) -> tuple[int, int]:
    """
    Zero-based coordinates of the primer hits *from the start of the sequence*.

    0   1   2   3  [4]  5   6   7   8   9   10  11 [12] 13  14  15  16
    a   t   c   c   a   c   g   t   t   t   g   a   g   c   t   a   a
                    |   |   |   |                   |   |   |   |
                    a   c   g   t                   g   c   t   a
                    F primer                        R primer (revcomp)
    Returns (4, 12)
    """
    f = f_hit - 1
    r = n - r_hit - len(r_primer) + 1
    return f, r


def make_alignment(f_pos: int, r_pos: int,
                   f_primer: str, r_primer: str, n: int) -> str:
    alignment = ''
    alignment += f_pos * '-'
    alignment += f_primer
    alignment += (r_pos - f_pos - len(f_primer)) * '-'
    alignment += r_primer
    alignment += (n - r_pos - len(r_primer)) * '-'
    return alignment


def run_primersearch(seq_file: str, primers_file: str,
                     mismatchpercent: int) -> str:
    ps = subprocess.run(
        [
            'primersearch',
            '-filter',
            '-seqall', seq_file,
            '-infile', primers_file,
            '-mismatchpercent', str(mismatchpercent)
        ],
        capture_output=True
    )
    ps.check_returncode()
    return ps.stdout.decode('utf-8')


def get_seq_index_by_name(alignment: Alignment, name: str) -> int:
    i = 0
    for record in alignment.sequences:
        if record.id == name:
            break
        i += 1
    else:
        raise ValueError(f'Entry {name} not found in file')
    return i


def align_all_primers(
        amplimers: list[dict], primers: dict[str, tuple[str, str]],
        alignment: Alignment, reference_seq_name: str,
) -> tuple[list[str], list[str]]:

    ref_i = get_seq_index_by_name(alignment, reference_seq_name)
    ref_n = len(alignment.sequences[ref_i])

    aligned_primers_names: list[str] = []
    aligned_primers_seqs: list[str] = []
    for amp in amplimers:
        f_ungapped, r_ungapped = convert_primersearch_to_pos(
            f_hit=amp['f'], r_hit=amp['r'],
            f_primer=primers[amp['name']][0], r_primer=primers[amp['name']][1],
            n=ref_n
        )
        f_a = np.where(alignment.indices[ref_i] == f_ungapped)[0][0]
        r_a = np.where(alignment.indices[ref_i] == r_ungapped)[0][0]
        aligned_primers = make_alignment(
            f_a, r_a,
            primers[amp['name']][0],
            primers[amp['name']][1],
            alignment.length
        )
        aligned_primers_seqs.append(aligned_primers)
        aligned_primers_names.append(amp['name'] + '_' + amp['id'])
    return aligned_primers_names, aligned_primers_seqs


def extend_alignment(
        alignment: Alignment, new_names: list[str], new_seqs: list[str],
) -> Alignment:
    seqs = [alignment[i] for i in range(len(alignment))]
    ids = [record.id for record in alignment.sequences]
    seqs.extend(new_seqs)
    seqs = [s.encode('ascii') for s in seqs]
    ids.extend(new_names)
    pseqs, pcoords = Alignment.parse_printed_alignment(seqs)
    ann_seqs = [SeqRecord(Seq(s), id=id) for s, id in zip(pseqs, ids)]
    return Alignment(ann_seqs, pcoords)


if __name__ == '__main__':
    args = parser.parse_args()
    primers: dict[str, tuple[str, str]] = parse_primers(args.infile)
    ps_result = run_primersearch(
        args.seqall, args.infile, args.mismatchpercent
    )
    if args.print_primersearch:
        print(ps_result, file=sys.stderr)
    amplimers: list[dict] = parse_primersearch(ps_result)
    input_seq: dict[str, str] = parse_usa(args.seqall)
    alignment: Alignment = Align.read(input_seq['file'], input_seq['format'])
    aligned_primers_names, aligned_primers_seqs = align_all_primers(
        amplimers, primers, alignment, input_seq['entry']
    )
    al = extend_alignment(
        alignment, aligned_primers_names, aligned_primers_seqs,
    )
    if args.outfile == '-':
        outfile = sys.stdout
    else:
        outfile = open(args.outfile, 'w')
    Align.write(al, outfile, 'fasta')
    outfile.close()

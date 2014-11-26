#!/usr/bin/env python

import argparse
import sys

parser = argparse.ArgumentParser(description='Given BLAST tabular output, find sequences that are very similar to each other (by default identical)')
parser.add_argument('--min_percent_id', default=100.0, type=float)
parser.add_argument('--display_submatches', default=False, action='store_true')
parser.add_argument('--display_contained_sequences', default=False, action='store_true', 
    help='Display matches where one sequences is wholly contained in another')
parser.add_argument('--display_full_matches', default=False)
parser.add_argument('--db_connect_string', help='DB connect string to query length from DB')
parser.add_argument('blast_input_file', type=argparse.FileType(), help='BLAST tabular output file (e.g. from blastall -m 8)')
parser.add_argument('output_file', type=argparse.FileType(), nargs='?', default=sys.stdout)
args = parser.parse_args()

if not (args.display_submatches or args.display_contained_sequences or args.display_full_matches):
    sys.stderr.write('Need to Displayy one of at least submatcehs, contained sequences or full matches')
    sys.exit(1)

session = None
if args.db_connect_string != None:
    from seabass_model import *
    session = get_session(args.db_connect_string)

if (args.display_contained_sequences or args.display_full_matches) and session == None:
    sys.stderr.write('If you want to check that matches are full length or are wholly contained, you need to specify a DB connect string (--db_connect_string)\n')
    sys.exit(1)

pairs_seen = set()
for line in args.blast_input_file:
    fields = line.rstrip().split()
    assert len(fields) == 12, "Invalid BLAST tabular format, expected 12 fields got {} on line {}".format(len(fields), line)
    (query_id, subject_id, percent_id, alignment_length, num_mismatches, num_gap_open, q_start, q_end,
        s_start, s_end, e_val, score) = fields
    if query_id == subject_id:
        # we're not interested in self to self matches :)
        continue
    pair = tuple(sorted((query_id, subject_id)))
    if pair in pairs_seen:
        continue
    (percent_id, e_val, score) = [float(x) for x in (percent_id, e_val, score)]
    (alignment_length, num_mismatches, num_gap_open, q_start, q_end, s_start, s_end) = [int(x) for x in 
        (alignment_length, num_mismatches, num_gap_open, q_start, q_end, s_start, s_end)]
    if percent_id >= args.min_percent_id:
        if session != None:
            q_gene_name = query_id.split('|')[1]
            q_gene_query = session.query(Gene).filter(Gene.name == q_gene_name)
            q_gene = q_gene_query.one()
            s_gene_name = subject_id.split('|')[1]
            s_gene_query = session.query(Gene).filter(Gene.name == s_gene_name)
            s_gene = s_gene_query.one()
        q_hit_length = q_end - q_start + 1
        s_hit_length = s_end - s_start + 1
        got_hit = False
        if not args.display_submatches:
            if (args.display_full_matches and
                q_hit_length == len(q_gene.protein) and s_hit_length == len(s_gene.protein)):
                    got_hit = True
            if (args.display_contained_sequences and ((q_hit_length == len(q_gene.protein) and s_hit_length < len(s_gene.protein)) ^
                (s_hit_length == len(s_gene.protein) and q_hit_length < len(q_gene.protein)))):
                # print q_hit_length, len(q_gene.protein), s_hit_length, len(s_gene.protein)
                # print (q_hit_length == len(q_gene.protein) and s_hit_length < len(s_gene.protein)), (s_hit_length == len(s_gene.protein) and q_hit_length < len(q_gene.protein))
                got_hit = True
        else:
            got_hit = True
        if got_hit:
            args.output_file.write(line)
            pairs_seen.add(pair)

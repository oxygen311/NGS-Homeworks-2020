import networkx as nx
import numpy as np

from os.path import splitext
from Bio import SeqIO
from Bio.Seq import Seq
from operator import itemgetter


class DeBruijinGraph(nx.DiGraph):
    def __init__(self, k):
        super().__init__()
        self.k = k

    def __add_read_one_strand(self, s):
        for i in range(len(s) - self.k):
            start_v = s[i:i + self.k]
            end_v = s[i + 1:i + self.k + 1]
            edge_s = s[i:i + self.k + 1]

            if self.has_edge(start_v, end_v):
                self[start_v][end_v]['coverage'] += 1
            else:
                self.add_edge(start_v, end_v, coverage=1, seq=edge_s)

    def add_reads_from_fasta(self, path):
        reads = SeqIO.parse(path, splitext(path)[1][1:])
        for read in reads:
            self.__add_read_one_strand(read.seq)
            self.__add_read_one_strand(read.seq.reverse_complement())

    def write_dot(self, path):
        out_g = nx.DiGraph()
        node_id, free_id = {}, 1
        for node in self.nodes():
            if self.in_degree(node) == 0 and self.out_degree(node) == 0: continue
            if node not in node_id.keys():
                node_id[node] = f'v{free_id} (forward)'
                node_id[node.reverse_complement()] = f'v{free_id} (rev_comp)'
                free_id += 1
        for start, end, data in self.edges(data=True):
            out_g.add_edge(node_id[start], node_id[end], label=f'"len:{len(data["seq"]) - self.k},cov:{round(data["coverage"], 1)}"')

        nx.drawing.nx_pydot.write_dot(out_g, path)

    def write_fasta(self, path):
        used_seqs = set()
        with open(path, 'w') as f:
            for i, (s, e, seq) in enumerate(self.edges(data='seq')):
                if not seq.reverse_complement() in used_seqs:
                    print('contig', i + 1, file=f)
                    print(seq, file=f)
                used_seqs.add(seq)

    def compress(self):
        def find_path(node):
            path = [node]

            # add backward edges
            cur_node = node
            while self.in_degree(cur_node) == 1:
                cur_node = list(self.in_edges(cur_node))[0][0]
                path.append(cur_node)
                if self.out_degree(cur_node) != 1: break
            path = path[::-1]

            # add forward edges
            cur_node = node
            while self.out_degree(cur_node) == 1:
                cur_node = list(self.out_edges(cur_node))[0][1]
                path.append(cur_node)
                if self.in_degree(cur_node) != 1: break

            return path

        def compress_path(path):
            seq, coverages = '', []

            for start, end in zip(path, path[1:]):
                data = self.get_edge_data(start, end)
                seq += data['seq'] if seq == '' else data['seq'][-1]
                coverages.append(data['coverage'])
                self.remove_edge(start, end)
                self.remove_edge(end.reverse_complement(), start.reverse_complement())
            self.add_edge(path[0], path[-1], coverage=np.mean(coverages), seq=seq)
            self.add_edge(path[-1].reverse_complement(), path[0].reverse_complement(),
                          coverage=np.mean(coverages), seq=seq.reverse_complement())

        for node in self.nodes():
            if self.in_degree(node) == 1 and self.out_degree(node) == 1:
                path = find_path(node)
                compress_path(path)

    def get_coverages(self):
        return list(map(itemgetter(2), self.edges(data='coverage')))

    def get_lengths(self):
        return list(map(len, map(itemgetter(2), self.edges(data='seq'))))

    def remove_weak_edges(self, only_tales, threshold_from_mean):
        coverage_threshold = np.mean(self.get_coverages()) * threshold_from_mean
        length_threshold = np.mean(self.get_lengths()) * threshold_from_mean
        for start, end, data in list(self.edges(data=True)):
            is_tail = self.in_degree(start) == 0 or self.out_degree(start) == 0 or \
                      self.in_degree(end) == 0   or self.out_degree(end)   == 0
            weak = data['coverage'] < coverage_threshold and len(data['seq']) < length_threshold
            if weak and (is_tail or not only_tales):
                self.remove_edge(start, end)
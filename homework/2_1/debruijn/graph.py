from collections import defaultdict
import numpy as np


class DeBruijnGraph:
    def __init__(self, k):
        self.k = k
        self.graph = defaultdict(list)
        self.coverage = defaultdict(int)
        self.edge_labels = dict()
        self.in_edges = defaultdict(list)

    def add_read(self, read):
        for i in range(len(read) - self.k + 1):
            kmer = read[i:i + self.k]
            prefix = kmer[:-1]
            suffix = kmer[1:]
            self.graph[prefix].append(suffix)
            self.coverage[kmer] += 1
            self.edge_labels[(prefix, suffix)] = kmer
            self.in_edges[suffix].append(prefix)

    def build_from_reads(self, reads):
        for read in reads:
            self.add_read(str(read.seq))

    def get_graph(self):
        return dict(self.graph)

    def get_coverage(self):
        return dict(self.coverage)

    def get_edge_labels(self):
        return self.edge_labels

    def compress(self):
        new_graph = defaultdict(list)
        new_edge_labels = dict()
        visited = set()

        def average_cov(path):
            total_len = 0
            weighted_cov = 0
            for i in range(len(path) - 1):
                u, v = path[i], path[i + 1]
                label = self.edge_labels[(u, v)]
                cov = self.coverage[label]
                l = len(label)
                total_len += l
                weighted_cov += cov * l
            return weighted_cov / total_len if total_len > 0 else 0

        for node in list(self.graph.keys()):
            if node in visited:
                continue
            if len(self.in_edges[node]) != 1:
                continue
            path = [node]
            while len(self.graph[path[-1]]) == 1:
                next_node = self.graph[path[-1]][0]
                if len(self.in_edges[next_node]) != 1 or next_node == node:
                    break
                path.append(next_node)
                visited.add(next_node)

            if len(path) > 1:
                compressed_label = path[0] + ''.join(p[-1] for p in path[1:])
                new_graph[path[0]].append(path[-1])
                new_edge_labels[(path[0], path[-1])] = compressed_label
                self.coverage[compressed_label] = average_cov(path)
            else:
                for v in self.graph[node]:
                    new_graph[node].append(v)
                    new_edge_labels[(node, v)] = self.edge_labels[(node, v)]

        self.graph = new_graph
        self.edge_labels = new_edge_labels
        new_in_edges = defaultdict(list)
        for u, vs in new_graph.items():
            for v in vs:
                new_in_edges[v].append(u)
        self.in_edges = new_in_edges

    def _reverse_edges(self, node):
        for u, vs in self.graph.items():
            if node in vs:
                yield u

    def remove_tips(self, quantile=0.3):
        scores = []
        tips = []

        for node in list(self.graph):
            in_deg = len(self.in_edges[node])
            if in_deg > 1:
                continue
            path = [node]
            score = 0
            visited_in_path = set(path)
            while len(self.graph[path[-1]]) == 1:
                nxt = self.graph[path[-1]][0]
                if nxt in visited_in_path:
                    break
                visited_in_path.add(nxt)
                edge = (path[-1], nxt)
                label = self.edge_labels.get(edge)
                if label is None:
                    break
                score += self.coverage[label] * len(label)
                path.append(nxt)
                if len(self.graph[nxt]) != 1:
                    break
            if len(path) > 1:
                scores.append(score)
                tips.append(path)

        if not scores:
            return

        threshold = np.quantile(scores, quantile)

        for path, score in zip(tips, scores):
            if score <= threshold:
                for i in range(len(path) - 1):
                    u, v = path[i], path[i + 1]
                    if v in self.graph[u]:
                        self.graph[u].remove(v)
                        self.edge_labels.pop((u, v), None)

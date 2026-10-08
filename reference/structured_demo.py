"""Replay a huge known tensor power without claiming a new rank improvement."""

import json

from .fields import FiniteField
from .schemes import strassen_scheme
from .structured_certificates import CertificateGraph


def main():
    graph = CertificateGraph(FiniteField(2))
    seed = graph.scheme(strassen_scheme(graph.field).as_tensor_scheme())
    root = graph.power(seed, 25)
    payload = json.dumps(graph.to_dict(root))
    replayed, checked = CertificateGraph.from_dict(json.loads(payload))
    node = replayed.node(checked)
    print(json.dumps({
        'description': '25th tensor power of a supplied known seven-term scheme',
        'axes': node.shape,
        'terms': node.terms,
        'nodes': len(replayed.nodes),
        'serialized_bytes': len(payload.encode()),
        'target_coefficient_000': replayed.tensor_coefficient(node.target, 0, 0, 0),
        'first_A_coefficient': replayed.scheme_coefficient(checked, 0, 0, 0),
        'support': replayed.scheme_support_summary(checked),
        'expanded': False,
        'new_rank_improvement': False,
    }, indent=2))


if __name__ == '__main__':
    main()

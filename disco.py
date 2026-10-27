#!/usr/bin/env python3
"""
DISCO / DISCO-R: decompose multi-copy gene-family trees into single-copy trees.

One script covers DISCO (Willson et al. 2022) and all eight variants studied in
Hakim et al. (wQFM-DISCO). The variants differ along three independent axes:

    --prune   disco | refined     DISCO, or DISCO-R (re-adds the part of the pruned
                                  subtree whose species are missing from the kept side)
    --side    smaller | larger | random
                                  which child of a duplication node is pruned
    costs     --dup-cost/--loss-cost, or --random-costs
                                  weights used when choosing the root

    --variant N (1..8) is a shortcut for the presets in variants.pdf:

    N  pruning  pruned subtree  dup     loss
    1  DISCO    smaller         1       1        (= plain DISCO)
    2  DISCO    larger          1       1
    3  DISCO    random          1       1
    4  DISCO-R  smaller         1       1
    5  DISCO    smaller         1       0.5
    6  DISCO-R  smaller         1       0.5
    7  DISCO    smaller         random  random
    8  DISCO-R  smaller         random  random

Species sets are stored as integer bitmasks, so set unions/intersections are
single machine-word-ish operations instead of Python set allocations.
"""
import argparse
import random

import treeswift
from treeswift import Node

__version__ = "2.0"

VARIANTS = {
    1: dict(prune="disco",   side="smaller", dup=1, loss=1,   random_costs=False),
    2: dict(prune="disco",   side="larger",  dup=1, loss=1,   random_costs=False),
    3: dict(prune="disco",   side="random",  dup=1, loss=1,   random_costs=False),
    4: dict(prune="refined", side="smaller", dup=1, loss=1,   random_costs=False),
    5: dict(prune="disco",   side="smaller", dup=1, loss=0.5, random_costs=False),
    6: dict(prune="refined", side="smaller", dup=1, loss=0.5, random_costs=False),
    7: dict(prune="disco",   side="smaller", dup=1, loss=1,   random_costs=True),
    8: dict(prune="refined", side="smaller", dup=1, loss=1,   random_costs=True),
}


def _popcount(x):
    return x.bit_count() if hasattr(x, "bit_count") else bin(x).count("1")


# --------------------------------------------------------------------------- #
# Species <-> bitmask
# --------------------------------------------------------------------------- #
class SpeciesIndex:
    """Maps leaf labels to single-bit species masks (cached per label)."""

    def __init__(self, gene_to_species=lambda x: x):
        self.gene_to_species = gene_to_species
        self.bit = {}        # species name -> mask
        self.names = []      # bit position -> species name
        self._leaf = {}      # leaf label -> mask

    def __call__(self, label):
        m = self._leaf.get(label)
        if m is None:
            sp = self.gene_to_species(label)
            m = self.bit.get(sp)
            if m is None:
                m = 1 << len(self.names)
                self.bit[sp] = m
                self.names.append(sp)
            self._leaf[label] = m
        return m

    def species(self, mask):
        out, i = [], 0
        while mask:
            if mask & 1:
                out.append(self.names[i])
            mask >>= 1
            i += 1
        return out


def _index(gene_to_species):
    return gene_to_species if isinstance(gene_to_species, SpeciesIndex) else SpeciesIndex(gene_to_species)


# --------------------------------------------------------------------------- #
# Root-scoring costs
# --------------------------------------------------------------------------- #
class Costs:
    """
    Duplication/loss cost of joining species sets a and b (bitmasks):
      disjoint                        -> 0           (speciation)
      a == b                          -> dup
      one side equals the union       -> dup + loss
      otherwise                       -> dup + 2*loss
    Random mode reproduces DISCO-7/8: randint(-2,10), randint(-2,15), randint(-2,20).
    """

    def __init__(self, dup=1, loss=1, random_costs=False, seed=None):
        if random_costs:
            try:
                import numpy
                self._randint = numpy.random.RandomState(seed).randint
            except ImportError:
                self._randint = random.Random(seed).randrange
            self.random = True
        else:
            self.random = False
            vals = [dup, dup + loss, dup + 2 * loss]
            self.c = [int(v) if float(v).is_integer() else v for v in vals]

    def __call__(self, total, a, b):
        if not a & b:
            return 0
        if total == a or total == b:
            k = 0 if a == b else 1
        else:
            k = 2
        if self.random:
            return self._randint(-2, (10, 15, 20)[k])
        return self.c[k]


DEFAULT_COSTS = Costs()


# --------------------------------------------------------------------------- #
# Tree helpers
# --------------------------------------------------------------------------- #
def unroot(tree):
    """Unroot a treeswift tree without contracting (A,B); to A;."""
    if tree.root is None:
        return tree
    if tree.root.num_children() == 2:
        left, right = tree.root.child_nodes()
        if not right.is_leaf():
            right.contract()
        elif not left.is_leaf():
            left.contract()
    tree.is_rooted = False
    return tree


def reroot_on_edge(tree, node):
    if not node.is_root():
        if getattr(node, "edge_length", None) in (None, 0):
            node.edge_length = 1
        tree.reroot(node, length=node.edge_length / 2)


def relabel(tree, gene_to_species=lambda x: x):
    for leaf in tree.traverse_postorder(internal=False):
        leaf.set_label(gene_to_species(leaf.get_label()))
    return tree


def _ensure_binary_root(tree):
    if tree.root.num_children() != 2:
        reroot_on_edge(tree, tree.root.child_nodes()[0])


# --------------------------------------------------------------------------- #
# Core algorithm
# --------------------------------------------------------------------------- #
def remove_in_paralogs(tree, gene_to_species=lambda x: x):
    """Collapse clades containing a single species. Returns number of leaves removed."""
    idx = _index(gene_to_species)
    _ensure_binary_root(tree)
    removed = 0
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.s = idx(node.get_label())
        else:
            s = 0
            for c in node.children:
                s |= c.s
            node.s = s
            if s & (s - 1) == 0:   # exactly one species
                for c in node.children[1:]:
                    node.remove_child(c)
                    removed += 1
    tree.suppress_unifurcations()
    return removed


def get_min_root(tree, gene_to_species=lambda x: x, verbose=False, costs=None):
    """
    ASTRAL-Pro style rooting: returns (best_node, min_score, ties), where the best
    root lies on the edge above best_node.
    """
    idx = _index(gene_to_species)
    score = costs or DEFAULT_COSTS

    if tree.root.num_children() == 0:
        tree.root.s = idx(tree.root.get_label())
        return tree.root, 0, []

    _ensure_binary_root(tree)
    tree.resolve_polytomies()

    # down pass
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.down = idx(node.get_label())
            node.d_score = 0
        else:
            if len(node.children) != 2:
                raise Exception("Vertex has more than 2 children")
            left, right = node.children
            node.down = left.down | right.down
            node.d_score = left.d_score + right.d_score + score(node.down, left.down, right.down)

    # children of the root
    root = tree.root
    left, right = root.children
    left.up, left.u_score = right.down, right.d_score
    right.up, right.u_score = left.down, left.d_score

    min_score = left.u_score + left.d_score + score(left.up | left.down, left.up, left.down)
    if not left.is_leaf():
        best_root = left
    elif not right.is_leaf():
        best_root = right
    else:                       # two-leaf tree: keep rooting
        best_root = root
    ties = [best_root]

    # up pass (preorder, internal nodes below the root's children)
    for node in tree.traverse_preorder(leaves=False):
        if node is root or node is left or node is right:
            continue
        parent = node.parent
        pc = parent.children
        other = pc[0] if pc[0] is not node else pc[1]
        node.up = parent.up | other.down
        node.u_score = parent.u_score + other.d_score + score(node.up, parent.up, other.down)
        total = node.u_score + node.d_score + score(node.up | node.down, node.up, node.down)
        if total == min_score:
            ties.append(node)
        if total < min_score:
            min_score, best_root, ties = total, node, [node]

    if verbose:
        print("Best root had score", min_score, "there were", len(ties), "ties.")
    return best_root, min_score, ties


def tag(tree, gene_to_species=lambda x: x):
    """Tag internal nodes 'S' (speciation) or 'D' (duplication) under current rooting."""
    idx = _index(gene_to_species)
    tree.suppress_unifurcations()
    tree.resolve_polytomies()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.s = idx(node.get_label())
            node.n_dup = 0
        else:
            left, right = node.children
            node.s = left.s | right.s
            node.n_dup = left.n_dup + right.n_dup
            if left.s & right.s:
                node.tag = "D"
                node.n_dup += 1
            else:
                node.tag = "S"
    tree.n_dup = tree.root.n_dup
    tree.species_index = idx


def _restricted_copy(node, excl):
    """
    Copy of the subtree at `node` keeping only leaves whose species are not in
    `excl`, with unifurcations suppressed (edge lengths summed exactly as
    treeswift's suppress_unifurcations would). Returns None if nothing survives.
    """
    def rec(n):
        # returns (copied_node, extra_edge_from_collapsed_ancestors, collapsed?) or None
        if not (n.cur & ~excl):
            return None
        if not n.children:
            c = Node(label=n.label, edge_length=n.edge_length)
            c.s = c.cur = n.cur
            return c, None, False
        kept = [r for r in map(rec, n.children) if r is not None]
        if len(kept) == 1:
            c, extra, _ = kept[0]
            if n.edge_length is not None:
                extra = n.edge_length if extra is None else extra + n.edge_length
            return c, extra, True
        c = Node(label=n.label, edge_length=n.edge_length)
        cur = 0
        # treeswift moves a collapsed child to the end of its parent's child list;
        # mirror that so the output newick is identical, not just the topology
        for k, extra, collapsed in [r for r in kept if not r[2]] + [r for r in kept if r[2]]:
            _apply_extra(k, extra)
            c.add_child(k)
            cur |= k.cur
        c.cur = cur
        return c, None, False

    r = rec(node)
    if r is None:
        return None
    c, extra, _ = r
    _apply_extra(c, extra)
    return c


def _apply_extra(node, extra):
    if extra is not None:
        node.edge_length = (0 if node.edge_length is None else node.edge_length) + extra


def _as_tree(root):
    t = treeswift.Tree()
    t.root = root
    root.parent = None
    return t


def decompose(tree, single_tree=False, prune="disco", side="smaller", rng=None):
    """
    Decompose a tagged tree (run `tag` first) at its duplication nodes.

    prune: 'disco'   - cut off the chosen child subtree
           'refined' - DISCO-R: also re-attach, under the duplication node, the part
                       of the cut subtree whose species are absent from the kept side
    side:  'smaller' | 'larger' | 'random' - which child is cut
    Returns the list of resulting trees (the main tree last).
    """
    refined = prune == "refined"
    if side == "random" and rng is None:
        rng = random.Random()
    out = []
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.cur = node.s
            continue
        if node.tag == "D":
            left, right = node.children
            if side == "smaller":
                delete = left if _popcount(left.s) < _popcount(right.s) else right
            elif side == "larger":
                delete = left if _popcount(left.s) > _popcount(right.s) else right
            else:
                delete = left if rng.random() < 0.5 else right
            keep = right if delete is left else left

            sd = _restricted_copy(delete, keep.cur) if refined else None
            node.remove_child(delete)
            if not single_tree:
                t = _as_tree(delete)
                t.suppress_unifurcations()
                out.append(t)
            if sd is not None:
                node.add_child(sd)
        cur = 0
        for c in node.children:
            cur |= c.cur
        node.cur = cur
    tree.suppress_unifurcations()   # duplication nodes left with one child
    out.append(tree)
    return out


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def main(args):
    if args.delimiter is not None:
        d, n = args.delimiter, args.nth_delimiter
        gene_to_species = lambda x: d.join(x.split(d)[:n])
    else:
        gene_to_species = lambda x: x
    idx = SpeciesIndex(gene_to_species)

    costs = Costs(args.dup_cost, args.loss_cost, args.random_costs, args.seed)
    side_rng = random.Random(args.seed)

    if args.output is None:
        stem, ext = args.input.rsplit(".", 1)
        output = stem + "-decomp." + ext
    else:
        output = args.output

    outgroup_file = args.input.rsplit(".", 1)[0] + "_outgroups.txt"
    og = open(outgroup_file, "w") if args.outgroups else None

    with open(args.input) as fi, open(output, "w") as fo:
        for i, line in enumerate(fi, 1):
            tree = treeswift.read_tree_newick(line)
            if isinstance(tree, list):
                assert len(tree) == 0, "Could not interpret {} on line {:d} as a tree".format(line, i)
                continue

            if args.remove_in_paralogs:
                n_par = remove_in_paralogs(tree, idx)

            root, score, ties = get_min_root(tree, idx, costs=costs)
            reroot_on_edge(tree, root)
            tag(tree, idx)

            if args.verbose:
                print("Tree ", i, ": Tree has ", _popcount(tree.root.s), " species.", sep="")
                if args.remove_in_paralogs:
                    print(n_par, "in-paralogs removed prior to rooting/scoring.")
                if _popcount(tree.root.s) < 2:
                    print("Uninformative")
                elif tree.n_dup == 0:
                    print("Single-Copy")
                else:
                    small = min(tree.root.children, key=lambda c: _popcount(c.s))
                    print("Best root had score ", score, " with ", tree.n_dup,
                          " non-terminal" if args.remove_in_paralogs else "",
                          " duplications; there were ", len(ties), " ties.\nOutgroup: {",
                          ",".join(sorted(map(str, idx.species(small.s)))), "}", sep="")

            if args.no_decomp:
                out = [tree]
            else:
                out = [t for t in decompose(tree, args.single_tree, args.prune, args.side, side_rng)
                       if t.num_nodes(internal=False) >= args.minimum]

            for t in out:
                if not args.no_decomp:
                    unroot(t)
                if not args.keep_labels:
                    relabel(t, gene_to_species)
                t.suppress_unifurcations()
                fo.write(t.newick() + "\n")

            if args.verbose:
                print("Decomposition strategy outputted", len(out), "tree(s) with minimum size", args.minimum, ".\n")

            if og is not None:
                og_tree = treeswift.read_tree_newick(line)
                root, score, ties = get_min_root(og_tree, idx, costs=costs)
                reroot_on_edge(og_tree, root)
                tag(og_tree, idx)
                if _popcount(og_tree.root.s) >= 2 and og_tree.n_dup >= 1:
                    og.write("Tree " + str(i) + ":\n")
                    for t in ties:
                        reroot_on_edge(og_tree, t)
                        tag(og_tree, idx)
                        small = min(og_tree.root.children, key=lambda c: _popcount(c.s))
                        og.write("{" + ",".join(sorted(map(str, idx.species(small.s)))) + "}\n")
    if og is not None:
        og.close()


def build_parser():
    p = argparse.ArgumentParser(
        description="DISCO / DISCO-R v" + __version__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Variants (--variant):\n" + "\n".join(
            "  {}: prune={prune}, side={side}, {}".format(
                k, "random costs" if v["random_costs"] else "dup={dup}, loss={loss}".format(**v), **v)
            for k, v in VARIANTS.items()))
    p.add_argument("-i", "--input", required=True, help="Input tree list file")
    p.add_argument("-o", "--output", help="Output tree list file")
    p.add_argument("-d", "--delimiter", help="Delimiter separating species name from rest of leaf label")
    p.add_argument("-n", "--nth-delimiter", type=int, help="Split on nth delimiter (only works with -d)")
    p.add_argument("-m", "--minimum", type=int, default=4, help="Minimum tree size outputted")
    p.add_argument("-v", "--verbose", action="store_true", help="Enables verbose output")
    p.add_argument("--keep-labels", action="store_true",
                   help="Keep original leaf labels instead of species labels (only relevant with -d)")
    p.add_argument("-s", "--single-tree", "--single_tree", dest="single_tree", action="store_true",
                   help="Only output single large tree")
    p.add_argument("--no-decomp", action="store_true", help="Outputs rooted trees without decomposition")
    p.add_argument("--outgroups", action="store_true", help="Output outgroups to file (including ties)")
    p.add_argument("--remove-in-paralogs", "--remove_in_paralogs", dest="remove_in_paralogs",
                   action="store_true", help="Remove in-paralogs before rooting/scoring tree")

    g = p.add_argument_group("variants")
    g.add_argument("-V", "--variant", type=int, choices=sorted(VARIANTS),
                   help="Preset from Hakim et al. (overrides the options below)")
    g.add_argument("-r", "--prune", choices=["disco", "refined"], default="disco",
                   help="'refined' = DISCO-R")
    g.add_argument("--side", choices=["smaller", "larger", "random"], default="smaller",
                   help="Which child subtree of a duplication is pruned")
    g.add_argument("--dup-cost", type=float, default=1, help="Duplication cost for rooting")
    g.add_argument("--loss-cost", type=float, default=1, help="Loss cost for rooting")
    g.add_argument("--random-costs", action="store_true", help="Random dup/loss costs (variants 7/8)")
    g.add_argument("--seed", type=int, help="Random seed (side=random, --random-costs)")
    return p


if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()
    if args.variant is not None:
        for k, v in VARIANTS[args.variant].items():
            setattr(args, k if k not in ("dup", "loss") else k + "_cost", v)
    if args.delimiter is None:
        if args.nth_delimiter is not None:
            parser.error("Cannot set -n without a delimiter")
        if args.keep_labels:
            parser.error("Cannot use --keep-labels without a delimiter")
    elif args.nth_delimiter is None:
        args.nth_delimiter = 1
    if args.single_tree and args.no_decomp:
        parser.error("Cannot combine --single-tree and --no-decomp")
    if not args.verbose and args.remove_in_paralogs:
        print("--remove-in-paralogs is meaningless without --verbose, as it does not change the optimal rooting. "
              "It may also slow the program.")
    main(args)

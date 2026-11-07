"""Two-track H-H coupling: Bally-Rablen for general H-H, DU4 for geminal CH2 2JHH.

Sources (both verified from the authors' own files, not from secondary summaries)
---------------------------------------------------------------------------------
Bally, T.; Rablen, P. R. J. Org. Chem. 2011, 76, 4818-4830, doi 10.1021/jo200513q.
  The author script `mk_input_file` (SI zip, jo200513q_si_002) writes:
      #n B3LYP/gen nmr=(FCOnly,ReadAtoms) int=nobasistransform geom=check
      <every non-H element> 0 / 6-31G(d,p) / ****
      H 0 / 8 s functions + 1 p / ****
      atoms=H
  and the paper states the FC terms "should be scaled by 0.9155" (p. 4827).
  SI Table S4a decodes the method as "H F H D N": mixed basis on H only, FC only,
  ReadAtoms/H, default SCF convergence and grid, integral=nobasistransform.

Kutateladze, A. G.; Mukhina, O. A. J. Org. Chem. 2014, 79, 8397-8406,
doi 10.1021/jo501781b.  SI S2-S4: geometry B3LYP/6-31G(d); FC only with
int=nobasistransform; basis C 4-31G, N/O/F/Li/Be/B 6-31G, Si/S/P/Cl/As/Se/Br/I 3-21G*,
H = optimised 8s1p; geminal CH2 fits per environment (SI S4, geminal panel).

Design note
-----------
Environment-specific scaling of geminal CH2 couplings is a *judgement call* (ring size,
pi-adjacency, extrapolation), so this module does NOT pick a class automatically at run
time.  It computes the two FC matrices, exports the raw per-conformer FC values and a
scaling worksheet listing every candidate formula; the class is chosen and recorded
during post-processing (by the analyst or the assistant), with the rationale written
into the same worksheet.
"""
ROUTE = "#p B3LYP/Gen 6D NMR=(FCOnly) Int=NoBasisTransform G09Defaults"

BALLY_FACTOR = 0.9155
BALLY_SOURCE = "https://doi.org/10.1021/jo200513q"
DU4_SOURCE = "https://doi.org/10.1021/jo501781b"

# Heavy atoms: the author script writes 6-31G(d,p) for every non-H element present.
BALLY_BASIS = {"C": "6-31G(d,p)", "O": "6-31G(d,p)"}
BALLY_DEFAULT = "6-31G(d,p)"
DU4_BASIS = {"C": "4-31G", "N": "6-31G", "O": "6-31G", "F": "6-31G",
             "Li": "6-31G", "Be": "6-31G", "B": "6-31G",
             "Si": "3-21G*", "S": "3-21G*", "P": "3-21G*", "Cl": "3-21G*",
             "As": "3-21G*", "Se": "3-21G*", "Br": "3-21G*", "I": "3-21G*"}

DU4_CLASSES = {
    "sp3_general": dict(slope=1.240, intercept_hz=0.57, fit="1.240 fc + 0.57",
                        n_points=23, rmsd_hz=0.13),
    "sp3_adjacent_sp2": dict(slope=1.072, intercept_hz=-1.66, fit="1.072 fc - 1.66",
                             n_points=10, rmsd_hz=0.09),
    "small_ring": dict(slope=1.249, intercept_hz=1.05, fit="1.249 fc + 1.05",
                       n_points=25, rmsd_hz=0.14),
    "sp2": dict(slope=1.337, intercept_hz=0.52, fit="1.337 fc + 0.52",
                n_points=20, rmsd_hz=0.24),
}

# Hydrogen basis exactly as printed in the author's mk_input_file: eight s primitives
# (four decontracted 6-31G plus four added tight s) and one p with exponent 1.1.
BALLY_H_BLOCK = [
    "S 1 1.00", "0.1873113696D+02 1.0",
    "S 1 1.00", "0.2825394365D+01 1.0",
    "S 1 1.00", "0.6401216923D+00 1.0",
    "S 1 1.00", "0.1612777588D+00 1.0",
    "S 1 1.00", "0.5619341088D+02 1.0",
    "S 1 1.00", "0.1685802326D+03 1.0",
    "S 1 1.00", "0.5057406979D+03 1.0",
    "S 1 1.00", "0.1517222094D+04 1.0",
    "P 1 1.00", "0.1100000000D+01 1.0",
]
DU4_H_BLOCK = [
    "S 1 1.00", "14.573150211 1.0", "S 1 1.00", "3.9686235317 1.0",
    "S 1 1.00", ".35279124965 1.0", "S 1 1.00", ".18072966380 1.0",
    "S 1 1.00", "240.39092912 1.0", "S 1 1.00", "385.73555302 1.0",
    "S 1 1.00", "1685.6845423 1.0", "S 1 1.00", "2387.1673902 1.0",
    "P 1 1.00", ".49109009605 1.0",
]


def basis_block(kind, elements):
    """Gen section; returns (text, notes)."""
    if kind not in ("bally", "du4"):
        raise ValueError("Unknown coupling basis kind " + str(kind))
    notes, lines = [], []
    for el in sorted(set(elements)):
        if el == "H":
            continue
        if kind == "bally":
            name = BALLY_BASIS.get(el, BALLY_DEFAULT)
            if el not in BALLY_BASIS:
                notes.append("%s uses the author-script default %s" % (el, name))
        else:
            name = DU4_BASIS.get(el)
            if name is None:
                raise ValueError("DU4 SI defines no basis for element " + el +
                                 "; add it to DU4_BASIS with a documented source")
        lines += ["%s 0" % el, name, "****"]
    h_block = BALLY_H_BLOCK if kind == "bally" else DU4_H_BLOCK
    lines += ["H 0"] + list(h_block) + ["****"]
    return "\n".join(lines) + "\n\n", notes


def typed_graph(mol2_path):
    """adjacency plus bond orders from the MOL2 bond block"""
    from pathlib import Path
    text = Path(mol2_path).read_text()
    if "@<TRIPOS>BOND" not in text:
        raise ValueError("Missing MOL2 bond connectivity: " + str(mol2_path))
    graph, order = {}, {}
    for line in text.split("@<TRIPOS>BOND", 1)[1].split("@<TRIPOS>", 1)[0].splitlines():
        p = line.split()
        if not p:
            continue
        a, b = int(p[1]), int(p[2])
        graph.setdefault(a, set()).add(b)
        graph.setdefault(b, set()).add(a)
        order[frozenset((a, b))] = p[3] if len(p) > 3 else "1"
    return graph, order


def ring_sizes(graph, start):
    """Ring sizes through `start`: shortest path between two of its neighbours that
    avoids `start`, plus the two bonds back to `start`."""
    sizes = set()
    nbrs = list(graph.get(start, ()))
    for i, a in enumerate(nbrs):
        for b in nbrs[i + 1:]:
            dist = {a: 0}
            frontier = [a]
            while frontier:
                nxt = []
                for u in frontier:
                    for v in graph.get(u, ()):
                        if v == start or v in dist:
                            continue
                        dist[v] = dist[u] + 1
                        nxt.append(v)
                frontier = nxt
            if b in dist:
                sizes.add(dist[b] + 2)
    return sizes


def describe_ch2(graph, order, geometry):
    """Advisory description of every CH2: ring sizes, pi-adjacency, formula candidates.
    No class is chosen here; the choice is recorded during post-processing."""
    out = []
    for i, atom in enumerate(geometry, 1):
        if atom[0] != "C":
            continue
        hs = sorted(j for j in graph.get(i, ()) if geometry[j - 1][0] == "H")
        if len(hs) != 2:
            continue
        heavy = [j for j in graph.get(i, ()) if geometry[j - 1][0] != "H"]
        sp2 = [j for j in heavy
               if any(order.get(frozenset((j, k))) in ("2", "ar") for k in graph.get(j, ()))]
        rings = sorted(ring_sizes(graph, i))
        if rings and min(rings) in (3, 4):
            hint = "small_ring"
            note = "3/4-membered ring: DU4 has no ring fit, values are extrapolative"
        elif sp2:
            hint = "sp3_adjacent_sp2"
            note = "adjacent to sp2 carbon/oxygen"
        else:
            hint = "sp3_general"
            note = "plain sp3 CH2"
        out.append(dict(carbon_atom=i, H_a=hs[0], H_b=hs[1],
                        ring_sizes=",".join(map(str, rings)),
                        adjacent_sp2_atoms=",".join(map(str, sp2)),
                        suggested_class=hint, suggestion_note=note,
                        candidates={k: dict(v) for k, v in DU4_CLASSES.items()}))
    return out


def j_from(fc_hz, cls):
    fit = DU4_CLASSES[cls]
    return fit["slope"] * fc_hz + fit["intercept_hz"]


def bally_j(fc_hz):
    return BALLY_FACTOR * fc_hz


def describe():
    return {"route": ROUTE,
            "bally": {"factor": BALLY_FACTOR,
                      "basis": "non-H 6-31G(d,p) (author default) + H 8s1p",
                      "source": BALLY_SOURCE,
                      "method_code": "H F H D N (SI Table S4a)",
                      "scope": "general_HH_single_linear_scale"},
            "du4": {"classes": DU4_CLASSES, "basis": DU4_BASIS,
                    "hydrogen_basis": "optimised 8s1p (SI S3)",
                    "source": DU4_SOURCE, "scope": "geminal_2JHH_of_CH2_only"},
            "note": "geminal class is chosen at post-processing, not automatically"}

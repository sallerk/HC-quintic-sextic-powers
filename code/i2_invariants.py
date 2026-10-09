"""i2_invariants.py -- Hodge classes on X_F^k (middle Kunneth part H^{2m}(X_F)^{(x)k}) for very general F, d = 5, 6.

Model (inputs, as in the note: Theorem A, the theorem of Cattani-Deligne-Kaplan and the sector model):
 * H^{2m}(X_F) = C h^m (+) (+)_c V_c^prim, and V_c^prim ~ (x)_s P_s(c_s) (x) L_c as Hodge structures up to the
   one-dimensional twist L_c of type (m - sigma, n0 + sigma - 1), sigma = sum_s c_s / d (reps in [0,d-1]),
   n0 = #{s: c_s = 0}; P(0) := trivial piece of dimension d-1 and type (0,0); P(g) of type (g-1, d-1-g).
   (Hodge numbers of this model are checked against the Jacobian ring below.)
 * Hodge classes of a very general member = vectors with finite monodromy orbits (CDK) of type (km,km), i.e.
   invariants of the identity component M of the Zariski closure of the monodromy. Per index s (independent):
     d=5: M_s = SL3 x SL3 acting on P(1)=std_A, P(4)=std_A^*, P(2)=std_B, P(3)=std_B^*, trivially on P(0);
     d=6: M_s = SL4 x Sp4 acting on P(2)=std, P(4)=std^*, P(3)=std(Sp4); trivially on P(0), P(1), P(5) (finite).
 * Invariant subspaces have a single Hodge type (checked by the formulas below), so a Galois orbit of sector tuples
   contributes rational Hodge classes iff every member has type (km,km).
Categories: PAIRING (only GL-type invariants and P(0)/h factors), CMFIN (d=6: uses P(+-1) factors, no det),
EXCEPTIONAL (some SL_n factor with #std != #dual, i.e. a determinant).
Output: i2_invariants.json and a summary.
"""
import itertools, json, math, sys
from collections import Counter, defaultdict
from functools import lru_cache


# ---------------- invariant dimensions via the Weyl character trick
def poly_mul(A, B):
    C = defaultdict(int)
    for ka, va in A.items():
        for kb, vb in B.items():
            C[tuple(x + y for x, y in zip(ka, kb))] += va * vb
    return {k: v for k, v in C.items() if v}


@lru_cache(maxsize=None)
def inv_SL(n, p, q):
    """dim of SL_n-invariants in V^{(x)p} (x) V*^{(x)q}: sum over t of coefficients at t(1,..,1) of
    chi * prod_{i<j} (1 - x_j/x_i) (GL_n weights)."""
    e = lambda i, s=1: tuple(s if j == i else 0 for j in range(n))
    std = {e(i): 1 for i in range(n)}
    dual = {e(i, -1): 1 for i in range(n)}
    P = {tuple([0] * n): 1}
    for _ in range(p):
        P = poly_mul(P, std)
    for _ in range(q):
        P = poly_mul(P, dual)
    for i in range(n):
        for j in range(i + 1, n):
            f = {tuple([0] * n): 1, tuple(1 if t == j else (-1 if t == i else 0) for t in range(n)): -1}
            P = poly_mul(P, f)
    return sum(v for k, v in P.items() if len(set(k)) == 1)


@lru_cache(maxsize=None)
def inv_Sp4(r):
    std = {(1, 0): 1, (-1, 0): 1, (0, 1): 1, (0, -1): 1}
    P = {(0, 0): 1}
    for _ in range(r):
        P = poly_mul(P, std)
    for a in [(1, -1), (1, 1), (2, 0), (0, 2)]:  # positive roots
        P = poly_mul(P, {(0, 0): 1, (-a[0], -a[1]): -1})
    return P.get((0, 0), 0)


# ---------------- per-index data
def per_s(d, labels):
    """labels: tuple of labels (ints mod d) of the P-factors at index s (copies that are h^m are omitted).
    Returns (inv_dim, (type_p, type_q) as Fractions*12 ints, category flags)."""
    c = Counter(labels)
    if d == 5:
        p1, q1, p2, q2, n0 = c[1], c[4], c[2], c[3], c[0]
        i = inv_SL(3, p1, q1) * inv_SL(3, p2, q2) * 4 ** n0
        if i == 0:
            return 0, None, None
        # types (x3 to stay integral): U_A (a=0): (q1, p1); U_B (a=1): ((p2+2q2)/3, (2p2+q2)/3)
        t = (3 * q1 + (p2 + 2 * q2), 3 * p1 + (2 * p2 + q2))
        exc = (p1 != q1) or (p2 != q2)
        return i, t, ("EXC" if exc else "PAIR")
    if d == 6:
        p, q, r, n1, n5, n0 = c[2], c[4], c[3], c[1], c[5], c[0]
        i = inv_SL(4, p, q) * inv_Sp4(r) * 4 ** (n1 + n5) * 5 ** n0
        if i == 0:
            return 0, None, None
        # x12: SL4 (a=1): ((p+3q)/4, (3p+q)/4) -> x12: (3(p+3q), 3(3p+q)); Sp4: (r/2, r/2) -> (6r, 6r);
        # finite: (n5, n1) -> (12 n5, 12 n1)
        t = (3 * (p + 3 * q) + 6 * r + 12 * n5, 3 * (3 * p + q) + 6 * r + 12 * n1)
        exc = p != q
        cat = "EXC" if exc else ("CMFIN" if n1 + n5 > 0 else "PAIR")
        return i, t, cat
    raise ValueError


SCALE = {5: 3, 6: 12}


def L_type(d, m, c):
    sigma = sum(c) // d
    n0 = sum(1 for x in c if x == 0)
    return (m - sigma, n0 + sigma - 1)


def model_hodge(d, c):
    """Hodge numbers of (x)_s P_s(c_s) (x) L_c (model), as dict p -> dim, p = second index."""
    m = len(c) - 1
    Lt = L_type(d, m, c)
    out = Counter({0: 1})
    for cs in c:
        if cs == 0:
            fac = {0: d - 1}
        else:
            fac = {0: cs - 1, 1: d - 1 - cs}  # (1,0)-part contributes q=0, (0,1)-part q=1
        new = Counter()
        for a, va in out.items():
            for b, vb in fac.items():
                if vb:
                    new[a + b] += va * vb
        out = new
    return {q + Lt[1]: v for q, v in out.items() if v}


def jac_hodge(d, c):
    m = len(c) - 1
    hilb = lambda e: min(e + 1, 2 * d - 3 - e) if 0 <= e <= 2 * d - 4 else 0
    opts = [[(e, hilb(e)) for e in range(2 * d - 3) if (e - (cs - 2)) % d == 0 and hilb(e)] for cs in c]
    out = Counter()
    for combo in itertools.product(*opts):
        e = sum(x[0] for x in combo)
        if (e + 2 * m + 2) % d:
            continue
        out[(e + 2 * m + 2) // d - 1] += math.prod(x[1] for x in combo)
    return dict(out)


def sectors(d, m):
    return [c for c in itertools.product(range(d), repeat=m + 1) if sum(c) % d == 0 and any(c)] + \
           [tuple([0] * (m + 1))]


def analyse(d, m, k, max_examples=6):
    secs = sectors(d, m)
    H = "h"
    items = secs + [H]
    units = [u for u in range(1, d) if math.gcd(u, d) == 1]
    S = SCALE[d]
    Ltypes = {c: L_type(d, m, c) for c in secs}
    # per multiset of sectors (copies unordered); evaluate dim and type; then Galois orbits of multisets
    res = {}
    for ms in itertools.combinations_with_replacement(range(len(items)), k):
        tup = [items[i] for i in ms]
        dim = 1
        tx, ty = 0, 0
        cats = set()
        for s in range(m + 1):
            labels = tuple(c[s] for c in tup if c != H)
            i, t, cat = per_s(d, labels)
            if i == 0:
                dim = 0
                break
            dim *= i
            tx += t[0]
            ty += t[1]
            cats.add(cat)
        if dim == 0:
            continue
        for c in tup:
            if c == H:
                tx += S * m
                ty += S * m
            else:
                tx += S * Ltypes[c][0]
                ty += S * Ltypes[c][1]
        nperm = math.factorial(k)
        for v in Counter(ms).values():
            nperm //= math.factorial(v)
        hodge = (tx == ty == S * k * m)
        cat = "EXC" if "EXC" in cats else ("CMFIN" if "CMFIN" in cats else "PAIR")
        res[ms] = (dim * nperm, hodge, cat)
    # Galois orbits
    idx = {c: i for i, c in enumerate(items)}
    def act(u, ms):
        return tuple(sorted(idx[H] if items[i] == H else idx[tuple((u * x) % d for x in items[i])] for i in ms))
    seen = set()
    tot = Counter()
    examples = defaultdict(list)
    for ms in res:
        if ms in seen:
            continue
        orb = {act(u, ms) for u in units}
        seen |= orb
        dims = [res.get(o, (0, False, None)) for o in orb]
        alldim = sum(x[0] for x in dims)
        cat = res[ms][2]
        tot["finite_orbit_dim"] += alldim
        if all(x[1] for x in dims) and all(o in res for o in orb):
            tot["hodge_" + cat] += alldim
            if cat != "PAIR" and len(examples[cat]) < max_examples:
                examples[cat].append([[("h" if items[i] == H else items[i]) for i in o] for o in sorted(orb)][:1])
        elif any(x[1] for x in dims):
            tot["partial_orbits_" + cat] += 1
    return dict(tot), {k_: v for k_, v in examples.items()}


def R(d, m):
    if d in (3, 4, 6):
        return sum(math.factorial(m + 1) // (math.factorial(k) ** 2 * math.factorial(m + 1 - 2 * k))
                   * (d - 2) ** (2 * k) * (d - 1) ** (m + 1 - 2 * k) for k in range((m + 1) // 2 + 1))
    return (d - 1) ** (m + 1)


def main():
    out = {"checks": {}, "results": {}}
    # invariant-theory sanity
    assert inv_SL(3, 1, 1) == 1 and inv_SL(3, 3, 0) == 1 and inv_SL(3, 2, 2) == 2 and inv_SL(3, 2, 0) == 0
    assert inv_SL(4, 4, 0) == 1 and inv_SL(4, 2, 2) == 2 and inv_SL(4, 3, 0) == 0
    assert inv_Sp4(2) == 1 and inv_Sp4(4) == 3 and inv_Sp4(1) == 0
    # model Hodge numbers vs Jacobian ring
    nbad = 0
    for d in (5, 6):
        for m in (1, 2, 3):
            for c in sectors(d, m):
                if any(c):
                    if model_hodge(d, c) != jac_hodge(d, c):
                        nbad += 1
    out["checks"]["model_vs_jacobian_mismatches"] = nbad
    print("model vs Jacobian ring mismatches:", nbad)
    plan = [(5, 1, 4), (5, 2, 4), (5, 3, 3), (6, 1, 4), (6, 2, 4), (6, 3, 3)]
    for d, m, kmax in plan:
        for k in range(1, kmax + 1):
            tot, ex = analyse(d, m, k)
            key = f"d={d},m={m},k={k}"
            out["results"][key] = {"totals": tot, "examples": ex}
            print(key, tot, ("examples: " + json.dumps(ex)) if ex else "")
            sys.stdout.flush()
        print(f"  k=1 check: R({d},{m}) = {R(d, m)} (+1 for h^m)")
    with open(sys.argv[1] if len(sys.argv) > 1 else "i2_invariants.json", "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()

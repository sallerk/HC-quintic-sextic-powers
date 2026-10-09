"""formulas.py -- checks of the closed formulas (exact integer arithmetic).

(G) Galois lemma (d = 5): for every multiset D in (Z/5)^x of size <= 10,
      sum_{g in D} (2<u g> - 5) = 0 for all units u   <=>   mult(g) = mult(-g) for all g.
(N) N_k = #{x in ((Z/5)^x)^k : sum x = 0} = (4^k + 4 (-1)^k)/5, k <= 8 (enumeration).
(P) dim of the Hodge classes in H^{2m}(X_F) (x) H^{2m}(X_F), F very general, d = 5:
      (1 + 4^{m+1})^2 + sum_{j=0}^{m-1} C(m+1, j) N_{m+1-j} 16^j   vs the two invariant counts (m = 1, 2, 3).
usage: python formulas.py counts2.json i2_invariants.json out.json
"""
import itertools, json, math, sys


def galois_check(maxsize=10):
    units = [1, 2, 3, 4]
    bad = 0
    n = 0
    for size in range(1, maxsize + 1):
        for D in itertools.combinations_with_replacement(units, size):
            n += 1
            cond = all(sum(2 * ((u * g) % 5) - 5 for g in D) == 0 for u in units)
            bal = all(D.count(g) == D.count((-g) % 5) for g in units)
            bad += cond != bal
    return n, bad


def N_enum(k):
    return sum(1 for x in itertools.product(range(1, 5), repeat=k) if sum(x) % 5 == 0)


def main():
    out = {}
    out["galois_multisets"], out["galois_failures"] = galois_check()
    out["N"] = {k: [N_enum(k), (4 ** k + 4 * (-1) ** k) // 5] for k in range(1, 9)}
    out["N_ok"] = all(a == b for a, b in out["N"].values())
    c2 = json.load(open(sys.argv[1]))
    i2 = json.load(open(sys.argv[2]))
    rows = {}
    for m in (1, 2, 3):
        f = (1 + 4 ** (m + 1)) ** 2 + sum(math.comb(m + 1, j) * (4 ** (m + 1 - j) + 4 * (-1) ** (m + 1 - j)) // 5 * 16 ** j
                                          for j in range(m))
        key = f"d=5,m={m},k=2"
        r2 = c2[key]["hodge_total"]
        ri = i2["results"][key]["totals"]
        ri = ri.get("hodge_PAIR", 0) + ri.get("hodge_EXC", 0) + ri.get("hodge_CMFIN", 0)
        rows[m] = dict(formula=f, counts2=r2, i2=ri, ok=(f == r2 == ri))
    out["XxX"] = rows
    print(json.dumps(out, indent=1))
    with open(sys.argv[3], "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()

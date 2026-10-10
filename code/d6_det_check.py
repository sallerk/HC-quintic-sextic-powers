"""d6_det_check.py -- d = 6: tally of the Hodge classes in sectors with determinants (D(c) nonempty).

Uses the model of counts2.py, imported from the same folder (monodromy SL4 x Sp4 on P(2), P(4) / P(3), trivial on
P(0), P(1), P(5)). For every sector tuple whose finite-orbit space is of Hodge type and has nonempty determinant
content D(c), it records D2, D4 (numbers of determinants on P(2), P(4)), n5 - n1 (number of entries 5 = -1 minus
number of entries 1) and the dimension. Checks: n5 - n1 = 2 (D2 - D4) for each such tuple (the balance (2) in the
proof of Theorem D(ii), which is the condition of Lemma 8.3), and the totals against the "det" column of Section 9
(2, 356262, 32).

Usage: python d6_det_check.py <out.json>
"""
import itertools, json, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import counts2 as c2  # noqa: E402


def tally(d, factors):
    ms = [len(f) - 1 for f in factors]
    n = sum(ms)
    labels = sorted(set(l for f in factors for l in f))
    secs = [[None] + c2.sectors(d, len(f)) for f in factors]
    pre = []
    for f, m, sl in zip(factors, ms, secs):
        rows = []
        for c in sl:
            if c is None:
                rows.append((m, tuple(() for _ in labels)))
            else:
                rows.append((m - sum(c) // d, tuple(tuple(c[s] for s, l in enumerate(f) if l == lab)
                                                    for lab in labels)))
        pre.append(rows)
    out = {}
    exc_total = 0
    balance_fail = 0
    for idx in itertools.product(*[range(len(sl)) for sl in secs]):
        dim, p, ok = 1, 0, True
        dets = []
        n1 = n5 = D2 = D4 = 0
        single = True
        for i, j in enumerate(idx):
            p += pre[i][j][0]
        for li in range(len(labels)):
            ent = []
            for i, j in enumerate(idx):
                ent.extend(pre[i][j][1][li])
            r = c2.label_info(d, tuple(sorted(ent)))
            if r[0] == 0:
                ok = False
                break
            dim *= r[0]
            if len(r[1]) != 1:
                single = False
            p += r[1][0][0]
            dets.extend(r[2])
            n1 += r[3]; n5 += r[4]; D2 += r[5]; D4 += r[6]
        if not ok or not single or p != n or not dets:
            continue
        exc_total += dim
        key = f"D2={D2},D4={D4},n5-n1={n5 - n1}"
        out[key] = out.get(key, 0) + dim
        if n5 - n1 != 2 * (D2 - D4):
            balance_fail += 1
    return {"det_total": exc_total, "balance_fail": balance_fail, "by_type": out}


def main():
    jobs = {
        "X_F^4, m=1": (6, [[0, 1]] * 4, 2),
        "X_F^4, m=2": (6, [[0, 1, 2]] * 4, 356262),
        "X(G,G,G,G,G',G'')": (6, [[0, 0, 0, 0, 1, 2]], 32),
        "X(G,G,G,G')": (6, [[0, 0, 0, 1]], 0),
        "X(G,G,G,G,G')": (6, [[0, 0, 0, 0, 1]], None),
    }
    res = {}
    for name, (d, f, expected) in jobs.items():
        t = time.time()
        r = tally(d, f)
        r["expected_det"] = expected
        r["match"] = (expected is None) or (r["det_total"] == expected)
        r["sec"] = round(time.time() - t, 1)
        res[name] = r
        print(name, r)
        sys.stdout.flush()
    with open(sys.argv[1], "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()

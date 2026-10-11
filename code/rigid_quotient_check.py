r"""Cohomology of the quotients of Section 6 of the note ("Rigid quotients").

Setting: nu >= 2, n >= 1 with nu | n + 2, G a binary form of degree n + 2 with distinct zeros, and the
curve C = {y^nu + G(u,v) = 0} in P(1, 1, (n+2)/nu), a cyclic cover of P^1 of degree nu, totally ramified
over the n + 2 zeros of G. The group mu_nu acts on y. Let mu_nu^{n,1} be the kernel of the product map
mu_nu^n -> mu_nu, and let calG = mu_nu^{n,1} x| S_n act on C^n (mu_nu^n factorwise, S_n by permuting the
factors).

The script computes, in exact arithmetic, the Betti and Hodge numbers of H^*(C^n, Q)^calG, the
cohomology of Omega_G = C^n / calG. It uses the Kunneth basis of H^*(C^n) built from a basis of H^*(C)
adapted to the characters of mu_nu and to the Hodge decomposition: H^0, H^2 (trivial character), and
for each character gamma = 1, ..., nu - 1 the eigenspace P(gamma) of H^1, of dimension n with
(1,0)-part of dimension (n+2) gamma / nu - 1. The invariants of mu_nu^{n,1} are spanned by the
Kunneth basis tensors on which mu_nu^n acts through a character trivial on mu_nu^{n,1}, i.e. whose
factors all have the same character. The dimension of the S_n-invariants of their span is the average
over S_n of the trace, where S_n permutes the tensor factors with the Koszul sign.

It compares the Euler characteristic of Omega_G with that of the normal cyclic cover of degree nu of
P^n branched along n + 2 hyperplanes in general position, nu e(P^n - A) + e(A), where A is the union
of the hyperplanes and e(A) is computed by inclusion-exclusion. By Proposition 6.2 of the note the two
varieties are isomorphic.

Cases: (nu, n) = (5, 3), the curves y^5 + G(u,v) = 0 of Theorem B, and (nu, n) = (3, 4), the curves
w^3 + G(u,v) = 0 of Section 9.

Usage: python rigid_quotient_check.py <out.json>
"""
import itertools
import json
import sys
from fractions import Fraction
from math import comb

CASES = [(5, 3), (3, 4)]


def invariant_cohomology(nu, n):
    """Return {(degree, p, q): dimension} for H^*(C^n)^calG."""
    k = n + 2  # number of branch points
    # basis of H^*(C): (degree, character, p, q)
    basis = [(0, 0, 0, 0), (2, 0, 1, 1)]
    for g in range(1, nu):
        p10 = Fraction(k * g, nu) - 1
        assert p10.denominator == 1
        p10 = int(p10)
        basis += [(1, g, 1, 0)] * p10 + [(1, g, 0, 1)] * (n - p10)
    nb = len(basis)
    assert sum(1 for b in basis if b[0] == 1) == (nu - 1) * n  # = 2 * genus

    def koszul_sign(sigma, t):
        """Sign of moving the tensor factors of t, factor i going to position sigma[i]."""
        odd = [i for i in range(n) if basis[t[i]][0] % 2]
        tgt = [sigma[i] for i in odd]
        inv = sum(1 for a in range(len(tgt)) for b in range(a + 1, len(tgt)) if tgt[a] > tgt[b])
        return -1 if inv % 2 else 1

    perms = list(itertools.permutations(range(n)))
    counts = {}
    for t in itertools.product(range(nb), repeat=n):
        if len({basis[i][1] for i in t}) != 1:
            continue  # not invariant under mu_nu^{n,1}
        key = tuple(sum(basis[i][j] for i in t) for j in (0, 2, 3))
        s = sum(koszul_sign(sigma, t) for sigma in perms if all(t[sigma[i]] == t[i] for i in range(n)))
        if s:
            counts[key] = counts.get(key, Fraction(0)) + Fraction(s, len(perms))
    out = {}
    for key, c in counts.items():
        assert c.denominator == 1
        if c:
            out[key] = int(c)
    return out


def euler_cover(nu, n):
    """nu e(P^n - A) + e(A) for the union A of n + 2 hyperplanes in general position in P^n."""
    k = n + 2
    e_a = sum((-1) ** (j + 1) * comb(k, j) * (n + 1 - j) for j in range(1, min(k, n) + 1))
    return nu * (n + 1 - e_a) + e_a, e_a


def main():
    out_path = sys.argv[1]
    results = []
    for nu, n in CASES:
        h = invariant_cohomology(nu, n)
        betti = [sum(v for (d, p, q), v in h.items() if d == deg) for deg in range(2 * n + 1)]
        euler = sum((-1) ** i * b for i, b in enumerate(betti))
        e_cov, e_a = euler_cover(nu, n)
        hodge = {str(deg): {f'{p},{q}': v for (d, p, q), v in sorted(h.items(), reverse=True) if d == deg}
                 for deg in range(2 * n + 1) if any(d == deg for (d, _, _) in h)}
        results.append({'nu': nu, 'n': n, 'betti': betti, 'hodge': hodge, 'euler': euler,
                        'euler_cover': e_cov, 'euler_A': e_a, 'agree': euler == e_cov})
        print(f'nu={nu}, n={n}: Betti numbers of Omega_G: {betti}, Euler characteristic {euler}')
        for deg, hs in hodge.items():
            print(f'  H^{deg}: ' + ', '.join(f'h^{pq}={v}' for pq, v in hs.items()))
        print(f'  cyclic cover of P^{n} branched along {n + 2} hyperplanes: e = {e_cov} (e(A) = {e_a}); '
              f'agree: {euler == e_cov}')
    summary = {'cases': results, 'all_agree': all(r['agree'] for r in results)}
    with open(out_path, 'w', encoding='utf-8', newline='\n') as fh:
        json.dump(summary, fh, indent=1)
        fh.write('\n')
    print(f'all agree: {summary["all_agree"]}')


if __name__ == '__main__':
    main()

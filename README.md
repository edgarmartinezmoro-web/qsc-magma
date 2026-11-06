# qsc-magma

Computational material for the paper

> E. Martínez-Moro and M. de los Ríos de Antonio,
> *Quantum synchronizable codes, polynomial and quasi-polynomial shifts* (in preparation).

The paper extends Fujiwara's construction of quantum synchronizable codes (QSCs) from cyclic
codes to constacyclic, polycyclic and quasi-polycyclic codes. This repository holds the searches
behind its section *Experimental data*, their raw output, and two scripts that rebuild every
reported code from the stored data and check its parameters.

## Contents

| Path | Content |
|---|---|
| `magma/qp_qsc_search.m` | Magma library: quasi-polycyclic codes, annihilator dual, Gram twist, polynomial orders and the search `RunGalindo`. |
| `magma/run*.m` | The search scripts that produced `results/`. |
| `magma/verify_examples.m` | Rebuilds and checks the codes of the index-1 tables, the binomial table and the `q = 5` table. |
| `magma/verify_q7_one.m` | Rebuilds and checks one code of the `q = 7` table (`N = 98`), chosen on the command line. |
| `results/res*.txt` | Raw output of the searches (format below). |
| `python/` | The two exhaustive searches done in Python (Tables `tab:consta` and `tab:exp:aly`). |

All Magma scripts load the library with `load "qp_qsc_search.m";`, so they must be run from
`magma/`.

## Requirements

- Magma. All results were obtained with Magma V2.29-7. The library uses only standard
  intrinsics such as `LinearCode`, `Dual`, `MinimumDistance` and `Factorization`.
- Python 3 with `numpy` and `sympy`, for `python/`.

## Verifying the codes of the paper

```
cd magma
magma -b verify_examples.m
for n in 1 2 3 4 5 6 7 8; do
  nohup nice -n 10 magma -b i:=$n verify_q7_one.m > verify_q7_one_$n.log 2>&1 &
done
```

For each code the scripts print a line that starts with `OK` or `MISMATCH`, and at the end the
number of mismatches, which should be 0. `verify_examples.m` checks 39 codes in about half a
minute. In `verify_q7_one.m` the index is 1 = `[[98,64,7]]`, 2 = `[[98,58,7]]`, 3 = `[[98,54,7]]`,
4 = `[[98,50,7]]`, 5 = `[[98,42,7]]`, 6 = `[[98,48,5]]`, 7 = `[[98,32,5]]`, 8 = `[[98,14,5]]`; the
binomial codes take a few seconds and the trinomial codes between about 6 and 14 hours each,
because of the exact distances. The results can be collected with
`grep -h "^OK\|MISMATCH" verify_q7_one_*.log`. With Magma V2.29-7 all 47 codes gave `OK`.

For the structured family `C = <(g, k g), (0, h)>` in `A^2`, with `A = F_q[x]/<f>`, the scripts
check

- the hypotheses `g | h` and `h^2 | f`;
- that the annihilator dual `C°` is contained in `C` (Lemma `lem:exp:annGalindo`);
- that `Dual(Psi(C) G)` is contained in `Psi(C)`, so that `(Psi(C), Psi(C) G)` is a Euclidean CSS
  pair (Theorem `theo:cssanh`);
- `K`, `dX = d(Psi(C))` and `dZ = d(Psi(C) G)`;
- the synchronization data: the best proper divisor `g1` of `g` with `ord(g/g1) >= 3`, its distance
  `d1 = d(<g1>)` and the window `ord(g/g1)`.

For the index-1 codes of Aydin, Liu and Yoshino, `verify_examples.m` checks the annihilator
condition `<f/h2> ⊆ <g1>`, the Euclidean CSS pair `(Psi(<g1>), Psi(<h2>) G)`, the values of `K`,
`d1` and `d2`, and the QSC families `t | f/(h2 g1)` with `ord(t)` and `K' = K - deg t`.

Here `Psi` is coefficient expansion and `G = I ⊗ (ε(x^{i+j}))` is the Hankel Gram matrix of the
annihilator form. The section headers printed by `verify_examples.m` use the labels of an earlier
draft: `sec:trinoex`, `tab:trinoqsc` and `tab:trinoqsc5` are now Section `sec:exp:index1` and
Tables `tab:exp:aly` and `tab:exp:aly5`; `tab:binomcap` is `tab:exp:binomcap`; `tab:trino25` is
`tab:exp:q5`; and `ex:qp50` is item 3 of Example `ex:exp:best`.

## Which file supports which part of the paper

Labels refer to the LaTeX source of the paper.

| Paper | Files |
|---|---|
| Table `tab:consta` (cyclic versus constacyclic codes) | `python/cc.py` |
| Tables `tab:exp:aly` and `tab:exp:aly5` (index 1, trinomial codes) | `magma/verify_examples.m`; constacyclic competitors with `python/comp.py` |
| Table `tab:exp:binomcap`, `q = 3` | `res1b_q3.txt`, `res1b_q3_m24.txt`, `res3_q3.txt` |
| Table `tab:exp:binomcap`, `q = 5` | `res0_test.txt` (`m = 10`), `res2b_q5_m25_QC.txt` (`m = 25`), `res3_q5.txt` (`m = 50`) |
| Table `tab:exp:binomcap`, `q = 7` | `res1b_q7_m14.txt`, `res1b_q7_m21.txt`, `res1b_q7_m28.txt`, `res1b_q7_m42.txt`, `res2b_q7_m49_QC.txt` |
| Table `tab:exp:q5`, trinomial columns | `res2b_q5_m25_tri.txt`, `res3_q5_tri25.txt`, `res4_q5_m25_tri.txt` (three independent runs) |
| Table `tab:exp:q5`, columns for `x^25 - 1` | `res2b_q5_m25_QC.txt` |
| Table `tab:exp:q7`, binomial rows | `res2b_q7_m49_QC.txt`; checked by `magma/verify_q7_one.m` |
| Table `tab:exp:q7`, trinomial rows | `res5_q7_m49_t1.txt`, `res5_q7_m49_t2.txt`, `res5_q7_m49_t3.txt` (`f = x^49 + x^7 + b`, `b = 1, 2, 3`); checked by `magma/verify_q7_one.m` |
| Example `ex:exp:best`, items 1 and 2 | `res5_q7_m49_t1.txt`, `res5_q7_m49_t3.txt` |
| Example `ex:exp:best`, item 3 | `res4_q5_m25_tri.txt` |

## Format of the result files

Each code found by `RunGalindo` and saved with `SaveGalindo` takes four lines:

```
[[N,K]] family d=.. dX=.. dZ=.. d1=.. win=.. rho=.. f=...
  g=...
  h=...
  k=...
```

`family` is `QC`, `QconstaC` or `QP(rho=r)`, `d = min(dX, dZ)`, `d1` is the distance of the marker
code `<g1>` in the synchronization block, `win` is `ord(g/g1)` (the window is
`a_l + a_r < win`), and `rho` is the twisted acyclicity of `f`. The code is
`C = <(g, k g), (0, h)>` in `(F_q[x]/<f>)^2`. In some files the variable of the polynomial ring is
printed as `$.1`; read it as `x`.

## Reproducing the searches

Run the scripts from `magma/`:

```
magma -b run0_test.m        # quick test; also gives q = 5, m = 10
magma -b run1b_consta.m     # quasi-cyclic and quasi-constacyclic codes, q = 3, 7
magma -b run2b_dml.m        # lengths of Du-Ma-Liu: N = 50 (q = 5) and N = 98 (q = 7)
magma -b run3_s2.m          # m = l p^2, and repeated-root trinomials with m = 25
magma -b run4_tri.m         # targeted search with repeated-root trinomials, m = 25
magma -b run5_q7_light.m    # repeated-root trinomials, m = 49, one output file per trinomial
```

The output files are written to `magma/`; we moved them to `results/` by hand. Some practical
points:

- The searches sample random multipliers `k`, and random pairs `(g, h)` when there are more than
  `maxpairs` of them. No seed was fixed, so a rerun gives equivalent but not identical output.
  Add `SetSeed(1);` at the top of a script to make it deterministic. The stored `f, g, h, k` rebuild
  every reported code exactly.
- `SaveGalindo` appends to its output file, so remove old `res*.txt` files before a rerun.
- Running times go from seconds (`run0_test.m`) to many hours (`run3_s2.m`, `run4_tri.m`,
  `run5_q7_light.m`). Each search runs on a single core.
- `run5_q7_light.m` loops over the 36 trinomials `x^49 + x^(7i) + b`. The paper uses the first
  three (`res5_q7_m49_t1.txt` to `t3.txt`); the other trinomials have not been run to the end.

## Python material

These scripts use `numpy` and `sympy` and do not need Magma. Run them from `python/`.

| Script | Purpose |
|---|---|
| `lib.py` | Polynomial arithmetic over prime fields, linear algebra mod `p` and minimum distance; used by `comp.py`. |
| `cc.py` | Cyclic versus constacyclic codes (Table `tab:consta`). Usage: `python3 cc.py q lambdas n0 n1`, for example `python3 cc.py 3 1,2 4 16`. |
| `comp.py` | Exhaustive search of cyclic and constacyclic competitors (all `λ`, all divisors of `x^n - λ`) for Table `tab:exp:aly`. Usage: `python3 comp.py n q K' D1 D2 minord`, for example `python3 comp.py 10 13 1 4 4 3`. |

## Citation

If you use this material, please cite the paper above.

## License

The code (`magma/` and `python/`) is licensed under the GNU General Public License v3.0; see
`LICENSE`. The data in `results/` are licensed under the Creative Commons Attribution 4.0
International License (CC BY 4.0); see `results/LICENSE.md`.

Copyright (c) 2026 Edgar Martínez-Moro and Miguel de los Ríos de Antonio.

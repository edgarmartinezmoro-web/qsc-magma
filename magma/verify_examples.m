// =====================================================================
//  verify_examples.m
//  Rebuilds every code reported in the paper from its stored data and checks
//  the claimed properties. Run it from the magma/ directory:
//
//      magma -b verify_examples.m
//
//  For each code it prints OK or MISMATCH; the summary at the end gives the
//  number of mismatches (it should be 0).
//
//  Checked properties
//   * Structured family (Lemma lem:exp:annGalindo):
//       C = <(g, k g), (0, h)> in A^2, A = F_q[x]/<f>
//       - hypotheses g | h, h^2 | f;
//       - annihilator dual C° contained in C;
//       - Dual(Psi(C) G) contained in Psi(C)       (Theorem theo:cssanh);
//       - K = 2 dim C - 2m, dX = d(Psi(C)), dZ = d(Psi(C) G);
//       - synchronization block: best proper divisor g1 of g with ord(g/g1) >= 3,
//         d1 = d(<g1>) and window ord(g/g1)  (same rule as the search).
//   * Codes of Aydin-Liu-Yoshino (index 1):
//       - g1 | f, h2 | f, annihilator condition <f/h2> in <g1> (i.e. g1 | f/h2);
//       - Dual(Psi(<h2>) G) contained in Psi(<g1>);
//       - K, d1 = d(Psi(<g1>)), d2 = d(Psi(<h2>) G);
//       - for each listed t: t | f/(h2 g1), ord(t) and K' = K - deg t.
// =====================================================================

load "qp_qsc_search.m";

function VerifyGalindo(label, f, g, h, k, exp)
    P := Parent(f); m := Degree(f); N := 2*m;
    hyp := IsDivisibleBy(h, g) and IsDivisibleBy(f, h^2);
    X  := QPCode([ [g, (k*g) mod f], [P!0, h] ], f, 2);
    X2 := Twist(X, f, 2);
    ann := AnnDual(X, f, 2) subset X;
    css := Dual(X2) subset X;
    K  := 2*Dimension(X) - N;
    dX := MinimumDistance(X);
    dZ := MinimumDistance(X2);
    found := false; d1 := 0; w1 := 0;
    for g1 in MonicDivisors(g) do
        if Degree(g1) eq Degree(g) then continue; end if;
        w := PolyOrder(g div g1);
        if w lt 3 then continue; end if;
        d := MinimumDistance(PolyCode(g1, f));
        if (not found) or d gt d1 or (d eq d1 and w gt w1) then
            d1 := d; w1 := w; found := true;
        end if;
    end for;
    got := < K, dX, dZ, d1, w1 >;
    ok := hyp and ann and css and (got eq exp);
    printf "%-6o %-32o [[%o,%o]]  dX=%o dZ=%o d1=%o window<%o  (hyp %o, C°<=C %o, CSS %o)\n",
           ok select "OK" else "MISMATCH", label, N, K, dX, dZ, d1, w1, hyp, ann, css;
    if not ok then printf "       expected <K,dX,dZ,d1,window> = %o, got %o\n", exp, got; end if;
    return ok;
end function;

function VerifyALY(label, q, fseq, g1seq, h2seq, exp, ts)
    F := GF(q); P := PolynomialRing(F);
    f := P!fseq; g1 := P!g1seq; h2 := P!h2seq; n := Degree(f);
    hyp := IsDivisibleBy(f, g1) and IsDivisibleBy(f, h2);
    u := f div h2;
    ann := IsDivisibleBy(u, g1);
    X1 := PolyCode(g1, f);
    X2 := Twist(PolyCode(h2, f), f, 1);
    css := Dual(X2) subset X1;
    K  := Dimension(X1) + Dimension(X2) - n;
    d1 := MinimumDistance(X1);
    d2 := MinimumDistance(X2);
    okt := true;
    for T in ts do
        t := P!T[1];
        if not (IsDivisibleBy(u div g1, t) and PolyOrder(t) eq T[2] and K - Degree(t) eq T[3]) then
            okt := false;
        end if;
    end for;
    got := < K, d1, d2 >;
    ok := hyp and ann and css and okt and (got eq exp);
    printf "%-6o %-32o K=%o d1=%o d2=%o  (hyp %o, ann %o, CSS %o, families %o)\n",
           ok select "OK" else "MISMATCH", label, K, d1, d2, hyp, ann, css, okt;
    if not ok then printf "       expected <K,d1,d2> = %o, got %o\n", exp, got; end if;
    return ok;
end function;

fails := 0;
t0 := Cputime();

print "\n== Codes of Aydin-Liu-Yoshino (Section sec:trinoex, Tables tab:trinoqsc and tab:trinoqsc5) ==";
fails +:= (VerifyALY("ALY [[10,2]]_13", 13, [12, 11, 0, 0, 0, 0, 0, 0, 0, 0, 1], [8, 10, 1, 1, 1], [12, 6, 10, 1, 1], <2,4,4>, [ <[2, 1],12,1> ]) select 0 else 1);
fails +:= (VerifyALY("ALY [[15,3]]_13", 13, [10, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [9, 8, 4, 6, 5, 2, 1], [5, 9, 6, 0, 12, 5, 1], <3,5,5>, []) select 0 else 1);
fails +:= (VerifyALY("ALY [[16,10]]_13", 13, [9, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 4, 9, 1], [6, 11, 11, 1], <10,3,3>, [ <[8, 1],4,9> ]) select 0 else 1);
fails +:= (VerifyALY("ALY [[18,10]]_13", 13, [12, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [8, 6, 8, 12, 1], [10, 6, 1, 5, 1], <10,4,4>, []) select 0 else 1);
fails +:= (VerifyALY("ALY [[24,16]]_13", 13, [9, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [12, 12, 2, 6, 1], [12, 1, 2, 7, 1], <16,4,4>, []) select 0 else 1);
fails +:= (VerifyALY("ALY [[27,17]]_7", 7, [4, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 1, 3, 1, 3, 1], [4, 2, 3, 4, 6, 1], <17,4,4>, [ <[2, 0, 0, 1],18,14>, <[2, 4, 3, 2, 5, 1],16806,12> ]) select 0 else 1);
fails +:= (VerifyALY("ALY [[24,18]]_17", 17, [14, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [12, 16, 6, 1], [5, 16, 11, 1], <18,3,3>, []) select 0 else 1);
fails +:= (VerifyALY("ALY [[25,18]]_5", 5, [3, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [2, 1, 3, 1], [4, 0, 0, 2, 1], <18,3,3>, []) select 0 else 1);
fails +:= (VerifyALY("ALY [[36,28]]_13", 13, [9, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [5, 0, 2, 1], [1, 12, 2, 6, 2, 1], <28,3,4>, []) select 0 else 1);
fails +:= (VerifyALY("ALY [[60,32]]_3", 3, [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [1, 1, 2, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1], [2, 0, 2, 2, 0, 0, 0, 2, 0, 2, 2, 1, 0, 2, 0, 2, 1], <32,5,6>, [ <[1, 2, 1],6,30>, <[1, 0, 2, 1],26,29>, <[1, 2, 0, 2, 1, 1],78,27>, <[1, 2, 1, 1, 0, 2, 1],182,26>, <[1, 0, 1, 2, 2, 0, 0, 0, 1],546,24> ]) select 0 else 1);
fails +:= (VerifyALY("ALY [[72,48]]_3", 3, [2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [2, 1, 0, 0, 1, 1, 2, 1, 2, 1, 1, 0, 1], [2, 2, 0, 0, 1, 2, 2, 2, 2, 2, 1, 0, 1], <48,5,5>, [ <[1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 1],484,38> ]) select 0 else 1);

// ------------------------------ q = 3 ------------------------------
F := GF(3); P<x> := PolynomialRing(F);

print "\n== Table tab:binomcap, q = 3 ==";
fails +:= (VerifyGalindo("binomcap q=3 m=6", x^6 + 1, x^2 + 1, x^2 + 1, x^5 + x^3 + x^2 + 2, <4,2,2,1,4>) select 0 else 1);
fails +:= (VerifyGalindo("binomcap q=3 m=12", x^12 + 1, x^4 + 1, x^4 + 1, x^10 + x^8 + x^7 + x^6 + x^2 + x, <8,2,2,2,8>) select 0 else 1);
fails +:= (VerifyGalindo("binomcap q=3 m=18", x^18 + 2, x^8 + 2*x^6 + 2*x^2 + 1, x^8 + 2*x^6 + 2*x^2 + 1, x^17 + 2*x^16 + 2*x^15 + 2*x^14 + x^11 + 2*x^10 + 2*x^9 + 2*x^8 + x^7 + 2*x^5 + x^4 + x^3 + x, <4,3,3,3,6>) select 0 else 1);
fails +:= (VerifyGalindo("binomcap q=3 m=24", x^24 + 1, x^8 + 1, x^8 + 1, x^22 + x^20 + 2*x^19 + x^18 + 2*x^17 + x^16 + x^15 + x^12 + 2*x^10 + 2*x^9 + 2*x^8 + 2*x^7 + 2*x^6 + x^5 + 2*x^4 + x^3 + x^2 + 1, <16,2,2,2,16>) select 0 else 1);
fails +:= (VerifyGalindo("binomcap q=3 m=36", x^36 + 1, x^16 + x^12 + x^4 + 1, x^16 + x^12 + x^4 + 1, 2*x^33 + 2*x^32 + 2*x^31 + 2*x^30 + 2*x^28 + 2*x^23 + x^21 + x^20 + x^19 + 2*x^18 + 2*x^17 + 2*x^15 + x^13 + x^10 + x^9 + 2*x^7 + 2*x^6 + 2*x^5 + 2*x^4 + 2*x^3, <8,3,3,3,72>) select 0 else 1);

// ------------------------------ q = 5 ------------------------------
F := GF(5); P<x> := PolynomialRing(F);

print "\n== Table tab:trino25 and Example ex:qp50 (repeated-root trinomials, m = 25) ==";
fails +:= (VerifyGalindo("trino25 K=10 [ex:qp50(b)]", x^25 + x^10 + 4, x^10 + 2*x^7 + 3*x^5 + x^4 + 3*x^2 + 1, x^10 + 2*x^7 + 3*x^5 + x^4 + 3*x^2 + 1, 3*x^24 + x^23 + x^20 + 2*x^16 + 3*x^15 + 3*x^13 + 4*x^11 + x^10 + 4*x^9 + 3*x^7 + 3*x^6 + 4*x^5 + 2*x^4 + x^2 + 3, <10,6,6,6,4>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=12 [ex:qp50(a)]", x^25 + 2*x^10 + 4, x^9 + 4*x^8 + x^7 + 3*x^6 + 2*x^5 + x^4 + 3*x^3 + 2*x^2 + 4*x + 1, x^10 + 4*x^7 + 3*x^5 + 4*x^4 + x^2 + 1, 3*x^24 + 3*x^23 + x^22 + 4*x^21 + x^19 + 4*x^18 + 4*x^17 + 2*x^16 + 4*x^14 + 2*x^13 + 3*x^10 + 4*x^7 + 3*x^5 + 3*x^4 + 2*x^3 + 4*x^2 + x + 3, <12,6,6,5,20>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=14", x^25 + x^10 + 3, x^9 + x^8 + x^7 + 3*x^6 + 3*x^5 + 4*x^4 + x + 1, x^9 + x^8 + x^7 + 3*x^6 + 3*x^5 + 4*x^4 + x + 1, 2*x^24 + 2*x^23 + x^19 + 3*x^18 + 4*x^17 + 3*x^13 + 4*x^12 + 3*x^11 + 3*x^10 + 3*x^9 + x^8 + 4*x^7 + 3*x^5 + 2*x^4 + 2*x^3 + 3*x^2 + 3*x + 1, <14,6,6,5,4>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=16", x^25 + x^10 + 3, x^8 + 2*x^7 + 3*x^6 + x^5 + 4*x^4 + 3*x^3 + 3*x^2 + 3*x + 4, x^9 + x^8 + x^7 + 3*x^6 + 3*x^5 + 4*x^4 + x + 1, 3*x^24 + 2*x^23 + 4*x^22 + x^21 + 2*x^20 + 2*x^19 + x^18 + x^17 + 2*x^15 + x^14 + 4*x^13 + 3*x^12 + 3*x^9 + x^8 + 3*x^6 + 4*x^4 + 2*x^3 + 2*x^2 + 2*x + 1, <16,6,6,5,4>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=18", x^25 + x^10 + 3, x^7 + 3*x^5 + 4*x^3 + 3*x + 2, x^9 + x^8 + x^7 + 3*x^6 + 3*x^5 + 4*x^4 + x + 1, 3*x^24 + 4*x^23 + x^22 + 3*x^21 + 4*x^19 + 2*x^18 + 4*x^17 + x^16 + 4*x^15 + x^14 + x^12 + 4*x^11 + 4*x^10 + x^9 + x^8 + 3*x^7 + x^6 + 2*x^5 + 2*x^4 + 3*x^3 + x^2 + 4*x, <18,6,6,3,124>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=20", x^25 + 2*x^10 + 4, x^6 + 4*x^5 + 3*x^4 + 2*x^3 + 2*x^2 + 3*x + 4, x^9 + 3*x^8 + 4*x^7 + x^6 + 3*x^5 + 2*x^4 + x + 3, 3*x^24 + 3*x^23 + 4*x^22 + 2*x^21 + x^20 + 2*x^17 + 4*x^15 + 4*x^12 + 4*x^11 + 4*x^9 + 4*x^8 + 3*x^5 + x^4 + 3*x^3 + 3*x^2 + 2*x + 4, <20,6,6,3,124>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=22", x^25 + x^10 + 3, x^5 + 3*x^4 + 3*x^3 + 4*x^2 + x + 4, x^9 + x^8 + x^7 + 3*x^6 + 3*x^5 + 4*x^4 + x + 1, x^24 + 3*x^23 + 2*x^22 + x^21 + 3*x^19 + 4*x^18 + x^17 + 3*x^15 + 2*x^14 + 4*x^13 + x^12 + 2*x^11 + 3*x^9 + 3*x^8 + 4*x^7 + 2*x^6 + 4*x^5 + 2*x + 2, <22,6,6,3,20>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=24", x^25 + x^10 + 3, x^4 + x^3 + x^2 + 2*x + 2, x^9 + x^8 + x^7 + 3*x^6 + 3*x^5 + 4*x^4 + x + 1, 3*x^24 + 2*x^23 + 2*x^22 + 3*x^21 + 3*x^19 + x^18 + 2*x^17 + 4*x^16 + 4*x^15 + 2*x^14 + 4*x^13 + x^12 + 4*x^11 + 3*x^10 + 4*x^8 + 3*x^7 + x^6 + 3*x^5 + 2*x^4 + 4*x^3 + x^2 + 4*x + 2, <24,6,6,3,4>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=26", x^25 + x^10 + 3, x^3 + 3*x^2 + 1, x^9 + x^8 + x^7 + 3*x^6 + 3*x^5 + 4*x^4 + x + 1, 4*x^24 + x^23 + x^22 + x^21 + 2*x^20 + 2*x^19 + 3*x^18 + x^17 + 2*x^16 + 4*x^15 + x^13 + x^10 + x^9 + x^7 + 4*x^5 + 4*x^4 + 2*x^3 + 4*x + 4, <26,6,6,2,20>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=28", x^25 + x^10 + 3, x^2 + 4*x + 4, x^9 + x^8 + x^7 + 3*x^6 + 3*x^5 + 4*x^4 + x + 1, 3*x^24 + 3*x^23 + 3*x^22 + 4*x^21 + 4*x^20 + 3*x^19 + x^18 + 2*x^16 + x^15 + 3*x^14 + 2*x^13 + 3*x^11 + x^10 + 3*x^9 + 4*x^8 + 2*x^7 + 3*x^6 + 3*x^5 + 4*x^4 + 2*x + 3, <28,6,6,2,4>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=30", x^25 + x^5 + 1, x^2 + x + 4, x^8 + 4*x^7 + 4*x^6 + 3*x^5 + x^4 + 3*x^2 + 2*x + 3, 3*x^24 + 4*x^22 + x^21 + 2*x^20 + x^19 + x^17 + x^16 + 3*x^15 + 2*x^14 + 4*x^13 + 4*x^11 + x^9 + 2*x^8 + x^7 + 3*x^4 + 3*x^3 + x + 2, <30,5,5,2,4>) select 0 else 1);
fails +:= (VerifyGalindo("trino25 K=32", x^25 + x^5 + 1, x + 3, x^8 + 4*x^7 + 4*x^6 + 3*x^5 + x^4 + 3*x^2 + 2*x + 3, x^24 + 3*x^22 + x^21 + 2*x^19 + x^18 + 3*x^16 + 2*x^15 + 3*x^14 + 2*x^13 + 3*x^12 + 3*x^11 + 3*x^10 + 3*x^7 + 3*x^5 + 2*x^4 + 3*x^2 + x + 2, <32,5,5,1,4>) select 0 else 1);

print "\n== Table tab:binomcap, q = 5 ==";
fails +:= (VerifyGalindo("binomcap q=5 m=10", x^10 + 4, x^4 + 3*x^2 + 1, x^4 + 3*x^2 + 1, 2*x^8 + 3*x^7 + 4*x^6 + x^4 + 3*x^3 + 4*x, <4,3,3,2,10>) select 0 else 1);
fails +:= (VerifyGalindo("binomcap q=5 m=25", x^25 + 4, x^12 + 3*x^11 + x^10 + 3*x^7 + 4*x^6 + 3*x^5 + x^2 + 3*x + 1, x^12 + 3*x^11 + x^10 + 3*x^7 + 4*x^6 + 3*x^5 + x^2 + 3*x + 1, 3*x^24 + 3*x^23 + 4*x^21 + 4*x^19 + 4*x^18 + 4*x^17 + x^16 + 4*x^15 + 4*x^14 + 3*x^13 + 2*x^12 + 2*x^11 + 2*x^10 + 4*x^9 + x^8 + 2*x^7 + 3*x^6 + 4*x^5 + 2*x^4 + 4*x^3 + 4*x^2 + x, <2,4,4,3,25>) select 0 else 1);
fails +:= (VerifyGalindo("binomcap q=5 m=50", x^50 + 4, x^20 + 4*x^19 + 3*x^18 + 2*x^17 + x^16 + 4*x^14 + 3*x^13 + 2*x^12 + x^11 + 3*x^10 + x^9 + 2*x^8 + 3*x^7 + 4*x^6 + x^4 + 2*x^3 + 3*x^2 + 4*x + 1, x^21 + 3*x^20 + 4*x^19 + 4*x^18 + 4*x^17 + 4*x^16 + 4*x^15 + 4*x^14 + 4*x^13 + 4*x^12 + 2*x^11 + 3*x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + 2*x + 4, 2*x^49 + 4*x^48 + 3*x^47 + 2*x^46 + 2*x^45 + 2*x^44 + 4*x^43 + 2*x^42 + 3*x^41 + 2*x^40 + x^37 + 2*x^36 + x^35 + 4*x^34 + 2*x^31 + 4*x^29 + 4*x^28 + 2*x^27 + x^26 + 4*x^23 + 3*x^22 + 3*x^21 + x^20 + x^19 + 3*x^18 + 4*x^17 + x^16 + x^15 + 4*x^14 + 4*x^13 + 3*x^11 + 4*x^9 + x^8 + 4*x^7 + x^6 + 2*x^5 + 3*x^4 + x^3 + 4, <18,4,4,4,50>) select 0 else 1);

// ------------------------------ q = 7 ------------------------------
F := GF(7); P<x> := PolynomialRing(F);

print "\n== Table tab:binomcap, q = 7 ==";
fails +:= (VerifyGalindo("binomcap q=7 m=14", x^14 + 6, x^6 + 4*x^4 + 3*x^2 + 6, x^6 + 4*x^4 + 3*x^2 + 6, 2*x^13 + 3*x^12 + 4*x^11 + 3*x^10 + 5*x^9 + 4*x^7 + 2*x^6 + 2*x^5 + x^4 + 3*x^3 + 4*x^2 + 4*x + 3, <4,4,4,4,14>) select 0 else 1);
fails +:= (VerifyGalindo("binomcap q=7 m=21", x^21 + 6, x^7 + x^6 + 6*x^5 + x^4 + 3*x^3 + x^2 + 5*x + 3, x^7 + x^6 + 6*x^5 + x^4 + 3*x^3 + x^2 + 5*x + 3, 6*x^20 + 4*x^19 + 2*x^18 + 5*x^17 + x^16 + x^15 + 5*x^14 + 6*x^13 + 5*x^12 + x^11 + 4*x^10 + 3*x^9 + 3*x^7 + 4*x^6 + 6*x^5 + 6*x^4 + 4*x^3 + 4*x^2 + x + 4, <14,4,4,4,21>) select 0 else 1);
fails +:= (VerifyGalindo("binomcap q=7 m=28", x^28 + 4, x^12 + 5*x^8 + 6*x^4 + 1, x^12 + 5*x^8 + 6*x^4 + 1, 5*x^27 + 4*x^26 + 5*x^25 + x^24 + 4*x^23 + x^21 + 4*x^20 + x^19 + 6*x^17 + 3*x^16 + 4*x^15 + 2*x^14 + 3*x^13 + 2*x^12 + 2*x^11 + 2*x^9 + 5*x^8 + 5*x^7 + x^6 + 2*x^5 + 5*x^4 + 3*x^3 + 3*x^2 + 2*x, <8,4,4,4,168>) select 0 else 1);
fails +:= (VerifyGalindo("binomcap q=7 m=42", x^42 + 3, x^18 + 2*x^12 + 6*x^6 + 6, x^18 + 2*x^12 + 6*x^6 + 6, 5*x^41 + 5*x^40 + 5*x^39 + 3*x^38 + 2*x^37 + 3*x^36 + x^35 + 3*x^34 + 5*x^33 + 4*x^32 + 6*x^31 + 6*x^29 + 4*x^28 + 6*x^27 + 5*x^26 + x^25 + 3*x^24 + x^23 + x^22 + 5*x^21 + x^20 + x^19 + x^16 + 5*x^15 + 3*x^14 + x^13 + 2*x^12 + 5*x^11 + 2*x^10 + 4*x^9 + 5*x^8 + 2*x^7 + x^6 + 6*x^4 + 2*x^3 + x^2 + 4*x + 4, <12,4,4,4,126>) select 0 else 1);
fails +:= (VerifyGalindo("binomcap q=7 m=49", x^49 + 6, x^21 + 4*x^14 + 3*x^7 + 6, x^22 + 6*x^21 + 4*x^15 + 3*x^14 + 3*x^8 + 4*x^7 + 6*x + 1, 6*x^48 + x^47 + 5*x^46 + 6*x^45 + 2*x^44 + 5*x^42 + 4*x^41 + 3*x^40 + 2*x^39 + 4*x^38 + x^37 + 6*x^36 + 6*x^35 + 5*x^33 + 5*x^32 + 4*x^30 + 2*x^29 + 6*x^27 + 4*x^26 + 4*x^25 + 5*x^24 + 6*x^23 + 6*x^22 + x^21 + 6*x^19 + 4*x^18 + x^16 + 2*x^15 + 3*x^14 + 2*x^13 + 6*x^12 + 2*x^11 + 6*x^10 + 4*x^9 + 4*x^8 + 6*x^7 + x^6 + x^5 + 6*x^4 + 3*x^3 + 4*x^2 + 5*x + 5, <12,5,5,4,7>) select 0 else 1);

print "\n== q = 7, m = 49: codes compared with Du-Ma-Liu in the text ==";
fails +:= (VerifyGalindo("dml98 K=30", x^49 + 6, x^10 + 4*x^9 + 3*x^8 + 6*x^7 + 6*x^3 + 3*x^2 + 4*x + 1, x^24 + 4*x^23 + 3*x^22 + 6*x^21 + 4*x^17 + 2*x^16 + 5*x^15 + 3*x^14 + 3*x^10 + 5*x^9 + 2*x^8 + 4*x^7 + 6*x^3 + 3*x^2 + 4*x + 1, 4*x^48 + 6*x^47 + 4*x^46 + 2*x^45 + 5*x^43 + 4*x^42 + 4*x^41 + 5*x^40 + 3*x^39 + 2*x^38 + 4*x^37 + 3*x^36 + x^35 + 5*x^34 + 4*x^33 + x^32 + 5*x^31 + x^30 + 4*x^28 + x^27 + x^26 + 3*x^25 + 5*x^24 + 5*x^23 + 5*x^22 + x^21 + 5*x^20 + 4*x^18 + 6*x^17 + 5*x^16 + 2*x^15 + x^13 + 6*x^11 + x^10 + x^9 + 3*x^8 + 3*x^7 + x^6 + 5*x^5 + 5*x^4 + 3*x^3 + 3*x^2 + 4*x, <30,5,5,3,7>) select 0 else 1);
fails +:= (VerifyGalindo("dml98 K=44", x^49 + 6, x^4 + 3*x^3 + 6*x^2 + 3*x + 1, x^23 + 5*x^22 + x^21 + 4*x^16 + 6*x^15 + 4*x^14 + 3*x^9 + x^8 + 3*x^7 + 6*x^2 + 2*x + 6, 4*x^48 + 2*x^45 + 4*x^44 + 3*x^42 + x^41 + 4*x^40 + 6*x^39 + 2*x^38 + 2*x^37 + 5*x^36 + 5*x^34 + 3*x^33 + 6*x^32 + 5*x^31 + 6*x^30 + 4*x^29 + 3*x^28 + 2*x^27 + 2*x^26 + 5*x^25 + 4*x^24 + x^23 + 2*x^22 + 4*x^21 + 6*x^19 + 2*x^18 + 4*x^17 + 4*x^15 + 5*x^14 + 5*x^12 + 3*x^11 + 4*x^10 + 2*x^9 + x^8 + 5*x^7 + 3*x^6 + x^4 + 3*x^3 + 3*x^2 + 4*x + 5, <44,5,5,2,7>) select 0 else 1);
fails +:= (VerifyGalindo("dml98 K=48", x^49 + 6, x^2 + 5*x + 1, x^23 + 5*x^22 + x^21 + 4*x^16 + 6*x^15 + 4*x^14 + 3*x^9 + x^8 + 3*x^7 + 6*x^2 + 2*x + 6, x^48 + 5*x^47 + 4*x^46 + 4*x^45 + 5*x^44 + 4*x^43 + 5*x^42 + 4*x^40 + 2*x^39 + 2*x^38 + 6*x^37 + 3*x^36 + 3*x^35 + x^34 + 3*x^33 + x^32 + 6*x^31 + 6*x^29 + 4*x^28 + 5*x^27 + 6*x^25 + x^24 + 3*x^23 + 4*x^20 + 4*x^19 + 5*x^18 + 6*x^15 + x^14 + 4*x^13 + 5*x^12 + 6*x^11 + x^10 + 4*x^9 + 2*x^8 + 4*x^6 + 2*x^5 + 5*x^4 + 6*x^3 + 5*x + 5, <48,5,5,1,7>) select 0 else 1);

printf "\n==============================================\n";
printf "Mismatches: %o     (CPU time %o s)\n", fails, Cputime(t0);
printf "==============================================\n";
quit;

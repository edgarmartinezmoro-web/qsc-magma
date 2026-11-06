// =====================================================================
//  verify_q7_one.m
//  Rebuilds one code of Table tab:exp:q7 of the paper and checks its
//  parameters, so that the eight codes can be checked in parallel. Run it from the
//  magma/ directory with the index of the code (1..8):
//
//      magma -b i:=3 verify_q7_one.m
//
//  1 = [[98,64,7]], 2 = [[98,58,7]], 3 = [[98,54,7]], 4 = [[98,50,7]],
//  5 = [[98,42,7]], 6 = [[98,48,5]], 7 = [[98,32,5]], 8 = [[98,14,5]].
//
//  The data f, g, h, k are copied from results/res5_q7_m49_t1.txt,
//  res5_q7_m49_t2.txt, res5_q7_m49_t3.txt (trinomial rows) and
//  results/res2b_q7_m49_QC.txt (binomial rows).
//
//  The checks are the same as in verify_examples.m:
//   * hypotheses g | h and h^2 | f;
//   * annihilator dual C° contained in C (Lemma lem:exp:annGalindo);
//   * Dual(Psi(C) G) contained in Psi(C) (Theorem theo:cssanh);
//   * K = 2 dim C - 2m, dX = d(Psi(C)), dZ = d(Psi(C) G) (exact distances);
//   * synchronization block: best proper divisor g1 of g with ord(g/g1) >= 3,
//     d1 = d(<g1>) and window ord(g/g1) (same rule as the search).
//
//  The exact distances of the five trinomial codes take between about 6 and
//  14 hours each with Magma V2.29-7; the binomial codes take a few seconds.
// =====================================================================

load "qp_qsc_search.m";

// Index of the code to check, given on the command line as i:=n (1..8).
idx := StringToInteger(i);

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

fails := 0;
t0 := Cputime();
F := GF(7); P<x> := PolynomialRing(F);


// source: results/res5_q7_m49_t1.txt, t2.txt, t3.txt
if idx eq 1 then
    fails +:= (VerifyGalindo("trinomial [[98,64,7]]", x^49 + x^7 + 3, x^4 + 6*x^3 + 4*x^2 + 6*x + 6, x^13 + 2*x^12 + 5*x^10 + 6*x^9 + 5*x^8 + 2*x^7 + 5*x^6 + 6*x^5 + 5*x^4 + 2*x^2 + 2*x + 3, 3*x^48 + 5*x^47 + 6*x^46 + 5*x^45 + 4*x^44 + x^43 + 3*x^41 + 5*x^40 + 3*x^39 + 4*x^38 + x^37 + 4*x^36 + 2*x^34 + 5*x^33 + 6*x^31 + x^30 + 5*x^29 + 3*x^28 + x^27 + x^26 + 6*x^25 + 2*x^24 + 3*x^23 + 4*x^22 + 5*x^21 + 3*x^18 + 4*x^17 + x^16 + x^15 + 2*x^14 + 5*x^13 + x^12 + 4*x^11 + 5*x^10 + 3*x^9 + 4*x^8 + 6*x^6 + 5*x^4 + 3*x^3 + 5*x + 6, <64,7,7,2,48>) select 0 else 1);
end if;
if idx eq 2 then
    fails +:= (VerifyGalindo("trinomial [[98,58,7]]", x^49 + x^7 + 1, x^5 + 6*x^4 + 2*x^3 + 4*x^2 + 5*x + 6, x^15 + 4*x^14 + 3*x^13 + 2*x^12 + 4*x^11 + x^9 + 2*x^7 + 5*x^6 + 6*x^5 + 4*x^4 + 2*x^2 + 4*x + 2, 5*x^48 + 3*x^47 + 5*x^46 + 3*x^44 + 3*x^42 + 6*x^41 + 5*x^40 + x^38 + 5*x^37 + x^36 + 2*x^35 + 3*x^34 + 6*x^33 + 4*x^32 + 5*x^31 + 2*x^30 + 5*x^29 + 4*x^28 + 4*x^27 + 5*x^26 + 6*x^25 + 2*x^23 + 2*x^22 + 5*x^21 + 3*x^20 + 5*x^19 + 2*x^18 + 5*x^17 + 4*x^16 + 6*x^15 + x^13 + x^12 + 5*x^11 + x^10 + x^9 + x^8 + 5*x^7 + 4*x^5 + 2*x^4 + 5*x^3 + 4*x^2 + 2*x + 6, <58,7,7,2,48>) select 0 else 1);
end if;
if idx eq 3 then
    fails +:= (VerifyGalindo("trinomial [[98,54,7]] (ex. 1)", x^49 + x^7 + 1, x^6 + 3*x^5 + x^4 + 4*x^3 + x^2 + 3*x + 5, x^16 + x^15 + x^14 + x^12 + 6*x^11 + x^10 + 3*x^8 + 4*x^7 + 2*x^6 + 3*x^5 + x^4 + 4*x^3 + 3*x^2 + 2*x + 5, 5*x^48 + 6*x^47 + 3*x^46 + x^45 + 2*x^44 + 2*x^43 + 4*x^42 + 4*x^41 + x^40 + x^39 + 2*x^38 + 4*x^37 + 6*x^36 + 6*x^35 + 4*x^34 + 6*x^33 + 2*x^32 + x^31 + 3*x^30 + 4*x^29 + 2*x^28 + 4*x^27 + 5*x^26 + 6*x^25 + 4*x^24 + 4*x^23 + 6*x^22 + 6*x^21 + 4*x^20 + 2*x^19 + 6*x^17 + 5*x^16 + 2*x^15 + 6*x^14 + x^12 + 2*x^11 + 3*x^10 + 6*x^9 + 3*x^8 + 5*x^7 + 2*x^6 + 6*x^5 + 3*x^4 + 4*x^3 + 5*x^2 + 6*x + 1, <54,7,7,3,16>) select 0 else 1);
end if;
if idx eq 4 then
    fails +:= (VerifyGalindo("trinomial [[98,50,7]]", x^49 + x^7 + 2, x^10 + 3*x^9 + 5*x^8 + 3*x^7 + 2*x^6 + 5*x^5 + 5*x^4 + 4*x^3 + 2*x^2 + x + 5, x^14 + 3*x^12 + x^11 + 6*x^10 + 4*x^8 + 3*x^6 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3, 4*x^48 + x^47 + 5*x^46 + 2*x^45 + 5*x^44 + 3*x^43 + 3*x^42 + x^41 + 3*x^40 + x^39 + x^38 + 6*x^37 + 2*x^35 + 2*x^34 + x^33 + 4*x^32 + 2*x^31 + 3*x^30 + x^29 + 4*x^28 + 4*x^26 + 4*x^24 + 6*x^22 + x^19 + x^17 + x^16 + 2*x^15 + 5*x^14 + 5*x^13 + 3*x^12 + 6*x^11 + 6*x^10 + 2*x^9 + 3*x^8 + x^7 + 2*x^6 + 6*x^5 + 6*x^4 + 5*x^3 + 2*x^2 + 2*x, <50,7,7,4,48>) select 0 else 1);
end if;
if idx eq 5 then
    fails +:= (VerifyGalindo("trinomial [[98,42,7]] (ex. 2)", x^49 + x^7 + 3, x^12 + 4*x^11 + 3*x^10 + 2*x^9 + x^8 + 3*x^7 + 5*x^6 + x^5 + 3*x^4 + x^3 + x^2 + x + 3, x^16 + 3*x^15 + 3*x^14 + x^12 + 4*x^11 + 2*x^9 + 6*x^8 + 2*x^6 + 2*x^5 + x^4 + x^3 + 5*x^2 + 5*x + 3, 2*x^48 + 2*x^47 + 4*x^46 + 5*x^44 + 3*x^43 + 5*x^42 + 3*x^41 + x^38 + 5*x^37 + 2*x^36 + 6*x^35 + 5*x^34 + 6*x^33 + x^31 + 5*x^30 + 3*x^28 + 3*x^27 + 2*x^26 + 2*x^25 + x^23 + 4*x^22 + 6*x^21 + 4*x^20 + x^19 + x^18 + 5*x^17 + 2*x^16 + 2*x^15 + 2*x^12 + x^10 + 3*x^9 + x^8 + 5*x^7 + 5*x^6 + 6*x^5 + 4*x^4 + 3*x^3 + x^2 + 3*x + 1, <42,7,7,5,48>) select 0 else 1);
end if;

// source: results/res2b_q7_m49_QC.txt
if idx eq 6 then
    fails +:= (VerifyGalindo("binomial [[98,48,5]]", x^49 + 6, x^2 + 5*x + 1, x^23 + 5*x^22 + x^21 + 4*x^16 + 6*x^15 + 4*x^14 + 3*x^9 + x^8 + 3*x^7 + 6*x^2 + 2*x + 6, x^48 + 5*x^47 + 4*x^46 + 4*x^45 + 5*x^44 + 4*x^43 + 5*x^42 + 4*x^40 + 2*x^39 + 2*x^38 + 6*x^37 + 3*x^36 + 3*x^35 + x^34 + 3*x^33 + x^32 + 6*x^31 + 6*x^29 + 4*x^28 + 5*x^27 + 6*x^25 + x^24 + 3*x^23 + 4*x^20 + 4*x^19 + 5*x^18 + 6*x^15 + x^14 + 4*x^13 + 5*x^12 + 6*x^11 + x^10 + 4*x^9 + 2*x^8 + 4*x^6 + 2*x^5 + 5*x^4 + 6*x^3 + 5*x + 5, <48,5,5,1,7>) select 0 else 1);
end if;
if idx eq 7 then
    fails +:= (VerifyGalindo("binomial [[98,32,5]]", x^49 + 6, x^10 + 4*x^9 + 3*x^8 + 6*x^7 + 6*x^3 + 3*x^2 + 4*x + 1, x^23 + 5*x^22 + x^21 + 4*x^16 + 6*x^15 + 4*x^14 + 3*x^9 + x^8 + 3*x^7 + 6*x^2 + 2*x + 6, 2*x^48 + x^45 + 6*x^44 + 2*x^43 + x^42 + 2*x^41 + 4*x^40 + 4*x^39 + 6*x^38 + 5*x^37 + 6*x^36 + 2*x^35 + 5*x^34 + 4*x^32 + 5*x^31 + 2*x^30 + x^29 + 4*x^28 + 3*x^27 + 3*x^26 + 4*x^25 + 4*x^24 + 3*x^23 + 6*x^22 + 4*x^21 + 4*x^20 + 2*x^19 + 3*x^18 + 2*x^17 + 2*x^16 + 2*x^15 + 3*x^14 + x^13 + x^12 + 5*x^11 + 5*x^10 + 5*x^9 + 4*x^8 + 4*x^7 + 2*x^6 + 2*x^4 + 3*x^3 + 4*x^2 + x + 2, <32,5,5,3,7>) select 0 else 1);
end if;
if idx eq 8 then
    fails +:= (VerifyGalindo("binomial [[98,14,5]]", x^49 + 6, x^20 + x^19 + x^18 + x^17 + x^16 + x^15 + x^14 + 5*x^13 + 5*x^12 + 5*x^11 + 5*x^10 + 5*x^9 + 5*x^8 + 5*x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1, x^22 + 6*x^21 + 4*x^15 + 3*x^14 + 3*x^8 + 4*x^7 + 6*x + 1, 2*x^48 + x^47 + 5*x^46 + x^43 + 6*x^42 + 6*x^41 + x^40 + 4*x^39 + 6*x^38 + 6*x^37 + 2*x^36 + x^35 + 6*x^34 + 2*x^33 + 3*x^32 + 6*x^31 + 6*x^30 + 5*x^29 + 5*x^28 + 3*x^27 + 2*x^25 + 2*x^24 + 4*x^23 + 2*x^21 + 5*x^20 + 3*x^18 + 5*x^17 + 3*x^16 + 2*x^15 + 5*x^14 + 3*x^13 + 6*x^11 + x^10 + 6*x^9 + 6*x^8 + 4*x^7 + 2*x^5 + x^4 + x^2 + 5*x + 3, <14,5,5,4,7>) select 0 else 1);
end if;

printf "\n==============================================\n";
printf "Mismatches: %o     (CPU time %o s)\n", fails, Cputime(t0);
printf "==============================================\n";
quit;

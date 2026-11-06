// Lengths m = l p^2 (roots of multiplicity p^2): the structured family is no longer
// limited to d <= 3, and some lambda are not m-th powers (new quasi-constacyclic codes).
load "qp_qsc_search.m";
B := RunGalindo(3, [18, 36], "consta", 30, 3);   PrintBest(B);  SaveGalindo(B, "res3_q3.txt");
B := RunGalindo(5, [50], "consta", 20, 3);       PrintBest(B);  SaveGalindo(B, "res3_q5.txt");
// Repeated-root trinomials at the length of the Du-Ma-Liu code [[50,4,4]]_5
B := RunGalindo(5, [25], "reptri", 20, 3);       PrintBest(B);  SaveGalindo(B, "res3_q5_tri25.txt");
print "END run3";
quit;

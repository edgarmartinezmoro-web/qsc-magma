// Quick test: it should finish in a few seconds with no assertion failures
load "qp_qsc_search.m";
B := RunGalindo(5, [10], "consta", 3, 3);
PrintBest(B);
SaveGalindo(B, "res0_test.txt");
print "OK test";
quit;

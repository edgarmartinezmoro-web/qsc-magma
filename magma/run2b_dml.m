// Reduced version of run2: only the Du-Ma-Liu lengths N = 50 and N = 98,
// with small samples. Each part saves its file as soon as it finishes.
load "qp_qsc_search.m";
B := RunGalindo(5, [25], "consta", 10, 3 : maxpairs := 100);  PrintBest(B);  SaveGalindo(B, "res2b_q5_m25_QC.txt");
B := RunGalindo(5, [25], "reptri", 3, 3 : maxpairs := 30);    PrintBest(B);  SaveGalindo(B, "res2b_q5_m25_tri.txt");
B := RunGalindo(7, [49], "consta", 5, 3 : maxpairs := 60);    PrintBest(B);  SaveGalindo(B, "res2b_q7_m49_QC.txt");
print "END run2b";
quit;

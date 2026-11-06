// Split version of run1 (not included): each call saves its results as soon as it finishes.
load "qp_qsc_search.m";
B := RunGalindo(7, [14], "consta", 20, 3);      PrintBest(B);  SaveGalindo(B, "res1b_q7_m14.txt");
B := RunGalindo(7, [21], "consta", 20, 3);      PrintBest(B);  SaveGalindo(B, "res1b_q7_m21.txt");
B := RunGalindo(3, [6, 12, 18], "consta", 20, 3); PrintBest(B);  SaveGalindo(B, "res1b_q3.txt");
B := RunGalindo(3, [24], "consta", 20, 3);      PrintBest(B);  SaveGalindo(B, "res1b_q3_m24.txt");
B := RunGalindo(7, [28], "consta", 10, 3 : maxpairs := 150);  PrintBest(B);  SaveGalindo(B, "res1b_q7_m28.txt");
B := RunGalindo(7, [42], "consta", 5, 3 : maxpairs := 60);    PrintBest(B);  SaveGalindo(B, "res1b_q7_m42.txt");
print "END run1b";
quit;

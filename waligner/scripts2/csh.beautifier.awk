/^#/{line++;printf("%4d %2d %2d\t",line,n1,n2);print;next;}
/#/{split ($0,aa,"#");$0=aa[1];}
/\/\//{split ($0,aa,"//");$0=aa[1];}
/ then/{if ($NF == "then") {n2++; ok=1 ;}}
/else if/{if (ok<1 && $NF == "then") {n2--; ok=1 ;}}
/endif/{if (ok<1 && $1 == "endif") {n2--; ok=1; }}
/end/{if (ok<1 && $1 == "end") {n1--; ok=1 ; }}
/foreach/{if (ok<1) {n1++; ok=1 ;}}
/ while/{if (ok<1){n1++; ok=1 ;}}

{line++;printf("%4d %2d %2d\t",line,n1,n2);print;ok=0;}

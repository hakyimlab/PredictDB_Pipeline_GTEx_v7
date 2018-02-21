#!/bin/bash

echo "file status" > o_summary.txt
echo "" > e_summary.txt
for p in joblogs/*.out; do
   echo -n "$p " >> o_summary.txt
   cat $p | tail -1 | \
awk '
{
   n = split($i, array, " / ")
   if (n == 0) {
     printf("not_ran\n")
   } else if (array[1] == array[2]) {
     printf("complete\n")
   } else {
     printf("incomplete\n")
   }
}
END { if (!NR) print "not ran" }
' >> o_summary.txt
done
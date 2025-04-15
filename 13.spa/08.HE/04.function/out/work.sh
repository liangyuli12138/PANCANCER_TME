cut -d "," -f 2 TLS_C1.diff.csv|head -n 101 |awk '!/name/ {print $1"\tTLS_C1"}' > fib.diff.list
cut -d "," -f 2 TLS_C2.diff.csv |head -n 101 |awk '!/name/ {print $1"\tTLS_C2"}' >> fib.diff.list
cut -d "," -f 2 TLS_C3.diff.csv |head -n 101 |awk '!/name/ {print $1"\tTLS_C3"}' >> fib.diff.list


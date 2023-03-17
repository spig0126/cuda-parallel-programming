start=77199
end=80000

echo "START add testing $start ~ $end!";

for (( c=$start; c<=$end; c++ ));
do ./add.out $c;
done;
echo;
echo "DONE add testing $start ~ $end!";
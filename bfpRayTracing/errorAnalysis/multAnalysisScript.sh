start=1001
end=5000

echo "START mult testing $start ~ $end!"

for (( c=$start; c<=$end; c++ ));
do ./mult.out $c;
done;
echo;
echo "DONE mult testing $start ~ $end!";
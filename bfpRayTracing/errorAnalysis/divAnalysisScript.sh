start=77199
end=80000

echo "START div testing $start ~ $end!";

for (( c=$start; c<=$end; c++ ));
do ./div.out $c;
done;
echo;
echo "DONE div testing $start ~ $end!";
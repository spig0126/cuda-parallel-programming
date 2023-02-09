# start=11220
start=6001
end=10000

echo "START div testing $start ~ $end!";

for (( c=$start; c<=$end; c++ ));
do ./div.out $c;
done;
echo;
echo "DONE div testing $start ~ $end!";
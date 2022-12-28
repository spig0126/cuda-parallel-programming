#include <fstream>
#include <iostream>
#include <math.h>
#include <limits>

using namespace std;

typedef numeric_limits<float> flt;

int main(){
    string file_path = "all_possible_nums.txt";
    ofstream write_file(file_path.data());

    //set maximum precision
    write_file.precision(flt::max_digits10 + 2);

    unsigned int MAX_NUM = (unsigned int)pow(2, 32) - 1;
    float *f;

    for (unsigned int n = 0; n <= MAX_NUM; n++)
    {
        f = reinterpret_cast<float*>(&n);
        write_file << *f << endl;
    }

    write_file.close();

    return 0;
}
#include <cstdio>
#include <iostream>
#include <fstream>

int main()
{
    std::ifstream fin("input.txt");
    std::ofstream fout("output.txt");
    char st[100];
    int N;
    fin>>N;
    for(int i = 0; i < N; ++i){
	fin.getline(st, sizeof(st));
	if(i % 30 == 0) fout<<st<<std::endl;
    }
    return 0;
}

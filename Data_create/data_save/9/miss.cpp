/*  ATTENTION!
 *  RAND_MAX is set to 32767 in some compiler, which will be too small for large scale data.
 *  (rand() % 2147483647 doesn't make sense)
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>

int main(int argc, char* argv[])
{
    srand(time(0));
    int N, A;
    int num, per, m, flag, missed;
    double x;
    std::ifstream fin("input.txt");
    fin>>N>>A;
    num = N*A;
    std::cout<<N<<" datas, "<<A<<" attributions"<<std::endl;
    std::cout<<"Input the percentage(0-100): ";
    if(argc > 1) {
	per = atoi(argv[1]);
	std::cout<<per<<std::endl;
    }
    else std::cin>>per;
    fin.close();
    char output[50];
    for(int ii = 0; ii < 10; ++ii){
	    sprintf(output, "%d%%-miss-%d.txt", per, ii);
	    fin.open("input.txt");
	    fin>>N>>A;
	    num = N * A;
	    std::ofstream fout(output);
	    m = per*num/100;
	    missed = 0;
	    fout<<N<<std::endl;
	    for(int i = 0; i < N; ++i){
				for(int j = 0; j < A; ++j){
			    int k = (rand()<<15)+rand();
			    fin>>x;
			    if((k%num)<(m-missed)) { fout<<"? "; missed++; }
			    else fout<<x<<' ';
			    num--;
				}
				fin>>flag;
				//num--;
				fout<<flag<<std::endl;
		  }
	    if(m == missed){
				std::cout<<"Done. "<<m<<" data missed"<<std::endl;
	    }
	    else{
				std::cout<<"Error. Only "<<missed<<" data missed (the number should be "<<m<<std::endl;
	    }
	    fin.close();
	    fout.close();
    }
    return 0;
}

#include <iostream>
using namespace std;

int main()
{
    double a, b;
    double x = 0, y = 0;
    int N = 0;
    while(cin>>a>>b){
	if(b < 0){
	    cout<<x/N<<' '<<y/N<<endl;
	    x = y = 0;
	    N = 0;
	}
	else{
	    x += a; y += b;
	    N++;
	}

    }
    return 0;
}

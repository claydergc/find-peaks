#include <iostream>
#include "Util.h"

using namespace std;

int main()
{
	float signal[14] = {0,1,1,1,1,1,1,5,1,1,1,1,7,1};
	vector<float> in(signal, signal + sizeof(signal) / sizeof(float));
	vector<int> out;

	findPeaks(in, out);

	cout<<"Maxima found:"<<endl;

	for(int i=0; i<out.size(); ++i)
		cout<<in[out[i]]<<" ";

	cout<<endl;
}

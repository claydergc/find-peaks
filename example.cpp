#include <iostream>
#include "PeakFinder.h"

int main()
{
	//float inArr[14] = {0,1,1,1,1,1,1,5,1,1,1,1,1,7};
	float inArr[4] =  { 1, 0, 0, 1 };

	std::vector<float> in(inArr, inArr + sizeof(inArr) / sizeof(float));
	std::vector<int> out;

    PeakFinder::findPeaks(in, out, false);

	if(out.size()==0)
	{
		std::cout<<"No peaks"<<std::endl;
		return 0;
	}

	std::cout<<"Maxima found:"<<std::endl;

	for(int i=0; i<out.size(); ++i)
		std::cout<<in[out[i]]<<" ";

	std::cout<<std::endl;

	return 0;
}

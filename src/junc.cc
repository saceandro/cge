#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

typedef std::pair<int,int> READ;
typedef std::vector<READ> READS;

int main()
{
  ifstream infilestream("reads.txt", std::ifstream::in);


  READS res;
  while(true)
    {
      READ re;
      infilestream >> re.first >> re.second;
      if (infilestream.eof()) break;
      res.push_back(re);
    }
  infilestream.close();

  for (READS::iterator it=res.begin(); it!=res.end(); ++it)
    cout << it->first << "\t" << it->second << endl;
  
  return 0;
}

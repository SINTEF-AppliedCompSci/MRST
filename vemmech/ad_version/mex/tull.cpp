#include <iostream>


using namespace std;
int main()
{
  double a = 0;
  for (size_t i = 0; i != 1000000000; ++i) 
    a = a+1;
  cout << a<< endl;
  return 0;
}

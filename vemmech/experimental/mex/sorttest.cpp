#include <algorithm>
#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <future>

using namespace std;

int main()
{
  int N = 2e8;

  // fill vector with random stuff
   std::uniform_real_distribution<double> unif(0,1);
   std::default_random_engine re;

   auto t0 = chrono::high_resolution_clock::now();
   vector<double> vec(N);
   std::generate(vec.begin(), vec.end(), [&]() {return unif(re);});
   auto t1 = chrono::high_resolution_clock::now();

   chrono::duration<double> diff = t1-t0;
   
   cout << "Time to generate indata: " << diff.count() << endl;

   // sorting the whole shebang
   auto cpy = vec;
   t0 = chrono::high_resolution_clock::now();
   sort(cpy.begin(), cpy.end());
   t1 = chrono::high_resolution_clock::now();
   diff = t1-t0;
   cout << "Time to sort in one pass: " << diff.count() << endl;
        
   t0 = chrono::high_resolution_clock::now();
   sort(cpy.begin(), cpy.end());
   t1 = chrono::high_resolution_clock::now();
   diff = t1-t0;
   cout << "Time to re-sort in one pass: " << diff.count() << endl;

   // for parallel thread at it
   int num_threads = 8;
   size_t chunklen = N/num_threads;
   // auto sort_job = [](vector<double>::iterator begin,
   //                    vector<double>::iterator end) { std::sort(begin, end); };

   t0 = chrono::high_resolution_clock::now();
   vector<future<void> > futures(num_threads);
   //cpy = vec; // reset to "unsorted state"
   for (int i = 0; i != num_threads; ++i) {
     cout << "initializing future " << i << endl;
     if (i != num_threads - 1)
       futures[i] = async([&] { sort(cpy.begin() + i * chunklen,
                                    cpy.begin() + (i+1) * chunklen);} );
     else
       futures[i] = async([&] { sort(cpy.begin() + i * chunklen, cpy.end()); });
     
   }
   // wait for all thread to finish
   for (int i = 0; i != num_threads; ++i) {
     futures[i].get();
   }
   futures.clear();
   t1 = chrono::high_resolution_clock::now();
   diff = t1 - t0;
   cout << "Threads finished in: " << diff.count() << endl;

   // final merge of thread results
   t0 = chrono::high_resolution_clock::now();
   for (int i = 1; i != num_threads; ++i) {

     auto cstart = cpy.begin() + i * chunklen;
     auto cend = (i==num_threads-1) ? cpy.end() : cstart + chunklen;
       
     std::inplace_merge(cpy.begin(), cstart, cend);
   }
   t1 = chrono::high_resolution_clock::now();
   diff = t1 - t0;
   cout << "Merge finished in: " << diff.count() << endl;
     
   
}


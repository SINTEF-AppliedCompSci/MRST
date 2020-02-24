#include <algorithm>
#include <vector>
#include <iostream>
#include <chrono>
#include <random>
#include <future>

using namespace std;

int main()
{
  int N = 4e8;

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
   int nte = 5;
   int num_threads = pow(2, nte);
   size_t chunklen = N/num_threads;
   // auto sort_job = [](vector<double>::iterator begin,
   //                    vector<double>::iterator end) { std::sort(begin, end); };

   t0 = chrono::high_resolution_clock::now();
   vector<future<void> > futures(num_threads);
   vector<vector<double>> chunks(num_threads);
   cpy = vec; // reset to "unsorted state"
   for (int i = 0; i != num_threads; ++i) {
     cout << "initializing future " << i << endl;
     if (i != num_threads - 1)
       chunks[i].insert(chunks[i].end(),
                        cpy.begin() + i * chunklen,
                        cpy.begin() + (i+1) * chunklen);
     else
       chunks[i].insert(chunks[i].end(),
                        cpy.begin() + i * chunklen, cpy.end());

     futures[i] = async([&chunks, i] { sort(chunks[i].begin(), chunks[i].end());});
   }
     
   cout << "All futures created." << endl;
   // wait for all thread to finish
   for (int i = 0; i != num_threads; ++i) {
     futures[i].get();
   }
   cout << "All futures finished." << endl;
   //futures.clear();
   t1 = chrono::high_resolution_clock::now();
   diff = t1 - t0;
   cout << "Threads finished in: " << diff.count() << endl;

   // final merge of thread results
   t0 = chrono::high_resolution_clock::now();

   bool merge_parallel = true;
   //bool merge_parallel = false;

   if (merge_parallel) {
   
     vector<future<void>> futures2(num_threads/2);
     for (int n = 0; n != nte; ++n) {
       cout << "Pass : " << n << endl;
       for (int i = 0; i != chunks.size(); i+=2) {
         futures2[i/2] = async([&chunks, i] {
             size_t mid = chunks[i].size();
             chunks[i].insert(chunks[i].end(),
                              chunks[i+1].begin(), chunks[i+1].end());
             std::inplace_merge(chunks[i].begin(),
                                chunks[i].begin() + mid,
                                chunks[i].end());
           });
       }
       // wait for all to finish
       vector<vector<double>> tmp(chunks.size()/2);
       for (int i = 0; i != chunks.size() / 2; ++i) {
         futures2[i].get();
         tmp[i].swap(chunks[2*i]);
       }
       
       chunks.swap(tmp);
       cout << "Number of chunks reduced to: " << chunks.size() << endl;
       
     }
   } else {
     vector<double> result(chunks[0]);
     for (int i = 1; i != chunks.size(); ++i)
       result.insert(result.end(), chunks[i].begin(), chunks[i].end());
     
     for (int i = 0; i != num_threads - 1; ++i) {
       const auto cstart = result.begin() + i * chunklen;
       const auto cend = (i == num_threads - 1) ? result.end() : cstart + chunklen;
       std::inplace_merge(result.begin(), cstart, cend);
     }
   }

   t1 = chrono::high_resolution_clock::now();
   diff = t1 - t0;
   cout << "Merge finished in: " << diff.count() << endl;
     
   
}


#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <array>
#include <unordered_map>
#include <numeric>
#include <functional>

using namespace std;

struct pairhash {
  size_t operator()(const pair<size_t, size_t>& arg) const {
    return arg.first;
    //return hash<unsigned int>()(arg.first) ^ hash<unsigned int>()(arg.second);
  }
};

struct Entry {
  size_t i1;
  size_t i2;
};

inline bool operator<(const Entry& e1, const Entry e2) {
  return (e1.i1 == e2.i1) ? (e1.i2 < e2.i2) : (e1.i1 < e2.i1);
}
inline bool operator==(const Entry& e1, const Entry e2) {
  return (e1.i1 == e2.i1) && (e1.i2 == e2.i2);
}


struct entryhash {
  size_t operator()(const Entry& e) const {
    return e.i1;
    //return arg.first;
    //return hash<unsigned int>()(arg.first) ^ hash<unsigned int>()(arg.second);
  }
};



// void accumarray(vector<size_t>& icol,
//                 vector<size_t>& irow,
//                 vector<double>& value)
// {
//   const auto N = icol.size();

  
  
//   vector<Entry> tmp;
//   tmp.reserve(icol.size());
//   transform(icol.begin(), icol.end(), irow.begin(), back_inserter(tmp),
//             [](size_t i, size_t j) {return Entry {i, j};});
//   sort(tmp.begin(), tmp.end());

  
  
// }

void accumarray(vector<size_t>& icol,
                vector<size_t>& irow,
                vector<double>& value)
{
  const auto N = icol.size();

  // unordered_map<pair<size_t, size_t>,
  //               double, pairhash> umap;
  unordered_map<Entry, double, entryhash> umap;

  
  vector<size_t>::const_iterator icp = icol.begin();
  vector<size_t>::const_iterator irp = irow.begin();
  vector<double>::const_iterator vp = value.begin();

  while (icp != icol.end()) {
    //const auto rc = make_pair(*irp++, *icp++);
    const auto rc = Entry {*irp++, *icp++};
    const auto v = *vp++;

    auto it = umap.find(rc);
    if (it == umap.end()) {
      umap.insert({rc, v});
    } else {
      it->second += v;
    }
  }

  cout << "Converting" << endl;
  // converting back
  size_t NNZ = umap.size();
  icol.resize(NNZ);
  irow.resize(NNZ);
  value.resize(NNZ);
  size_t counter = 0;

  for (auto it : umap) {
    irow[counter] = it.first.i1;
    icol[counter] = it.first.i2;
    value[counter++] = it.second;
  }
}

int main() {

  cout << "Reading input data" << endl;
  vector<size_t> irow;
  ifstream is("ic");
  
  const auto start_read = chrono::high_resolution_clock::now();
  while (!is.eof()) {
    double x;
    is >> x;
    irow.push_back(x);
  }
  const auto end_read = chrono::high_resolution_clock::now();
  cout << "read time: " << chrono::duration_cast<chrono::milliseconds>(end_read - start_read).count() << endl;

  vector<size_t> icol(irow.size());
  vector<double> value(irow.size(), 1.0);
  iota(icol.begin(), icol.end(), 0);
  
  cout << "Now accumulating..." << endl;
  const auto start_accum = chrono::high_resolution_clock::now();
  accumarray(irow, icol, value);
  const auto end_accum = chrono::high_resolution_clock::now();  
  cout << "Finished" << endl;
  cout << "Accumulation time: " << chrono::duration_cast<chrono::milliseconds>(end_accum - start_accum).count() << endl;

}





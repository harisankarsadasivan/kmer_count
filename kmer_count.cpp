#include <bits/stdc++.h>
#include <cstdlib>
#include <future>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

using namespace std;

//modify this for increased asynchronous parallelism from thread pools in histogram computation
#define NUMBER_OF_THREADS 3

//helper function for data distribution in multi-threading
//helps split large input data by returning vector of equally spaced anchors indices
vector<int> get_indices_for_worker_threads(size_t size, int no_of_threads) {
  vector<int> anchor_idx{};
  int w = size / no_of_threads;
  int idx = 0;
  while (idx < size) {
    anchor_idx.push_back(idx);
    idx += w;
  }
  if (idx != size - 1)
    anchor_idx.push_back(size - 1);

  return anchor_idx;
}

//returns map of kmer and it's count.
map<string, int> get_histogram(string &data, size_t size, int no_of_threads,
                               int KMER_SIZE) {
  //store histogram of kmers
  map<string, int> histogram{};

  auto indexes = get_indices_for_worker_threads(size, no_of_threads);
  vector<future<map<string, int>>> futures{};

  // loop to start threads
  for (auto i = 0; i < indexes.size() - 1; i++) {
    int head = indexes[i];
    int tail = indexes[i + 1];

    //multi-threaded asynchronous function call
    future<map<string, int>> future =
        async(launch::async, [data, head, tail, size, KMER_SIZE]() {
          string kmer_key = "";
          map<string, int> sub_histogram{};

          for (auto i = head; (i < tail) && (i < (size - KMER_SIZE + 1)); ++i) {
            for (auto j = 0; j < KMER_SIZE; j++) {

              kmer_key += char(data[i + j]);
            }

            sub_histogram[kmer_key] += 1u;
            kmer_key = "";
          }
          //return sub-histogram per thread
          return sub_histogram;
        });

    futures.push_back(std::move(future));
  }

  // loop to aggregate threads results
  for (auto i = 0; i < no_of_threads; ++i) {
    auto result = futures[i].get();

    for (auto item : result) {

      histogram[item.first] += item.second;
    }

  }
  return histogram;
}

int main(int argc, char *argv[]) {
  string piped_input;

  // read from cin input buffer
  cin >> piped_input;

  // set n,k(kmer_size) and other variables: 
  //prev remembers previous value of kmer while printing out
  //n_count is a counter
  int n = atoi(argv[2]), prev = 0, n_count = 1, KMER_SIZE = atoi(argv[1]);

  // find histogram-- key is kmer and value is frequency of kmer
  // we do this for compute using asynchsronous parallelization
  auto histogram =
      get_histogram(piped_input, piped_input.length(), NUMBER_OF_THREADS, KMER_SIZE);

  // create a multimap where key is frequency and value is a kmer
  multimap<int, string, greater<int>> kmer_table;
  for (auto item : histogram)
    kmer_table.insert(pair<int, string>(item.second, item.first));

  // print output: iterate over all kmers in multimap and print kmers that have same frequency
  for (auto i = kmer_table.begin(); i != kmer_table.end(); i++) {

    if (i->first != prev) {
      //exit condition
      if (n_count == n + 1) {
        cout << "\n";
        break;
      }
      cout << "\n";
      cout << i->first << "x\t" << i->second;
      n_count++;
      prev = i->first;
    } else
      cout << "\t" << i->second;
  }
  cout << "\n";

  return 0;
}
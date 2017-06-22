#ifndef RaggedArray_inc
#define RaggedArray_inc

#include<vector>

using namespace std;

typedef vector<int>::iterator v_iterator;

struct RaggedArray {
  vector<int> data; //originally x
  vector<int> pos; //pos
  int p;
  RaggedArray() {p=0;}; //constructor initialising p to 0
  // functions
  void new_row() { pos.push_back(data.size()); }; //inserts size of data into pos
  void push_back(int new_data_pt) { 
    data.push_back(new_data_pt);
    p=max(p, new_data_pt);
  }; //insert integer to data
  void push_back(v_iterator new_data_pts_begin, v_iterator new_data_pts_end) {
    new_row();
    data.insert(data.end(), new_data_pts_begin, new_data_pts_end);
    p=max(*max_element(new_data_pts_begin, new_data_pts_end), p);
  }; //insert size of data into pos, then insert to end of data a sequence of integers determined by the pointers new_data_pts_begin and new_data_pts_end
  v_iterator begin(int i) {
    return data.begin() + pos[i];
  }; //pointer corresponding to beginning of ith observation in data (starting from observation 0)
  v_iterator end(unsigned int i) {
    if (i < pos.size() - 1) {
      return data.begin() + pos[i+1];
    } else {
      return data.end();
    }
  }; //pointer corresponding to end of ith observation in data (starting from observation 0)
  int nrow() {
    return pos.size();
  }; //number of observations
  int ncol() {
    return p+1;
  }; //number of variables
};

#endif

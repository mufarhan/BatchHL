#ifndef HGHWAY_LABELING_H_
#define HGHWAY_LABELING_H_

#include <sys/time.h>
#include <iostream>
#include <thread>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <fstream>

#include "two_layer_queue.h"

//
// NOTE: Currently only unweighted and undirected graphs are supported.
//

class HighwayLabelling {
public:
  // Constructs an index from a graph, given as a list of edges.
  HighwayLabelling(std::string filename, int k);
  HighwayLabelling();
  ~HighwayLabelling();

  void ConstructHighwayLabelling(int i, int topk[]);
  void BuildIndex(int topk[]);

  void UpdateLabelling(std::string filename, int m, int p);
  void BHL_Plus(std::vector<std::pair<std::string, std::pair<int, int> > > updates);
  void BHL_Plus_Parallel(int i, std::vector<std::pair<std::string, std::pair<int, int> > > updates);
  void BHL(std::vector<std::pair<std::string, std::pair<int, int> > > updates);
  void BHL_Parallel(int i, std::vector<std::pair<std::string, std::pair<int, int> > > updates);
  bool prunable(int i, int u, uint8_t *temp, uint8_t *A);

  void deallocate();

  void SelectLandmarks_HD(int topk[]);
  int LabellingSize();

  uint8_t query(int r, int v);
  uint8_t min(uint8_t a, uint8_t b);
  void RemoveLandmarks(int topk[]);
  uint8_t HC_UB_naive(int s, int t);
  void QueryDistance(std::string pairs, std::string output);

  void storeLabelling(std::string filename);
  void loadLabelling_Full(std::string filename, int topk[]);
  void loadLabelling_Pruned(std::string filename);

  void PrintLabelingStatistics(std::string s);
  void PrintQueryStatistics();

private:
  int V;  // total number of vertices
  long E; // total number of edges
  int K; // total number of landmarks

  uint8_t **distances, **distances_1;
  uint8_t **highway, **highway_1;
  uint8_t **vertices;
  uint8_t *C;
  std::vector<std::vector<int> > adj;
  std::map<int, uint8_t> landmarks;

  double GetCurrentTimeSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
  }

  long GetCurrentTimeMicroSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t) 1e6 + tv.tv_usec;
  }

  long GetCurrentTimeMilliSec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * 1000LL + tv.tv_usec / 1000;
  }

  // Statistics
  double time_, time_querying_sec_;
  long time_querying_microsec_, time_querying_millisec_;
};

HighwayLabelling::HighwayLabelling() { }

HighwayLabelling::~HighwayLabelling() { }

void HighwayLabelling::deallocate() {

  for(int i = 0; i < V; i++) {
    delete [] distances[i];
    delete [] distances_1[i];
  }
  delete [] distances;
  delete [] distances_1;

  for(int i = 0; i < K; i++) {
    delete [] highway[i];
    delete [] highway_1[i];
  }
  delete [] highway;
  delete [] highway_1;
}

HighwayLabelling::HighwayLabelling(std::string filename, int k) {
  K = k; V = 0; E = 0;

  std::ifstream ifs(filename);
  if (ifs.is_open()) {
    ifs >> V >> E;

    adj.reserve(V);

    int v, w, deg;
    while (ifs >> v >> deg) {
      adj[v].reserve(deg);
      for (int i = 0; i < deg; i++) {
        ifs >> w;
        adj[v].push_back(w);
      }
    }
    ifs.close();

    std::cout << "V : " << V << " E : " << E << std::endl << std::endl;
  } else
      std::cout << "Unable to open file" << std::endl;
}

int HighwayLabelling::LabellingSize() {
  long size = 0;
  for (int i = 0; i < V; i++) {
    for (int j = 0; j < K; j++) {
      if(distances[i][j] != 111)
        size++;
    }
  }
  
  return (V + 2 *size) / (1024 * 1024);
}

void HighwayLabelling::ConstructHighwayLabelling(int i, int topk[]) {

  uint8_t *P = new uint8_t[V];
  for(int j = 0; j < V; j++)
    P[j] = 111;

  std::queue<int> que[2];
  que[0].push(topk[i]); que[0].push(-1);
  distances[topk[i]][i] = 0; distances_1[topk[i]][i] = 0; P[topk[i]] = 0; int use = 0;
  while (!que[0].empty()) {
    int u = que[use].front();
    que[use].pop();

    if(u == -1) {
      use = 1 - use;
      que[use].push(-1);
      continue;
    }

    for (int w : adj[u]) {
      if (P[w] == 111) {
        P[w] = P[u] + 1;
        if(use == 1 || landmarks.count(w) > 0)
          que[1].push(w);
        else {
          que[0].push(w);
          distances[w][i] = P[w];
          distances_1[w][i] = P[w];
        }
      }
    }
  }

  for(int j = 0; j < K; j++) {
    if(P[topk[j]] != 111) {
      highway[i][j] = P[topk[j]];
      highway[j][i] = P[topk[j]];

      highway_1[i][j] = P[topk[j]];
      highway_1[j][i] = P[topk[j]];
    }
  }

  delete [] P;
}

void HighwayLabelling::BuildIndex(int topk[]) {

  for(int i = 0; i < K; i++)
    landmarks[topk[i]] = i;

  // Initialization
  distances = new uint8_t*[V];
  distances_1 = new uint8_t*[V];
  for(int i = 0; i < V; i++) {
    distances[i] = new uint8_t[K];
    distances_1[i] = new uint8_t[K];
    for(int j = 0; j < K; j++) {
      distances[i][j] = 111;
      distances_1[i][j] = 111;
    }
  }

  highway = new uint8_t*[K];
  highway_1 = new uint8_t*[K];
  for(int i = 0; i < K; i++) {
    highway[i] = new uint8_t[K];
    highway_1[i] = new uint8_t[K];
  }

  // Start computing Highway Labelling (HL)
  time_ = -GetCurrentTimeSec();
  for (int i = 0; i < K; i++)
    ConstructHighwayLabelling(i, topk);
  time_ += GetCurrentTimeSec();

  std::cout << "Construction Time (sec.): " << time_ << " Labelling Size: " << LabellingSize() << " MB" << std::endl;
}

void HighwayLabelling::UpdateLabelling(std::string filename, int m, int p) {
  std::ifstream ifs(filename); int a, b; std::string op;
  std::vector<std::pair<std::string, std::pair<int, int> > > updates;

  while(ifs >> op >> a >> b) {
    if(op == "EI") {
      adj[a].push_back(b);
      adj[b].push_back(a);
    } else if(op == "ED") {
      adj[a].erase(std::remove(adj[a].begin(), adj[a].end(), b), adj[a].end());
      adj[b].erase(std::remove(adj[b].begin(), adj[b].end(), a), adj[b].end());
    }
    updates.push_back(std::make_pair(op, std::make_pair(a, b)));
  }
  ifs.close();

  time_ = -GetCurrentTimeSec();
  if(p == 0) {
    if(m == 0)
      BHL_Plus(updates);
    else if (m == 1)
      BHL(updates);
  } else if(p == 1) {
    std::thread t[K];
    if(m == 0) {
      for (a = 0; a < K; a++)
        t[a] = std::thread(&HighwayLabelling::BHL_Plus_Parallel, this, a, updates);

      for (a = 0; a < K; a++)
        t[a].join();
    } else if(m == 1) {
      for (a = 0; a < K; a++)
        t[a] = std::thread(&HighwayLabelling::BHL_Parallel, this, a, updates);

      for (a = 0; a < K; a++)
        t[a].join();
    }
  }
  time_ += GetCurrentTimeSec();
  std::cout << "Batch Update Time (sec.): " << time_ << " Updated Labelling Size: " << LabellingSize() << " MB" << std::endl;

  for(a = 0; a < V; a++) {
    for(b = 0; b < K; b++)
      distances_1[a][b] = distances[a][b];
  }

  for(a = 0; a < K; a++) {
    for(b = 0; b < K; b++) {
      highway_1[a][b] = highway[a][b];
    }
  }
}

/** BHL^+ **/
void HighwayLabelling::BHL_Plus(std::vector<std::pair<std::string, std::pair<int, int> > > updates) {
  uint8_t *A = new uint8_t[V];
  uint8_t *temp = new uint8_t[V];
  int *AFF_VERTS = new int[V]; int c; uint8_t br, beta;

  for(int i = 0; i < V; i++)
    A[i] = 111;

  for(int i = 0; i < K; i++) {

    // computing affected vertices
    std::priority_queue<std::pair<uint8_t, int>, std::vector<std::pair<uint8_t, int> >, std::greater<std::pair<uint8_t, int> > > que;
    std::queue<std::pair<uint8_t, int> > que1; c = 0;
    for(std::pair<std::string, std::pair<int, int> > iter : updates) {

      temp[iter.second.first] = query(i, iter.second.first);
      temp[iter.second.second] = query(i, iter.second.second);

      if(temp[iter.second.first] > temp[iter.second.second]) {
        if(iter.first == "EI") {
          br = 4*(temp[iter.second.second] + 1) + (distances_1[iter.second.second][i] == 111 || landmarks.count(iter.second.first) > 0?0:2) + 1; beta = 4*temp[iter.second.first] + (distances_1[iter.second.first][i] == 111?0:2);
          if(br <= beta)
            que.push(std::make_pair(br, iter.second.first));
        } else if(iter.first == "ED") {
          br = 4*(temp[iter.second.second] + 1) + (distances_1[iter.second.second][i] == 111 || landmarks.count(iter.second.first) > 0?0:2); beta = 4*temp[iter.second.first] + (distances_1[iter.second.first][i] == 111?0:2);
          if(br <= beta)
            que.push(std::make_pair(br, iter.second.first));
        }
      } else if(temp[iter.second.first] < temp[iter.second.second]) {
        if(iter.first == "EI") {
          br = 4*(temp[iter.second.first] + 1) + (distances_1[iter.second.first][i] == 111 || landmarks.count(iter.second.second) > 0?0:2) + 1; beta = 4*temp[iter.second.second] + (distances_1[iter.second.second][i] == 111?0:2);
          if(br <= beta)
            que.push(std::make_pair(br, iter.second.second));
        } else if(iter.first == "ED") {
          br = 4*(temp[iter.second.first] + 1) + (distances_1[iter.second.first][i] == 111 || landmarks.count(iter.second.second) > 0?0:2); beta = 4*temp[iter.second.second] + (distances_1[iter.second.second][i] == 111?0:2);
          if(br <= beta)
            que.push(std::make_pair(br, iter.second.second));
        }
      }
    }

    while (!que1.empty() || !que.empty()) {

      std::pair<uint8_t, int> p;
      if(!que1.empty() && !que.empty()) {
        if(que.top().first < que1.front().first) {
          p = que.top(); que.pop();
        } else {
          p = que1.front(); que1.pop();
        }
      } else if(!que1.empty()) {
        p = que1.front(); que1.pop();
      } else {
        p = que.top(); que.pop();
      }

      if(A[p.second] == 111) {

        for (int w : adj[p.second]) {
          temp[w] = query(i, w);

          br = p.first + 4 | ((p.first + 4 & (1 << 1)) == 0 || landmarks.count(w) > 0?0:2); beta = 4*temp[w] + (distances_1[w][i] == 111?0:2);
          if(br <= beta)
            que1.push(std::make_pair(br, w));
        }

        A[p.second] = p.first;
        AFF_VERTS[c] = p.second; c++;
      }
    }

    // computing boundary vertices
    std::vector<std::pair<uint8_t, int> > V_aff;
    for(int j = 0; j < c; j++) {
      temp[AFF_VERTS[j]] = 111;
      for (int w : adj[AFF_VERTS[j]]) {
        if(A[w] == 111)
          temp[AFF_VERTS[j]] = min(temp[AFF_VERTS[j]], temp[w] + 1);
      }

      if(temp[AFF_VERTS[j]] != 111)
        V_aff.push_back(std::make_pair(temp[AFF_VERTS[j]], AFF_VERTS[j]));
      distances[AFF_VERTS[j]][i] = 111;
    }

    // updating the labelling
    if(V_aff.size() > 0) {
      std::sort(V_aff.begin(), V_aff.end());

      std::queue<int> quee[2]; int use = 0;
      uint8_t d = V_aff[0].first; int x = 0;
      while(x < V_aff.size() && V_aff[x].first == d) {
        if(prunable(i, V_aff[x].second, temp, A)){
          quee[1].push(V_aff[x].second);
        } else {
          distances[V_aff[x].second][i] = temp[V_aff[x].second];
          quee[0].push(V_aff[x].second);
        }
        A[V_aff[x].second] = 111;
        x++;
      }

      quee[use].push(-1);
      while (!quee[0].empty() || x < V_aff.size()) {
        int u = quee[use].front();
        quee[use].pop();

        if(u == -1) {
          if(use == 0) { d++;
            while(x < V_aff.size() && V_aff[x].first == d) {
              if(A[V_aff[x].second] != 111) {
                if(prunable(i, V_aff[x].second, temp, A)) {
                  quee[1].push(V_aff[x].second);
                } else {
                  distances[V_aff[x].second][i] = temp[V_aff[x].second];
                  quee[0].push(V_aff[x].second);
                }
                A[V_aff[x].second] = 111;
              }
              x++;
            }
          }
          use = 1 - use;
          quee[use].push(-1);
          continue;
        }

        for (int w : adj[u]) {
          if(A[w] != 111) {
            temp[w] = temp[u] + 1;
            if(use == 1 || prunable(i, w, temp, A)) {
              quee[1].push(w);
            } else {
              distances[w][i] = temp[w];
              quee[0].push(w);
            }
            A[w] = 111;
          }
        }
      }
    }

    if(i == K - 1)
      continue;

    for(int j = 0; j < c; j++) {
      if(A[AFF_VERTS[j]] != 111)
        A[AFF_VERTS[j]] = 111;
    }
  }

  delete [] A;
  delete [] temp;
  delete [] AFF_VERTS;
}

/** BHL^p **/
void HighwayLabelling::BHL_Plus_Parallel(int i, std::vector<std::pair<std::string, std::pair<int, int> > > updates) {
  uint8_t *A = new uint8_t[V];
  uint8_t *temp = new uint8_t[V];
  int *AFF_VERTS = new int[V]; int c; uint8_t br, beta;

  for(int j = 0; j < V; j++)
    A[j] = 111;

  // computing affected vertices
  std::priority_queue<std::pair<uint8_t, int>, std::vector<std::pair<uint8_t, int> >, std::greater<std::pair<uint8_t, int> > > que;
  std::queue<std::pair<uint8_t, int> > que1; c = 0;
  for(std::pair<std::string, std::pair<int, int> > iter : updates) {

    temp[iter.second.first] = query(i, iter.second.first);
    temp[iter.second.second] = query(i, iter.second.second);

    if(temp[iter.second.first] > temp[iter.second.second]) {
      if(iter.first == "EI") {
        br = 4*(temp[iter.second.second] + 1) + 2 + (distances_1[iter.second.second][i] == 111 || landmarks.count(iter.second.first) > 0?0:1); beta = 4*temp[iter.second.first] + (distances_1[iter.second.first][i] == 111?0:2);
        if(br <= beta)
          que.push(std::make_pair(br, iter.second.first));
      } else if(iter.first == "ED") {
        br = 4*(temp[iter.second.second] + 1) + (distances_1[iter.second.second][i] == 111 || landmarks.count(iter.second.first) > 0?0:1); beta = 4*temp[iter.second.first] + (distances_1[iter.second.first][i] == 111?0:2);
        if(br <= beta)
          que.push(std::make_pair(br, iter.second.first));
      }
    } else if(temp[iter.second.first] < temp[iter.second.second]) {
      if(iter.first == "EI") {
        br = 4*(temp[iter.second.first] + 1) + 2 + (distances_1[iter.second.first][i] == 111 || landmarks.count(iter.second.second) > 0?0:1); beta = 4*temp[iter.second.second] + (distances_1[iter.second.second][i] == 111?0:2);
        if(br <= beta)
          que.push(std::make_pair(br, iter.second.second));
      } else if(iter.first == "ED") {
        br = 4*(temp[iter.second.first] + 1) + (distances_1[iter.second.first][i] == 111 || landmarks.count(iter.second.second) > 0?0:1); beta = 4*temp[iter.second.second] + (distances_1[iter.second.second][i] == 111?0:2);
        if(br <= beta)
          que.push(std::make_pair(br, iter.second.second));
      }
    }
  }

  while (!que1.empty() || !que.empty()) {

    std::pair<uint8_t, int> p;
    if(!que1.empty() && !que.empty()) {
      if(que.top().first < que1.front().first) {
        p = que.top(); que.pop();
      } else {
        p = que1.front(); que1.pop();
      }
    } else if(!que1.empty()) {
      p = que1.front(); que1.pop();
    } else {
      p = que.top(); que.pop();
    }

    if(A[p.second] == 111) {

      for (int w : adj[p.second]) {
        temp[w] = query(i, w);

        br = p.first + 4 | ((p.first + 4 & 1) == 0 || landmarks.count(w) > 0?0:1); beta = 4*temp[w] + (distances_1[w][i] == 111?0:2);
        if(br <= beta)
          que1.push(std::make_pair(br, w));
      }

      A[p.second] = p.first;
      AFF_VERTS[c] = p.second; c++;
    }
  }

  // computing boundary vertices
  std::vector<std::pair<uint8_t, int> > V_aff;
  for(int j = 0; j < c; j++) {
    temp[AFF_VERTS[j]] = 111;
    for (int w : adj[AFF_VERTS[j]]) {
      if(A[w] == 111)
        temp[AFF_VERTS[j]] = min(temp[AFF_VERTS[j]], temp[w] + 1);
    }

    if(temp[AFF_VERTS[j]] != 111)
      V_aff.push_back(std::make_pair(temp[AFF_VERTS[j]], AFF_VERTS[j]));
    distances[AFF_VERTS[j]][i] = 111;
  }

  // updating the labelling
  if(V_aff.size() > 0) {
    std::sort(V_aff.begin(), V_aff.end());

    std::queue<int> quee[2]; int use = 0;
    uint8_t d = V_aff[0].first; int x = 0;
    while(x < V_aff.size() && V_aff[x].first == d) {
      if(prunable(i, V_aff[x].second, temp, A)){
        quee[1].push(V_aff[x].second);
      } else {
        distances[V_aff[x].second][i] = temp[V_aff[x].second];
        quee[0].push(V_aff[x].second);
      }
      A[V_aff[x].second] = 111;
      x++;
    }

    quee[use].push(-1);
    while (!quee[0].empty() || x < V_aff.size()) {
      int u = quee[use].front();
      quee[use].pop();

      if(u == -1) {
        if(use == 0) { d++;
          while(x < V_aff.size() && V_aff[x].first == d) {
            if(A[V_aff[x].second] != 111) {
              if(prunable(i, V_aff[x].second, temp, A)) {
                quee[1].push(V_aff[x].second);
              } else {
                distances[V_aff[x].second][i] = temp[V_aff[x].second];
                quee[0].push(V_aff[x].second);
              }
              A[V_aff[x].second] = 111;
            }
            x++;
          }
        }
        use = 1 - use;
        quee[use].push(-1);
        continue;
      }

      for (int w : adj[u]) {
        if(A[w] != 111) {
          temp[w] = temp[u] + 1;
          if(use == 1 || prunable(i, w, temp, A)) {
            quee[1].push(w);
          } else {
            distances[w][i] = temp[w];
            quee[0].push(w);
          }
          A[w] = 111;
        }
      }
    }
  }

  delete [] A;
  delete [] temp;
  delete [] AFF_VERTS;
}

void HighwayLabelling::BHL(std::vector<std::pair<std::string, std::pair<int, int> > > updates) {
  uint8_t *A = new uint8_t[V];
  uint8_t *temp = new uint8_t[V];
  int *AFF_VERTS = new int[V]; int c;

  for(int i = 0; i < V; i++)
    A[i] = 111;

  for(int i = 0; i < K; i++) {

    // computing affected vertices
    std::priority_queue<std::pair<uint8_t, int>, std::vector<std::pair<uint8_t, int> >, std::greater<std::pair<uint8_t, int> > > que;
    std::queue<std::pair<uint8_t, int> > que1; c = 0;
    for(std::pair<std::string, std::pair<int, int> > iter : updates) {

      temp[iter.second.first] = query(i, iter.second.first);
      temp[iter.second.second] = query(i, iter.second.second);

      if(temp[iter.second.first] > temp[iter.second.second]) {
        if(iter.first == "EI")
          que.push(std::make_pair(temp[iter.second.second] + 1, iter.second.first));
        else if(iter.first == "ED")
          que.push(std::make_pair(temp[iter.second.first], iter.second.first));
      } else if(temp[iter.second.first] < temp[iter.second.second]) {
        if(iter.first == "EI")
          que.push(std::make_pair(temp[iter.second.first] + 1, iter.second.second));
        else if(iter.first == "ED")
          que.push(std::make_pair(temp[iter.second.second], iter.second.second));
      }
    }

    while (!que1.empty() || !que.empty()) {

      std::pair<uint8_t, int> p;
      if(!que.empty() && !que1.empty()) {
        if(que.top().first < que1.front().first) {
          p = que.top(); que.pop();
        } else {
          p = que1.front(); que1.pop();
        }
      } else if(!que1.empty()) {
        p = que1.front(); que1.pop();
      } else {
        p = que.top(); que.pop();
      }

      if(A[p.second] == 111) {

        for (int w : adj[p.second]) {
          temp[w] = query(i, w);
          if(p.first + 1 <= temp[w]) {
            que1.push(std::make_pair(p.first + 1, w));
          }
        }

        A[p.second] = p.first;
        AFF_VERTS[c] = p.second; c++;
      }
    }

    // computing boundary vertices
    std::vector<std::pair<uint8_t, int> > V_aff;
    for(int j = 0; j < c; j++) {
      temp[AFF_VERTS[j]] = 111;
      for (int w : adj[AFF_VERTS[j]]) {
        if(A[w] == 111)
          temp[AFF_VERTS[j]] = min(temp[AFF_VERTS[j]], temp[w]+1);
      }

      if(temp[AFF_VERTS[j]] != 111)
        V_aff.push_back(std::make_pair(temp[AFF_VERTS[j]], AFF_VERTS[j]));
      distances[AFF_VERTS[j]][i] = 111;
    }

    // updating the labelling
    if(V_aff.size() > 0) {
      std::sort(V_aff.begin(), V_aff.end());

      std::queue<int> quee[2]; int use = 0;
      uint8_t d = V_aff[0].first; int x = 0;
      while(x < V_aff.size() && V_aff[x].first == d) {
        if(prunable(i, V_aff[x].second, temp, A)) {
          quee[1].push(V_aff[x].second);
        } else {
          distances[V_aff[x].second][i] = temp[V_aff[x].second];
          quee[0].push(V_aff[x].second);
        }
        A[V_aff[x].second] = 111;
        x++;
      }

      quee[use].push(-1);
      while (!quee[0].empty() || x < V_aff.size()) {
        int u = quee[use].front();
        quee[use].pop();

        if(u == -1) {
          if(use == 0) { d++;
            while(x < V_aff.size() && V_aff[x].first == d) {
              if(A[V_aff[x].second] != 111) {
                if(prunable(i, V_aff[x].second, temp, A)) {
                  quee[1].push(V_aff[x].second);
                } else {
                  distances[V_aff[x].second][i] = temp[V_aff[x].second];
                  quee[0].push(V_aff[x].second);
                }
                A[V_aff[x].second] = 111;
              }
              x++;
            }
          }
          use = 1 - use;
          quee[use].push(-1);
          continue;
        }

        for (int w : adj[u]) {
          if(A[w] != 111) {
            temp[w] = temp[u] + 1;
            if(use == 1 || prunable(i, w, temp, A)) {
              quee[1].push(w);
            } else {
              distances[w][i] = temp[w];
              quee[0].push(w);
            }
            A[w] = 111;
          }
        }
      }
    }

    if(i == K - 1)
      continue;

    for(int j = 0; j < c; j++) {
      if(A[AFF_VERTS[j]] != 111)
        A[AFF_VERTS[j]] = 111;
    }
  }

  delete [] A;
  delete [] temp;
  delete [] AFF_VERTS;
}

void HighwayLabelling::BHL_Parallel(int i, std::vector<std::pair<std::string, std::pair<int, int> > > updates) {
  uint8_t *A = new uint8_t[V];
  uint8_t *temp = new uint8_t[V];
  int *AFF_VERTS = new int[V]; int c;

  for(int j = 0; j < V; j++)
    A[j] = 111;

  // computing affected vertices
  std::priority_queue<std::pair<uint8_t, int>, std::vector<std::pair<uint8_t, int> >, std::greater<std::pair<uint8_t, int> > > que;
  std::queue<std::pair<uint8_t, int> > que1; c = 0;
  for(std::pair<std::string, std::pair<int, int> > iter : updates) {

    temp[iter.second.first] = query(i, iter.second.first);
    temp[iter.second.second] = query(i, iter.second.second);

    if(temp[iter.second.first] > temp[iter.second.second]) {
      if(iter.first == "EI")
        que.push(std::make_pair(temp[iter.second.second] + 1, iter.second.first));
      else if(iter.first == "ED")
        que.push(std::make_pair(temp[iter.second.first], iter.second.first));
    } else if(temp[iter.second.first] < temp[iter.second.second]) {
      if(iter.first == "EI")
        que.push(std::make_pair(temp[iter.second.first] + 1, iter.second.second));
      else if(iter.first == "ED")
        que.push(std::make_pair(temp[iter.second.second], iter.second.second));
    }
  }

  while (!que1.empty() || !que.empty()) {

    std::pair<uint8_t, int> p;
    if(!que.empty() && !que1.empty()) {
      if(que.top().first < que1.front().first) {
        p = que.top(); que.pop();
      } else {
        p = que1.front(); que1.pop();
      }
    } else if(!que1.empty()) {
      p = que1.front(); que1.pop();
    } else {
      p = que.top(); que.pop();
    }

    if(A[p.second] == 111) {

      for (int w : adj[p.second]) {
        temp[w] = query(i, w);
        if(p.first + 1 <= temp[w])
          que1.push(std::make_pair(p.first + 1, w));
      }

      A[p.second] = p.first;
      AFF_VERTS[c] = p.second; c++;
    }
  }

  // computing boundary vertices
  std::vector<std::pair<uint8_t, int> > V_aff;
  for(int j = 0; j < c; j++) {
    temp[AFF_VERTS[j]] = 111;
    for (int w : adj[AFF_VERTS[j]]) {
      if(A[w] == 111)
        temp[AFF_VERTS[j]] = min(temp[AFF_VERTS[j]], temp[w]+1);
    }

    if(temp[AFF_VERTS[j]] != 111)
      V_aff.push_back(std::make_pair(temp[AFF_VERTS[j]], AFF_VERTS[j]));
    distances[AFF_VERTS[j]][i] = 111;
  }

  // updating the labelling
  if(V_aff.size() > 0) {
    std::sort(V_aff.begin(), V_aff.end());

    std::queue<int> quee[2]; int use = 0;
    uint8_t d = V_aff[0].first; int x = 0;
    while(x < V_aff.size() && V_aff[x].first == d) {
      if(prunable(i, V_aff[x].second, temp, A)) {
        quee[1].push(V_aff[x].second);
      } else {
        distances[V_aff[x].second][i] = temp[V_aff[x].second];
        quee[0].push(V_aff[x].second);
      }
      A[V_aff[x].second] = 111;
      x++;
    }

    quee[use].push(-1);
    while (!quee[0].empty() || x < V_aff.size()) {
      int u = quee[use].front();
      quee[use].pop();

      if(u == -1) {
        if(use == 0) { d++;
          while(x < V_aff.size() && V_aff[x].first == d) {
            if(A[V_aff[x].second] != 111) {
              if(prunable(i, V_aff[x].second, temp, A)) {
                quee[1].push(V_aff[x].second);
              } else {
                distances[V_aff[x].second][i] = temp[V_aff[x].second];
                quee[0].push(V_aff[x].second);
              }
              A[V_aff[x].second] = 111;
            }
            x++;
          }
        }
        use = 1 - use;
        quee[use].push(-1);
        continue;
      }

      for (int w : adj[u]) {
        if(A[w] != 111) {
          temp[w] = temp[u] + 1;
          if(use == 1 || prunable(i, w, temp, A)) {
            quee[1].push(w);
          } else {
            distances[w][i] = temp[w];
            quee[0].push(w);
          }
          A[w] = 111;
        }
      }
    }
  }

  delete [] A;
  delete [] temp;
  delete [] AFF_VERTS;
}

bool HighwayLabelling::prunable(int i, int u, uint8_t *temp, uint8_t *A) {

  if(landmarks.count(u) > 0) {
    highway[i][landmarks[u]] = temp[u];
    return true;
  } else {
    for (int w : adj[u]) {
      if(A[w] == 111) {
        if(temp[w] == temp[u] - 1 && distances[w][i] == 111)
        return true;
      }
    }
  }
  return false;
}

uint8_t HighwayLabelling::query(int r, int v) {

  uint8_t m = 111;
  for(int i = 0; i < K; i++) {
    m = min(m, distances_1[v][i] + highway_1[r][i]);
  }
  return m;
}

uint8_t HighwayLabelling::min(uint8_t a, uint8_t b) {
  return (a < b) ? a : b;
}

void HighwayLabelling::SelectLandmarks_HD(int topk[]) {
  std::vector<std::pair<int, int> > deg(V); long sum = 0;
  for (int v = 0; v < V; v++) {
    deg[v] = std::make_pair(adj[v].size(), v);
    sum = sum + adj[v].size();
  }
  std::sort(deg.rbegin(), deg.rend());

  for (int v = 0; v < K; v++)
    topk[v] = deg[v].second;
}

void HighwayLabelling::RemoveLandmarks(int top[]) {
  for(int i = 0; i < K; i++) {
    for (int v : adj[landmarks[i]]) {
      adj[v].erase(std::remove(adj[v].begin(), adj[v].end(), landmarks[i]));
      adj[v].shrink_to_fit();
    }
    adj[landmarks[i]].clear();
    adj[landmarks[i]].shrink_to_fit();
  }
}

uint8_t HighwayLabelling::HC_UB_naive(int s, int t) {

  uint8_t m = 111; int i, j;
  for(i = 0; i < C[s]; i++) {
    for (j = 0; j < C[t]; j++)
      m = min(m, distances[s][i] + highway[vertices[s][i]][vertices[t][j]] + distances[t][j]);
  }
  return m;
}

void HighwayLabelling::QueryDistance(std::string pairs, std::string output) {
  std::vector<TwoLayerQueue> qque; std::vector<uint8_t> qdist[2];

  qdist[0].resize(V, 111); qdist[1].resize(V, 111);
  qque.push_back(TwoLayerQueue(V)); qque.push_back(TwoLayerQueue(V));

  time_querying_millisec_ = 0; int s = 0, t = 0; int total = 0;
  std::ifstream ifs(pairs); std::ofstream ofs(output);
  while(ifs >> s >> t) { total++;

    double a = -GetCurrentTimeMilliSec();

    uint8_t dist_upper = HC_UB_naive(s, t);
    uint8_t res = dist_upper, dis[2] = {0, 0};
    for (int dir = 0; dir < 2; dir++){
      int v = dir == 0 ? s : t;
      qque[dir].clear();
      qque[dir].push(v);
      qque[dir].next();
      qdist[dir][v] = 0;
    }

    while (!qque[0].empty() && !qque[1].empty()) {
      int use = (qque[0].size() <= qque[1].size()) ? 0 : 1;
      dis[use]++;

      if (dis[0] + dis[1] == dist_upper) {
        res = dis[0] + dis[1];
        goto LOOP_END;
      }

      while (!qque[use].empty()) {

        int v = qque[use].front();
        qque[use].pop();

        for (int w : adj[v]) {

          uint8_t &src_d = qdist[    use][w];
          uint8_t &dst_d = qdist[1 - use][w];
          if (src_d != 111) continue;
          if (dst_d != 111) {
            res = qdist[use][v] + 1 + dst_d;
            goto LOOP_END;
          } else {
            qque[use].push(w);
            qdist[use][w] = qdist[use][v] + 1;
          }
	}
      }
      qque[use].next();
    }
    LOOP_END:

    a += GetCurrentTimeMilliSec();
    time_querying_millisec_ += a;

    for (int dir = 0; dir < 2; dir++) {
      for (int v : qque[dir]) {
        qdist[dir][v] = 111;
      }
      qque[dir].clear();
    }
    ofs << s << " " << t << " " << (int) min(res, dist_upper) << "\n";
  }
  ifs.close();
  ofs.close();

  std::cout << "Average Query Time (ms) : " << (double) time_querying_millisec_ / total << std::endl;

  /** freeing memory **/
  for(int i = 0; i < V; i++) {
    delete [] distances[i];
    delete [] vertices[i];
  }
  delete [] distances;
  delete [] vertices;
  delete [] C;

  for(int i = 0; i < K; i++)
    delete [] highway[i];
  delete [] highway;
}

void HighwayLabelling::loadLabelling_Full(std::string filename, int topk[]) {

  for(int i = 0; i < K; i++)
    landmarks[topk[i]] = i;

  std::ifstream ifs(std::string(filename) + std::string("_index"));

  distances = new uint8_t*[V]; distances_1 = new uint8_t*[V];
  for(int i = 0; i < V; i++) {
    distances[i] = new uint8_t[K]; distances_1[i] = new uint8_t[K];
    for(int j = 0; j < K; j++) {
      distances[i][j] = 111;
      distances_1[i][j] = 111;
    }
  }

  uint8_t C, idx;
  for(int i = 0; i < V; i++) {
    ifs.read((char*)&C, sizeof(C));
    for(uint8_t j = 0; j < C; j++) {
      ifs.read((char*)&idx, sizeof(idx));
      ifs.read((char*)&distances[i][idx], sizeof(distances[i][idx]));
      distances_1[i][idx] = distances[i][idx];
    }
  }
  ifs.close();

  ifs.open(std::string(filename) + std::string("_highway"));
  highway = new uint8_t*[K]; highway_1 = new uint8_t*[K];
  for(uint8_t i = 0; i < K; i++) {
    highway[i] = new uint8_t[K]; highway_1[i] = new uint8_t[K];
    for(uint8_t j = 0; j < K; j++) {
      ifs.read((char*)&highway[i][j], sizeof(highway[i][j]));
      highway_1[i][j] = highway[i][j];
    }
  }
  ifs.close();
}

void HighwayLabelling::loadLabelling_Pruned(std::string filename) {
  std::ifstream ifs(std::string(filename) + std::string("_index"));

  C = new uint8_t[V];
  vertices = new uint8_t*[V];
  distances = new uint8_t*[V];

  for(int i = 0; i < V; i++) {
    ifs.read((char*)&C[i], sizeof(C[i]));
    vertices[i] = new uint8_t[C[i]];
    distances[i] = new uint8_t[C[i]];
    for(uint8_t j = 0; j < C[i]; j++) {
      ifs.read((char*)&vertices[i][j], sizeof(vertices[i][j]));
      ifs.read((char*)&distances[i][j], sizeof(distances[i][j]));
    }
  }
  ifs.close();

  ifs.open(std::string(filename) + std::string("_highway"));
  highway = new uint8_t*[K];
  for(uint8_t i = 0; i < K; i++) {
    highway[i] = new uint8_t[K];
    for(uint8_t j = 0; j < K; j++)
      ifs.read((char*)&highway[i][j], sizeof(highway[i][j]));
  }
  ifs.close();
}

void HighwayLabelling::storeLabelling(std::string filename) {
  std::ofstream ofs(std::string(filename) + std::string("_index"));

  for(int i = 0; i < V; i++) {
    uint8_t C = 0;
    for(int j = 0; j < K; j++) {
      if(distances[i][j] != 111)
        C++;
    }

    ofs.write((char*)&C, sizeof(C));
    for(uint8_t j = 0; j < K; j++) {
      if(distances[i][j] != 111) {
        ofs.write((char*)&j, sizeof(j));
        ofs.write((char*)&distances[i][j], sizeof(distances[i][j]));
      }
    }
  }
  ofs.close();

  ofs.open(std::string(filename) + std::string("_highway"));
  for(int i = 0; i < K; i++) {
    for(int j = 0; j < K; j++) {
      ofs.write((char*)&highway[i][j], sizeof(highway[i][j]));
    }
  }
  ofs.close();
}
#endif  // PRUNED_LANDMARK_LABELING_H_

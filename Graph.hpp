#ifndef GRAPH_HPP
#define GRAPH_HPP
#include<bits/stdc++.h>
#include <sys/time.h>
#include <omp.h>
using namespace std;
#define mp(a,b) make_pair(a,b)
using value=float;
int global_update_nodes;
int loop_times;

vector<pair<unsigned, unsigned> > changeEdge;
struct timeval supdate_tm, eupdate_tm;

class Graph{
  public:
  Graph(){
    V=0;E=0;
  }

  struct Item{
    int node;
    value cb;
    Item()=default;
    Item(int n,value c):node(n), cb(c){};
    bool operator < (const struct Item a) const{
      if(this->cb!=a.cb)
        return this->cb>a.cb;
      else if(this->cb==a.cb){
        return this->node>a.node;
      }
    }
  };
  void ReadGraph(const char*);//file name ,mode
  void ConstructDirected();//file name ,mode
  int BaseBSearch(unsigned);
  void DeleteSearch(unsigned);
  void UpdateSMap(unsigned, unsigned, unsigned);
  //void UpdateSMap_atomic(unsigned, unsigned, unsigned, mutex*);
  void Update(unsigned, bool*);
  void EgoBWCal(unsigned, bool*);
  int OptBSearch(unsigned, value);
  void InsertUptCB(unsigned, unsigned, unsigned, set<Item>&, set<Item>&);
  void DeleteUptCB(unsigned, unsigned, unsigned, set<Item>&, set<Item>&);
  void LocalUptSmap(unsigned , unsigned, vector<unsigned>);
  void InsertUptTop(unsigned, unsigned, unsigned, set<Item>&, set<Item>&);
  void DeleteUptTop(unsigned, unsigned, unsigned, set<Item>&, set<Item>&);
  void BelongCompute(unsigned u, unsigned ,unsigned, vector<unsigned>&, bool*, set<Item>&, set<Item>&, bool); 
  void ExcludeCompute(unsigned u, unsigned ,unsigned, vector<unsigned>&, bool* , set<Item>& , set<Item>&, bool );
  void EBWEnumCpt(unsigned);
  double rand_delete(int dnum, unsigned K);
  double rand_insert(int dnum, unsigned K);
  double rand_topdelete(int dnum, unsigned K);
  double rand_topinsert(int dnum, unsigned K);
  void init_HtopR(unsigned K);
  void init_HtopR(const char* filename, unsigned k);
  void uplocalHR_top(unsigned u, unsigned K, set<Item>&, set<Item> &);
  void Destroy();
  void ResetH(int, int);
  void LocalUptSmapTop_bef(int, int, vector<unsigned>&);
  value LocalUptSmapTop(bool, unsigned, unsigned , unsigned , vector<unsigned>& );
  void updateH_CB(int, int);
  void random_generate(int);

  
  void EBWEnumCpt(const char*, int);
  void UpdateSMap_atomic(unsigned u, unsigned v, unsigned w, mutex* map_mutex);
  void ParallelEdge(const char* ,int);
  pair<int*, int> GraphColoring();

  inline void read_changes(const char* filename){
    changeEdge.clear();
    ifstream ifs;
    ifs.open(filename);
    unsigned u, v;
    while(ifs>>u>>v){
      changeEdge.emplace_back(u,v);
    }
  }

  value compute_cb_from_update(unsigned node, unsigned u, unsigned v){
    global_update_nodes++;
    if(node==u||node==v){
      unsigned another=(node==u?v:u);

      if(is_delete){
        for(auto pr:SMap[node]){
          if(pr.first.first==another||pr.first.second==another){
            Cb[node]=Cb[node]-1.0/(pr.second+1);
          }
          else{
            Cb[node]=Cb[node]-1.0/(pr.second+1)+1.0/pr.second;
          }
        }
      }
      else if(!is_delete){
        for(auto pr :SMap[node]){
          if(pr.first.first==another||pr.first.second==another){
            Cb[node]=Cb[node]+1.0/(pr.second+1);
          }
          else{
            Cb[node]=Cb[node]-1.0/pr.second+1.0/(pr.second+1);
          }
        }
      }
    }
    else{
      if(is_delete){
        for(auto pr:SMap[node]){
          if((pr.first.first==u&&pr.first.second==v)||(pr.first.first==v&&pr.first.second==u))
          {
            Cb[node]=Cb[node]+1.0/(pr.second+1);
          }
          else{
            Cb[node]=Cb[node]-1.0/(pr.second+1)+1.0/pr.second;
          }
        }
      }
      else if(!is_delete){
        for(auto pr:SMap[node]){
          if((pr.first.first==v&&pr.first.second==u)||(pr.first.first==u&&pr.first.second==v)){
            //cout << "pair(" << pr.first.first << ", " << pr.first.second << ") --- " << pr.second << endl;
            Cb[node]=Cb[node]-1.0/(pr.second+1);
          }
          else{
            //cout << "pair(" << pr.first.first << ", " << pr.first.second << ") --- " << pr.second << endl;
            Cb[node]=Cb[node]-1.0/pr.second+1.0/(pr.second+1);
          }
        }
      }
      else printf("error\n");
    }
    return Cb[node];
  }

  value compute(unsigned u){
    value res=Degree[u]*(Degree[u]-1)/2;
    for(auto i:Edge_map[u]){
      for(auto j:Edge_map[u]){
        if(i>=j) continue;
        if(Edge_map[i].find(j)!=Edge_map[i].end()) {
            res--;
          continue;
        }
        int cnt=0;
        for(auto k:Edge_map[u]){
          if(i==k||j==k) continue;
          if(Edge_map[i].find(k)!=Edge_map[i].end()&&
            Edge_map[j].find(k)!=Edge_map[j].end())
            cnt++;
        }
        if(cnt){
          res--;
          res+=1.0/(cnt+1);
        }
      }
    }
    return res;
  }
  value check(unsigned u, ofstream* ofs){
    value res=Degree[u]*(Degree[u]-1)/2;
    for(auto i:Edge_map[u]){
      for(auto j:Edge_map[u]){
        if(i>=j) continue;
        if(Edge_map[i].find(j)!=Edge_map[i].end()) {
            res--;
          continue;
        }
        int cnt=0;
        for(auto k:Edge_map[u]){
          if(i==k||j==k) continue;
          if(Edge_map[i].find(k)!=Edge_map[i].end()&&
            Edge_map[j].find(k)!=Edge_map[j].end())
            cnt++;
        }
        if(cnt){
          res--;
          res+=1.0/(cnt+1);
        }
      }
    }
    *ofs<<u<<" "<<res<<endl;
  }
  void OutputRes(vector<unsigned> R){
    printf("-------------results-------------\n");
    ofstream ofs("output.txt");
    for(auto i:R){
      cout<<i<<" "<<Cb[i]<<"\n";
    } 

  }
  void OutputRes(set<Item> R, bool outputtofile=false){
    if(outputtofile){
      ofstream ofs("output_r.txt");
      for(auto i:R)
        ofs<<i.node<<" "<< Cb[i.node] <<endl;
      ofs.close();
    }
    else{
      for(auto i:R)
        cout<<i.node<<" "<< Cb[i.node] <<endl;
    }
  }

  inline double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (tv.tv_usec / 1e6);
  }

  int getv(){

  	return V;
  }

  int gete(){

  	return E;
  }
  

  private:
  bool is_delete;
  unsigned V, E;
  vector<set<unsigned> > Edge_map;
  vector<set<unsigned> > Edge_map_directed;
  vector<map<pair<unsigned,unsigned>, unsigned>> SMap;
  map<unsigned, vector<unsigned>> record;
  unsigned* Degree;
  unsigned* count;
  //for ego-betweenness
  unsigned K;
  bool * vis;
  value* bound;
  value* Cb;
  value* Cb_origin;//do not modify it

  bool* double_vis;
  bool* u_vis;
  bool* v_vis;
  bool* my_vis;

  set<Item> Hset;//do not modify it
  set<Item> Rset;//do not modify it
  //set<Item> LHset;


  vector<bool*> global_vis;

};

void Graph::Destroy(){
    delete[] Degree;
    delete[] count;
    delete[] vis;
    delete[] bound;
    delete[] Cb;
 }

void Graph::ReadGraph(const char* filename){
  ifstream ifs;
  ifs.open(filename, ifstream::in);
  if(!ifs.is_open()){
    cout<<"error opening file\n";
  }
  unsigned a, b;
  vector<pair<unsigned, unsigned>>edges;
  V=0;
  while(ifs>>a>>b){
    edges.push_back(mp(a,b));
    if(a>V) V=a;
    if(b>V) V=b;
  }
  V++;
  Edge_map.resize(V);
  for(auto pr:edges){
    Edge_map[pr.first].insert(pr.second);
    Edge_map[pr.second].insert(pr.first);
  }
  Degree=new unsigned[V];
  for(int i=0; i<V; i++)
    Degree[i]=Edge_map[i].size();
  printf("V=%d\n", V);

  double_vis=(bool*) calloc(V,sizeof(bool));
  u_vis=(bool*) calloc(V,sizeof(bool));
  v_vis=(bool*) calloc(V,sizeof(bool));
  my_vis=(bool*) calloc(V,sizeof(bool));
}

void Graph::ConstructDirected(){
  //degree larger than myself (break ties by larger id)
  Edge_map_directed.resize(V);
  for(int i=0; i<V; i++){
    for(auto nei:Edge_map[i]){
      if(Degree[i]>Degree[nei]||(Degree[nei]==Degree[i]&&i>nei))
        Edge_map_directed[i].insert(nei);
    }
  } 
}
/**updatesmap twice ver.**/
int Graph::BaseBSearch(unsigned K){
  this->K=K;
  bound=new value[V];
  Cb=new value[V];
  SMap.resize(V);
  for(int i=0; i<V; i++){
    bound[i]=(Degree[i]*(Degree[i]-1))/2;
    Cb[i]=bound[i];
  }
  int call_compute=0;
  //vector<unsigned> R;
  set<Item> larger_first_R;
  ConstructDirected();
  set<Item> sorted_order;
  for(int i=0; i<V; i++){
    Item item(i, Cb[i]);
    sorted_order.insert(item);
  }
//start seaching
  bool* IsNeighbor=new bool[V];
  for(int i=0; i<V; i++)
    IsNeighbor[i]=false;
  
  while(!sorted_order.empty()){//the updated ones are moved
    auto cur=sorted_order.begin();
    int curnode=cur->node;

    if(larger_first_R.size()==K){
      if(larger_first_R.rbegin()->cb>=bound[curnode])
        break;
    }

    call_compute++;
    for(auto nei:Edge_map_directed[curnode]){
      IsNeighbor[nei]=true;
    }

    for(auto v:Edge_map_directed[curnode]){
      for(auto w:Edge_map_directed[v]){
        if(IsNeighbor[w]){
          #ifdef DEBUG
          cout<<"cur="<<curnode<<" v="<<v<<" w="<<w<<endl;
          #endif
          UpdateSMap(curnode,v,w);
          UpdateSMap(v,curnode,w);
          //UpdateSM(ap(w,curnode,v);
          if(curnode<v)
            SMap[w][mp(curnode,v)]=0;
          else
            SMap[w][mp(v, curnode)]=0;
          
        }
      }
    }

    for(auto nei:Edge_map_directed[curnode]){
      IsNeighbor[nei]=false;
    }
    
    #ifdef DEBUG
    printf("cur=%d: \n", curnode);
    #endif

    for(auto i:SMap[curnode]){
      if(i.first.first>i.first.second) continue;
      Cb[curnode]--;
      #ifdef DEBUG
      printf("pair<%d %d>=%d\n", i.first.first, i.first.second, i.second);
      #endif
      if(i.second){
        Cb[curnode]+=1.0/(i.second+1);
      }
    }


    if(larger_first_R.size()<K){
      Item item_(curnode, Cb[curnode]);
      //R.push_back(curnode);
      larger_first_R.insert(item_);
    }
    else{//>min
      value min_cb=larger_first_R.rbegin()->cb;
      if(Cb[curnode]>min_cb){
        larger_first_R.erase(--larger_first_R.end());
        Item item_(curnode,Cb[curnode]);
        larger_first_R.insert(item_);
      }
    }
    sorted_order.erase(cur);
  }
  //printf("tot compute time=%d\n", call_compute);
  
  OutputRes(larger_first_R);
  //Destroy();
  return call_compute;
}

void Graph::UpdateSMap(unsigned u, unsigned v, unsigned w){
  /*if(Edge_map[u].size()==2){
    SMap[u][mp(v,w)]=0;
    return ;
  }
  */
  for(auto x:Edge_map[u]){
    bool x_v=Edge_map[x].find(v)!=Edge_map[x].end();//true if in E
    bool x_w=Edge_map[x].find(w)!=Edge_map[x].end();
    if(x!=v&&x_v){
      auto it =SMap[u].find(mp(x, v));
      if(it==SMap[u].end())
        SMap[u][mp(x,v)]=0;
    }
    if(x!=w&&x_w){
      auto it =SMap[u].find(mp(x, w));
      if(it==SMap[u].end())
        SMap[u][mp(x,w)]=0;
    }
    if(x==v||x==w) continue;
    if(x!=v&&x_v&&!x_w){
      auto it=SMap[u].find(mp(x, w));
      if(it==SMap[u].end()||it->second)
        SMap[u][mp(x,w)]++;
    }
    if(!x_v&&x_w){
      auto it=SMap[u].find(mp(x, v));
      if(it==SMap[u].end()||it->second)
        SMap[u][mp(x,v)]++;
        SMap[w][mp(x,v)]++;
    }
  }
}


// void Graph::UpdateSMap_atomic(unsigned u, unsigned v, unsigned w, mutex* map_mutex){
//   map_mutex[u].lock();
//   for(auto x:Edge_map[u]){
//     bool x_v=Edge_map[x].find(v)!=Edge_map[x].end();//true if in E
//     bool x_w=Edge_map[x].find(w)!=Edge_map[x].end();
//     if(x!=v&&x_v){
//       auto it =SMap[u].find(mp(x, v));
//       if(it==SMap[u].end())
//         SMap[u][mp(x,v)]=0;
//     }
//     if(x!=w&&x_w){
//       auto it =SMap[u].find(mp(x, w));
//       if(it==SMap[u].end())
//         SMap[u][mp(x,w)]=0;
//     }
//     if(x==v||x==w) continue;
//     if(x!=v&&x_v&&!x_w){
//       auto it=SMap[u].find(mp(x, w));
//       if(it==SMap[u].end()||it->second)
//         SMap[u][mp(x,w)]++;
//     }
//     if(!x_v&&x_w){
//       map_mutex[w].lock();
//       auto it=SMap[u].find(mp(x, v));
//       if(it==SMap[u].end()||it->second)
//         SMap[u][mp(x,v)]++;
//         SMap[w][mp(x,v)]++;
//       map_mutex[w].unlock();
//     }
//   }
//   map_mutex[u].unlock();
// }

void Graph::DeleteSearch(unsigned K){
  this->K=K;
  bound=new value[V];
  Cb=new value[V];
  vis=new bool[V];
  SMap.resize(V);
  for(int i=0; i<V; i++){
    bound[i]=(Degree[i]*(Degree[i]-1))/2;
    Cb[i]=bound[i];
    vis[i]=0;
  }
  count=new unsigned[V];

  vector<unsigned> R;
  set<Item> sorted_order;
  for(int i=0; i<V; i++){
    Item item(i, Cb[i]);
    sorted_order.insert(item);
  }
//start seaching
  bool* IsDeleted=new bool[V];
  for(int i=0; i<V; i++){
    IsDeleted[i]=false;
  }

  
  while(!sorted_order.empty()){//the updated ones are moved
    auto cur=sorted_order.begin();
    int curnode=cur->node;
    if(R.size()==K){
      bool flag=true;
      for(auto v:R){
        if(Cb[v]<bound[curnode]){
          flag=false;break;
        }
      }
      if(flag) break;
    }

    #ifdef DEBUG
    printf("cur=%d: \n", curnode);
    #endif

    Update(curnode,IsDeleted);


    for(auto i:SMap[curnode]){
      Cb[curnode]--;
      if(curnode==24110)
        printf("pair<%d %d>=%d\n", i.first.first, i.first.second, i.second);
      #ifdef DEBUG
      printf("pair<%d %d>=%d\n", i.first.first, i.first.second, i.second);
      #endif
      if(i.second){
        Cb[curnode]+=1.0/(i.second+1);
      }
    }

    //printf("Cb[%d]=%lf degree=%d\n", curnode, Cb[curnode], Degree[curnode]);

    if(R.size()<K){
      R.push_back(curnode);
    }
    else{//>min
      value min_cb=Cb[R[0]];
      unsigned min_id=0;
      for(int i=0; i<R.size(); i++){
        if(Cb[R[i]]<min_cb) {
          min_cb=Cb[R[i]];
          min_id=i;
        }
      }
      if(Cb[curnode]>min_cb){
        R[min_id]=curnode;
      }
    }
    sorted_order.erase(cur);
    IsDeleted[curnode]=true;
  }
  OutputRes(R);
}

void Graph::Update(unsigned u, bool* deleted){
  vector<unsigned> neighbors;
  vector<unsigned> deleted_neighbors;
  for(auto i:Edge_map[u]){
    //count[i]=0;
    if(!deleted[i]){
      neighbors.push_back(i);
    }
    else{
      deleted_neighbors.push_back(i);
    }
    record[i].clear();
  }
  for(auto pr: SMap[u]){
    if(!pr.second){
      record[pr.first.first].push_back(pr.first.second);
      record[pr.first.second].push_back(pr.first.first);
    }
  }

  //compute the special case
  //deleted nodes are in ascending order
  for(int pos_i=0; pos_i<deleted_neighbors.size(); pos_i++){
    unsigned i=deleted_neighbors[pos_i];
    for(auto pr:record[i]){
      vis[pr]=true;
    }

    for(int pos_j=pos_i+1; pos_j<deleted_neighbors.size(); pos_j++){
      unsigned j=deleted_neighbors[pos_j];
      if(vis[j]==true) continue;
      for(auto pr_:record[j]){
        if(vis[pr_]&&!deleted[pr_]){
          SMap[u][mp(i,j)]++;
          SMap[pr_][mp(i,j)]++;
        }
      }
    }
    for(auto pr:record[i]){
      vis[pr]=false;
    }
  }
  for(int idx_i=0; idx_i<neighbors.size(); idx_i++){
    unsigned i=neighbors[idx_i];
    for(int idx_j=idx_i+1; idx_j<neighbors.size(); idx_j++){
      unsigned j=neighbors[idx_j];
      if(Edge_map[i].find(j)!=Edge_map[i].end()){
        //decide order:first< second
        SMap[u][mp(i,j)]=0;
        if(u<j) SMap[i][mp(u,j)]=0;
        else SMap[i][mp(j, u)]=0;
        if(u<i) SMap[j][mp(u,i)]=0;
        else SMap[j][mp(i, u)]=0;
        //
        for(auto n:record[i]){
          if(j<n){
            auto it=SMap[u].find(mp(j,n));
            if(it==SMap[u].end()||it->second){
              SMap[u][mp(j,n)]++;
              SMap[i][mp(j,n)]++;
            }
          } 
          else{
            auto it=SMap[u].find(mp(n, j));
            if(it==SMap[u].end()||it->second){
              SMap[u][mp(n,j)]++;
              SMap[i][mp(n,j)]++;

            }
          }
        }
        for(auto n:record[j]){
          if(i<n){
            auto it=SMap[u].find(mp(i,n));
            if(it==SMap[u].end()||it->second){
              SMap[u][mp(i,n)]++;
              SMap[j][mp(i,n)]++;

            }
          } 
          else{
            auto it=SMap[u].find(mp(n, i));
            if(it==SMap[u].end()||it->second){
              SMap[u][mp(n, i)]++;
              SMap[j][mp(n, i)]++;

            }
          }
        }
        record[i].push_back(j);
        record[j].push_back(i);
      }
    }
  }

}

void Graph::EgoBWCal(unsigned v, bool*deleted){
  vector<unsigned> Neighbors;
  vector<unsigned> deleted_neighbors;
  for(auto i:Edge_map[v]){
    //count[i]=0;
    if(!deleted[i]){
      Neighbors.push_back(i);
    }
    else
      deleted_neighbors.push_back(i);
    record[i].clear();
  }
  for(auto pr: SMap[v]){
    if(!pr.second){
      record[pr.first.first].push_back(pr.first.second);
      record[pr.first.second].push_back(pr.first.first);
    }
  }

  for(int pos_i=0; pos_i<deleted_neighbors.size(); pos_i++){
    unsigned i=deleted_neighbors[pos_i];
    for(auto pr:record[i]){
      vis[pr]=true;
    }
    for(int pos_j=pos_i+1; pos_j<deleted_neighbors.size(); pos_j++){
      unsigned j=deleted_neighbors[pos_j];
      if(vis[j]==true) continue;
      for(auto pr_:record[j]){
        if(vis[pr_]&&!deleted[pr_]){
          SMap[v][mp(i,j)]++;
          SMap[pr_][mp(i,j)]++;
        }
      }
    }
    for(auto pr:record[i]){
      vis[pr]=false;
    }
  }

  for(int idx_i=0; idx_i<Neighbors.size(); idx_i++){
    int i=Neighbors[idx_i];
    for(int idx_j=idx_i; idx_j<Neighbors.size(); idx_j++){
      int j=Neighbors[idx_j];
      if(idx_i==idx_j) continue;
      for(int idx_k=0; idx_k<Neighbors.size(); idx_k++){
        if(idx_k==idx_i||idx_k==idx_j) continue;
        int k=Neighbors[idx_k];
        //if(v==25758) printf("i=%d, j=%d, k=%d\n", i, j, k);
        bool i_j=(Edge_map[i].find(j)!=Edge_map[i].end());
        bool i_k=(Edge_map[i].find(k)!=Edge_map[i].end());
        bool j_k=(Edge_map[j].find(k)!=Edge_map[j].end());
        if(i_j)  {
          if(i<j) SMap[v][mp(i,j)]=0;
          else SMap[v][mp(j,i)]=0;
          if(v<j) SMap[i][mp(v,j)]=0;
          else SMap[i][mp(j,v)]=0;
          if(v<i) SMap[j][mp(v,i)]=0;
          else SMap[j][mp(i,v)]=0;
        }
        if(j_k){
          if(j<k)  SMap[v][mp(j,k)]=0;
          else SMap[v][mp(k,j)]=0;
          if(v<k) SMap[j][mp(v,k)]=0;
          else SMap[j][mp(k,v)]=0;
          if(v<j) SMap[k][mp(v,j)]=0;
          else SMap[k][mp(j,v)]=0;
        }
        if(i_k){
          if(i<k) SMap[v][mp(i,k)]=0;
          else SMap[v][mp(k,i)]=0;
          if(v<k) SMap[i][mp(v, k)]=0;
          else SMap[i][mp(k,v)]=0;
          if(v<i) SMap[k][mp(v,i)]=0;
          else SMap[k][mp(i,v)]=0;
        }
        if(i_j&&i_k&&!j_k){
          bool flag=false;
          if(j>k){
            swap(j,k);flag=true;
          }
          if(SMap[i].find(mp(j,k))==SMap[i].end()){
            SMap[i][mp(j,k)]=0;
          }
          if(SMap[v].find(mp(j,k))==SMap[v].end()){
            SMap[v][mp(j,k)]=0;
          }
          SMap[i][mp(j,k)]++;SMap[v][mp(j,k)]++;
          if(flag)
            swap(j,k);
        }
        if(i_j&&j_k&&!i_k){
          bool flag=false;
          if(i>k){
            swap(i,k);flag=true;
          }
          if(SMap[j].find(mp(i,k))==SMap[j].end())
            SMap[j][mp(i, k)]=0;
          if(SMap[v].find(mp(i,k))==SMap[v].end())
            SMap[v][mp(i, k)]=0;
          SMap[j][mp(i, k)]++;SMap[v][mp(i, k)]++;
          if(flag) swap(i, k);
        }
        if(i_k&&j_k&&!i_j){
          bool flag=false;
          if(i>j){
            swap(i, j); flag=true;
          }
          if(SMap[k].find(mp(i,j))==SMap[k].end())
            SMap[k][mp(i,j)]=0;
          if(SMap[v].find(mp(i,j))==SMap[v].end())
            SMap[v][mp(i,j)]=0;
          SMap[k][mp(i,j)];SMap[v][mp(i,j)]++;
          if(flag) swap(i, j);
        }
      }
    }
  }
}

int Graph::OptBSearch(unsigned K, value theta){
  this->K=K;
  bound=new value[V];
  Cb=new value[V];
  Cb_origin=new value[V];
  vis=new bool[V];
  SMap.resize(V);
  for(int i=0; i<V; i++){
    bound[i]=(Degree[i]*(Degree[i]-1))/2;
    Cb[i]=bound[i];
    vis[i]=0;
  }
  int call_compute=0;
  //vector<unsigned> R;
  set<Item> large_first_R;
  
  //ConstructDirected();
  set<Item> sorted_order;
  for(int i=0; i<V; i++){
    Item item(i, Cb[i]);
    sorted_order.insert(item);
  }
  //start seaching
  bool* IsDeleted=new bool[V];
  for(int i=0; i<V; i++){
    IsDeleted[i]=false;
  }
  unsigned pre_fetch=V+1;
  while(!sorted_order.empty()){//the updated ones are moved
    auto cur=sorted_order.begin();
    sorted_order.erase(cur);
    int curnode=cur->node;
    value curcb=cur->cb;
    //compute new_bound
    value new_bound=Degree[curnode]*(Degree[curnode]-1)/2;
    for(auto pr:SMap[curnode]){
      new_bound--;
      if(pr.second)
        new_bound+=1.0/(pr.second+1);
    }
 // printf("curnode=%u, new bound=%lf, curcb=%lf, deg=%d\n",curnode, new_bound,curcb,Degree[curnode]);
    value min_cb;
    unsigned min_id;
    
    if(theta*new_bound < curcb){
      if(large_first_R.size()<K){
        Item item(curnode, new_bound);
        sorted_order.insert(item);
        //continue;
      }
      else{
        min_cb=large_first_R.rbegin()->cb;
        min_id=large_first_R.rbegin()->node;
        if(new_bound>min_cb){
          Item item(curnode, new_bound);
          sorted_order.insert(item);
          // printf("bool=%d\n",sorted_order.find(item)!=sorted_order.end());
          //continue;
        }
      }
      continue;
    }
    if(!large_first_R.empty())
      min_cb=large_first_R.rbegin()->cb;
    if(large_first_R.size()==K&&curcb<=min_cb) break;
    //exact search
    call_compute++;
    Update(curnode, IsDeleted);
    #ifdef DEBUG
    printf("cur=%d: \n", curnode);
    #endif
    for(auto i:SMap[curnode]){
      if(i.first.first>i.first.second) continue;

      Cb[curnode]--;
      #ifdef DEBUG
      printf("pair<%d %d>=%d\n", i.first.first, i.first.second, i.second);
      #endif
      if(i.second){
        Cb[curnode]+=1.0/(i.second+1);
      }
    }
    if(large_first_R.size()<K){
      Item new_item(curnode, Cb[curnode]);
      large_first_R.insert(new_item);
     // printf("add\n");
    }
    else{//>min
      // min_cb=Cb[R[0]];
      // min_id=0;
      min_cb=large_first_R.rbegin()->cb;
      min_id=large_first_R.rbegin()->node;
      // for(int i=0; i<R.size(); i++){
      //   if(Cb[R[i]]<min_cb) {
      //     min_cb=Cb[R[i]];
      //     min_id=i;
      //   }
      // }
      if(Cb[curnode]>min_cb){
        large_first_R.erase(--large_first_R.end());
        Item new_item(curnode, Cb[curnode]);
        large_first_R.insert(new_item);
       // printf("add\n");
      }
    }
    IsDeleted[curnode]=true;
  }
  //printf("tot compute time=%d\n", call_compute);
  //OutputRes(large_first_R);
  Hset = large_first_R;
  //LHset = large_first_R;
  for(int i=0; i<V; i++){
    Cb_origin[i]=Cb[i];
  }
  //Destroy();
  return call_compute;
}

void Graph::InsertUptCB(unsigned add_u, unsigned add_v, unsigned K, set<Item>& H, set<Item>& R){
  
  SMap.resize(V);
  for(int i=0; i<V; i++) SMap[i].clear();
  //insert(u, v) to graph
  Edge_map[add_u].insert(add_v);
  Edge_map[add_v].insert(add_u);
  Degree[add_u]++;
  Degree[add_v]++;

  //find intersection of N(u) and N(v)
  vector<unsigned> vec_L;
  for(int i=0; i<V; i++) vis[i]=false;
  for(auto i:Edge_map[add_v])
    vis[i]=true;
  for(auto i:Edge_map[add_u]){
    if(vis[i])
      vec_L.push_back(i);
  }
  for(auto i:Edge_map[add_v])
    vis[i]=false;
    /*
  printf("modify %d %d",add_u, add_v);
  for(auto i:vec_L)
    printf(" %d", i);
  printf("\n");
  */
  // printf("common neighbors: ");
  // for(auto i:vec_L){
  //   printf("%u ", i);
  // }
  // printf("\n");
  //update smap
  LocalUptSmap(add_u, add_v, vec_L);
  // printf("pr 350: ");
  // for(auto pr: SMap[350]){
  //   printf("<%u, %u>=%d\n", pr.first.first, pr.first.second, pr.second);
  // }
  for(int i=0; i<2; i++){
    unsigned u, v;
    if(i==0) u=add_u,v=add_v;
    else u=add_v, v=add_u; 
    global_update_nodes++;
    for(auto pr :SMap[u]){
      if(pr.first.first==v||pr.first.second==v){
        //cout << "pair(" << pr.first.first << ", " << pr.first.second << ") --- " << pr.second << endl;
        Cb[u]=Cb[u]+1.0/(pr.second+1);
      }
      else{
      	//cout << "pair(" << pr.first.first << ", " << pr.first.second << ") --- " << pr.second << endl;
        Cb[u]=Cb[u]-1.0/pr.second+1.0/(pr.second+1);
      }
    }
    uplocalHR_top(u, K ,H, R );
  }
  //printf("Cb[%u] =%lf, Cb[%u]=%lf\n", add_u, Cb[add_u], add_v, Cb[add_v]);
  for(auto x: vec_L){
    bool flag=false;
    global_update_nodes++;
    for(auto pr:SMap[x]){
      //if(x==1028) printf("(%d %d) %d\n", pr.first.first, pr.first.second, pr.second);
      if((pr.first.first==add_v&&pr.first.second==add_u)||(pr.first.first==add_u&&pr.first.second==add_v)){
        //cout << "pair(" << pr.first.first << ", " << pr.first.second << ") --- " << pr.second << endl;
        Cb[x]=Cb[x]-1.0/(pr.second+1);
        flag=true;
      }
      else{
      	//cout << "pair(" << pr.first.first << ", " << pr.first.second << ") --- " << pr.second << endl;
        Cb[x]=Cb[x]-1.0/pr.second+1.0/(pr.second+1);
      }
      //if(x==1028) printf("origin=%lf, new=%lf\n", Cb_origin[1028], Cb[1028]);
    }
    uplocalHR_top(x, K, H, R);
    //printf("Cb[%u]= %lf\n", x, Cb[x]);
  }
  //printf("Cb_origin[145]=%lf, Cb[145]=%lf\n", Cb_origin[145], Cb[145]);
  OutputRes(R, true);
}

void Graph::DeleteUptCB(unsigned delete_u, unsigned delete_v, unsigned K, set<Item>& H, set<Item>& R){
  SMap.resize(V);
  for(int i=0; i<V; i++) SMap[i].clear();

  //find intersection of N(u) and N(v)
  vector<unsigned> vec_L;
  for(int i=0; i<V; i++) vis[i]=false;
  for(auto i:Edge_map[delete_v])
    vis[i]=true;
  for(auto i:Edge_map[delete_u]){
    if(vis[i])
      vec_L.push_back(i);
  }
  for(auto i:Edge_map[delete_v])
    vis[i]=false;

  //update smap
  LocalUptSmap(delete_u, delete_v, vec_L);


  for(int i=0; i<2; i++){
    global_update_nodes++;
    unsigned u, v;
    if(i==0) u=delete_u, v=delete_v;
    else u=delete_v, v=delete_u;
    for(auto pr:SMap[u]){
      if(pr.first.first==v||pr.first.second==v){
        //cout << "&&&" << pr.second << endl;
        Cb[u]=Cb[u]-1.0/(pr.second+1);
      }
      else{
      	//cout << "&&&" << pr.second << endl;
        Cb[u]=Cb[u]-1.0/(pr.second+1)+1.0/pr.second;
      }
    }
    uplocalHR_top(u, K, H ,R);
  }
  //printf("Cb[%u] =%lf, Cb[%u]=%lf\n", delete_u, Cb[delete_u], delete_v, Cb[delete_v]);
  //cout << "---" << endl;
  for(auto x:vec_L){
    global_update_nodes++;
    for(auto pr:SMap[x]){
      if((pr.first.first==delete_u&&pr.first.second==delete_v)||(pr.first.first==delete_v&&pr.first.second==delete_u))
      {
        //cout << "&&&" << pr.second << endl;
        Cb[x]=Cb[x]+1.0/(pr.second+1);
      }
      else{
      	//cout << "&&&" << pr.second << endl;
        Cb[x]=Cb[x]-1.0/(pr.second+1)+1.0/pr.second;
      }
    }
    uplocalHR_top(x, K, H, R);
    //printf("Cb[%u]= %lf\n", x, Cb[x]);
  }
  //cout << "ok" << endl;
  //delete u and v from graph
  Edge_map[delete_u].erase(delete_v);
  Edge_map[delete_v].erase(delete_u);
  Degree[delete_u]--;
  Degree[delete_v]--;

  OutputRes(R, true);
  //cout << "ok" << endl;
      
}

void Graph::LocalUptSmap(unsigned add_u, unsigned add_v, vector<unsigned> L){
  if(add_u>add_v)
    swap(add_u,add_v);
  double exec_time=0;
  exec_time-=get_time();
  unsigned* local_vis= new unsigned[V];
  for(int i=0; i<V; i++) local_vis[i]=0;
  for(auto l: L){
    SMap[l][mp(add_u, add_v)]=0;
  }
  
  for(auto nei:Edge_map[add_v])
    local_vis[nei]=1;
  for(auto nei:Edge_map[add_u]){
    if(local_vis[nei]) continue;
    if(add_v<nei)
      SMap[add_u][mp(add_v,nei)]=0;
    else if(add_v>nei)
      SMap[add_u][mp(nei, add_v)]=0;
  }
  for(auto nei:Edge_map[add_v])
    local_vis[nei]=0;
  
  for(auto nei:Edge_map[add_u])
    local_vis[nei]=1;
  for(auto nei:Edge_map[add_v]){
    if(local_vis[nei]) continue;
    if(add_u<nei)
      SMap[add_v][mp(add_u, nei)]=0;
    else if(add_u>nei)
      SMap[add_v][mp(nei, add_u)]=0;
  }
  for(auto nei:Edge_map[add_u])
    local_vis[nei]=0;

  for(auto p: L){
    for(auto nei: Edge_map[p])
      vis[nei]=true;
    for(int ii=0; ii<2; ii++){
      unsigned u, v;
      if(ii==0) u=add_u, v=add_v;
      else u=add_v, v=add_u;
      for(auto x:Edge_map[u]){
        if(vis[x]){
          if(Edge_map[x].find(v)==Edge_map[x].end()&&x!=v){
            if(x<v){
              SMap[u][mp(x,v)]++;
            }
            else{
              SMap[u][mp(v, x)]++;
            }
            
            vector<unsigned> tool_vec={x,v,p};
            unsigned max_degree=0;
            int max_idx=-1;
            for(auto i:tool_vec){
              if(Degree[i]>max_degree){
                max_degree=Degree[i];
                max_idx=i;
              }
            }        
            assert(max_idx!=-1);

            for(auto nei: Edge_map[max_idx]){
              local_vis[nei]=1;
            }
            for(auto i:tool_vec){
              if(i!=max_idx){
                for(auto nei:Edge_map[i]){
                  if(local_vis[nei])
                    local_vis[nei]++;
                }
              }
            }

            int sum=0;
            for(auto nei:Edge_map[max_idx]){
              if(local_vis[nei]==3) sum++;
              local_vis[nei]=0;
            }
            if(x<v)
              SMap[p][mp(x,v)]=sum;
            else
              SMap[p][mp(v,x)]=sum;

          }
        }
      }
    }
    for(auto q:L){
      //N[q] intersection N[p]
      if(Edge_map[p].find(q)==Edge_map[p].end()&&p>q){

        for(auto nei:Edge_map[q]){
          if(vis[nei]==true){
            local_vis[nei]=2;
          }
        }

        int sum=0;
        for(auto nei:Edge_map[add_u]){
          if(local_vis[nei])
            sum++;
        }
        SMap[add_u][mp(q,p)]=sum;
        sum=0;
        for(auto nei:Edge_map[add_v]){
          if(local_vis[nei])
            sum++;
        }
        SMap[add_v][mp(q,p)]=sum;

        for(auto nei:Edge_map[q]){
          if(vis[nei]==true){
            local_vis[nei]=0;
          }
        }

      }

      if(Edge_map[p].find(q)!=Edge_map[p].end()&&p<q){
        SMap[p][mp(add_u, add_v)]++;
        SMap[q][mp(add_u, add_v)]++;
      }
    }

    for(auto nei: Edge_map[p])
      vis[nei]=false;
  }

  delete[] local_vis;

  exec_time += get_time();
 // printf("exec_time for updating =%lf(s)\n", exec_time);
}

void Graph::InsertUptTop(unsigned add_u, unsigned add_v, unsigned K, set<Item>& H, set<Item>& R){
  //strcut H
  //bool* flag=new bool[V];
  is_delete=false;
  SMap.resize(V);
  for(int i=0; i<V; i++) SMap[i].clear();

  Edge_map[add_u].insert(add_v);
  Edge_map[add_v].insert(add_u);
  Degree[add_u]++;
  Degree[add_v]++;
///modify
  //find intersection of N(u) and N(v)
  vector<unsigned> vec_L;
  for(int i=0; i<V; i++) vis[i]=false;
  for(auto i:Edge_map[add_v])
    v_vis[i]=true;
  for(auto i:Edge_map[add_u]){
    if(v_vis[i])
      vec_L.push_back(i);
  }
  for(auto i:Edge_map[add_u])
    u_vis[i]=true;
  //update smap
  //printf("u=%d, v=%d\n", add_u, add_v);
  LocalUptSmapTop_bef(add_u, add_v, vec_L);

  bool* flag= (bool*)calloc(V, sizeof(bool));

  //vis= new bool[V];
  for(int i=0; i<V; i++) vis[i]=false;

  for(auto i:R) vis[i.node]=true;

  for(int i=0; i<2; i++){
    unsigned u;
    if(i==0) u=add_u;
    else u=add_v;

    if(vis[u]){// belong to R
      BelongCompute(u, add_u, add_v, vec_L, flag, R, H, 1);
    }
    else{
      value upperbound=Degree[u]*(Degree[u]-1)/2.0;
      if(upperbound>R.rbegin()->cb){
        value new_cb=LocalUptSmapTop(1, u, add_u, add_v, vec_L);
        //printf("new cb of %d is %lf\n", u, new_cb);
        for(auto it=H.begin(); it!=H.end(); it++){
          if(it->node==u){
            H.erase(it);break;
          }
        }
        H.insert(Item(u, new_cb));
        flag[u]=false;
        if(new_cb>R.rbegin()->cb){
          vis[R.rbegin()->node]=false;
          vis[u]=true;
          R.erase(*R.rbegin());
          R.insert(Item(u, new_cb));
        }
      }
      else
        flag[u]=true;
    }

  }
  for(auto x:vec_L){
    if(vis[x]){
      BelongCompute(x, add_u, add_v, vec_L, flag, R, H, 1);
    }
    else{
      flag[x]=true;
    }
  }

  //OutputRes(R_set);
  for(auto i:Edge_map[add_v])
    v_vis[i]=false;
  for(auto i:Edge_map[add_u])
    u_vis[i]=false;
    //printf("Cb_origin[145]=%lf, Cb[145]=%lf\n", Cb_origin[145], Cb[145]);

  OutputRes(R, true);
}

void Graph::DeleteUptTop(unsigned delete_u, unsigned delete_v, unsigned K, set<Item>& H, set<Item>& R){
  //strcut H
  //bool* flag=new bool[V];
  is_delete=true;
  SMap.resize(V);
  for(int i=0; i<V; i++) SMap[i].clear();
///modify
  //find intersection of N(u) and N(v)
  vector<unsigned> vec_L;
  for(int i=0; i<V; i++) vis[i]=false;
  for(auto i:Edge_map[delete_v])
    v_vis[i]=true;
  for(auto i:Edge_map[delete_u]){
    if(v_vis[i])
      vec_L.push_back(i);
  }
  for(auto i:Edge_map[delete_u])
    u_vis[i]=true;
  //update smap

  LocalUptSmapTop_bef(delete_u, delete_v, vec_L);

  bool* flag= (bool*)calloc(V, sizeof(bool));

  for(auto i:R){
    vis[i.node]=true;
  } 

  for(int i=0; i<2; i++){
    unsigned u;
    if(i==0) u=delete_u;
    else u=delete_v; 

    if(vis[u]){
      //printf("%u belong\n", u);
      BelongCompute(u, delete_u,  delete_v, vec_L,flag, R, H, 0);
    }
    else{
      //printf("%u exclude\n", u);
      ExcludeCompute(u, delete_u,  delete_v, vec_L, flag, R, H, 0);
    }
  }

  for(auto x:vec_L){
    if(!vis[x]){
      ExcludeCompute(x, delete_u, delete_v, vec_L, flag, R, H, 0);
    }
    else{
      flag[x]=true;
    }
  }
  Edge_map[delete_u].erase(delete_v);
  Edge_map[delete_v].erase(delete_u);
  Degree[delete_v]--;
  Degree[delete_u]--;
  for(auto i:Edge_map[delete_v])
    v_vis[i]=false;
  for(auto i:Edge_map[delete_u])
    u_vis[i]=false;
  OutputRes(R, true);
}

void Graph::BelongCompute(unsigned u, unsigned node_u, unsigned node_v, vector<unsigned>& vec_L, bool* flag, set<Item>& R_set, set<Item>& H, bool type){

  vector<Item> new_add;
  vector<Item> new_delete;
  value cb_new=LocalUptSmapTop(type, u, node_u, node_v, vec_L);

  //printf("new cb of %u is %lf\n", u, cb_new);
  new_add.push_back(Item(u, cb_new));
  flag[u]=false;

//here
  for(auto it=H.begin(); it!=H.end(); it++){
    if(it->node==u){
      H.erase(it);break;
    }
  }
  for(auto it: R_set){
    if(it.node==u){
      R_set.erase(it);break;
    }
  }
  R_set.insert(Item(u, cb_new));
  value min_val=V*(V-1)/2.0;
  for(auto it=R_set.rbegin(); it!=R_set.rend(); it++){
    if(it->node!=u&&it->cb<min_val){
      min_val=it->cb;break;
    }
  }
  assert(min_val!=V*(V-1)/2.0);
  //printf("node=%d, cb[new]=%lf, min_val=%lf\n",u, Cb[u], min_val);
  if(cb_new<min_val){
    //sort H;
    //printf("yes\n");
    for(auto item_:H){
      loop_times++;
      //printf("item=%d\n", item_.node);
      if(!vis[item_.node]){
        if(flag[item_.node]==false&&cb_new<item_.cb){ 
          //erase add_u
          auto it=R_set.begin();
          for(; it!=R_set.end(); it++){
            if(it->node==u){
              break;
            }
          }
          vis[it->node]=false;
          vis[item_.node]=true;
          R_set.erase(it);
          R_set.insert(item_);
          break;
        }
        else
        {
          if(flag[item_.node]==true){
            cb_new=LocalUptSmapTop(type, u, node_u, node_v, vec_L);
            flag[item_.node]=false;
            new_add.push_back(Item(item_.node, cb_new));
            new_delete.push_back(item_);
          }
        }
        
      }
    }
  }
  for(auto item_:new_delete)
    H.erase(item_);
  for(auto item_:new_add){
    H.insert(item_);
  }
}

void Graph::ExcludeCompute(unsigned u, unsigned node_u, unsigned node_v, vector<unsigned>& vec_L, bool* flag, set<Item>& R_set, set<Item>& H, bool type){
  value upperbound=Degree[u]*(Degree[u]-1)/2.0;
  if(upperbound>R_set.rbegin()->cb){
    value new_cb=LocalUptSmapTop(type, u, node_u, node_v, vec_L);
      //if(u==2552) printf("modify %lf\n", compute(2552));

    for(auto it=H.begin(); it!=H.end(); it++){
      if(it->node==u){
        H.erase(it);break;
      }
    }
    H.insert(Item(u, new_cb));
    flag[u]=false;
    while(flag[R_set.rbegin()->node]==true){
      unsigned node_=R_set.rbegin()->node;
      value new_cb_=LocalUptSmapTop(type, u, node_u, node_v, vec_L);
      flag[node_]=false;
      for(auto ii: H){
        if(ii.node==node_){
          H.erase(ii); break;
        }
      }
      H.insert(Item(node_, new_cb_));
      for(auto ii:R_set){
        if(ii.node=node_){
          R_set.erase(ii); break;
        }
      }
      R_set.insert(Item(node_, new_cb_));
    }
    if(new_cb>R_set.rbegin()->cb){
      vis[R_set.rbegin()->node]=false;
      vis[u]=true;
      R_set.erase(*R_set.rbegin());
      R_set.insert(Item(u, new_cb));
    }
  }
  else{
    flag[u]=true;
  }
}                  


void Graph::EBWEnumCpt(unsigned k){
  SMap.resize(V);
  Cb=new value[V];
  ConstructDirected();
  mutex* map_mutex= new mutex[V];

  #pragma omp parallel for 
  for(int i=0; i<V; i++){
    bool* vis = new bool[V];
    for(int j=0; j<V; j++)
      vis[j]=false;
    for(auto v: Edge_map_directed[i])
      vis[v]=true;
    for(auto v: Edge_map_directed[i]){
      for(auto w:Edge_map_directed[v]){
        if(vis[w]==true){
          UpdateSMap_atomic(i, v, w, map_mutex);
          UpdateSMap_atomic(v, i, w, map_mutex);
          map_mutex[w].lock();
          if(i<v) SMap[w][mp(i,v)]=0;
          else SMap[w][mp(v, i)]=0;
          map_mutex[w].unlock();
        }
      }
    }
    for(auto v: Edge_map_directed[i])
      vis[v]=false;
  }
  #pragma omp parallel for
  for(int curnode=0; curnode<V; curnode++){
    for(auto i:SMap[curnode]){
      if(i.first.first>i.first.second) continue;
      Cb[curnode--];
      if(i.second){
        Cb[curnode]+=1.0/(i.second+1);
      }
    }
  }
  
}

void Graph::random_generate(int dnum){
  srand((unsigned) time(0));
    //for delete
  int cnt = 0;
  while(cnt < dnum){
      unsigned u = (rand() % (V-1-0+1))+ 0;
      unsigned v;
      if(Edge_map[u].size() == 0){
        continue;
      }
      unsigned vid = (rand() % (Edge_map[u].size()-1-0+1))+0;
      set<unsigned>::iterator itedge = Edge_map[u].begin();
      int incnt = 0;
      while(itedge != Edge_map[u].end()){
        if(incnt == vid){
          v = *itedge;
          break;
        }
        incnt++;
        itedge++;
      }
      if(itedge == Edge_map[u].end()){
        continue;
      }
      //cout << "delete edge (" << u << ", " << v << ")" << endl; 
      changeEdge.emplace_back(u, v);
      cnt++;    
  }
  ofstream ofs("delete_edges.txt");
  for(auto pr:changeEdge){
    ofs<<pr.first<<" "<<pr.second<<endl;
  }
  ofs.close();
//for insert
  cnt=0;
  ofs.open("insert_edges.txt");
  set<pair<unsigned, unsigned> >avoid_dup;
  while(cnt<dnum){
      unsigned u = (rand() % (V-1-0+1))+ 0;
      unsigned v= (rand() % (V-1-0+1))+ 0;
      while(Edge_map[u].find(v)!=Edge_map[u].end()){
        v=(rand() % (V-1-0+1))+ 0;
      }
      if(avoid_dup.find(make_pair(u, v))==avoid_dup.end()&&avoid_dup.find(make_pair(v,u))==avoid_dup.end()){
        avoid_dup.insert(make_pair(u,v));
      }
      ofs<<u<<" "<<v<<endl;
      cnt++;
  }
  printf("write done\n");
  

}

double Graph::rand_delete(int dnum, unsigned K){
    double time = 0.0;
    read_changes("delete_edges.txt");
	  gettimeofday(&supdate_tm, NULL);
    for(int i = 0; i < dnum; i++){

        int u = changeEdge[i].first;
        int v = changeEdge[i].second;
        //cout << "insert edge (" << u << ", " << v << ")" << endl; 
        set<Item> H=Hset;
        for(int j=0; j<V; j++)
          Cb[j]=Cb_origin[j];
        set<Item> R=Rset;
        DeleteUptCB(u, v, K, H, R);
        Edge_map[u].insert(v);
        Edge_map[v].insert(u);
        Degree[u]++;
        Degree[v]++;
    } 
    gettimeofday(&eupdate_tm, NULL);
    time += double(eupdate_tm.tv_sec - supdate_tm.tv_sec) + double(eupdate_tm.tv_usec - supdate_tm.tv_usec) / 1000000.0f;
   	printf("global uptdate nodes =%d\n", global_update_nodes);
     // for(int i = 0; i < V; i++){
    //   cout<<i<<" "<<Cb[i]<<"\n";
    // } 
    return time/dnum;
}

double Graph::rand_insert(int dnum, unsigned K){
    double time = 0.0;
    read_changes("insert_edges.txt");
    gettimeofday(&supdate_tm, NULL);
    global_update_nodes=0;
    for(int i = 0; i < dnum; i++){

        int u = changeEdge[i].first;
        int v = changeEdge[i].second;
        //cout << "insert edge (" << u << ", " << v << ")" << endl; 
        set<Item> H=Hset;
        for(int j=0; j<V; j++)
          Cb[j]=Cb_origin[j];
        set<Item>R= Rset;
        InsertUptCB(u, v, K, H, R);
        Edge_map[u].erase(v);
        Edge_map[v].erase(u);
        Degree[u]--;
        Degree[v]--;
    }
    gettimeofday(&eupdate_tm, NULL);
    time = double(eupdate_tm.tv_sec - supdate_tm.tv_sec) + double(eupdate_tm.tv_usec - supdate_tm.tv_usec) / 1000000.0f;
    printf("global uptdate nodes =%d\n", global_update_nodes);
    // for(int i = 0; i < V; i++){
    //   cout<<i<<" "<<Cb[i]<<"\n";
    // } 
    return time/dnum;
}


double Graph::rand_topdelete(int dnum, unsigned K){
    // changeEdge.clear();
    read_changes("delete_edges.txt");
    double time = 0.0;
    gettimeofday(&supdate_tm, NULL);
    global_update_nodes=0;
    loop_times=0;

	  for(int i = 0; i < dnum; i++){ 
	    int u = changeEdge[i].first;
	    int v = changeEdge[i].second;
	    //cout << "insert edge (" << u << ", " << v << ")" << endl; 
      set<Item> H=Hset;
      for(int j=0; j<V; j++)
        Cb[j]=Cb_origin[j];
      set<Item>R= Rset;

	    DeleteUptTop(u, v, K, H, R);
      Edge_map[u].insert(v);
      Edge_map[v].insert(u);
      Degree[u]++;
      Degree[v]++;
	 } 
    gettimeofday(&eupdate_tm, NULL);
    time += double(eupdate_tm.tv_sec - supdate_tm.tv_sec) + double(eupdate_tm.tv_usec - supdate_tm.tv_usec) / 1000000.0f;
    printf("global uptdate nodes =%d\n", global_update_nodes);
  // gettimeofday(&eupdate_tm, NULL);
   //time += double(eupdate_tm.tv_sec - supdate_tm.tv_sec) + double(eupdate_tm.tv_usec - supdate_tm.tv_usec) / 1000000.0f;

   	// for(int i = 0; i < V; i++){
    //   cout<<i<<" "<<Cb[i]<<"\n";
    // } 
  return time/dnum;
}

double Graph::rand_topinsert(int dnum, unsigned K){
    double time = 0.0;
    read_changes("insert_edges.txt");
    //gettimeofday(&supdate_tm, NULL);
    gettimeofday(&supdate_tm, NULL);
    global_update_nodes=0;
    loop_times=0;
    for(int i = 0; i < dnum; i++){
        //printf("new round.....\n");

        int u = changeEdge[i].first;
        int v = changeEdge[i].second;
        //cout << "insert edge (" << u << ", " << v << ")" << endl; 
        set<Item> H=Hset;
        for(int j=0; j<V; j++)
          Cb[j]=Cb_origin[j];
        set<Item>R= Rset;
        InsertUptTop(u, v, K, H, R);   
        Edge_map[u].erase(v);
        Edge_map[v].erase(u);
        Degree[u]--;
        Degree[v]--;
    }
    gettimeofday(&eupdate_tm, NULL);
    time += double(eupdate_tm.tv_sec - supdate_tm.tv_sec) + double(eupdate_tm.tv_usec - supdate_tm.tv_usec) / 1000000.0f;
    printf("global uptdate nodes =%d, loop times=%d\n", global_update_nodes, loop_times);
    // for(int i = 0; i < V; i++){
    //   cout<<i<<" "<<Cb[i]<<"\n";
    // } 
    return time/dnum;
}

void Graph::init_HtopR(unsigned k){

	set<Item>::iterator it = Hset.begin();
  Rset.clear();
  K=k;

	int cnt = 0;
	while(it != Hset.end()){
		if(cnt == k){
			break;
		}
		cnt++;
		Rset.insert(*it);
		it++;
	}
  //for(int i=0; i<V; i++)
  //  printf("cb[%d]=%lf\n", i, Cb_origin[i]);
	// cout << "hset---------" << endl;
	// OutputRes(Hset);
	// cout << "rset---------" << endl;
	// OutputRes(Rset);
    return;
}

void Graph::init_HtopR(const char* filename, unsigned k){
	
	Cb=new value[V];
  	Cb_origin=new value[V];
  	vis=new bool[V];

	ifstream ifs;
    ifs.open(filename);
    unsigned u;
    float cbu;

    while(ifs>>u>>cbu){

    	Item new_item(u, cbu);
        Hset.insert(new_item);
        Cb[u] = cbu;
    }

	set<Item>::iterator it = Hset.begin();
 	Rset.clear();
  	K=k;

	int cnt = 0;
	while(it != Hset.end()){
		if(cnt == k){
			break;
		}
		cnt++;
		Rset.insert(*it);
		it++;
	}
	for(int i=0; i<V; i++){
    	Cb_origin[i]=Cb[i];
  	}

  //for(int i=0; i<V; i++)
  //  printf("cb[%d]=%lf\n", i, Cb_origin[i]);
	// cout << "hset---------" << endl;
	// OutputRes(Hset);
	// cout << "rset---------" << endl;
	// OutputRes(Rset);
    return;
}

void Graph::
ResetH(int u, int v){
  //function: modify values in Hset, Rset and Cb
  vector<int> wait_to_rest;
  wait_to_rest.push_back(u);
  wait_to_rest.push_back(v);
  bool* nei_vis=(bool*) calloc(V, sizeof(bool));
  for(int nei:Edge_map[u]){//get common neighbor setL
    if(Edge_map[v].find(nei)!=Edge_map[v].end()){
      wait_to_rest.push_back(nei);
      nei_vis[nei]=true;
    }
  }
  for(auto it=Hset.begin(); it!=Hset.end(); it++){
    if(nei_vis[it->node]){
      Hset.erase(it++);
    }
  }
  for(int i=0; i<wait_to_rest.size(); i++){
    value new_cb=compute(wait_to_rest[i]);
    Cb[wait_to_rest[i]]=new_cb;
    Hset.insert(Item(wait_to_rest[i], new_cb));
  }
  Rset.clear();
  int cnt=0;
  for(auto it:Hset){
    if(cnt>=K) break;
    Rset.insert(it);
  }
}

void Graph::LocalUptSmapTop_bef(int u, int v, vector<unsigned>& vec_L){
  if(u<v){
    for(auto l:vec_L){
      SMap[l][mp(u, v)]=0;
    }
  }
  else{
    for(auto l:vec_L){
      SMap[l][mp(v, u)]=0;
    }
  }
  for(auto nei:Edge_map[u]){
    if(v_vis[nei]) continue;
    if(v<nei)
      SMap[u][mp(v,nei)]=0;
    else if(v>nei)
      SMap[u][mp(nei, v)]=0;
  }
  for(auto nei:Edge_map[v]){
    if(u_vis[nei]) continue;
    if(u<nei)
      SMap[v][mp(u, nei)]=0;
    else if(u>nei)
      SMap[v][mp(nei, u)]=0;
  }
}

value Graph::LocalUptSmapTop(bool type, unsigned node, unsigned u, unsigned v, vector<unsigned>& L){//type==0->delete, type==1->insert
  //blabla
  global_update_nodes++;
  if(node==u||node==v){
    int another=(node==u?v:u);
    bool* my_vis_=(node==u?u_vis:v_vis);
    //for(auto nei:Edge_map[node]) my_vis[nei]=true;
    for(auto p:L){
      for(auto x:Edge_map[p]){
        if(my_vis_[x]){//common neighbor of p and node
          if(x!=another && Edge_map[x].find(another)==Edge_map[x].end()){
            if(x<another)
              SMap[node][mp(x, another)]++;
            else
              SMap[node][mp(another, x)]++;
          }
        }
      }
      for(auto q:L){
        if(Edge_map[p].find(q)==Edge_map[p].end()&&p>q){
          for(auto nei:Edge_map[q]){
            if(my_vis_[nei]) double_vis[nei]=true;
          }
          int sum=0;
          for(auto nei:Edge_map[p]){
            if(double_vis[nei]) sum++;
          }
          SMap[node][mp(q,p)]+=sum;
          // for(auto nei:Edge_map[q]){
          //   if(my_vis[nei]) double_vis[nei]=false;
          }
        }
      }
    
    //for(auto nei:Edge_map[node]) my_vis[nei]=false;
    if(type==0){
      for(auto pr:SMap[node]){
        if(pr.first.first==another||pr.first.second==another){
          Cb[node]=Cb[node]-1.0/(pr.second+1);
        }
        else{
          Cb[node]=Cb[node]-1.0/(pr.second+1)+1.0/pr.second;
        }
      }
    }
    else if(type==1){
      for(auto pr :SMap[node]){
        if(pr.first.first==another||pr.first.second==another){
          Cb[node]=Cb[node]+1.0/(pr.second+1);
        }
        else{
          Cb[node]=Cb[node]-1.0/pr.second+1.0/(pr.second+1);
        }
      }
    }
    else{
      printf("error input!\n");
    }

  }
  else{
    for(auto x:Edge_map[node]) my_vis[x]=true;
    for(auto x:Edge_map[node]){
      if(u_vis[x]){
        if(x!=v&&!v_vis[x]){
          int cnt=0;
          for(auto y:Edge_map[x]){
            if(v_vis[y]&&my_vis[y]){
              cnt++;
            }
          }
         // if(node==1028) printf("1<%d> ", y);
          if(x<v)
            SMap[node][mp(x,v)]+=cnt;
          else SMap[node][mp(v,x)]+=cnt;
        }
      }
    }
    for(auto x:Edge_map[node]){
      if(v_vis[x]){
        if(x!=u&&!u_vis[x]){
          int cnt=0;
          for(auto y:Edge_map[x]){
            if(u_vis[y]&&my_vis[y]){
              cnt++;
            }
          }
           //if(node==1028) printf("1<%d> ", y);
          if(x<u)
            SMap[node][mp(x,u)]+=cnt;
          else SMap[node][mp(u,x)]+=cnt;
        }
      }
    }
    int cnt=0;
    for(auto q:L){
      if(my_vis[q]) cnt++;
    }
    if(u<v)
      SMap[node][mp(u,v)]+=cnt;
    else SMap[node][mp(v,u)]+=cnt;
    for(auto x:Edge_map[node]) my_vis[x]=false;

    if(type==0){
      for(auto pr:SMap[node]){
        if((pr.first.first==u&&pr.first.second==v)||(pr.first.first==v&&pr.first.second==u))
          Cb[node]+=1.0/(pr.second+1);
        else
          Cb[node]=Cb[node]-1.0/(pr.second+1)+1.0/pr.second;
      }
    }
    else if(type==1){
      for(auto pr:SMap[node]){
        if((pr.first.first==u&&pr.first.second==v)||(pr.first.first==v&&pr.first.second==u))
          Cb[node]-=1.0/(pr.second+1);
        else
          Cb[node]=Cb[node]-1.0/pr.second+1.0/(pr.second+1);
      }
    }
    else printf("error\n");
  }

  return Cb[node];
}

void Graph::uplocalHR_top(unsigned u, unsigned K, set<Item>& H, set<Item>& R){

	for(auto it = H.begin(); it != H.end(); it++){
      if(it->node==u){
        H.erase(it);break;
      }
    }

    H.insert(Item(u, Cb[u]));
    
    R.clear();
    set<Item>::iterator it = H.begin();
  	int cnt = 0;
  	while(it != H.end()){
  		if(cnt == K){
  			break;
  		}
  		cnt++;
  		R.insert(*it);
  		it++;
  	}
	// cout << "hset---------" << endl;
	// OutputRes(Hset);
	// cout << "rset---------" << endl;
	// OutputRes(Rset);
    return;
}


void Graph::EBWEnumCpt( const char* output, int thread_num=0){
  double exec_time = 0;//计算迭代时间
  exec_time -= get_time();

  if(thread_num==0)
    thread_num=omp_get_num_procs()*2;
  printf("thread_num=%d\n", thread_num);
  omp_set_num_threads(thread_num);
  global_vis.resize(thread_num);
  for(int i=0; i<thread_num; i++){
    global_vis[i]= (bool*) calloc(V, sizeof(bool));
  }

  SMap.resize(V);
  value* Cb=new value[V];
  #pragma omp parallel for
  for(int i=0; i<V; i++){
    Cb[i]=(Degree[i]*(Degree[i]-1))/2;
  }
  ConstructDirected();
  mutex* map_mutex= new mutex[V];

  //printf("start...\n");
  #pragma omp parallel for schedule(dynamic)
  for(int i=0; i<V; i++){
    int thread_id=omp_get_thread_num();
    bool* vis = global_vis[thread_id];
    
    for(auto v: Edge_map_directed[i])
      vis[v]=true;
      
    for(auto v: Edge_map_directed[i]){
      for(auto w:Edge_map_directed[v]){
        if(vis[w]==true){
          UpdateSMap_atomic(i, v, w, map_mutex);
          UpdateSMap_atomic(v, i, w, map_mutex);
          
          map_mutex[w].lock();
          if(i<v) SMap[w][mp(i, v)]=0;
          else SMap[w][mp(v, i)]=0;
          map_mutex[w].unlock();
          
        }
      }
    }
    
    for(auto v: Edge_map_directed[i])
      vis[v]=false;
      
  }

  #pragma omp parallel for schedule(dynamic)
  for(int curnode=0; curnode<V; curnode++){
    for(auto i:SMap[curnode]){
      if(i.first.first>i.first.second) continue;
      Cb[curnode]--;
      //if(curnode==10)
      //printf("pair<%d %d>=%d\n", i.first.first, i.first.second, i.second);
      if(i.second){
        Cb[curnode]+=1.0/(i.second+1);
      }
      
    }
  }
  delete[] map_mutex;

  exec_time+=get_time();
  printf("exec_time for computing =%lf(s)\n", exec_time);
  
  ofstream ofs(output);
  for(int i=0; i<V; i++){
    ofs<<i<<" "<<Cb[i]<<endl;
  }
}

void Graph::UpdateSMap_atomic(unsigned u, unsigned v, unsigned w, mutex* map_mutex){
  vector<tuple<unsigned, unsigned, unsigned> > wait_list;
  for(auto x:Edge_map[u]){
    bool x_v=Edge_map[x].find(v)!=Edge_map[x].end();//true if in E
    bool x_w=Edge_map[x].find(w)!=Edge_map[x].end();
    if(x!=v&&x_v){
      wait_list.push_back(make_tuple(x, v, 0));
    }
    if(x!=w&&x_w){
      wait_list.push_back(make_tuple(x, w, 0));
    }
    if(x==v||x==w) continue;
    if(x!=v&&x_v&&!x_w){
      wait_list.push_back(make_tuple(x, w, 1));
    }
    if(!x_v&&x_w){
      wait_list.push_back(make_tuple(x, v, 1));
      map_mutex[w].lock();
      SMap[w][mp(x,v)]++;
      map_mutex[w].unlock();
    }
  }
  map_mutex[u].lock();
  for(auto p:wait_list){
    if(get<2>(p)==0)
      SMap[u][mp(get<0>(p),get<1>(p))]=0;
    else
      SMap[u][mp(get<0>(p),get<1>(p))]++;
  }
  map_mutex[u].unlock();
}

void Graph::ParallelEdge( const char* output, int thread_num=0){
  double exec_time = 0;//计算迭代时间
  exec_time -= get_time();

  if(thread_num==0)
    thread_num=omp_get_num_procs()*2;
  printf("thread_num=%d\n", thread_num);
  omp_set_num_threads(thread_num);
  global_vis.resize(thread_num);
  for(int i=0; i<thread_num; i++){
    global_vis[i]= (bool*) calloc(V, sizeof(bool));
  }

  auto res=GraphColoring();
  int* color=res.first;
  SMap.resize(V);
  mutex* mu= new mutex[V];
  value* Cb = new value[V];
  #pragma omp parallel for
  for(int i=0; i<V; i++){
    Cb[i]=(Degree[i]*(Degree[i]-1))/2;
  }

  #pragma omp parallel for schedule(dynamic)
  for(int i=0; i<V; i++){
    int thread_id=omp_get_thread_num();
    bool* vis=global_vis[thread_id];
    for(auto j:Edge_map[i]) vis[j]=true;
    vector<unsigned> list; 
    for(auto j:Edge_map[i]){
      list.clear();
      if(j>i) break;//j<i
      for(auto it=Edge_map[j].begin(); it!=Edge_map[j].end(); it++){
        unsigned k=*it;
      //for(auto k:Edge_map[j]){
        if(vis[k]){
          list.push_back(k);
          if(k<i&&k<j){
            mu[i].lock();
            SMap[i][mp(k, j)]=0;
            mu[i].unlock();
            mu[j].lock();
            SMap[j][mp(k, i)]=0;
            mu[j].unlock();
            mu[k].lock();
            SMap[k][mp(j, i)]=0;
            mu[k].unlock();
          }
        }
      }


      for(int p=0; p<list.size(); p++){
        unsigned node_p=list[p];
        for(int q=p+1; q<list.size(); q++){
          unsigned node_q=list[q];
          /*
          if(color[node_p]==color[node_q]){//no edge
            mu[i].lock();
            SMap[i][mp(node_p, node_q)]++;
            mu[i].unlock();
            mu[j].lock();
            SMap[j][mp(node_p, node_q)]++;
            mu[j].unlock();
          }
          */
          //else{
            if(Edge_map[node_p].find(node_q)==Edge_map[node_p].end()){
              mu[i].lock();
              SMap[i][mp(node_p, node_q)]++;
              mu[i].unlock();
              mu[j].lock();
              SMap[j][mp(node_p, node_q)]++;
              mu[j].unlock();
            }
          //}
        }

      }
    }

    for(auto j:Edge_map[i]) vis[j]=false;
  }

  //printf("start computing...\n");
  #pragma omp parallel for schedule(dynamic)
  for(unsigned curnode=0; curnode<V; curnode++){
    //printf("node=%d\n", curnode);
    for(auto i:SMap[curnode]){
      if(i.first.first>i.first.second) {
        //printf("first=%d, second=%d\n", i.first.first, i.first.second);
        continue;
      }
      Cb[curnode]--;
      //#ifdef DEBUG
      //if(curnode==101)
      //printf("pair<%d %d>=%d\n", i.first.first, i.first.second, i.second);
      //#endif
      if(i.second){
        Cb[curnode]+=1.0/(i.second+1);
      }
    }
  }
  exec_time += get_time();
  printf("exec_time for computing =%lf(s)\n", exec_time);

  ofstream ofs(output);
  for(int i=0; i<V; i++)
    ofs<<i<<" "<<Cb[i]<<endl;
}


pair<int*, int> Graph::GraphColoring(){
  
  double exec_time = 0;//计算迭代时间
  exec_time -= get_time();

  int* color= new int[V];
  printf("start coloring...\n");
  bool* cvis=new bool[V];
  int* head = new int[V];
  int* nxt = new int[V];
  int max_degree=0;

  #pragma omp parallel for
  for(int i=0; i<V; i++)  head[i] = V;
  #pragma omp parallel for
  for(int i=0; i<V; i++)  cvis[i] = 0;
    #pragma omp parallel for
  for(int i=0; i<V; i++)  color[i] = V;//V means illegal

  for(int i=0; i<V; i++){
    nxt[i]=head[Degree[i]];
    head[Degree[i]]=i;
    if(Degree[i] > max_degree) max_degree=Degree[i];
  }

  //decide the color of vertices
  int max_color = 0;
  for(int ii=max_degree; ii>=1; ii--){
    for(int jj=head[ii]; jj!=V; jj=nxt[jj]){
      int u = jj;
      for(auto j : Edge_map[u]) {
        int c = color[j];
        if(c != V) {
            cvis[c] = 1;
        }
      }

      for(int j = 0;;j ++){
        if(!cvis[j]) {
          color[u] = j;
          if(j > max_color) max_color = j;
          break;
        }
      }
      for(auto j : Edge_map[u]) {
        int c = color[j];
        if(c != V) cvis[c] = 0;
      }
    }
  }
    max_color++;
    delete[] cvis;
    delete[] head;
    delete[] nxt;
    //compute d_attribute[node]

    printf("max color=%d\n", max_color);
    exec_time += get_time();
    printf("exec_time for coloring =%lf(s)\n", exec_time);
    return make_pair(color,max_color);
}


#endif
#include<bits/stdc++.h>
#include"Graph.hpp"
using namespace std;
int main(int argc, char** argv){

  printf("ego-betweenness\n");
  Graph* graph=new Graph();
  //load graph first
  graph->ReadGraph(argv[1]);
  printf("read success, data = %s\n", argv[1]);

  int ext_cnt = 0;

  if(strcmp(argv[2],"base")==0)
    ext_cnt = graph->BaseBSearch(atoi(argv[3]));
  else if(strcmp(argv[2],"opt")==0){
    if(argc<5) {
      printf("./main graph algorithm k fac\n");
      return 0;
    }
    ext_cnt = graph->OptBSearch(atoi(argv[3]),atof(argv[4]));
  }

  //graph update
  graph->init_HtopR(argv[2], atoi(argv[3]));
  printf("init success, k = %d, hfile = %s\n", atoi(argv[3]), argv[2]);
  graph->random_generate(atoi(argv[4]));
  printf("generate edges success, update edge num = %d\n", atoi(argv[4]));
  double avg_dtime = graph-> rand_delete(atoi(argv[4]), atoi(argv[3]));
  printf("avg_dtime=%lf(s)\n", avg_dtime);

  double avg_itime = graph->rand_insert(atoi(argv[4]), atoi(argv[3]));
  printf("avg_itime=%lf(s)\n", avg_itime);

  double avg_topdtime = graph->rand_topdelete(atoi(argv[4]), atoi(argv[3]));
  printf("avg_topdtime=%lf(s)\n", avg_topdtime);

  double avg_topitime = graph->rand_topinsert(atoi(argv[4]), atoi(argv[3]));
  printf("avg_topitime=%lf(s)\n", avg_topitime);


  //node parallel version
  //graph->EBWEnumCpt(argv[4],atoi(argv[3]));
  //edge parallel version
  //graph->ParallelEdge(argv[4], atoi(argv[3]));

  
  //graph->EBWEnumCpt(1);
  graph->Destroy();
  delete graph;
  return 0;
  
}
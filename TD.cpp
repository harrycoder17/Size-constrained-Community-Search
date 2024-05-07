//author, Haicheng Guo (Email, haicheng_guo@126.com)
// data : 2021-10-8
// function : output  trussness of edge
// argv[0]: run program ;  argv[1] : input file  ; argv[2]: output file
// output file  format
// node number    edge number
// h[i]  h[i+1]
//  vID   trussness   ...  h[i+1] - h[i]


#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<omp.h>
#include<cstring>
#include<unordered_map>
#include<map>



using namespace std;

#define num_of_thread 10

struct G_ds{
    unsigned *h;
    unsigned *e;
    unsigned sz;
    unsigned es;
};

struct Edge{
    unsigned u;
    unsigned v;
};

vector<Edge> ID2ends;



inline char nc(){
    ///*
  static char buf[1000000],*p1=buf,*p2=buf;
  if (p1==p2) { p2=(p1=buf)+fread(buf,1,1000000,stdin); if (p1==p2) return EOF; }
  return *p1++;
    //*/return getchar();
}

inline int read(unsigned &x){
  char c=nc();unsigned b=1;
  if (c=='\377') return 0;
  for (;!(c>='0' && c<='9');c=nc()) if (c=='\377') return 0;   //Ϊɶ����ȥ��
  for (x=0;c>='0' && c<='9';x=x*10+c-'0',c=nc()); x*=b;
  return 1;
}

inline void InFile(G_ds &G,char* path){
    unsigned pre    = 999;//�����name�Ǵ�1��ʼ��
    unsigned ptr_node = 0;
    unsigned ptr_edge = 0;
	FILE* fp=freopen(path,"r",stdin);

    unsigned ver1,ver2;
    read(ver1);read(ver2);

	G.e =(unsigned *)malloc((ver2+5)*sizeof(unsigned ));
    G.h =(unsigned*)malloc((ver1+5)*sizeof(unsigned));

    //unsigned  *ver_to_ID =(unsigned *)malloc(memory_node_id*sizeof(unsigned ));
    unsigned total = ver2;
    while(total --){
       // cout<<total<<endl;
        unsigned ver1 = 0,ver2 = 0;
        int t=read(ver1);
        if (t==0) break;
        t=read(ver2);//t=read(ver3);
        if(pre  != ver1){
            G.h[ptr_node++] = ptr_edge;
            //ver_to_ID[ver1] = ptr_node++;
            G.e[ptr_edge++] = ver2;
            pre = ver1;
	   // cout<<ver1<<" "<<ver2<<endl;
        }else{
            G.e[ptr_edge++] = ver2;
	   // cout<<ver1<<" "<<ver2<<endl;
        }

    }//while

        G.sz = ptr_node;
        G.es = ptr_edge;
        G.h[ptr_node] = ptr_edge;

       /***  name to ID***/
    // #pragma omp parallel for schedule(static) num_threads(num_of_thread)
    // for (int i = 0; i < G.es; i++)  {
    //      G.e[i] = ver_to_ID[G.e[i]];
    // }

//    for(unsigned i=0;i<G.sz;++i){
//            cout<< "i = " << i<<endl;
//        for(unsigned j=G.h[i];j<G.h[i+1];++j){
//            cout<<" "<<G.e[j];
//        }
//        cout<<endl;
//    }

   // free(ver_to_ID);
	fclose(stdin);
}


inline void core_parallel(G_ds &g,vector<int>& deg,unsigned todo,unsigned n) {

    int level = 1;
    while (todo > 0){
    #pragma omp parallel num_threads(num_of_thread)
    {
        //unsigned thread_num = omp_get_num_threads();
        // Size of cache line
        const unsigned BUFFER_SIZE = n / num_of_thread + 1;
        //unsigned buff[BUFFER_SIZE];
        unsigned * buff = new unsigned[BUFFER_SIZE]();
        unsigned index = 0;

        #pragma omp for schedule(static)
        for (long i = 0; i < n; i++) {

            if (deg[i] == level) {
                buff[index] = i;
                index++;
            }
        }

/******************************scan pro **********************/
        unsigned start = 0;
        while(index > start){

            unsigned v = buff[start];
            ++start;
            //Check the adj list of vertex v

            for (unsigned i=g.h[v];i<g.h[v+1];i++)
            {

               unsigned u =  g.e[i];
                if (deg[u] > level) {
                    int du =  __sync_fetch_and_sub(&deg[u], 1);
                    if(du <= level )  __sync_fetch_and_add(&deg[u],1);
                    if (du == (level+1)) {
                        buff[index] = u;
                        index++;
                    }
                }
            }//for-end
        }//while-end
        #pragma omp barrier
        __sync_fetch_and_sub(&todo,index);
        delete [] buff;
    }//  #pragma omp parallel
    ++level;
    }//while-end
    return;
}

//ParK to compute k core decompsitions in parallel
inline void ParK(G_ds &g, vector<int> &deg){


    unsigned todo = g.sz;
    unsigned n = todo;
   #pragma omp parallel for schedule(static) num_threads(num_of_thread)
    for (long i = 0; i < n; i++)  {
        //deg[i] = g[i].size();
        deg[i] = g.h[i+1]-g.h[i];
    }
    core_parallel(g,deg,todo,n);

//    for(unsigned i=0;i<g.sz;++i){
//        cout<<deg[i]<<" ";
//    }
}

inline void calculateSupport_parallel(G_ds& G,vector<long> &ID2sup,unsigned * edge2id,vector<unordered_map<int,int>>& edgeid) {//calculate the support number on edges and return the max sup
    unsigned ID=0;


    //#pragma omp parallel for
    for(unsigned i=0;i<G.sz;i++){//iterate through vertices, i is the id of the first vertex

        unordered_map<int,int> temp;
	    //for(NeighborsWithHeads ::iterator it=a.begin();it!=a.end();it++){//iterate through neighbors of vi
        for( unsigned j=G.h[i];j<G.h[i+1];++j){

            unsigned it_first = G.e[j];
            if(it_first <= i){
                continue;
            }

            Edge tempEdge;
            tempEdge.u=i;
            tempEdge.v= it_first;

            temp.insert(pair<int,int>(it_first,ID));
            ID2ends.emplace_back(tempEdge);


            edge2id[j] = ID++;
        }

        edgeid.emplace_back(temp);
    }

    ID2sup.resize(ID);
#pragma omp parallel num_threads(num_of_thread)
    {

    bool * thash = new bool[G.sz]();
    #pragma omp for schedule(static)
	for (unsigned i=0;i<ID;i++){//iterate through vertices, i is the id of the first vertex

        unsigned a = ID2ends[i].u; //nbs of the first vertex
        unsigned b = ID2ends[i].v;
        unsigned a_size = G.h[a+1]-G.h[a];
        unsigned b_size = G.h[b+1]-G.h[b];



        unsigned ptr;
        unsigned ptr_hash;

        if(a_size > b_size){
            for(unsigned j = G.h[a];j<G.h[a+1];++j){
                thash[G.e[j]] = true;
            }
            ptr= b;
            ptr_hash = a;
        }
        else{
            for(unsigned j = G.h[b];j<G.h[b+1];++j){
                thash[G.e[j]] = true;
            }
            ptr_hash = b;
            ptr = a;
        }

        for(unsigned j = G.h[ptr];j<G.h[ptr+1];++j){
                if(thash[G.e[j]]) {  ++ID2sup[i] ;}
        }

        for(unsigned j = G.h[ptr_hash];j<G.h[ptr_hash+1];++j){
                thash[G.e[j]] = false ;
        }

	}//for
	      delete [] thash;
    }//pragma end
}


inline void SortSwap(unsigned i,unsigned ID,vector<unsigned>& posi, vector<unsigned>& ver,vector<long>& ID2sup, vector<unsigned> & bin){
    
    if(ID2sup[ID]>i){
        unsigned du = ID2sup[ID];
        unsigned pu = posi[ID];
        unsigned pm = bin[du];
        unsigned m  = ver[pm];

        if(ID!=m){
            posi[ID] = pm;
            posi[m]  = pu;
            ver[pu]  = m;
            ver[pm]  = ID;
        }
        ++bin[du];
        --ID2sup[ID];
    }else {//放到bin[1]
        unsigned du = ID2sup[ID];
        unsigned pu = posi[ID];
        unsigned pm = bin[du];
        unsigned m  = ver[pm];

        posi[ID] = pm;
        posi[m]  = pu;
        ver[pu]  = m;
        ver[pm]  = ID;
    }

    return;
}



//void TrussDecomp(vector<khash_t(32)*> Sub, khash_t(32)* G2Sid, vector<unsigned> &ID2sup, vector<vector<unsigned>>& ID2end){

void k_truss(G_ds G,vector<long>& ID2sup,unsigned &result0, unsigned &result1,vector<unordered_map<int,int> >edgeid){  

    unsigned num_sup = *max_element(ID2sup.begin(),ID2sup.end());
    unsigned num_ver = ID2sup.size();
    vector<bool>  flag(num_ver,false);

    vector<unsigned> bin(num_sup+1,0);    //下标是度，值为个数/起点
    for(unsigned i=0;i<num_ver;++i){
        ++bin[ID2sup[i]];
    }

    unsigned start = 0; //起始位置下标为0
    for(unsigned i=0;i<=num_sup;++i){ // 经函数，值从个数变为度起点(初始起点是0)
        unsigned num = bin[i];
        bin[i] = start;
        start += num;
    }

    vector<unsigned>posi(num_ver);//下标是顶点ID，值是按度排序位置
    vector<unsigned>ver(num_ver);//下标是按度排序位置，值是顶点ID
    for(unsigned i=0;i<num_ver;++i){
        posi[i] = bin[ID2sup[i]];
        ver[posi[i]] = i;
        ++bin[ID2sup[i]];
    }
    for(unsigned i=num_sup;i>0;--i){//bin恢复，值为度的起点
        bin[i] = bin[i-1];
    }
    bin[0] = 0;

    bool * thash = new bool[G.sz]();
//////////////////////////////////////////////////////////////////////////////////////////////////////
    for(unsigned i=0;i<=num_sup;++i){       //sup >= 0
        unsigned end_bin =  i!=num_sup?bin[i+1]:num_ver;
        unsigned sum = 0;

        while(bin[i]<end_bin){
            unsigned cur_id = ver[bin[i]];
            ++bin[i];
            ++sum;
            flag[cur_id] = true;
            
            // teminal condition
            if(bin[i] == num_ver ){

                bin.resize(0);
                ver.resize(0);
                flag.resize(0);
                posi.resize(0);
                result0 = *max_element(ID2sup.begin(),ID2sup.end()) ;
                result1 = 0;
                for(unsigned val:ID2sup) {
                    if(val == result0) ++result1;
                }
                delete [] thash;
                result0 += 2;
                //ResultSub(Sub,G2Sid,ID2sup);
                return;
            }//end

            unsigned u = ID2ends[cur_id].u;
            unsigned v = ID2ends[cur_id].v;

            unsigned g ;  //greater size
            unsigned l ;  //less    size
            
            if(G.h[u+1]-G.h[u] > G.h[v+1]-G.h[v]){
                g = u;
                l = v;
            }else{
                g = v;
                l = u;
            }
            for(unsigned m = G.h[g];m<G.h[g+1];++m){
                    thash[G.e[m]] = true;
            }
            
            for(unsigned k=G.h[l];k<G.h[l+1];++k){
                
                unsigned  w = G.e[k];
                if(!thash[w]) continue;


                //(u,w)                           
                unsigned ID;
                if(u<w)  ID = edgeid[u][w];
                else     ID = edgeid[w][u];

                unsigned ID1;
                if(v<w)  ID1 = edgeid[v][w];
                else     ID1 = edgeid[w][v];

                if(flag[ID]||flag[ID1])  continue;

                //(u,w) (v,w)
                SortSwap(i,ID,posi,ver,ID2sup,bin);
                SortSwap(i,ID1,posi,ver,ID2sup,bin);

            }//for-end
            end_bin =  i!=num_sup?bin[i+1]:num_ver;

            for(unsigned m = G.h[g];m<G.h[g+1];++m){
                    thash[G.e[m]] = false;
            }
        }//while-end

       //cout<<endl;
    }//for-end

    return;

}

void output(G_ds G,   vector<unordered_map<int,int>>  edgeid , vector<long> ID2sup ,char* path){

    ofstream out( path);

    out<< G.sz<<" "<<G.es<<endl;

    for(unsigned i = 0; i< G.sz;++i){

        unsigned beginn = G.h[i];
        unsigned endd = G.h[i+1];

       out<<beginn<<" "<<endd<<endl;
        //int num = beginn - endd;
        multimap<unsigned,unsigned, greater<unsigned>> temp;
        for(unsigned j=  beginn ;j< endd ; ++j){

            unsigned id = G.e[j];

            unsigned te1 = i <id ? i: id;
            unsigned te2 = i<id ?  id:i;

            unsigned truss = ID2sup[edgeid[te1].find(te2)->second] + 2;

            out<<id<<" "<<truss<<endl;

        }


    }

   out.close();

}

int main(int argc, char** argv){

double time0 =omp_get_wtime() ;

    G_ds G;
    InFile(G,argv[1]);

double time1 =omp_get_wtime() ;
cout<<"read time"<< time1 - time0 <<endl;


    vector<unordered_map<int,int>> edgeid;
    vector<long> ID2sup;

    unsigned *edge2id = (unsigned *)malloc( G.es * sizeof(unsigned) );
    calculateSupport_parallel(G,ID2sup,edge2id,edgeid);/////////////////////


    unsigned result0,result1;
    k_truss(G ,ID2sup,result0, result1,edgeid);     ////////////                                      /////////////////
    

   output(G,edgeid,ID2sup,argv[2]);
   cout<<" output  success"<<endl;
return 0;

}

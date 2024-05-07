//author, Haicheng Guo (Email, haicheng_guo@126.com)
//data, 2022-4-10
//function, return the subgraph satisfying the constraint of STCS problem
//parameter, argv[1] input file (e.g., truss.txt); argv[2] query vertex; argv[3] size lower limit; argv[4] size upper limit; argv[6] runnig time limit; 

/*input file format
total node number    total edge number
h[i]  h[i+1]
v_id truss  { total pairs h[i+1] - h[i]}
note that you can modify the format based on your situlation and the graph struct classs "G_ds" in the following.
*/


#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<cstring>
#include<unordered_map>
#include<map>
#include<unordered_set>
#include<stdlib.h>
#include<cstring>
#include<omp.h>
#include<malloc.h>
#include<set>
#include<list>
#include<queue>
#include "khash.h"

using namespace std;

// gragh struct
struct G_ds{
    unsigned *h;//hash from u_id to the address of the first neighbor of u_id in arry e
    unsigned *e;//store all nighbors of vertex u
    unsigned *t;//hash from edge_id to its trussness
    unsigned hs;//node number
    unsigned es;//edge number
};
G_ds G;

// domination struct
struct DOM{
    vector<unsigned> toge;   // domination by eath other
    vector<unsigned> withdraw;
    unsigned tru;            // kl
};

struct NEIDEG{
    unsigned val;
    unsigned deg;
    unsigned tru;
};

// kash 
KHASH_MAP_INIT_INT(32,unsigned)
khiter_t k;

// global variable
vector<khash_t(32)*> khash;

unordered_set<unsigned> NeiDeg; //store the 1-hop of query
vector<NEIDEG> key;



unsigned query;  // query vertex
unsigned ku;     // upper bound of optimal k
unsigned kl;     // lower bound of optimal k
unsigned size_l; // lower bound of size  
unsigned size_u; // upper bound of size 
unsigned trq;


// record conhesive score
double r_t = 0;
double r_a = 0;
double r_d = 0;
double r_c = 0;
double r_s = 0;


char * path_out;

vector<unsigned> C;      // the partial solution set for heuristic algorihtm
vector<unsigned> R;      // the candidate set for heuristic algorihtm
map<unsigned,unsigned> tr;
//unordered_set<unsigned> CC;
map<unsigned,unsigned> CC;   // the partial solution set for exact algorihtm
map<unsigned,unsigned> Result;
unordered_set<unsigned> Cstar;  //stores vertex whose trussness is no bigger kl in C
map<unsigned,unsigned> RR;  // the candidate set for exact algorihtm
unordered_map<unsigned,unsigned> RRH;
vector<unsigned> VerSup;   // store the vertex support 
//unordered_map<unsigned,unsigned> RR_LB;
map<unsigned,DOM> Dom;

bool IsLess = false;
bool IsSuc = false;
bool OverTime = false;
unordered_set<unsigned> H;  // store the best community currently
vector<unsigned> H_cur; //store kl 


///////////////exact
double pro_start_time   = 0;
double exact_start_time = 0; 
double maxtime          = 1000;   // runtime limit

////////////test
double sum_suc = 0;
double sum_le  = 0;
/************  input ********/
inline char nc(){
    static char buf[1000000],*p1=buf,*p2=buf;
    if (p1==p2) { p2=(p1=buf)+fread(buf,1,1000000,stdin); if (p1==p2) return EOF; }
    return *p1++;
    //*/return getchar();
}
inline int read(unsigned &x){
    char c=nc();unsigned b=1;
    if (c=='\377') return 0;
    for (;!(c>='0' && c<='9');c=nc()) if (c=='\377') return 0;   
    for (x=0;c>='0' && c<='9';x=x*10+c-'0',c=nc()); x*=b;
    return 1;
}
void Infile(char * Inpath){
    
    double time_st = omp_get_wtime();

    FILE* fp=freopen(Inpath,"r",stdin);


    unsigned ver1 = 0,ver2 = 0;
    int t=read(ver1);//node number  
    t=read(ver2);// edge nunmber

    G.h =(unsigned*)malloc((ver1+1)*sizeof(unsigned)); // h[ver1+1] store address of end
    G.e =(unsigned *)malloc((ver2)*sizeof(unsigned ));
    G.t =(unsigned *)malloc((ver2)*sizeof(unsigned ));

    G.hs = ver1;
    G.es = ver2;
    G.h[ver1] = ver2;

    unsigned i = 0;

    while(1){

        t=read(ver1);
        if (t==0) break;
        t=read(ver2);

        G.h[i++] = ver1;

        unsigned beginn = ver1;
        unsigned endd   = ver2;

        khash_t(32) * temp = kh_init(32);
        int temp_ret;

        for(unsigned j = beginn;j< endd;++j){
            t=read(ver1);
            t=read(ver2);

            G.e[j] = ver1;
            G.t[j] = ver2;

            k = kh_put(32,temp,ver1,&temp_ret);
            kh_value(temp,k) = ver2;

        }

        khash.emplace_back(temp);
        //kh_destroy(32,temp);

    }//while


     cout<<Inpath<<" node number: "<<G.hs<<" edge number:"<<G.es<<endl;


	fclose(stdin);

    cout<<"read time "<<omp_get_wtime() - time_st<<endl;
    

}


// return vertex trussness of vertex
inline unsigned FindT(unsigned ver){
    if(tr.find(ver)!=tr.end()) return tr.find(ver)->second;
    unsigned t = 0;
    for(unsigned i=G.h[ver];i<G.h[ver+1];++i){
        t = G.t[i] > t ? G.t[i]:t;
    }
    tr.insert(make_pair(ver,t));
    return t;
}


inline void QueueClear(queue<unsigned>& q) {
    queue<unsigned> emp;
    swap(emp, q);
}
inline unsigned EdgeIDinSub(unsigned u,unsigned w,vector<khash_t(32)*> Sub,khash_t(32)* G2Sid){

    unsigned u2Sub = kh_value(G2Sid,kh_get(32,G2Sid,u)); 
    unsigned w2Sub = kh_value(G2Sid,kh_get(32,G2Sid,w));

    if(u<w)  return kh_value(Sub[u2Sub],kh_get(32,Sub[u2Sub],w));
    else     return kh_value(Sub[w2Sub],kh_get(32,Sub[w2Sub],u));
}

/************  triangle functions ********/
void TriUpdate(unsigned id, vector<khash_t(32)*>& Sub,khash_t(32)* G2Sid, vector<unsigned> &ID2sup, vector<vector<unsigned>>& ID2end){

    unordered_set<unsigned> SupUpdate;  //store the id of edge needed updating sup 

    khash_t(32) * temp = kh_init(32);
    Sub.emplace_back(temp);
    //kh_destroy(32,temp);

    int G2Sid_ret;
    khiter_t k;

    // Update G --> S id inserting
    k = kh_put(32,G2Sid,id,&G2Sid_ret);
    kh_value(G2Sid,k) = Sub.size()-1;

    if(id == query) return;

    for(unsigned i=G.h[id];i<G.h[id+1];++i){
        
        unsigned v = G.e[i];

        if(kh_get(32,G2Sid,v) == kh_end(G2Sid)) continue;
        unsigned ptr_l, ptr_s;
        if(G.h[id+1] - G.h[id] > G.h[v+1] - G.h[v]){
            ptr_l = id;
            ptr_s = v;
        }
        else{
            ptr_l = v;
            ptr_s = id;
        }

        //Update new edge as id inserted 
        // id v
        vector<unsigned> temp(2);
        temp[0] = id; temp[1] = v;
        ID2end.emplace_back(temp);
        temp.resize(0);

        unsigned te1, te2;
        if(v<id){
            te1 = kh_value(G2Sid,kh_get(32,G2Sid,v));
            te2 = id;
        }
        else{
            te1 = kh_value(G2Sid,kh_get(32,G2Sid,id));
            te2 = v;
        }

        k = kh_put(32,Sub[te1],te2,&G2Sid_ret);
        kh_value(Sub[te1],k) = ID2end.size() - 1;

        unsigned sum_tri = 0;

        for(unsigned j = G.h[ptr_s];j<G.h[ptr_s+1];++j){
            unsigned w = G.e[j];
            if(kh_get(32,G2Sid,w) == kh_end(G2Sid)) continue;

            if(kh_get(32,khash[ptr_l],w) != kh_end(khash[ptr_l])){          
                // v w
                SupUpdate.insert(EdgeIDinSub(v,w,Sub,G2Sid));

                ++ sum_tri;
            }
        }

        ID2sup.emplace_back(sum_tri);
    }
    // Update the sup of edge existing in sub
    for(unsigned val:SupUpdate){
        ++ID2sup[val];
    }

    return;
}

// triangle updating in C
void TriDele(unsigned id, vector<khash_t(32)*>& Sub, khash_t(32)* G2Sid, vector<unsigned> &ID2sup, vector<vector<unsigned>>& ID2end){
   

    unordered_set<unsigned> SupUpex;
    unsigned sum_edge_del = 0;

    for(unsigned i=G.h[id];i<G.h[id+1];++i){
        
        unsigned v = G.e[i];

        if(kh_get(32,G2Sid,v) == kh_end(G2Sid)) continue;

        // id v
        ++sum_edge_del;   // delete edge(id,v)

        unsigned te1, te2;
        if(v<id){
            te1 = kh_value(G2Sid,kh_get(32,G2Sid,v));
            te2 = id;
            kh_del(32,Sub[te1],kh_get(32,Sub[te1],te2));
        }
        else{
            //....
        }
        //delete (id,v) in Sub
        
        //prepare for calculating triangle
        unsigned ptr_l, ptr_s;
        if(G.h[id+1] - G.h[id] > G.h[v+1] - G.h[v]){
            ptr_l = id;
            ptr_s = v;
        }
        else{
            ptr_l = v;
            ptr_s = id;
        }

        for(unsigned j = G.h[ptr_s];j<G.h[ptr_s+1];++j){
            unsigned w = G.e[j];
            if(kh_get(32,G2Sid,w) == kh_end(G2Sid)) continue;


            if(kh_get(32,khash[ptr_l],w) != kh_end(khash[ptr_l])){
                
                // v w
                unsigned te1, te2;
                if(v<w){
                    te1 = kh_value(G2Sid,kh_get(32,G2Sid,v));
                    te2 = w;
                }
                else{
                    te1 = kh_value(G2Sid,kh_get(32,G2Sid,w));
                    te2 = v;
                }
                SupUpex.insert(kh_value(Sub[te1],kh_get(32,Sub[te1],te2)));
            }
        }

    }

    while(sum_edge_del--){
        ID2end.pop_back();
        ID2sup.pop_back();
    }
    for(unsigned val:SupUpex){
        --ID2sup[val];
    }

    kh_del(32,G2Sid,kh_get(32,G2Sid,id));
    Sub.pop_back();

    return;
}

inline void TriSubBuild(vector<unsigned>& C,khash_t(32)* G2Sid,vector<khash_t(32)*>& Sub,unsigned short* State,
 vector<unsigned> & ID2sup,vector<vector<unsigned>> & ID2end){
    
    unsigned id_sup = 0;
    unsigned id_ver = 0;

    int G2Sid_ret;
    int cur_ret;
    khint32_t k;

    for(unsigned val:C){

        //new id in Sub for triangel computing 
        k = kh_put(32,G2Sid,val,&G2Sid_ret);
        kh_value(G2Sid,k) = id_ver++;

        // store the neighbor of the val 
        khash_t(32) * cur = kh_init(32);

        for(unsigned i=G.h[val];i<G.h[val+1];++i){
            unsigned v = G.e[i];
            if(v > val && State[v] == 1){

                vector<unsigned> temp(2);
                temp[0] = val; temp[1] = v;
                ID2end.emplace_back(temp);
                temp.resize(0);

                k = kh_put(32,cur,v,&cur_ret);
                kh_value(cur,k) = id_sup++;

            }
        }

        Sub.emplace_back(cur);
    }

    //store the value of support
    ID2sup.resize(id_sup);

}

// triangle count for heu
void TriCalculate(vector<unsigned> & C,unsigned short *State, vector<khash_t(32)*>& Sub, khash_t(32)* G2Sid,
 vector<unsigned> & ID2sup,vector<vector<unsigned>> & ID2end){
    
    TriSubBuild(C,G2Sid,Sub,State,ID2sup,ID2end);

    unsigned id_sup = 0;
    for(auto val:C){
    
        for(unsigned i=G.h[val];i<G.h[val+1];++i){
            unsigned v = G.e[i];
        
            if(v > val && State[v] == 1){
           
                unsigned g,l; //greater less
                if(G.h[val+1]-G.h[val] > G.h[v+1]-G.h[v]){
                    g = val;
                    l = v;    
                }
                else{
                    g = v;
                    l = val;
                }

                //find common vertex
                unsigned sum_sup = 0;
                for(unsigned j=G.h[l];j<G.h[l+1];++j){
                    unsigned w = G.e[j];
                    if(kh_get(32,G2Sid,w)==kh_end(G2Sid)) continue;
                    if(kh_get(32,khash[g],w)==kh_end(khash[g])) continue;
                    ++sum_sup;
                }
                ID2sup[id_sup++] = sum_sup;
             
            }

        }
        
    }

    return;
}


/************  truss decomption ********/
inline bool ExactSub(unsigned sup,vector<khash_t(32)*> Sub, khash_t(32)* G2Sid,vector<unsigned> ID2sup){

//     cout<<"Exact Sub"<<sup<<endl;

    unsigned sum_H = 1;
    queue<unsigned> Q;
    Q.push(query);
    bool * inqueue = new bool [G.hs]{};
    inqueue[query] = true;

    while(!Q.empty()){
        unsigned u = Q.front();
        Q.pop();

        for(unsigned i=G.h[u];i<G.h[u+1];++i){  

            unsigned w = G.e[i];

            if(kh_get(32,G2Sid,w)==kh_end(G2Sid)) continue;

            if(VerSup[kh_value(G2Sid,kh_get(32,G2Sid,w))] >= sup && !inqueue[w]) {
                Q.push(w);inqueue[w]=true;++sum_H;
            }
            
            // if(ID2sup[EdgeIDinSub(u,w,Sub,G2Sid)] >= sup && !inqueue[w]) {Q.push(w);inqueue[w]=true;++sum_H;}
        }
    }

    if(sum_H >= size_l){
        H.clear();
        for(unsigned i=0;i<G.hs;++i){
            if(inqueue[i]) {H.insert(i);}
        }

        delete[] inqueue;
        return true;
    }

    delete [] inqueue;
    return false;
 
}
// true --> update H
inline bool ResultSub(vector<khash_t(32)*> Sub, khash_t(32)* G2Sid, vector<unsigned> &ID2sup){

    VerSup.clear();
    VerSup.resize(kh_size(G2Sid));   // hash  origal id --->  its vertex trussness
 
    // compute vertex support of every vertex   new update
    
    for (unsigned val :C){
        
        unsigned versup = 0;
        for(unsigned i=G.h[val];i<G.h[val+1];++i){
            unsigned v = G.e[i];
            if(kh_get(32,G2Sid,v)==kh_end(G2Sid)) continue;

            unsigned ID = EdgeIDinSub(val,v,Sub,G2Sid);

            versup = versup > ID2sup[ID] ? versup : ID2sup[ID];
        }

        VerSup[kh_value(G2Sid,kh_get(32,G2Sid,val))] = versup;

    }


    // find the maximal trussness including query
    unsigned up_cur = VerSup[kh_value(G2Sid,kh_get(32,G2Sid,query))];

    //computing the size of up_cur-truss subgraph in Sub
    

    up_cur = up_cur < ku-2?up_cur:ku-2;  // not that support

    //fast check 

    for(int i = up_cur; i>=0; --i ){
        // if(kl > 1 && i <= kl -2) break;
        if(ExactSub(i,Sub,G2Sid,ID2sup)) {kl=i+2; return true;} 
    }


    return false;
}

inline bool ResultSubExa(vector<khash_t(32)*> Sub, khash_t(32)* G2Sid, vector<unsigned> &ID2sup){

    VerSup.clear();
    VerSup.resize(kh_size(G2Sid));   // hash  origal id --->  its vertex trussness
 
    // compute vertex support of every vertex   new update
    
    for (map<unsigned,unsigned>:: iterator it = CC.begin();it != CC.end();++it){
        
        unsigned versup = 0;
        unsigned val = it->first;
        for(unsigned i=G.h[val];i<G.h[val+1];++i){
            unsigned v = G.e[i];
            if(kh_get(32,G2Sid,v)==kh_end(G2Sid)) continue;

            unsigned ID = EdgeIDinSub(val,v,Sub,G2Sid);

            versup = versup > ID2sup[ID] ? versup : ID2sup[ID];
        }

        VerSup[kh_value(G2Sid,kh_get(32,G2Sid,val))] = versup;

    }


    // find the maximal trussness including query
    unsigned up_cur = VerSup[kh_value(G2Sid,kh_get(32,G2Sid,query))];

    //computing the size of up_cur-truss subgraph in Sub
    

    up_cur = up_cur < ku-2?up_cur:ku-2;  // not that support

    for(int i = up_cur; i>=0; --i ){
        // if(kl > 1 && i <= kl -2) break;
        if(ExactSub(i,Sub,G2Sid,ID2sup)) {kl=i+2; return true;} 
    }


    return false;
}

// sub-operation of binsort of truss decomp
inline void SortSwap(unsigned i,unsigned ID,vector<unsigned>& posi, vector<unsigned>& ver,vector<unsigned>& ID2sup, vector<unsigned> & bin){
    
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
    }else {//bin[1]
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
void TrussDecomp(vector<khash_t(32)*> Sub, khash_t(32)* G2Sid, vector<unsigned> &ID2supp, vector<vector<unsigned>>& ID2end){

    vector<unsigned> ID2sup(ID2supp);

    unsigned num_sup = *max_element(ID2sup.begin(),ID2sup.end());
    unsigned num_ver = ID2sup.size();
    vector<bool>  flag(num_ver,false);

    // binsort in the following
    vector<unsigned> bin(num_sup+1,0);    

    for(unsigned i=0;i<num_ver;++i){
        ++bin[ID2sup[i]];
    }


    unsigned start = 0; 
    for(unsigned i=0;i<=num_sup;++i){ 
        unsigned num = bin[i];
        bin[i] = start;
        start += num;
    }

    vector<unsigned>posi(num_ver);
    vector<unsigned>ver(num_ver);
    for(unsigned i=0;i<num_ver;++i){
        posi[i] = bin[ID2sup[i]];
        ver[posi[i]] = i;
        ++bin[ID2sup[i]];
    }
    for(unsigned i=num_sup;i>0;--i){
        bin[i] = bin[i-1];
    }
    bin[0] = 0;


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
            if(bin[i] == num_ver || ID2sup[cur_id]>=ku-2){

                bin.resize(0);
                ver.resize(0);
                flag.resize(0);
                posi.resize(0);
                // cout<<"result extra begin"<<endl;

                ResultSub(Sub,G2Sid,ID2sup);
                //cout<<"result extra end"<<endl;
                return;
            }//end

            unsigned u = ID2end[cur_id][0];
            unsigned v = ID2end[cur_id][1];

            unsigned u2Sub = kh_value(G2Sid,kh_get(32,G2Sid,u));
            unsigned v2Sub = kh_value(G2Sid,kh_get(32,G2Sid,v));

            unsigned g ;  //greater size
            unsigned l ;  //less    size
            
            if(G.h[u+1]-G.h[u] > G.h[v+1]-G.h[v]){
                g = u;
                l = v;
            }else{
                g = v;
                l = u;
            }
            
            for(unsigned k=G.h[l];k<G.h[l+1];++k){
                
                unsigned  w = G.e[k];
                if(kh_get(32,G2Sid,w) == kh_end(G2Sid) || kh_get(32,khash[g],w) == kh_end(khash[g])) continue;

                unsigned w2Sub = kh_value(G2Sid,kh_get(32,G2Sid,w));

                //(u,w)                           
                unsigned ID;
                if(u<w)  ID = kh_value(Sub[u2Sub],kh_get(32,Sub[u2Sub],w));
                else     ID = kh_value(Sub[w2Sub],kh_get(32,Sub[w2Sub],u));

                unsigned ID1;
                if(v<w)  ID1 = kh_value(Sub[v2Sub],kh_get(32,Sub[v2Sub],w));
                else     ID1 = kh_value(Sub[w2Sub],kh_get(32,Sub[w2Sub],v));

                if(flag[ID]||flag[ID1])  continue;

                //(u,w) (v,w)
                SortSwap(i,ID,posi,ver,ID2sup,bin);
                SortSwap(i,ID1,posi,ver,ID2sup,bin);

            }//for-end
            end_bin =  i!=num_sup?bin[i+1]:num_ver;
        }//while-end

       //cout<<endl;
    }//for-end

    return;

}


/************  find ku ********/
// return conditions: 0, [l,h]; 1. >h; >=2 = kl;
inline int ExtraK(unsigned ku_cur){

    unsigned sum_H = 1;
    queue<unsigned> Q;
    Q.push(query);
    bool * inqueue = new bool [G.hs]{};
    inqueue[query] = true;

    while(!Q.empty()){
        unsigned u = Q.front();
        Q.pop();

        for(unsigned i=G.h[u]; i<G.h[u+1];++i){
            unsigned v = G.e[i];
            //cout<<v<<" v-truss: "<<kh_value(khash[u],kh_get(32,khash[u],v))<<endl;
            if(FindT(v)>=ku_cur){  // new update
                if(!inqueue[v]){
                    Q.push(v);
                    ++sum_H;
                    inqueue[v] = true;
                    if(sum_H > size_u) break;
                }
            }
            //else break;
        } 
        if(sum_H > size_u) break;  
    }

    if(sum_H >= size_l && sum_H != 1){
        if(sum_H <= size_u){  //[l,h]
 
            H.clear();
            for(unsigned i=0;i<G.hs;++i){if(inqueue[i]) H.insert(i);}
            ku = ku_cur;
            delete [] inqueue;
            return 0;
        }
        else{               //>h
            ku = ku_cur;
            delete [] inqueue;
            return 1;
        }
    }

    //then check ku =? tru(query) - x;

    unsigned tru = 0;
    H_cur.clear();
    for(unsigned i=0;i<G.hs;++i){
        if(inqueue[i]){
            for(unsigned j= G.h[i];j<G.h[i+1];++j){
                unsigned tru_v = FindT(G.e[j]);         // new update
                if(!inqueue[G.e[j]] && tru_v < ku_cur && tru_v > tru) tru = tru_v;
            }
            H_cur.emplace_back(i);
        }
    }
    delete [] inqueue;

    if(tru == 0) return -1;
    return tru;
}
inline int Findku(){

    unsigned flag = trq > size_u ?size_u:trq;
    while(1){
        flag = ExtraK(flag);
        if(flag == 0 ) return 0;
        if(flag == 1 ) return 1;

        if(flag == -1) {return -1;}
    }  
}


/************  find clique ********/
inline bool comparison(NEIDEG &a, NEIDEG &b){// bigger
        return a.deg > b.deg;
}

inline void FindCliRe (vector<unsigned> & Cli){

    //unordered_set<unsigned> NeiDeg;
    for(unsigned i=G.h[query];i<G.h[query+1];++i){
        if(G.t[i] < ku) continue;
        NeiDeg.insert(G.e[i]); 
    }    

    for( unordered_set<unsigned>::iterator it=NeiDeg.begin();it!=NeiDeg.end();++it){
        unsigned sum_cur = 0;
        for(unsigned i=G.h[*it];i<G.h[*it+1];++i){
            if(NeiDeg.find(G.e[i])!=NeiDeg.end()) ++ sum_cur;
            //if(G.t[i]>=ku) ++ sum_cur;
        } 
        NEIDEG temp;
        temp.val = *it; temp.deg = sum_cur; temp.tru = FindT(*it);//G.t[G.h[*it]];
        key.emplace_back(temp);
    }

    sort(key.begin(),key.end(),comparison);

    Cli.emplace_back(query);

    for(NEIDEG val:key){
        bool flag = true;
        unsigned w = val.val;
    
        for(unsigned i=1;i<Cli.size();++i){
            if(kh_get(32,khash[Cli[i]],w) == kh_end(khash[Cli[i]])) {flag = false;break;}
        }
        if(flag) Cli.emplace_back(w);
        if(Cli.size() == size_u) break;
    }
    
    return;

}

inline void FindCli(vector<unsigned> & Cli){


    // first check the neighbor of query
    unsigned sum = 1;
    unsigned ku_or = trq;//G.t[G.h[query]];
    for(unsigned i=G.h[query];i<G.h[query+1];++i){
        if(G.t[i] == ku_or) {++sum;}
    }

    if(sum == ku_or){   //note that sum  = 1
        H.clear();
        H.insert(query); 
        sum = 1;

        for(unsigned i=G.h[query];i<G.h[query+1];++i){
            if(G.t[i]==ku_or) {Cli.emplace_back(G.e[i]);}
            if(++sum >= ku) break;
        }  
        Cli.emplace_back(query);        
        return;   
    }


    // then general
    FindCliRe(Cli);
    return;
}


unsigned ScoreHeu(){

    int score_best = -1;
    unsigned ver_best = 0;
    unsigned truss_best = 0;
   
    for(unordered_map<unsigned,unsigned>::iterator it=RRH.begin();it!=RRH.end();++it){

        int score_cur = it->second;
        unsigned ver_cur = it->first;
        unsigned truss_cur = FindT(ver_cur);//G.t[G.h[ver_cur]];

        if(score_best < score_cur ){
            ver_best = ver_cur;
            truss_best = truss_cur; 
            score_best = score_cur;
        }
        else if(score_best ==  score_cur && truss_best < truss_cur){
            ver_best = ver_cur;
            truss_best = truss_cur; 
            score_best = score_cur;

        }
 
    }
    return ver_best;
}

inline void RUpdate(unsigned val,unsigned short * State){
   
    // update R due to the vertex is selected form R
    RRH.erase(val);     

   // update R due to C grows
    State[val] = 1;
    for(unsigned i = G.h[val]; i<G.h[val+1];++i){
        unsigned v = G.e[i];
        if(State[v]!=1 && FindT(v)>=ku) {   // new update
            unordered_map<unsigned,unsigned>::iterator it = RRH.find(v);
            if(it == RRH.end()){  
                RRH.insert(make_pair(v,0)); 
                State[v] = 2;

                if(G.t[i]>=ku) ++RRH[v];
            }
            else{if(G.t[i]>=ku) ++ it->second;}
        }
    }
    
   
    return;   
}

// Update R due to add
// R and state
inline void Raddupdate(unsigned val,unsigned short * State){
  
    for(unsigned i = G.h[val]; i<G.h[val+1];++i){
        unsigned v = G.e[i];
        if(State[v]==0 && FindT(v) > kl ) {
            map<unsigned,unsigned>:: iterator it = RR.find(v);     // new update, note that the score of RR is counducted in ExaScore
            if(it == RR.end()){RR.insert(make_pair(v,0));State[v] = 2;}
            //else{++ it->second;}
        }
    }

    return; 
}
//Update R due to delete
// R and state
inline void Rdelupdate(unsigned val,unsigned short * State){
    
    // for(unsigned i = G.h[val]; i<G.h[val+1];++i){
    //     unsigned v = G.e[i];
    //     if(State[v]==2 && G.t[i]> kl_cur ) {
    //        // unordered_map<unsigned,unsigned>:: iterator it = RR.find(v);
    //        // if(it!=RR.end() && it->second > 0) --it->second;
    //     }
    // }  
    return; 

}

// computer the cohesive score of H
void ScoreSta(unsigned short *State, vector<khash_t(32)*> Sub, khash_t(32)* G2Sid, 
    vector<unsigned>& ID2sup){

    double   num_edge     = 0;
    double   num_trianle  = 0;
    double   num_triple   = 0;
    double   max_degree   = 999999;
    double   average      = 0;
    double   density      = 0;
    double   cluster      = 0;


    // compute vertex support of every vertex   new update

    queue<unsigned> Q;
    Q.push(query);

    bool *inqueue = new bool[G.hs]{};
    inqueue[query] = true;

    while(!Q.empty()){
        unsigned u = Q.front();
        Q.pop();

        unsigned num_deg = 0;

        for(unsigned i = G.h[u];i<G.h[u+1];++i){ 
    
            unsigned w = G.e[i];
            if(H.find(w)==H.end()) continue;
    

            // unsigned ID = EdgeIDinSub(u,w,Sub,G2Sid);

            if(VerSup[kh_value(G2Sid,kh_get(32,G2Sid,w))]>= kl-2){

                ++num_deg;
                unsigned sum_tri = 0;
             
                for(unsigned j = G.h[u];j<G.h[u+1];++j){   
                  
                    unsigned v = G.e[j];
                    if(H.find(v)==H.end()) continue;

                    if(kh_get(32,khash[w],v)!= kh_end(khash[w])){    
                        //unsigned v = it->first;

                        if(VerSup[kh_value(G2Sid,kh_get(32,G2Sid,v))] >= kl-2){           
                            ++ sum_tri;
                        }
                              
                    }
                }
                num_trianle += sum_tri;

                if(!inqueue[w]) { Q.push(w);inqueue[w]=true;}
                
            }

            
        }

        max_degree = max_degree>num_deg?num_deg:max_degree;
        num_edge += num_deg;
        num_triple += num_deg*(num_deg-1)/2;;
                
    }

    num_edge /= 2;
    num_trianle /= 6;

    average = 2*num_edge/H.size();
    density = 2*num_edge/(H.size()*(H.size()-1));

    if(num_triple!=0 )cluster = 3*num_trianle/num_triple;
    else              cluster = 0;

  
    if(kl > r_t ||(kl == r_t && r_a <= average && r_d <= density && r_c <= cluster)){
        r_a = average;
        r_d = density;
        r_c = cluster;
        r_t = kl;
        r_s = H.size();
    }
    

    return ;

}

// for return directly
void ScoFin(){

    if(H.size() == 1) {
        // ofstream out(path_out,ios::app);
        // // out<<0<<" "<<0<<" ";
        // // out<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0;
        // out<<0<<" "<<0<<" "<<0;
        //out.close();
        return;
    }

    vector<unsigned> ID2sup; vector<vector<unsigned>> ID2end;
    vector<khash_t(32)*> Sub; khash_t(32)* G2Sid = kh_init(32);;
    unsigned short * InC = new unsigned short[G.hs]{};
    for(unsigned val: H){
        InC[val] = 1;
    }


    TriUpdate(query,Sub,G2Sid,ID2sup,ID2end);
    for(unsigned val : H){
        if(val == query) continue;
        TriUpdate(val,Sub,G2Sid,ID2sup,ID2end);
    }


    TrussDecomp(Sub,G2Sid,ID2sup,ID2end);

    ScoreSta(InC,Sub,G2Sid,ID2sup);

    delete [] InC;
    ID2sup.resize(0);
    ID2end.resize(0);
    Sub.resize(0);
    kh_destroy(32,G2Sid);

    return;
}

inline void CliStart(){
    if(ku > 2 && (ku == trq || trq >= size_u)){
        vector<unsigned> Cli;
        FindCli(Cli);
 
        if(Cli.size()==ku && Cli.size() >= size_l){
            H.clear();
            for(unsigned val:Cli){H.insert(val);}
            kl = ku;
            r_a = H.size()-1;
            r_d = 1;
            r_c = 1;
            //ScoFin();
            return;
        }

        C = Cli;
    }else if (H_cur.size() > 2){
        C = H_cur; 
    }
    else{
        C.emplace_back(query);
    }


}

/////////////////all maxmial clique including query///////////////////
inline void BK(unordered_set<unsigned> &M,unordered_set<unsigned> &P,unordered_set<unsigned> &X,vector<vector<unsigned>>& CliResult,bool * hop){
    if(omp_get_wtime() - pro_start_time >= maxtime) {OverTime = true; return;}
    if(P.size()==0 && X.size()==0){
        vector<unsigned> temp;
        for(unsigned val:M) temp.emplace_back(val);
        CliResult.emplace_back(temp);
        return;
    }

   if(P.size()==0) return;
    //find pivot vertex has the most neighbor in P add X
    unsigned pivot_ver = *P.begin();
    unsigned pivot_deg = 0;
   
    for(unsigned val:P){
        unsigned cur_deg = 0;
        for(unsigned i=G.h[val];i<G.h[val+1];++i){
            if(hop[G.e[i]]) ++cur_deg;
        }
        
        if(cur_deg > pivot_deg){
            pivot_ver = val;
            pivot_deg = cur_deg;
        }
    }

    for(unsigned val:X){
        unsigned cur_deg = 0;
        for(unsigned i=G.h[val];i<G.h[val+1];++i){
            if(hop[G.e[i]]) ++cur_deg;
        }
        
        if(cur_deg > pivot_deg){
            pivot_ver = val;
            pivot_deg = cur_deg;
        }
    }

    unordered_set<unsigned> cycle;
    for(auto val:P) cycle.insert(val);
    for(unsigned i=G.h[pivot_ver];i<G.h[pivot_ver+1];++i){
            if(hop[G.e[i]]) cycle.erase(G.e[i]);
    }

    for(unsigned val:cycle){
      
        unordered_set<unsigned> M1,P1,X1;
        for(auto i:M) M1.insert(i);
        M1.insert(val);

        for(unsigned i=G.h[val];i<G.h[val+1];++i){
            if(hop[G.e[i]] && P.find(G.e[i])!=P.end()) P1.insert(G.e[i]);
        }   

        for(unsigned i=G.h[val];i<G.h[val+1];++i){
            if(hop[G.e[i]] && X.find(G.e[i])!=X.end()) X1.insert(G.e[i]);
        }

        BK(M1,P1,X1,CliResult,hop);
        M1.clear();P1.clear();X1.clear();
        P.erase(val);
        X.insert(val);
    }


    //Nu.clear();

}

inline vector<vector<unsigned>> BKStart(unordered_set<unsigned>& CliInit){

    vector<vector<unsigned>> CliResult;
    unordered_set<unsigned> M,P,X;
    bool * hop = new bool[G.hs]{};
    for(unsigned val : CliInit){
        P.insert(val);
        hop[val] = true;
    }

    BK(M,P,X,CliResult,hop);

    M.clear();P.clear();X.clear();
   
    delete [] hop;
   
    return CliResult;
}


void Heuristic(){

    //CliStart();
 
    // init  for store sub
    vector<unsigned> ID2sup; vector<vector<unsigned>> ID2end;
    vector<khash_t(32)*> Sub; khash_t(32)* G2Sid = kh_init(32);;
    unsigned short * State = new unsigned short [G.hs]{};

    //update C and new R
    for(unsigned val : C){
        State[val] = 1;
    }
    for(unsigned val : C){
        RUpdate(val,State);
    }

    while(C.size()< size_u){
        unsigned ver_cand = ScoreHeu();
        //cout<<"scorebest "<<ver_cand<<endl;
        C.emplace_back(ver_cand);     
        RUpdate(ver_cand,State);          
    } 

    TriCalculate(C,State,Sub,G2Sid,ID2sup,ID2end);

    TrussDecomp(Sub,G2Sid,ID2sup,ID2end);

    double temp = omp_get_wtime();
    ScoreSta(State,Sub,G2Sid,ID2sup);
    pro_start_time += (omp_get_wtime() - temp);


    delete [] State;
    ID2sup.resize(0);
    ID2end.resize(0);
    Sub.resize(0);
    kh_destroy(32,G2Sid);

    return;

}

inline void CliHeu(){

    if(ku == 3) return;
     //find clique 
    unordered_set<unsigned> CliInit;          // used for Find clique
    vector<vector<unsigned>>BKresult;         //return result of clique
    unordered_set<unsigned> same;

    unsigned best_ver;
    bool flag_ver = true;
    for(NEIDEG val: key){
        if(same.find(val.val)==same.end()) {best_ver = val.val;flag_ver = false;break;}
    }
    if(flag_ver) return;

    for(unsigned i=G.h[best_ver];i<G.h[best_ver+1];++i){
        unsigned v = G.e[i];
        if(NeiDeg.find(v)!=NeiDeg.end()) CliInit.insert(v);
    }
    

    while(CliInit.size()>0){
        BKresult.clear();
        unordered_set<unsigned>:: iterator it = CliInit.begin();
        if(CliInit.size()==1){
            vector<unsigned> cur;cur.emplace_back(*it);
            BKresult.emplace_back(cur);
        }else if(CliInit.size()==2){
            unsigned u = *it;
            unsigned v = *(++it);
            if(kh_get(32,khash[u],v)!=kh_end(khash[u])){
                vector<unsigned> cur;cur.emplace_back(u);cur.emplace_back(v);
                BKresult.emplace_back(cur);
            }
            else{
                vector<unsigned> cur;cur.emplace_back(u);
                vector<unsigned> cur1;cur1.emplace_back(v);
                BKresult.emplace_back(cur);BKresult.emplace_back(cur1);
            }
        }else{
            BKresult = BKStart(CliInit);
        }
        
        if(OverTime) return;

        same.insert(best_ver);

        if(BKresult.size()>0){
            unsigned id=0,best_s=0,best_d=0;

            for(unsigned i=0;i<BKresult.size();++i){
                unsigned cur_size = BKresult[i].size();
                unsigned cur_deg = 0;

                for(unsigned val:BKresult[i]){
                    same.insert(val);
                    for(unsigned j=G.h[val];j<G.h[val+1];++j){
                        if(G.t[j]>=ku) ++ cur_deg;
                    }

                } 

                if(cur_size > best_s ){
                    best_s = cur_size;
                    best_d = cur_deg;
                    id = i;
                }
                else if(cur_size == best_s && cur_deg > best_d){
                    best_s = cur_size;
                    best_d = cur_deg;
                    id = i;
                }
            }

            C.clear();
            // H.clear();
            RRH.clear();
            for(unsigned val:BKresult[id]) {C.emplace_back(val);}
           
            C.emplace_back(query);C.emplace_back(best_ver);

            Heuristic();

            if(kl==ku) return;

            //diversity 
            if(BKresult.size()>1){
                set<unsigned> cli_cur;
                for(unsigned val:BKresult[id]) {cli_cur.insert(val);}
                id = 0;best_d=0;
                for(unsigned i=0;i<BKresult.size();++i){
                    unsigned cur_d = 0;
                    for(unsigned val:BKresult[i]){
                        if(cli_cur.find(val)==cli_cur.end()) ++cur_d;
                    }

                    if(cur_d > best_d){
                        best_d = cur_d;
                        id = i;
                    }

                }

                C.clear();
                // H.clear();
                RRH.clear();
                for(unsigned val:BKresult[id]) {C.emplace_back(val);}
               
                C.emplace_back(query);C.emplace_back(best_ver);

                Heuristic();
      
                if(kl==ku) return;

            }


        }

        flag_ver = true;
        for(NEIDEG val: key){
            if(same.find(val.val)==same.end()) {best_ver = val.val;flag_ver = false;break;}
        }
        if(flag_ver) return;  
        
        CliInit.clear();
        for(unsigned i=G.h[best_ver];i<G.h[best_ver+1];++i){
            unsigned v = G.e[i];
            if(NeiDeg.find(v)!=NeiDeg.end()) CliInit.insert(v);
        }

    }
}



// Update R Cstar if kl Updated
inline void UpdateRCstar(unsigned short * State){

    Cstar.clear();
    for(map<unsigned,unsigned>::iterator iter = CC.begin();iter!=CC.end();++iter){
        if(FindT(iter->first) <= kl) Cstar.insert(iter->first);
    }
    

    // queue<unsigned> Q;
    // Q.push(query);
    // bool * inqueue = new bool [G.hs]{};
    // inqueue[query] = true;

    // while(!Q.empty()){
    //     unsigned u = Q.front();
    //     Q.pop();

    //     for(unsigned i=G.h[u];i<G.h[u+1];++i){  

    //         unsigned w = G.e[i];

    //         if(CC.find(w)==CC.end()) continue;
    //         if(kh_value(khash[w],kh_get(32,khash[w],u)) >kl && !inqueue[w]) {Q.push(w);inqueue[w]=true;}
    //     }
    // }

    // for(pair<unsigned,unsigned> val:CC){
    //     if(!inqueue[val.first]) Cstar.insert(val.first);
    // }

    // if(Cstar.size() == 0){
    //     delete [] inqueue;
    //     return;
    // }

    // for(map<unsigned,unsigned>::iterator iter = CC.begin();iter!=CC.end();++iter){
    //     unsigned sum_deg = 0;
    //     for(unsigned i=G.h[iter->first];i!=G.h[iter->first+1];++i){
    //         map<unsigned,unsigned>::iterator it = CC.find(G.e[i]);
    //         if(it!=CC.end()&&G.t[i]>kl){
    //             ++ sum_deg;
    //         }
    //     }
    //     iter->second = sum_deg;
    // }

    // delete [] inqueue;
    return;
}

// judge whether is there vertex  belongs to Cstar in C
// parameter flag = true shows Cstar can delete ver
// return true shows it is in C
bool JudCstar(bool flag,const vector<unsigned>& Ver = vector<unsigned>()){
    if(flag){
        for(unsigned val:Ver) if(Cstar.find(val)!=Cstar.end()) {Cstar.erase(val);}
    }
    
    for(pair<unsigned,unsigned> val:CC) if(Cstar.find(val.first)!=Cstar.end()) return true;
    return false;
}

// score fun in exact
// return 0 shows there is no vertex satisfies strategies
inline unsigned ScoreExact(vector<unsigned>&VerBest, unsigned short * State){
    unsigned ver_truss,ver_cur;
    unsigned best_ver = 0,best_truss = 0;
    unsigned flag_re = false;
    unsigned budget = size_u - CC.size() - 1;
    int best_score = -1,ver_score;   // new update
    vector<unsigned> temp;
   
    for(map<unsigned,unsigned>:: iterator it = RR.begin();it!=RR.end();++it){
        ver_cur   = it->first;
        ver_score = 0;
        ver_truss = FindT(ver_cur);
        unsigned sum_score = 0;  //for < ku
        unsigned del_flag = false;
        for(unsigned i = G.h[ver_cur]; i<G.h[ver_cur+1];++i){
            unsigned v = G.e[i];
            if(State[v]==1 && G.t[i]>kl) {
                ++ ver_score;
            } else if(State[v]==1 && ver_truss > kl){
                del_flag = true;
            }
        }
 
        // it->second = sum_score;
    
        if(ver_score==0 && del_flag) temp.emplace_back(ver_cur);  
        
        //reduce1   
        if(budget + ver_score < kl) continue;           // reudce re1

        if( ver_score > best_score){
            flag_re = true;
            best_truss = ver_truss;
            best_score = ver_score;
            best_ver   = ver_cur;
        }else if(ver_score == best_score && ver_truss > best_truss){
            best_truss = ver_truss;
            best_score = ver_score;
            best_ver   = ver_cur;
        }

    }

    for(auto val:temp) {State[val] = 0;RR.erase(val);}
    // store best_ver
    if(!flag_re) {return 0;}

    VerBest.emplace_back(best_ver);
    //reduce re2
    if(best_truss == kl+1){
          
        unsigned sum  = 0;
        bool flag_equ = true;
        for(unsigned i = G.h[best_ver]; i<G.h[best_ver+1];++i){   // all the neighbor are lawful
            if(G.t[i]>kl){
                if(State[G.e[i]]==3){flag_equ = false;break;}
                ++sum;
            }
        }

        if(sum == kl && flag_equ){
            for(unsigned i = G.h[best_ver]; i<G.h[best_ver+1];++i){
                if( G.t[i]>kl && State[G.e[i]]!= 1) { VerBest.emplace_back(G.e[i]);}
            }
        }
    }

    return 1;

}

inline void Domination(unsigned val){
 
    map<unsigned,DOM>:: iterator it = Dom.find(val);
    if(it == Dom.end()){     // need computing 
        DOM temp;
        Dom.insert(make_pair(val,temp));
        for(unsigned i = G.h[val];i<G.h[val+1];++i){
            
            unsigned v = G.e[i];
            if(G.h[val+1]-G.h[val] > G.h[v+1]-G.h[v]){  //withdraw
                bool flag_cur = true;
                for(unsigned j = G.h[v];j<G.h[v+1];++j){
            
                    if(kh_get(32,khash[val],G.e[j]) == kh_end(khash[val]) && G.e[j]!=val){flag_cur = false;break;}
                }

                if(flag_cur)  temp.withdraw.emplace_back(v);
            }
            else if(G.h[val+1]-G.h[val] == G.h[v+1]-G.h[v]){                                      //toge
                bool flag_cur = true;
                for(unsigned j = G.h[val];j<G.h[val+1];++j){
                    if(kh_get(32,khash[v],G.e[j]) == kh_end(khash[v])&& G.e[j]!=v){flag_cur = false;break;}
                }

                if(flag_cur) temp.toge.emplace_back(v);
            }               
        }
    
    }
    else{
        // return;
    }
    
}

inline void RaddVerBest(unsigned ver, unsigned short * State,unsigned kl_cur){
    
    if(FindT(ver)<=kl_cur) return;

    // for(unsigned i = G.h[ver]; i<G.h[ver+1];++i){
    //     unsigned v = G.e[i];
    //     if(State[v]==1 && G.t[i]>kl_cur ) {
    //         ++ cur_sum;        
    //     }
    // }
    RR.insert(make_pair(ver,0));
    State[ver] = 2;

    return;
}

// Truss update with edge insert designed by TCP paper
// the function is the subset of TrussInsert
// supcur stores current trussness-2 of every edge
void TrussInsertSub(vector<unsigned> & supcur, khash_t(32)* G2Sid, vector<khash_t(32)*>& Sub, vector<vector<unsigned>>& ID2end, unsigned u, unsigned v){
    
    // compute k1 k2
    unsigned ptr_l, ptr_s;
    if(G.h[u+1] - G.h[u] > G.h[v+1] - G.h[v]){
        ptr_l = u;
        ptr_s = v;
    }
    else{
        ptr_l = v;
        ptr_s = u;
    }

    multimap<unsigned,vector<unsigned>,greater<unsigned>> SortSup;  // sort sup of neighbor's edges
   
    for(unsigned j = G.h[ptr_s];j<G.h[ptr_s+1];++j){
        unsigned w = G.e[j];
        if(kh_get(32,G2Sid,w) == kh_end(G2Sid)) continue;


        if(kh_get(32,khash[ptr_l],w) != kh_end(khash[ptr_l])){
                
            // v w
            unsigned te1, te2;
            if(v<w){
                te1 = kh_value(G2Sid,kh_get(32,G2Sid,v));
                te2 = w;
            }
            else{
                te1 = kh_value(G2Sid,kh_get(32,G2Sid,w));
                te2 = v;
            }

            k = kh_get(32,Sub[te1],te2);
            if(k == kh_end(Sub[te1])) continue;

            unsigned id1 = kh_value(Sub[te1],k);
            unsigned sup1 = supcur[id1];
            
            if(u<w){
                te1 = kh_value(G2Sid,kh_get(32,G2Sid,u));
                te2 = w;
            }
            else{
                te1 = kh_value(G2Sid,kh_get(32,G2Sid,w));
                te2 = u;
            }

            k = kh_get(32,Sub[te1],te2);
            if(k == kh_end(Sub[te1])) continue;

            unsigned id2 = kh_value(Sub[te1],k);
            unsigned sup2 = supcur[id2];
            
            sup1 = sup1 < sup2 ? sup1: sup2;
            vector<unsigned> temp(2);
            temp[0] = id1; temp[1] = id2;
            SortSup.insert(make_pair(sup1,temp));
            temp.resize(0);
        }     
    }// iterate neighbor edge of insterting edge e

    unsigned k1 = 0,k2 = 0;
    unsigned sum_n = 0;
    bool flag = true;
    for(multimap<unsigned,vector<unsigned>,greater<unsigned>>:: iterator it = SortSup.begin();it!=SortSup.end();++it){
        ++ sum_n;
        unsigned temp = it->first;
        if(flag && sum_n >= temp) {
            k1 = temp;
            flag = false;
        }
        if(sum_n >= temp + 1){
            k2 = temp + 1;
            break;
        }
    }

    // ps k1, case is <2>
    k1 = (flag == true && SortSup.size() != 0)? SortSup.size() : k1;

    // ps k2, case is <1>
    k2 = k2 < k1 ? k1 : k2;

    //update sup of insert edge
    supcur.emplace_back(k1);

    
    /***************************************************/
    if(k2 == 0 || (--SortSup.end())->first >= k2) return; // kmax = k2 - 1 = -1 
    else --k2;          // k2 --> kmax
    
    // Initialize the vector that stroing edges with possibly changing
    vector<map<unsigned,unsigned>> LK (k2+1);

    //iterate neighbor's edges of inster edge, SortSup stores them
    multimap<unsigned,vector<unsigned>,greater<unsigned>>:: iterator iter = --SortSup.end();
    for(unsigned i = 0; i< SortSup.size();++i){
        unsigned kk = iter->first;
        if(kk<= k2){
            unsigned id1 = iter->second[0];
            unsigned id2 = iter->second[1];
            if(supcur[id1] == kk) LK[kk].insert(make_pair(id1,0));
            if(supcur[id2] == kk) LK[kk].insert(make_pair(id2,0));

        }else{
            break;
        }
        --iter;
    }

    /******************************************************/
    // search edge with connecting by ajjacent triangles k
    for(int i = k2;i>=0;--i){
        if(LK[i].size()==0) continue;

        queue<unsigned> Q;
        for(map<unsigned,unsigned>::iterator iter =  LK[i].begin(); iter!=LK[i].end();++iter){
            Q.push(iter->first);
            //cout<<iter->first<<endl;
        }

        //cout<<"Q size:"<<Q.size()<<endl;
        unordered_set<unsigned> NEXT;

        while(!Q.empty()){
            unsigned idcur = Q.back();
            Q.pop();

            unsigned u1 = ID2end[idcur][0];
            unsigned v1 = ID2end[idcur][1];

            unsigned ptr_l, ptr_s;
            if(G.h[u1+1] - G.h[u1] > G.h[v1+1] - G.h[v1]){
                ptr_l = u1;
                ptr_s = v1;
            }
            else{
                ptr_l = v1;
                ptr_s = u1;
            }

            unsigned sum_s = 0;

            for(unsigned j = G.h[ptr_s];j<G.h[ptr_s+1];++j){
                unsigned w1 = G.e[j];
    
                if(kh_get(32,G2Sid,w1) == kh_end(G2Sid)) continue;

                if(kh_get(32,khash[ptr_l],w1) != kh_end(khash[ptr_l])){
                        
                    // v w
                    unsigned te1, te2;
                    if(v1<w1){
                        te1 = kh_value(G2Sid,kh_get(32,G2Sid,v1));
                        te2 = w1;
                    }
                    else{
                        te1 = kh_value(G2Sid,kh_get(32,G2Sid,w1));
                        te2 = v1;
                    }

                    //cout<<"        point1"<<endl;
                    k = kh_get(32,Sub[te1],te2);
                    if(k == kh_end(Sub[te1])) continue;
                    unsigned id1 = kh_value(Sub[te1],k);
                    
                    if(u1<w1){
                        te1 = kh_value(G2Sid,kh_get(32,G2Sid,u1));
                        te2 = w1;
                    }
                    else{
                        te1 = kh_value(G2Sid,kh_get(32,G2Sid,w1));
                        te2 = u1;
                    }

                    //cout<<"        point2"<<endl;
                    k = kh_get(32,Sub[te1],te2);
                    if(k == kh_end(Sub[te1])) continue;
                    unsigned id2 = kh_value(Sub[te1],k);
                    
                    //cout<<"        point3"<<endl;
                    if(supcur[id1] < i || supcur[id2] < i) continue;
                    ++ sum_s;
                    //cout<<"id1:"<<id1<<" val:"<<supcur[id1]<<endl;
                    //cout<<"id2:"<<id1<<" val:"<<supcur[id2]<<endl;
                    if(supcur[id1] == i && LK[i].find(id1)== LK[i].end()) {LK[i].insert(make_pair(id1,0)); Q.push(id1);}
                    if(supcur[id2] == i && LK[i].find(id2)== LK[i].end()) {LK[i].insert(make_pair(id2,0)); Q.push(id2);}

                }     
            }// 

            LK[i][idcur] = sum_s;
            if(sum_s <= i) NEXT.insert(idcur);
                
        }// Q
        //cout<<"LK size:"<<LK[i].size()<<endl;
        QueueClear(Q);

        /****************************************/
        while (!NEXT.empty()){  
            unsigned idcur = *NEXT.begin();
            NEXT.erase(idcur);
            LK[i].erase(idcur);
            
            unsigned u1 = ID2end[idcur][0];
            unsigned v1 = ID2end[idcur][1];

            unsigned ptr_l, ptr_s;
                if(G.h[u1+1] - G.h[u1] > G.h[v1+1] - G.h[v1]){
                    ptr_l = u1;
                    ptr_s = v1;
                }
                else{
                    ptr_l = v1;
                    ptr_s = u1;
                }

                for(unsigned j = G.h[ptr_s];j<G.h[ptr_s+1];++j){
                    unsigned w1 = G.e[j];
                    if(kh_get(32,G2Sid,w1) == kh_end(G2Sid)) continue;

                    if(kh_get(32,khash[ptr_l],w1) != kh_end(khash[ptr_l])){
                            
                        // v w
                        unsigned te1, te2;
                        if(v1<w1){
                            te1 = kh_value(G2Sid,kh_get(32,G2Sid,v1));
                            te2 = w1;
                        }
                        else{
                            te1 = kh_value(G2Sid,kh_get(32,G2Sid,w1));
                            te2 = v1;
                        }

                        k = kh_get(32,Sub[te1],te2);
                        if(k == kh_end(Sub[te1])) continue;
                        unsigned id1 = kh_value(Sub[te1],k);
                        
                        if(u1<w1){
                            te1 = kh_value(G2Sid,kh_get(32,G2Sid,u1));
                            te2 = w1;
                        }
                        else{
                            te1 = kh_value(G2Sid,kh_get(32,G2Sid,w1));
                            te2 = u1;
                        }

                        k = kh_get(32,Sub[te1],te2);
                        if(k == kh_end(Sub[te1])) continue;
                        unsigned id2 = kh_value(Sub[te1],k);

                        if(supcur[id1] < i || supcur[id2] < i) continue;

                        if(supcur[id1] == i && LK[i].find(id1) == LK[i].end()) continue;
                        if(supcur[id2] == i && LK[i].find(id2) == LK[i].end()) continue;

                        if(LK[i].find(id1) == LK[i].end()){
                            --LK[i][id1];
                            if(LK[i][id1] <= i) NEXT.insert(id1);
                        }

                        if(LK[i].find(id2) == LK[i].end()){
                            --LK[i][id2];
                            if(LK[i][id2] <= i) NEXT.insert(id2);
                        }

                    }     
                }// 

        }

        
        for(map<unsigned,unsigned>::iterator iter =  LK[i].begin(); iter!=LK[i].end();++iter){
            supcur[iter->first] = i + 1;
            //cout<<"update:"<<supcur[iter->first]<<endl;
        }
        //cout<<"end is "<<i<<endl;
    }// for supmax -> 0 (or kmax -> 2)

    LK.resize(0);
}

// truss maintenance
// Sup stores edges (trussenss-2) value at each C
// Sub stores hash (new id u, old id v)--> id of edge
void TrussInsert(unsigned id, vector<vector<unsigned>>& Sup, vector<khash_t(32)*>& Sub,khash_t(32)* G2Sid, vector<vector<unsigned>>& ID2end){
    
    // update G2Sid;
    int G2Sid_ret;
    khiter_t k;

    //id = query 
    if(id == query){
        k = kh_put(32,G2Sid,id,&G2Sid_ret);
        kh_value(G2Sid,k) = 0;

        khash_t(32) * temp = kh_init(32);
        Sub.emplace_back(temp);

        return;
    }

    // update sub 
    khash_t(32) * temp = kh_init(32);
    Sub.emplace_back(temp);
    
    k = kh_put(32,G2Sid,id,&G2Sid_ret);
    kh_value(G2Sid,k) = Sub.size()-1;

    // insert all edges between id and C

    // supportness for current operation 
    // Sup[0] only has one edge
    vector<unsigned> supcur;
    if(CC.size() > 2){
        for(unsigned val:Sup[CC.size()-3]){
            supcur.emplace_back(val);
        }
    }  

    for(unsigned i=G.h[id];i<G.h[id+1];++i){
        unsigned v = G.e[i];
        if(kh_get(32,G2Sid,v) == kh_end(G2Sid)) continue;
        
        // update ID2end sub
        vector<unsigned> temp(2);
        temp[0] = id; temp[1] = v;
        ID2end.emplace_back(temp);
        temp.resize(0);

        // update Sub
        unsigned te1, te2;
        if(v<id){
            te1 = kh_value(G2Sid,kh_get(32,G2Sid,v));
            te2 = id;
        }
        else{
            te1 = kh_value(G2Sid,kh_get(32,G2Sid,id));
            te2 = v;
        }

        k = kh_put(32,Sub[te1],te2,&G2Sid_ret);
        kh_value(Sub[te1],k) = ID2end.size() - 1;

        // the key to update trussenss in C;

        //cout<<"**** sub in "<<id<<" "<<v<<endl;
        TrussInsertSub(supcur,G2Sid,Sub,ID2end,id,v);
       // if(supcur.size()>0)     for(auto val:supcur) cout<<val<<" "; cout<<endl;
        //cout<<" sub out "<<endl;
    }

    // update Sup
    Sup.emplace_back(supcur);
    supcur.resize(0);
}

// truss maintenance
// Sup stores edges (trussenss-2) value at each C
// Sub stores hash (new id u, old id v)--> id of edge
void TrussDel(unsigned id, vector<vector<unsigned>>& Sup, vector<khash_t(32)*>& Sub,khash_t(32)* G2Sid, vector<vector<unsigned>>& ID2end){
    
    // if id == query, the state never happen
    
    // update Sub G2SId
    Sub.pop_back();
    kh_del(32,G2Sid,kh_get(32,G2Sid,id));

    // update Sub ID2end
    for(unsigned i=G.h[id];i<G.h[id+1];++i){
        
        unsigned v = G.e[i];

        if(kh_get(32,G2Sid,v) == kh_end(G2Sid)) continue;

        unsigned te1, te2;
        if(v<id){
            te1 = kh_value(G2Sid,kh_get(32,G2Sid,v));
            te2 = id;
            kh_del(32,Sub[te1],kh_get(32,Sub[te1],te2));
        }

        ID2end.pop_back();
    }
    //update Sup
    Sup.pop_back();

}
inline bool LBOpt(unsigned short * State, vector<khash_t(32)*>& Sub,khash_t(32)* G2Sid, 
vector<unsigned> &ID2sup, vector<vector<unsigned>>& ID2end){
    //First step, initialize
    int budget_max = int(size_u) - int(CC.size());
    int cost_max = 0;

    vector<vector<int>> Candidate (CC.size(),vector<int>(2,0)) ;   // Index<Budget,Cost>
    vector<int> bmax(CC.size());       // store bmax(u) for updating budget
    

    // second budget and cost
    for(map<unsigned,unsigned>::iterator it_CC = CC.begin(); it_CC != CC.end();++it_CC){

        unsigned u = it_CC->first;
        unsigned ID_Candidate = kh_value(G2Sid,kh_get(32,G2Sid,u));
        

        // budget and vertex support
        unsigned sum_neigh = 0;
        unsigned SupU = 0;
        for(unsigned i= G.h[u];i<G.h[u+1];++i){
            unsigned v = G.e[i];
          
            if(State[v] == 2){
                ++ sum_neigh;
            }
            else if(State[v] == 1){
                unsigned SupC = ID2sup[EdgeIDinSub(u,v,Sub,G2Sid)];
                SupU = SupU >  SupC ? SupU : SupC;
            }
        }

        Candidate[ID_Candidate][0] = sum_neigh < budget_max ? sum_neigh:budget_max;
        Candidate[ID_Candidate][1] = (kl-1) <= SupU ? 0 : kl-1-SupU;
        cost_max = cost_max > Candidate[ID_Candidate][1] ? cost_max : Candidate[ID_Candidate][1];
        bmax[ID_Candidate] = Candidate[ID_Candidate][0];
        // budget < cost, return false;
        if( Candidate[ID_Candidate][0] < Candidate[ID_Candidate][1]) return false;

    }

    //fast judget operation
    if(*min_element(bmax.begin(),bmax.end()) >= 2*cost_max) return true;

    // third step, updating budget
    for(map<unsigned,unsigned>::iterator it = CC.begin(); it != CC.end();++it){

        unsigned u = it->first; 
        unsigned ID_Candidate = kh_value(G2Sid,kh_get(32,G2Sid,u));   // use new id of G2Sid

        // cost(u) = 0
        if(Candidate[ID_Candidate][1] == 0) continue;

        // find all neighbors of u in R
        set<unsigned> NeighOut;
        for(unsigned i= G.h[u];i<G.h[u+1];++i){
            unsigned v = G.e[i]; 
            if(State[v] == 2){
                NeighOut.insert(v);
            }
        }

        // conduct the operation of updating budget
        map<unsigned,unsigned>::iterator it_cur =  it;  // it_cur = it shows current point
        ++ it_cur;
        for(; it_cur != CC.end();++it_cur){

            unsigned ver = it_cur->first;
            unsigned id_ver = kh_value(G2Sid,kh_get(32,G2Sid,ver));

            if(bmax[id_ver] >= Candidate[id_ver][1] + Candidate[ID_Candidate][1]) continue;

            bool flag = true;
            if(G.h[ver+1]-G.h[ver] > NeighOut.size()){
                for(set<unsigned>::iterator it_set = NeighOut.begin();it_set != NeighOut.end();++ it_set){
                    if(kh_get(32,khash[ver],*it_set) != kh_end(khash[ver])){
                        flag = false;
                        break;
                    }
                }
            }
            else{
                for(unsigned j = G.h[ver];j<G.h[ver+1];++j){
                    if(NeighOut.find(G.e[j])!=NeighOut.end()){
                        flag = false;
                        break;
                    }
                }
            }
            
            if(flag){
                // updating budget
                if(budget_max -Candidate[ID_Candidate][1] < Candidate[id_ver][1]) return false;   // new updata bmax[id_ver]
            }
        } 

    }// second step end 

    return true;   
}



// Update C and R for exact
// ture shows add VerBest to C;
// false shows delete verbest from C;
inline void UpdateCR(vector<vector<unsigned>>& Sup,unsigned short * State, vector<khash_t(32)*>& Sub,khash_t(32)* G2Sid, 
vector<vector<unsigned>>& ID2end,bool flag,unsigned ver){

    if(flag){
        //C
        //cout<<"cr "<<ver<<endl;
        // unsigned sum_deg = 0;
        // for(unsigned i=G.h[ver];i!=G.h[ver+1];++i){
        //     map<unsigned,unsigned>::iterator it = CC.find(G.e[i]);
        //     if(it!=CC.end()&&G.t[i]>kl){
        //         ++ it->second;
        //         ++ sum_deg;
        //     }
        // }

        CC.insert(make_pair(ver,0));
        //TriUpdate(ver,Sub,G2Sid,ID2sup,ID2end);
        //cout<<"CC: ";
        //for(auto val:CC) cout<<val.first<<" "; cout<<endl;
        //cout<<"in insert"<<endl;
        TrussInsert(ver,Sup,Sub,G2Sid,ID2end);
        //cout<<"out insert"<<endl;  
        State[ver] = 1;
        RR.erase(ver);

        //R
        Raddupdate(ver,State);
    }
    else{

        // for(unsigned i=G.h[ver];i!=G.h[ver+1];++i){
        //     unsigned v = G.e[i];
        //     map<unsigned,unsigned>::iterator it = CC.find(v);
        //     if(it!=CC.end()&&G.t[i]>kl){
        //         -- it->second;
        //     }
        // }
        CC.erase(ver);
       // cout<<"in del"<<endl;   
        TrussDel(ver,Sup,Sub,G2Sid,ID2end);
       // cout<<"out del"<<endl;  
        //TriDele(ver,Sub,G2Sid,ID2sup,ID2end);
        State[ver] = 3;
    
        Rdelupdate(ver,State);
    }

    return;
}
struct comp//compare
{
    bool operator()(const NEIDEG &a, const NEIDEG &b)  //greater sort
    {
        return a.deg >= b.deg;
    }
};

inline void LB3Updata(set<NEIDEG,comp> & Candidata,set<unsigned>& VerDel){
    //cout<<"size "<<VerDel.size()<<endl;
    set<NEIDEG,comp>::iterator it = Candidata.begin();
    vector<set<NEIDEG,comp>::iterator> loc;
    for(;it!=Candidata.end();++it){
        if(VerDel.find((*it).val)!=VerDel.end()) {  
               
            loc.emplace_back(it);  
                  
        }
    }

    for(set<NEIDEG,comp>::iterator iter: loc){
        if((*iter).deg > 1){    
            NEIDEG temp;
            temp.deg = (*iter).deg - 1;
            temp.val = (*iter).val; 
            Candidata.erase(iter);
            Candidata.insert(temp);
        }
        else{
            Candidata.erase(iter);
        }

    }
    loc.resize(0);

}

inline bool LowerBound3(set<NEIDEG,comp>& Candidata, unsigned sum_gkl){

    
    map<unsigned,unsigned> CC_temp(CC);
    
    set<unsigned> RR_used;

    unsigned budget = size_u - CC.size();

    set<NEIDEG,comp>::iterator it = Candidata.begin();

    unsigned sum_new  = 0;   // new vertex ?=kl due to RR

    for(unsigned i=0;i<budget;++i){
    
        if( sum_gkl+sum_new == CC.size() || it == Candidata.end() ) break;
            unsigned v = (*it).val;
            
            if((*it).deg == 0) break;
            RR_used.insert(v);
            for(unsigned j=G.h[v];j<G.h[v+1];++j){

                map<unsigned,unsigned>::iterator iter = CC_temp.find(G.e[j]);
            
                if( iter ==CC_temp.end() || G.t[j]<=kl) continue;
    
                ++ iter->second;

                if(iter->second == kl){
                    ++ sum_new;
                    unsigned w = iter->first;
                    set<unsigned> VerDel;

                    for(unsigned k=G.h[w];k<G.h[w+1];++k){
                        unsigned u = G.e[k];
                        if(G.t[k]>kl && RR.find(u)!=RR.end() && RR_used.find(u) == RR_used.end()){  // in Candidaata not in RR_used
                            VerDel.insert(u);
                        }
                    }
            
                    LB3Updata(Candidata,VerDel);
                    VerDel.clear();
                }
            }
        ++it;
    }


    if(sum_new + sum_gkl >= CC.size()) {return false;}

    return true;
}

inline bool LowerBound2(){
    
    set<NEIDEG,comp> Candidata;

    for(pair<unsigned,unsigned> val:RR){
        NEIDEG temp;
    
        if(val.second == 0) continue;
        temp.val = val.first;temp.deg = val.second;
        Candidata.insert(temp);
    }

    unsigned budget = size_u - CC.size();

    set<NEIDEG,comp>::iterator it = Candidata.begin();

    unsigned sum_new  = 0;   // new vertex ?=kl due to RR
    unsigned sum_gkl  = 0;
    
    for(pair<unsigned,unsigned> val:CC){
        if(val.second>=kl) ++sum_gkl;
    }

    float sum_lack = 0;
    for(pair<unsigned,unsigned> val:CC){
        if(val.second < kl) sum_lack += (kl - val.second) ;
    }

    map<unsigned,unsigned> CC_temp(CC);
    
    for(unsigned i=0;i<budget;++i){
        if( sum_lack  <= 0 || it == Candidata.end()) break;
        unsigned v = (*it).val;

        sum_lack -= (*it).deg ;

        ++it;
    }


    if(sum_lack <=0) return false;
    else true;


    Candidata.clear();

    return false;

}
// true,terminate
inline bool LowerBound1 (){  // general 
 
    unsigned sum_lack = 0;
    for(pair<unsigned,unsigned> val:CC){
        if(val.second < kl) sum_lack += (kl - val.second) ;
    }

    if(size_u >= (sum_lack + CC.size())) {return false;}


    if(LowerBound2()) return true;

    return false;
}

void TrussExaDecomp(vector<khash_t(32)*>& Sub,khash_t(32)* G2Sid,vector<unsigned>& ID2sup){
    ResultSubExa(Sub,G2Sid,ID2sup);
}
void Exact(vector<vector<unsigned>>& Sup,unsigned short * State, vector<khash_t(32)*>& Sub,khash_t(32)* G2Sid, vector<vector<unsigned>>& ID2end){

    // over time budget
    exact_start_time = omp_get_wtime();
    if(exact_start_time - pro_start_time >= maxtime){
        OverTime = true;
        return;
    }

    // find the best resutl
    if(kl == ku) return;

    //truss decomp
    if( (CC.size() == size_u) || (RR.size() == 0 && CC.size()>= size_l && CC.size()!=1)){
        unsigned kl_temp = kl;

        // TrussDecomp(Sub,G2Sid,ID2sup,ID2end);
        TrussExaDecomp(Sub,G2Sid,Sup[Sup.size()-1]);

        if(kl == ku) return;
        if(kl > kl_temp){
            UpdateRCstar(State);             // Update R and C;
            
            double temp = omp_get_wtime();
            ScoreSta(State,Sub,G2Sid,Sup[Sup.size()-1]);
            pro_start_time += (omp_get_wtime() - temp);

        }
        else{
            kl = kl_temp;
        }
    
    }
    

    // end condition
    if(CC.size()==size_u || RR.size()==0) return;

    if(JudCstar(false)) return;
    // select the highest vertex 


    unsigned kl_cur = kl;
    vector<unsigned> VerBest;
    if(ScoreExact(VerBest,State)==0) { return;} 

  
    if(CC.size()*(kl+1) > (size_u - CC.size())&& Sup.size() >0 && !LBOpt(State,Sub,G2Sid,Sup[Sup.size()-1],ID2end)) return;
    
    // C add verbest
    for(unsigned val:VerBest) UpdateCR(Sup,State,Sub,G2Sid,ID2end,true,val); 

    exact_start_time = omp_get_wtime();
    if(exact_start_time - pro_start_time >= maxtime){
        OverTime = true;
        return;
    }

    Exact(Sup,State,Sub,G2Sid,ID2end);
    if(kl==ku || OverTime) return;

    // C delete verbest
    // at the same time, recover C and R
    if(JudCstar(true,VerBest)) {

        for(unsigned val:VerBest) UpdateCR(Sup,State,Sub,G2Sid,ID2end,false,val); 


        for(unsigned val:VerBest)   RaddVerBest(val,State,kl_cur);

        
        return;
    }

    for(unsigned val:VerBest) UpdateCR(Sup,State,Sub,G2Sid,ID2end,false,val); 



    Exact(Sup,State,Sub,G2Sid,ID2end);


    if(kl==ku || OverTime) return;
    for(unsigned val:VerBest)   RaddVerBest(val,State,kl_cur);


    return;  
}


void FrameWork(){

    pro_start_time = omp_get_wtime();

    kl = 1;
    ku = 2;
    H.clear();
    H_cur.clear();
    C.clear();
    R.clear();
    IsLess = false;
    IsSuc = false;
    OverTime = false;
    r_t = 0;
    r_a = 0;
    r_d = 0;
    r_c = 0;
    r_s = 0;
    trq = FindT(query);

    // obtain ku
    int flag = Findku();
    if(flag == -1){IsLess = true;H.clear();H.insert(query);ScoFin();return;}
    if(flag == 0) {kl = ku;IsSuc = true; C.clear(); for(auto val:H) C.emplace_back(val); ScoFin();return;}        //[l,h]

    // heurstic algorithm
    CliStart();

    if(kl == ku) return;

    Heuristic();

    if(kl == ku) return;

    if(H_cur.size() == 0) {CliHeu();}

    //exact algorithm
    if(kl != ku){

        // H.clear();
        C.clear();
        CC.clear();
        RR.clear();

        //vector<unsigned> ID2sup; 
        vector<vector<unsigned>> ID2end;
        vector<khash_t(32)*> Sub; khash_t(32)* G2Sid = kh_init(32);
        vector<vector<unsigned>> Sup;
        // 1, in C; 2, in R; 3, forbidden
        unsigned short *State = new unsigned short[G.hs]();
        UpdateCR(Sup,State,Sub,G2Sid,ID2end,true,query);
        Exact(Sup,State,Sub,G2Sid,ID2end);
       
    }
    
    return;
}


int main(int argc,char ** argv){

    Infile(argv[1]);

    // path_out = argv[2];

    // ofstream out;
    // out.open(path_out,ios::app);

    query = atoi(argv[2]);   // query vertex
    size_l =atoi(argv[3]) ;  // the lower bound of size
    size_u = atoi(argv[4]);  // the upper bound of size
    maxtime = atoi(argv[5]); // the time limit

    if(size_u < size_l || query > G.hs) { cout<<"Input is error!"<<endl; return 0;}

    FrameWork();   // the main func

    double end_time = omp_get_wtime();
    cout<<endl<<"INPUT: query "<<query<<"; size constraint ("<<size_l<<","<<size_u<<"); time limit "<<maxtime<<endl;

    if(!IsLess)      cout<<"OUTPUT: the upper bound of min-trussness "<<ku<<"; current min-trussness "<<kl<<"; "<<endl<<"average density:"<<r_a<<"; internal density:"<<r_d<<"; clustering:"<<r_c<<"; ";
    else             cout<<"error! Note that size constraint is not suitable"<<"; ";

    cout<<"running time: "<<end_time - pro_start_time<<endl;

    cout<<endl<<"Returned subgraph:"<<endl;
    for(auto val:H) cout<<val<<" "; cout<<endl;

//    out.close();

    return 0;
}

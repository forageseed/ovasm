#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>

// Streaming exact DNA marker extraction. No assembly/recruitment classification.
struct Node { std::array<int,4> next{{-1,-1,-1,-1}}; int fail=0; unsigned mask=0; };
int main(int argc,char**argv){
    if(argc!=4)return 2;
    std::ios::sync_with_stdio(false);std::cin.tie(nullptr);
    std::array<unsigned char,256> code;code.fill(4);
    for(int i=0;i<4;i++){code[(unsigned char)"ACGT"[i]]=i;code[(unsigned char)"acgt"[i]]=i;}
    std::vector<Node> nodes(1);std::ifstream seeds(argv[1]);std::string group,word;
    if(!seeds)throw std::runtime_error("cannot open seeds");
    while(seeds>>group>>word){
        unsigned mask=group=="left"?1:group=="right"?2:group=="control"?4:0;
        if(!mask||word.empty())throw std::runtime_error("invalid seed");
        for(int strand=0;strand<2;strand++){
            std::string seq=word;
            if(strand){seq.clear();for(auto i=word.rbegin();i!=word.rend();++i){int c=code[(unsigned char)*i];if(c>3)throw std::runtime_error("bad seed base");seq+="TGCA"[c];}}
            int state=0;
            for(char b:seq){int c=code[(unsigned char)b];if(c>3)throw std::runtime_error("bad seed base");int next=nodes[state].next[c];if(next<0){next=nodes.size();nodes[state].next[c]=next;nodes.emplace_back();}state=next;}
            nodes[state].mask|=mask;
        }
    }
    std::queue<int> queue;
    for(int c=0;c<4;c++){int next=nodes[0].next[c];if(next<0)nodes[0].next[c]=0;else queue.push(next);}
    while(!queue.empty()){
        int state=queue.front();queue.pop();
        for(int c=0;c<4;c++){
            int next=nodes[state].next[c];
            if(next<0)nodes[state].next[c]=nodes[nodes[state].fail].next[c];
            else{nodes[next].fail=nodes[nodes[state].fail].next[c];nodes[next].mask|=nodes[nodes[next].fail].mask;queue.push(next);}
        }
    }
    std::ofstream hits(argv[2]),stats(argv[3]);if(!hits||!stats)throw std::runtime_error("output unavailable");
    std::string h,s,p,q;uint64_t records=0,bases=0,selected=0;std::array<uint64_t,8> counts{};
    while(std::getline(std::cin,h)){
        if(!std::getline(std::cin,s)||!std::getline(std::cin,p)||!std::getline(std::cin,q))throw std::runtime_error("truncated FASTQ");
        if(h.empty()||h[0]!='@'||p.empty()||p[0]!='+'||s.size()!=q.size())throw std::runtime_error("invalid FASTQ");
        int state=0;unsigned flags=0;
        for(unsigned char ch:s){int c=code[ch];if(c>3){state=0;continue;}state=nodes[state].next[c];flags|=nodes[state].mask;}
        records++;bases+=s.size();counts[flags]++;
        if(flags&3){hits<<h<<'\n'<<s<<'\n'<<p<<'\n'<<q<<'\n';selected++;}
        if(records%100000==0)std::cerr<<"progress records="<<records<<" bases="<<bases<<" selected="<<selected<<'\n';
    }
    if(!std::cin.eof())throw std::runtime_error("FASTQ stream read error");
    hits.flush();if(!hits)throw std::runtime_error("output write failure");
    stats<<"records\t"<<records<<"\nbases\t"<<bases<<"\nselected\t"<<selected<<'\n';
    for(int i=0;i<8;i++)stats<<"flags_"<<i<<'\t'<<counts[i]<<'\n';
    std::cerr<<"complete records="<<records<<" bases="<<bases<<" selected="<<selected<<'\n';
}

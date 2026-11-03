// Runtime zlib is installed; development header is absent. Stable public C ABI.
extern "C" {
typedef void* gzFile;
gzFile gzopen(const char*,const char*);
int gzbuffer(gzFile,unsigned);
char* gzgets(gzFile,char*,int);
const char* gzerror(gzFile,int*);
int gzclose(gzFile);
}
constexpr int Z_OK=0, Z_STREAM_END=1;
#include <cstdio>
#include <cstdint>
#include <string>
#include <sstream>
#include <array>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cerrno>
#include <cstring>
#include <chrono>
bool line(gzFile f,std::string& s) {
 s.clear(); char b[65536];
 for(;;) { char* r=gzgets(f,b,sizeof(b)); if(!r){int e;gzerror(f,&e);if(e!=Z_OK&&e!=Z_STREAM_END)throw std::runtime_error("gzip read failure");if(s.empty())return false;throw std::runtime_error("unterminated FASTQ line");}s+=b;if(!s.empty()&&s.back()=='\n')return true; }
}
std::string name(const std::string& h){if(h.empty()||h[0]!='@')throw std::runtime_error("invalid header");auto n=h.substr(1,h.find_first_of(" \t\r\n",1)-1);if(n.size()>2&&n[n.size()-2]=='/'&&(n.back()=='1'||n.back()=='2'))n.resize(n.size()-2);return n;}
uint64_t hash(const std::string& s){uint64_t h=14695981039346656037ULL^20261103ULL;for(unsigned char c:s){h^=c;h*=1099511628211ULL;}h^=h>>30;h*=0xbf58476d1ce4e5b9ULL;h^=h>>27;h*=0x94d049bb133111ebULL;return h^(h>>31);}
size_t bases(const std::string&s){size_t n=s.size();while(n&&(s[n-1]=='\n'||s[n-1]=='\r'))--n;return n;}
int main(int argc,char**argv){try{if(argc!=6)throw std::runtime_error("usage: pan_pair_sample R1.fastq.gz R2.fastq.gz OUT1.fastq OUT2.fastq REPORT.json; fixed one-in-five pair-name hash, salt=20261103; outputs must not exist");for(int i=3;i<6;i++)if(std::filesystem::exists(argv[i]))throw std::runtime_error(std::string("output already exists: ")+argv[i]);gzFile a=gzopen(argv[1],"rb"),b=gzopen(argv[2],"rb");if(!a||!b)throw std::runtime_error("open input");gzbuffer(a,1<<20);gzbuffer(b,1<<20);FILE* x=fopen(argv[3],"wx"),*y=fopen(argv[4],"wx");if(!x||!y)throw std::runtime_error("cannot exclusively create output");setvbuf(x,nullptr,_IOFBF,1<<20);setvbuf(y,nullptr,_IOFBF,1<<20);uint64_t n=0,k=0,ib=0,ob=0;std::array<std::string,4>ra,rb;auto start=std::chrono::steady_clock::now();for(;;){bool aa=line(a,ra[0]),bb=line(b,rb[0]);if(aa!=bb)throw std::runtime_error("different mate counts");if(!aa)break;for(int i=1;i<4;i++)if(!line(a,ra[i])||!line(b,rb[i]))throw std::runtime_error("truncated pair");auto key=name(ra[0]);if(key!=name(rb[0]))throw std::runtime_error("pair name mismatch at "+std::to_string(n));if(ra[2][0]!='+'||rb[2][0]!='+'||bases(ra[1])!=bases(ra[3])||bases(rb[1])!=bases(rb[3]))throw std::runtime_error("invalid FASTQ pair");n++;ib+=bases(ra[1])+bases(rb[1]);if(hash(key)%5==0){k++;ob+=bases(ra[1])+bases(rb[1]);for(int i=0;i<4;i++){if(fwrite(ra[i].data(),1,ra[i].size(),x)!=ra[i].size()||fwrite(rb[i].data(),1,rb[i].size(),y)!=rb[i].size())throw std::runtime_error("output write failure");}}if(n%10000000==0)std::cerr<<"pairs "<<n<<" kept "<<k<<"\n";}if(fclose(x)!=0||fclose(y)!=0||gzclose(a)!=Z_OK||gzclose(b)!=Z_OK)throw std::runtime_error("close/write failure");double sec=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();FILE* report=fopen(argv[5],"wx");if(!report)throw std::runtime_error("cannot exclusively create report");std::ostringstream j;j<<"{\"pairs_scanned\":"<<n<<",\"pairs_kept\":"<<k<<",\"bases_scanned\":"<<ib<<",\"bases_kept\":"<<ob<<",\"sampling_seconds\":"<<sec<<",\"salt\":20261103,\"rule\":\"FNV1a+mix64(normalized_pair_name) mod 5 == 0\",\"pair_names_checked_exactly\":true,\"read_lines_unchanged\":true}\n";auto report_text=j.str();if(fwrite(report_text.data(),1,report_text.size(),report)!=report_text.size()||fclose(report)!=0)throw std::runtime_error("report write failure");return 0;}catch(const std::exception&e){std::cerr<<e.what()<<"\nOutputs created during a failed run may be partial; they are retained, never overwritten. A successful report plus exit 0 is required.\n";return 1;}}

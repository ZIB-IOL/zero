// Exact fill-only refinement audit, grouped by parent cells.
// For each fine offset and query phase, record exactly the neighbor parent
// offset and neighbor phase. Bit sets compress these finite relations; every
// fine-cell pair tested by the flat checker is retained, without sampling.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

struct Pair { int32_t i,j; };
struct Paint { int32_t i,j,c; };
static_assert(sizeof(Pair)==8 && sizeof(Paint)==12);
template<class T> std::vector<T> read(const char* path) {
    std::ifstream f(path,std::ios::binary|std::ios::ate);
    if (!f) throw std::runtime_error("missing input");
    auto bytes=f.tellg();
    if(bytes<0 || bytes%sizeof(T)) throw std::runtime_error("malformed input");
    std::vector<T> data(static_cast<size_t>(bytes)/sizeof(T));
    f.seekg(0);
    if(!f.read(reinterpret_cast<char*>(data.data()),bytes)) throw std::runtime_error("short input");
    return data;
}

struct Relation {
    Pair q;
    uint64_t active=0;
    std::array<uint64_t,64> neighbors{};
};
struct Block {
    Pair q;
    std::array<uint64_t,5> colors{};
};

int main(int argc,char** argv) {
  try {
    if(argc!=9) throw std::runtime_error("parent rows columns fh fw offsets paint report");
    int64_t hh=std::stoll(argv[2]),ww=std::stoll(argv[3]);
    int64_t ffh=std::stoll(argv[4]),ffw=std::stoll(argv[5]);
    if(hh<=0 || ww<=0 || ffh<=0 || ffw<=0 || ffh>64 || ffw>64 || ffh*ffw>64
       || hh*ffh>INT32_MAX/2 || ww*ffw>INT32_MAX/2)
        throw std::runtime_error("unsupported dimensions");
    int32_t h=hh,w=ww,fh=ffh,fw=ffw,H=h*fh,W=w*fw,phases=fh*fw;
    auto parent=read<uint8_t>(argv[1]);
    auto offsets=read<Pair>(argv[6]);
    auto paint=read<Paint>(argv[7]);
    if(parent.size()!=static_cast<size_t>(h)*w) throw std::runtime_error("wrong parent size");
    for(auto c:parent) if(c>5) throw std::runtime_error("invalid parent label");
    std::vector<Relation> relations;
    std::unordered_map<int64_t,int32_t> index;
    for(auto d:offsets) {
        if(d.i<0 || d.i>=H || d.j<0 || d.j>=W) throw std::runtime_error("invalid offset");
        for(int32_t s=0;s<fh;++s) for(int32_t t=0;t<fw;++t) {
            int32_t ni=s+d.i,nj=t+d.j;
            int32_t qi=ni/fh,qj=nj/fw;
            if(qi>=h) qi-=h;
            if(qj>=w) qj-=w;
            int64_t key=static_cast<int64_t>(qi)*w+qj;
            auto inserted=index.emplace(key,relations.size());
            if(inserted.second) { Relation rel; rel.q={qi,qj}; relations.push_back(rel); }
            auto &rel=relations[inserted.first->second];
            int32_t p=s*fw+t,np=(ni%fh)*fw+(nj%fw);
            rel.active |= uint64_t(1)<<p;
            rel.neighbors[p] |= uint64_t(1)<<np;
        }
    }
    // Sparse blocked lookup from parent coordinates to newly painted blocks.
    int32_t bh=(h+31)/32,bw=(w+31)/32;
    std::vector<int32_t> pages(static_cast<size_t>(bh)*bw,-1);
    int32_t page_count=0;
    for(auto x:paint) {
        if(x.i<0 || x.i>=H || x.j<0 || x.j>=W || x.c<0 || x.c>=5)
            throw std::runtime_error("invalid paint");
        int32_t i=x.i/fh,j=x.j/fw;
        auto &page=pages[static_cast<size_t>(i>>5)*bw+(j>>5)];
        if(page<0) page=page_count++;
    }
    std::vector<int32_t> cells(static_cast<size_t>(page_count)*1024,-1);
    std::vector<Block> blocks;
    for(auto x:paint) {
        int32_t i=x.i/fh,j=x.j/fw;
        if(parent[static_cast<size_t>(i)*w+j]!=5) throw std::runtime_error("not a fill-only extension");
        auto page=pages[static_cast<size_t>(i>>5)*bw+(j>>5)];
        auto &entry=cells[static_cast<size_t>(page)*1024+(i&31)*32+(j&31)];
        if(entry<0) {entry=blocks.size(); Block b; b.q={i,j};blocks.push_back(b);}
        auto &b=blocks[entry];
        uint64_t bit=uint64_t(1)<<((x.i%fh)*fw+x.j%fw);
        for(auto colored:b.colors) if(colored&bit) throw std::runtime_error("duplicate paint");
        b.colors[x.c] |= bit;
    }
    uint64_t conflicts=0;
    for(size_t u=0;u<blocks.size();++u) {
        const auto &a=blocks[u];
        for(const auto &rel:relations) {
            int32_t i=a.q.i+rel.q.i,j=a.q.j+rel.q.j;
            if(i>=h) i-=h;
            if(j>=w) j-=w;
            auto c=parent[static_cast<size_t>(i)*w+j];
            if(c<5) {
                uint64_t hits=a.colors[c]&rel.active;
                while(hits) {
                    auto p=__builtin_ctzll(hits); hits &= hits-1;
                    conflicts += __builtin_popcountll(rel.neighbors[p]);
                }
            } else {
                auto page=pages[static_cast<size_t>(i>>5)*bw+(j>>5)];
                if(page<0) continue;
                auto v=cells[static_cast<size_t>(page)*1024+(i&31)*32+(j&31)];
                if(v<0) continue;
                const auto &b=blocks[v];
                for(int c=0;c<5;++c) if(b.colors[c]) {
                    uint64_t hits=a.colors[c]&rel.active;
                    while(hits) {
                        auto p=__builtin_ctzll(hits); hits &= hits-1;
                        conflicts += __builtin_popcountll(rel.neighbors[p]&b.colors[c]);
                    }
                }
            }
        }
        if(u%10000==0) std::cout<<"grouped fill audit "<<u<<"/"<<blocks.size()
                              <<" conflict witnesses "<<conflicts<<std::endl;
    }
    std::ofstream report(argv[8]);
    report<<"{\"checked\":"<<paint.size()<<",\"conflict_witnesses\":"<<conflicts
          <<",\"parent_blocks\":"<<blocks.size()<<",\"parent_offsets\":"<<relations.size()<<"}\n";
    if(!report) throw std::runtime_error("cannot write report");
    return 0;
  } catch(const std::exception &e) {std::cerr<<e.what()<<std::endl;return 1;}
}

// Exact integer label checks for all changed real cells and all mask offsets.
// Geometry is separately generated with exact Python integer arithmetic.
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

struct Pair { int32_t i, j; };
static_assert(sizeof(Pair)==8);

std::vector<Pair> read_pairs(const char* path) {
    std::ifstream f(path,std::ios::binary|std::ios::ate);
    if (!f) throw std::runtime_error("cannot read pair file");
    auto bytes=f.tellg();
    if (bytes<0 || bytes%8) throw std::runtime_error("malformed pair file");
    std::vector<Pair> result(static_cast<size_t>(bytes)/8);
    f.seekg(0);
    if (!f.read(reinterpret_cast<char*>(result.data()),bytes))
        throw std::runtime_error("short pair file");
    return result;
}

int main(int argc,char** argv) {
  try {
    if (argc!=7) throw std::runtime_error("grid rows columns offsets changed report");
    const int64_t hh=std::stoll(argv[2]),ww=std::stoll(argv[3]);
    if (hh<=0 || ww<=0 || hh>INT32_MAX/2 || ww>INT32_MAX/2)
        throw std::runtime_error("unsupported dimensions");
    const int32_t h=hh,w=ww;
    const size_t size=static_cast<size_t>(h)*w;
    int fd=open(argv[1],O_RDONLY);
    if (fd<0) throw std::runtime_error("cannot read grid");
    struct stat st;
    if (fstat(fd,&st) || static_cast<uint64_t>(st.st_size)!=size)
        throw std::runtime_error("grid length mismatch");
    auto grid=static_cast<const uint8_t*>(mmap(nullptr,size,PROT_READ,MAP_PRIVATE,fd,0));
    if (grid==MAP_FAILED) throw std::runtime_error("cannot map grid");
    auto offsets=read_pairs(argv[4]);
    auto changed=read_pairs(argv[5]);
    for (auto q:offsets)
        if (q.i<0 || q.i>=h || q.j<0 || q.j>=w)
            throw std::runtime_error("offset outside normalized torus");
    uint64_t conflicts=0;
    size_t checked=0;
    for (auto q:changed) {
        if (q.i<0 || q.i>=h || q.j<0 || q.j>=w)
            throw std::runtime_error("changed coordinate outside grid");
        const auto color=grid[static_cast<size_t>(q.i)*w+q.j];
        if (color>=5) throw std::runtime_error("changed endpoint is not real colored");
        for (auto offset:offsets) {
            int32_t i=q.i+offset.i,j=q.j+offset.j;
            if (i>=h) i-=h;
            if (j>=w) j-=w;
            conflicts += grid[static_cast<size_t>(i)*w+j]==color;
        }
        ++checked;
        if (checked%20000==0)
            std::cout<<"native extension audit "<<checked<<"/"<<changed.size()
                     <<" conflict witnesses "<<conflicts<<std::endl;
    }
    std::ofstream report(argv[6]);
    report<<"{\"checked\":"<<checked<<",\"conflict_witnesses\":"<<conflicts<<"}\n";
    if (!report) throw std::runtime_error("cannot write report");
    munmap(const_cast<uint8_t*>(grid),size);
    close(fd);
    return 0;
  } catch (const std::exception& e) {
    std::cerr<<e.what()<<std::endl;
    return 1;
  }
}

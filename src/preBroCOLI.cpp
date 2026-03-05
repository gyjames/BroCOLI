#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <thread>
#include <numeric>
#include <cstring>
#include "WFA2-lib/bindings/cpp/WFAligner.hpp"
#include "edlib.h"
#include <array>
#include <cstdint>
#include <sys/stat.h>

using namespace std;
using Clock = std::chrono::high_resolution_clock;
using microsec = std::chrono::microseconds;  

constexpr int magicseqXKMER = 3;
constexpr int magicseqALLKMER = 9;

constexpr int XKMER10x3v3 = 7;

constexpr int visiumXKMER = 11;
constexpr int visiumALLKMER = 11;


int guagua = 0;
long long total_First = 0;
long long total_BarcodeX = 0;
long long total_BarcodeY = 0;
long long total_BarcodeZ = 0;
long long total_Final = 0;

int chushi = 0;
int erfen = 0;
int bitwei = 0;

int polyTOK = 0;
int noJG = 0;
int noAlign = 0;
int noX = 0;
int noY = 0;
int noZ = 0;
int noLen = 0;
int noFinal = 0;
int Ambiguous = 0;
int total6 = 0;
int jige = 0;

// std::mutex bigMutex;

// Append .1 to version for dev code, remove for release
// e.g. 1.00.1 (dev) goes to 1.01 (release)
const static string VERSION="1.0.0";

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  std::cerr << "usage: preBroCOLI [options] [reads_input]\n\n";
  std::cerr << "  reads_input: a .fastq or .fasta file. Will read from stdin if empty.\n\n";
  std::cerr << "  options: \n";
  std::cerr << "     -k known_list   Either 1) a text file of expected barcodes in the first column,\n";
  std::cerr << "                     one row per barcode, or 2) a comma separate string of barcodes.\n";
  std::cerr << "                     Without this option, preBroCOLI will search and report possible barcodes.\n";
  std::cerr << "                     The generated list can be used for known_list in subsequent runs.\n";
  std::cerr << "     -i true/false   Replace read ID with barcodes+UMI, remove search strings\n";
  std::cerr << "                     including flanking sequenence and split read if multiple\n";
  std::cerr << "                     barcodes found (default: true).\n";
  std::cerr << "     -s true/false   Sort reads into separate files by barcode (default: false)\n";
  std::cerr << "     -c true/false   Add a _C suffix to the read identifier of any chimeric reads\n";
  std::cerr << "                     (default: false). For instance if,\n";
  std::cerr << "                       @BC_UMI#READID_+1of2\n";
  std::cerr << "                     is chimeric, it will become:\n";
  std::cerr << "                       @BC_UMI#READID_+1of2_C\n";
  std::cerr << "     -n prefix       Prefix for output filenames.\n";
  std::cerr << "     -e N            Maximum edit distance to barcode (default 2).\n";
  std::cerr << "     -f N            Maximum edit distance to primer+polyT (default 8).\n";
  std::cerr << "     -p N            Number of threads (default: 1).\n\n";
  std::cerr << "     Specifying adaptor / barcode structure : \n";
  std::cerr << "     -x sequence Append flanking sequence to search for\n";
  std::cerr << "     -b sequence Append the barcode pattern to search for\n";
  std::cerr << "     -u sequence Append the UMI pattern to search for\n";
  std::cerr << "     Notes:\n";
  std::cerr << "          The order of these options matters\n";
  std::cerr << "          ? - can be used as a wildcard\n";
  std::cerr << "     When no search pattern x,b,u option is provided, the following default pattern is used: \n";
  std::cerr << "          primer: CTACACGACGCTCTTCCGATCT\n";
  std::cerr << "          barcode: ????????????????\n";
  std::cerr << "          UMI: ????????????\n";
  std::cerr << "          polyT: TTTTTTTTT\n";
  std::cerr << "     which is the same as providing: \n";
  std::cerr << "         -x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ???????????? -x TTTTTTTTT\n\n";

  std::cerr << "  Predefined search schemes:\n";
  std::cerr << "\n  -h  Print this usage information.\n\n";
  std::cerr << "Have a different barcode scheme you would like preBroCOLI to work with? Post a request at:\n" ;
  std::cerr << "https://github.com/gyjames/BroCOLI\n" ;

  std::cerr << "If you use preBroCOLI in your research, please cite our paper:\n" ;
  std::cerr << "Paper \n";
  std::cerr << endl;
}


bool file_exists(const std::string& path) {
  struct stat st;

  if (stat(path.c_str(), &st) != 0) {
    // std::cerr << "[ERROR] stat failed for " << path << ": " << std::strerror(errno) << std::endl;
    std::cerr << "You first run, stat failed for " << path << ": " << std::strerror(errno) << std::endl;
    return false;
  }

  if (!S_ISREG(st.st_mode)) {
    std::cerr << "[ERROR] path is not a regular file: " << path << std::endl;
    return false;
  }

  std::ifstream f(path.c_str());
  if (!f.good()) {
    std::cerr << "[ERROR] file exists but is not readable: " << path << std::endl;
    return false;
  }

  return true;
}

static const std::array<unsigned char, 256> lut = []{
    std::array<unsigned char, 256> table;
    table.fill(0xFF); // 初始化为非法值
    table['A'] = 0; table['a'] = 0;
    table['C'] = 1; table['c'] = 1;
    table['G'] = 2; table['g'] = 2;
    table['T'] = 3; table['t'] = 3;
    return table;
}();

// 将 12bp 字符串编码为 uint32_t
inline bool encode_dna_2bit_n(const char* seq, int len, uint32_t& out) {
    out = 0;
    for (int i = 0; i < len; ++i) {
        unsigned char val = lut[static_cast<unsigned char>(seq[i])];
        if (val == 0xFF) return false;
        out = (out << 2) | val;
    }
    return true;
}

struct UmiIndex {
    // 1. 数据源：不再存 string，直接引用或拷贝 uint32_t
    // 这里选择拷贝一份，因为后续为了 Cache 友好，独立的连续内存更好
    std::vector<uint32_t> valid_umis; 

    // 2. K-mer 查找表
    // LUT[kmer_int_val] = { index_in_valid_umis, ... }
    std::vector<std::vector<int>> kmer_lut;

    int k = 5;
    int umi_len = 12;
    // 解码;
    static std::string decode32(uint32_t val, int len) {
        std::string res(len, 'N');
        static const char dec[] = "ACGT";
        for (int i = len - 1; i >= 0; --i) {
            res[i] = dec[val & 3];
            val >>= 2;
        }
        return res;
    }
    // 构建函数：输入是已经编码好的整数向量
    void build(const std::vector<uint32_t>& umis, int k_val = 5, int len_val = 12) {
        k = k_val;
        umi_len = len_val;
        valid_umis = umis; // 拷贝数据 (vector copy is fast enough)

        // 初始化查找表：4^k 大小
        int num_kmers = 1 << (2 * k); 
        kmer_lut.clear();
        kmer_lut.resize(num_kmers);

        // 掩码：用于提取 k-mer (2*k 个 1)
        // 比如 k=5, mask = 1111111111 (binary) = 1023
        uint32_t mask = (1U << (2 * k)) - 1;

        // 遍历所有合法的 UMI
        for (int i = 0; i < valid_umis.size(); ++i) {
            uint32_t u_val = valid_umis[i];

            // 滑动窗口提取 K-mer (纯位运算，无 string 构造)
            // UMI 12bp = 24 bits. 高位是序列开头。
            for (int p = 0; p <= umi_len - k; ++p) {
                // 计算右移位数：
                // 比如 p=0 (取最左边), 我们要把高位移到最低位
                // 总位数 = umi_len * 2
                // 需要保留的结束位 = (umi_len - p - k) * 2
                int shift = 2 * (umi_len - p - k);
                uint32_t kmer_val = (u_val >> shift) & mask;
                // 记录：这个 k-mer 出现在第 i 个 UMI 中
                kmer_lut[kmer_val].push_back(i);
            }
        }
    }
};


void get_magicseq_inform(const std::string& whitelistPath, const std::string& barcodePathX, const std::string& barcodePathY, 
              const std::string& barcodePathZ, std::unordered_map<std::string, std::vector<uint32_t>>& WhiteList, 
              std::unordered_set<std::string>& BarcodesX, std::unordered_set<std::string>& BarcodesY, 
              std::unordered_set<std::string>& BarcodesZ) {

  if (file_exists(whitelistPath)) {
    std::ifstream file; 
    file.open(whitelistPath);
    std::cerr << "[magicseq] BarcodeX BarcodeY BarcodeZ and UMI from: " << whitelistPath << "\n";
    std::string line;
    std::string bc; bc.reserve(24);
    uint32_t umi_val = 0;

    while (getline(file, line)) {
      if (line.size() < 24) continue;             
      BarcodesX.emplace(line, 0, 8);
      BarcodesY.emplace(line, 8, 8);
      BarcodesZ.emplace(line, 16, 8);
      bc = line.substr(0, 24);
      std::vector<uint32_t>& umi_list = WhiteList[bc];
      if (line.size() > 36) {
        uint32_t umi12;
        if (encode_dna_2bit_n(line.data() + 25, 12, umi_val)) {
          umi_list.push_back(umi_val);
        }
      }
    }
    file.close();

    size_t total_umis = 0;
    for (auto& pair : WhiteList) {
      std::vector<uint32_t>& vec = pair.second;
      if (!vec.empty()) {
          std::sort(vec.begin(), vec.end());
          auto last = std::unique(vec.begin(), vec.end());
          vec.erase(last, vec.end());
          vec.shrink_to_fit(); 
      }
      total_umis += vec.size();
    }        

    std::cerr << "Number of known whitelist barcodes: " << WhiteList.size() << "\n"; 
    std::cerr << "Total unique UMIs loaded: " << total_umis << "\n";      
    if (WhiteList.empty()) {
      std::cerr << "The WhiteList file does not provide any barcodes.\n";     
      print_usage();
      exit(1);
    }

  } else {
    std::cerr << "[INFO] whitelist file does not exist: " << whitelistPath << "\n";
    if (file_exists(barcodePathX) and file_exists(barcodePathX) and file_exists(barcodePathX)) {
      std::cerr << "[INFO] All the barcode information has been provided. " << "\n";
      std::ifstream file; 
      file.open(barcodePathX);
      std::cerr << "BarcodeX from: " << barcodePathX << "\n";
      std::string line;
      while (getline(file, line)) {
        std::istringstream line_stream(line);
        std::string token;
        while (line_stream >> token) {
          if (token.size() > 7) {
            BarcodesX.insert(token);
          }
        }
      }
      file.close();
      std::cerr << "Number of known barcodesX: " << BarcodesX.size() << "\n";
      if (BarcodesX.empty()) {
        print_usage();
        exit(1);
      }

      file.clear();
      file.open(barcodePathY);
      std::cerr << "BarcodeY from: " << barcodePathY << "\n";
      while (getline(file, line)) {
        std::istringstream line_stream(line);
        std::string token;
        while (line_stream >> token) {
          if (token.size() > 7) {
            BarcodesY.insert(token);
          }
        }
      }
      file.close();
      std::cerr << "Number of known barcodesY: " << BarcodesY.size() << "\n";
      if (BarcodesY.empty()) {
        print_usage();
        exit(1); 
      }

      file.clear();
      file.open(barcodePathZ);
      std::cerr << "BarcodeZ from: " << barcodePathZ << "\n";
      while (getline(file, line)) {
        std::istringstream line_stream(line);
        std::string token;
        while (line_stream >> token) {
          if (token.size() > 7) {
            BarcodesZ.insert(token);
          }
        }
      }
      file.close();
      std::cerr << "Number of known barcodesZ: " << BarcodesZ.size() << "\n";
      if (BarcodesZ.empty()) {
        print_usage();
        exit(1);
      }
      std::cerr << "[INFO] Extracted all the X, Y, and Z candidates." << std::endl;

    } else {
      std::cerr << "[ERROR] There are no barcode informations at all." << std::endl;
      exit(1);
    }
  }
}


void get_10x3v3_inform(const std::string& whitelistPath, const std::string& barcodePathX, std::unordered_map<std::string, std::vector<uint32_t>>& WhiteList, std::unordered_set<std::string>& BarcodesX) {
  if (file_exists(whitelistPath)) {
    std::ifstream file; 
    file.open(whitelistPath);
    std::cerr << "[10x3v3] Barcode and UMI from: " << whitelistPath << "\n";
    std::string line;
    std::string bc; bc.reserve(16);
    uint32_t umi_val = 0;

    while (getline(file, line)) {
      if (line.size() < 16) continue;          
      bc = line.substr(0, 16);   
      BarcodesX.emplace(bc);
      std::vector<uint32_t>& umi_list = WhiteList[bc];
      if (line.size() > 28) {
        uint32_t umi12;
        if (encode_dna_2bit_n(line.data() + 17, 12, umi_val)) {
          umi_list.push_back(umi_val);
        }
      }
    }
    file.close();

    size_t total_umis = 0;
    for (auto& pair : WhiteList) {
      std::vector<uint32_t>& vec = pair.second;
      if (!vec.empty()) {
          std::sort(vec.begin(), vec.end());
          auto last = std::unique(vec.begin(), vec.end());
          vec.erase(last, vec.end());
          vec.shrink_to_fit(); 
      }
      total_umis += vec.size();
    }        

    std::cerr << "Number of known whitelist barcodes: " << WhiteList.size() << "\n"; 
    std::cerr << "Total unique UMIs loaded: " << total_umis << "\n";      
    if (WhiteList.empty()) {
      std::cerr << "The WhiteList file does not provide any barcodes.\n";     
      print_usage();
      exit(1);
    }
  } else {
    std::cerr << "[INFO] whitelist file does not exist: " << whitelistPath << "\n";
    if (file_exists(barcodePathX)) {
      std::cerr << "[INFO] All the barcode information has been provided. " << "\n";
      std::ifstream file; 
      file.open(barcodePathX);
      std::cerr << "BarcodeX from: " << barcodePathX << "\n";
      std::string line;
      while (getline(file, line)) {
        std::istringstream line_stream(line);
        std::string token;
        while (line_stream >> token) {
          if (token.size() > 15) {
            size_t pos = token.rfind('-');
            if (pos != std::string::npos) {
              bool all_digit = true;
              for (size_t i = pos + 1; i < token.size(); ++i) {
                  if (!isdigit(token[i])) {
                      all_digit = false;
                      break;
                  }
              }
              if (all_digit) token = token.substr(0, pos);               
            }
            if (token.size() > 15) BarcodesX.insert(token);
          }
        }
      }
      file.close();
      std::cerr << "Number of known barcodesX: " << BarcodesX.size() << "\n";
      if (BarcodesX.empty()) {
        print_usage();
        exit(1);
      }
    }
  }
}

void get_visium_inform(const std::string& whitelistPath, const std::string& barcodePathX, std::unordered_map<std::string, std::vector<uint32_t>>& WhiteList, std::unordered_set<std::string>& BarcodesX) {
  if (file_exists(whitelistPath)) {
    std::ifstream file; 
    file.open(whitelistPath);
    std::cerr << "[visium] Barcode and UMI from: " << whitelistPath << "\n";
    std::string line;
    std::string bc; bc.reserve(16);
    uint32_t umi_val = 0;

    while (getline(file, line)) {
      if (line.size() < 16) continue;          
      bc = line.substr(0, 16);   
      BarcodesX.emplace(bc);
      std::vector<uint32_t>& umi_list = WhiteList[bc];
      if (line.size() > 28) {
        uint32_t umi12;
        if (encode_dna_2bit_n(line.data() + 17, 12, umi_val)) {
          umi_list.push_back(umi_val);
        }
      }
    }
    file.close();

    size_t total_umis = 0;
    for (auto& pair : WhiteList) {
      std::vector<uint32_t>& vec = pair.second;
      if (!vec.empty()) {
          std::sort(vec.begin(), vec.end());
          auto last = std::unique(vec.begin(), vec.end());
          vec.erase(last, vec.end());
          vec.shrink_to_fit(); 
      }
      total_umis += vec.size();
    }        

    std::cerr << "Number of known whitelist barcodes: " << WhiteList.size() << "\n"; 
    std::cerr << "Total unique UMIs loaded: " << total_umis << "\n";      
    if (WhiteList.empty()) {
      std::cerr << "The WhiteList file does not provide any barcodes.\n";     
      print_usage();
      exit(1);
    }
  } else {
    std::cerr << "[INFO] whitelist file does not exist: " << whitelistPath << "\n";
    if (file_exists(barcodePathX)) {
      std::cerr << "[INFO] All the barcode information has been provided. " << "\n";
      std::ifstream file; 
      file.open(barcodePathX);
      std::cerr << "BarcodeX from: " << barcodePathX << "\n";
      std::string line;
      while (getline(file, line)) {
        std::istringstream line_stream(line);
        std::string token;
        while (line_stream >> token) {
          if (token.size() > 15) {
            size_t pos = token.rfind('-');
            if (pos != std::string::npos) {
              bool all_digit = true;
              for (size_t i = pos + 1; i < token.size(); ++i) {
                  if (!isdigit(token[i])) {
                      all_digit = false;
                      break;
                  }
              }
              if (all_digit) token = token.substr(0, pos);               
            }
            if (token.size() > 15) BarcodesX.insert(token);
          }
        }
      }
      file.close();
      std::cerr << "Number of known barcodesX: " << BarcodesX.size() << "\n";
      if (BarcodesX.empty()) {
        print_usage();
        exit(1);
      }
    }
  }
}


struct BarcodeUMIindex
{
  std::unordered_map<std::string, std::vector<std::string>> barcodeX_kmer_index;
  std::unordered_map<std::string, std::vector<std::string>> barcodeY_kmer_index;
  std::unordered_map<std::string, std::vector<std::string>> barcodeZ_kmer_index;
  std::unordered_map<std::string, std::vector<std::string>> whiteList_kmer_index;
  std::unordered_map<std::string, UmiIndex> barcode_UMI_kmer_index;
};


void get_magicseq_barcode_umi_index(std::unordered_map<std::string, std::vector<uint32_t>>& WhiteList, 
                  std::unordered_set<std::string>& BarcodesX, std::unordered_set<std::string>& BarcodesY,
                  std::unordered_set<std::string>& BarcodesZ, BarcodeUMIindex& XYZUindex) {

  std::string kmer; kmer.resize(magicseqXKMER);
  std::string bc_ext; bc_ext.resize(10);

  std::cerr << "Building [magicseq] K-mer indices for Barcode X Y Z.\n";

  for (const auto& bc:BarcodesX) {
    bc_ext.clear();
    bc_ext.append("T"); bc_ext += bc; bc_ext.append("C"); 
    for (int p = 0; p + magicseqXKMER <= bc_ext.size(); ++p) {
      std::memcpy(&kmer[0], bc_ext.data() + p, magicseqXKMER);
      XYZUindex.barcodeX_kmer_index[kmer].emplace_back(bc);
    }
  }

  for (const auto& bc:BarcodesY) {
    bc_ext.clear();
    bc_ext.append("A"); bc_ext += bc; bc_ext.append("T"); 
    for (int p = 0; p + magicseqXKMER <= bc_ext.size(); ++p) {
      std::memcpy(&kmer[0], bc_ext.data() + p, magicseqXKMER);
      XYZUindex.barcodeY_kmer_index[kmer].emplace_back(bc);
    }
  }  
  
  for (const auto& bc:BarcodesZ) {
      for (int p = 0; p + magicseqXKMER <= bc.size(); ++p) {
        kmer = bc.substr(p, magicseqXKMER);
        XYZUindex.barcodeZ_kmer_index[kmer].emplace_back(bc);
      }
  }

  kmer.resize(magicseqALLKMER);
  for (const auto& wl: WhiteList) {
    for (int p = 0; p + magicseqALLKMER <= wl.first.size(); ++p) {
      std::memcpy(&kmer[0], wl.first.data() + p, magicseqALLKMER);
      XYZUindex.whiteList_kmer_index[kmer].emplace_back(wl.first);
    }
  }
  
  std::cerr << "K-mer index build complete.\n";

  std::cerr << "Building [magicseq] UMI K-mer for " << WhiteList.size() << " barcodes...\n";
  
  for (const auto& pair : WhiteList) {
      const std::string& bc_key = pair.first;
      const std::vector<uint32_t>& umi_vec = pair.second;
      UmiIndex& idx = XYZUindex.barcode_UMI_kmer_index[bc_key];
      idx.build(umi_vec, 4, 12);
  }
  std::cerr << "UMI K-mer index build complete.\n";

}


void get_visium_barcode_umi_index(std::unordered_map<std::string, std::vector<uint32_t>>& WhiteList, 
                  std::unordered_set<std::string>& BarcodesX, BarcodeUMIindex& XYZUindex) {
  
  std::string kmer; kmer.resize(visiumXKMER);
  std::string bc_ext; bc_ext.resize(21);

  std::cerr << "Building [visium] K-mer indices for Barcode.\n";

  for (const auto& bc:BarcodesX) {
    bc_ext.clear();
    // bc_ext.append("TCTTCCGATCT");  
    bc_ext.append("CCGATCT"); 
    // bc_ext.append("TTCCGATCT");  
    // bc_ext.append("CCGATCT");       
    // bc_ext.append("GATCT"); 
    bc_ext += bc;
    for (int p = 0; p + visiumXKMER <= bc_ext.size(); ++p) {
      std::memcpy(&kmer[0], bc_ext.data() + p, visiumXKMER);
      XYZUindex.barcodeX_kmer_index[kmer].emplace_back(bc);
    }
  }

  kmer.resize(visiumALLKMER);
  for (const auto& wl: WhiteList) {
    for (int p = 0; p + visiumALLKMER <= wl.first.size(); ++p) {
      std::memcpy(&kmer[0], wl.first.data() + p, visiumALLKMER);
      XYZUindex.whiteList_kmer_index[kmer].emplace_back(wl.first);
    }
  }
  
  std::cerr << "K-mer index build complete.\n";

  std::cerr << "Building [visium] UMI K-mer for " << WhiteList.size() << " barcodes...\n";
  for (const auto& pair : WhiteList) {
      const std::string& bc_key = pair.first;
      const std::vector<uint32_t>& umi_vec = pair.second;
      UmiIndex& idx = XYZUindex.barcode_UMI_kmer_index[bc_key];
      idx.build(umi_vec, 4, 12);
  }
  std::cerr << "UMI K-mer index build complete.\n"; 

}




void get_seq_information(const std::string SeqMode, const std::string& whitelistPath, const std::string& barcodePathX, const std::string& barcodePathY, 
                const std::string& barcodePathZ, std::unordered_map<std::string, std::vector<uint32_t>>& WhiteList, 
                std::unordered_set<std::string>& BarcodesX, std::unordered_set<std::string>& BarcodesY, 
                std::unordered_set<std::string>& BarcodesZ, BarcodeUMIindex& XYZUindex) {

  if (SeqMode == "magicseq") {
    std::cerr << "Your data is [MAGIC-seq] dataset.\n";
    get_magicseq_inform(whitelistPath, barcodePathX, barcodePathY, barcodePathZ, WhiteList, BarcodesX, BarcodesY, BarcodesZ);    
    XYZUindex.whiteList_kmer_index.reserve(WhiteList.size());
    get_magicseq_barcode_umi_index(WhiteList, BarcodesX, BarcodesY, BarcodesZ, XYZUindex);
    return;

  } else if (SeqMode == "10x3v3") {
    std::cerr << "Your data is [10x3v3] dataset.\n";
    get_10x3v3_inform(whitelistPath, barcodePathX, WhiteList, BarcodesX);
    XYZUindex.whiteList_kmer_index.reserve(WhiteList.size());


  } else if (SeqMode == "visium") {
    std::cerr << "Your data is [visium] dataset.\n";
    get_visium_inform(whitelistPath, barcodePathX, WhiteList, BarcodesX);
    XYZUindex.whiteList_kmer_index.reserve(WhiteList.size());
    get_visium_barcode_umi_index(WhiteList, BarcodesX, XYZUindex);
    return;

  } else if (SeqMode == "HD") {
    std::cerr << "Your data is [HD] dataset.\n";


  } else if (SeqMode == "stereoseq") {
    std::cerr << "Your data is [stereoseq] dataset.\n";


  } else {
    std::cerr << "The -q parameter has the following types: [magicseq] [10x3v3] [visium] [HD] [stereoseq] \n";
  }
}


void get_seq_search_pattern(const std::string SeqMode, std::vector<std::pair<std::string, std::string>>& Pattern) {

  if (SeqMode == "magicseq") {
    std::cerr << "Using MAGIC-seq search pattern: \n";

    Pattern = {{"primer", "CTACACGACGCTCTTCCGATCT"},
              {"BCX", std::string(8, '?')},
              {"LinkerXY", "CAGTCATGTCATGAGCTA"},
              {"BCY", std::string(8, '?')},
              {"LinkerYZ", "TGATGCGACACTGATCGA"},
              {"BCZ", std::string(8, '?')},
              {"UMI", std::string(12, '?')},
              {"polyT", "TTTTTTTTTT"} };
    
    for (auto& i : Pattern) std::cerr << i.first << ": " << i.second << "\n";

  } else if (SeqMode == "10x3v3") {
    std::cerr << "Using 10x3v3 search pattern: \n";

    Pattern = {{"primer", "CTACACGACGCTCTTCCGATCT"},
              {"BC", std::string(16, '?')},
              {"UMI", std::string(12, '?')},
              {"polyT", "TTTTTTTTTT"} };

    for (auto& i : Pattern) std::cerr << i.first << ": " << i.second << "\n";

  } else if (SeqMode == "visium") {
    std::cerr << "Using visium search pattern: \n";

    Pattern = {{"primer", "CTACACGACGCTCTTCCGATCT"},
              {"BC", std::string(16, '?')},
              {"UMI", std::string(12, '?')},
              {"polyT", "TTTTTTTTTT"} };

    for (auto& i : Pattern) std::cerr << i.first << ": " << i.second << "\n";

  } else if (SeqMode == "HD") {
    std::cerr << "Using HD search pattern: \n";

  } else if (SeqMode == "stereoseq") {
    std::cerr << "Using stereoseq search pattern: \n";

  } else {
    std::cerr << "preBroCOLI only has the following search patterns: [magicseq] [10x3v3] [visium] [HD] [stereoseq]\n";    
  }

}



// complement nucleotides - used to reverse complement string
char complement(char& c){
  switch(c) {
  case 'A' : return 'T';
  case 'T' : return 'A';
  case 'G' : return 'C';
  case 'C' : return 'G';
  default: return 'N';
  }
}

//Inplace reverse complement
void reverse_complement(string & seq){
   reverse(seq.begin(), seq.end());
   transform(seq.begin(), seq.end(), seq.begin(), complement);
}

//Holds the found barcode and associated information
struct Barcode {
  string barcode;
  string umi;
  int editd;
  int flank_editd;
  int flank_start;
  int flank_end;
  bool unambiguous;
} ;

struct SearchResult {
  string read_id;
  string qual_scores;
  string line;
  string rev_line;
  vector<Barcode> vec_bc_for;
  vector<Barcode> vec_bc_rev;
  int count;
  bool chimeric;
};

// reads的barcode, 白名单的barcode; 0; 2
unsigned int edit_distance(const std::string& s1, const std::string& s2, unsigned int &end, int max_editd){

  std::size_t len1 = s1.size() + 1, len2 = s2.size() + 1;
  const char * s1_c = s1.c_str(); const char * s2_c = s2.c_str();

  vector<unsigned int> dist_holder(len1*len2);

  dist_holder[0]=0; //[0][0]
  for(unsigned int j = 1; j < len2; ++j) dist_holder[j] = j; //[0][j];
  for(unsigned int i = 1; i < len1; ++i) dist_holder[i*len2] = 0; //[i][0];

  int best=len2;
  end=len1-1;

  //loop over the distance matrix elements and calculate running distance
  for (unsigned int j = 1; j < len2; ++j) {
    bool any_below_threshold = false; // flag used for early exit
    for (unsigned int i = 1; i < len1; ++i) {
      int sub = (s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1; // are the bases the same?

      // if yes, no need to increment distance
      if (sub == 0) {
        dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + (j - 1)];
      } else {
        // clang-format off
        dist_holder[i*len2+j] = std::min({ //j+i*len2  //set[i][j]
          dist_holder[(i-1)*len2+j] + 1, //[i-1][j]
          dist_holder[i*len2+(j-1)] + 1, //[i][j-1]
          dist_holder[(i-1)*len2+(j-1)] + 1}); // ((s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1) });
      }
      if (dist_holder[i * len2 + j] <= max_editd) {
        any_below_threshold = true;
      }

      // if this is the last row in j
      if (j == (len2 - 1) && dist_holder[i * len2 + j] < best) {
        // check if this is the best running score
        best = dist_holder[i * len2 + j];
        end = i; // update the end position of alignment
      }
    }

    // early exit to save time.
    if(!any_below_threshold) {
      return(100);
    }
  }

  return best; // return edit distance
}

int edit_distance_ptr(const char* s1, int n1, const char* s2, int n2) {
    // 如果 n1, n2 很小 (<=12), 可以直接用栈数组做 DP 表，避免 vector<vector>
    int dp[13][13]; 
    for (int i = 0; i <= n1; i++) dp[i][0] = i;
    for (int j = 0; j <= n2; j++) dp[0][j] = j;
    for (int i = 1; i <= n1; i++) {
        for (int j = 1; j <= n2; j++) {
            int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
            dp[i][j] = std::min({
                dp[i - 1][j] + 1,       // deletion
                dp[i][j - 1] + 1,       // insertion
                dp[i - 1][j - 1] + cost // substitution
            });
        }
    }
    return dp[n1][n2];
}

// 在已经匹配到barcode的位置后, 从测序读段read中提取UMI序列;
// extract UMI from the read after barcode matching
std::string get_umi(const std::string &seq,
                    const std::vector<std::pair<std::string, std::string>> &search_pattern,
                    const std::vector<int> &read_to_subpatterns,
                    const int umi_index, UmiIndex &umi_map, 
                    const std::string &SEQID) { 
  if (umi_index < 0) return "";   
  const size_t umi_start = read_to_subpatterns[umi_index];
  const size_t umi_length = search_pattern[umi_index].second.size();

  std::string UMI_init; UMI_init.reserve(umi_length);  
  if (umi_start >= seq.size()) {
      return std::string(umi_length, 'N');
  }

  std::string umi_pad = "";
  if (seq.size() < umi_start + umi_length) {
      umi_pad = std::string(search_pattern[umi_index].second.length() - umi_length, 'N');
  }

  UMI_init = seq.substr(umi_start, umi_length) + umi_pad;

  // return UMI_init;

  if (umi_map.valid_umis.empty()) {
      return UMI_init;
  }  
  
  uint32_t exact_match_val = 0;
  bool can_encode = true;

  // 1. 尝试将字符串转为整数 (复用全局 lut 表)
  // 这里手写循环是为了极致性能，避免函数调用开销
  for (int i = 0; i < umi_length; ++i) {
      unsigned char base = lut[static_cast<uint8_t>(UMI_init[i])];
      if (base == 0xFF) { 
          can_encode = false; // 含有 N, 无法精确匹配
          break; 
      }
      exact_match_val = (exact_match_val << 2) | base;
  }

  // 2. 二分查找
  if (can_encode) {
      // std::binary_search 依赖于 vector 是有序的
      if (std::binary_search(umi_map.valid_umis.begin(), 
                             umi_map.valid_umis.end(), 
                             exact_match_val)) {
          // std::cerr << "直接在里面:" << UMI_init << std::endl;
          // std::cerr << std::endl;
          erfen = erfen + 1;
          return UMI_init; 
      }
  }

  static thread_local std::vector<int> votes;
  if (votes.size() < umi_map.valid_umis.size()) {
      votes.resize(umi_map.valid_umis.size());
  }  
  // 只清零我们需要用到的部分,节省时间;
  std::memset(votes.data(), 0, umi_map.valid_umis.size() * sizeof(int));
  int k = umi_map.k; 
  int len = UMI_init.size();
  if (len < k) return UMI_init;
  // --- 3. 投票过程 (Zero-Copy) ---;
  for (int p = 0; p <= len - k; ++p) {
      uint32_t kmer_val = 0;
      bool valid_kmer = true;
      for (int j = 0; j < k; ++j) {
          // 注意 1：必须强转为 unsigned char 也就是 uint8_t 作为下标，防止负数下标越界
          unsigned char base = lut[static_cast<uint8_t>(UMI_init[p + j])];
          // 注意 2：你的非法值是 0xFF，所以这里不能判断 < 0，而要判断是否等于 0xFF
          if (base == 0xFF) { 
              valid_kmer = false; 
              break; // 遇到 N 或非法字符，跳过这个 k-mer
          } 
          kmer_val = (kmer_val << 2) | base;
      }
      if (valid_kmer) {
          // 如果 k-mer 合法，再去查索引
          // 注意：需确保 kmer_val 不会超出 lut 范围，但在 2-bit 编码下，
          // 5-mer 最大值是 1023，数组肯定是够大的
          if (kmer_val < umi_map.kmer_lut.size()) {
              const auto& candidates = umi_map.kmer_lut[kmer_val];
              for (int candidate_id : candidates) {
                  votes[candidate_id]++;
              }
          }
      }
  }
  // --- 4. 寻找最佳匹配 ---
  int best_id = -1;
  int min_dist = 999;
  int threshold_dist = 3;

  for (int i = 0; i < umi_map.valid_umis.size(); ++i) {
      if (votes[i] > 0) {
          // 1. 在栈上开辟一个小缓存，不涉及任何内存分配
          char buffer[16]; 
          
          // 2. 现场解码 (On-the-fly decoding)
          uint32_t val = umi_map.valid_umis[i];
          static const char dec[] = "ACGT";
          // 从后往前解，或者根据你的编码逻辑从前往后
          // 假设高位是字符串开头：
          for (int j = 0; j < umi_length; ++j) {
              int shift = 2 * (umi_length - 1 - j);
              buffer[j] = dec[(val >> shift) & 3];
          }
          // 注意：buffer 不是以 \0 结尾的 C-string，处理时要注意长度
          
          // 3. 计算编辑距离
          // 这里你需要一个支持 (char*, int, char*, int) 参数的 edit_distance 函数
          // 或者简单的把 string 的 data() 传进去
          int d = edit_distance_ptr(UMI_init.data(), UMI_init.size(), buffer, umi_length);
          
          if (d < min_dist) {
              min_dist = d;
              best_id = i;
          }
      }
  }

  // --- 5. 返回结果 ---
  if (best_id != -1 && min_dist <= threshold_dist) {
      // std::cerr << "解码得到UMI_final:" << UmiIndex::decode32(umi_map.valid_umis[best_id], umi_length) << std::endl;
      // std::cerr << std::endl;
      bitwei = bitwei + 1;
      return UmiIndex::decode32(umi_map.valid_umis[best_id], umi_length);
  }
  // std::cerr << "readID:" << SEQID << std::endl;
  // std::cerr << "UMI_init:" << UMI_init << "  ";
  // std::cerr << "还是用的初始的:" << UMI_init << std::endl;
  // std::cerr << "序列:" << seq.substr(0,120) << std::endl;
  // std::cerr << std::endl;

  chushi = chushi + 1;
  return UMI_init;
}


size_t find_polyT_start_allow_errors(const std::string& s, size_t min_T = 10, size_t max_nonT = 1) {
    // l:当前窗口左端点; cntT:窗口内T的个数; cntNonT:窗口内非T的个数;
    size_t l = 0;
    size_t cntT = 0, cntNonT = 0;
    // bestL:目前找到的最佳poly-T起点; 最佳窗口里的T数量;
    size_t bestL = std::string::npos;
    size_t bestT = 0;
    // r是窗口右端点; 
    for (size_t r = 0; r < s.size(); ++r) {
        // 更新窗口内的计数;
        if (s[r] == 'T') ++cntT;
        else ++cntNonT;
        // 如果非T超标, 就移动左端点缩窗口;
        while (cntNonT > max_nonT) {
            if (s[l] == 'T') --cntT;
            else --cntNonT;
            ++l;
        }
        // 如果非T超标, 就移动左端点缩窗口;
        size_t l2 = l;
        size_t cntT2 = cntT;
        size_t cntNonT2 = cntNonT;
        while (l2 <= r && s[l2] != 'T') {
            if (s[l2] == 'T') --cntT2;
            else --cntNonT2;
            ++l2;
        }
        // 判断这个"起点是T的候选窗口"是否合格, 并更新最优;
        if (l2 <= r && cntNonT2 <= max_nonT && cntT2 >= min_T) {
            if (cntT2 > bestT || (cntT2 == bestT && l2 < bestL)) {
                bestT = cntT2;
                bestL = l2;
            }
        }
    }
    if (bestL == std::string::npos) return 0;
    else return bestL;
}



Barcode get_magicseq_barcode(const string& readID, string & seq, unordered_set<string> *known_barcodesX,
        unordered_set<string> *known_barcodesY, unordered_set<string> *known_barcodesZ,
		    int flank_max_editd,
        const std::vector<std::pair<std::string, std::string>> &search_pattern,
        std::unordered_map<std::string, std::vector<uint32_t>> *WHITELIST,
        BarcodeUMIindex& buindex) {

  Clock::time_point t1, t2;

  Barcode barcode;
  barcode.editd=-1000; barcode.flank_editd=100; barcode.unambiguous = false;

  EdlibEqualityPair additionalEqualities[32] = {
    {'R', 'A'}, {'R', 'G'},
    {'K', 'G'}, {'K', 'T'},
    {'S', 'G'}, {'S', 'C'},
    {'Y', 'C'}, {'Y', 'T'},
    {'M', 'A'}, {'M', 'C'},
    {'W', 'A'}, {'W', 'T'},
    {'B', 'C'}, {'B', 'G'}, {'B', 'T'},
    {'H', 'A'}, {'H', 'C'}, {'H', 'T'},
    {'?', 'A'}, {'?', 'C'}, {'?', 'G'}, {'?', 'T'},
    {'D', 'A'}, {'D', 'G'}, {'D', 'T'},
    {'V', 'A'}, {'V', 'C'}, {'V', 'G'},
    {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
  };
  // 初始化edlib比对配置; EDLIB_MODE_HW 半全局比对模式(适合找子串); EDLIB_TASK_PATH 要求返回路径, 以便之后推算子模式在read中的位置; IUPAC模糊匹配表; 表示表中共有28对匹配规则;
  EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 32};

  // 将所有的序列拼接起来;
  std::string search_string; search_string.reserve(105);
  std::vector<long unsigned int> subpattern_lengths; subpattern_lengths.reserve(search_pattern.size());
  for (const auto &pair : search_pattern) {
    search_string += pair.second;
    subpattern_lengths.push_back(pair.second.length());
  }

  // 如果read太短就跳过这个read, 返回一个空的Barcode structure;
  if (seq.length() < search_string.length()) {
    noLen = noLen + 1;
    return barcode;
  }

  t1 = Clock::now();
  EdlibAlignResult result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
  t2 = Clock::now();
  total_First += std::chrono::duration_cast<microsec>(t2 - t1).count();

  if (result.status != EDLIB_STATUS_OK | result.numLocations == 0) {
    noAlign = noAlign + 1;
    edlibFreeAlignResult(result);
    return(barcode); 
  } 

  barcode.flank_editd = result.editDistance;
  barcode.flank_start = result.startLocations[0];
  barcode.flank_end = result.endLocations[0];

  // 依次计算累加的总长度; std::partial_sum是 C++ <numeric> 头文件里的一个标准算法; 把输入序列中每个元素依次累加, 输出累加到当前位置的总和;
  std::vector<long unsigned int> subpattern_ends;
  subpattern_ends.resize(subpattern_lengths.size());
  std::partial_sum(subpattern_lengths.begin(), subpattern_lengths.end(), subpattern_ends.begin());

  // 记录每个子pattern (primer/BC/UMI/polyA) 在read中的起止坐标; 首先只有序列刚开始的位置;
  vector<int> read_to_subpatterns;
  read_to_subpatterns.reserve(subpattern_ends.size() + 1);
  read_to_subpatterns.emplace_back(barcode.flank_start); 

  int i_read = barcode.flank_start;
  int i_pattern = 0;
  int i_subpattern = 0;

  // walk through edlib aligment; 现在target是每条reads的序列, query是搜索模式拼接序列;
  // 0 for match; 1 for insertion to target; 2 for insertion to query; 3 for mismatch;
  std::vector<unsigned char> alignment_vector(result.alignment, result.alignment + result.alignmentLength);
  for (const auto& value : alignment_vector) {
    if (value != 1) i_read++;
    if (value != 2) i_pattern++;
    if (i_pattern >= subpattern_ends[i_subpattern]) {
      read_to_subpatterns.emplace_back(i_read);
      i_subpattern++;
    }
  }

  edlibFreeAlignResult(result);

  int bcx_index = 1;
  int bcy_index = 3;
  int bcz_index = 5;
  int umi_index = 6;

  int edit_XY = 10; // 初始为10;
  int edit_Z = edit_XY - 2; // 初始为8;

  // 从reads中barcode开始的位置提取; 提取barcode长度的序列;
  std::string exact_umi;
  std::string BarcodeX, BarcodeY, BarcodeZ;
  std::vector<std::string> BarcodeVecStrX, BarcodeVecStrY, BarcodeVecStrZ;
  
  t1 = Clock::now();
  // barcode X 
  std::string exact_bcx = seq.substr(read_to_subpatterns[bcx_index], read_to_subpatterns[bcx_index+1]-read_to_subpatterns[bcx_index]);
  int barcodeX_editd = 100;
  bool barcodeX_unambiguous = false;
  if (known_barcodesX->size()==0 || (known_barcodesX->find(exact_bcx) != known_barcodesX->end())) {
    barcodeX_editd = 0;
    barcodeX_unambiguous = true;
    BarcodeX = exact_bcx;
    BarcodeVecStrX.push_back(BarcodeX);
  } else {
    unordered_map<string, int> candidate_votes;  
    std::string bcx_ext;
    bcx_ext.reserve(exact_bcx.size()+2);  
    bcx_ext.append("T");          
    bcx_ext += exact_bcx;            
    bcx_ext.append("C");  

    for (int p = 0; p + magicseqXKMER <= bcx_ext.size(); ++p) {
        string kmer = bcx_ext.substr(p, magicseqXKMER);
        auto it = buindex.barcodeX_kmer_index.find(kmer);
        if (it != buindex.barcodeX_kmer_index.end()) {
            for (string id : it->second)
                candidate_votes[id]++;
        }
    }
    vector<pair<string,int>> cands(candidate_votes.begin(), candidate_votes.end());
    sort(cands.begin(), cands.end(), [](const std::pair<string,int>& a, const std::pair<string,int>& b){
        return a.second > b.second;
    });
    if (cands.size() > 8) cands.resize(8);

    int left_bound = max(read_to_subpatterns[bcx_index-1], 0);
    int max_length = read_to_subpatterns[bcx_index+2] - read_to_subpatterns[bcx_index-1];
    std::string barcode_seq = seq.substr(left_bound, max_length);

    unsigned int editDistance, endDistance;
    search_string.reserve(60);
    for (const auto& cands_barcode:cands) {
      search_string.clear();
      search_string.append(search_pattern[bcx_index-1].second); search_string.append(cands_barcode.first); search_string.append(search_pattern[bcx_index+1].second);
      editDistance = edit_distance(barcode_seq, search_string, endDistance, edit_XY);

      if (editDistance == barcodeX_editd) {
        barcodeX_unambiguous = false;
        if (editDistance < edit_XY) {
          BarcodeVecStrX.push_back(cands_barcode.first);
        }
      } else if (editDistance < barcodeX_editd && editDistance < edit_XY) {
        BarcodeVecStrX.clear();
        BarcodeVecStrX.push_back(cands_barcode.first);
        barcodeX_unambiguous = true;
        barcodeX_editd = editDistance; 
        BarcodeX = cands_barcode.first;
      }
    }

    if (BarcodeVecStrX.size() == 0) {
      left_bound = max(read_to_subpatterns[bcx_index], 0);
      max_length = read_to_subpatterns[bcx_index+1] - read_to_subpatterns[bcx_index];
      barcode_seq = seq.substr(left_bound, max_length);
      for (const auto& cands_barcode:cands) {
        search_string.clear();
        search_string.append(cands_barcode.first);
        editDistance = edit_distance(barcode_seq, search_string, endDistance, 2);
        if (editDistance < 3) {
          BarcodeVecStrX.push_back(cands_barcode.first);
          barcodeX_unambiguous = true;
          barcodeX_editd = editDistance; 
          BarcodeX = cands_barcode.first;
        }
      }
    }
  }

  t2 = Clock::now();
  total_BarcodeX += std::chrono::duration_cast<microsec>(t2 - t1).count();

  if (BarcodeVecStrX.size() == 0 and BarcodeX.size() == 0) {
    noX = noX + 1;
    return barcode;
  }
  if (BarcodeVecStrX.size() == 1) exact_bcx = BarcodeVecStrX[0];
  // exact_bcx = BarcodeVecStrX[0];

  // barcode Y;
  t1 = Clock::now();
  BarcodeVecStrY.clear();
  // std::string exact_bcy = seq.substr(read_to_subpatterns[bcy_index], search_pattern[bcy_index].second.length());
  std::string exact_bcy = seq.substr(read_to_subpatterns[bcy_index], read_to_subpatterns[bcy_index+1]-read_to_subpatterns[bcy_index]);
  int barcodeY_editd = 100;
  bool barcodeY_unambiguous = false;
  if (known_barcodesY->size()==0 || (known_barcodesY->find(exact_bcy) != known_barcodesY->end())) {
    barcodeY_editd = 0;
    barcodeY_unambiguous = true;
    BarcodeY = exact_bcy;
    BarcodeVecStrY.push_back(BarcodeY);
  } else {
    unordered_map<string, int> candidate_votes;  
    std::string bcy_ext;
    bcy_ext.reserve(exact_bcy.size()+2);  
    bcy_ext.append("A");          
    bcy_ext += exact_bcy;            
    bcy_ext.append("T");  

    for (int p = 0; p + magicseqXKMER <= bcy_ext.size(); ++p) {
        string kmer = bcy_ext.substr(p, magicseqXKMER);
        auto it = buindex.barcodeY_kmer_index.find(kmer);
        if (it != buindex.barcodeY_kmer_index.end()) {
            for (string id : it->second)
                candidate_votes[id]++;
        }
    }

    vector<pair<string,int>> cands(candidate_votes.begin(), candidate_votes.end());
    sort(cands.begin(), cands.end(), [](const std::pair<string,int>& a, const std::pair<string,int>& b){
        return a.second > b.second;
    });
    if (cands.size() > 8) cands.resize(8);

    // linkerXY+barcodeY+linkerYZ
    int left_bound = max(read_to_subpatterns[bcy_index-1], 0);
    int max_length = read_to_subpatterns[bcy_index+2] - read_to_subpatterns[bcy_index-1];
    std::string barcode_seq = seq.substr(left_bound, max_length); 

    // 遍历所有的白名单中的barcode;
    unsigned int editDistance, endDistance;
    search_string.reserve(60);
    for (const auto& cands_barcode:cands) {
      search_string.clear();
      search_string.append(search_pattern[bcy_index-1].second); search_string.append(cands_barcode.first); search_string.append(search_pattern[bcy_index+1].second);
      editDistance = edit_distance(barcode_seq, search_string, endDistance, edit_XY);
      if (editDistance == barcodeY_editd) {
        barcodeY_unambiguous = false;
        if (editDistance < edit_XY) {
          BarcodeVecStrY.push_back(cands_barcode.first);
        }
      } else if (editDistance < barcodeY_editd && editDistance < edit_XY) {
        BarcodeVecStrY.clear();
        BarcodeVecStrY.push_back(cands_barcode.first);
        barcodeY_unambiguous = true;
        barcodeY_editd = editDistance; 
        BarcodeY = cands_barcode.first;
      }
    }

    if (BarcodeVecStrY.size() == 0) {
      left_bound = max(read_to_subpatterns[bcy_index], 0);
      max_length = read_to_subpatterns[bcy_index+1] - read_to_subpatterns[bcy_index];
      barcode_seq = seq.substr(left_bound, max_length);
      for (const auto& cands_barcode:cands) {
        search_string.clear();
        search_string.append(cands_barcode.first);
        editDistance = edit_distance(barcode_seq, search_string, endDistance, 2);
        if (editDistance < 3) {
          BarcodeVecStrY.push_back(cands_barcode.first);
          barcodeY_unambiguous = true;
          barcodeY_editd = editDistance; 
          BarcodeY = cands_barcode.first;
        }
      }
    }
  }
  t2 = Clock::now();
  total_BarcodeY += std::chrono::duration_cast<microsec>(t2 - t1).count();

  if (BarcodeVecStrY.size() == 0 and BarcodeY.size() == 0) {
    noY = noY + 1;
    return barcode;
  }
  if (BarcodeVecStrY.size() == 1) exact_bcy = BarcodeVecStrY[0];
  // exact_bcy = BarcodeVecStrY[0];

  // barcode Z;
  t1 = Clock::now();
  BarcodeVecStrZ.clear();
  // std::string exact_bcz = seq.substr(read_to_subpatterns[bcz_index], search_pattern[bcz_index].second.length());
  std::string exact_bcz = seq.substr(read_to_subpatterns[bcz_index], read_to_subpatterns[bcz_index+1]-read_to_subpatterns[bcz_index]);
  int barcodeZ_editd = 100;
  bool barcodeZ_unambiguous = false;
  if (known_barcodesZ->size()==0 || (known_barcodesZ->find(exact_bcz) != known_barcodesZ->end())){
    barcodeZ_editd = 0;
    barcodeZ_unambiguous = true;
    BarcodeZ = exact_bcz;
    BarcodeVecStrZ.push_back(BarcodeZ);
  } else {
    // linkYZ+barcodeZ+umi+ployT;
    int left_bound = max(read_to_subpatterns[bcz_index-1], 0);
    int max_length = read_to_subpatterns[bcz_index+1] - read_to_subpatterns[bcz_index-1];
    std::string barcode_seq = seq.substr(left_bound, max_length);
    // 遍历所有的白名单中的barcode;
    unordered_set<string>::iterator known_barcodes_itr = known_barcodesZ->begin();
    unsigned int editDistance, endDistance;
    search_string.reserve(60);
    for(; known_barcodes_itr != known_barcodesZ->end(); known_barcodes_itr++){
      search_string.clear();
      search_string.append(search_pattern[bcz_index-1].second); search_string.append(*known_barcodes_itr);
      editDistance = edit_distance(barcode_seq, search_string, endDistance, edit_Z);
      if (editDistance == barcodeZ_editd) {
        barcodeZ_unambiguous = false;
        if (editDistance < edit_Z) {
           BarcodeVecStrZ.push_back(*known_barcodes_itr);
        }        
      } else if (editDistance < barcodeZ_editd && editDistance < edit_Z) {
        BarcodeVecStrZ.clear();
        BarcodeVecStrZ.push_back(*known_barcodes_itr);        
        barcodeZ_unambiguous = true;
        barcodeZ_editd = editDistance;
        BarcodeZ = *known_barcodes_itr;
      }
    }

    if (BarcodeVecStrZ.size() == 0) {
      left_bound = max(read_to_subpatterns[bcz_index], 0);
      max_length = read_to_subpatterns[bcz_index+1] - read_to_subpatterns[bcz_index];
      barcode_seq = seq.substr(left_bound, max_length);
      for(; known_barcodes_itr != known_barcodesZ->end(); known_barcodes_itr++) {
        search_string.clear();
        search_string.append(*known_barcodes_itr);
        editDistance = edit_distance(barcode_seq, search_string, endDistance, 2);
        if (editDistance < 3) {
          BarcodeVecStrZ.push_back(*known_barcodes_itr);
          barcodeZ_unambiguous = true;
          barcodeZ_editd = editDistance; 
          BarcodeZ = *known_barcodes_itr;
        }
      }      
    }
  }
  t2 = Clock::now();
  total_BarcodeZ += std::chrono::duration_cast<microsec>(t2 - t1).count();

  if (BarcodeVecStrZ.size() == 0 and BarcodeZ.size() == 0) {
    noZ = noZ + 1;
    return barcode;
  }  
  if (BarcodeVecStrZ.size() == 1) exact_bcz = BarcodeVecStrZ[0];
  // exact_bcz = BarcodeVecStrZ[0];

  search_string.reserve(100);
  std::string Final_Barcode; Final_Barcode.reserve(24);
  t1 = Clock::now();
  if (WHITELIST->size() > 0) {
    std::string whiteListString;
    whiteListString.reserve(24);
    std::unordered_set<std::string> Final24;
    for (const auto& x:BarcodeVecStrX) {
      for (const auto& y:BarcodeVecStrY) {
        for (const auto& z:BarcodeVecStrZ) {
          Final24.clear();
          whiteListString.clear();
          whiteListString.append(x); whiteListString.append(y); whiteListString.append(z);
          for (int p = 0; p + magicseqALLKMER <= whiteListString.size(); ++p) {
            string kmer = whiteListString.substr(p, magicseqALLKMER);
            auto it = buindex.whiteList_kmer_index.find(kmer);
            if (it != buindex.whiteList_kmer_index.end()) {
              for (string id : it->second) {
                Final24.insert(id);
              }
            }
          }
        }
      }
    }

    wfa::WFAlignerGapAffine aligner(4,6,2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    std::string linshi;
    int left_bound = max(read_to_subpatterns[bcx_index-1], 0);
    int max_length = read_to_subpatterns[bcz_index+1] - read_to_subpatterns[bcx_index-1];
    std::string barcode_seq = seq.substr(left_bound, max_length);

    for(const auto& fin:Final24) {
      search_string.clear();
      linshi = fin.substr(0, 8);
      search_string.append(search_pattern[0].second); search_string.append(linshi); 
      linshi = fin.substr(8, 8);
      search_string.append(search_pattern[2].second); search_string.append(linshi); 
      linshi = fin.substr(16, 8);
      search_string.append(search_pattern[4].second); search_string.append(linshi);
      aligner.alignEnd2End(search_string, barcode_seq);
      if (aligner.getAlignmentScore() > barcode.editd) {
        Final_Barcode = fin;
        barcode.editd = aligner.getAlignmentScore();
      }
    }
  } else {
    for (const auto& x:BarcodeVecStrX) {
      for (const auto& y:BarcodeVecStrY) {
        for (const auto& z:BarcodeVecStrZ) {
          search_string.clear();
          search_string.append(search_pattern[0].second); search_string.append(x);
          search_string.append(search_pattern[2].second); search_string.append(y);
          search_string.append(search_pattern[4].second); search_string.append(z);
          search_string.append("NNNNNNNNNNNN"); search_string.append("TTTTTTTTTT");
          EdlibAlignResult result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
          if (result.editDistance < barcode.editd) {
            BarcodeX = x; BarcodeY = y; BarcodeZ = z;
            barcode.editd = result.editDistance;
          }
        }
      }
    }
  }
  t2 = Clock::now();
  total_Final += std::chrono::duration_cast<microsec>(t2 - t1).count();

  if (!Final_Barcode.empty()) {
    search_string.clear();
    search_string.append(exact_bcx); search_string.append(exact_bcy); search_string.append(exact_bcz);
    unsigned int endDistance;
    barcode.editd = edit_distance(Final_Barcode, search_string, endDistance, 8);
    exact_umi = get_umi(seq, search_pattern, read_to_subpatterns, umi_index, buindex.barcode_UMI_kmer_index.at(Final_Barcode), readID);
    barcode.unambiguous = true;
    barcode.umi = exact_umi;
    barcode.barcode = Final_Barcode;

    if (barcode.editd > 6) {
      jige = jige + 1;
    }

  } else {
    noFinal = noFinal+1;
  }

  return barcode;
}



Barcode get_visium_barcode(const std::string& readID, std::string& seq, unordered_set<string> *known_barcodesX,
		    int flank_max_editd, const std::vector<std::pair<std::string, std::string>> &search_pattern,
        std::unordered_map<std::string, std::vector<uint32_t>> *WHITELIST,
        BarcodeUMIindex& buindex) {
  
  Clock::time_point t1, t2;

  Barcode barcode;
  barcode.editd=-1000; barcode.flank_editd=100; barcode.unambiguous = false;

  EdlibEqualityPair additionalEqualities[32] = {
    {'R', 'A'}, {'R', 'G'},
    {'K', 'G'}, {'K', 'T'},
    {'S', 'G'}, {'S', 'C'},
    {'Y', 'C'}, {'Y', 'T'},
    {'M', 'A'}, {'M', 'C'},
    {'W', 'A'}, {'W', 'T'},
    {'B', 'C'}, {'B', 'G'}, {'B', 'T'},
    {'H', 'A'}, {'H', 'C'}, {'H', 'T'},
    {'?', 'A'}, {'?', 'C'}, {'?', 'G'}, {'?', 'T'},
    {'D', 'A'}, {'D', 'G'}, {'D', 'T'},
    {'V', 'A'}, {'V', 'C'}, {'V', 'G'},
    {'N', 'A'}, {'N', 'C'}, {'N', 'G'}, {'N', 'T'}
  };
  // 初始化edlib比对配置; EDLIB_MODE_HW 半全局比对模式;
  EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 32};

  // 将所有的序列拼接起来;
  std::string search_string; search_string.reserve(60);
  std::vector<long unsigned int> subpattern_lengths; subpattern_lengths.reserve(search_pattern.size());
  for (const auto &pair : search_pattern) {
    search_string += pair.second;
    subpattern_lengths.push_back(pair.second.length());
  }

  if (seq.length() < search_string.length()) return barcode;

  t1 = Clock::now();
  EdlibAlignResult result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
  t2 = Clock::now();
  total_First += std::chrono::duration_cast<microsec>(t2 - t1).count();

  if (result.status != EDLIB_STATUS_OK | result.numLocations == 0) {
    noAlign = noAlign + 1;
    edlibFreeAlignResult(result);
    return(barcode); 
  } 

  barcode.flank_editd = result.editDistance;
  barcode.flank_start = result.startLocations[0];
  barcode.flank_end = result.endLocations[0];

  std::vector<long unsigned int> subpattern_ends; subpattern_ends.resize(subpattern_lengths.size());
  std::partial_sum(subpattern_lengths.begin(), subpattern_lengths.end(), subpattern_ends.begin());

  vector<int> read_to_subpatterns; 
  read_to_subpatterns.reserve(subpattern_ends.size() + 1);
  read_to_subpatterns.emplace_back(barcode.flank_start); 

  int i_read = barcode.flank_start;
  int i_pattern = 0;
  int i_subpattern = 0;

  std::vector<unsigned char> alignment_vector(result.alignment, result.alignment + result.alignmentLength);
  for (const auto& value : alignment_vector) {
    if (value != 1) {
      i_read++;
    }
    if (value != 2) {
      i_pattern++;
    }
    if (i_pattern >= subpattern_ends[i_subpattern]) {
      read_to_subpatterns.emplace_back(i_read);
      i_subpattern++;
    }
  }

  edlibFreeAlignResult(result);

  int bc_index = 1;
  int umi_index = 2;
  int edit_X = 5;

  std::string exact_umi;
  std::string BarcodeX;
  std::vector<std::string> BarcodeVecStrX;
  
  t1 = Clock::now();
  // barcode
  std::string exact_bcx = seq.substr(read_to_subpatterns[bc_index], read_to_subpatterns[bc_index+1]-read_to_subpatterns[bc_index]);
  int barcodeX_editd = 100;
  bool barcodeX_unambiguous = false;
  if (known_barcodesX->size()==0 || (known_barcodesX->find(exact_bcx) != known_barcodesX->end())) {
    barcodeX_editd = 0;
    barcodeX_unambiguous = true;
    BarcodeX = exact_bcx;
    BarcodeVecStrX.push_back(BarcodeX);
  } else {
    unordered_map<string, int> candidate_votes;  
    std::string bcx_ext;
    bcx_ext.reserve(exact_bcx.size() + 15);
    bcx_ext.append("CCGATCT"); //"TCTTCCGATCT"
    bcx_ext += exact_bcx;

    for (int p = 0; p + visiumXKMER <= bcx_ext.size(); ++p) {
        std::string kmer = bcx_ext.substr(p, visiumXKMER);
        auto it = buindex.barcodeX_kmer_index.find(kmer);
        if (it != buindex.barcodeX_kmer_index.end()) {
            for (std::string id : it->second) candidate_votes[id]++;
        }
    }

    std::vector<std::pair<std::string,int>> cands(candidate_votes.begin(), candidate_votes.end());
    sort(cands.begin(), cands.end(), [](const std::pair<string,int>& a, const std::pair<std::string,int>& b){
        return a.second > b.second;
    });

    if (cands.size() > 10) cands.resize(10);

    for (const auto& aaa:cands) {
      BarcodeVecStrX.push_back(aaa.first);
    }
  }

  t2 = Clock::now();
  total_BarcodeX += std::chrono::duration_cast<microsec>(t2 - t1).count();

  int left_bound = max(read_to_subpatterns[bc_index] - 5, 0);
  int max_length = search_pattern[bc_index].second.length() + 2 * 5;
  std::string barcode_seq = seq.substr(left_bound, max_length);

  if (BarcodeVecStrX.empty() and BarcodeX.size() == 0) {
    noX = noX + 1;
    unordered_set<string>::iterator known_barcodes_itr = known_barcodesX->begin();
    unsigned int editDistance, endDistance;
    for(; known_barcodes_itr != known_barcodesX->end(); known_barcodes_itr++) {
      search_string = *known_barcodes_itr;
      editDistance = edit_distance(barcode_seq, search_string, endDistance, 3);
      if (editDistance == barcode.editd) {
        barcode.unambiguous = false;
      } else if (editDistance < barcode.editd && editDistance <= 3) {
        BarcodeVecStrX.clear();
        barcode.unambiguous = true;
        barcode.editd = editDistance;
        barcode.barcode = *known_barcodes_itr;
        BarcodeVecStrX.push_back(barcode.barcode);
        if (editDistance == 0) break;
      }
    } 
  }

  if (BarcodeVecStrX.size() == 1) exact_bcx = BarcodeVecStrX[0];  
  
  t1 = Clock::now();

  std::string Final_Barcode; Final_Barcode.reserve(16);
  if (WHITELIST->size() > 0) {
    wfa::WFAlignerGapAffine aligner(4,5,1, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    for(const auto& fin:BarcodeVecStrX) {
      aligner.alignEndsFree(fin, 0, 0, barcode_seq, 1, 1);
      if (aligner.getAlignmentScore() > barcode.editd) {
        Final_Barcode = fin;
        barcode.editd = aligner.getAlignmentScore();
      }
    }

  }
  t2 = Clock::now();
  total_Final += std::chrono::duration_cast<microsec>(t2 - t1).count();

  if (!Final_Barcode.empty()) {
    unsigned int endDistance;
    barcode.editd = edit_distance(barcode_seq, Final_Barcode, endDistance, edit_X);
    exact_umi = get_umi(seq, search_pattern, read_to_subpatterns, umi_index, buindex.barcode_UMI_kmer_index.at(Final_Barcode), readID);
    barcode.unambiguous = true;
    barcode.umi = exact_umi;
    barcode.barcode = Final_Barcode;
  }

  // std::cerr << "barcode.barcode:" << barcode.barcode << std::endl;
  return barcode;
}


// 用来在一条测序读段sequence中递归地查找所有可能存在的barcode;
std::vector<Barcode> big_barcode_search(const std::string& ModeSeq, const string& seqID, string& sequence, 
                                  unordered_set<string>& known_barcodes_X, 
                                  unordered_set<string>& known_barcodes_Y, 
                                  unordered_set<string>& known_barcodes_Z,
                                  int max_flank_editd,
                                  const std::vector<std::pair<std::string, std::string>> &search_pattern,
                                  std::unordered_map<std::string, std::vector<uint32_t>>& WHITELIST,
                                  BarcodeUMIindex& buindex) {

  std::vector<Barcode> return_vec;

  Barcode result;

  if (ModeSeq == "magicseq") {
    result = get_magicseq_barcode(seqID, sequence, &known_barcodes_X, &known_barcodes_Y, &known_barcodes_Z,
                                  max_flank_editd, search_pattern, &WHITELIST, buindex);
    if (result.editd <= 2 && result.unambiguous) return_vec.push_back(result);

    if (!result.unambiguous) Ambiguous = Ambiguous + 1;

    if(result.editd > 2) total6 = total6 + 1;

  } else if (ModeSeq == "10x3v3") {


  } else if (ModeSeq == "visium") {
    result = get_visium_barcode(seqID, sequence, &known_barcodes_X,
                                  max_flank_editd, search_pattern, &WHITELIST, buindex);

    if (result.editd <= 2 && result.unambiguous) return_vec.push_back(result);

  }

  // 如果已经找到至少一个barcode, 则继续在剩下的序列里搜索其他可能的barcode;
  if(return_vec.size() > 0) {
    string masked_sequence = sequence;
    for(int i = 0; i < return_vec.size(); i++) {
      int flank_length = return_vec.at(i).flank_end-return_vec.at(i).flank_start;
      masked_sequence.replace(return_vec.at(i).flank_start, flank_length, string(flank_length,'X'));
    }
    vector<Barcode> masked_res;
    masked_res = big_barcode_search(ModeSeq, seqID, masked_sequence, known_barcodes_X, 
                                    known_barcodes_Y, known_barcodes_Z, 
                                    max_flank_editd, search_pattern, WHITELIST, buindex);
    return_vec.insert(return_vec.end(), masked_res.begin(), masked_res.end());
  }

  return(return_vec);
}


// 把字符串参数(如 "true"、"0"、"f"等)解析成bool值;
// utility function to check true/false input options
bool get_bool_opt_arg(string value){
  transform(value.begin(), value.end(), value.begin(), ::tolower);
  if (value.compare("true")==0 || value.compare("t")==0 || value.compare("1")==0){
    return true;
  } else if (value.compare("false")==0 || value.compare("f")==0 || value.compare("0")==0){
    return false;
  } else {
    cerr << "Unknown argument to boolean option\n";
    print_usage();
    exit(1);
  }
}

// print information about barcodes
void print_stats(string read_id, vector<Barcode> & vec_bc, ostream & out_stream){
  for(int b = 0; b < vec_bc.size(); b++) {
    out_stream << read_id << "\t"
	       << vec_bc.at(b).barcode << "\t"
	       << vec_bc.at(b).flank_editd << "\t"
	       << vec_bc.at(b).editd << "\t"
	       << vec_bc.at(b).umi << "\n";
  }
}

void print_line(string id, string read, string quals, bool is_fastq, ostream & out_stream){
  //output to the new read file
    if(is_fastq)
      out_stream << "@" << id << "\n";
    else
      out_stream << ">" << id << "\n";
    out_stream << read << "\n";
    if(is_fastq) {
      out_stream << "+" << id << "\n";
      out_stream << quals << "\n";
    }
}

//print fastq or fasta lines..
void print_read(string read_id, string read, string qual,
		vector<Barcode> & vec_bc, string prefix,
		bool split, unordered_set<string> & found_barcodes,
		bool trim_barcodes, bool chimeric, bool is_fastq) {
    auto vec_size = vec_bc.size();
    // loop over the barcodes found... usually will just be one.
    // 遍历所有检测到的barcodes, 通常一个read只有一个barcode;
    for (int b = 0; b < vec_size; b++) {
      // format the new read id. Using FLAMES format.
      stringstream ss;
      ss << (b + 1) << "of" << vec_size;
      if (chimeric) {
        ss << "_" << "C";
      }

      string barcode = vec_bc.at(b).barcode;
      // also add the proper FASTQ way: \tCB:Z:barcode\tUB:Z:umi
      string new_read_id =
          barcode + "_" + vec_bc.at(b).umi + "#" + read_id + ss.str() + "\tCB:Z:" + barcode + "\tUB:Z:" + vec_bc[b].umi;

      // work out the start and end base in case multiple barcodes
      // note: read_start+1, see issue #63
      int read_start = vec_bc.at(b).flank_end + 1;
      int read_length = read.length() - read_start;

      for (int f = 0; f < vec_size; f++) {
        int temp_read_length = vec_bc.at(f).flank_start - read_start;
        if (temp_read_length >= 0 && temp_read_length < read_length)
          read_length = temp_read_length;
      }
      string qual_new = ""; // don't trim the quality scores if it's a fasta file

      if (qual != "") {
        if((read_start+read_length)>(qual.length())) {
          std::cerr << "WARNING: sequence and quality lengths diff for read: " << read_id << ". Ignoring read." << endl;
          return;
        }
        qual_new = qual.substr(read_start, read_length);
      }
      string read_new = read.substr(read_start, read_length);

      if (b == 0 && !trim_barcodes) { // override if read shouldn't be cut
        new_read_id = read_id;
        read_new = read;
        qual_new = qual;
        b = vec_size; // force loop to exit after this iteration
      }

      // Skip reads that become empty after trimming
      if (read_new.length() == 0) {
        continue;
      }
      if (split) { 
        string outname = prefix + "_" + barcode + ".";
        if (qual == "")
          outname += "fasta";
        else
          outname += "fastq";
        ofstream outstream;
        if (found_barcodes.insert(barcode).second)
          outstream.open(outname); // override file if this is the first read
                                   // for the barcode
        else
          outstream.open(outname, ofstream::app); // remove file if this is the
                                                  // first read for the barcode
        print_line(new_read_id, read_new, qual_new, is_fastq, outstream);
        outstream.close();
      } else {
        print_line(new_read_id, read_new, qual_new, is_fastq, std::cout);
      }
    }
}

// 分线程负责对每个read进行条形码(barcode)匹配搜索;
// 对输入的每条read(DNA 序列), 在正向和反向两条链上搜索已知barcode, 记录找到的barcode列表,匹配数量,以及是否是嵌合序列;
// separated out from main so that this can be run with threads
void search_read(const std::string& ModeSeq, std::vector<SearchResult>& reads,
                unordered_set<string>& known_barcodes, unordered_set<string>& known_barcodes_Y, 
                unordered_set<string>& known_barcodes_Z, int flank_edit_distance,
                const std::vector<std::pair<std::string, std::string>>& search_pattern,
                std::unordered_map<std::string, std::vector<uint32_t>>& WHITELIST,
                BarcodeUMIindex& buindex) {

  for (int r = 0; r < reads.size(); r++) {

    auto forward_reads = big_barcode_search(ModeSeq, reads[r].read_id, reads[r].line, known_barcodes, known_barcodes_Y, 
                                            known_barcodes_Z, flank_edit_distance,
                                            search_pattern, WHITELIST, buindex);
    // get reverse complement
    reads[r].rev_line = reads[r].line;
    reverse_complement(reads[r].rev_line);

    auto reverse_reads = big_barcode_search(ModeSeq, reads[r].read_id, reads[r].rev_line, known_barcodes, known_barcodes_Y,
                                            known_barcodes_Z, flank_edit_distance,
                                            search_pattern, WHITELIST, buindex);
    reads[r].vec_bc_for = forward_reads;
    reads[r].vec_bc_rev = reverse_reads;
    reads[r].count = forward_reads.size() + reverse_reads.size();
    reads[r].chimeric = forward_reads.size() && reverse_reads.size();

    if (forward_reads.size() == 0) noJG = noJG + 1;

  }
}


// MAIN is here!!
int main(int argc, char **argv) {

  std::ios_base::sync_with_stdio(false);

  std::cerr << R"(
   ____               ____   ___  _      ___ 
  | __ )  _ __ ___   / ___| / _ \| |    |_ _|
  |  _ \ | '__/ _ \ | |    | | | | |     | | 
  | |_) || | | (_) || |___ | |_| | |___  | | 
  |____/ |_|  \___/  \____| \___/|_____||___|
  )" << "\n" << "preBroCOLI version" << VERSION << "\n";

  int expected_cells = 0;       //(d)
  int edit_distance = 8;        //(e) //暂时没用到-20251218;
  int flank_edit_distance = 20; //(f)

  std::string out_stat_filename = "reads_barcodes.txt";
  std::string out_bc_filename = "barcodes_counts.txt";
  std::string out_filename_prefix = "preBroCOLI"; //(n)

  bool split_file_by_barcode = false; //(s)
  bool remove_barcodes = true;        //(r)
  bool print_chimeric = false;        //(c)

  std::string input_whitelist_filename;
  std::string input_barcodeX_filename;
  std::string input_barcodeY_filename;
  std::string input_barcodeZ_filename;
  std::string MODE_str;

  // Set of known barcodes
  std::unordered_set<std::string> known_barcodesX; known_barcodesX.reserve(3000);
  std::unordered_set<std::string> known_barcodesY; known_barcodesY.reserve(3000);
  std::unordered_set<std::string> known_barcodesZ; known_barcodesZ.reserve(3000);
  std::unordered_set<std::string> found_barcodes; found_barcodes.reserve(3000);
  std::unordered_map<std::string, std::vector<uint32_t>> whiteList; whiteList.reserve(3000);
  
  int n_threads = 1;
  int c;
  int params = 1;

  // 将指针位置为argv到argv+argc的所有元素构建到vector中;
  std::vector<char *> myArgs(argv, argv + argc);

  // 从命令行参数中解析选项; myArgs.size()就是argc; myArgs.data()返回的就是一个char**, 像argv一样使用;
  while ((c = getopt(myArgs.size(), myArgs.data(),
                     "q:x:y:z:w:u:i:e:f:n:s:hp:c:")) != EOF) {
    switch (c) {

      case 'q': {
        MODE_str = optarg;
        if (MODE_str == "magicseq") flank_edit_distance=20;
        if (MODE_str == "10x3v3") flank_edit_distance=8;
        if (MODE_str == "visium") flank_edit_distance=8;
        params += 2;
        break;
      }
      case 'x': { 
        input_barcodeX_filename = optarg;
        params += 2;
        break;
      }

      case 'y': { 
        input_barcodeY_filename = optarg;
        params += 2;
        break;
      }

      case 'z': { 
        input_barcodeZ_filename = optarg;
        params += 2;
        break;
      }

      case 'w': {
        input_whitelist_filename = optarg;
        params += 2;
        break; 
      }

      case 'i': {
        // 原始read:[flank_left][barcode][UMI][flank_right][real RNA sequence]
        // 默认是true; (1)把原始的read名称替换成由barcode+UMI组成的新ID;
        // (2)移除匹配到的flanking序列(包括前后两端的); (3)如果一个read中检测到多个barcode, 则将该read拆分成多个片段;
        remove_barcodes = get_bool_opt_arg(optarg);
        std::cerr << "Setting read IDs to be replaced: " << remove_barcodes << "\n";
        params += 2;
        break;
      }

      case 'e': {
        edit_distance = atoi(optarg);
        std::cerr << "Setting max barcode edit distance to " << edit_distance << "\n";
        params += 2;
        break;
      }

      case 'f': {
        flank_edit_distance = atoi(optarg);
        std::cerr << "Setting max flanking sequence edit distance to " << flank_edit_distance << "\n";
        params += 2;
        break;
      }

      case 'h': {
        print_usage();
        exit(0);
      }

      case 'n': {
        out_filename_prefix = optarg;
        std::cerr << "Setting output filename prefix to: " << out_filename_prefix << "\n";
        params += 2;
        break;
      }

      case 's': {
        split_file_by_barcode = get_bool_opt_arg(optarg);
        std::cerr << "Split read output into separate files by barcode: " << split_file_by_barcode << "\n";
        int max_split_bc = 50;
        if (known_barcodesX.size() > max_split_bc) {
          std::cerr << "Too many barcodes to split into separate files: " << known_barcodesX.size() << "> " << max_split_bc << "\n";
          split_file_by_barcode = false;
        }
        params += 2;
        break;
      }

      case 'c': {
        // 如果reds是嵌合体, 后面加上"_C"标记; 默认是false;
        print_chimeric = get_bool_opt_arg(optarg);
        params += 2;
        break;
      }

      case 'p': {
        n_threads = std::atoi(optarg);
        std::cerr << "Setting number of threads to " << n_threads << endl;
        params += 2;
        break;
      }

      case '?':
        std::cerr << "Unknown option.. stopping" << endl;
        print_usage();
        exit(1);
    }
  }

  BarcodeUMIindex BUindex;

  get_seq_information(MODE_str, input_whitelist_filename, input_barcodeX_filename, input_barcodeY_filename, input_barcodeZ_filename, 
                    whiteList, known_barcodesX, known_barcodesY, known_barcodesZ, BUindex);

  std::vector<std::pair<std::string, std::string>> search_pattern;
  get_seq_search_pattern(MODE_str, search_pattern);  

  std::istream *in;
  std::ifstream reads_ifs;

  if (params >= myArgs.size()) {
    std::cerr << "No filename given... getting reads from stdin...\n";
    in = &std::cin;
  } else {
    std::string reads_file = myArgs[params];
    reads_ifs.open(reads_file);
    if (!(reads_ifs.good())) {
      std::cerr << "Unable to open file " << reads_file << "\n";
      print_usage();
      exit(1);
    }
    in = &reads_ifs;
  }

  int bc_count = 0;
  int r_count = 0;
  int multi_bc_count = 0;

  std::ofstream out_stat_file;
  out_stat_filename = out_filename_prefix + "_" + out_stat_filename; // preBroCOLI_reads_barcodes.txt
  out_bc_filename = out_filename_prefix + "_" + out_bc_filename; // preBroCOLI_barcodes_counts.txt;

  if (known_barcodesX.size() > 0) {
    if (file_exists(out_stat_filename)){
      std::cerr << "File " << out_stat_filename << " already exists, overwriting.\n";
    } else {
      std::cerr << "File " << out_stat_filename << ", we can build it.\n";
    }
    out_stat_file.open(out_stat_filename, std::ios::trunc);
    out_stat_file << "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI\n";
  }

  std::cerr << "Searching for barcodes...\n";
  bool is_fastq = true;
  std::unordered_map<std::string, int> barcode_counts;
  std::string read_id_line;

  if (getline(*in, read_id_line)) { 
    if (read_id_line[0] == '>') { 
      is_fastq = false;
    } else if (read_id_line[0] == '@') {
    } else {
      std::cerr << "Unknown read format... exiting\n";
      exit(1);
    }
  }

  std::string line;
  while (getline(*in, line)) { 
    const int buffer_size = 2000;
    // 存储; 根据线程; 每个线程2000条reads;
    std::vector<std::vector<SearchResult>> sr_v(n_threads);
    for (int i = 0; i < n_threads; i++) {
      sr_v[i] = std::vector<SearchResult>(buffer_size);
    }

    std::vector<std::thread> threads(n_threads);

    // 首先遍历 线程×2000条reads; 存储信息然后处理;
    for (int t = 0; t < n_threads; t++) {
      for (int b = 0; b < buffer_size; b++) {

        SearchResult &sr = sr_v[t][b]; // 每个线程的第t个线程, 第b条reads;

        sr.line = line;
        std::string read_id;

        std::istringstream line_stream(read_id_line);
        // 字符串变为字符串流, 首先读取空格前的第一个部分, 也就是read的名称;
	      line_stream >> sr.read_id;
        // 删除第一个字符@;
        sr.read_id.erase(0, 1);

        if (!is_fastq) {
          std::string buffer_string;
          while (getline(*in, buffer_string) && buffer_string[0] != '>')
            sr.line += buffer_string;
          read_id_line = buffer_string;
        } else {
          // fastq格式;
          for (int s = 0; s < 2; s++) {
            getline(*in, sr.qual_scores); //第四行; 质量分数;
          }
          // 下一条read的名称行;
          getline(*in, read_id_line);
        }

        r_count++;
        if ( (r_count < 100000 && r_count % 10000 == 0) || (r_count % 100000 == 0)) {
          std::cerr << r_count / ((double)1000000) << " million reads processed..\n";
        }

        // this is quite ugly, must be a better way to do this..
        if (b == buffer_size - 1 && t == n_threads - 1) {
          // 如果当前循环的位置, 正好是最后一个线程(t), 并且这个线程的最后一个read(b), 就直接break;
          break; 
        } else if (!getline(*in, line)) {
          // 如果读取下一行失败(即文件结束了), 说明这批是最后一批reads;
          sr_v[t].resize(b + 1);
          threads[t] = std::thread(search_read, ref(MODE_str), ref(sr_v[t]),
                                   ref(known_barcodesX), ref(known_barcodesY),
                                   ref(known_barcodesZ), flank_edit_distance,
                                   ref(search_pattern), ref(whiteList),
                                   ref(BUindex) );
          for (int t2 = t + 1; t2 < n_threads; t2++) {
            sr_v[t2].resize(0);
          }
          goto print_result;
        }
      }
      threads[t] = std::thread(search_read, ref(MODE_str), ref(sr_v[t]), ref(known_barcodesX), ref(known_barcodesY),
                                ref(known_barcodesZ), flank_edit_distance, ref(search_pattern), ref(whiteList),
                                ref(BUindex) );
    }

  print_result:
    for (int t = 0; t < sr_v.size(); t++) {
      if (sr_v[t].size() > 0) threads[t].join();
        
      for (int r = 0; r < sr_v[t].size(); r++) {

        for (int b = 0; b < sr_v[t][r].vec_bc_for.size(); b++) {
          barcode_counts[sr_v[t][r].vec_bc_for.at(b).barcode]++;
        }
    
        for (int b = 0; b < sr_v[t][r].vec_bc_rev.size(); b++) {
          barcode_counts[sr_v[t][r].vec_bc_rev.at(b).barcode]++;
        }
        // 含有至少一个有效barcode的read数;
        if (sr_v[t][r].count > 0) {
          bc_count++;
        }
        // 嵌合(chimeric)的read数;
        if (sr_v[t][r].chimeric) {
          multi_bc_count++;
        }

        if (known_barcodesX.size() != 0) { 
          // if we are just looking for all possible barcodes don't output reads etc.
          print_stats(sr_v[t][r].read_id, sr_v[t][r].vec_bc_for, out_stat_file);
          print_stats(sr_v[t][r].read_id, sr_v[t][r].vec_bc_rev, out_stat_file);

          print_read(
            sr_v[t][r].read_id + "_+",
            sr_v[t][r].line,
            sr_v[t][r].qual_scores,
            sr_v[t][r].vec_bc_for,
            out_filename_prefix,
            split_file_by_barcode,
            found_barcodes,
            remove_barcodes,
            print_chimeric && sr_v[t][r].chimeric, // include chimeric information if requested
            is_fastq
          );

          // case we just want to print read once if multiple bc found.
          if (remove_barcodes || sr_v[t][r].vec_bc_for.size() == 0) {
            // for a chimeric read, we first need to reverse the quality scores
            reverse(sr_v[t][r].qual_scores.begin(), sr_v[t][r].qual_scores.end());

            print_read(
              sr_v[t][r].read_id + "_-",
              sr_v[t][r].rev_line,
              sr_v[t][r].qual_scores,
              sr_v[t][r].vec_bc_rev,
              out_filename_prefix,
              split_file_by_barcode,
              found_barcodes,
              remove_barcodes,
              print_chimeric && sr_v[t][r].chimeric,
              is_fastq
            );
          }
        }
      }
    }

  }
  reads_ifs.close();

  cerr << "Number of reads processed: " << r_count << "\n";
  cerr << "Number of reads where a barcode was found: " << bc_count << "\n";
  cerr << "Number of reads where more than one barcode was found: " << multi_bc_count << "\n";
  cerr << "All done!" << endl;
  cerr << "If you like BroCOLI, please cite us! " << endl;
  cerr << "First total: " << total_First/1000000 << " s\n";
  cerr << "BarcodeX total: " << total_BarcodeX/1000000 << " s\n";
  cerr << "BarcodeY total: " << total_BarcodeY/1000000 << " s\n";
  cerr << "BarcodeZ total: " << total_BarcodeZ/1000000 << " s\n";
  cerr << "Final total: " << total_Final/1000000 << " s\n";
  cerr << "直接找不到结构的: " << noJG << " \n";
  cerr << "长度不够的:" << noLen << "\n";
  cerr << "polyT这里不能找到index:" << polyTOK << "\n";
  cerr << "不能用编辑距离比对上找到结构:" << noAlign << "\n";
  cerr << "不能找到X:" << noX << "\n";
  cerr << "不能找到Y:" << noY << "\n";
  cerr << "不能找到Z:" << noZ << "\n";
  cerr << "不能找到Final:" << noFinal << "\n";
  cerr << "Ambiguous:" << Ambiguous << "\n";
  cerr << "总体小于>6:" << total6 << "\n";
  cerr << "------------------------------------" << "\n";
  cerr << "用的初始的UMI:" << chushi << "\n";
  cerr << "用的二分查找的UMI:" << erfen << "\n";
  cerr << "用的bit位的UMI:" << bitwei << "\n";
  cerr << "guagua:" << guagua << "\n";
  

  if (known_barcodesX.size() > 0) {
    out_stat_file.close();
    return (0);
  }

  if (barcode_counts.size() == 0)
    return (0);

  typedef std::pair<std::string, int> pair;
  std::vector<pair> bc_vec;

  copy(barcode_counts.begin(), barcode_counts.end(),
       back_inserter<vector<pair>>(bc_vec));

  std::sort(bc_vec.begin(), bc_vec.end(), [](const pair &l, const pair &r) {
    if (l.second != r.second)
      return l.second > r.second;
    return l.first < r.first;
  });

  std::vector<int> hist(bc_vec[0].second);
  std::ofstream out_bc_file;

  out_bc_file.open(out_bc_filename);

  for (auto const &bc_pair : bc_vec) {
    out_bc_file << bc_pair.first << "\t" << bc_pair.second << "\n";
    hist[bc_pair.second - 1]++;
  }

  out_bc_file.close();

  std::cout << "Reads\tBarcodes" << "\n";
  for (int i = hist.size() - 1; i >= 0; i--) std::cout << i + 1 << "\t" << hist[i] << "\n";

}


#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <list>
#include <numeric>
#include <cmath>
#include <limits>
#include <Eigen/Dense>
#include <getopt.h>
#include <sys/stat.h>   
#include <sys/types.h>  
#include <dirent.h>
#include <thread>
#include <mutex>
#include <iomanip>
#include <condition_variable>
#include <queue>
#include <functional>
#include <atomic>
#include "thread_pool_executor.hpp"
#include <chrono>
#include <edlib.h>

std::mutex updatedGtfMutex;
std::mutex traceMutex;
std::mutex bigMutex;

struct Split_Result{
    std::string read_name;
    std::vector<std::string> tokens;
    int read_length = 0;
    int token_count = 0;

    std::string read_barcode;
    std::string read_umi;
};


inline void get_string_split_sc_fast(const std::string& s, char delimiter, 
                                     const std::string& umi_tag_fourth, const std::string& barcode_tag_fourth,
                                     Split_Result& r) {

    r.read_name.clear(); r.read_umi.clear(); r.read_barcode.clear();
    r.read_length = 0; r.token_count = 0;
    
    if (r.tokens.size() < 6) r.tokens.resize(6);

    const char* base = s.data();
    const size_t n = s.size();
    size_t start = 0;
    int field_idx = 0; 
    bool is_flexiplex = false;

    const size_t umi_tag_len = umi_tag_fourth.size();
    const size_t barcode_tag_len = barcode_tag_fourth.size();

    for (size_t i = 0; i <= n; ++i) {

        if (i == n || base[i] == delimiter) {
            const char* p = base + start;
            const size_t len = i - start;

            // --- 第 1 列 (field 0): read_name 及 Flexiplex 预检查 ---
            if (field_idx == 0) {
                r.read_name.assign(p, len);

                // Flexiplex 逻辑检查: 包含 "-1_"
                if (r.read_name.find("#") != std::string::npos) {
                    size_t first_underscore = r.read_name.find('_');
                    if (first_underscore != std::string::npos) {
                        std::string barcode_part = r.read_name.substr(0, first_underscore);

                        if (barcode_part.size() >= 2 && 
                            barcode_part.compare(barcode_part.size() - 2, 2, "-1") == 0) {
                             barcode_part.resize(barcode_part.size() - 2);
                        }
                        r.read_barcode = barcode_part;

                        size_t hash_pos = r.read_name.find('#', first_underscore + 1);
                        if (hash_pos != std::string::npos) {
                            r.read_umi = r.read_name.substr(first_underscore + 1, hash_pos - (first_underscore + 1));
                        }
                    }
                }

                // 如果成功提取了 Flexiplex 信息，更新 flag 并重组 read_name
                if (!r.read_umi.empty() && !r.read_barcode.empty()) {
                    is_flexiplex = true;
                    r.read_name = r.read_barcode + "-" + r.read_umi;
                }
            }
            // --- 第 3 到 7 列 (field 2-6): 提取 tokens ---
            // 对应原逻辑: 2 < number < 8 (即 number 为 3,4,5,6,7)
            else if (field_idx >= 1 && field_idx <= 5) {
                if (r.token_count < (int)r.tokens.size()) {
                    r.tokens[r.token_count].assign(p, len);
                    ++r.token_count;
                }
            }
            // --- 第 11 列 (field 10): read_length ---
            else if (field_idx == 10) {
                r.read_length = (int)len; // 注意：原代码用的 token.size() 即字符串长度
            }
            // --- 其他列 (Sicelore 标签检查) ---
            else {
                // 如果是 Flexiplex 模式，原逻辑在 number > 11 后直接返回(忽略后续)
                if (is_flexiplex) {
                    if (field_idx > 10) { 
                        // 相当于原代码的 return read_result;
                        // 为了保持单循环结构，我们直接 break 循环
                        break; 
                    }
                } 
                else {
                    // 非 Flexiplex 模式，检查剩余列的 Tag
                    if (len > 11) {
                        // 检查 UMI Tag 前缀
                        if (len >= umi_tag_len && std::strncmp(p, umi_tag_fourth.c_str(), umi_tag_len) == 0) {
                            // r.read_umi.assign(p, len);
                            r.read_umi.assign(p + umi_tag_len, len - umi_tag_len);
                        }
                        // 检查 Barcode Tag 前缀
                        else if (len >= barcode_tag_len && std::strncmp(p, barcode_tag_fourth.c_str(), barcode_tag_len) == 0) {
                            // r.read_barcode.assign(p, len);
                            r.read_barcode.assign(p + barcode_tag_len, len - barcode_tag_len);
                        }
                    }
                }
            }
            // 准备下一次迭代
            start = i + 1;
            ++field_idx;
        }
    }
    // 循环结束后，如果是 Sicelore 模式且找到了 UMI/Barcode，重组 read_name
    if (!is_flexiplex && !r.read_umi.empty() && !r.read_barcode.empty()) {
        r.read_name = r.read_barcode + "-" + r.read_umi;
    }
}


int IntervalIntersection(std::vector<std::array<int,2>> firstList, std::vector<std::array<int,2>> secondList) {
    std::vector<std::array<int,2>> res; //定义返回的结果;
    int result = 0;
    int i = 0;
    int j = 0;

    while(i < firstList.size() && j < secondList.size()){
        int start = std::max(firstList[i][0], secondList[j][0]);
        int end = std::min(firstList[i][1], secondList[j][1]);
        if(start <= end) {
            res.push_back({start, end});
            result = result + (end-start+1);       
        }
        if(firstList[i][1] < secondList[j][1]){
            i++;
        }else{
            j++;
        }
    }
    return result;
}


int IntervalMerge(std::vector<std::array<int,2>> firstList, std::vector<std::array<int,2>> secondList) {

    std::vector<std::array<int,2>> intervals;
    intervals.resize(firstList.size()+secondList.size());
    merge(firstList.begin(), firstList.end(), secondList.begin(), secondList.end(), intervals.begin());
    
    sort(intervals.begin(), intervals.end()); 
    std::vector<std::array<int,2>> merged; 
    int result = 0;
    
    for (int i=0; i<intervals.size();) {  
        int left = intervals[i][0];
        int right = intervals[i][1];
        while (i+1 < intervals.size() && right >= intervals[i+1][0]){
            right = std::max(right, intervals[i+1][1]);
            i++;
        }
        merged.push_back({left, right});
        result = result + (right-left+1);
        i++;
    }
    return result;
}


int findNonOverlappingIntervals(const std::vector<std::array<int, 2>>& intervals1, const std::vector<std::array<int, 2>>& intervals2) {
    std::vector<std::array<int, 2>> result;
    int finally_number = 0;

    for (const auto& interval1 : intervals1) {
        bool hasOverlap = false;
        for (const auto& interval2 : intervals2) {
            if (!(interval1[1] < interval2[0] || interval1[0] > interval2[1])) {

                hasOverlap = true;
                break;
            }
        }
        if (!hasOverlap) {
            result.push_back(interval1);
            finally_number = finally_number + (interval1[1]-interval1[0]+1);
        }
    }
    for (const auto& interval2 : intervals2) {
        bool hasOverlap = false;

        for (const auto& interval1 : intervals1) {
            if (!(interval2[1] < interval1[0] || interval2[0] > interval1[1])) {

                hasOverlap = true;
                break;
            }
        }
        if (!hasOverlap) {
            result.push_back(interval2);
            finally_number = finally_number + (interval2[1]-interval2[0]+1);
        }
    }
    return finally_number;
}

int IntervalMinDistance(const std::vector<std::array<int,2>>& firstList,
                        std::array<int,2> secondList) {
    int left_min  = INT_MAX;
    int right_min = INT_MAX;
    for (const auto& a : firstList) {
        left_min  = std::min(left_min,  std::abs(secondList[0] - a[0]));
        right_min = std::min(right_min, std::abs(secondList[1] - a[1]));
    }
    return std::min(left_min, right_min);
}

std::vector<std::array<int,2>> mergeIntervals(std::vector<std::array<int,2>>& intervals) {
    if (intervals.empty()) return {};
    std::sort(intervals.begin(), intervals.end(),
              [](const std::array<int,2>& a, const std::array<int,2>& b) { return a[0] < b[0]; });
    std::vector<std::array<int,2>> merged;
    merged.push_back(intervals[0]);

    for (size_t i = 1; i < intervals.size(); ++i) {
        if (intervals[i][0] <= merged.back()[1]) {
            merged.back()[1] = std::max(merged.back()[1], intervals[i][1]);
        } else {
            merged.push_back(intervals[i]);
        }
    }
    return merged;
}


int get_read_match_gene_length(std::vector<std::array<int,2>> CIGAR_vector){

    int Mlength = 0;
    for (auto each_array : CIGAR_vector) {
        Mlength = Mlength + (each_array[1] - each_array[0]);
    }
    return Mlength;
}


struct Read_intervals_and_Mlength{
    std::vector<std::array<int,2>> ReadIntervals;
    int ReadMatchLength;
};

Read_intervals_and_Mlength get_read_intervals_fast(const std::string& CIGARvalue, const std::string& initpos){

    Read_intervals_and_Mlength out;
    out.ReadMatchLength = 0;
    out.ReadIntervals.reserve(CIGARvalue.size() / 2 + 1);

    int pos = 0; 
    for (char ch : initpos) pos = pos * 10 + (ch - '0');

    int number = 0;
    for (char c : CIGARvalue) {
        if (c >= '0' && c <= '9') {
            number = number * 10 + (c - '0');
        } else {
            switch (c) {
            case 'M': case '=': case 'X':
                out.ReadIntervals.emplace_back(std::array<int,2> {pos, pos + number});
                pos += number;
                out.ReadMatchLength += number;
                break;
            case 'D': case 'N':
                pos += number;
                break;
            case 'I': case 'S':
                out.ReadMatchLength += number;
                break;
            case 'H': case 'P': break;
            }
            number = 0;
        }
    }
    return out;
}



struct Reads_Clusters {

    std::streampos lastPos;
    std::streampos newPos;
    std::map<std::string, std::vector<std::array<int,2>>> Mymap;
    std::unordered_map<std::string, std::string> ReadsBarcode;
    std::unordered_map<std::string, std::string> ReadsUMI;
    std::map<std::string, int> Mylen;
    std::map<std::string, int> MyFlag;
    std::string SetRef_name;
    std::array<int,2> ClusterCoverage;
    int mapqless1 = 0;
    int sameumi = 0;
    int mappingdiff = 0;
    int surveyNum = 0;    
}; 


Reads_Clusters get_each_cluster_reads_sc(std::ifstream& samfile, std::streampos CurrentPos, std::streampos EndPos, const int& MAPQ,
                                        const std::string& umi_tag_third, const std::string& barcode_tag_third) {

    Reads_Clusters NewCluster = {};
    std::map<std::string, std::vector<std::array<int,2>>> read_informs;
    std::map<std::string, int> read_len;
    std::map<std::string, int> read_flag;

    std::unordered_map<std::string, int> BarcodeUMI2reads;
    std::unordered_map<std::string, std::string> read_of_Barcode;
    std::unordered_map<std::string, std::string> read_of_UMI;

    std::string line;
    std::string now_gene;
    std::string last_chr;

    samfile.seekg(CurrentPos, std::ios::beg);
    std::streampos earlyPos = CurrentPos;
    int early_begin_pos {0};
    int early_end_pos {0};
    int read_mapq {0};
    int Flag {0};

    Split_Result This_Line;
    This_Line.tokens.resize(6);
    Read_intervals_and_Mlength CIGAR_interval {};

	while (getline(samfile, line))
	{
        earlyPos = CurrentPos; 
        CurrentPos = samfile.tellg();

        if (line[0] != '@') {
            NewCluster.surveyNum = NewCluster.surveyNum + 1;
            
            get_string_split_sc_fast(line, '\t', umi_tag_third, barcode_tag_third, This_Line); 

            if (This_Line.read_barcode.size() == 0 || This_Line.read_umi.size() == 0) continue;

            int read_mapq = std::stoi(This_Line.tokens[3]);
            int Flag = std::stoi(This_Line.tokens[0]);

            if (read_mapq >= MAPQ and (Flag & 0x900)==0) {

                if ( (BarcodeUMI2reads.find(This_Line.read_name) == BarcodeUMI2reads.end()) || (BarcodeUMI2reads[This_Line.read_name] < This_Line.read_length) ) {
                    
                    BarcodeUMI2reads[This_Line.read_name] = This_Line.read_length;         
                    CIGAR_interval = get_read_intervals_fast(This_Line.tokens[4], This_Line.tokens[2]);
                    if (This_Line.read_length == CIGAR_interval.ReadMatchLength) {

                        now_gene = This_Line.tokens[1];
                        int read_in_gene_begin_pos = std::stoi(This_Line.tokens[2]);
                        std::array<int,2> lastEle = CIGAR_interval.ReadIntervals.back();
                        int read_in_gene_end_pos = lastEle[1] - 1;

                        if (early_begin_pos == 0) {
                            early_begin_pos = read_in_gene_begin_pos;
                            early_end_pos = read_in_gene_end_pos;
                            last_chr = now_gene;
                        } else {
                            if (now_gene != last_chr){
                                NewCluster.lastPos = earlyPos;
                                NewCluster.newPos = CurrentPos;
                                NewCluster.Mymap = read_informs;
                                NewCluster.Mylen = read_len;
                                NewCluster.MyFlag = read_flag;
                                NewCluster.ClusterCoverage[0] = early_begin_pos;
                                NewCluster.ClusterCoverage[1] = early_end_pos;
                                NewCluster.SetRef_name = last_chr;
                                NewCluster.ReadsBarcode = read_of_Barcode;
                                NewCluster.ReadsUMI = read_of_UMI;
                                NewCluster.surveyNum = NewCluster.surveyNum - 1;
                                break;
                            } else{
                                if ((read_in_gene_begin_pos - early_end_pos > 0) || (early_begin_pos - read_in_gene_end_pos > 0)) {
                                    NewCluster.lastPos = earlyPos;
                                    NewCluster.newPos = CurrentPos;
                                    NewCluster.Mymap = read_informs;
                                    NewCluster.Mylen = read_len;
                                    NewCluster.MyFlag = read_flag;
                                    NewCluster.ClusterCoverage[0] = early_begin_pos;
                                    NewCluster.ClusterCoverage[1] = early_end_pos;
                                    NewCluster.SetRef_name = last_chr;
                                    NewCluster.ReadsBarcode = read_of_Barcode;
                                    NewCluster.ReadsUMI = read_of_UMI;
                                    NewCluster.surveyNum = NewCluster.surveyNum - 1;
                                    break;
                                } else {
                                    early_begin_pos = (read_in_gene_begin_pos>early_begin_pos)?early_begin_pos:read_in_gene_begin_pos;
                                    early_end_pos = (read_in_gene_end_pos>early_end_pos)?read_in_gene_end_pos:early_end_pos;
                                    last_chr = now_gene;
                                }
                            }
                        }
                        read_informs[This_Line.read_name] = CIGAR_interval.ReadIntervals;
                        read_of_Barcode[This_Line.read_name] = This_Line.read_barcode;
                        read_of_UMI[This_Line.read_name] = This_Line.read_umi;
                        read_len[This_Line.read_name] = This_Line.read_length;

                        if (Flag & 16) {
                            read_flag[This_Line.read_name] = 0; // -
                        } else {
                            read_flag[This_Line.read_name] = 1; // +
                        }

                    } else {
                        NewCluster.mappingdiff = NewCluster.mappingdiff + 1;
                    }                                     
                } else {
                    //UMI相同的reads个数;
                    NewCluster.sameumi = NewCluster.sameumi + 1;
                }
            } else {
                //mapping quality小于1;
                NewCluster.mapqless1 = NewCluster.mapqless1 + 1;
            }
        }
	}
    if (CurrentPos >= EndPos) {
        NewCluster.lastPos = earlyPos;
        NewCluster.newPos = CurrentPos;
        NewCluster.Mymap = read_informs;
        NewCluster.Mylen = read_len;
        NewCluster.MyFlag = read_flag;
        NewCluster.ClusterCoverage[0] = early_begin_pos;
        NewCluster.ClusterCoverage[1] = early_end_pos;
        NewCluster.SetRef_name = last_chr;
        NewCluster.ReadsBarcode = read_of_Barcode;
        NewCluster.ReadsUMI = read_of_UMI;
    }
    return NewCluster;
}


struct SpliceJs {

    std::unordered_map<std::string, std::vector<std::array<int,2>>> reads_SJs;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> reads_SJs_left;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> reads_SJs_right;
    std::unordered_map<std::string, std::array<int,2>> reads_begin_end;
    std::unordered_map<std::string, std::array<int,2>> reads_single_exon;

}; 


SpliceJs get_reads_allSJs(std::map<std::string, std::vector<std::array<int,2>>>& RCs, int SJD) {

    SpliceJs NewSJs = {};
    std::unordered_map<std::string, std::vector<std::array<int,2>>> AllRead_SJs;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> AllRead_SJs_left;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> AllRead_SJs_right;
    std::unordered_map<std::string, std::array<int,2>> AllRead_begin_end;
    std::unordered_map<std::string, std::array<int,2>> AllRead_SingleExon;

    for (const auto &pair : RCs) {

        std::string Rname = pair.first; 
        std::vector<std::array<int,2>> RMs = pair.second;
        std::vector<std::array<int,2>> ReadSJs;
        std::vector<std::array<int,2>> ReadSJs_left;
        std::vector<std::array<int,2>> ReadSJs_right;
        std::array<int,2> Read_SingleExon;
        std::array<int,2> Read_begin_end;

        for (auto it = RMs.begin(); it!= RMs.end()-1; it++) {
            
            if (it == RMs.begin()){
                Read_begin_end[0] = (*(it))[0];
            }
            if (it== RMs.end()-2){
                Read_begin_end[1] = (*(it+1))[1];
            }
            int later = (*(it+1))[0];
            int former = (*(it))[1];
            if (later - former > SJD) {

                std::array<int,2> eachSJ;
                eachSJ[0] = former;
                eachSJ[1] = later - 1;
                ReadSJs.push_back(eachSJ);

                std::array<int,2> eachSJ_left;
                eachSJ_left[0] = (*(it))[0];
                eachSJ_left[1] = former - 1;
                ReadSJs_left.push_back(eachSJ_left);

                std::array<int,2> eachSJ_right;
                eachSJ_right[0] = later;
                eachSJ_right[1] = (*(it+1))[1];
                ReadSJs_right.push_back(eachSJ_right);
            }
        }

        if (ReadSJs.size() == 0) {
            Read_SingleExon[0] = (*RMs.begin())[0];
            Read_SingleExon[1] = (*(RMs.end()-1))[1] - 1;
            AllRead_SingleExon[Rname] = Read_SingleExon;
        } else {

            AllRead_SJs[Rname] = ReadSJs;
            AllRead_SJs_left[Rname] = ReadSJs_left;
            AllRead_SJs_right[Rname] = ReadSJs_right;
            AllRead_begin_end[Rname] = Read_begin_end;
        }
    }

    NewSJs.reads_SJs = AllRead_SJs;
    NewSJs.reads_SJs_left = AllRead_SJs_left;
    NewSJs.reads_SJs_right = AllRead_SJs_right;
    NewSJs.reads_begin_end = AllRead_begin_end;
    NewSJs.reads_single_exon = AllRead_SingleExon;
    return NewSJs;
}



struct GTF
{
    std::map<std::string, std::map<std::string, std::vector<std::array<int,2>>>> GTF_transcript;
    std::map<std::string, std::map<std::string, std::array<int,2>>> GTF_gene; 
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> GTF_transcript_strand;
    std::map<std::string, std::map<std::string, std::string>> GTF_gene_strand;
    std::map<std::string, std::map<std::string, std::vector<std::string>>> GTF_gene2transcript;
};

struct unGTF
{
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::array<int,2>>>> GTF_transcript;
    std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> GTF_gene; 
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> GTF_transcript_strand;
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> GTF_gene_strand;
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> GTF_gene2transcript;
};

unGTF get_gtf_annotation(std::string& GTFFile_name) {
    GTF GTFAll_Info;

    if (!GTFFile_name.empty()) {

        std::cerr << "***** Now open the gtf file: " << GTFFile_name << "! *****" << std::endl;        

        std::map<std::string, std::vector<std::array<int,2>>> ChrEach;
        std::unordered_map<std::string, std::string> ChrTranscriptStrand;
        std::map<std::string, std::array<int,2>> GeneEach;      
        std::vector<std::array<int,2>> GTFAnno_SJs;  
        std::array<int,2> annoExon;
        std::ifstream GTFFile; 

        GTFFile.open(GTFFile_name);
        if (!GTFFile.is_open())	{
            std::cerr << "We can't open GTF file !" << GTFFile_name << std::endl;
            exit(EXIT_FAILURE);
        }        
        std::string line;
        std::string now_chr_name = "Han";
        std::string early_chr_name = "Han";
        std::string early_gene_transcript_name = "WHH";
        std::string now_gene_transcript_name = "WHH";

        std::string early_Exonstrand = "Wei";
        std::string now_Exonstrand = "Wei";

        while (std::getline(GTFFile, line)){
          
            if (line[0] != '#'){

                std::istringstream lineiss(line);

                std::vector<std::string> tokens;

                std::string token;

                while (std::getline(lineiss, token, '\t')) {
                    tokens.push_back(token);
                } 
                
                if (tokens[2] == "exon"){
                    early_chr_name = now_chr_name;
                    now_chr_name = tokens[0]; 
                    early_Exonstrand = now_Exonstrand;
                    now_Exonstrand = tokens[6]; 

                    std::string AllEndT = tokens.back();

                    std::istringstream EndTT(AllEndT);
                    std::string Endtoken;
                    std::vector<std::string> Endtokens;

                    while (std::getline(EndTT, Endtoken, '"')){
                        Endtokens.push_back(Endtoken);
                    }
                    early_gene_transcript_name = now_gene_transcript_name;
                    now_gene_transcript_name = Endtokens[1]+'|'+Endtokens[3];

                    if ((GTFAnno_SJs.size() != 0) && (early_gene_transcript_name != now_gene_transcript_name)){
                        if ((GTFAnno_SJs.size()>1) && (GTFAnno_SJs[0][0] > GTFAnno_SJs[1][0])){ 
                            std::reverse(GTFAnno_SJs.begin(),GTFAnno_SJs.end());
                        } 
                        ChrEach[early_gene_transcript_name] = GTFAnno_SJs;
                        ChrTranscriptStrand[early_gene_transcript_name] = early_Exonstrand;
                        GTFAnno_SJs.clear();                        
                    }

                    if ((now_chr_name != early_chr_name) && (ChrEach.size() != 0)){
                        GTFAll_Info.GTF_transcript[early_chr_name] = ChrEach;
                        GTFAll_Info.GTF_transcript_strand[early_chr_name] = ChrTranscriptStrand;
                        ChrEach.clear();
                        ChrTranscriptStrand.clear();                     
                    }
                    annoExon[0] = std::stoi(tokens[3]);
                    annoExon[1] = std::stoi(tokens[4]);
                    GTFAnno_SJs.push_back(annoExon);
                }
            }
        }

        if (GTFFile.eof()) {

            if ((GTFAnno_SJs.size()>1) && (GTFAnno_SJs[0][0] > GTFAnno_SJs[1][0])) {
                std::reverse(GTFAnno_SJs.begin(), GTFAnno_SJs.end());
            }            
            ChrEach[now_gene_transcript_name] = GTFAnno_SJs;
            GTFAll_Info.GTF_transcript[now_chr_name] = ChrEach;
            ChrTranscriptStrand[now_gene_transcript_name] = now_Exonstrand;
            GTFAll_Info.GTF_transcript_strand[now_chr_name] = ChrTranscriptStrand;
            std::cerr << "***** The GTF file has been read to the end ! *****" << std::endl;
        } 
        else if (GTFFile.fail()) {
            std::cerr << "File FALSE !" << std::endl;
        }
        else {std::cerr << "unkown reason" << std::endl;}

        GTFFile.close();

        std::vector<int> All_gene_begin;
        std::vector<int> All_gene_end;
        std::vector<std::string> tx_name;
        size_t pos;
        std::string now_gene_name;
        std::string ago_gene_name;
        std::string now_chr;
        std::string ago_chr;
        std::string Strand;
     
        for (const auto& eachChr:GTFAll_Info.GTF_transcript){

            now_chr = eachChr.first;
            GTFAll_Info.GTF_gene[now_chr] = {};
            GTFAll_Info.GTF_gene_strand[now_chr] = {};
            GTFAll_Info.GTF_gene2transcript[now_chr] = {};

            if (All_gene_begin.size() != 0 && now_chr != ago_chr) {
                auto min_it = std::min_element(All_gene_begin.begin(), All_gene_begin.end());
                auto max_it = std::max_element(All_gene_end.begin(), All_gene_end.end());
                GTFAll_Info.GTF_gene[ago_chr][ago_gene_name][0] = *min_it;
                GTFAll_Info.GTF_gene[ago_chr][ago_gene_name][1] = *max_it;
                GTFAll_Info.GTF_gene2transcript[ago_chr][ago_gene_name] = tx_name;
                All_gene_begin.clear();
                All_gene_end.clear();
                tx_name.clear();
                ago_gene_name = "";
            }

            for (const auto& eachGene:eachChr.second) {
                size_t pos = eachGene.first.find('|');
                now_gene_name = eachGene.first.substr(0, pos);
                Strand = GTFAll_Info.GTF_transcript_strand[now_chr][eachGene.first];

                if (ago_gene_name == now_gene_name) {
                    All_gene_begin.push_back(eachGene.second[0][0]);
                    All_gene_end.push_back(eachGene.second[eachGene.second.size()-1][1]);
                    tx_name.push_back(eachGene.first);
                } else {
                    if (ago_gene_name.size() != 0) {
                        auto min_it = std::min_element(All_gene_begin.begin(), All_gene_begin.end());
                        auto max_it = std::max_element(All_gene_end.begin(), All_gene_end.end());
                        GTFAll_Info.GTF_gene[eachChr.first][ago_gene_name][0] = *min_it;
                        GTFAll_Info.GTF_gene[eachChr.first][ago_gene_name][1] = *max_it;
                        GTFAll_Info.GTF_gene2transcript[eachChr.first][ago_gene_name] = tx_name;
                    }
                    All_gene_begin.clear();
                    All_gene_end.clear();
                    tx_name.clear();
                    All_gene_begin.push_back(eachGene.second[0][0]);
                    All_gene_end.push_back(eachGene.second[eachGene.second.size()-1][1]);
                    tx_name.push_back(eachGene.first);
                    GTFAll_Info.GTF_gene_strand[now_chr][now_gene_name] = Strand;
                }
                ago_gene_name = now_gene_name;
            }
            ago_chr = now_chr;
        }

        if (!All_gene_begin.empty() && !ago_gene_name.empty()) {
            auto min_it = std::min_element(All_gene_begin.begin(), All_gene_begin.end());
            auto max_it = std::max_element(All_gene_end.begin(), All_gene_end.end());
            GTFAll_Info.GTF_gene[ago_chr][ago_gene_name][0] = *min_it;
            GTFAll_Info.GTF_gene[ago_chr][ago_gene_name][1] = *max_it;
            GTFAll_Info.GTF_gene2transcript[ago_chr][ago_gene_name] = tx_name;
        }

    } else {
        std::cerr << "***** No GTF files are put into the program! *****" << std::endl;         
    }

    // 转为unordered_map;
    unGTF unGTFAll_Info;
    unGTFAll_Info.GTF_transcript.reserve(GTFAll_Info.GTF_transcript.size());
    for (const auto& outer_pair : GTFAll_Info.GTF_transcript) {
        const std::string& key1 = outer_pair.first;
        const auto& inner_map = outer_pair.second;
        for (const auto& inner_pair : inner_map) {
            const std::string& key2 = inner_pair.first;
            const auto& value = inner_pair.second;
            unGTFAll_Info.GTF_transcript[key1][key2] = value;
        }
    }

    unGTFAll_Info.GTF_gene.reserve(GTFAll_Info.GTF_gene.size());
    for (const auto& outer_pair : GTFAll_Info.GTF_gene) {
        const std::string& key1 = outer_pair.first;
        const auto& inner_map = outer_pair.second;
        for (const auto& inner_pair : inner_map) {
            const std::string& key2 = inner_pair.first;
            const auto& value = inner_pair.second;
            unGTFAll_Info.GTF_gene[key1][key2] = value;
        }
    }    

    unGTFAll_Info.GTF_gene_strand.reserve(GTFAll_Info.GTF_gene_strand.size());
    for (const auto& outer_pair : GTFAll_Info.GTF_gene_strand) {
        const std::string& key1 = outer_pair.first;
        const auto& inner_map = outer_pair.second;
        for (const auto& inner_pair : inner_map) {
            const std::string& key2 = inner_pair.first;
            const auto& value = inner_pair.second;
            unGTFAll_Info.GTF_gene_strand[key1][key2] = value;
        }
    }    

    unGTFAll_Info.GTF_gene2transcript.reserve(GTFAll_Info.GTF_gene2transcript.size());
    for (const auto& outer_pair : GTFAll_Info.GTF_gene2transcript) {
        const std::string& key1 = outer_pair.first;
        const auto& inner_map = outer_pair.second;
        for (const auto& inner_pair : inner_map) {
            const std::string& key2 = inner_pair.first;
            const auto& value = inner_pair.second;
            unGTFAll_Info.GTF_gene2transcript[key1][key2] = value;
        }
    }      

    unGTFAll_Info.GTF_transcript_strand = std::move(GTFAll_Info.GTF_transcript_strand);
    return unGTFAll_Info;
}



struct GTFsj {
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::array<int,2>>>> mSJs;
    std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> mSJsBE;
    std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> SE;
};


GTFsj get_SJs_SE(std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::array<int,2>>>>& Known_Exon){

    GTFsj Known_SJ_SE;

    if (Known_Exon.size() != 0){

        std::unordered_map<std::string, std::vector<std::array<int,2>>> Every_transcript_All;
        std::vector<std::array<int,2>> Every_transcript;
        std::string transaction_name;
        std::string chr_name;
        
        std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::array<int,2>>>> chr_Every_SJs_all;
        std::unordered_map<std::string, std::vector<std::array<int,2>>> Every_SJs_All;
        std::vector<std::array<int,2>> Every_SJs;
        std::array<int,2> SSJs;

        std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> chr_Every_SE_all;
        std::unordered_map<std::string, std::array<int,2>> Every_SE_All;
        std::array<int,2> Every_SE;

        std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> chr_Every_SJ_begin_end_all;
        std::unordered_map<std::string, std::array<int,2>> Every_SJ_begin_end_all;
        std::array<int,2> Every_SJ_begin_end;

        for (const auto& pair : Known_Exon) {

            chr_name = pair.first;
            Every_transcript_All = pair.second;

            Every_SJs_All.clear();
            Every_SE_All.clear();
            Every_SJ_begin_end_all.clear();

            for (auto it = Every_transcript_All.begin(); it != Every_transcript_All.end(); ++it) {

                transaction_name = it->first;
                Every_transcript = it->second;

                if (Every_transcript.size() == 1){

                    Every_SE[0] = (Every_transcript.front())[0];
                    Every_SE[1] = (Every_transcript.front())[1];
                    Every_SE_All[transaction_name] = Every_SE;
                } else {

                    Every_SJs.clear();
                    for (auto exonit = Every_transcript.begin(); exonit != Every_transcript.end() - 1; ++exonit) {
  
                        if (exonit == Every_transcript.begin()) {
                            Every_SJ_begin_end[0] = (*exonit)[0];
                        }
                        if (exonit == Every_transcript.end() - 2) {
                            Every_SJ_begin_end[1] = ((*(exonit+1))[1]);
                        }
                        SSJs[0] = ((*exonit)[1])+1;
                        SSJs[1] = ((*(exonit+1))[0])-1;
                        Every_SJs.push_back(SSJs);
                    }
                    Every_SJs_All[transaction_name] = Every_SJs;
                    Every_SJ_begin_end_all[transaction_name] = Every_SJ_begin_end;
                }
            }

            if (Every_SJs_All.size() != 0) {
                chr_Every_SJs_all[chr_name] = Every_SJs_All;
            }
            if (Every_SE_All.size() != 0) {
                chr_Every_SE_all[chr_name] = Every_SE_All;
            }
            if (Every_SJ_begin_end_all.size() != 0) {
                chr_Every_SJ_begin_end_all[chr_name] = Every_SJ_begin_end_all;
            }
        }

        Known_SJ_SE.mSJs = chr_Every_SJs_all;
        Known_SJ_SE.SE = chr_Every_SE_all;
        Known_SJ_SE.mSJsBE = chr_Every_SJ_begin_end_all;
    }
    return Known_SJ_SE;
}





std::unordered_map<std::string, std::string> Read_fasta_file(std::string& fastafile_name){

    std::unordered_map<std::string, std::string> chrgeneString;

    std::vector<std::string> VecStrings;

    std::cout << "***** Now open Fasta file: " << fastafile_name << "! *****" << std::endl;

    std::ifstream FastaFile; 

	FastaFile.open(fastafile_name);

	if (!FastaFile.is_open()) {
		std::cout << "We can't open Fasta file !" << fastafile_name << std::endl;
		exit(EXIT_FAILURE);
	}

    std::string ChrGeneName;
    std::string concatenatedString;
    std::stringstream accumulatedStringStream;
    std::string line;
    std::string token;
    
    while (std::getline(FastaFile, line))
    {      
        if (line[0] == '>'){
            //第一行;
            if (!ChrGeneName.empty()){
                concatenatedString = accumulatedStringStream.str();
                chrgeneString[ChrGeneName] = concatenatedString;
                accumulatedStringStream.str("");
            }
            std::istringstream iss(line);
            int count = 0;
            while (std::getline(iss, token, ' ')){
                count++;

                if (count == 1){
                    ChrGeneName = token.substr(1);
                }
            }

        } else {

            accumulatedStringStream << line;
        }
    }

    if (!ChrGeneName.empty()) {
        concatenatedString = accumulatedStringStream.str();
        chrgeneString[ChrGeneName] = concatenatedString;
    }

    if (FastaFile.eof()) {

        std::cout << "***** The Fasta file has been read to the end ! *****" << std::endl;
    } 
    else if (FastaFile.fail()) {
        std::cout << "File FALSE !" << std::endl;
    }
    else {std::cout << "unkown reason" << std::endl;}

    //关闭文件;
    FastaFile.close();
    return chrgeneString;
}


bool ifSjNoError(const std::vector<std::array<int,2>>& leftSj, const std::vector<std::array<int,2>>& rightSj, int& Sj_Number) {
    bool flag;
    std::array<int,2> thisSj_left = leftSj[Sj_Number];
    std::array<int,2> thisSj_right = rightSj[Sj_Number];
    if ((thisSj_left[1]-thisSj_left[0]>10) && (thisSj_right[1]-thisSj_right[0]>10)) {
        flag = true;
    } else {
        flag = false;
    }
    return flag;
}


int ifSjSignal(std::string& Fasta_chr, const std::array<int,2>& Sj_Array, int& itsFlag) {
    int flag = 0;
    std::string signal1 = Fasta_chr.substr(Sj_Array[0]-1, 2);
    std::string signal2 = Fasta_chr.substr(Sj_Array[1]-2, 2);

    std::transform(signal1.begin(), signal1.end(), signal1.begin(), ::toupper);
    std::transform(signal2.begin(), signal2.end(), signal2.begin(), ::toupper);

    if (itsFlag == 1) {
        if (signal1 == "GT" && signal2 == "AG"){
            flag = 1;
        } else if (signal1 == "GC" && signal2 == "AG") {
            flag = 1;
        } else if (signal1 == "AT" && signal2 == "AC") {
            flag = 1;
        } else {
            flag = 0;
        }
    } else {
        if (signal1 == "CT" && signal2 == "AC") {
            flag = 2;
        } else if (signal1 == "CT" && signal2 == "GC") {
            flag = 2;
        } else if (signal1 == "GT" && signal2 == "AT") {
            flag = 2;
        } else {
            flag = 0;
        }
    }                
    return flag;
}



static struct option long_options[] = {
    {"sam", required_argument, 0, 's'},
    {"fasta", required_argument, 0, 'f'},
    {"gtf", required_argument, 0, 'g'},
    {"output", required_argument, 0, 'o'},
    {"SJDistance", required_argument, 0, 'j'},  
    {"support", required_argument, 0, 'n'},   
    {"single_exon_boundary", required_argument, 0, 'e'},
    {"graph_distance", required_argument, 0, 'd'},  
    {"thread", required_argument, 0, 't'}, 
    {"umi", required_argument, 0, 'u'},
    {"barcode", required_argument, 0, 'b'},
    {"output_min_read_count", required_argument, 0, 'r'}, 
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}  
};

bool directoryExists(const std::string& path) {
    struct stat info;

    if (stat(path.c_str(), &info) != 0) {

        return false;
    } else if (info.st_mode & S_IFDIR) {

        return true;
    } else {

        return false;
    }
}

bool createDirectory(const std::string& path) {

    if (mkdir(path.c_str(), 0755) == 0) {

        return true;
    } else {

        std::cerr << "Error creating directory: " << path << strerror(errno) << std::endl;
        return false;
    }
}

std::string joinPath(const std::string& directory, const std::string& filename) {
    
    if (directory.back() == '/') {
        return directory + filename;  
    } else {
        return directory + "/" + filename;  
    }
}


void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << std::endl;
    std::cout << "Required parameter: [-s <sam_file>] [-f <fasta_file>] [-o <output_file>] " << std::endl;
    std::cout << "Optional parameters: [-g <gtf_file>] [-j <SJDistance>] [-n <SJ_support_read_number>] [-d <Graph_distance>] [-t <Thread>] [-r <output_min_read_count>]" << std::endl;
    std::cout << "Help: [-h <help>]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -s, --sam                   SAM file path. Details of the 3 input methods can be found on GitHub. (required)" << std::endl;
    std::cout << "  -f, --fasta                 FASTA file path. FASTA file requires the chromosome names to match the GTF file. (required)" << std::endl;
    std::cout << "  -g, --gtf                   input annotation file in GTF format. (optional, Recommendation provided)" << std::endl;
    std::cout << "  -o, --output                output folder path. (required)" << std::endl;
    // std::cout << "  -m, --mode                  Sequencing method (default: 0(cDNA), !0(direct RNA))" << std::endl;
    std::cout << "  -j, --SJDistance            the minimum distance determined as intron. (optional, default:18)" << std::endl;
    std::cout << "  -n, --support               min perfect read count for all splice junctions of novel isoform. (optional, default:2)" << std::endl;
    std::cout << "  -e, --single_exon_boundary  belongs to the isoform scope of a single exon. (optional, default:60)" << std::endl;
    std::cout << "  -d, --graph_distance        the distance threshold for constructing the isoform candidate distance graph. (optional, default:60)" << std::endl;
    std::cout << "  -t, --thread                thread number (optional, default:8)." << std::endl;
    std::cout << "  -u, --umi                   umi_tag (optional, default:U8:Z:)." << std::endl;
    std::cout << "  -b, --barcode               barcode_tag (optional, default:BC:Z:)." << std::endl;
    std::cout << "  -r, --output_min_read_count the minimum number of transcripts in the output results. (optional, int, default:1)" << std::endl;
    std::cout << "  -h, --help                  show this help information." << std::endl;
}

std::vector<std::string> check_catalog_exist(const std::string& output_path) {
    std::vector<std::string> outputFileVector;

    if (!directoryExists(output_path)) {
        std::cout << "Directory does not exist, creating it..." << std::endl;

        if (createDirectory(output_path)) {
            std::cout << "Directory created successfully: " << output_path << std::endl;
        } else {
            std::cerr << "Failed to create directory: " << output_path << std::endl;
            exit(EXIT_FAILURE);
        }
    } else {
        std::cout << "Directory already exists: " << output_path << std::endl;
    } 

    std::string updatedGtfPath = joinPath(output_path, "updated_annotations.gtf");
    outputFileVector.push_back(updatedGtfPath);
    std::ofstream gtf_file(updatedGtfPath, std::ios::trunc);
    gtf_file.close();

    std::string TracePath = joinPath(output_path, "compatible_isoform.tsv");
    outputFileVector.push_back(TracePath);
    std::ofstream TraceIsoform(TracePath, std::ios::trunc);
    if (TraceIsoform.is_open()){
        TraceIsoform << "read_id" << '\t' << "category" << '\t' << "isoform_id" << '\t' << "gene_id" << '\t' << "barcode_id" << '\t' << "file" << '\n';
    }
    TraceIsoform.close();   
    return outputFileVector;
}

static inline bool ends_with(const std::string& s, const std::string& suf) {
    return s.size() >= suf.size() &&
           s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

static inline std::string trim_copy(std::string x) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    x.erase(x.begin(), std::find_if(x.begin(), x.end(), not_space));
    x.erase(std::find_if(x.rbegin(), x.rend(), not_space).base(), x.end());
    return x;
}

std::vector<std::string> traverse_sam_file(const std::string& sam_file_path, const std::string& output_path){
    std::vector<std::string> sam_file_vector;
    struct stat sam_stat;

    if (stat(sam_file_path.c_str(), &sam_stat) != 0) {
        std::cerr << "* Not a valid file or folder! * " << sam_file_path << " : " << std::strerror(errno) << "\n";
        return sam_file_vector;
    }

    if (S_ISREG(sam_stat.st_mode)) {
        if (ends_with(sam_file_path, ".txt") || ends_with(sam_file_path, ".tsv")) {
            std::cerr << "* This is a txt/tsv file! * << " << sam_file_path << std::endl;
            std::ifstream infile(sam_file_path);
            if (!infile) {
                std::cerr << "The file cannot be opened: " << sam_file_path << " : " << std::strerror(errno) << "\n";
                return sam_file_vector;
            }       

            std::string line;
            while (std::getline(infile, line)) {
                line = trim_copy(line);
                if (line.empty()) continue;
                if (!line.empty() && line[0] == '#') continue;
                sam_file_vector.push_back(line);
            }
            infile.close();
            if (sam_file_vector.empty()) {
                std::cerr << "^-^ There are 0 sam files in total. ^-^\n";
            } else {
                std::cerr << "^-^ There are " << sam_file_vector.size() << " sam files in total. ^-^\n";
            }
        } else if (ends_with(sam_file_path, ".sam")) {
            std::cerr << "* Only one sam file is entered! * " << sam_file_path << "\n";
            sam_file_vector.push_back(sam_file_path);
        } else {
            std::cerr << "* Not a valid file type (expect .sam/.txt/.tsv)! * " << sam_file_path << "\n";
            return sam_file_vector;
        }
        
    } else if (S_ISDIR(sam_stat.st_mode)) {
        std::cerr << "* A folder was entered! * " << sam_file_path << "\n";
        DIR* dir = opendir(sam_file_path.c_str());
        if (!dir) {
            std::cerr << "opendir failed: " << sam_file_path << " : " << std::strerror(errno) << "\n";
            return sam_file_vector;
        }   
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            if (!entry->d_name || entry->d_name[0] == '.') continue;
            std::string fileName = entry->d_name;
            if (ends_with(fileName, ".sam")) {
                sam_file_vector.push_back(fileName);
            }
        }
        closedir(dir);
        std::cerr << "^-^ There are " << sam_file_vector.size() << " sam files in total. ^-^\n";
    } else {
        std::cerr << "* Not a valid file or folder! * " << sam_file_path << "\n";
        return sam_file_vector;
    }

    std::string File_explain_path = joinPath(output_path, "file_explain.txt");
    std::ofstream Explain_file(File_explain_path, std::ios::trunc);    
    if (!Explain_file) {
        std::cerr << "Cannot write explain file: " << File_explain_path << " : " << std::strerror(errno) << "\n";
        return sam_file_vector;
    }   

    Explain_file << "File" << '\t' << "File_Path" << '\n';
    for (int i = 0; i < sam_file_vector.size(); i++) {
        std::cout << "SAM File " << i << " : " << sam_file_vector[i] << std::endl;
        Explain_file << i << '\t' << sam_file_vector[i] << '\n';
    }
    Explain_file.close();
    return sam_file_vector;
}


std::streampos findNextLineStart(std::ifstream& f, std::streampos pos) {
    if (pos <= std::streampos(0)) return std::streampos(0);
    f.clear();
    // 已经在行首;
    f.seekg(pos - std::streamoff(1));
    char c = '\0';
    if (f.get(c)) {
        if (c == '\n') {
            return pos;
        }
    } else {
        f.clear();
        return std::streampos(0);
    } 
    f.clear();
    // 跳过残行
    f.seekg(pos);
    std::string dummy;
    std::getline(f, dummy);          
    return f.tellg();                     
}


void processChunk(const std::string& one_sam_file_path, const std::streampos& start, const std::streampos& end, 
                const std::string& output_path, const int file_i, const int SJ_Distance, 
                std::unordered_map<std::string, std::string>& FastaRef, const int& mapq,
                const std::string& umi_tag_second, const std::string& barcode_tag_second){

    int Group_index = 0;
    std::string chrchr;
    Reads_Clusters each_cluster_informs;
    SpliceJs each_read_SJs_informs;
    std::array<int,2> readBeginEnd;
    
    bool thisSjNoError = 0;
    int thisSjSignal = 0;
    int SjNumber = 0; 

    std::ostringstream oss;
    oss << "Read_" << std::setw(4) << std::setfill('0') << file_i << ".txt";
    std::string File_name = oss.str();
    std::string ReadInformPath = joinPath(output_path, File_name);

    std::ofstream ReadInform(ReadInformPath.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
    if (!ReadInform.is_open()) {
        std::cerr << "Cannot open output file: " << ReadInformPath << std::endl;
        return;
    }

    std::ifstream samfile(one_sam_file_path.c_str());
    if (!samfile.is_open()) {
        std::cerr << "Cannot open SAM file: " << one_sam_file_path << std::endl;
        return;
    }
    samfile.seekg(start);

    std::streampos Current_Position = samfile.tellg();
    std::streampos Last_Position = Current_Position;

    // 大缓冲区：1MB 初始容量（根据内存可调）
    std::string outBuf;
    outBuf.reserve(10 * 1024 * 1024);  // 10 MB
    const std::size_t FLUSH_THRESHOLD = 10 * 1024 * 1024;  // 缓冲区达到 10MB 就写入磁盘

    while (Current_Position < end) {
        Group_index++;
        each_cluster_informs = get_each_cluster_reads_sc(samfile, Last_Position, end, mapq, umi_tag_second, barcode_tag_second);
        chrchr = each_cluster_informs.SetRef_name;
        Last_Position = each_cluster_informs.lastPos;
        Current_Position = each_cluster_informs.newPos;
        each_read_SJs_informs = get_reads_allSJs(each_cluster_informs.Mymap, SJ_Distance);  
                                                                                            
        for (const auto& eachRead:each_read_SJs_informs.reads_SJs) {
            readBeginEnd = each_read_SJs_informs.reads_begin_end[eachRead.first];
            ReadInform << eachRead.first << '\t' << chrchr << '\t' << Group_index << '\t'
                        << each_cluster_informs.ClusterCoverage[0] << '\t' 
                        << each_cluster_informs.ClusterCoverage[1] << '\t' << readBeginEnd[0] << '\t' 
                        << readBeginEnd[1] << '\t' << each_cluster_informs.Mylen[eachRead.first] << '\t'
                        << each_cluster_informs.ReadsBarcode[eachRead.first] << '\t'
                        << each_cluster_informs.ReadsUMI[eachRead.first] << '\t'
                        << each_cluster_informs.MyFlag[eachRead.first];
            
            SjNumber = 0;
            for (const auto& Sj:eachRead.second) {
                ReadInform << '\t' << Sj[0] << '\t' << Sj[1];

                thisSjNoError = ifSjNoError(each_read_SJs_informs.reads_SJs_left[eachRead.first], each_read_SJs_informs.reads_SJs_right[eachRead.first], SjNumber); 
                if (thisSjNoError == 1) {

                    thisSjSignal = ifSjSignal(FastaRef[chrchr], Sj, each_cluster_informs.MyFlag[eachRead.first]);

                    if (thisSjSignal == 1) {
                        ReadInform << '\t' << 1;
                    } else if (thisSjSignal == 2) {
                        ReadInform << '\t' << 2;
                    } else {
                        ReadInform << '\t' << 0;
                    }
                } else {
                    ReadInform << '\t' << 0;
                }
                SjNumber++;                    
            }
            ReadInform << '\n';
        }
        
        for (const auto& eachRead:each_read_SJs_informs.reads_single_exon) {
            readBeginEnd = eachRead.second;
            ReadInform << eachRead.first << '\t' << chrchr << '\t' << Group_index << '\t'
                        << each_cluster_informs.ClusterCoverage[0] << '\t' 
                        << each_cluster_informs.ClusterCoverage[1] << '\t' << readBeginEnd[0] << '\t' 
                        << readBeginEnd[1] << '\t' << each_cluster_informs.Mylen[eachRead.first] << '\t'
                        << each_cluster_informs.ReadsBarcode[eachRead.first] << '\t'
                        << each_cluster_informs.ReadsUMI[eachRead.first] << '\t'
                        << each_cluster_informs.MyFlag[eachRead.first] << '\n';         
        } 
    }
    samfile.close();
    ReadInform.close();
    // std::cout << "^-^ Thread: " << std::this_thread::get_id() << " has completed processing! ^-^" << std::endl;
    std::cerr << "^-^ Thread: " << file_i << " has completed processing! ^-^" << std::endl;
}



Split_Result get_line_split(const std::string& s, char delimiter){

    Split_Result read_result = {};

    std::istringstream iss(s);

    std::string token;
  
    int number = 0;
    while (std::getline(iss, token, delimiter)){
        number++;
        if (number == 1) {
            read_result.read_name = token;
        } else {
            read_result.tokens.push_back(token);
        }
    }
    return read_result;
}


Split_Result get_sj_split(const std::string& s, char delimiter){

    Split_Result sj_result = {};

    std::istringstream iss(s);

    std::string token;

    while (std::getline(iss, token, delimiter)){
        sj_result.tokens.push_back(token);
    }
    return sj_result;
}



struct FileSplit {

    std::map<std::string, std::vector<std::array<int,2>>> chr_coverage; 
    std::map<std::array<int,2>, std::array<std::streampos,2>> coverage2pos; 

    std::vector<std::streampos> reads_pointer; 
    std::vector<int> group_reads_number; 
    // 输出的文件;
    std::string readtxt_path; 
    std::map<std::string, std::set<std::string>> file_barcodeSet;
    int FileNo;
};

template <typename T>
std::vector<std::size_t> sort_indexes_e(std::vector<T> &v)
{
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
    return idx;
}



FileSplit Merge_Read_Small_Files(const std::string& SmallFilePath, const int& SAMFileNumber) {
    FileSplit Chunk_Bang;

    std::vector<std::string> Read_x_vec;

    DIR* dir = opendir(SmallFilePath.c_str());
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_name[0] == 'R') { 
            Read_x_vec.push_back(entry->d_name);
        } 
    }
    std::sort(Read_x_vec.begin(), Read_x_vec.end()); 
    closedir(dir);

    std::streampos Readpos;
    std::string File_total_name = SmallFilePath + "/All_Read.txt"; 
    std::ofstream AllReadInform(File_total_name, std::ios::trunc);  
    Chunk_Bang.readtxt_path = File_total_name;
    Readpos = AllReadInform.tellp();
    Chunk_Bang.reads_pointer.push_back(Readpos);

    std::string line;
    std::string chr_name;
    std::string last_chr;
    std::string earlyGroup = "+";
    std::string thisGroup = "+";
    int Group_index = 0;
    int Group_low;
    int Group_high;

    std::set<std::string> cellSet;

    std::map<std::string, std::vector<std::string>> Group_temporary;
    std::array<int,2> this_coverage;

    for (int file_number = 0; file_number < Read_x_vec.size(); file_number++) {

        std::string fileName = SmallFilePath + "/" + Read_x_vec[file_number];
        std::ifstream SmallSamFile(fileName);
        earlyGroup = "+";
        thisGroup = "+";

        while (getline(SmallSamFile, line)) {
            Split_Result This_Line = get_line_split(line, '\t');
            cellSet.insert(This_Line.tokens[7]);
            thisGroup = This_Line.tokens[1];
            this_coverage = {std::stoi(This_Line.tokens[2]), std::stoi(This_Line.tokens[3])};

            if (thisGroup == earlyGroup) {
                auto it = Group_temporary.find(This_Line.read_name);
                if (it != Group_temporary.end()) {
                    if (std::stoi(This_Line.tokens[6]) > std::stoi(Group_temporary[This_Line.read_name][6])){
                        Group_temporary[This_Line.read_name] = This_Line.tokens;
                    }
                } else {
                    Group_temporary[This_Line.read_name] = This_Line.tokens;
                }
            } else { 

                if (Group_temporary.size() != 0) {

                    if (last_chr == This_Line.tokens[0]) {

                        if (Group_low <= this_coverage[1] && this_coverage[0] <= Group_high) {

                            auto it = Group_temporary.find(This_Line.read_name);
                            if (it != Group_temporary.end()){
                                if (std::stoi(This_Line.tokens[6]) > std::stoi(Group_temporary[This_Line.read_name][6])) {
                                    Group_temporary[This_Line.read_name] = This_Line.tokens;
                                }
                            } else {
                                Group_temporary[This_Line.read_name] = This_Line.tokens;
                            }
                            Group_high = this_coverage[1]>Group_high ? this_coverage[1] : Group_high;
                            Group_low = this_coverage[0]<Group_low ? this_coverage[0] : Group_low;
                        } else {

                            Group_index++;
                            for (const auto& Read:Group_temporary){
                                AllReadInform << SAMFileNumber << '\t' << Read.first;
                                for (int j = 0; j < Read.second.size(); j++) {
                                    if (j == 1) {
                                        AllReadInform << '\t' << Group_index;
                                    } else if (j == 2) {
                                        AllReadInform << '\t' << Group_low;
                                    } else if (j == 3) {
                                        AllReadInform << '\t' << Group_high;
                                    } else {
                                        AllReadInform << '\t' << Read.second[j];
                                    }
                                }
                                chr_name = Read.second[0];
                                AllReadInform << '\n';
                            }

                            Readpos = AllReadInform.tellp();
                            Chunk_Bang.coverage2pos.insert({{Group_low, Group_high}, {Chunk_Bang.reads_pointer.back(),Readpos}});
                            Chunk_Bang.reads_pointer.push_back(Readpos); //这个Group结束;
                            Chunk_Bang.reads_pointer.push_back(Readpos); //下个Group开始;
                            Chunk_Bang.group_reads_number.push_back(Group_temporary.size());

                            auto it = Chunk_Bang.chr_coverage.find(chr_name);
                            if (it != Chunk_Bang.chr_coverage.end()){
                                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                            } else {
                                Chunk_Bang.chr_coverage[chr_name] = {};
                                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                            }
                            
                            Group_temporary.clear();
                            Group_temporary[This_Line.read_name] = This_Line.tokens;
                            Group_high = this_coverage[1];
                            Group_low = this_coverage[0];
                        }
                    } else {
                        Group_index++;
                        for (const auto& Read:Group_temporary){
                            AllReadInform << SAMFileNumber << '\t' << Read.first;
                            for (int j = 0; j < Read.second.size(); j++) {
                                if (j == 1) {
                                    AllReadInform << '\t' << Group_index;
                                } else if (j == 2) {
                                    AllReadInform << '\t' << Group_low;
                                } else if (j == 3) {
                                    AllReadInform << '\t' << Group_high;
                                } else {
                                    AllReadInform << '\t' << Read.second[j];
                                }
                            }
                            chr_name = Read.second[0];
                            AllReadInform << '\n';
                        }

                        Readpos = AllReadInform.tellp();
                        Chunk_Bang.coverage2pos.insert({{Group_low, Group_high}, {Chunk_Bang.reads_pointer.back(),Readpos}});
                        Chunk_Bang.reads_pointer.push_back(Readpos); //这个Group结束;
                        Chunk_Bang.reads_pointer.push_back(Readpos); //下个Group开始;     
                        Chunk_Bang.group_reads_number.push_back(Group_temporary.size());                     
                        
                        auto it = Chunk_Bang.chr_coverage.find(chr_name);
                        if (it != Chunk_Bang.chr_coverage.end()){
                            Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                        } else {
                            Chunk_Bang.chr_coverage[chr_name] = {};
                            Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                        }
                   
                        Group_temporary.clear();
                        Group_temporary[This_Line.read_name] = This_Line.tokens;    
                        Group_high = this_coverage[1];
                        Group_low = this_coverage[0];
                                                                  
                    }
                } else {
                    Group_temporary[This_Line.read_name] = This_Line.tokens;
                    Group_high = this_coverage[1];
                    Group_low = this_coverage[0];
                }
            }
            earlyGroup = thisGroup; 
            last_chr = This_Line.tokens[0];
        }

        if (file_number == (Read_x_vec.size()-1) && SmallSamFile.eof()){
            Group_index++;
            for (const auto& Read:Group_temporary){
                AllReadInform << SAMFileNumber << '\t' << Read.first;
                for (int j = 0; j < Read.second.size(); j++) {
                    if (j == 1) {
                        AllReadInform << '\t' << Group_index;
                    } else if (j == 2) {
                        AllReadInform << '\t' << Group_low;
                    } else if (j == 3) {
                        AllReadInform << '\t' << Group_high;
                    } else {
                        AllReadInform << '\t' << Read.second[j];
                    }
                }
                chr_name = Read.second[0];
                AllReadInform << '\n';
            }
            Readpos = AllReadInform.tellp();
            Chunk_Bang.coverage2pos.insert({{Group_low, Group_high}, {Chunk_Bang.reads_pointer.back(),Readpos}});
            Chunk_Bang.reads_pointer.push_back(Readpos);
            Chunk_Bang.group_reads_number.push_back(Group_temporary.size());
            
            auto it = Chunk_Bang.chr_coverage.find(chr_name);
            if (it != Chunk_Bang.chr_coverage.end()) {
                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
            } else {
                Chunk_Bang.chr_coverage[chr_name] = {};
                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
            }
        }
        SmallSamFile.close();
        std::remove(fileName.c_str());
    }
    AllReadInform.close();
    Chunk_Bang.FileNo = 1;
    Chunk_Bang.file_barcodeSet[std::to_string(SAMFileNumber)] = cellSet;
    return Chunk_Bang;
}


std::map<std::string, std::vector<std::array<int,2>>> Merge_Read_Interval(std::map<int, std::map<std::string, std::vector<std::array<int,2>>>>& FileCoverage) {

    std::map<std::string, std::vector<std::array<int,2>>> chr_range_temp;
    std::vector<std::array<int,2>> raw_intervals;
    std::vector<std::array<int,2>> sorted_intervals;
    if (FileCoverage.size() > 0) {

        for (const auto& eachFile:FileCoverage) {
            if (eachFile.second.size() > 0) {

                for (const auto& eachChr:eachFile.second) {
                    auto itchr = chr_range_temp.find(eachChr.first);
                    if (itchr != chr_range_temp.end()) {
                        chr_range_temp[eachChr.first].insert(chr_range_temp[eachChr.first].end(), eachChr.second.begin(), eachChr.second.end());
                    } else {
                        chr_range_temp[eachChr.first] = {};
                        chr_range_temp[eachChr.first].insert(chr_range_temp[eachChr.first].end(), eachChr.second.begin(), eachChr.second.end());
                    }
                }
            }
        }
    }

    std::map<std::string, std::vector<std::array<int,2>>> ChrInterval;
    if (chr_range_temp.size() > 0) {
        for (const auto& eachChr:chr_range_temp) {
            raw_intervals = eachChr.second;
            sort(raw_intervals.begin(), raw_intervals.end());
            sorted_intervals.clear();

            for (int i=0; i<raw_intervals.size();) {
                int left = raw_intervals[i][0];
                int right = raw_intervals[i][1];
                while (i+1 < raw_intervals.size() && right >= raw_intervals[i+1][0]){
                    right = std::max(right, raw_intervals[i+1][1]);
                    i++;
                }
                sorted_intervals.push_back({left, right});
                i++;
            }
            ChrInterval[eachChr.first] = sorted_intervals;
        }
    }
    return ChrInterval;
}


void Merge_Read_Big_Files_Group(const int& NEWNEW,
                                const std::string& output_path, 
                                const std::array<int,2>& Regin_BE,
                                std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>>& file_pointers,
                                const std::string& chrchr,
                                std::vector<std::streampos>& group_read_pointer,
                                std::vector<int>& group_read_number,
                                std::ofstream& BigFile) {

    std::streampos ThisFile_ThisGroup_pointer;
    std::streampos end;
    std::string line;
    std::map<int, std::vector<std::vector<std::string>>> Group_store;

    std::map<int, std::array<std::streampos,2>> This_Coverage_pointers = file_pointers[chrchr][Regin_BE];

    if (This_Coverage_pointers.size() > 0) {
        for (auto it = This_Coverage_pointers.begin(); it != This_Coverage_pointers.end(); ++it) {
            int i = it->first;
            Group_store[i] = {};
            const std::array<std::streampos, 2>& positions = it->second;
            // 遍历每个文件;
            std::string ReadFilePath = output_path + "/sam_" + std::to_string(i) + "/All_Read.txt";
            // 读取文件;
            std::ifstream samfile_one(ReadFilePath);
            // 文件跳入到这里;
            ThisFile_ThisGroup_pointer = positions[0];
            end = positions[1];
            samfile_one.seekg(ThisFile_ThisGroup_pointer);
            while (getline(samfile_one, line)) {
                Split_Result This_Line = get_sj_split(line, '\t');
                Group_store[i].push_back(This_Line.tokens);
                ThisFile_ThisGroup_pointer = samfile_one.tellg();
                if (ThisFile_ThisGroup_pointer < end) {
                    continue;
                } else {
                    break;
                }
            }
            samfile_one.close();
        }
    }

    int count = 0;
    int line_count = -1;
    if (BigFile.is_open()) 
    {
        {
            std::unique_lock<std::mutex> lock(bigMutex);
            for (const auto& eachFile:Group_store) {
                for (const auto& eachLine:eachFile.second) {
                    count++;
                    line_count = -1;
                    for (const auto& eachpart:eachLine) {
                        line_count = line_count + 1;
                        if (line_count == 0) {
                            BigFile << eachpart;
                            continue;
                        } else if (line_count == 3) {
                            BigFile << '\t' << NEWNEW;
                            continue;
                        } else if (line_count == 4) {
                            BigFile << '\t' << Regin_BE[0];
                            continue;
                        } else if (line_count == 5) {
                            BigFile << '\t' << Regin_BE[1];
                            continue;
                        } else {
                            BigFile << '\t' << eachpart;
                        }
                    }
                    BigFile << '\n';
                }
            }
            group_read_number.push_back(count);
            ThisFile_ThisGroup_pointer = BigFile.tellp();
            group_read_pointer.push_back(ThisFile_ThisGroup_pointer); //这个cluster结束;
            group_read_pointer.push_back(ThisFile_ThisGroup_pointer); //下个cluster开始;
        }
    } else {
        std::cout << "^-^ *** Merge large files failed *** ^-^" << std::endl;
    }
}



std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>> get_pointers(
        std::map<std::string, std::vector<std::array<int,2>>>& chrcoverage,
        std::map<int, std::map<std::string, std::vector<std::array<int,2>>>>& FileCoverage,
        std::map<int, std::map<std::array<int,2>, std::array<std::streampos,2>>>& groups_pointer) {

    std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>> chr_pointer_coverage;

    for (const auto& EachChr:chrcoverage) {
        chr_pointer_coverage[EachChr.first] = {};
        for (const auto& Region:EachChr.second) {
            chr_pointer_coverage[EachChr.first][Region] = {};
        }
    }

    std::array<std::streampos,2> thisFile_group_be;

    for (const auto& EachFile:FileCoverage) {

        int FileNumber = EachFile.first;
        std::map<std::array<int,2>, std::array<std::streampos,2>> thisFile_group_pointer = groups_pointer[FileNumber];

        for (const auto& EachChr:EachFile.second) {
            int Rhino = 0; //每个大区域索引;
            int Rhino_small = 0;
            std::vector<std::streampos> temp_pointer; 
            std::string ChrName = EachChr.first;
            std::vector<std::array<int,2>> Big_Rhino = chrcoverage[ChrName]; 

            for (const auto& small:EachChr.second) {

                Rhino_small = Rhino_small + 1;

                for (int i = Rhino; i < Big_Rhino.size(); i++) {
                    if (small[0] >= Big_Rhino[i][0] && small[1] <= Big_Rhino[i][1]) {

                        temp_pointer.push_back(thisFile_group_pointer[small][0]);
                        temp_pointer.push_back(thisFile_group_pointer[small][1]);

                        if (Rhino_small == EachChr.second.size()) {
                            thisFile_group_be[0] = temp_pointer[0];
                            thisFile_group_be[1] = temp_pointer[temp_pointer.size()-1];
                            chr_pointer_coverage[ChrName][Big_Rhino[i]].insert({FileNumber, thisFile_group_be});
                            temp_pointer.clear();
                        }
                        Rhino = i;
                        break;
                    } else if (small[0] > Big_Rhino[i][1]) {

                        if (temp_pointer.size() > 0) {
                            thisFile_group_be[0] = temp_pointer[0];
                            thisFile_group_be[1] = temp_pointer[temp_pointer.size()-1];
                            temp_pointer.clear();
                            chr_pointer_coverage[ChrName][Big_Rhino[i]].insert({FileNumber, thisFile_group_be});
                        }
                    } else {
                        continue;
                    }

                }
            }
        }
    }
    return chr_pointer_coverage;
}


FileSplit thread_all_read_sam_files(const std::string& sam_file_path, 
                                    std::vector<std::string>& samfilevec,
                                    const int& numThreads, 
                                    const std::string& outputPath, 
                                    const int& SJDistance, 
                                    std::unordered_map<std::string, std::string>& fasta_ref,
                                    const int& MAPQ, 
                                    const std::string& umi_tag_first,
                                    const std::string& barcode_tag_first) {
    
    std::vector<std::string> sam_file_vec;
    if (samfilevec.size() == 1) {
        sam_file_vec.push_back(sam_file_path);
    } else {
        for (const auto& eachFile:samfilevec) {
            if (eachFile.find('/') != std::string::npos) {
                sam_file_vec.push_back(eachFile);
            } else {
                std::string chunkFilePath = joinPath(sam_file_path, eachFile);
                sam_file_vec.push_back(chunkFilePath);
            }
        }
    }

    FileSplit BigBang;
    std::map<int, std::map<std::string, std::vector<std::array<int,2>>>> File_chr_coverage;
    std::map<int, std::map<std::array<int,2>, std::array<std::streampos,2>>> File_group_pointer;
    std::vector<std::streampos> startposVec;
    std::vector<std::streampos> endposVec;
    std::map<std::string, std::set<std::string>> All_File_All_Barcode;

    if (sam_file_vec.size() != 0) {
        ThreadPool Bigfilepool(numThreads);
        std::vector<std::future<void>> myJobs;
        for (int samFileNumber = 0; samFileNumber < sam_file_vec.size(); ++samFileNumber) {
            
            startposVec.clear();
            endposVec.clear();

            std::cerr << "*** Start processing SAM File " << samFileNumber << " ***" << std::endl;
            std::cerr << "***** " << sam_file_vec[samFileNumber] << " *****" << std::endl;
            
            std::string chunkFilePath = "sam_" + std::to_string(samFileNumber);
            chunkFilePath = joinPath(outputPath, chunkFilePath);

            if (!directoryExists(chunkFilePath)) {
                createDirectory(chunkFilePath);      
            }
            struct stat statBuf;
            if (stat(sam_file_vec[samFileNumber].c_str(), &statBuf) != 0){
                std::cerr << "We can't get file size !" << sam_file_vec[samFileNumber] << std::endl;
                exit(EXIT_FAILURE);
            }
            std::streampos fileSize = statBuf.st_size;

            std::streampos chunkSize = fileSize / numThreads;

            std::ifstream one_sam(sam_file_vec[samFileNumber]);

            if (!one_sam.is_open())	{
                std::cerr << "We can't open SAM file !" << sam_file_vec[samFileNumber] << std::endl;
                exit(EXIT_FAILURE);
            }

            for (int i = 0; i < numThreads; ++i) {
                std::streampos startPos = i * chunkSize;
                std::streampos endPos = (i == numThreads - 1) ? (fileSize) : static_cast<std::streampos>((i + 1) * chunkSize);

                if (i > 0) {
                    startPos = findNextLineStart(one_sam, startPos);
                }
                startposVec.push_back(startPos);
                endposVec.push_back(endPos);
            }
            one_sam.close();

            myJobs.clear();
            for (int i = 0; i < numThreads; ++i) {
                myJobs.emplace_back(Bigfilepool.enqueue([&,i]() {
                    processChunk(
                        sam_file_vec[samFileNumber],
                        startposVec[i], 
                        endposVec[i],
                        chunkFilePath, 
                        i, 
                        SJDistance, 
                        std::ref(fasta_ref),
                        MAPQ,
                        umi_tag_first,
                        barcode_tag_first
                    );
                }));
            }
            for (auto& future : myJobs) {
                future.get(); 
            }
            std::cerr << "^-^ [" << sam_file_vec[samFileNumber] << "] All threads are finished generating small files! ^-^" << std::endl;
            std::cerr << "^-^ Start of merge small files ! ^-^" << std::endl;
            BigBang = Merge_Read_Small_Files(chunkFilePath, samFileNumber);
            std::cout << "^-^ End of merge small files ! ^-^" << std::endl;

            if (sam_file_vec.size() > 1) {
                File_chr_coverage[samFileNumber] = BigBang.chr_coverage;
                File_group_pointer[samFileNumber] = BigBang.coverage2pos;
                All_File_All_Barcode[std::to_string(samFileNumber)] = BigBang.file_barcodeSet[std::to_string(samFileNumber)];
            }
            
        }
        if (sam_file_vec.size() > 1) {
            std::cout << "^-^ The number of sam files is greater than one. Large files need to be merged. ^-^" << std::endl;

            std::map<std::string, std::vector<std::array<int,2>>> ChrCoverage = Merge_Read_Interval(File_chr_coverage);

            std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>> ChrCoverage_SmallPointer = 
                                                                        get_pointers(ChrCoverage, File_chr_coverage, File_group_pointer);

            std::vector<std::streampos> Allreads_pointer;
            std::vector<int> Allreads_group_number;            
            // 输出的文件名;
            std::string FinallyFile = joinPath(outputPath, "MergedRead.txt");
            std::ofstream FinallyReadInform(FinallyFile, std::ios::trunc); 
            std::streampos Readpos = FinallyReadInform.tellp();
            Allreads_pointer.push_back(Readpos);

            myJobs.clear();
            int new_group = 0;

            std::cout << "^-^ Start of merge large File ! ^-^" << std::endl;
            for (const auto& eachChr:ChrCoverage_SmallPointer) {
                std::string ChrName = eachChr.first;
                // std::cout << ChrName << std::endl;
                for (const auto& eachRegion:eachChr.second) {
                    std::array<int,2> ourRegion = eachRegion.first;
                    new_group = new_group + 1;
                    myJobs.emplace_back(Bigfilepool.enqueue([new_group, ourRegion, ChrName, &outputPath, &ChrCoverage_SmallPointer, &Allreads_pointer, &Allreads_group_number, &FinallyReadInform]() {
                        Merge_Read_Big_Files_Group(
                            new_group,
                            outputPath,
                            ourRegion,
                            ChrCoverage_SmallPointer,
                            ChrName,
                            Allreads_pointer,
                            Allreads_group_number,
                            FinallyReadInform
                        );
                    }));
                }
            }
            for (auto& future : myJobs) {
                future.get();
            }
            std::cout << "^-^ End of merge large File ! ^-^" << std::endl;

            BigBang.readtxt_path = FinallyFile;
            BigBang.group_reads_number = Allreads_group_number;
            BigBang.reads_pointer = Allreads_pointer;
            BigBang.FileNo = sam_file_vec.size();
            BigBang.file_barcodeSet = All_File_All_Barcode;
            FinallyReadInform.close();
        }
    } else {
        std::cerr << "* There are no files in the folder! *" << strerror(errno) << std::endl;
    }
    return BigBang;
}



std::vector<std::unique_ptr<std::ofstream>> write_quantification_files(const std::string& outputPath, 
                                std::map<std::string, std::set<std::string>>& eachFileBarcode){
    std::vector<std::unique_ptr<std::ofstream>> geneisoformfiles;

    for (int filenumber = 0; filenumber < eachFileBarcode.size(); filenumber++) {
        std::string thisFile = std::to_string(filenumber);
        std::string isoformfile = "counts_transcript_" + thisFile + ".txt";
        isoformfile = joinPath(outputPath, isoformfile);
        std::unique_ptr<std::ofstream> Q_output_Transcript(new std::ofstream(isoformfile, std::ios::trunc));
        *Q_output_Transcript << "transcript_id" << '\t' << "gene_id";

        for (const auto& barcode:eachFileBarcode[thisFile]) {
            *Q_output_Transcript << '\t' << barcode;
        }
        *Q_output_Transcript << '\n';
        // 将 unique_ptr 转移到 vector 中
        geneisoformfiles.push_back(std::move(Q_output_Transcript));
    }
    return geneisoformfiles;
}


std::vector<std::unique_ptr<std::ofstream>> write_quantification_Gene_files(const std::string& outputPath, 
                                std::map<std::string, std::set<std::string>>& eachFileBarcode){
    std::vector<std::unique_ptr<std::ofstream>> geneisoformfiles;

    for (int filenumber = 0; filenumber < eachFileBarcode.size(); filenumber++) {
        std::string thisFile = std::to_string(filenumber);
        std::string isoformfile = "counts_gene_" + thisFile + ".txt";
        isoformfile = joinPath(outputPath, isoformfile);
        std::unique_ptr<std::ofstream> Q_output_Transcript(new std::ofstream(isoformfile, std::ios::trunc));
        *Q_output_Transcript << "gene_id";

        for (const auto& barcode:eachFileBarcode[thisFile]) {
            *Q_output_Transcript << '\t' << barcode;
        }
        *Q_output_Transcript << '\n';
        geneisoformfiles.push_back(std::move(Q_output_Transcript));
    }
    return geneisoformfiles;
}



struct Line_Split{
    std::string read_file; //0
    std::string read_name; //1
    std::string chr_name; //2
    std::string Group_index; //3
    std::array<int,2> Group_coverage; //4,5
    std::array<int,2> read_coverage; //6,7
    int read_length; // 这个就是追踪文件的时候可以画Coverage; //8
    std::string read_barcode; // 9
    std::string read_umi; // 10
    int read_strand; //11
    std::vector<std::array<int,2>> read_SJ;
    std::vector<int> read_sj_quality;
    
};


Line_Split MakeLineSplit(const std::string& s, char delimiter){

    Line_Split read_result = {};
    // 字符串流;
    std::istringstream iss(s);

    std::string token;

    int number = 0;

    std::vector<int> secondProcess;
    std::array<int,2> Sj;
    while (std::getline(iss, token, delimiter)){
        number = number + 1;
        if (number == 1) {
            read_result.read_file = token; 
        } else if (number == 2) {
            read_result.read_name = token;
        } else if (number == 3) {
            read_result.chr_name = token;
        } else if (number == 4) {
            read_result.Group_index = token;
        } else if (number == 5) {
            read_result.Group_coverage[0] = std::stoi(token);
        } else if (number == 6) {
            read_result.Group_coverage[1] = std::stoi(token);
        } else if (number == 7) {
            read_result.read_coverage[0] = std::stoi(token);
        } else if (number == 8) {
            read_result.read_coverage[1] = std::stoi(token);
        } else if (number == 9) {
            read_result.read_length = std::stoi(token);
        } else if (number == 10) {
            read_result.read_barcode = token; 
        } else if (number == 11) {
            read_result.read_umi = token;
        } else if (number == 12) {
            read_result.read_strand = std::stoi(token);
        } else if (number > 12) {
            secondProcess.push_back(std::stoi(token));
        }
    }

    for (size_t i = 0; i < secondProcess.size(); i += 3) {
        Sj = {secondProcess[i], secondProcess[i+1]};
        read_result.read_SJ.push_back(Sj);
        read_result.read_sj_quality.push_back(secondProcess[i+2]);
    }
    return read_result;
}



struct GroupInformation {
    std::string GroupIndex;
    std::string chrName; 
    std::array<int,2> GroupCoverage; 
    std::unordered_map<std::string, std::string> GroupReadBarcodes;
    std::unordered_map<std::string, std::string> GroupReadUMIs;
    std::unordered_map<std::string, std::string> GroupReadFiles; 
    std::unordered_map<std::string, std::vector<std::array<int,2>>> GroupReadSjs; 
    std::unordered_map<std::string, std::array<int,2>> GroupReadCoverage; 
    std::map<std::array<int,2>, int> GroupSjs;
    std::unordered_map<std::string, std::vector<int>> GroupSigns; 
    std::unordered_map<std::string, std::array<int,2>> GroupSingleExon; 
    std::unordered_map<std::string, int> GroupReadStrands;
};

GroupInformation knowGroupInformation(std::streampos& startpos, 
                                    std::streampos& endpos, 
                                    const std::string& sam_file_path, 
                                    const int& Sj_Support_Read_Number) {
    GroupInformation groupinformation;
    std::ifstream FinalSamFile(sam_file_path);
    FinalSamFile.seekg(startpos);

    std::streampos Current_Position = FinalSamFile.tellg();
    std::streampos Last_Position = Current_Position;
    std::string line;

    std::map<std::array<int,2>, std::array<int,2>> Temp_Sjs;

    while (getline(FinalSamFile, line)) {

        Last_Position = Current_Position;
        Current_Position = FinalSamFile.tellg();
        if (Last_Position < endpos || Last_Position == FinalSamFile.eof()) {

            Line_Split Line_result = MakeLineSplit(line, '\t');

            groupinformation.GroupIndex = Line_result.Group_index;
            groupinformation.chrName = Line_result.chr_name;
            groupinformation.GroupCoverage = Line_result.Group_coverage;
            groupinformation.GroupReadFiles[Line_result.read_name] = Line_result.read_file;

            groupinformation.GroupReadCoverage[Line_result.read_name] = Line_result.read_coverage;
            groupinformation.GroupReadBarcodes[Line_result.read_name] = Line_result.read_barcode;
            groupinformation.GroupReadUMIs[Line_result.read_name] = Line_result.read_umi;
            groupinformation.GroupReadStrands[Line_result.read_name] = Line_result.read_strand;

            if (Line_result.read_SJ.size() > 0) {
                groupinformation.GroupReadSjs[Line_result.read_name] = Line_result.read_SJ;
                groupinformation.GroupSigns[Line_result.read_name] = Line_result.read_sj_quality;

                for (int i = 0; i < Line_result.read_sj_quality.size(); i++) {
                    if (Line_result.read_sj_quality[i] > 0) {
                        std::array<int,2> Temp_array = Line_result.read_SJ[i];
                        auto it = Temp_Sjs.find(Temp_array);
                        if (it != Temp_Sjs.end()) {
                            Temp_Sjs[Temp_array][0] = Temp_Sjs[Temp_array][0] + 1;
                        } else {
                            Temp_Sjs[Temp_array] = {1, Line_result.read_sj_quality[i]};
                        }
                    }
                }
            } else {
                groupinformation.GroupSingleExon[Line_result.read_name] = Line_result.read_coverage;
            }
        } else {
            break;
        }
    }
    FinalSamFile.close();
    
    for (const auto& Sj:Temp_Sjs) {
        if (Sj.second[0] > Sj_Support_Read_Number) {
            groupinformation.GroupSjs[Sj.first] = Sj.second[1];
        }
    }
    return groupinformation;
}



struct GroupAnnotation
{   
    std::unordered_map<std::string, std::array<int,2>> group_genes; // 所有的基因;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> group_me_transcripts; //所有的多exon转录本; exon;
    std::unordered_map<std::string, std::array<int,2>> group_se_transcripts; //所有的单exon转录本;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> group_me_SJs; //所有的多exon转录本; SJ;
};

GroupAnnotation get_group_single_exon_gene_annotation (
                                std::unordered_map<std::string, std::array<int,2>>& AnnoGene, 
                                std::array<int,2>& Group_Coverage, 
                                std::unordered_map<std::string, std::vector<std::string>>& gene2tx,
                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& tx2exon,
                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& tx2sj,
                                std::string& CHR ) {
    GroupAnnotation thisGroupAnno; 
    thisGroupAnno.group_genes.reserve(AnnoGene.size() / 20);
    for (const auto& eachGene:AnnoGene) {
        if ( !( (eachGene.second[1] < Group_Coverage[0]) or (Group_Coverage[1] < eachGene.second[0]) ) ) {
            thisGroupAnno.group_genes.emplace(eachGene.first, eachGene.second);
        }
    }

    if (!thisGroupAnno.group_genes.empty()) {

        thisGroupAnno.group_me_transcripts.reserve(thisGroupAnno.group_genes.size() * 5);
        thisGroupAnno.group_me_SJs.reserve(thisGroupAnno.group_genes.size() * 5);
        thisGroupAnno.group_se_transcripts.reserve(thisGroupAnno.group_genes.size() * 3);

        for (const auto& eachGene:thisGroupAnno.group_genes) {
            auto it_gene_tx = gene2tx.find(eachGene.first);
            if (it_gene_tx == gene2tx.end()) continue;
            const auto& thisGene_txname = it_gene_tx->second;  // 引用;
            for (const auto& eachTx:thisGene_txname) {
                auto it_tx_exon = tx2exon.find(eachTx);
                if (it_tx_exon == tx2exon.end()) continue; 
                const auto& thistx_exon = it_tx_exon->second;
                if (thistx_exon.empty()) continue;
                if (thistx_exon.size() > 1) { //多exon的;
                    thisGroupAnno.group_me_transcripts.emplace(eachTx, thistx_exon);
                    auto it_tx_sj = tx2sj.find(eachTx);
                    if (it_tx_sj != tx2sj.end()) {
                        thisGroupAnno.group_me_SJs.emplace(eachTx, it_tx_sj->second);
                    }
                } else {
                    thisGroupAnno.group_se_transcripts.emplace(eachTx, thistx_exon[0]);
                }
            }
        }
    }

    return thisGroupAnno;
}



void write_single_exon_gtf_trace_sc(int& FileNo,
                                    std::unordered_map<std::string, std::string>& groupreadfiles,
                                    std::unordered_map<std::string, std::vector<std::string>>& singleexonwithreads,
                                    std::unordered_map<std::string, std::array<int,2>>& singleexongroupannotation,
                                    std::ofstream& Updated_Files, std::ofstream& Trace, 
                                    std::vector<std::unique_ptr<std::ofstream>>& IsoformFilePath,
                                    std::vector<std::mutex>& IsoformMutexes,
                                    std::unordered_map<std::string, std::string>& GTF_Transcript_Strand,
                                    std::unordered_map<std::string, std::string>& groupreadbarcodes,
                                    std::string& chrname,
                                    std::map<std::string, std::set<std::string>>& AllFileBarcode) {

    std::string updateBuffer;
    std::vector<std::string> isoformBuffers(FileNo);

    std::unordered_map<std::string, std::vector<std::string>> me_transcript_se_reads;
    me_transcript_se_reads.reserve(singleexonwithreads.size());

    std::vector<std::vector<std::string>> AllBarcodeVec(FileNo);
    for (int f = 0; f < FileNo; f++) {
        const std::string key = std::to_string(f);
        AllBarcodeVec[f].assign(AllFileBarcode[key].begin(), AllFileBarcode[key].end());
    }

    for (const auto& kv : singleexonwithreads) {
        const std::string& transcript = kv.first;
        const std::vector<std::string>& reads = kv.second;

        auto itAnno = singleexongroupannotation.find(transcript);

        if (itAnno == singleexongroupannotation.end()) {
            if (!reads.empty()) me_transcript_se_reads.emplace(transcript, reads);
            continue;
        }
        if (reads.empty()) continue;

        size_t pos = transcript.find('|');
        std::string first_part, second_part;
        if (pos != std::string::npos) {
            first_part  = transcript.substr(0, pos);
            second_part = transcript.substr(pos + 1);
        } else {
            first_part = second_part = transcript;
        }

        const auto& exonRange = itAnno->second;
        const std::string& strand = GTF_Transcript_Strand[transcript];

        updateBuffer.append(chrname); updateBuffer.push_back('\t');
        updateBuffer.append("annotated_isoform\ttranscript\t");
        updateBuffer.append(std::to_string(exonRange[0])); updateBuffer.push_back('\t');
        updateBuffer.append(std::to_string(exonRange[1]));
        updateBuffer.append("\t.\t"); updateBuffer.append(strand);
        updateBuffer.append("\t.\tgene_id \""); updateBuffer.append(first_part);
        updateBuffer.append("\"; transcript_id \""); updateBuffer.append(second_part);
        updateBuffer.append("\";\n");

        updateBuffer.append(chrname); updateBuffer.push_back('\t');
        updateBuffer.append("annotated_isoform\texon\t");
        updateBuffer.append(std::to_string(exonRange[0])); updateBuffer.push_back('\t');
        updateBuffer.append(std::to_string(exonRange[1]));
        updateBuffer.append("\t.\t"); updateBuffer.append(strand);
        updateBuffer.append("\t.\tgene_id \""); updateBuffer.append(first_part);
        updateBuffer.append("\"; transcript_id \""); updateBuffer.append(second_part);
        updateBuffer.append("\"; exon_number \"0\";\n");

        std::vector<std::vector<int>> counts(FileNo);
        for (int f = 0; f < FileNo; f++)
            counts[f].assign(AllBarcodeVec[f].size(), 0);

        for (const auto& read : reads) {
            int fileID = std::stoi(groupreadfiles[read]);
            const std::string& barcode = groupreadbarcodes[read];

            auto& barVec = AllBarcodeVec[fileID];
            for (size_t bi = 0; bi < barVec.size(); bi++) {
                if (barVec[bi] == barcode) {
                    counts[fileID][bi] += 1;
                    break;
                }
            }
        }

        for (int f = 0; f < FileNo; f++) {
            isoformBuffers[f].append(second_part);
            isoformBuffers[f].push_back('\t');
            isoformBuffers[f].append(first_part);
            for (auto c : counts[f]) {
                isoformBuffers[f].push_back('\t');
                isoformBuffers[f].append(std::to_string(c));
            }
            isoformBuffers[f].push_back('\n');
        }
    }

    if (Updated_Files.is_open()) {
        std::unique_lock<std::mutex> lock(updatedGtfMutex);
        Updated_Files << updateBuffer;
    }

    for (int f = 0; f < FileNo; f++) {
        std::unique_lock<std::mutex> lock(IsoformMutexes[f]);
        *(IsoformFilePath[f]) << isoformBuffers[f];

    }
    singleexonwithreads = std::move(me_transcript_se_reads);

}


struct ReadNode {
    std::string read_id;
    std::string umi;
    size_t umi_len;
    bool is_duplicate;
    // 使用 std::move 优化构造函数
    ReadNode(std::string id, std::string u) 
        : read_id(std::move(id)), umi(std::move(u)), is_duplicate(false) {
        umi_len = umi.length();
    }
};

template <typename KeyType>
void delete_umi_replicate(std::unordered_map<KeyType, std::vector<std::string>>& groupcluster,
                        const std::unordered_map<std::string, std::string>& groupreadbarcodes,
                        const std::unordered_map<std::string, std::string>& groupreadumi) {
    
    // 移除 result_cluster，不再需要双倍内存
    EdlibAlignConfig align_config = edlibNewAlignConfig(3, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);

    // 使用迭代器遍历，方便删除元素
    auto it_cluster = groupcluster.begin();
    while (it_cluster != groupcluster.end()) {
        const std::vector<std::string>& read_list = it_cluster->second;

        // 1. 构建 Barcode 分组
        std::unordered_map<std::string, std::vector<ReadNode>> barcode_groups;
        // 预估大小，避免扩容
        if (read_list.size() > 0) barcode_groups.reserve(read_list.size()); 

        for (const auto& read_id : read_list) {
            auto it_bc = groupreadbarcodes.find(read_id);
            auto it_umi = groupreadumi.find(read_id);
            if (it_bc != groupreadbarcodes.end() && it_umi != groupreadumi.end()) {
                barcode_groups[it_bc->second].emplace_back(read_id, it_umi->second);
            }
        }

        std::vector<std::string> kept_reads_in_cluster;
        kept_reads_in_cluster.reserve(read_list.size());

        // 2. 遍历每个Barcode组进行聚类;
        for (auto& bg_pair : barcode_groups) {
            std::vector<ReadNode>& nodes = bg_pair.second;
            const size_t num_reads = nodes.size();
            // Case 1: 单条read, 直接保留;
            if (num_reads == 1) {
                kept_reads_in_cluster.push_back(std::move(nodes[0].read_id));
                continue;
            }
            // Case 2:多条reads, 贪婪聚类;
            for (size_t i = 0; i < num_reads; ++i) {
                if (nodes[i].is_duplicate) continue;

                // 保留中心节点
                kept_reads_in_cluster.push_back(std::move(nodes[i].read_id));
                
                const char* center_umi_ptr = nodes[i].umi.c_str();
                const int center_len = static_cast<int>(nodes[i].umi_len);

                for (size_t j = i + 1; j < num_reads; ++j) {
                    if (nodes[j].is_duplicate) continue;
                    const int other_len = static_cast<int>(nodes[j].umi_len);
                    if (std::abs(center_len - other_len) > 3) continue;
                    // Edlib 比对
                    EdlibAlignResult alignResult = edlibAlign(
                        center_umi_ptr, center_len,
                        nodes[j].umi.c_str(), other_len,
                        align_config
                    );

                    if (alignResult.status == EDLIB_STATUS_OK && alignResult.editDistance != -1) nodes[j].is_duplicate = true;
                    edlibFreeAlignResult(alignResult);
                }
            }
        }

        if (!kept_reads_in_cluster.empty()) {
            // 将去重后的reads移动回原来的map value中, 旧的vector内容被释放, 新的内容接管, 内存波动极小;
            it_cluster->second = std::move(kept_reads_in_cluster);
            ++it_cluster;
        } else {
            it_cluster = groupcluster.erase(it_cluster);
        }
    }
}



struct SE_belong2_genetranscript {
    std::unordered_map<int, std::unordered_map<std::string, std::vector<std::string>>> file_SE_reads_gene_number;
    std::unordered_map<std::string, std::vector<std::string>> Transcript_with_SE_reads;
};


SE_belong2_genetranscript get_group_singleexon_reads_2gene(
                                    GroupAnnotation& thisGroupAnnotation,
                                    std::unordered_map<std::string, std::array<int,2>>& GroupSingleExonReads,
                                    std::unordered_map<std::string, std::string>& groupreadfiles,
                                    std::unordered_map<std::string, std::vector<std::string>>& gene2tx,
                                    std::unordered_map<std::string, std::vector<std::array<int,2>>>& tx2exon,
                                    int& FileNo, int& Edge, std::string& CHR,
                                    std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>>& AlreadyFSM,
                                    std::ofstream& Trace,
                                    std::unordered_map<std::string, std::string>& groupreadbarcodes,
                                    std::unordered_map<std::string, std::string>& groupreadumis) {
    
    SE_belong2_genetranscript SingleReads;
    
    const auto& group_genes = thisGroupAnnotation.group_genes;
    const auto& group_me_transcripts = thisGroupAnnotation.group_me_transcripts;
    const auto& group_se_transcripts = thisGroupAnnotation.group_se_transcripts;

    std::unordered_map<std::string, std::vector<std::string>> SingleReads_inGenes;
    SingleReads_inGenes.reserve(group_genes.size());

    for (const auto& EachGene:group_genes) SingleReads_inGenes.emplace(EachGene.first, std::vector<std::string>{});

    std::unordered_map<std::string, std::vector<std::string>> SingleReads_inTranscripts;
    SingleReads_inTranscripts.reserve(group_me_transcripts.size() + group_se_transcripts.size());

    for (const auto& EachTx:group_me_transcripts) SingleReads_inTranscripts.emplace(EachTx.first, std::vector<std::string>{});
    for (const auto& EachTx:group_se_transcripts) SingleReads_inTranscripts.emplace(EachTx.first, std::vector<std::string>{});

    std::unordered_map<std::string, std::vector<std::array<int,2>>> genes2exons; genes2exons.reserve(group_genes.size());
    for (const auto& g : group_genes) {
        auto it_gene = gene2tx.find(g.first);
        if (it_gene == gene2tx.end()) continue;
        std::vector<std::array<int,2>> exons;
        for (const std::string& tx : it_gene->second) {
            auto it_tx = tx2exon.find(tx);
            if (it_tx == tx2exon.end()) continue;
            exons.insert(exons.end(), it_tx->second.begin(), it_tx->second.end());
        }
        exons = mergeIntervals(exons);
        genes2exons.emplace(g.first, std::move(exons));
    }    

    std::vector<std::string> Temp_Gene; Temp_Gene.reserve(32);
    std::vector<std::array<int,2>> tmpReadVec(1);

    std::vector<int> edgeDist;
    
    for (const auto& EachSing:GroupSingleExonReads) {

        const auto& read_name = EachSing.first;
        const auto& read_exon = EachSing.second;
        Temp_Gene.clear();
        tmpReadVec[0] = read_exon;

        for (const auto& EachGene:group_genes) {
            if ( read_exon[0] > EachGene.second[1] + Edge || read_exon[1] < EachGene.second[0] - Edge ) continue;
            if (IntervalIntersection(genes2exons[EachGene.first], tmpReadVec) > 0) Temp_Gene.push_back(EachGene.first);       
        }

        if (Temp_Gene.size() == 0) continue;
        
        if (Temp_Gene.size() == 1) {
            const auto& gene = Temp_Gene[0];
            auto& SRiG_gene = SingleReads_inGenes[gene];
            SRiG_gene.push_back(read_name);

        } else {
            // 如果有多个基因, 看与哪个基因重复的更多;
            edgeDist.clear();
            edgeDist.reserve(Temp_Gene.size());
            
            std::vector<int> Temp_Dist; Temp_Dist.reserve(Temp_Gene.size());
            std::vector<int> Temp_exon; Temp_exon.reserve(Temp_Gene.size());
            
            for (const auto& gene_id : Temp_Gene) {
                auto it_gene_exons = genes2exons.find(gene_id);
                auto it_gene_txs = gene2tx.find(gene_id);
                if (it_gene_exons == genes2exons.end() || it_gene_txs == gene2tx.end()) {
                    Temp_Dist.push_back(0);
                    Temp_exon.push_back(0);
                    edgeDist.push_back(INT_MAX);
                    continue;
                }
                const auto& exons = it_gene_exons->second;
                const auto& transcripts = it_gene_txs->second;

                int min_dist = INT_MAX;
                int intersection = IntervalIntersection(exons, tmpReadVec);
                Temp_Dist.push_back(intersection);
                Temp_exon.push_back(exons.size());

                for (const auto& EachTx:transcripts) {
                    auto it_tx = tx2exon.find(EachTx);
                    if (it_tx == tx2exon.end()) continue;
                    int d = IntervalMinDistance(it_tx->second, read_exon);
                    if (d < min_dist) min_dist = d;
                    if (min_dist == 0) break; 
                }
                edgeDist.push_back(min_dist);
            }

            std::string MaxGeneName;         
            // 1.首先看看基因边界;
            std::vector<int> Next_Gene_Part;
            int min_gene_interval = *std::min_element(edgeDist.begin(), edgeDist.end());
            for (int i = 0; i < Temp_Gene.size(); i++) {
                if (edgeDist[i] - min_gene_interval < 10) {
                    Next_Gene_Part.push_back(i);
                }
            }
            if (Next_Gene_Part.size() == 1) {
                MaxGeneName = Temp_Gene[Next_Gene_Part[0]];
            } else {
                // 2.然后看看总体长度;
                std::vector<int> Next2_Gene_Part;
                std::vector<int> Temp_Dist2;
                for (int i = 0; i < Next_Gene_Part.size(); i++) {
                    Temp_Dist2.push_back(Temp_Dist[Next_Gene_Part[i]]);
                }
                int max_Gene_intersect = *std::max_element(Temp_Dist2.begin(), Temp_Dist2.end());
                for (int j = 0; j < Next_Gene_Part.size(); j++) {
                    if (max_Gene_intersect - Temp_Dist[Next_Gene_Part[j]] < 10) {
                        Next2_Gene_Part.push_back(j);
                    }
                }
                if (Next2_Gene_Part.size() == 1) {
                    MaxGeneName = Temp_Gene[Next2_Gene_Part[0]];
                } else {
                    std::vector<int> Next3_Gene_Part;
                    for (int j = 0; j < Next2_Gene_Part.size(); j++) {
                        if (Temp_exon[Next2_Gene_Part[j]] == 1) {
                            Next3_Gene_Part.push_back(Next2_Gene_Part[j]);
                        }
                    }
                    if (Next3_Gene_Part.size() == 1) {
                        MaxGeneName = Temp_Gene[Next3_Gene_Part[0]];
                    } else {
                        MaxGeneName = Temp_Gene[Next2_Gene_Part[0]];
                    }
                }
            }
            auto& SRiG_gene = SingleReads_inGenes[MaxGeneName];
            SRiG_gene.push_back(read_name);
        }
    }

    delete_umi_replicate(SingleReads_inGenes, groupreadbarcodes, groupreadumis);

    std::vector<int> Tx_Temp_Dist; 
    std::vector<int> Tx_Temp_exon; 
    std::vector<int> Tx_Temp_Interval;
    std::vector<std::string> Part_Transcript;

    for (const auto& EachGene:SingleReads_inGenes) {
        if (EachGene.second.size() > 0) {
            const auto& gene = EachGene.first;
            const auto& Temp_Transcript = gene2tx[gene];

            for (const auto& read_name:EachGene.second) {
                auto itread = GroupSingleExonReads.find(read_name);
                if (itread == GroupSingleExonReads.end()) continue;
                const auto& read_exon = itread->second;
                if (Temp_Transcript.size() == 1) {
                    // 只有一个转录本;
                    SingleReads_inTranscripts[Temp_Transcript[0]].push_back(read_name);
                } else {
                    // 这一个基因有多个转录本; 首先是交集筛选;
                    Tx_Temp_Dist.clear();
                    Tx_Temp_exon.clear();
                    Part_Transcript.clear();
                    Tx_Temp_Interval.clear();
                    for (const auto& EachTx:Temp_Transcript) {
                        auto it_tx = tx2exon.find(EachTx);
                        if (it_tx == tx2exon.end()) continue;
                        const auto& SJ_Vector = it_tx->second;
                        int intersection = IntervalIntersection(SJ_Vector, tmpReadVec);
                        if (intersection > 0) {
                            Part_Transcript.push_back(EachTx);
                            Tx_Temp_Dist.push_back(intersection);
                            Tx_Temp_exon.push_back(SJ_Vector.size());
                            Tx_Temp_Interval.push_back(IntervalMinDistance(SJ_Vector, read_exon));
                        }
                    }
                    // 如果这多个转录本跟reads都没有交集;
                    if (Part_Transcript.size() == 0) continue;
                    if (Part_Transcript.size() == 1) {
                        // 如果只有一个转录本跟reads有交集;
                        if (Part_Transcript[0].size() != 0) {
                            SingleReads_inTranscripts[Part_Transcript[0]].push_back(read_name);
                        }                    
                    } else {
                        // 如果多个转录本跟reads有交集;
                        std::string MaxTranscriptName;
                        auto it = std::find(Tx_Temp_exon.begin(), Tx_Temp_exon.end(), 1);     
                        // 这里面有单exon的转录本;
                        if (it != Tx_Temp_exon.end()) {

                            // 1.首先看看边界
                            std::vector<int> Next_Part;
                            int min_interval = *std::min_element(Tx_Temp_Interval.begin(), Tx_Temp_Interval.end());
                            for (int i = 0; i < Part_Transcript.size(); i++) {
                                if (Tx_Temp_Interval[i] - min_interval < 10 ) {
                                    Next_Part.push_back(i);
                                }
                            }

                            if (Next_Part.size() == 1) {
                                if (Part_Transcript[Next_Part[0]].size() != 0) {
                                    SingleReads_inTranscripts[Part_Transcript[Next_Part[0]]].push_back(read_name);
                                }
                            
                            } else {
                                // 2.总体长度;
                                std::vector<int> Next2_Part;
                                std::vector<int> Dist2;
                                for (int i = 0; i < Next_Part.size(); i++) {
                                    Dist2.push_back(Tx_Temp_Dist[Next_Part[i]]);
                                }
                                int max_intersect = *std::max_element(Dist2.begin(), Dist2.end());

                                for (int j = 0; j < Next_Part.size(); j++) {
                                    if (max_intersect - Tx_Temp_Dist[Next_Part[j]] < 10) {
                                        Next2_Part.push_back(Next_Part[j]);
                                    }
                                }

                                if (Next2_Part.size() == 1) {
                                    if (Part_Transcript[Next2_Part[0]].size() != 0) {
                                        SingleReads_inTranscripts[Part_Transcript[Next2_Part[0]]].push_back(read_name);
                                    }

                                } else {
                                    // 3.是不是单exon的;
                                    std::vector<int> Next3_Part;
                                    for (int j = 0; j < Next2_Part.size(); j++) {
                                        if (Tx_Temp_exon[Next2_Part[j]] == 1) {
                                            Next3_Part.push_back(Next2_Part[j]);
                                        }
                                    }
                                    if (Next3_Part.size() == 1) {
                                        if (Part_Transcript[Next3_Part[0]].size() != 0) {
                                            SingleReads_inTranscripts[Part_Transcript[Next3_Part[0]]].push_back(read_name);
                                        }
                                        
                                    } else {
                                        // 看FSM;
                                        size_t positions = 0;
                                        int Thiscount = 0;
                                        for (size_t i = 0; i < Next3_Part.size(); i++) {
                                            if (AlreadyFSM.find(Part_Transcript[Next3_Part[i]]) != AlreadyFSM.end()) {
                                                if (AlreadyFSM[Part_Transcript[Next3_Part[i]]].second.size() > Thiscount) {
                                                    Thiscount = AlreadyFSM[Part_Transcript[Next3_Part[i]]].second.size();
                                                    positions = Next3_Part[i];
                                                }
                                            }
                                        }
                                        MaxTranscriptName = Part_Transcript[positions];                                 
                                        if (MaxTranscriptName.size() != 0) {
                                            SingleReads_inTranscripts[MaxTranscriptName].push_back(read_name);
                                        }
                                    }
                                }
                            }

                        // 这里面没有单exon的转录本;
                        } else {

                            // 整个基因没有单个exon的情况; 就要看与其他转录本的交集了;
                            // 1.首先看看边界
                            std::vector<int> Next_Part;
                            int min_interval = *std::min_element(Tx_Temp_Interval.begin(), Tx_Temp_Interval.end());
                            for (int i = 0; i < Part_Transcript.size(); i++) {
                                if ( Tx_Temp_Interval[i] - min_interval < 10 ) {
                                    Next_Part.push_back(i);
                                }
                            } 
                            if (Next_Part.size() == 1) {
                                if (Part_Transcript[Next_Part[0]].size() != 0) {
                                    SingleReads_inTranscripts[Part_Transcript[Next_Part[0]]].push_back(read_name);
                                }
                            
                            } else {
                                // 2.总体长度;
                                std::vector<int> Next2_Part;
                                std::vector<int> Dist2;
                                for (int i = 0; i < Next_Part.size(); i++) {
                                    Dist2.push_back(Tx_Temp_Dist[Next_Part[i]]);
                                }
                                int max_intersect = *std::max_element(Dist2.begin(), Dist2.end());

                                for (int j = 0; j < Next_Part.size(); j++) {
                                    if (max_intersect - Tx_Temp_Dist[Next_Part[j]] < 10) {
                                        Next2_Part.push_back(Next_Part[j]);
                                    }
                                }
                                if (Next2_Part.size() == 1) {
                                    if (Part_Transcript[Next2_Part[0]].size() != 0) {
                                        SingleReads_inTranscripts[Part_Transcript[Next2_Part[0]]].push_back(read_name);
                                    }
                                } else {
                                    // 看FSM;
                                    size_t positions = 0;
                                    int Thiscount = 0;
                                    for (size_t i = 0; i < Next2_Part.size(); i++) {
                                        if (AlreadyFSM.find(Part_Transcript[Next2_Part[i]]) != AlreadyFSM.end()) {
                                            if (AlreadyFSM[Part_Transcript[Next2_Part[i]]].second.size() > Thiscount) {
                                                Thiscount = AlreadyFSM[Part_Transcript[Next2_Part[i]]].second.size();
                                                positions = Next2_Part[i];
                                            }
                                        }
                                    }                                
                                    MaxTranscriptName = Part_Transcript[positions];                              
                                    if (MaxTranscriptName.size() != 0) {
                                        SingleReads_inTranscripts[MaxTranscriptName].push_back(read_name);
                                    }
                                }
                            }               
                        }
                    }
                }

            }
        }
    }

    std::unordered_map<int, std::unordered_map<std::string, std::vector<std::string>>> FileSingleGeneNumber;
    FileSingleGeneNumber.reserve(FileNo);
    for (int i = 0; i < FileNo; i++) {
        FileSingleGeneNumber.emplace(i, std::unordered_map<std::string, std::vector<std::string>>{});
    }

    if (FileNo > 1) {
        std::map<int, std::vector<std::string>> thisGene_FileCounts;
        if (SingleReads_inGenes.size() > 0) {
            for (const auto& EachGene:SingleReads_inGenes) {
                if (EachGene.second.size() > 0) {
                    thisGene_FileCounts.clear();
                    for (int i = 0; i < FileNo; i++) {
                        thisGene_FileCounts[i] = {};
                    }
                    // 开始数数;
                    for (const auto& EachRead:EachGene.second) {
                        int thisReadFile = std::stoi(groupreadfiles[EachRead]);
                        thisGene_FileCounts[thisReadFile].push_back(EachRead);
                    }
                    // 最终数量结果写入;
                    for (const auto& EachFile:thisGene_FileCounts) {
                        FileSingleGeneNumber[EachFile.first][EachGene.first] = EachFile.second; 
                    }
                }
            }

            std::string first_part, second_part;
            for (const auto& EachTx:SingleReads_inTranscripts) {
                std::string buffer;
                buffer.reserve(EachTx.second.size() * 80);  
                size_t pos = EachTx.first.find('|');
                if (pos != std::string::npos) {
                    first_part = EachTx.first.substr(0, pos);
                    second_part = EachTx.first.substr(pos+1);
                }
                for (const auto& eachR:EachTx.second) {
                    buffer += eachR;
                    buffer += '\t';
                    buffer += "SE\t";
                    buffer += second_part;
                    buffer += '\t';
                    buffer += first_part;
                    buffer += '\t';
                    buffer += groupreadfiles[eachR];
                    buffer += '\n';
                }
                if (Trace.is_open()) {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << buffer;
                }
            }            
        }

    } else {
        if (SingleReads_inGenes.size() > 0) {
            for (const auto& EachGene:SingleReads_inGenes) {          
                FileSingleGeneNumber[0][EachGene.first] = EachGene.second;
            }

            std::string first_part, second_part;
            for (const auto& EachTx:SingleReads_inTranscripts) {
                std::string buffer;
                buffer.reserve(EachTx.second.size() * 80);  
                size_t pos = EachTx.first.find('|');
                if (pos != std::string::npos) {
                    first_part = EachTx.first.substr(0, pos);
                    second_part = EachTx.first.substr(pos+1);
                }
                for (const auto& eachR:EachTx.second) {
                    buffer += eachR;
                    buffer += '\t';
                    buffer += "SE\t";
                    buffer += second_part;
                    buffer += '\t';
                    buffer += first_part;
                    buffer += '\t';
                    buffer += "0";
                    buffer += '\n';
                }
                if (Trace.is_open()) {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << buffer;
                }
            }
        }
    }
    SingleReads.file_SE_reads_gene_number = std::move(FileSingleGeneNumber);
    SingleReads.Transcript_with_SE_reads = std::move(SingleReads_inTranscripts);

    return SingleReads;
}


void merge_high_sj(std::unordered_map<std::string, std::vector<std::array<int,2>>>& group_annotations, std::unordered_map<std::string, std::string>& transcript_strand, std::map<std::array<int,2>, int>& Grouphighsjs) {

    if (group_annotations.size() != 0) {
        for (const auto& isoform:group_annotations) {
            for (const auto& sj:isoform.second) {
                auto it = Grouphighsjs.find(sj);
                if (it == Grouphighsjs.end()) {
                    if (transcript_strand[isoform.first] == "+") {
                        Grouphighsjs[sj] = 1;
                    } else {
                        Grouphighsjs[sj] = 2;
                    }
                }
            }
        }
    } else {
        return;
    }
}


std::unordered_map<std::size_t, std::vector<std::string>> classifyReadsVec(std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs) {
    std::unordered_map<std::size_t, std::vector<std::string>> ClustersReads;
    std::map<std::vector<std::array<int,2>>, std::size_t> readsj2sizet;

    std::size_t clusterNumber = 0;
    for (const auto& pair : groupreadsjs) {
        const std::string& vecname = pair.first;
        const std::vector<std::array<int,2>>& vec = pair.second;
        auto itsj = readsj2sizet.find(vec);
        if (itsj != readsj2sizet.end()) {
            std::size_t exist_cluster_id = itsj->second;
            ClustersReads[exist_cluster_id].push_back(vecname);
        } else {
            // 新的cluster;
            readsj2sizet[vec] = clusterNumber;
            ClustersReads[clusterNumber].push_back(vecname);
            clusterNumber++;
        }
    }
    return ClustersReads;
}


std::unordered_map<std::size_t, std::array<int,2>> get_every_cluster_begin_end(std::unordered_map<std::size_t, std::vector<std::string>>& ClusterReads, 
                                                                               std::unordered_map<std::string, std::array<int,2>>& AllReadBE){
    std::unordered_map<std::size_t, std::array<int,2>> CluserBeginEnd;
    std::vector<int> ExonbeginVec;
    std::vector<int> ExonendVec;
    for (const auto& Tuan:ClusterReads){
        ExonbeginVec.clear();
        ExonendVec.clear();
        CluserBeginEnd[Tuan.first] = {};
        for (const auto& oneread:Tuan.second){
            ExonbeginVec.push_back(AllReadBE[oneread][0]);
            ExonendVec.push_back(AllReadBE[oneread][1]);
        }
        //分别对两个变量取值;
        auto min_it = std::min_element(ExonbeginVec.begin(), ExonbeginVec.end());
        auto max_it = std::max_element(ExonendVec.begin(), ExonendVec.end());
        CluserBeginEnd[Tuan.first][0] = *min_it;
        CluserBeginEnd[Tuan.first][1] = *max_it;
    }
    return CluserBeginEnd;
}


struct ReferenceCluster
{
    std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>> FSM;
    std::map<std::size_t, std::vector<std::string>> ISM; //保存FSM的ISM和Others的ISM;
    std::unordered_map<std::size_t, std::vector<std::string>> Others;
};

ReferenceCluster get_FSM_and_others_sc(std::unordered_map<std::size_t, std::vector<std::string>>& ClustersReads, 
                                   std::unordered_map<std::string, std::vector<std::array<int,2>>>& GroupAnno, 
                                   std::unordered_map<std::string, std::array<int,2>>& GroupAnnoBE, 
                                   std::unordered_map<std::string, std::array<int,2>>& AllReadBE, 
                                   std::unordered_map<std::string, std::vector<std::array<int,2>>>& AllSJs,
                                   std::unordered_map<std::string, std::string>& GroupReadFile, 
                                   std::unordered_map<std::string, std::string>& GroupReadBarcode,
                                   std::ofstream& Trace) {
    ReferenceCluster FSMISMO = {};

    std::vector<std::array<int,2>> each_read_SJs; 
    std::vector<std::string> statis_FSM; 
    std::vector<std::string> statis_ISM; 
    std::vector<int> ISM_flag; ISM_flag.reserve(30);

    std::string FSM_name; 
    std::vector<int> minFSM; minFSM.reserve(30); 

    std::string first_part, second_part;

    for (auto it = ClustersReads.begin(); it != ClustersReads.end(); ++it) {
        const auto& read_list = it->second;
        const auto& first_read = read_list[0];
        const auto& each_read_SJs_ref = AllSJs.at(first_read);
        each_read_SJs = each_read_SJs_ref;
        
        statis_FSM.clear(); 
        statis_ISM.clear();
        for (const auto& pair : GroupAnno) {
            const auto& ref_sj = pair.second;

            if (ref_sj.size() == each_read_SJs.size() && ref_sj == each_read_SJs) {
                statis_FSM.push_back(pair.first);
                continue;
            } else {
                if (ref_sj.size() > each_read_SJs.size()) {
                    ISM_flag.clear();
                    ISM_flag.reserve(each_read_SJs.size());
                    for (const auto& element : each_read_SJs) {
                        auto it_pos = std::find(ref_sj.begin(), ref_sj.end(), element);
                        if (it_pos == ref_sj.end()){
                            ISM_flag.push_back(10000);  
                            break; 
                        } else {
                            ISM_flag.push_back(it_pos - ref_sj.begin());
                        }
                    }
                    if (!ISM_flag.empty() &&
                        ISM_flag.front() != 10000 &&
                        ISM_flag.back() != 10000 &&
                        (ISM_flag.back() - ISM_flag.front() == (int)ISM_flag.size() - 1))
                    {
                        statis_ISM.push_back(pair.first);
                    }                    
                }
            }
        }

        if (statis_FSM.size() == 1) {
            FSM_name = statis_FSM[0];
            auto& fsm_entry = FSMISMO.FSM[FSM_name];
            fsm_entry.first = each_read_SJs;
            fsm_entry.second.insert(fsm_entry.second.end(),
                                    read_list.begin(), read_list.end());
            
            size_t pos = FSM_name.find('|');
            if (pos != std::string::npos){
                first_part = FSM_name.substr(0, pos);
                second_part = FSM_name.substr(pos + 1);
            }

            if (Trace.is_open()) {
                const auto& reads = ClustersReads[it->first];
                std::string buffer;
                buffer.reserve(reads.size() * 80);
                for (const auto& EachRead : reads) {
                    buffer += EachRead;
                    buffer += '\t';
                    buffer += "FSM\t";
                    buffer += second_part;
                    buffer += '\t';
                    buffer += first_part;
                    buffer += '\t';
                    buffer += GroupReadBarcode[EachRead];
                    buffer += '\t';
                    buffer += GroupReadFile[EachRead];
                    buffer += '\n';
                }
                {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << buffer;
                }
            } 

        } else if (!statis_FSM.empty()){
            std::string buffer;
            buffer.reserve(read_list.size() * 80);

            for (const auto& refx : statis_FSM) {
                auto& entry = FSMISMO.FSM[refx];
                entry.first = each_read_SJs;
                entry.second.clear();
            }

            for (const auto& readx : read_list) { 
                minFSM.clear();
                minFSM.reserve(statis_FSM.size());

                for (const auto& refx : statis_FSM){
                    minFSM.push_back(
                        std::abs(AllReadBE[readx][0] - GroupAnnoBE[refx][0]) +
                        std::abs(AllReadBE[readx][1] - GroupAnnoBE[refx][1])
                    );
                }

                int minmin = std::min_element(minFSM.begin(), minFSM.end()) - minFSM.begin();

                FSM_name = statis_FSM[minmin];
                FSMISMO.FSM[FSM_name].second.push_back(readx); 

                size_t pos = FSM_name.find('|');
                if (pos != std::string::npos){
                    first_part = FSM_name.substr(0, pos);
                    second_part = FSM_name.substr(pos + 1);
                }
                buffer += readx;
                buffer += '\t';
                buffer += "FSM\t";
                buffer += second_part;
                buffer += '\t';
                buffer += first_part;
                buffer += '\t';
                buffer += GroupReadBarcode[readx];
                buffer += '\t';
                buffer += GroupReadFile[readx];
                buffer += '\n';                
            } 

            if (Trace.is_open()) {
                std::unique_lock<std::mutex> lock(traceMutex);
                Trace << buffer;
            }                  

        } else { 

            if (statis_ISM.empty()){
                FSMISMO.Others[it->first] = read_list;
            } else {

                FSMISMO.ISM[it->first] = read_list;
                for (const auto& name : statis_ISM) {
                    auto itFSM = FSMISMO.FSM.find(name);
                    if (itFSM == FSMISMO.FSM.end()) {
                        auto& entry = FSMISMO.FSM[name];
                        entry.first = GroupAnno[name];
                        entry.second.clear();
                    }
                }
            }
        }
    }
    return FSMISMO;
}



struct HighLowClusters{
    std::map<size_t, std::vector<std::string>> HighConClusters;
    std::unordered_map<size_t, std::vector<std::string>> LowConClusters;
    std::unordered_map<size_t, std::string> HighStrand;
};

HighLowClusters get_HighLow_clusters(std::unordered_map<std::size_t, std::vector<std::string>>& Others_Clusters, 
                                    std::unordered_map<std::string, std::vector<std::array<int,2>>>& AllSJs, 
                                    std::unordered_map<std::string, std::vector<int>>& readsigns,
                                    const int& sjsupportnumbers){
    HighLowClusters HighLowC = {};
    std::map<std::size_t, std::vector<std::string>> onlyHighC;
    std::unordered_map<size_t, std::string> onlyHighStrand;

    std::map<std::size_t, std::vector<std::string>> onlyHighC_copy;
    std::unordered_map<size_t, std::string> onlyHighStrand_copy;

    std::vector<std::string> each_cluster; 
    std::size_t hashname; 
    std::vector<int> eachreadsign;
    
    int SumCount = 0;
    int StrandCount = 0;
    int pos1 = 0;
    int pos2 = 0;
    std::vector<int> results;
    int SjNumber = 0;
    bool StrandFlag = 0;

    for (auto it = Others_Clusters.begin(); it != Others_Clusters.end(); ++it) {
        hashname = it->first; 
        each_cluster = it->second; 
        if (each_cluster.size() > 1){
            SumCount = 0;
            StrandCount = 0;
            SjNumber = readsigns[each_cluster[0]].size();
            results.assign(SjNumber, 0);
            StrandFlag = 0;

            for (const auto& each:each_cluster){
                eachreadsign = readsigns[each];
                for (int i = 0; i < eachreadsign.size(); i++) {
                    if (eachreadsign[i] == 1) {
                        results[i] = results[i] + eachreadsign[i];
                        if (StrandCount == 2) {
                            StrandFlag = 1;
                        }
                        StrandCount = 1;
                    } else if (eachreadsign[i] == 2) {
                        results[i] = results[i] + 1;
                        if (StrandCount == 1) {
                            StrandFlag = 1;
                        }
                        StrandCount = 2;
                    }
                }
            }

            for (const auto& eachBit:results) {
                if (eachBit > sjsupportnumbers) {
                    SumCount = SumCount + 1;
                }
            }

            if (SumCount == SjNumber) {
                if (StrandFlag) {
                    // 链不一致;
                    HighLowC.LowConClusters[hashname] = Others_Clusters[hashname];               
                } else {
                    if (StrandCount == 1) {
                        onlyHighStrand[hashname] = "+";
                        onlyHighC[hashname] = Others_Clusters[hashname];
                    } else {
                        onlyHighStrand[hashname] = "-";
                        onlyHighC[hashname] = Others_Clusters[hashname];
                    }
                }
            } else {
                HighLowC.LowConClusters[hashname] = Others_Clusters[hashname];
                
            }
        } else {
            HighLowC.LowConClusters[hashname] = Others_Clusters[hashname];
        }
    }

    onlyHighC_copy = onlyHighC;
    onlyHighStrand_copy = onlyHighStrand;

    std::vector<std::array<int,2>> SJs1;
    std::vector<std::array<int,2>> SJs2;
    for (const auto& clus1: onlyHighC) {
        SJs1 = AllSJs[clus1.second[0]];
        for (const auto& clus2: onlyHighC){
            SJs2 = AllSJs[clus2.second[0]];
            if (SJs2.size()>SJs1.size()){
                if (std::includes(SJs2.begin(), SJs2.end(), SJs1.begin(), SJs1.end())){

                    pos1 = std::find(SJs2.begin(), SJs2.end(), SJs1[0]) - SJs2.begin();
                    pos2 = std::find(SJs2.begin(), SJs2.end(), SJs1[SJs1.size()-1]) - SJs2.begin();
                    if (pos2-pos1 == SJs1.size()-1){
                        HighLowC.LowConClusters[clus1.first] = clus1.second;
                        onlyHighC_copy.erase(clus1.first);
                        onlyHighStrand_copy.erase(clus1.first);
                    }
                }
            }
        }
    }
    for (const auto& eachHighCluster:onlyHighC_copy){
        HighLowC.HighConClusters[eachHighCluster.first] = eachHighCluster.second;
        HighLowC.HighStrand[eachHighCluster.first] = onlyHighStrand_copy[eachHighCluster.first];
    }
    return HighLowC;
}



void get_filtered_FSM(std::unordered_map<std::size_t, std::vector<std::string>>& LowReads, 
                      std::unordered_map<std::string, std::vector<std::array<int,2>>>& GroupAnno, 
                      std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>>& AllFSM, 
                      std::unordered_map<std::string, std::vector<std::array<int,2>>>& AllSJs, 
                      std::unordered_map<std::string, std::string>& GroupReadFile,
                      std::unordered_map<std::string, std::string>& GroupReadBarcode,
                      std::ofstream& Trace) {

    std::unordered_map<std::string, std::unordered_map<std::size_t, std::vector<std::string>>> RecycleAnnoName_value;
    std::set<std::size_t> DeleteItem;

    int RecycleLength = 0;
    int count = 0;
    std::string minAnnoName;
    int dis = 10000;
    int FL = 0;
    std::unordered_map<std::string, int> potential_FSM;
    std::vector<std::array<int,2>> LowReadSJVec;

    std::string first_part;
    std::string second_part;

    for (const auto& each_Readcluster:LowReads){
        LowReadSJVec = AllSJs[each_Readcluster.second[0]];

        potential_FSM.clear();
        for (const auto& each_anno:GroupAnno){

            if (each_anno.second.size() == LowReadSJVec.size()){
                RecycleLength = IntervalMerge(each_anno.second, LowReadSJVec) - IntervalIntersection(each_anno.second, LowReadSJVec) - findNonOverlappingIntervals(each_anno.second, LowReadSJVec);
                if (RecycleLength < 60 and RecycleLength > 0) {
                    potential_FSM[each_anno.first] = RecycleLength;
                }
            }
        }

        if (potential_FSM.size() == 1) {

            DeleteItem.insert(each_Readcluster.first);

            size_t pos = (*(potential_FSM.begin())).first.find('|');
            if (pos != std::string::npos){
                first_part = (*(potential_FSM.begin())).first.substr(0, pos);
                second_part = (*(potential_FSM.begin())).first.substr(pos + 1);
            }

            if (Trace.is_open()) {
                std::string buffer;
                buffer.reserve(each_Readcluster.second.size() * 80);
                for (const auto& EachRead : each_Readcluster.second) {
                    const auto& readFile = GroupReadFile[EachRead];
                    const auto& readBarcode = GroupReadBarcode[EachRead];
                    buffer += EachRead;
                    buffer += '\t';
                    buffer += "simFSM\t";
                    buffer += second_part;
                    buffer += '\t';
                    buffer += first_part;
                    buffer += '\t';
                    buffer += readBarcode;
                    buffer += '\t';
                    buffer += readFile;
                    buffer += '\n';
                }
                {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << buffer;                    
                }
            }                        
            if (RecycleAnnoName_value.find((*potential_FSM.begin()).first) != RecycleAnnoName_value.end()){
                RecycleAnnoName_value[(*potential_FSM.begin()).first][each_Readcluster.first] = each_Readcluster.second;
            } else{
                RecycleAnnoName_value[(*potential_FSM.begin()).first] = {};
                RecycleAnnoName_value[(*potential_FSM.begin()).first][each_Readcluster.first] = each_Readcluster.second;
            }
        } else if (potential_FSM.size() > 1) {
            DeleteItem.insert(each_Readcluster.first);
            dis = 100000;
            FL = 0;
            for (const auto& Po:potential_FSM){
                if (AllFSM.find(Po.first) != AllFSM.end()) {
                    if (AllFSM[Po.first].second.size() > FL) {
                        FL = AllFSM[Po.first].second.size();
                        minAnnoName = Po.first;
                    }
                }
            }
            if (FL == 0) {
                for (const auto& Po:potential_FSM) {
                    if (Po.second < dis) {
                        dis = Po.second;
                        minAnnoName = Po.first;
                    }
                }
            }
            size_t pos = minAnnoName.find('|');
            if (pos != std::string::npos){
                first_part = minAnnoName.substr(0, pos);
                second_part = minAnnoName.substr(pos + 1);
            }            

            if (Trace.is_open()) {
                std::string buffer;
                buffer.reserve(each_Readcluster.second.size() * 70);
                for (const auto& EachRead : each_Readcluster.second) {
                    const std::string& file = GroupReadFile[EachRead]; 
                    const auto& readBarcode = GroupReadBarcode[EachRead];
                    buffer += EachRead;
                    buffer += '\t';
                    buffer += "simFSM\t";
                    buffer += second_part;
                    buffer += '\t';
                    buffer += first_part;
                    buffer += '\t';
                    buffer += readBarcode;
                    buffer += '\t';
                    buffer += file;
                    buffer += '\n';
                }            
                {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << buffer;
                }
            }   

            if (RecycleAnnoName_value.find(minAnnoName) != RecycleAnnoName_value.end()) {
                RecycleAnnoName_value[minAnnoName][each_Readcluster.first] = each_Readcluster.second;
            } else{
                RecycleAnnoName_value[minAnnoName] = {};
                RecycleAnnoName_value[minAnnoName][each_Readcluster.first] = each_Readcluster.second;
            }
        }
    }   

    std::vector<std::string> newName;

    for (const auto& newFSM:RecycleAnnoName_value) {
        newName.clear();
        for (const auto& allCluster:newFSM.second) { 
            newName.insert(newName.end(), allCluster.second.begin(), allCluster.second.end());
        }

        if (AllFSM.find(newFSM.first) != AllFSM.end()){
            AllFSM[newFSM.first].second.insert(AllFSM[newFSM.first].second.end(), newName.begin(), newName.end());
        } else{
            AllFSM[newFSM.first].first = GroupAnno[newFSM.first];
            AllFSM[newFSM.first].second.insert(AllFSM[newFSM.first].second.end(), newName.begin(), newName.end());
        }
    }

    for (const auto& eachSJ:DeleteItem){
        auto it = LowReads.find(eachSJ);
        if (it != LowReads.end()){
            LowReads.erase(it);
        }
    }
}


struct SpliceChainClass
{
    std::unordered_map<std::size_t, std::array<int,2>> ClusterCoverage;
    std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>> FSM;
    std::map<std::size_t, std::vector<std::string>> ISM;
    std::map<size_t, std::vector<std::string>> HighConClusters;
    std::unordered_map<size_t, std::vector<std::string>> LowConClusters;
    std::unordered_map<size_t, std::string> HighStrand;
};

SpliceChainClass generate_splice_chain_class( 
                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs, 
                                std::unordered_map<std::string, std::array<int,2>>& groupreadcoverage, 
                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupannotations, 
                                std::unordered_map<std::string, std::vector<int>> groupreadsigns,
                                std::unordered_map<std::string, std::array<int,2>>& AnnoCoverage,
                                std::unordered_map<std::string, std::string>& groupreadfiles,
                                std::unordered_map<std::string, std::string>& groupreadbarcodes,
                                std::unordered_map<std::string, std::string>& groupreadumis,
                                std::ofstream& traceFilePath,
                                const int& Sj_Support_Number) {
    SpliceChainClass SCC;
    ReferenceCluster FsmIsmOthers;
    HighLowClusters Others2HighLow;
    std::unordered_map<std::size_t, std::vector<std::string>> groupCluster = classifyReadsVec(groupreadsjs);
    
    // std::unordered_map<std::size_t, std::vector<std::string>> New_groupCluster = delete_umi_replicate(groupCluster, groupreadbarcodes, groupreadumis);
    delete_umi_replicate(groupCluster, groupreadbarcodes, groupreadumis);    

    SCC.ClusterCoverage = get_every_cluster_begin_end(groupCluster, groupreadcoverage);
    FsmIsmOthers = get_FSM_and_others_sc(groupCluster, groupannotations, AnnoCoverage, groupreadcoverage, groupreadsjs, groupreadfiles, groupreadbarcodes, traceFilePath);
    Others2HighLow = get_HighLow_clusters(FsmIsmOthers.Others, groupreadsjs, groupreadsigns, Sj_Support_Number);

    get_filtered_FSM(Others2HighLow.LowConClusters, groupannotations, FsmIsmOthers.FSM, groupreadsjs, groupreadfiles, groupreadbarcodes, traceFilePath);

    SCC.FSM = FsmIsmOthers.FSM;
    SCC.ISM = FsmIsmOthers.ISM;
    SCC.HighConClusters = Others2HighLow.HighConClusters;
    SCC.HighStrand = Others2HighLow.HighStrand;
    SCC.LowConClusters = Others2HighLow.LowConClusters;

    return SCC;
}


struct DistanceInform
{
    std::unordered_map<int, std::set<int>> NodeDistanceSet;
    std::unordered_map<int, std::vector<std::array<int,2>>> Index2Unknown;
    std::unordered_map<int, std::string> Index2Anno;
    std::unordered_map<int, int> Index2Count;
    std::unordered_map<int, std::size_t> Index2hashname;
    std::unordered_map<int, std::string> Index2novelname; 
};


DistanceInform get_distance_matrix(std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupannotations, 
                                   std::map<std::size_t, std::vector<std::string>>& HighClusters, 
                                   std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs, int DistanceG, 
                                   std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>>& AllFSM){
    
    DistanceInform Graph_Init_inform;
    
    for (int i = 0; i < groupannotations.size()+HighClusters.size(); ++i){
        Graph_Init_inform.NodeDistanceSet[i] = {};
    }

    std::vector<std::array<int,2>> SJVec1;
    std::vector<std::array<int,2>> SJVec2;
    int LineElement = 0; //距离;
    int NameNumber = 0;
    int HighNumber = groupannotations.size();

    if (groupannotations.size() != 0){
        for (const auto& EachAnno : groupannotations){
            Graph_Init_inform.Index2Anno[NameNumber] = EachAnno.first;

            auto it = AllFSM.find(EachAnno.first);

            if (it != AllFSM.end()) {
                Graph_Init_inform.Index2Count[NameNumber] = (it->second).second.size();
            } else{
                Graph_Init_inform.Index2Count[NameNumber] = 0;
            }

            if (HighClusters.size() != 0){
                HighNumber = groupannotations.size();
                for (const auto& eachHC : HighClusters){
                    SJVec1 = groupreadsjs[eachHC.second[0]];
                    LineElement = IntervalMerge(EachAnno.second, SJVec1) - IntervalIntersection(EachAnno.second, SJVec1) - findNonOverlappingIntervals(EachAnno.second, SJVec1);

                    if ((LineElement < DistanceG+1) && (LineElement > 0)){
                        Graph_Init_inform.NodeDistanceSet[NameNumber].insert(HighNumber);
                        Graph_Init_inform.NodeDistanceSet[HighNumber].insert(NameNumber);
                    }
                    HighNumber++;
                }
            }
            NameNumber++;
        }
    }

    if (HighClusters.size() != 0) {
        for (const auto& HC1 : HighClusters){
            SJVec1 = groupreadsjs[HC1.second[0]];
            Graph_Init_inform.Index2Unknown[NameNumber] = SJVec1;
            Graph_Init_inform.Index2Count[NameNumber] = HC1.second.size();
            Graph_Init_inform.Index2hashname[NameNumber] = HC1.first;
            HighNumber = groupannotations.size();
            for (const auto& HC2 : HighClusters){
                SJVec2 = groupreadsjs[HC2.second[0]];
                LineElement = IntervalMerge(SJVec1, SJVec2) - IntervalIntersection(SJVec1, SJVec2) - findNonOverlappingIntervals(SJVec1, SJVec2);
                if ((LineElement < DistanceG+1) && (LineElement > 0)){
                    Graph_Init_inform.NodeDistanceSet[NameNumber].insert(HighNumber);
                }
                HighNumber++;
            }
            NameNumber++;
        }
    }
    return Graph_Init_inform;
}


int find_max_neighborhood_node(std::set<int>& fireNode, std::unordered_map<int, std::set<int>>& evidence){
    int node_max_nei = 0;
    int max_number = 0;
    for (const auto& node:fireNode){
        if (evidence[node].size() > max_number){
            node_max_nei = node;
            max_number = evidence[node].size();
        }
    }
    return node_max_nei;
}


std::vector<std::vector<int>> Find_Maximal_Cliques(std::unordered_map<int, std::set<int>>& adj){
    std::vector<std::vector<int>> CliqueVec;
    if (adj.empty()) return CliqueVec;
    //初始化Q;
    std::vector<int> Q = {};
    std::set<int> cand;
    std::set<int> subg;
    std::set<int> ext_u;
    std::vector<std::tuple<std::set<int>, std::set<int>, std::set<int>>> stack;

    int q = 0;
    int max_node = 0;
    std::set<int> adj_q;
    std::set<int> subg_q;
    std::set<int> cand_q;

    //adj已经有了;
    for (const auto& eachnode:adj) cand.insert(eachnode.first);
    subg = cand;
    Q.push_back(-1);    
    max_node = find_max_neighborhood_node(subg, adj);
    std::set_difference(cand.begin(), cand.end(),
                        adj[max_node].begin(), adj[max_node].end(),
                        std::inserter(ext_u, ext_u.begin()));

    while (true) {
        if (!ext_u.empty()) {
            q = *ext_u.begin(); 
            ext_u.erase(q);
            cand.erase(q);
            Q.back() = q;
            adj_q = adj[q];
            subg_q.clear();
            std::set_intersection(subg.begin(), subg.end(),
                                adj_q.begin(), adj_q.end(),
                                std::inserter(subg_q, subg_q.begin())); //subg_q
            if (subg_q.empty())
            {
                CliqueVec.push_back(Q);
            } 
            else
            {   
                cand_q.clear();
                std::set_intersection(cand.begin(), cand.end(),
                                adj_q.begin(), adj_q.end(),
                                std::inserter(cand_q, cand_q.begin())); //cand_q
            
                if (!cand_q.empty())
                {
                    stack.push_back(std::make_tuple(subg, cand, ext_u));
                    Q.push_back(-1); //添加一个占位符;
                    subg = subg_q;
                    cand = cand_q;
                    max_node = find_max_neighborhood_node(subg, adj);
                    ext_u.clear();
                    std::set_difference(cand.begin(), cand.end(),
                        adj[max_node].begin(), adj[max_node].end(),
                        std::inserter(ext_u, ext_u.begin()));
                } 
            }
        }
        else
        {   
            if (stack.size() != 0 && Q.size() != 0)
            {
                Q.pop_back();
                std::tie(subg, cand, ext_u) = stack.back();
                stack.pop_back();
            } 
            else
            {
                break;
            }
        }
    }
    return CliqueVec;
}


bool compareBySize(const std::vector<int>& a, const std::vector<int>& b) {
    
    if (a.empty() && b.empty()) {
        return false; 
    }
    if (a.empty()) {
        return true; 
    }
    if (b.empty()) {
        return false; 
    }
    auto AminElement = std::min_element(a.begin(), a.end());
    auto BminElement = std::min_element(b.begin(), b.end());
    // 比较最小元素
    if (*AminElement != *BminElement) {
        return *AminElement < *BminElement; 
    }
    return a.size() < b.size();
}



struct DetectionResults
{
    std::set<int> NodesKnownTrue;
    std::set<int> TrueNodeSet;
    std::set<int> FalseNodeSet;
};


DetectionResults Transcript_Detection(std::vector<std::vector<int>>& CqVec, 
                                      std::unordered_map<int, std::string>& Number2Anno, 
                                      std::unordered_map<int, std::vector<std::array<int,2>>>& Number2Unknown, 
                                      std::unordered_map<int, int>& Number2Count) {

    DetectionResults Node_Results = {};
    if (CqVec.empty()) return Node_Results;

    //函数用到的中间变量;
    std::set<int> subgraph_known_true_nodes;
    std::set<int> subgraph_remain_nodes;
    int standard_read_number = std::numeric_limits<int>::max();
    int known_true_number = 0;
    int max_node_value = 0;

    for (const auto& Anno : Number2Anno){
        Node_Results.NodesKnownTrue.insert(Anno.first);
    }

    if (CqVec.size() > 1) {
        std::sort(CqVec.begin(), CqVec.end(), compareBySize);
    }
    //遍历每个团;
    for (const auto& eachClique : CqVec){
        standard_read_number = std::numeric_limits<int>::max();

        std::set<int> CliqueSet(eachClique.begin(), eachClique.end());
        subgraph_known_true_nodes.clear();
        std::set_intersection(CliqueSet.begin(), CliqueSet.end(), Node_Results.NodesKnownTrue.begin(), Node_Results.NodesKnownTrue.end(),
                        std::inserter(subgraph_known_true_nodes, subgraph_known_true_nodes.begin()));
        subgraph_remain_nodes.clear();
        std::set_difference(CliqueSet.begin(), CliqueSet.end(), subgraph_known_true_nodes.begin(), subgraph_known_true_nodes.end(),
                    std::inserter(subgraph_remain_nodes, subgraph_remain_nodes.begin()));

        if (subgraph_known_true_nodes.size() != 0){

            if (subgraph_known_true_nodes.size() == 1){
                standard_read_number = Number2Count[(*subgraph_known_true_nodes.begin())];
            } else {

                for(const auto& each_True_node : subgraph_known_true_nodes){
                    if (Number2Count[each_True_node] < standard_read_number) {
                        standard_read_number = Number2Count[each_True_node];
                    }
                }
            }

            if (subgraph_remain_nodes.size() != 0){ 

                for (const auto& eachknownnode:subgraph_known_true_nodes){
                    Node_Results.TrueNodeSet.insert(subgraph_known_true_nodes.begin(), subgraph_known_true_nodes.end());
                }

                for (const auto& each_remain_node : subgraph_remain_nodes){
                    
                    if (standard_read_number == 0){
                        Node_Results.FalseNodeSet.insert(each_remain_node);
                    } else {
                        if (Number2Count[each_remain_node] < standard_read_number){
                            Node_Results.FalseNodeSet.insert(each_remain_node);
                        } else {
                            Node_Results.TrueNodeSet.insert(each_remain_node);
                        }
                    }
                }
            } else { 

                for (const auto& eachknownnode:subgraph_known_true_nodes){
                    if (Number2Count[eachknownnode] > 0){
                        Node_Results.TrueNodeSet.insert(eachknownnode);
                    }
                } 
            }

        } else {
            if (subgraph_remain_nodes.size() == 1){
                if (Number2Count[*subgraph_remain_nodes.begin()] > 0) {
                    Node_Results.TrueNodeSet.insert(*subgraph_remain_nodes.begin());
                }
            } else{
                max_node_value = 0;
                for (const auto& each_remain_node : subgraph_remain_nodes){
                    if (max_node_value < Number2Count[each_remain_node]){
                        max_node_value = Number2Count[each_remain_node];
                    }
                }
                for (const auto& each_remain_node : subgraph_remain_nodes){
                    if (max_node_value == Number2Count[each_remain_node]){
                        Node_Results.TrueNodeSet.insert(each_remain_node);
                    }
                    else
                    {
                        Node_Results.FalseNodeSet.insert(each_remain_node);
                    }
                }   
            }
        }
    }
    for (const auto& eachFalse:Node_Results.FalseNodeSet){
        auto iter = Node_Results.TrueNodeSet.find(eachFalse);
        if (iter != Node_Results.TrueNodeSet.end()) {
            Node_Results.TrueNodeSet.erase(eachFalse);
        }
    }
    
    return Node_Results;
}


struct Solvent_sc {
    std::map<std::string, std::vector<std::string>> File_FSM;
    std::map<std::size_t, std::vector<std::string>> File_ISM;
    std::map<std::size_t, std::vector<std::string>> File_HighConClusters;
    std::map<std::size_t, std::vector<std::string>> File_LowConClusters;
};

Solvent_sc get_Solvent_FsmIsmHigh(SpliceChainClass& FsmIsmHigh, int& FileNumber, GroupInformation& groupinformations) {
    Solvent_sc FileSpliceChains;

    std::string thisFileString = std::to_string(FileNumber);
    
    if (FsmIsmHigh.FSM.size() > 0) {
        for (const auto& eachCluster:FsmIsmHigh.FSM) {
            std::string AnnoName = eachCluster.first;
            FileSpliceChains.File_FSM[AnnoName] = {};
            for (const auto& eachRead:eachCluster.second.second) {
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                if (read_file == thisFileString) {
                    FileSpliceChains.File_FSM[AnnoName].push_back(eachRead);
                }
            }        
        }
    } 

    if (FsmIsmHigh.ISM.size() > 0) {
        for (const auto& eachCluster:FsmIsmHigh.ISM) {
            std::size_t NonName = eachCluster.first;
            FileSpliceChains.File_ISM[NonName] = {};
            for (const auto& eachRead:eachCluster.second) {
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                if (read_file == thisFileString) {
                    FileSpliceChains.File_ISM[NonName].push_back(eachRead);
                }
            }
        }
    }

    if (FsmIsmHigh.HighConClusters.size() > 0) {
        for (const auto& eachCluster:FsmIsmHigh.HighConClusters) {
            std::size_t NonName = eachCluster.first;
            FileSpliceChains.File_HighConClusters[NonName] = {};
            for (const auto& eachRead:eachCluster.second) {
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                if (read_file == thisFileString) {
                    FileSpliceChains.File_HighConClusters[NonName].push_back(eachRead);
                }
            }      
        }
    }

    if (FsmIsmHigh.LowConClusters.size() > 0) {
        for (const auto& eachCluster:FsmIsmHigh.LowConClusters) {
            std::size_t NonName = eachCluster.first;
            FileSpliceChains.File_LowConClusters[NonName] = {};
            for (const auto& eachRead:eachCluster.second) {
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                if (read_file == thisFileString) {
                    FileSpliceChains.File_LowConClusters[NonName].push_back(eachRead);
                }
            }
        }
    }

    return FileSpliceChains;
}


std::string get_co_sj_gene(std::vector<std::array<int,2>>& novel_isoform, 
                            std::unordered_map<std::string, std::array<int,2>>& gene_set,
                            std::unordered_map<std::string, std::vector<std::string>>& gene2tx,
                            std::unordered_map<std::string, std::vector<std::array<int,2>>>& tx2sj) {
    //
    std::vector<std::string> TempGene;
    std::vector<int> TempDist;
    std::set<std::array<int,2>> SJ_set; 
    std::vector<std::string> ALL_tx;
    std::set<std::array<int, 2>> novel_isoform_set(novel_isoform.begin(), novel_isoform.end());
    std::set<std::array<int, 2>> intersection;

    // 对每个基因循环;
    for (const auto& A_Gene:gene_set) {
        SJ_set.clear();
        intersection.clear();
        TempGene.push_back(A_Gene.first);
        ALL_tx = gene2tx[A_Gene.first];
        for (const auto& A_Tx:ALL_tx) {
            // 这个转录本在里面;
            if (tx2sj.find(A_Tx) != tx2sj.end()) {
                SJ_set.insert(tx2sj[A_Tx].begin(), tx2sj[A_Tx].end());
            }
        }
        std::set_intersection(novel_isoform_set.begin(), novel_isoform_set.end(), SJ_set.begin(), SJ_set.end(),
                          std::inserter(intersection, intersection.begin()));
        TempDist.push_back(intersection.size());
    }
    // 判断哪个基因多;
    auto max_it = std::max_element(TempDist.begin(), TempDist.end());
    std::string MaxGeneName;
    if (max_it != TempDist.end()) {
        if (*max_it > 0) {
            MaxGeneName = TempGene[std::distance(TempDist.begin(), max_it)];
        } else {
            MaxGeneName = "";
        }
    } else {
        MaxGeneName = "";
    }    
    return MaxGeneName;
}


struct OutputInformation
{
    std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, double>> Transcript_Annotations;
    std::unordered_map<std::string, std::string> transcript2gene;
};

OutputInformation Write_Detection_Transcript2gtf_AllFiles(std::ofstream& Updated_Files, std::ofstream& Trace,
                                                DetectionResults& noderesults, DistanceInform& Disinform, 
                                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupannotations, 
                                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& Annoexon, 
                                                std::unordered_map<std::size_t, std::array<int, 2>>& groupreadcoverage, 
                                                std::unordered_map<std::string, std::array<int,2>>& groupallgene, 
                                                std::unordered_map<std::string, std::array<int,2>>& Annogenecovergae, 
                                                std::unordered_map<std::string, std::string>& GTF_Transcript_Strand, 
                                                std::unordered_map<size_t, std::string>& High_Strand, 
                                                std::string& chrname, std::string& group_size, 
                                                std::map<size_t, std::vector<std::string>>& HighClusters,
                                                std::unordered_map<std::string, std::string>& groupreadfiles,
                                                std::unordered_map<std::string, std::string>& groupreadbarcodes,
                                                std::map<std::string, std::set<std::string>>& allfilebarcodeset,
                                                std::unordered_map<std::string, std::vector<std::string>>& gtf_gene2tx,
                                                std::unordered_map<std::string, std::string>& gtf_gene_std) {
    OutputInformation FinalAnnotations;
    if (noderesults.TrueNodeSet.empty()) return FinalAnnotations;
    std::string itsname;
    std::string itsbarcode;
    std::vector<std::array<int,2>> itssj;
    std::vector<std::array<int,2>> itsexon;
    std::vector<std::string> TempGene;
    std::string TempGeneName;
    std::vector<int> TempDist;
    int novel_count = 0;
    std::size_t nameya;
    std::array<int,2> BE;
    std::string first_part;
    std::string second_part;

    if (Updated_Files.is_open()) {

        for (const auto aaa:noderesults.TrueNodeSet) {

            if (aaa < Disinform.Index2Anno.size()) {

                itsname = Disinform.Index2Anno[aaa];
                itssj = groupannotations[itsname];
                itsexon = Annoexon[itsname];
                FinalAnnotations.Transcript_Annotations[itsname].first = itssj;
                FinalAnnotations.Transcript_Annotations[itsname].second = Disinform.Index2Count[aaa];

                size_t pos = itsname.find('|');
                if (pos != std::string::npos) {
                    first_part = itsname.substr(0, pos);
                    second_part = itsname.substr(pos + 1);
                }
                {
                    if (!itsexon.empty()) {
                        std::unique_lock<std::mutex> lock(updatedGtfMutex);
                        Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "transcript" << '\t' << itsexon[0][0] << '\t' << itsexon[itsexon.size()-1][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[itsname] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\";" << '\n';
                        for (int i = 0; i < itsexon.size(); i++) {
                            Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "exon" << '\t' << itsexon[i][0] << '\t' << itsexon[i][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[itsname] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\"; exon_number \"" << i << "\";" << '\n';
                        }
                    }
                }
                FinalAnnotations.transcript2gene[itsname] = first_part;

            } else {
                novel_count = novel_count + 1;
                itsname = chrname + "-novel-" + group_size + "-" + std::to_string(aaa) + "-" + std::to_string(novel_count);
                Disinform.Index2novelname[aaa] = itsname;
                itssj = Disinform.Index2Unknown[aaa];
                FinalAnnotations.Transcript_Annotations[itsname].first = itssj;
                FinalAnnotations.Transcript_Annotations[itsname].second = Disinform.Index2Count[aaa];
                nameya = Disinform.Index2hashname[aaa];
                BE = groupreadcoverage[nameya];
                
                TempGene.clear();
                if (groupallgene.size() != 0) {
                    TempGeneName = get_co_sj_gene(itssj, groupallgene, gtf_gene2tx, groupannotations);

                    if (TempGeneName.size() != 0) {
                        if (BE[0] > Annogenecovergae[TempGeneName][0] and BE[1] < Annogenecovergae[TempGeneName][1]) {
                            first_part = TempGeneName;
                        } else {
                            first_part = "NA";
                        }
    
                    } else {
                        for (const auto& every_gene:groupallgene) {
                            if (BE[0] > Annogenecovergae[every_gene.first][0] and BE[1] < Annogenecovergae[every_gene.first][1]){
                                TempGene.push_back(every_gene.first);
                            }
                        }
                        if (TempGene.size() == 0){
                            first_part = "NA";
                        } else if (TempGene.size() == 1) {
                            first_part = TempGene[0];
                        } else {
                            TempDist.clear();
                            for (const auto& each_gene:TempGene){
                                int a = abs(BE[1] - Annogenecovergae[each_gene][1]);
                                int b = abs(Annogenecovergae[each_gene][0] - BE[0]);
                                TempDist.push_back((a>b)?a:b);
                            }
                            auto min_it = std::min_element(TempDist.begin(), TempDist.end());
                            first_part = TempGene[std::distance(TempDist.begin(), min_it)];  
                        }
                    }

                } else {
                    first_part = "NA";
                }

                FinalAnnotations.transcript2gene[itsname] = first_part;
                {   
                    if (first_part != "NA") {
                        second_part = gtf_gene_std[first_part];
                    } else {
                        second_part = High_Strand[nameya];
                    }
                    if (!itssj.empty()) {
                        std::unique_lock<std::mutex> lock(updatedGtfMutex);
                        Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "transcript" << '\t' << BE[0] << '\t' << BE[1] << '\t' << "." << '\t' << second_part << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\";" << '\n';
                        for (int i = 0; i < itssj.size()+1; i++){
                            if (i == 0) {
                                Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << BE[0] << '\t' << itssj[i][0]-1 << '\t' << "." << '\t' << second_part << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" <<  '\n';
                            } else if (i == itssj.size()) {
                                Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << itssj[i-1][1]+1 << '\t' << BE[1] << '\t' << "." << '\t' << second_part << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" << '\n';
                            } else {
                                Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << itssj[i-1][1]+1 << '\t' << itssj[i][0]-1 << '\t' << "." << '\t' << second_part << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" << '\n';
                            }
                        }
                    }
                }

                if (Trace.is_open()) {
                    const auto& read_list = HighClusters[nameya];
                    std::string buffer;
                    buffer.reserve(read_list.size() * 80);  // 预分配，加速
                    for (const auto& EachRead : read_list) {
                        buffer += EachRead;
                        buffer += '\t';
                        buffer += "novel\t";
                        buffer += itsname;
                        buffer += '\t';
                        buffer += first_part;
                        buffer += '\t';
                        buffer += groupreadbarcodes[EachRead];
                        buffer += '\t';
                        buffer += groupreadfiles[EachRead];
                        buffer += '\n';
                    }
                    {
                        std::unique_lock<std::mutex> lock(traceMutex);
                        Trace << buffer;
                    }
                }
            }
        }
    }
    return FinalAnnotations;
}


void get_transcript_init (
    std::unordered_map<std::string, std::vector<double>>& File_TranscriptNumber,
    DetectionResults& noderesults, DistanceInform& Disinform, std::string& chrname, std::string& group_size, 
    int& FileNumber, std::unordered_map<std::string, std::string>& groupreadfiles,
    std::unordered_map<std::string, std::string>& groupreadbarcodes,
    Solvent_sc& filesolvent, const std::set<std::string>& itBarSet, std::unordered_map<std::string,int>& transcript2ID) {
    
    std::string fileIndex = std::to_string(FileNumber);
    int novel_count = 0;

    for (const auto aaa:noderesults.TrueNodeSet) {
        std::string itsname; 
        const std::vector<std::string>* reads = nullptr;

        if (aaa < Disinform.Index2Anno.size()) {
            itsname = Disinform.Index2Anno[aaa];
            auto it = filesolvent.File_FSM.find(itsname);
            if (it != filesolvent.File_FSM.end()) reads = &it->second;
        } else {
            ++novel_count;
            itsname = chrname + "-novel-" + group_size + "-" + std::to_string(aaa) + "-" + std::to_string(novel_count);
            Disinform.Index2novelname[aaa] = itsname;

            auto nameya = Disinform.Index2hashname[aaa];
            auto it = filesolvent.File_HighConClusters.find(nameya);
            if (it != filesolvent.File_HighConClusters.end()) reads = &it->second;
        }

        if (!reads || reads->empty()) continue;

        for (const auto& eachRead : *reads) {
            auto itBC = groupreadbarcodes.find(eachRead);
            if (itBC == groupreadbarcodes.end()) continue;

            const std::string& bc = itBC->second;
            auto itOut = File_TranscriptNumber.find(bc);
            if (itOut == File_TranscriptNumber.end()) continue;

            auto itIndex = transcript2ID.find(itsname);
            if (itIndex == transcript2ID.end()) continue;
            itOut->second[itIndex->second] += 1.0;
        }        
    }
}


struct Barcode_sc {
    std::unordered_map<std::string, std::map<std::string, int>> Barcode_FSM;
    std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>> Barcode_ISM;
    std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>> Barcode_HighConClusters;
    std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>> Barcode_LowConClusters;
};

Barcode_sc get_Barcode_FsmIsmHigh(std::map<std::string, std::vector<std::string>>& thisfilefsm,
                            std::map<std::size_t, std::vector<std::string>>& thisfileism,
                            std::map<std::size_t, std::vector<std::string>>& thisfileHC,
                            std::map<std::size_t, std::vector<std::string>>& thisfileLC,
                            std::unordered_map<std::string, std::string>& groupreadbarcodes) {
    // 输出;
    Barcode_sc barcodesc;
    std::unordered_map<std::string, std::map<std::string, int>> barcodefsm;
    std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>> barcodeism;
    std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>> barcodehc;
    std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>> barcodelc;
    
    // FSM
    for (std::map<std::string, std::vector<std::string>>::const_iterator it = thisfilefsm.begin(); it != thisfilefsm.end(); ++it) {
        const std::string& FSMname = it->first;
        const std::vector<std::string>& reads = it->second;
        for (std::vector<std::string>::const_iterator r = reads.begin(); r != reads.end(); ++r) {
            std::unordered_map<std::string, std::string>::const_iterator itb = groupreadbarcodes.find(*r);
            if (itb == groupreadbarcodes.end()) continue;
            const std::string& barcode = itb->second;
            ++barcodefsm[barcode][FSMname];
        }
    }

    // ISM
    for (std::map<std::size_t, std::vector<std::string>>::const_iterator it = thisfileism.begin(); it != thisfileism.end(); ++it) {
        std::size_t clusterId = it->first;
        const std::vector<std::string>& reads = it->second;
        for (std::vector<std::string>::const_iterator r = reads.begin(); r != reads.end(); ++r) {
            std::unordered_map<std::string, std::string>::const_iterator itb = groupreadbarcodes.find(*r);
            if (itb == groupreadbarcodes.end()) continue;
            const std::string& barcode = itb->second;
            barcodeism[barcode][clusterId].push_back(*r);
        }
    }

    // HC
    for (std::map<std::size_t, std::vector<std::string>>::const_iterator it = thisfileHC.begin(); it != thisfileHC.end(); ++it) {
        std::size_t clusterId = it->first;
        const std::vector<std::string>& reads = it->second;
        for (std::vector<std::string>::const_iterator r = reads.begin(); r != reads.end(); ++r) {
            std::unordered_map<std::string, std::string>::const_iterator itb = groupreadbarcodes.find(*r);
            if (itb == groupreadbarcodes.end()) continue;
            const std::string& barcode = itb->second;
            barcodehc[barcode][clusterId].push_back(*r);
        }
    }

    // LC
    for (std::map<std::size_t, std::vector<std::string>>::const_iterator it = thisfileLC.begin(); it != thisfileLC.end(); ++it) {
        std::size_t clusterId = it->first;
        const std::vector<std::string>& reads = it->second;
        for (std::vector<std::string>::const_iterator r = reads.begin(); r != reads.end(); ++r) {
            std::unordered_map<std::string, std::string>::const_iterator itb = groupreadbarcodes.find(*r);
            if (itb == groupreadbarcodes.end()) continue;
            const std::string& barcode = itb->second;
            barcodelc[barcode][clusterId].push_back(*r);
        }
    }

    barcodesc.Barcode_FSM = barcodefsm;
    barcodesc.Barcode_ISM = barcodeism;
    barcodesc.Barcode_HighConClusters = barcodehc;
    barcodesc.Barcode_LowConClusters = barcodelc;
    return barcodesc;
}




struct FileBarcodeFalseNode
{
    std::set<int> ThisFlaseNode;
    std::unordered_map<int, std::vector<std::string>> FasleNode2Count;
};

void get_File_Barcode_False_node (std::set<int>& FalseNodeSet, 
                            std::unordered_map<int, std::size_t>& Node2hashname, 
                            std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>>& barcodehc, 
                            const std::string& Barcode, FileBarcodeFalseNode& ThisFile_FalseNumber) {
    ThisFile_FalseNumber.ThisFlaseNode.clear();
    ThisFile_FalseNumber.FasleNode2Count.clear();

    std::map<size_t, std::vector<std::string>> this_file_HC;
    auto it = barcodehc.find(Barcode);
    if (it == barcodehc.end()) {
        return;
    } else {
        this_file_HC = it->second;
    }
    for (const auto& node:FalseNodeSet) {
        const auto & node_name = Node2hashname[node];
        auto which = this_file_HC.find(node_name);
        if (which != this_file_HC.end()) {
            ThisFile_FalseNumber.ThisFlaseNode.insert(node);
            ThisFile_FalseNumber.FasleNode2Count[node] = this_file_HC[node_name];
        }
    }
    return;
}


int whether_isoform_part_is_not(std::vector<std::array<int,2>> AnnoSJ, std::vector<std::array<int,2>> ISMSJ){
    int ISM_Flag = 0;
    auto it1 = std::find(AnnoSJ.begin(), AnnoSJ.end(), ISMSJ[0]);
    if (it1 != AnnoSJ.end()){
        int it1_index = it1 - AnnoSJ.begin();
        auto it2 = std::find(AnnoSJ.begin(), AnnoSJ.end(), ISMSJ[ISMSJ.size()-1]);
        if (it2 != AnnoSJ.end()){
            int it2_index = it2 - AnnoSJ.begin();
            if (it2_index - it1_index == ISMSJ.size() - 1){
                ISM_Flag = 1;
            }
        }
    } 
    return ISM_Flag;
}

std::string concatenateSet(const std::set<std::string>& stringSet) {
    if (stringSet.empty()) {
        return "";
    }
    std::ostringstream oss;
    std::string oss_string;
    for (const auto& str : stringSet) {
        oss << str << ',';
    }
    oss_string = oss.str();
    oss_string.pop_back();
    return oss_string;
}

struct IndicateFire
{
    Eigen::MatrixXd Indicate_Matrix;
    Eigen::VectorXd Cluster_Number;
    std::vector<std::string> Order_Transcript_Name_Vector;
};

void Quantification_initialization_Barcodes(std::map<std::size_t, std::vector<std::string>>& groupISM, 
                                            OutputInformation& FinallyAnnotations, 
                                            std::set<int>& truenodeset,
                                            std::set<int>& falsenodeset,
                                            std::unordered_map<int, std::vector<std::string>>& faslenode2count, 
                                            DistanceInform& Disinform, 
                                            std::string& group_size, 
                                            std::map<size_t, std::vector<std::string>>& HighClusters, 
                                            std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs,
                                            std::unordered_map<std::string, std::string>& groupreadfiles,
                                            std::unordered_map<std::string, std::string>& groupreadbarcodes,
                                            std::ofstream& Trace, IndicateFire& OutputResults) {
    OutputResults.Indicate_Matrix.resize(0, 0);
    OutputResults.Cluster_Number.resize(0);

    Eigen::MatrixXd known_ISM_matrix;
    std::vector<std::array<int,2>> AreadSJs;
    
    int RowIndex = -1; int ColIndex = -1;
    int flag = 0;

    const auto& Order_Transcript_Name = OutputResults.Order_Transcript_Name_Vector;
    const int Transcript_Number = Order_Transcript_Name.size();

    Eigen::RowVectorXd EachClusterRowVector(Transcript_Number);
    Eigen::VectorXd known_ISM_cluster_numbers;

    //ISM的追踪文件;
    std::vector<std::string> ism_gtf_name;
    std::set<std::string> first_part_set;
    std::set<std::string> second_part_set;
    std::string first_part;
    std::string second_part;
    std::vector<int> ColIndexVec;
    double max_ratio = 0;
    double this_ratio = 0;

    if (!groupISM.empty()) {
        RowIndex = -1;
        for (const auto& eachISM:groupISM) {
            AreadSJs = groupreadsjs[eachISM.second[0]];
            EachClusterRowVector.setZero();
            ColIndex = -1;
            ism_gtf_name.clear();
            first_part_set.clear();
            second_part_set.clear();  
            max_ratio = 0;
            this_ratio = 0;          

            for (const auto& eachTransctipt:Order_Transcript_Name) {
                ColIndex = ColIndex + 1;
                if (AreadSJs.size() < FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size()){
                    if ((AreadSJs[0][0] >= FinallyAnnotations.Transcript_Annotations[eachTransctipt].first[0][0]) && (AreadSJs[AreadSJs.size()-1][1] <= FinallyAnnotations.Transcript_Annotations[eachTransctipt].first[(FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size()-1)][1])){
                        flag = whether_isoform_part_is_not(FinallyAnnotations.Transcript_Annotations[eachTransctipt].first, AreadSJs);
                        if (flag == 1) {
                            this_ratio = static_cast<double>(AreadSJs.size())/FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size();
                            
                            if (max_ratio < this_ratio) {
                                max_ratio = this_ratio;
                            }                            
                            EachClusterRowVector(ColIndex) = 1;
                            ism_gtf_name.push_back(eachTransctipt);
                            ColIndexVec.push_back(ColIndex);
                        } else {
                            EachClusterRowVector(ColIndex) = 0;
                        }
                    } else{
                        EachClusterRowVector(ColIndex) = 0;
                    }
                } else {
                    EachClusterRowVector(ColIndex) = 0;
                }
            }

            if (ColIndexVec.size() > 1) {
                if (max_ratio > 0.4) {
                    for (const auto& ColIndex:ColIndexVec) {
                        std::string TransName = Order_Transcript_Name[ColIndex];
                        this_ratio = static_cast<double>(AreadSJs.size())/FinallyAnnotations.Transcript_Annotations[TransName].first.size();
                        if (this_ratio < max_ratio) {
                            EachClusterRowVector[ColIndex] = 0;
                            ism_gtf_name.erase(std::remove(ism_gtf_name.begin(), ism_gtf_name.end(), TransName), ism_gtf_name.end());
                        }
                    }
                }
            } 
            if (ism_gtf_name.size() != 0) {

                RowIndex = RowIndex + 1;
                known_ISM_matrix.conservativeResize(RowIndex+1, Transcript_Number);
                known_ISM_matrix.row(RowIndex) = EachClusterRowVector;
                known_ISM_cluster_numbers.conservativeResize(RowIndex+1, 1);
                known_ISM_cluster_numbers(RowIndex, 0) = eachISM.second.size();

                for (const auto& EachTranscript:ism_gtf_name){
                    size_t pos = EachTranscript.find('|');
                    if (pos != std::string::npos) {
                        first_part_set.insert(EachTranscript.substr(0, pos));
                        second_part_set.insert(EachTranscript.substr(pos + 1));
                    } else {
                        second_part_set.insert(EachTranscript);
                        const auto geneName = FinallyAnnotations.transcript2gene[EachTranscript];
                        first_part_set.insert(geneName);                        
                    }
                }

                first_part = concatenateSet(first_part_set);
                second_part = concatenateSet(second_part_set);
                
                if (Trace.is_open()) {
                    const auto& reads = eachISM.second; 
                    std::string buffer;
                    buffer.reserve(reads.size() * 80); 
                    for (const auto& EachRead:reads) {
                        buffer += EachRead;
                        buffer += '\t';
                        buffer += "ISM\t";
                        buffer += second_part;
                        buffer += '\t';
                        buffer += first_part;
                        buffer += '\t';
                        buffer += groupreadbarcodes.at(EachRead);
                        buffer += '\t';
                        buffer += groupreadfiles.at(EachRead);
                        buffer += '\n';
                    }
                    {
                        std::unique_lock<std::mutex> lock(traceMutex);
                        Trace << buffer; 
                    }
                }
            } else {}
        }
    }

    std::set<int> MergeFalseSet(falsenodeset.begin(), falsenodeset.end());
    Eigen::MatrixXd False_novel_candidate_matrix(MergeFalseSet.size(), Transcript_Number);

    std::set<int> FalseNode_NeighborSet;
    std::set<int> FalseNode_TrueSet;
    std::string itsname;
    Eigen::RowVectorXd FalseNode_RowVector(Transcript_Number);
    RowIndex = -1;
    int countC = 0;
    
    if (MergeFalseSet.size() != 0) {
        std::set<int> MergedTrueSet(truenodeset.begin(), truenodeset.end());
        
        for (const auto& FalseNode:MergeFalseSet) {
   
            RowIndex = RowIndex + 1;        
            FalseNode_NeighborSet.clear();
            FalseNode_TrueSet.clear();            
            FalseNode_NeighborSet = Disinform.NodeDistanceSet[FalseNode];

            std::set_intersection(MergedTrueSet.begin(), MergedTrueSet.end(), FalseNode_NeighborSet.begin(), FalseNode_NeighborSet.end(),
            std::inserter(FalseNode_TrueSet, FalseNode_TrueSet.begin()));

            FalseNode_RowVector.setZero();
            ism_gtf_name.clear();

            if (FalseNode_TrueSet.size() != 0) {

                for(const auto& EachTNode:FalseNode_TrueSet){
                    if (EachTNode < Disinform.Index2Anno.size()){
                        itsname = Disinform.Index2Anno[EachTNode];
                        auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                        FalseNode_RowVector[itwhere] = 1;
                        ism_gtf_name.push_back(itsname);
                    } else {
                        itsname = Disinform.Index2novelname[EachTNode];
                        auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                        FalseNode_RowVector[itwhere] = 1;
                        ism_gtf_name.push_back(itsname);
                    }
                }

                False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
            } 
            else
            {   
                countC = 0;
                for(const auto& neinode:FalseNode_NeighborSet) {
                    if (Disinform.NodeDistanceSet[neinode].size() != 0){
                        for (const auto& node:Disinform.NodeDistanceSet[neinode]){
                            if (node < Disinform.Index2Anno.size()) {
                                itsname = Disinform.Index2Anno[node];
                                auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                                FalseNode_RowVector[itwhere] = 1;
                                ism_gtf_name.push_back(itsname);
                                countC++;
                            } else {
                                if (truenodeset.find(node) != truenodeset.end()) {
                                    itsname = Disinform.Index2novelname[node];
                                    auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                                    FalseNode_RowVector[itwhere] = 1;
                                    ism_gtf_name.push_back(itsname);
                                    countC++;
                                }
                            }
                        }
                    }
                }
                if (countC != 0){
                    False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
                } else {

                    std::vector<int> houX = {};
                    int distance = 0;
                    for (const auto& isofrorm:Order_Transcript_Name){
                        distance = IntervalMerge(FinallyAnnotations.Transcript_Annotations[isofrorm].first, Disinform.Index2Unknown[FalseNode]) - IntervalIntersection(FinallyAnnotations.Transcript_Annotations[isofrorm].first, Disinform.Index2Unknown[FalseNode]);
                        houX.push_back(distance);
                    }

                    auto min_it = std::min_element(houX.begin(), houX.end());
                    int min_index = std::distance(houX.begin(), min_it);
                    FalseNode_RowVector[min_index] = 1;
                    False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
                    ism_gtf_name.push_back(Order_Transcript_Name[min_index]);
                }
            }

            first_part_set.clear();
            second_part_set.clear();

            for (const auto& EachTranscript:ism_gtf_name){
                size_t pos = EachTranscript.find('|');
                if (pos != std::string::npos) {
                    first_part_set.insert(EachTranscript.substr(0, pos));
                    second_part_set.insert(EachTranscript.substr(pos + 1));
                } else {
                    first_part_set.insert(FinallyAnnotations.transcript2gene[EachTranscript]);
                    second_part_set.insert(EachTranscript);
                }
            }
            first_part = concatenateSet(first_part_set);
            second_part = concatenateSet(second_part_set);

            if (Trace.is_open()) {
                const auto& reads = HighClusters[Disinform.Index2hashname[FalseNode]];
                std::string buffer;
                buffer.reserve(reads.size() * 80);      
                for (const auto& EachRead:reads) {
                    buffer += EachRead;
                    buffer += '\t';
                    buffer += "simnovel\t";
                    buffer += second_part;
                    buffer += '\t';
                    buffer += first_part;
                    buffer += '\t';
                    buffer += groupreadbarcodes.at(EachRead);
                    buffer += '\t';
                    buffer += groupreadfiles.at(EachRead);
                    buffer += '\n';
                }
                {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << buffer;
                }
            }            
        }
    }

    Eigen::VectorXd False_novel_candidate_cluster_numbers(MergeFalseSet.size());
    RowIndex = -1;
    for (const auto& FalseNode:MergeFalseSet){
        RowIndex = RowIndex + 1;
        False_novel_candidate_cluster_numbers[RowIndex] = faslenode2count[FalseNode].size();
    }

    if (known_ISM_matrix.rows() != 0){
        if (False_novel_candidate_matrix.rows() != 0){
            OutputResults.Indicate_Matrix.resize(known_ISM_matrix.rows()+False_novel_candidate_matrix.rows(), known_ISM_matrix.cols());
            OutputResults.Cluster_Number.resize(known_ISM_matrix.rows()+False_novel_candidate_matrix.rows(), known_ISM_cluster_numbers.cols());
            OutputResults.Indicate_Matrix << known_ISM_matrix,
                                             False_novel_candidate_matrix;
            OutputResults.Cluster_Number << known_ISM_cluster_numbers,
                                            False_novel_candidate_cluster_numbers;
        }
        else
        {
            OutputResults.Indicate_Matrix = known_ISM_matrix;
            OutputResults.Cluster_Number = known_ISM_cluster_numbers;
        }
    } else {
        if (False_novel_candidate_matrix.rows() != 0){
            OutputResults.Indicate_Matrix = False_novel_candidate_matrix;
            OutputResults.Cluster_Number = False_novel_candidate_cluster_numbers;
        }
    }
    return;
}


void EM_Alg_Barcodes(std::unordered_map<std::string, std::vector<double>>& filetranscriptnumber, 
                        IndicateFire& Indicate_Number, std::unordered_map<std::string, std::string>& tx2ge,
                        const std::string& Barcode, std::unordered_map<std::string, std::unordered_map<std::string, double>>& filekgenenumber,
                        std::unordered_map<std::string,int>& transcript2ID) {
    
    std::vector<double> Order_Transcript_Number(Indicate_Number.Order_Transcript_Name_Vector.size(), 0.0);
    auto& current_transcript_map = filetranscriptnumber[Barcode]; 
    int index = 0;

    for (const auto& eachIso:Indicate_Number.Order_Transcript_Name_Vector) {
        auto it = transcript2ID.find(eachIso);
        if ( it != transcript2ID.end()) {
            Order_Transcript_Number[index] = current_transcript_map[it->second];
        }
        index++;
    }

    Eigen::VectorXd Order_Transcript_Vector = Eigen::Map<Eigen::VectorXd>(Order_Transcript_Number.data(), Order_Transcript_Number.size());

    if (Indicate_Number.Indicate_Matrix.rows() != 0 && Indicate_Number.Indicate_Matrix.cols() != 0) {
        Eigen::VectorXd P_Col_init0 = Eigen::VectorXd::Constant(
            Indicate_Number.Order_Transcript_Name_Vector.size(),
            1.0 / Indicate_Number.Order_Transcript_Name_Vector.size() 
        );

        Eigen::MatrixXd Z(Indicate_Number.Indicate_Matrix.rows(), Indicate_Number.Indicate_Matrix.cols());
        Eigen::VectorXd P1(P_Col_init0.size());
        Eigen::MatrixXd AnnoN;
        int CountCyc = 1; 
        double sum_abs_diff;

        do {
            Z = Indicate_Number.Indicate_Matrix.array().rowwise() * P_Col_init0.transpose().array();
            Z.array().colwise() /= (Indicate_Number.Indicate_Matrix * P_Col_init0).array();
            AnnoN = ((Indicate_Number.Cluster_Number.transpose()) * Z).transpose().reshaped(Indicate_Number.Indicate_Matrix.cols(),1);
            P1 = AnnoN + Order_Transcript_Vector;
            P1 = P1 / P1.sum();
            sum_abs_diff = (P1 - P_Col_init0).cwiseAbs().sum();
            P_Col_init0 = P1;
            CountCyc = CountCyc + 1;
        } while ((sum_abs_diff > 5e-2) and (CountCyc < 10));

        std::string its_name; double em_count;
        
        for (int i = 0; i < AnnoN.rows(); i++){
            its_name = Indicate_Number.Order_Transcript_Name_Vector[i];
            em_count = AnnoN(i, 0);
            if (em_count > 0.0) {
                auto it = transcript2ID.find(its_name);
                if (it != transcript2ID.end()) current_transcript_map[it->second] += em_count;
            }
        }
    }

    std::string geneName; std::string second_part;
    filekgenenumber[Barcode] = {};
    auto& thisGeneBarcode = filekgenenumber[Barcode];

    for (const auto& eachAnno:transcript2ID) {
        size_t pos = eachAnno.first.find('|');
        if (pos != std::string::npos) {
            second_part = eachAnno.first.substr(pos + 1);
        } else {
            second_part = eachAnno.first;
        }
        geneName = tx2ge[eachAnno.first];
        if (geneName != "NA") {
            auto it = thisGeneBarcode.find(geneName);
            if (it != thisGeneBarcode.end()) {
                it->second = it->second + current_transcript_map[eachAnno.second];
            } else {
                thisGeneBarcode[geneName] = current_transcript_map[eachAnno.second];
            }
        }
    } 
}   


void write_transcript_file(std::unordered_map<std::string, std::vector<double>>& filetranscriptnumber,
                            OutputInformation& FinallyAnnotations, const std::set<std::string>& barcodeset, int& filenumber,
                            std::vector<std::unique_ptr<std::ofstream>>& AllIsoformFilePath,
                            std::vector<std::mutex>& AllIsoformMutexes, int& rc_threshold,
                            std::unordered_map<std::string, std::vector<std::string>>& me_transcript_se_reads,
                            std::unordered_map<std::string, std::string>& GroupReadBarcodes,
                            std::unordered_map<std::string, std::string>& GroupReadFiles,
                            std::ofstream& Updated_Files, std::string& chrname,
                            std::unordered_map<std::string, std::vector<std::array<int,2>>>& rawgtf_isoform,
                            std::unordered_map<std::string, std::string>& rawgtf_strand,
                            std::unordered_map<std::string,int>& transcript2ID) {
    
    std::string FILENUMBER = std::to_string(filenumber);
    std::unordered_map<std::string, std::unordered_map<std::string, double>> transcriptbarcodenumber;

    if (FinallyAnnotations.Transcript_Annotations.size() > 0 or me_transcript_se_reads.size() > 0) {
        for (const auto& eachAnno : FinallyAnnotations.Transcript_Annotations) {
            transcriptbarcodenumber[eachAnno.first];
        }

        for (const auto& barcode_pair : filetranscriptnumber) {
            const std::string& barcode = barcode_pair.first;
            const auto& barcodeCount = barcode_pair.second;
            
            for (const auto& tx_pair : transcript2ID) {
                const std::string& tx_name = tx_pair.first;
                double count = barcodeCount[tx_pair.second];
                transcriptbarcodenumber[tx_name][barcode] = count; 
            }
        }

        std::ostringstream gtf_buffer; 
        bool has_gtf_content = false;
        
        for (const auto& eachAnno:me_transcript_se_reads) {
            const std::string AnnoName = eachAnno.first;
            const auto& reads = eachAnno.second;
            auto it_tx = transcriptbarcodenumber.find(AnnoName);

            if (it_tx != transcriptbarcodenumber.end()) {
                for (const auto& eachRead:reads) {
                    if (GroupReadFiles[eachRead] == FILENUMBER) {
                        std::string BARCODE = GroupReadBarcodes[eachRead];
                        it_tx->second[BARCODE] += 1.0;
                    }
                }
            } else {
                size_t pos = AnnoName.find('|');
                std::string first_part, second_part;
                if (pos != std::string::npos) {
                    first_part = AnnoName.substr(0, pos);
                    second_part = AnnoName.substr(pos + 1);
                }
                if (Updated_Files.is_open()) {
                    has_gtf_content = true;
                    std::vector<std::array<int,2>>& thisVec = rawgtf_isoform[AnnoName];
                    const std::string& strand = rawgtf_strand[AnnoName];
                    gtf_buffer << chrname << "\tannotated_isoform\ttranscript\t" 
                            << thisVec[0][0] << "\t" << thisVec.back()[1] << "\t.\t" 
                            << strand << "\t.\tgene_id \"" << first_part 
                            << "\"; transcript_id \"" << second_part << "\";\n";                    
                    for (size_t i = 0; i < thisVec.size(); i++) {
                        gtf_buffer << chrname << "\tannotated_isoform\texon\t" 
                                << thisVec[i][0] << "\t" << thisVec[i][1] << "\t.\t" 
                                << strand << "\t.\tgene_id \"" << first_part 
                                << "\"; transcript_id \"" << second_part 
                                << "\"; exon_number \"" << i << "\";\n";
                    }                    
                }

                auto& current_counts = transcriptbarcodenumber[AnnoName]; 
                for (const auto& eachRead : reads) {
                    if (GroupReadFiles[eachRead] == FILENUMBER) {
                        std::string BARCODE = GroupReadBarcodes[eachRead];
                        current_counts[BARCODE] += 1.0;
                    }
                }       
            }
        }

        if (has_gtf_content && Updated_Files.is_open()) {
            std::unique_lock<std::mutex> lock(updatedGtfMutex);
            Updated_Files << gtf_buffer.str();
        }

        if (!transcriptbarcodenumber.empty()) {
            std::ostringstream main_file_buffer;
            main_file_buffer.str().reserve(transcriptbarcodenumber.size() * 1024); 

            std::string transcript_name, gene_name;
            for (const auto& eachAnno:transcriptbarcodenumber) {
                const std::string& full_name = eachAnno.first;
                const auto& counts_map = eachAnno.second;

                size_t pos = full_name.find('|');
                if (pos != std::string::npos) {
                    transcript_name = full_name.substr(pos + 1);
                    gene_name = full_name.substr(0, pos);
                } else {
                    transcript_name = full_name;
                    gene_name = "NA";
                }

                main_file_buffer << transcript_name << '\t' << gene_name;
                for (const auto& eachBar : barcodeset) {
                    double val = 0.0;
                    auto it_val = counts_map.find(eachBar);
                    if (it_val != counts_map.end()) {
                        val = it_val->second;
                    }
                    if (val < rc_threshold) {
                        main_file_buffer << '\t' << 0;
                    } else {
                        main_file_buffer << '\t' << val;
                    }
                }
                main_file_buffer << '\n';
            }
            {
                std::unique_lock<std::mutex> lock(AllIsoformMutexes[filenumber]);
                *(AllIsoformFilePath[filenumber]) << main_file_buffer.str();
            }
        }
    }
}


void write_gene_file(std::unordered_map<std::string, std::unordered_map<std::string, double>>& filegenenumber,
                    std::unordered_map<std::string, std::vector<std::string>>& segenenumber,
                    std::unordered_map<std::string, std::string>& read2barcode,
                    const std::set<std::string>& barcodeset, int& filenumber,
                    std::vector<std::unique_ptr<std::ofstream>>& AllGeneFilePath,
                    std::vector<std::mutex>& AllGeneMutexes, int& rc_threshold) {

    std::vector<std::string> barcode_list(barcodeset.begin(), barcodeset.end());
    std::unordered_map<std::string, int> barcode2idx;
    barcode2idx.reserve(barcode_list.size());
    for (int i = 0; i < barcode_list.size(); i++) {
        barcode2idx[barcode_list[i]] = i;
    }
    const int B = barcode_list.size();

    std::unordered_map<std::string, std::vector<int>> SE_gene_counts;
    SE_gene_counts.reserve(segenenumber.size());
    for (auto& kv : segenenumber) {
        const std::string& gene = kv.first;
        std::vector<int> counts(B, 0);
        for (const std::string& read : kv.second) {
            auto it = read2barcode.find(read);
            if (it != read2barcode.end()) {
                int bi = barcode2idx[it->second];
                counts[bi]++;
            }
        }
        SE_gene_counts[gene] = std::move(counts);
    }

    std::unordered_set<std::string> geneSet;
    if (!filegenenumber.empty()) {
        const auto& firstInner = filegenenumber.begin()->second;
        for (auto& g : firstInner) geneSet.insert(g.first);
    }
    for (auto& kv : SE_gene_counts) geneSet.insert(kv.first);

    std::unordered_map<std::string, std::vector<double>> matrix;
    matrix.reserve(geneSet.size());
    for (const std::string& gene : geneSet) {
        std::vector<double> counts(B, 0);

        if (!filegenenumber.empty()) {
            for (auto& bar : barcodeset) {
                auto it = filegenenumber[bar].find(gene);
                if (it != filegenenumber[bar].end()) {
                    counts[barcode2idx[bar]] = it->second;
                }
            }
        }
        auto it = SE_gene_counts.find(gene);
        if (it != SE_gene_counts.end()) {
            for (int i = 0; i < B; i++) {
                counts[i] += it->second[i];
            }
        }
        matrix[gene] = std::move(counts);
    }

    std::string out;
    out.reserve(geneSet.size() * (20 + B * 6));
    for (auto& kv : matrix) {
        out += kv.first;
        const std::vector<double>& vec = kv.second;
        for (double c : vec) {
            if (c < rc_threshold) c = 0;
            out += '\t';
            out += std::to_string((int)(c+0.5)); // ★取整了, 结果加和有点区别;
            // out += std::to_string(c);
        }
        out += '\n';
    }
    {
        std::unique_lock<std::mutex> lock(AllGeneMutexes[filenumber]);
        *(AllGeneFilePath[filenumber]) << out;
    }

}


void DetectQuant(GroupAnnotation& groupanno, SpliceChainClass& splicechainclass, GroupInformation& groupinform, 
                const int& groupdistance, unGTF& gtfexon, GTFsj& gtfsjs,
                std::ofstream& gtfFilePath, std::ofstream& traceFilePath, int& FileNo,
                std::map<std::string, std::set<std::string>>& AllFile_BarcodeSet,
                std::vector<std::unique_ptr<std::ofstream>>& IsoformFilePath,
                std::vector<std::mutex>& IsoformMutexes,
                std::vector<std::unique_ptr<std::ofstream>>& GeneFilePath,
                std::vector<std::mutex>& GeneMutexes,
                int& RC_threshold, SE_belong2_genetranscript& groupSeReads) {
   
    std::string chrchr = groupinform.chrName;
    std::string groupnumber = groupinform.GroupIndex;
    DistanceInform DMatrix_GraphNode = get_distance_matrix(groupanno.group_me_transcripts, 
                                                           splicechainclass.HighConClusters, 
                                                           groupinform.GroupReadSjs, 
                                                           groupdistance, splicechainclass.FSM);

    std::vector<std::vector<int>> CliquesVector = Find_Maximal_Cliques(DMatrix_GraphNode.NodeDistanceSet);

    DetectionResults NodeResults = Transcript_Detection(CliquesVector, DMatrix_GraphNode.Index2Anno, 
                                                        DMatrix_GraphNode.Index2Unknown, DMatrix_GraphNode.Index2Count);
    CliquesVector.clear();

    OutputInformation Finally_Annotations = Write_Detection_Transcript2gtf_AllFiles(gtfFilePath, traceFilePath, 
                                            NodeResults, DMatrix_GraphNode, 
                                            gtfsjs.mSJs[chrchr], gtfexon.GTF_transcript[chrchr], 
                                            splicechainclass.ClusterCoverage, groupanno.group_genes, 
                                            gtfexon.GTF_gene[chrchr], 
                                            gtfexon.GTF_transcript_strand[chrchr], 
                                            splicechainclass.HighStrand, chrchr, groupnumber, 
                                            splicechainclass.HighConClusters,
                                            groupinform.GroupReadFiles,
                                            groupinform.GroupReadBarcodes,
                                            AllFile_BarcodeSet,
                                            gtfexon.GTF_gene2transcript[chrchr],
                                            gtfexon.GTF_gene_strand[chrchr]); 
    splicechainclass.HighStrand.clear();

    const size_t T = Finally_Annotations.Transcript_Annotations.size();
    std::unordered_map<std::string, int> transcript_to_id; transcript_to_id.reserve(T);
    int idx = 0;
    for(const auto& t : Finally_Annotations.Transcript_Annotations) {
        transcript_to_id[t.first] = idx++;
    }

    if (FileNo > 1) {

        std::unordered_map<std::string, std::vector<double>> File_k_TranscriptNumber;
        std::unordered_map<std::string, std::unordered_map<std::string, double>> File_K_GeneNumber;
        FileBarcodeFalseNode This_File_False_Node; This_File_False_Node.FasleNode2Count.reserve(20);
        IndicateFire InitFirefly; InitFirefly.Order_Transcript_Name_Vector.reserve(T);
        for(const auto& t : Finally_Annotations.Transcript_Annotations) {
            InitFirefly.Order_Transcript_Name_Vector.push_back(t.first);
        }

        for (int k = 0; k < FileNo; k++) {
            File_k_TranscriptNumber.clear(); File_K_GeneNumber.clear();
            const std::string FileStringK = std::to_string(k);
            const auto& ALLBarcode = AllFile_BarcodeSet[FileStringK];
            File_k_TranscriptNumber.reserve(ALLBarcode.size());
            for (const auto& bc : ALLBarcode) {
                File_k_TranscriptNumber.emplace(bc, std::vector<double>(T, 0.0));
            }         

            if (!Finally_Annotations.Transcript_Annotations.empty()) {
                Solvent_sc SpliceChainSolvent = get_Solvent_FsmIsmHigh(splicechainclass, k, groupinform);
                get_transcript_init(File_k_TranscriptNumber, NodeResults, DMatrix_GraphNode, chrchr, groupnumber, k, 
                                    groupinform.GroupReadFiles, groupinform.GroupReadBarcodes,
                                    SpliceChainSolvent, ALLBarcode, transcript_to_id);
                
                Barcode_sc ThisFileBarcodeSpliceChain = get_Barcode_FsmIsmHigh(SpliceChainSolvent.File_FSM,
                                                                    SpliceChainSolvent.File_ISM,
                                                                    SpliceChainSolvent.File_HighConClusters,
                                                                    SpliceChainSolvent.File_LowConClusters,
                                                                    groupinform.GroupReadBarcodes);

                for (const auto& ThisBarcode:ALLBarcode) {
                    get_File_Barcode_False_node(NodeResults.FalseNodeSet,
                                                DMatrix_GraphNode.Index2hashname,
                                                ThisFileBarcodeSpliceChain.Barcode_HighConClusters,
                                                ThisBarcode, This_File_False_Node);

                    Quantification_initialization_Barcodes(ThisFileBarcodeSpliceChain.Barcode_ISM[ThisBarcode], 
                                                Finally_Annotations, NodeResults.TrueNodeSet,
                                                This_File_False_Node.ThisFlaseNode, 
                                                This_File_False_Node.FasleNode2Count,
                                                DMatrix_GraphNode, groupnumber, 
                                                ThisFileBarcodeSpliceChain.Barcode_HighConClusters[ThisBarcode], 
                                                groupinform.GroupReadSjs, 
                                                groupinform.GroupReadFiles,
                                                groupinform.GroupReadBarcodes,
                                                traceFilePath, InitFirefly);

                    EM_Alg_Barcodes(File_k_TranscriptNumber, InitFirefly, Finally_Annotations.transcript2gene, ThisBarcode, File_K_GeneNumber, transcript_to_id);                                       
                }
            }
            write_transcript_file(File_k_TranscriptNumber, Finally_Annotations, ALLBarcode, 
                                k, IsoformFilePath, IsoformMutexes, RC_threshold,
                                groupSeReads.Transcript_with_SE_reads, groupinform.GroupReadBarcodes, groupinform.GroupReadFiles,
                                gtfFilePath, chrchr, gtfexon.GTF_transcript[chrchr], gtfexon.GTF_transcript_strand[chrchr], transcript_to_id);

            write_gene_file(File_K_GeneNumber, groupSeReads.file_SE_reads_gene_number[k], 
                            groupinform.GroupReadBarcodes, ALLBarcode, k, 
                            GeneFilePath, GeneMutexes, RC_threshold);
        }

    } else {
        int k = 0;
        const auto& ALLBarcode = AllFile_BarcodeSet["0"];
        std::unordered_map<std::string, std::vector<double>> File_k_TranscriptNumber;
        File_k_TranscriptNumber.reserve(ALLBarcode.size());
        for (const auto& bc : ALLBarcode) {
            File_k_TranscriptNumber.emplace(bc, std::vector<double>(T, 0.0));
        }
        std::unordered_map<std::string, std::unordered_map<std::string, double>> File_K_GeneNumber;

        if (!Finally_Annotations.Transcript_Annotations.empty()) {
            Solvent_sc SpliceChainSolvent = get_Solvent_FsmIsmHigh(splicechainclass, k, groupinform);
            get_transcript_init(File_k_TranscriptNumber, NodeResults, DMatrix_GraphNode, chrchr, groupnumber, k, 
                                groupinform.GroupReadFiles, groupinform.GroupReadBarcodes,
                                SpliceChainSolvent, ALLBarcode, transcript_to_id);

            Barcode_sc ThisFileBarcodeSpliceChain = get_Barcode_FsmIsmHigh(SpliceChainSolvent.File_FSM,
                                                                SpliceChainSolvent.File_ISM,
                                                                SpliceChainSolvent.File_HighConClusters,
                                                                SpliceChainSolvent.File_LowConClusters,
                                                                groupinform.GroupReadBarcodes);

            FileBarcodeFalseNode This_File_False_Node; This_File_False_Node.FasleNode2Count.reserve(20);
            IndicateFire InitFirefly; InitFirefly.Order_Transcript_Name_Vector.reserve(T);
            for(const auto& t : Finally_Annotations.Transcript_Annotations) {
                InitFirefly.Order_Transcript_Name_Vector.push_back(t.first);
            }
            
            for (const auto& ThisBarcode:ALLBarcode) {
                get_File_Barcode_False_node(NodeResults.FalseNodeSet,
                                            DMatrix_GraphNode.Index2hashname,
                                            ThisFileBarcodeSpliceChain.Barcode_HighConClusters,
                                            ThisBarcode, This_File_False_Node); 

                Quantification_initialization_Barcodes(ThisFileBarcodeSpliceChain.Barcode_ISM[ThisBarcode], 
                                                        Finally_Annotations, NodeResults.TrueNodeSet,
                                                        This_File_False_Node.ThisFlaseNode, 
                                                        This_File_False_Node.FasleNode2Count,
                                                        DMatrix_GraphNode, groupnumber, 
                                                        ThisFileBarcodeSpliceChain.Barcode_HighConClusters[ThisBarcode], 
                                                        groupinform.GroupReadSjs, 
                                                        groupinform.GroupReadFiles,
                                                        groupinform.GroupReadBarcodes,
                                                        traceFilePath, InitFirefly);

                EM_Alg_Barcodes(File_k_TranscriptNumber, InitFirefly, Finally_Annotations.transcript2gene, ThisBarcode, File_K_GeneNumber, transcript_to_id);
            }
        }

        write_transcript_file(File_k_TranscriptNumber, Finally_Annotations, ALLBarcode, 
                            k, IsoformFilePath, IsoformMutexes, RC_threshold,
                            groupSeReads.Transcript_with_SE_reads, groupinform.GroupReadBarcodes, groupinform.GroupReadFiles,
                            gtfFilePath, chrchr, gtfexon.GTF_transcript[chrchr], gtfexon.GTF_transcript_strand[chrchr], transcript_to_id);

        write_gene_file(File_K_GeneNumber, groupSeReads.file_SE_reads_gene_number[k], groupinform.GroupReadBarcodes, 
                        ALLBarcode, k, GeneFilePath, GeneMutexes, RC_threshold);

    }
    
}


void processGroup(std::streampos& start, std::streampos& end, 
                  const std::string& sam_file_path, 
                  const int& Sj_supportReadNumber, unGTF& gtf_full, 
                  GTFsj& gtf_splice, const int& GraphDis,
                  std::ofstream& updatedgtffile, std::ofstream& tracefile,
                  std::vector<std::unique_ptr<std::ofstream>>& isoformfilePath,
                  std::vector<std::mutex>& isoformMutexes,
                  std::vector<std::unique_ptr<std::ofstream>>& genefilePath,
                  std::vector<std::mutex>& geneMutexes,
                  int& fileno, std::string& outputPath,
                  std::map<std::string, std::set<std::string>>& filebarcodeSet,
                  int& singleEdge, int& ReadCount_threshold) {
    
    GroupInformation group_information = knowGroupInformation(start, end, sam_file_path, Sj_supportReadNumber);
    std::string chrchr = group_information.chrName;
    std::array<int,2> groupcoverage = group_information.GroupCoverage;
    
    SE_belong2_genetranscript group_se_reads;
    GroupAnnotation group_annotation;
    SpliceChainClass spliceclass;
    
    if ( (!group_information.GroupSingleExon.empty()) or (!group_information.GroupReadSjs.empty()) ) {
        group_annotation = get_group_single_exon_gene_annotation(gtf_full.GTF_gene[chrchr], groupcoverage,
                                gtf_full.GTF_gene2transcript[chrchr], gtf_full.GTF_transcript[chrchr], gtf_splice.mSJs[chrchr], chrchr);
    }

    if (!group_information.GroupReadSjs.empty()) {
        spliceclass = generate_splice_chain_class(group_information.GroupReadSjs, group_information.GroupReadCoverage, 
                                                group_annotation.group_me_SJs, group_information.GroupSigns, 
                                                gtf_splice.mSJsBE[chrchr], group_information.GroupReadFiles, group_information.GroupReadBarcodes, 
                                                group_information.GroupReadUMIs, tracefile, Sj_supportReadNumber);
        group_information.GroupSigns.clear();
        group_information.GroupReadCoverage.clear();
    }

    if (!group_information.GroupSingleExon.empty()) {
        group_se_reads = get_group_singleexon_reads_2gene(group_annotation, group_information.GroupSingleExon, 
                                        group_information.GroupReadFiles, gtf_full.GTF_gene2transcript[chrchr], 
                                        gtf_full.GTF_transcript[chrchr], fileno, singleEdge, chrchr, spliceclass.FSM, tracefile,
                                        group_information.GroupReadBarcodes, group_information.GroupReadUMIs);

        write_single_exon_gtf_trace_sc(fileno, group_information.GroupReadFiles, group_se_reads.Transcript_with_SE_reads, 
                                    group_annotation.group_se_transcripts, updatedgtffile, tracefile, isoformfilePath, isoformMutexes, 
                                    gtf_full.GTF_transcript_strand[chrchr], group_information.GroupReadBarcodes, 
                                    chrchr, filebarcodeSet);    
    }
    group_information.GroupSingleExon.clear();

    DetectQuant(group_annotation, spliceclass, group_information, GraphDis, gtf_full, gtf_splice, 
                updatedgtffile, tracefile, fileno,
                filebarcodeSet, isoformfilePath, isoformMutexes, genefilePath, geneMutexes,
                ReadCount_threshold, group_se_reads);

    spliceclass.FSM.clear(); spliceclass.HighConClusters.clear();spliceclass.ISM.clear();
    spliceclass.LowConClusters.clear();spliceclass.HighStrand.clear();

    if (group_information.GroupReadSjs.size() + group_information.GroupSingleExon.size() > 10000) {
        std::unique_lock<std::mutex> lock(bigMutex);
        std::cerr << (group_information.GroupReadSjs.size() + group_information.GroupSingleExon.size()) / ((double)1000000) << " M reads processed..\n";
    }
}



int main(int argc, char* argv[])
{   
    std::ios_base::sync_with_stdio(false);
    int option_index = 0;
    int c;

    std::string samfile_name;
    std::string fastafile_name;
    std::string gtffile_name;
    std::string output_file_name;

    int SJDistance = 18;
    int SJ_support_read_number = 2;
    int Graph_distance = 60;
    int Thread = 10;
    int Read_count = 0;
    int single_exon_edge = 500;
    std::string umi = "U8:Z:";
    std::string barcode = "BC:Z:";
    int mapq = 0;

    std::cerr << R"(
     ____               ____   ___  _      ___ 
    | __ )  _ __ ___   / ___| / _ \| |    |_ _|
    |  _ \ | '__/ _ \ | |    | | | | |     | | 
    | |_) || | | (_) || |___ | |_| | |___  | | 
    |____/ |_|  \___/  \____| \___/|_____||___|
    )" << std::endl;

    std::cerr << "         BroCOLI  Version: 1.0.0" << std::endl;

    while ((c = getopt_long(argc, argv, "s:f:g:o:j:n:m:e:d:t:u:b:r:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 's':
                samfile_name = optarg;
                break;
            case 'f':
                fastafile_name = optarg;
                break;
            case 'g':
                gtffile_name = optarg;
                break;
            case 'o':
                output_file_name = optarg;
                break;
            case 'j':
                SJDistance = std::stoi(optarg);  
                break;
            case 'n':
                SJ_support_read_number = std::stoi(optarg);   
                break;
            case 'm':
                mapq = std::stoi(optarg);   
                break;
            case 'e':
                single_exon_edge = std::stoi(optarg);   
                break;
            case 'd':
                Graph_distance = std::stoi(optarg);   
                break;
            case 't':
                Thread = std::stoi(optarg);   
                break;
            case 'u':
                umi = optarg;   
                break;
            case 'b':
                barcode = optarg;  
                break;
            case 'r':
                Read_count = std::stoi(optarg);   
                break;            
            case 'h':
                print_usage(argv[0]);
                exit(EXIT_FAILURE);              
            case '?':
                std::cerr << "Invalid option" << std::endl;
                exit(EXIT_FAILURE);
            default:
                break;
        }
    }

    if (samfile_name.empty() || fastafile_name.empty() || output_file_name.empty()) {
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> sam_file_vec = traverse_sam_file(samfile_name, output_file_name);

    std::vector<std::string> outputFileVec = check_catalog_exist(output_file_name);
    std::ofstream gtf_file(outputFileVec[0], std::ios::app);
    std::ofstream trace_file(outputFileVec[1], std::ios::app);

    std::cerr << "*****" << std::endl;
    std::cerr << "Input file: " << samfile_name << std::endl;
    std::cerr << "FASTA file: " << fastafile_name << std::endl;
    std::cerr << "GTF file: " << gtffile_name << std::endl;
    std::cerr << "Output file: " << output_file_name << std::endl;
    std::cerr << "Single exon boundary : " << single_exon_edge << std::endl;
    std::cerr << "SJ Distance: " << SJDistance << std::endl;
    std::cerr << "SJ support read number: " << SJ_support_read_number << std::endl;
    std::cerr << "MAPQ: " << mapq << std::endl;
    std::cerr << "Graph distance: " << Graph_distance << std::endl;
    std::cerr << "Thread: " << Thread << std::endl;
    std::cerr << "Output min read count: " << Read_count << std::endl;
    std::cerr << "*****" << std::endl;

    std::cerr << "*** " << "Read and process the files ...... " << std::endl;
    std::unordered_map<std::string, std::string> Fasta = Read_fasta_file(fastafile_name);

    auto start = std::chrono::high_resolution_clock::now();
    FileSplit BroCOLIfile = thread_all_read_sam_files(samfile_name, sam_file_vec, Thread, output_file_name, SJDistance, Fasta, mapq, umi, barcode);
    Fasta.clear();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cerr << "Read file information cost time = " << diff.count() << " s\n";     

    std::vector<std::unique_ptr<std::ofstream>> BroCOLIQuantfile = write_quantification_files(output_file_name, BroCOLIfile.file_barcodeSet);
    std::vector<std::mutex> TranscriptMutexes(BroCOLIQuantfile.size());
    std::vector<std::unique_ptr<std::ofstream>> BroCOLIGenefile = write_quantification_Gene_files(output_file_name, BroCOLIfile.file_barcodeSet);
    std::vector<std::mutex> GeneMutexes(BroCOLIGenefile.size());

    unGTF GTF_full = get_gtf_annotation(gtffile_name);
    GTFsj GTF_Splice = get_SJs_SE(GTF_full.GTF_transcript);
    std::vector<std::size_t> Group_idx = sort_indexes_e(BroCOLIfile.group_reads_number);
    std::cerr << "*** File processing completed! ***\n";

    ThreadPool BroCOLIpool(Thread);
    std::vector<std::future<void>> futures;

    std::cerr << "*** " << "Transcript identification and quantification ......\n";
    std::cerr << "*** " << "Only output the completion status of the larger clusters ......\n";
    std::cerr << "*** " << "Clusters with less than 10k reads will not be output ......\n";
    start = std::chrono::high_resolution_clock::now();
    for (const auto& i:Group_idx) {
        futures.emplace_back(BroCOLIpool.enqueue([&, i]() { 
            processGroup(
                BroCOLIfile.reads_pointer[2*i], 
                BroCOLIfile.reads_pointer[2*i+1], 
                BroCOLIfile.readtxt_path, 
                SJ_support_read_number, 
                std::ref(GTF_full), 
                std::ref(GTF_Splice),
                Graph_distance,
                gtf_file,
                trace_file,
                BroCOLIQuantfile,
                TranscriptMutexes,
                BroCOLIGenefile,
                GeneMutexes,
                BroCOLIfile.FileNo,
                std::ref(output_file_name),
                std::ref(BroCOLIfile.file_barcodeSet),
                single_exon_edge,
                Read_count);
            }));
    }
    
    for (auto& future : futures) {
        future.get();  
    }
    end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cerr << "Identification and quantification cost time = " << diff.count() << " s\n";

    gtf_file.close();
    trace_file.close();
    std::cerr << "*** BroCOLI quantification has been successfully completed! ***\n";


    return 0;
}



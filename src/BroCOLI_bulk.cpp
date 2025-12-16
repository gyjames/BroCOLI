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

std::mutex updatedGtfMutex;
std::mutex isoformCountMutex;
std::mutex geneCountMutex;
std::mutex traceMutex;
std::mutex bigMutex;

struct Split_Result{
    std::string read_name;
    std::vector<std::string> tokens;
    int read_length;
};


Split_Result get_string_split_fast(const std::string& s) {
    Split_Result r;
    r.read_length = 0;

    int field = 0;
    size_t start = 0;
    size_t i = 0;

    for (; i <= s.size(); ++i) {
        if (i == s.size() || s[i] == '\t') {
            const char* p = s.data() + start;
            size_t len = i - start;
            if (field == 0) {
                r.read_name.assign(p, len);
            } else if (field < 7) {
                r.tokens.emplace_back(p, len);
            } else if (field == 9) {
                r.read_length = len; 
            }
            field++;
            start = i + 1;
        }
    }
    return r;
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

    int pos = std::stoi(initpos);
    int number = 0;
    for (char c : CIGARvalue) {
        if (c >= '0' && c <= '9') {
            number = number * 10 + (c - '0');
        } else {
            switch (c) {
            case 'M':
                out.ReadIntervals.push_back({pos, pos + number});
                pos += number;
                out.ReadMatchLength += number;
                break;
            case 'D':
            case 'N':
                pos += number;
                break;
            case 'I':
            case 'S':
                out.ReadMatchLength += number;
                break;
            case 'H':
            case 'P':
                break;
            }
            number = 0;
        }
    }
    return out;
}



int IntervalIntersection(const std::vector<std::array<int,2>>& firstList, const std::vector<std::array<int,2>>& secondList) {
    std::vector<std::array<int,2>> res; //定义返回的结果;
    int result = 0;
    int i = 0;
    int j = 0;

    while(i < firstList.size() && j < secondList.size()){
        int start = std::max(firstList[i][0], secondList[j][0]);
        int end = std::min(firstList[i][1], secondList[j][1]);
        if(start <= end){
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




int IntervalMerge(const std::vector<std::array<int,2>>& firstList, const std::vector<std::array<int,2>>& secondList) {
    //将两个序列集合合并为一个;
    std::vector<std::array<int,2>> intervals;
    intervals.resize(firstList.size()+secondList.size());
    merge(firstList.begin(), firstList.end(), secondList.begin(), secondList.end(), intervals.begin());
    
    sort(intervals.begin(), intervals.end()); //对所有的区间从小到大排序;
    std::vector<std::array<int,2>> merged; //初始化输出结果;
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
                // Overlapping, don't include this interval
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
                // Overlapping, don't include this interval
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


struct Reads_Clusters {
    std::streampos lastPos;
    std::streampos newPos;
    std::map<std::string, std::vector<std::array<int,2>>> Mymap;
    std::map<std::string, int> Mylen;
    std::map<std::string, int> MyFlag;
    std::string SetRef_name;
    std::array<int,2> ClusterCoverage;

    int MAPQ = 0;
    int mapqless1 = 0;
    int mappingdiff = 0;
    int surveyNum = 0;    
}; 


Reads_Clusters get_each_cluster_reads(std::ifstream& samfile, std::streampos CurrentPos, std::streampos EndPos, const int& MAPQ) {

    Reads_Clusters NewCluster {};
    std::map<std::string, std::vector<std::array<int,2>>> read_informs;
    std::map<std::string, int> read_len;
    std::map<std::string, int> read_flag; // 1 + ; 0 -
    std::string line;
    std::string now_gene;
    std::string last_chr;

    samfile.seekg(CurrentPos, std::ios::beg);
    std::streampos earlyPos = CurrentPos;
    int early_begin_pos = 0;
    int early_end_pos = 0;

    Split_Result This_Line {};
    Read_intervals_and_Mlength CIGAR_interval {};

	while (getline(samfile, line))
	{
        earlyPos = CurrentPos; 
        CurrentPos = samfile.tellg();

        if (line[0] != '@'){
            NewCluster.surveyNum = NewCluster.surveyNum + 1;
            // This_Line = get_string_split(line, '\t');         
            This_Line = get_string_split_fast(line);  

            int read_mapq = std::stoi(This_Line.tokens[3]) ;
            if (read_mapq >= MAPQ) {
                
                CIGAR_interval = get_read_intervals_fast(This_Line.tokens[4], This_Line.tokens[2]);
                if (This_Line.read_length == CIGAR_interval.ReadMatchLength){

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
                            NewCluster.surveyNum = NewCluster.surveyNum - 1;
                            break;
                        } else{

                            if ((read_in_gene_begin_pos - early_end_pos > 0) || (early_begin_pos - read_in_gene_end_pos > 0)){

                                NewCluster.lastPos = earlyPos;
                                NewCluster.newPos = CurrentPos;
                                NewCluster.Mymap = read_informs;
                                NewCluster.Mylen = read_len;
                                NewCluster.MyFlag = read_flag;
                                NewCluster.ClusterCoverage[0] = early_begin_pos;
                                NewCluster.ClusterCoverage[1] = early_end_pos;
                                NewCluster.SetRef_name = last_chr;
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
                    read_len[This_Line.read_name] = This_Line.read_length;
                    int Flag = std::stoi(This_Line.tokens[0]);
                    if (Flag & 16) {
                        read_flag[This_Line.read_name] = 0; // -
                    } else {
                        read_flag[This_Line.read_name] = 1; // +
                    }

                } else {
                    NewCluster.mappingdiff = NewCluster.mappingdiff + 1;
                }
                
            } else {
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





struct GTFsj{
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

    std::cerr << "***** Now open Fasta file: " << fastafile_name << "! *****" << std::endl;

    std::ifstream FastaFile; 

	FastaFile.open(fastafile_name);

	if (!FastaFile.is_open()) {
		std::cerr << "We can't open Fasta file !" << fastafile_name << std::endl;
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

        std::cerr << "***** The Fasta file has been read to the end ! *****" << std::endl;
    } 
    else if (FastaFile.fail()) {
        std::cerr << "File FALSE !" << std::endl;
    }
    else {std::cerr << "unkown reason" << std::endl;}

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
    size_t pos1 = Sj_Array[0] - 1;
    if (pos1 + 2 > Fasta_chr.size()) {
        std::cerr << "ERROR: signal1 substring out of range! " << std::endl;
        std::cerr << "ERROR perhaps it is because the reference used during the mapping process does not match the given reference! " << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::string signal1 = Fasta_chr.substr(pos1, 2);


    size_t pos2 = Sj_Array[1] - 2;
    if (pos2 + 2 > Fasta_chr.size()) {
        std::cerr << "ERROR: signal2 substring out of range!" << std::endl;
        std::cerr << "ERROR perhaps it is because the reference used during the mapping process does not match the given reference! " << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::string signal2 = Fasta_chr.substr(pos2, 2);

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
    {"mode", required_argument, 0, 'm'},
    {"SJDistance", required_argument, 0, 'j'},  
    {"support", required_argument, 0, 'n'},   
    {"single_exon_boundary", required_argument, 0, 'e'},
    {"graph_distance", required_argument, 0, 'd'},  
    {"thread", required_argument, 0, 't'}, 
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
    std::cout << "  -m, --mapq                  Mapping quality." << std::endl;
    std::cout << "  -j, --SJDistance            the minimum distance determined as intron. (optional, default:18)" << std::endl;
    std::cout << "  -n, --support               min perfect read count for all splice junctions of novel isoform. (optional, default:2)" << std::endl;
    std::cout << "  -e, --single_exon_boundary  belongs to the isoform scope of a single exon. (optional, default:60)" << std::endl;
    std::cout << "  -d, --graph_distance        the distance threshold for constructing the isoform candidate distance graph. (optional, default:60)" << std::endl;
    std::cout << "  -t, --thread                thread number. (optional, default:8)" << std::endl;
    std::cout << "  -r, --output_min_read_count the minimum number of transcripts in the output results. (optional, int, default:1)" << std::endl;
    std::cout << "  -h, --help                  show this help information." << std::endl;
}

std::vector<std::string> check_catalog_exist(const std::string& output_path,
                                            std::vector<std::string>& samFileVec) {
    std::vector<std::string> outputFileVector;

    if (!directoryExists(output_path)) {
        std::cerr << "Directory does not exist, creating it..." << std::endl;

        if (createDirectory(output_path)) {
            std::cerr << "Directory created successfully: " << output_path << std::endl;
        } else {
            std::cerr << "Failed to create directory: " << output_path << std::endl;
            exit(EXIT_FAILURE);
        }
    } else {
        std::cerr << "Directory already exists: " << output_path << std::endl;
    } 


    std::string updatedGtfPath = joinPath(output_path, "updated_annotitions.gtf");
    outputFileVector.push_back(updatedGtfPath);
    std::ofstream gtf_file(updatedGtfPath, std::ios::trunc);
    gtf_file.close();

    std::string quantTranscriptPath = joinPath(output_path, "counts_transcript.txt");
    outputFileVector.push_back(quantTranscriptPath);
    std::ofstream Q_output_Transcript(quantTranscriptPath, std::ios::trunc);
    
    std::string quantGenePath = joinPath(output_path, "counts_gene.txt");
    outputFileVector.push_back(quantGenePath);
    std::ofstream Q_output_Gene(quantGenePath, std::ios::trunc);
    
    std::string TracePath = joinPath(output_path, "compatible_isoform.tsv");
    outputFileVector.push_back(TracePath);
    std::ofstream TraceIsoform(TracePath, std::ios::trunc);

    if (samFileVec.size() == 1){
        if (Q_output_Transcript.is_open()){
            Q_output_Transcript << "transcript_id" << '\t' << "gene_id" << '\t' << samFileVec[0] << '\n';
        }
        if (Q_output_Gene.is_open()){
            Q_output_Gene << "gene_id" << '\t' << samFileVec[0] << '\n';
        }
        if (TraceIsoform.is_open()){
            TraceIsoform << "read_id" << '\t' << "category" << '\t' << "isoform_id" << '\t' << "gene_id" << '\t' << "file" << '\n';
        }
    } else {
        if (Q_output_Transcript.is_open()){
            Q_output_Transcript << "transcript_id" << '\t' << "gene_id";
        }
        if (Q_output_Gene.is_open()){
            Q_output_Gene << "gene_id";
        }        
        if (TraceIsoform.is_open()){
            TraceIsoform << "read_id" << '\t' << "category" << '\t' << "isoform_id" << '\t' << "gene_id" << '\t' << "file" << '\n';
        }
        for (const auto& eachFile:samFileVec) {
            Q_output_Transcript << '\t' << eachFile;
            Q_output_Gene << '\t' << eachFile;
        }
        Q_output_Transcript << '\n';
        Q_output_Gene << '\n';
    }
    Q_output_Transcript.close();
    Q_output_Gene.close(); 
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
            std::cerr << "* This is a txt/tsv file! * " << sam_file_path << "\n";
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
        while ((entry = readdir(dir)) == readdir(dir)) {
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
        std::cerr << "SAM File " << i << " : " << sam_file_vector[i] << std::endl;
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
                  std::unordered_map<std::string, std::string>& FastaRef, const int& mapq) {

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

    // 打开输出文件 —— 保持文本模式；我们会用 write() 写整个 buffer
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
        each_cluster_informs = get_each_cluster_reads(samfile, Last_Position, end, mapq);
        chrchr = each_cluster_informs.SetRef_name;
        Last_Position = each_cluster_informs.lastPos;
        Current_Position = each_cluster_informs.newPos;
        each_read_SJs_informs = get_reads_allSJs(each_cluster_informs.Mymap, SJ_Distance);

        // 缓存常用引用，减少 map 查找
        auto &reads_SJs = each_read_SJs_informs.reads_SJs;
        auto &reads_begin_end = each_read_SJs_informs.reads_begin_end;
        auto &reads_SJs_left = each_read_SJs_informs.reads_SJs_left;
        auto &reads_SJs_right = each_read_SJs_informs.reads_SJs_right;
        auto &reads_single_exon = each_read_SJs_informs.reads_single_exon;

        auto &cluster_Mylen = each_cluster_informs.Mylen;
        auto &cluster_MyFlag = each_cluster_informs.MyFlag;
        const std::array<int,2> &clusterCoverage = each_cluster_informs.ClusterCoverage;

        // 遍历多外显子 reads（含 SJ）
        for (const auto& eachRead : reads_SJs) {
            const std::string &read_id = eachRead.first;
            // 尽量用 find 一次取得 Mylen / MyFlag，避免用 operator[]（可能插入）
            int read_len = 0;
            int read_flag = 0;
            auto it_len = cluster_Mylen.find(read_id);
            if (it_len != cluster_Mylen.end()) read_len = it_len->second;
            auto it_flag = cluster_MyFlag.find(read_id);
            if (it_flag != cluster_MyFlag.end()) read_flag = it_flag->second;
            auto it_be = reads_begin_end.find(read_id);
            if (it_be != reads_begin_end.end()) {
                readBeginEnd = it_be->second;
            } else {
                readBeginEnd = {0,0};
            }
            // 拼接行（不要频繁写流）
            outBuf.append(read_id);
            outBuf.push_back('\t');
            outBuf.append(chrchr);
            outBuf.push_back('\t');
            outBuf.append(std::to_string(Group_index));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(clusterCoverage[0]));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(clusterCoverage[1]));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(readBeginEnd[0]));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(readBeginEnd[1]));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(read_len));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(read_flag));

            SjNumber = 0;
            for (const auto& Sj : eachRead.second) {
                outBuf.push_back('\t');
                outBuf.append(std::to_string(Sj[0]));
                outBuf.push_back('\t');
                outBuf.append(std::to_string(Sj[1]));

                // 判断 splice 信号（仍保持原逻辑）
                bool sjOk = false;
                auto it_left = reads_SJs_left.find(read_id);
                auto it_right = reads_SJs_right.find(read_id);
                if (it_left != reads_SJs_left.end() && it_right != reads_SJs_right.end()) {
                    sjOk = ifSjNoError(it_left->second, it_right->second, SjNumber);
                }
                if (sjOk) {
                    // 这里 access FastaRef[chrchr] 可能比较昂贵，先用 find
                    auto itF = FastaRef.find(chrchr);
                    if (itF != FastaRef.end()) {
                        thisSjSignal = ifSjSignal(itF->second, Sj, read_flag);
                    } else {
                        thisSjSignal = 0;
                    }

                    if (thisSjSignal == 1) {
                        outBuf.push_back('\t');
                        outBuf.append("1");
                    } else if (thisSjSignal == 2) {
                        outBuf.push_back('\t');
                        outBuf.append("2");
                    } else {
                        outBuf.push_back('\t');
                        outBuf.append("0");
                    }
                } else {
                    outBuf.push_back('\t');
                    outBuf.append("0");
                }
                SjNumber++;
            }
            outBuf.push_back('\n');
            // 达到阈值就一次性写磁盘并清空 buffer
            if (outBuf.size() >= FLUSH_THRESHOLD) {
                ReadInform.write(outBuf.data(), static_cast<std::streamsize>(outBuf.size()));
                outBuf.clear();
            }
        }

        // single exon reads
        for (const auto& eachRead : reads_single_exon) {
            const std::string &read_id = eachRead.first;
            int read_len = 0;
            int read_flag = 0;
            auto it_len = cluster_Mylen.find(read_id);
            if (it_len != cluster_Mylen.end()) read_len = it_len->second;
            auto it_flag = cluster_MyFlag.find(read_id);
            if (it_flag != cluster_MyFlag.end()) read_flag = it_flag->second;
            readBeginEnd = eachRead.second;
            outBuf.append(read_id);
            outBuf.push_back('\t');
            outBuf.append(chrchr);
            outBuf.push_back('\t');
            outBuf.append(std::to_string(Group_index));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(clusterCoverage[0]));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(clusterCoverage[1]));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(readBeginEnd[0]));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(readBeginEnd[1]));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(read_len));
            outBuf.push_back('\t');
            outBuf.append(std::to_string(read_flag));
            outBuf.push_back('\n');
            if (outBuf.size() >= FLUSH_THRESHOLD) {
                ReadInform.write(outBuf.data(), static_cast<std::streamsize>(outBuf.size()));
                outBuf.clear();
            }
        }
    } // while
    // flush 剩余数据
    if (!outBuf.empty()) {
        ReadInform.write(outBuf.data(), static_cast<std::streamsize>(outBuf.size()));
        outBuf.clear();
    }
    samfile.close();
    ReadInform.close();
    std::cerr << "^-^ One of the threads for splitting the file has completed its processing! ^-^" << std::endl;
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

    std::string readtxt_path; 
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

    std::map<std::string, std::vector<std::string>> Group_temporary;
    std::array<int,2> this_coverage;

    for (int file_number = 0; file_number < Read_x_vec.size(); file_number++){
        std::string fileName = SmallFilePath + "/" + Read_x_vec[file_number];
        std::ifstream SmallSamFile(fileName);
        earlyGroup = "+";
        thisGroup = "+";

        while (getline(SmallSamFile, line)) {
            Split_Result This_Line = get_line_split(line, '\t');
            thisGroup = This_Line.tokens[1];
            this_coverage = {std::stoi(This_Line.tokens[2]), std::stoi(This_Line.tokens[3])};

            if (thisGroup == earlyGroup) {

                auto it = Group_temporary.find(This_Line.read_name);
                if (it != Group_temporary.end()){
                    if (std::stoi(This_Line.tokens[6]) > std::stoi(Group_temporary[This_Line.read_name][6])){
                        Group_temporary[This_Line.read_name] = This_Line.tokens;
                    }
                } else {
                    Group_temporary[This_Line.read_name] = This_Line.tokens;
                }
            } else { 

                if (Group_temporary.size() != 0) {

                    if (last_chr == This_Line.tokens[0]) {

                        if (Group_low <= this_coverage[1] && this_coverage[0] <= Group_high){

                            auto it = Group_temporary.find(This_Line.read_name);
                            if (it != Group_temporary.end()){
                                if (std::stoi(This_Line.tokens[6]) > std::stoi(Group_temporary[This_Line.read_name][6])){
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
                            Chunk_Bang.reads_pointer.push_back(Readpos); 
                            Chunk_Bang.reads_pointer.push_back(Readpos); 
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
                        Chunk_Bang.reads_pointer.push_back(Readpos); 
                        Chunk_Bang.reads_pointer.push_back(Readpos);  
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
    return Chunk_Bang;
}


std::map<std::string, std::vector<std::array<int,2>>> Merge_Read_Interval(std::map<int, std::map<std::string, std::vector<std::array<int,2>>>>& FileCoverage) {

    std::map<std::string, std::vector<std::array<int,2>>> chr_range_temp;
    std::vector<std::array<int,2>> raw_intervals;
    std::vector<std::array<int,2>> sorted_intervals;
    if (FileCoverage.size() != 0) {

        for (const auto& eachFile:FileCoverage) {
            if (eachFile.second.size() != 0) {

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
    if (chr_range_temp.size() != 0) {
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

            std::string ReadFilePath = output_path + "/sam_" + std::to_string(i) + "/All_Read.txt";

            std::ifstream samfile_one(ReadFilePath);

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
            group_read_pointer.push_back(ThisFile_ThisGroup_pointer); 
            group_read_pointer.push_back(ThisFile_ThisGroup_pointer); 
        }
    } else {
        std::cerr << "^-^ *** Merge large files failed *** ^-^" << std::endl;
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
            int Rhino = 0; 
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
                                    const int& MAPQ) {
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
    if (sam_file_vec.size() != 0) {

        ThreadPool Bigfilepool(numThreads);
        std::vector<std::future<void>> myJobs;
        for (int samFileNumber = 0; samFileNumber < sam_file_vec.size(); ++samFileNumber) {
            
            startposVec.clear();
            endposVec.clear();

            std::cerr << "*** " << "Start processing SAM File " << samFileNumber << " ***" << std::endl;
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
                        MAPQ
                    );
                }));
            }
            for (auto& future : myJobs) {
                future.get();
            }
            std::cerr << "^-^ [" << sam_file_vec[samFileNumber] << "] All threads are finished generating small files! ^-^" << std::endl;
            std::cerr << "^-^ Start of merge small files ! ^-^" << std::endl;
            BigBang = Merge_Read_Small_Files(chunkFilePath, samFileNumber);
            std::cerr << "^-^ End of merge small files ! ^-^" << std::endl;
            
            if (sam_file_vec.size() > 1) {
                File_chr_coverage[samFileNumber] = BigBang.chr_coverage;
                File_group_pointer[samFileNumber] = BigBang.coverage2pos;

            }

        }
        if (sam_file_vec.size() > 1) {
            std::cerr << "^-^ The number of sam files is greater than one. Large files need to be merged. ^-^" << std::endl;
            std::map<std::string, std::vector<std::array<int,2>>> ChrCoverage = Merge_Read_Interval(File_chr_coverage);
            std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>> ChrCoverage_SmallPointer = 
                                                                        get_pointers(ChrCoverage, File_chr_coverage, File_group_pointer);
            std::vector<std::streampos> Allreads_pointer;
            std::vector<int> Allreads_group_number;            

            std::string FinallyFile = joinPath(outputPath, "MergedRead.txt");
            std::ofstream FinallyReadInform(FinallyFile, std::ios::trunc); 
            std::streampos Readpos = FinallyReadInform.tellp();
            Allreads_pointer.push_back(Readpos);
            myJobs.clear();
            int new_group = 0;

            std::cerr << "^-^ Start of merge large File ! ^-^" << std::endl;

            for (const auto& eachChr:ChrCoverage_SmallPointer) {
                std::string ChrName = eachChr.first;
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

            std::cerr << "^-^ End of merge large File ! ^-^" << std::endl;

            BigBang.readtxt_path = FinallyFile;
            BigBang.group_reads_number = Allreads_group_number;
            BigBang.reads_pointer = Allreads_pointer;
            BigBang.FileNo = sam_file_vec.size();

            FinallyReadInform.close();
        }
    } else {
        std::cerr << "* There are no files in the folder! *" << strerror(errno) << std::endl;
    }
    return BigBang;
}



struct Line_Split{
    std::string read_file; //0
    std::string read_name; //1
    std::string chr_name; //2
    std::string Group_index; //3
    std::array<int,2> Group_coverage; //4,5
    std::array<int,2> read_coverage; //6,7
    int read_length; // 这个就是追踪文件的时候可以画Coverage; //8
    int read_strand; //9
    std::vector<std::array<int,2>> read_SJ;
    std::vector<int> read_sj_quality;
};

Line_Split MakeLineSplit(const std::string& s, char delimiter){
    // 初始化返回值结构体;
    Line_Split read_result = {};
    // 字符串流;
    std::istringstream iss(s);
    // 定义字符串流中的每一段字符串;
    std::string token;
    int number = 0;
    // 用于二次处理;
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
            read_result.read_strand = std::stoi(token);
        } else if (number > 10) {
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
    std::string chrName; //染色体;
    std::array<int,2> GroupCoverage; //Group范围;
    std::unordered_map<std::string, std::string> GroupReadFiles; //Group每条reads的所属文件; (公共, single exon 与 Multi exon都有);
    std::unordered_map<std::string, std::vector<std::array<int,2>>> GroupReadSjs; //Group每条reads的SJs;
    std::unordered_map<std::string, std::array<int,2>> GroupReadCoverage; //Group每条reads的开头结尾;
    std::map<std::array<int,2>, int> GroupSjs; //Group中根据Reads得到的高置信SJ;
    std::unordered_map<std::string, std::vector<int>> GroupSigns; //Group中所有reads的信号;
    std::unordered_map<std::string, std::array<int,2>> GroupSingleExon; // single exon的reads;
    std::unordered_map<std::string, int> GroupReadStrands;
};

GroupInformation knowGroupInformation(std::streampos& startpos, std::streampos& endpos, const std::string& sam_file_path, const int& Sj_Support_Read_Number) {
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



struct SE_belong2_genetranscript {
    std::unordered_map<int, std::unordered_map<std::string, double>> file_SE_reads_gene_number;
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
                                    std::ofstream& Trace) {
    SE_belong2_genetranscript SingleReads;
    const auto& group_genes = thisGroupAnnotation.group_genes;
    const auto& group_me_transcripts = thisGroupAnnotation.group_me_transcripts;
    const auto& group_se_transcripts = thisGroupAnnotation.group_se_transcripts;

    // gene-reads;
    std::unordered_map<std::string, std::vector<std::string>> SingleReads_inGenes;
    SingleReads_inGenes.reserve(group_genes.size());

    for (const auto& EachGene:group_genes) {
        SingleReads_inGenes.emplace(EachGene.first, std::vector<std::string>{});
    }

    std::unordered_map<std::string, std::vector<std::string>> SingleReads_inTranscripts;
    SingleReads_inTranscripts.reserve(group_me_transcripts.size() + group_se_transcripts.size());

    for (const auto& EachTx:group_me_transcripts) {
        SingleReads_inTranscripts.emplace(EachTx.first, std::vector<std::string>{});
    }
    for (const auto& EachTx:group_se_transcripts) {
        SingleReads_inTranscripts.emplace(EachTx.first, std::vector<std::string>{});
    }

    // 生成一个gene:{exons}的unordered_map的列表;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> genes2exons;
    genes2exons.reserve(group_genes.size());
    
    // 遍历所有基因
    for (const auto& g : group_genes) {
        auto it_gene = gene2tx.find(g.first);
        if (it_gene == gene2tx.end()) continue;
        std::vector<std::array<int,2>> exons;
        for (const std::string& tx : it_gene->second) {
            auto it_tx = tx2exon.find(tx);
            if (it_tx == tx2exon.end()) continue;
            exons.insert(exons.end(), it_tx->second.begin(), it_tx->second.end());
        }
        // 赋值到结果map
        exons = mergeIntervals(exons);
        genes2exons.emplace(g.first, std::move(exons));
    }    

    std::vector<std::string> Temp_Gene; Temp_Gene.reserve(32);
    std::vector<std::array<int,2>> tmpReadVec(1);
    std::vector<std::string> Part_Transcript;

    std::vector<int> Tx_Temp_Dist; 
    std::vector<int> Tx_Temp_exon; 
    std::vector<int> Tx_Temp_Interval;
    std::vector<int> edgeDist;
    
    for (const auto& EachSing:GroupSingleExonReads) {

        const auto& read_name = EachSing.first;
        const auto& read_exon = EachSing.second;
        Temp_Gene.clear();
        tmpReadVec[0] = read_exon;

        for (const auto& EachGene:group_genes) {
            if ( read_exon[0] > EachGene.second[1] + Edge || read_exon[1] < EachGene.second[0] - Edge ) { continue; }
            if (IntervalIntersection(genes2exons[EachGene.first], tmpReadVec) > 0) { Temp_Gene.push_back(EachGene.first); }            
        }

        int changdu = 0;
        if (Temp_Gene.size() == 0) { continue; } 
        
        Tx_Temp_Dist.clear();
        Tx_Temp_exon.clear();
        Part_Transcript.clear();
        Tx_Temp_Interval.clear();
        
        if (Temp_Gene.size() == 1) {
            const auto& gene = Temp_Gene[0];
            SingleReads_inGenes[gene].push_back(read_name);

            // 从这里继续分这个基因的reads给isoform;
            const auto& Temp_Transcript = gene2tx[gene];

            if (Temp_Transcript.size() == 1) {
                // 只有一个转录本;
                SingleReads_inTranscripts[Temp_Transcript[0]].push_back(read_name);
            } else {
                // 这一个基因有多个转录本; 首先是交集筛选;
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
                if (Part_Transcript.size() == 0) { continue; } 
                
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

        } else {
            // 如果有多个基因, 看与哪个基因重复的更多;
            edgeDist.clear();
            edgeDist.reserve(Temp_Gene.size());

            std::vector<int> Temp_Dist;
            std::vector<int> Temp_exon;
            Temp_Dist.reserve(Temp_Gene.size());
            Temp_exon.reserve(Temp_Gene.size());

            for (const auto& EachGene:Temp_Gene) {
                int MinDist = INT_MAX;
                const auto& exons = genes2exons[EachGene];
                int intersection = IntervalIntersection(exons, tmpReadVec);
                Temp_Dist.push_back(intersection);
                Temp_exon.push_back(exons.size());
                for (const auto& EachTx:gene2tx[EachGene]) {
                    MinDist = std::min(MinDist, IntervalMinDistance(tx2exon[EachTx], read_exon));
                }
                edgeDist.push_back(MinDist);
                // std::cerr << EachGene << " intersection:" << intersection << " exonNumber:" << SJ_Vector.size() << " MinDist:" << MinDist << std::endl;
            }

            // 首先看看有没有单个exon的基因; 但是有多个单exon的基因呢?
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

            if (MaxGeneName.size() != 0) {
                SingleReads_inGenes[MaxGeneName].push_back(read_name);
                // 开始转录本;
                const auto& Temp_Transcript = gene2tx[MaxGeneName];

                if (Temp_Transcript.size() == 1) {
                    // 只有一个转录本;
                    SingleReads_inTranscripts[Temp_Transcript[0]].push_back(read_name);
                } else {
                    // 这一个基因有多个转录本; 首先是交集筛选;
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
                    if (Part_Transcript.size() == 0) { continue; } 
                    
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
    genes2exons.clear();

    std::unordered_map<int, std::unordered_map<std::string, double>> FileSingleGeneNumber;
    FileSingleGeneNumber.reserve(FileNo);
    for (int i = 0; i < FileNo; i++) {
        FileSingleGeneNumber.emplace(i, std::unordered_map<std::string, double>{});
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
                        FileSingleGeneNumber[EachFile.first][EachGene.first] = static_cast<double>(EachFile.second.size()); 
                    }
                }
            }
        }

    } else {
        if (SingleReads_inGenes.size() > 0) {
            for (const auto& EachGene:SingleReads_inGenes) {   
                FileSingleGeneNumber[0][EachGene.first] = static_cast<double>(EachGene.second.size());
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


void write_single_exon_gtf_trace(int& FileNo,
    std::unordered_map<std::string, std::string>& groupreadfiles,
    std::unordered_map<std::string, std::vector<std::string>>& singleexonwithreads,
    std::unordered_map<std::string, std::array<int,2>>& singleexongroupannotation,
    std::ofstream& Updated_Files, std::ofstream& isoformFilePath,
    std::unordered_map<std::string, std::string>& GTF_Transcript_Strand,
    std::string& chrname) {

    std::unordered_map<std::string, std::vector<double>> FileSingleExonNumber;
    FileSingleExonNumber.reserve(singleexonwithreads.size());

    std::unordered_map<std::string, std::vector<std::string>> remainReads;
    remainReads.reserve(singleexonwithreads.size());

    std::string updateBuffer;
    std::string isoformBuffer;

    for (const auto& eachAnno:singleexonwithreads) {
        const std::string& transcript = eachAnno.first;
        auto itAnno = singleexongroupannotation.find(transcript);

        if (itAnno == singleexongroupannotation.end()) {
            if (!eachAnno.second.empty()) remainReads.emplace(transcript, eachAnno.second);
            continue;
        }
        if (eachAnno.second.empty()) continue;
        const auto& exonRange = itAnno->second;

        const auto strandIt = GTF_Transcript_Strand.find(transcript);
        const std::string strand = (strandIt != GTF_Transcript_Strand.end()) ? strandIt->second : ".";

        const size_t pos = eachAnno.first.find('|');
        std::string first_part, second_part;
        if (pos != std::string::npos) {
            first_part  = transcript.substr(0, pos);
            second_part = transcript.substr(pos + 1);
        } else {
            first_part = second_part = transcript;
        }

        updateBuffer.append(chrname); updateBuffer.push_back('\t');
        updateBuffer.append("annotated_isoform\ttranscript\t");
        updateBuffer.append(std::to_string(exonRange[0])); updateBuffer.push_back('\t');
        updateBuffer.append(std::to_string(exonRange[1])); updateBuffer.append("\t.\t");
        updateBuffer.append(strand); updateBuffer.append("\t.\tgene_id \"");
        updateBuffer.append(first_part); updateBuffer.append("\"; transcript_id \"");
        updateBuffer.append(second_part); updateBuffer.append("\";\n");

        updateBuffer.append(chrname); updateBuffer.push_back('\t');
        updateBuffer.append("annotated_isoform\texon\t");
        updateBuffer.append(std::to_string(exonRange[0])); updateBuffer.push_back('\t');
        updateBuffer.append(std::to_string(exonRange[1])); updateBuffer.append("\t.\t");
        updateBuffer.append(strand); updateBuffer.append("\t.\tgene_id \"");
        updateBuffer.append(first_part); 
        updateBuffer.append("\"; transcript_id \"");
        updateBuffer.append(second_part);
        updateBuffer.append("\"; exon_number \"0\";\n");                

        std::vector<double> fileCounts(FileNo, 0.0);

        for (const auto& read : eachAnno.second) {
            auto itGF = groupreadfiles.find(read);
            const std::string& fileID = itGF->second;
            int idx = std::stoi(fileID); // FileNo 由上层保证正确
            fileCounts[idx] += 1;
        }

        FileSingleExonNumber.emplace(transcript, std::move(fileCounts));
    }

    // 批量输出文件（仅上锁一次）
    if (Updated_Files.is_open()) {
        std::unique_lock<std::mutex> lock(updatedGtfMutex);
        Updated_Files << updateBuffer;
    }   
    
    if (isoformFilePath.is_open()) {
        for (const auto& kv : FileSingleExonNumber) {
            const std::string& transcript = kv.first;
            const std::vector<double>& counts = kv.second;
            size_t pos = transcript.find('|');
            std::string first_part, second_part;
            if (pos != std::string::npos) {
                first_part  = transcript.substr(0, pos);
                second_part = transcript.substr(pos + 1);
            } else {
                first_part = second_part = transcript;
            }
            isoformBuffer.append(second_part);
            isoformBuffer.push_back('\t');
            isoformBuffer.append(first_part);
            for (double x : counts) {
                isoformBuffer.push_back('\t');
                isoformBuffer.append(std::to_string(static_cast<long long>(x)));
            }
            isoformBuffer.push_back('\n');
        }
        {
            std::unique_lock<std::mutex> lock(isoformCountMutex);
            isoformFilePath << isoformBuffer;
        }
    }    

    singleexonwithreads = std::move(remainReads);
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

    std::vector<std::array<int,2>> vec;
    std::string vecname;
    std::size_t cluster = 0;
    std::size_t clusterNumber = 0;
    for (const auto& pair : groupreadsjs) {
        vecname = pair.first;
        vec = pair.second;
        auto itsj = readsj2sizet.find(vec);
        if (itsj != readsj2sizet.end()) {
            cluster = readsj2sizet[vec];
            ClustersReads[cluster].push_back(vecname);
        } else {
            // 新的cluster;
            readsj2sizet[vec] = clusterNumber;
            ClustersReads[clusterNumber].push_back(vecname);
            clusterNumber = clusterNumber + 1;
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


ReferenceCluster get_FSM_and_others(std::unordered_map<std::size_t, std::vector<std::string>>& ClustersReads, 
                                   std::unordered_map<std::string, std::vector<std::array<int,2>>>& GroupAnno, 
                                   std::unordered_map<std::string, std::array<int,2>>& GroupAnnoBE, 
                                   std::unordered_map<std::string, std::array<int,2>>& AllReadBE, 
                                   std::unordered_map<std::string, std::vector<std::array<int,2>>>& AllSJs,
                                   std::unordered_map<std::string, std::string>& GroupReadFile, 
                                   std::ofstream& Trace) {
    ReferenceCluster FSMISMO = {};

    std::vector<std::array<int,2>> each_read_SJs; 
    std::vector<std::string> statis_FSM; 
    std::vector<std::string> statis_ISM; 
    std::vector<int> ISM_flag; 
    ISM_flag.reserve(30);

    std::string FSM_name; 
    std::vector<int> minFSM; 
    minFSM.reserve(30); 

    std::string first_part;
    std::string second_part; 

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
                        if (it_pos == ref_sj.end()) {
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
            if (pos != std::string::npos) {
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
                    buffer += GroupReadFile[EachRead];
                    buffer += '\n';
                }
                {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << buffer;
                }
            }

        } else if (!statis_FSM.empty()) {
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
                int min_index = std::min_element(minFSM.begin(), minFSM.end()) - minFSM.begin();

                FSM_name = statis_FSM[min_index];
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
            } else{
                //这种的情况是ISM类, 不用管是哪种ISM; //如果后续改识别ISM类型的时候, 更改这里; ★★★★★★★
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
    //所有的需要进一步进行识别的Cluster, 以及他们对应的数量;
    std::map<size_t, std::vector<std::string>> HighConClusters;
    //淘汰的Cluster及其他们对应的数量;
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
    
    //中间用到的变量;
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
            if (SJs2.size()>SJs1.size()) {
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
                // if (RecycleLength < 20*LowReadSJVec.size() and RecycleLength > 0) {
                if (RecycleLength < 60 and RecycleLength > 0) {    
                    potential_FSM[each_anno.first] = RecycleLength;
                }
            }
        }

        if (potential_FSM.size() == 1){

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
                    buffer += EachRead;
                    buffer += '\t';
                    buffer += "simFSM\t";
                    buffer += second_part;
                    buffer += '\t';
                    buffer += first_part;
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
            // 当前已经有全长的;
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
                    buffer += EachRead;
                    buffer += '\t';
                    buffer += "simFSM\t";
                    buffer += second_part;
                    buffer += '\t';
                    buffer += first_part;
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

    for (const auto& newFSM:RecycleAnnoName_value){
        newName.clear();
        for (const auto& allCluster:newFSM.second){ 
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
                                std::ofstream& traceFilePath,
                                const int& Sj_Support_Number, std::string& CHR) {
    SpliceChainClass SCC;
    ReferenceCluster FsmIsmOthers;
    HighLowClusters Others2HighLow;
    std::unordered_map<std::size_t, std::vector<std::string>> groupCluster = classifyReadsVec(groupreadsjs);

    SCC.ClusterCoverage = get_every_cluster_begin_end(groupCluster, groupreadcoverage);
    FsmIsmOthers = get_FSM_and_others(groupCluster, groupannotations, AnnoCoverage, groupreadcoverage, groupreadsjs, groupreadfiles, traceFilePath); //readTraceMutex;
    Others2HighLow = get_HighLow_clusters(FsmIsmOthers.Others, groupreadsjs, groupreadsigns, Sj_Support_Number);  

    get_filtered_FSM(Others2HighLow.LowConClusters, groupannotations, FsmIsmOthers.FSM, groupreadsjs, groupreadfiles, traceFilePath); //readTraceMutex;     

    // 将有用结果输出;
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

    int Merge_L;
    int Inter_L;
    int Nonoverlap_L;

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
                    Merge_L = IntervalMerge(EachAnno.second, SJVec1);
                    Inter_L = IntervalIntersection(EachAnno.second, SJVec1);
                    Nonoverlap_L = findNonOverlappingIntervals(EachAnno.second, SJVec1);
                    LineElement = Merge_L - Inter_L - Nonoverlap_L;
                    // LineElement = IntervalMerge(EachAnno.second, SJVec1) - IntervalIntersection(EachAnno.second, SJVec1);
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
                Merge_L = IntervalMerge(SJVec1, SJVec2);
                Inter_L = IntervalIntersection(SJVec1, SJVec2);
                Nonoverlap_L = findNonOverlappingIntervals(SJVec1, SJVec2);
                LineElement = Merge_L - Inter_L - Nonoverlap_L;
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

    if (adj.size() == 0){
        return CliqueVec;
    } else{
        //adj已经有了;
        for (const auto& eachnode:adj) cand.insert(eachnode.first);
        subg = cand;
        Q.push_back(-1);
        //首先找到邻居节点最多的节点;        
        max_node = find_max_neighborhood_node(subg, adj);

        std::set_difference(cand.begin(), cand.end(),
                            adj[max_node].begin(), adj[max_node].end(),
                            std::inserter(ext_u, ext_u.begin()));

        while (true)
        {
            if (!ext_u.empty()) 
            {
                q = *ext_u.begin(); //获取集合中的第一个元素
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

    std::set<int> subgraph_known_true_nodes;
    std::set<int> subgraph_remain_nodes;
    int standard_read_number = std::numeric_limits<int>::max();
    int known_true_number = 0;
    int max_node_value = 0;

    for (const auto& Anno : Number2Anno){
        Node_Results.NodesKnownTrue.insert(Anno.first);
    }

    if (CqVec.size() != 0){

        if (CqVec.size() > 1) {
            std::sort(CqVec.begin(), CqVec.end(), compareBySize);
        }

        for (const auto& eachClique : CqVec) {
            standard_read_number = std::numeric_limits<int>::max();

            std::set<int> CliqueSet(eachClique.begin(), eachClique.end());
            subgraph_known_true_nodes.clear();
            std::set_intersection(CliqueSet.begin(), CliqueSet.end(), Node_Results.NodesKnownTrue.begin(), Node_Results.NodesKnownTrue.end(),
                          std::inserter(subgraph_known_true_nodes, subgraph_known_true_nodes.begin()));
            
            subgraph_remain_nodes.clear();
            std::set_difference(CliqueSet.begin(), CliqueSet.end(), subgraph_known_true_nodes.begin(), subgraph_known_true_nodes.end(),
                        std::inserter(subgraph_remain_nodes, subgraph_remain_nodes.begin()));

            if (subgraph_known_true_nodes.size() != 0) {
                if (subgraph_known_true_nodes.size() == 1) {
                    standard_read_number = Number2Count[(*subgraph_known_true_nodes.begin())];
                } else {

                    for(const auto& each_True_node : subgraph_known_true_nodes) {

                        if (Number2Count[each_True_node] < standard_read_number) {
                            standard_read_number = Number2Count[each_True_node];
                        }
                    }
                }

                if (subgraph_remain_nodes.size() != 0) { 

                    for (const auto& eachknownnode:subgraph_known_true_nodes) {
                        Node_Results.TrueNodeSet.insert(subgraph_known_true_nodes.begin(), subgraph_known_true_nodes.end());
                    }

                    for (const auto& each_remain_node : subgraph_remain_nodes){
                        if (standard_read_number == 0) {
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
                    for (const auto& eachknownnode:subgraph_known_true_nodes) {
                        // 十分关键; 对于速度; ☆★☆;☆★☆;☆★☆;☆★☆;☆★☆
                        if (Number2Count[eachknownnode] > 0) {
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
            if (iter != Node_Results.TrueNodeSet.end()){
                Node_Results.TrueNodeSet.erase(eachFalse);
            }
        }
    }
    return Node_Results;
}


struct Solvent {
    std::map<std::string, std::unordered_map<std::string, int>> File_FSM;
    std::map<std::string, std::map<std::size_t, std::vector<std::string>>> File_ISM;
    std::map<std::string, std::map<std::size_t, std::vector<std::string>>> File_HighConClusters;
    std::map<std::string, std::map<std::size_t, std::vector<std::string>>> File_LowConClusters;
};

Solvent get_Solvent_FsmIsmHigh(SpliceChainClass& FsmIsmHigh, int& FileNumber, GroupInformation& groupinformations) {
    Solvent FileSpliceChains;
    std::unordered_map<std::string, int> thisTranscript_File2Count; // FSM的临时变量;
    std::map<std::string, std::vector<std::string>> thisTranscript_File2Count_ISM; // ISM和High的临时变量;
    std::string Fi;

    if (FsmIsmHigh.FSM.size() > 0) {

        for (int i = 0; i < FileNumber; i++) {
            Fi = std::to_string(i);
            FileSpliceChains.File_FSM[Fi] = {};
        }        

        for (const auto& eachCluster:FsmIsmHigh.FSM) {
            thisTranscript_File2Count.clear();
            for (int i = 0; i < FileNumber; i++) {
                Fi = std::to_string(i);
                thisTranscript_File2Count[Fi] = 0;
            }
            std::string AnnoName = eachCluster.first;
 
            for (const auto& eachRead:eachCluster.second.second) {

                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                thisTranscript_File2Count[read_file] = thisTranscript_File2Count[read_file] + 1;
            }

            for (const auto& eachFile:thisTranscript_File2Count) {
                if (eachFile.second > 0) {
                    FileSpliceChains.File_FSM[eachFile.first][AnnoName] = eachFile.second;
                }
            }            
        }
    } // 结束;

    if (FsmIsmHigh.ISM.size() > 0) {
        for (int i = 0; i < FileNumber; i++) {
            FileSpliceChains.File_ISM[std::to_string(i)] = {};
        }
        for (const auto& eachCluster:FsmIsmHigh.ISM) {
            thisTranscript_File2Count_ISM.clear();
            for (int i = 0; i < FileNumber; i++) {
                thisTranscript_File2Count_ISM[std::to_string(i)] = {};
            }
            std::size_t NonName = eachCluster.first;
            // 对这个ISM所属的每条reads进行循环;
            for (const auto& eachRead:eachCluster.second) {
                // 这条read所在的文件;
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                thisTranscript_File2Count_ISM[read_file].push_back(eachRead);
            }
            // 所有的reads循环结束;
            // 结果填写到ISM中;
            for (const auto& eachFile:thisTranscript_File2Count_ISM) {
                if (eachFile.second.size() > 0) {
                    FileSpliceChains.File_ISM[eachFile.first][NonName] = eachFile.second;
                }
            }
        }
    }

    // 拆分HighConfidenceClusters;
    if (FsmIsmHigh.HighConClusters.size() > 0) {
        // 首先初始化;
        for (int i = 0; i < FileNumber; i++) {
            FileSpliceChains.File_HighConClusters[std::to_string(i)] = {};
        }
        // 对每个HighConCluster循环;
        for (const auto& eachCluster:FsmIsmHigh.HighConClusters) {
            thisTranscript_File2Count_ISM.clear();
            for (int i = 0; i < FileNumber; i++) {
                thisTranscript_File2Count_ISM[std::to_string(i)] = {};
            }
            std::size_t NonName = eachCluster.first;
            // 对这个HighConCluster所属的每条reads进行循环;
            for (const auto& eachRead:eachCluster.second) {
                // 这条read所在的文件;
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                thisTranscript_File2Count_ISM[read_file].push_back(eachRead);
            }
            // 所有的reads循环结束;
            // 结果填写到HighConCluster中;
            for (const auto& eachFile:thisTranscript_File2Count_ISM) {
                if (eachFile.second.size() > 0) {
                    FileSpliceChains.File_HighConClusters[eachFile.first][NonName] = eachFile.second;
                }
            }             
        }
    }

    // 拆分LowConfidenceClusters;
    if (FsmIsmHigh.LowConClusters.size() > 0) {
        // 首先初始化;
        for (int i = 0; i < FileNumber; i++) {
            FileSpliceChains.File_LowConClusters[std::to_string(i)] = {};
        }
        // 对每个LowConCluster循环;
        for (const auto& eachCluster:FsmIsmHigh.LowConClusters) {
            thisTranscript_File2Count_ISM.clear();
            for (int i = 0; i < FileNumber; i++) {
                thisTranscript_File2Count_ISM[std::to_string(i)] = {};
            }
            std::size_t NonName = eachCluster.first;
            // 对这个LowConCluster所属的每条reads进行循环;            
            for (const auto& eachRead:eachCluster.second) {
                // 这条read所在的文件;
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                thisTranscript_File2Count_ISM[read_file].push_back(eachRead);                
            }
            // 所有的reads循环结束;
            // 结果填写到LowConCluster中;
            for (const auto& eachFile:thisTranscript_File2Count_ISM) {
                if (eachFile.second.size() > 0) {
                    FileSpliceChains.File_LowConClusters[eachFile.first][NonName] = eachFile.second;
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

    std::vector<std::string> TempGene;
    std::vector<int> TempDist;
    std::set<std::array<int,2>> SJ_set; 
    std::vector<std::string> ALL_tx;
    // std::set<std::array<int, 2>> intersection;
    int intersection{0};

    // 对每个基因循环;
    for (const auto& A_Gene:gene_set) {
        SJ_set.clear();
        // intersection.clear();
        TempGene.push_back(A_Gene.first);
        ALL_tx = gene2tx[A_Gene.first];
        for (const auto& A_Tx:ALL_tx) {
            // 这个转录本在里面;
            if (tx2sj.find(A_Tx) != tx2sj.end()) {
                SJ_set.insert(tx2sj[A_Tx].begin(), tx2sj[A_Tx].end());
            }
        }
        std::vector<std::array<int, 2>> SJ_vec(SJ_set.begin(), SJ_set.end());
        intersection = IntervalIntersection(SJ_vec, novel_isoform);
        TempDist.push_back(intersection);
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
    std::unordered_map<std::string, std::unordered_map<std::string, double>> File_TranscriptNumber;
};

//18. 将识别的结果输出到txt;
OutputInformation Write_Detection_Transcript2gtf(std::ofstream& Updated_Files, std::ofstream& Trace,
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
                                                std::unordered_map<std::string, std::vector<std::string>>& gtf_gene2tx) {
    // std::mutex& GtfMutex, std::mutex& TraceMutex;
    OutputInformation FinalAnnotations;
    std::string itsname;
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

    if (Updated_Files.is_open()){
        // 第一次
        for (const auto aaa:noderesults.TrueNodeSet) {
            if (aaa < Disinform.Index2Anno.size()){
                //有注释的部分;
                itsname = Disinform.Index2Anno[aaa];
                itssj = groupannotations[itsname];
                itsexon = Annoexon[itsname];
                FinalAnnotations.Transcript_Annotations[itsname].first = itssj;
                FinalAnnotations.Transcript_Annotations[itsname].second = Disinform.Index2Count[aaa];
                
                size_t pos = itsname.find('|');
                if (pos != std::string::npos) {
                    // 根据 "|" 将字符串分割为两部分
                    first_part = itsname.substr(0, pos);
                    second_part = itsname.substr(pos + 1);
                }
                {
                    std::unique_lock<std::mutex> lock(updatedGtfMutex);
                    Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "transcript" << '\t' << itsexon[0][0] << '\t' << itsexon[itsexon.size()-1][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[itsname] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\";" << '\n';
                    for (int i = 0; i < itsexon.size(); i++){
                        Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "exon" << '\t' << itsexon[i][0] << '\t' << itsexon[i][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[itsname] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\"; exon_number \"" << i << "\";" << '\n';
                    }
                }
                FinalAnnotations.transcript2gene[itsname] = first_part;

            } else {
                novel_count = novel_count + 1;
                //没有注释的部分;
                itsname = chrname + "-novel-" + group_size + "-" + std::to_string(aaa) + "-" + std::to_string(novel_count);
                Disinform.Index2novelname[aaa] = itsname;
                itssj = Disinform.Index2Unknown[aaa];
                FinalAnnotations.Transcript_Annotations[itsname].first = itssj;
                FinalAnnotations.Transcript_Annotations[itsname].second = Disinform.Index2Count[aaa];
                nameya = Disinform.Index2hashname[aaa];
                std::array<int,2> exonBE = groupreadcoverage[nameya];
                BE[0] = itssj.front()[0];
                BE[1] = itssj.back()[1];

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
                            if (BE[0] > Annogenecovergae[every_gene.first][0] and BE[1] < Annogenecovergae[every_gene.first][1]) {
                                TempGene.push_back(every_gene.first);
                            }
                        }
                        if (TempGene.size() == 0) {
                            first_part = "NA";
                        } else if (TempGene.size() == 1) {
                            if (BE[0] > Annogenecovergae[TempGene[0]][0] and BE[1] < Annogenecovergae[TempGene[0]][1]) {
                                first_part = TempGene[0];
                            } else {
                                first_part = "NA";
                            }
                        } else {
                            TempDist.clear();
                            for (const auto& each_gene:TempGene){
                                int a = abs(exonBE[1] - Annogenecovergae[each_gene][1]);
                                int b = abs(Annogenecovergae[each_gene][0] - exonBE[0]);
                                TempDist.push_back((a<b)?a:b);
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
                    std::unique_lock<std::mutex> lock(updatedGtfMutex);
                    Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "transcript" << '\t' << BE[0] << '\t' << BE[1] << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\";" << '\n';
                    for (int i = 0; i < itssj.size()+1; i++) {
                        if (i == 0) {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << BE[0] << '\t' << itssj[i][0]-1 << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" <<  '\n';
                        } else if (i == itssj.size()) {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << itssj[i-1][1]+1 << '\t' << BE[1] << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" << '\n';
                        } else {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << itssj[i-1][1]+1 << '\t' << itssj[i][0]-1 << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" << '\n';
                        }
                    }
                }
                
                if (Trace.is_open()) {
                    const auto& read_list = HighClusters[nameya];
                    std::string buffer;
                    buffer.reserve(read_list.size() * 80);              
                    for (const auto& EachRead : read_list) {
                        const std::string& file = groupreadfiles[EachRead];
                        buffer += EachRead;
                        buffer += '\t';
                        buffer += "novel\t";
                        buffer += itsname;
                        buffer += '\t';
                        buffer += first_part;
                        buffer += '\t';
                        buffer += file;
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

OutputInformation Write_Detection_Transcript2gtf_MultiFiles(std::ofstream& Updated_Files, std::ofstream& Trace,
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
                                                int& FileNumber, std::unordered_map<std::string, std::string>& groupreadfiles,
                                                Solvent& filesolvent,
                                                std::unordered_map<std::string, std::vector<std::string>>& gtf_gene2tx){
    OutputInformation FinalAnnotations;
    std::string itsname;
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

    for (int i = 0; i < FileNumber; i++) {
        FinalAnnotations.File_TranscriptNumber[std::to_string(i)] = {};
    }

    if (Updated_Files.is_open()) {
        // 第一次
        for (const auto aaa:noderesults.TrueNodeSet){
            if (aaa < Disinform.Index2Anno.size()){
                //有注释的部分;
                itsname = Disinform.Index2Anno[aaa];
                itssj = groupannotations[itsname];
                itsexon = Annoexon[itsname];
                FinalAnnotations.Transcript_Annotations[itsname].first = itssj;
                FinalAnnotations.Transcript_Annotations[itsname].second = Disinform.Index2Count[aaa];

                for (int i = 0; i < FileNumber; i++) {
                    auto is1 = filesolvent.File_FSM[std::to_string(i)].find(itsname);
                    if (is1 != filesolvent.File_FSM[std::to_string(i)].end()) {
                        FinalAnnotations.File_TranscriptNumber[std::to_string(i)][itsname] = filesolvent.File_FSM[std::to_string(i)][itsname];
                    } else {
                        FinalAnnotations.File_TranscriptNumber[std::to_string(i)][itsname] = 0;
                    }
                }

                size_t pos = itsname.find('|');
                if (pos != std::string::npos) {
                    // 根据 "|" 将字符串分割为两部分
                    first_part = itsname.substr(0, pos);
                    second_part = itsname.substr(pos + 1);
                }
                {
                    std::unique_lock<std::mutex> lock(updatedGtfMutex);
                    Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "transcript" << '\t' << itsexon[0][0] << '\t' << itsexon[itsexon.size()-1][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[itsname] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\";" << '\n';
                    for (int i = 0; i < itsexon.size(); i++){
                        Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "exon" << '\t' << itsexon[i][0] << '\t' << itsexon[i][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[itsname] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\"; exon_number \"" << i << "\";" << '\n';
                    }
                }
                FinalAnnotations.transcript2gene[itsname] = first_part;

            } else {
                // 这个是novel了;
                novel_count = novel_count + 1;
                //没有注释的部分;
                itsname = chrname + "-novel-" + group_size + "-" + std::to_string(aaa) + "-" + std::to_string(novel_count);
                Disinform.Index2novelname[aaa] = itsname;
                itssj = Disinform.Index2Unknown[aaa];
                FinalAnnotations.Transcript_Annotations[itsname].first = itssj;
                FinalAnnotations.Transcript_Annotations[itsname].second = Disinform.Index2Count[aaa];
                nameya = Disinform.Index2hashname[aaa];
                std::array<int,2> exonBE = groupreadcoverage[nameya];
                BE[0] = itssj.front()[0];
                BE[1] = itssj.back()[1];
                
                for (int i = 0; i < FileNumber; i++) {
                    auto is1 = filesolvent.File_HighConClusters[std::to_string(i)].find(nameya);
                    if (is1 != filesolvent.File_HighConClusters[std::to_string(i)].end()) {
                        FinalAnnotations.File_TranscriptNumber[std::to_string(i)][itsname] = filesolvent.File_HighConClusters[std::to_string(i)][nameya].size();
                    } else {
                        FinalAnnotations.File_TranscriptNumber[std::to_string(i)][itsname] = 0;
                    }
                }
                
                TempGene.clear();
                if (groupallgene.size() != 0) {
                    TempGeneName = get_co_sj_gene(itssj, groupallgene, gtf_gene2tx, groupannotations);

                    if (TempGeneName.size() != 0) {
                        if (BE[0] > Annogenecovergae[TempGeneName][0] and BE[1] < Annogenecovergae[TempGeneName][1])  {
                            first_part = TempGeneName;
                        } else {
                            first_part = "NA";
                        }
                    } else {
                        for (const auto& every_gene:groupallgene) {
                            if ( BE[0] > Annogenecovergae[every_gene.first][0] and BE[1] < Annogenecovergae[every_gene.first][1] ) {
                                TempGene.push_back(every_gene.first);
                            }
                        }
                        if (TempGene.size() == 0) {
                            first_part = "NA";
                        } else if (TempGene.size() == 1) {
                            if ( BE[0] > Annogenecovergae[TempGene[0]][0] and BE[1] < Annogenecovergae[TempGene[0]][1] ) {
                                first_part = TempGene[0];
                            } else {
                                first_part = "NA";
                            }
                        } else {
                            TempDist.clear();
                            for (const auto& each_gene:TempGene) {
                                int a = abs(exonBE[1] - Annogenecovergae[each_gene][1]);
                                int b = abs(Annogenecovergae[each_gene][0] - exonBE[0]);
                                TempDist.push_back((a<b)?a:b);
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
                    std::unique_lock<std::mutex> lock(updatedGtfMutex);
                    Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "transcript" << '\t' << BE[0] << '\t' << BE[1] << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\";" << '\n';
                    for (int i = 0; i < itssj.size()+1; i++){
                        if (i == 0) {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << BE[0] << '\t' << itssj[i][0]-1 << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" <<  '\n';
                        } else if (i == itssj.size()) {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << itssj[i-1][1]+1 << '\t' << BE[1] << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" << '\n';
                        } else {
                            Updated_Files << chrname << '\t' << "novel_isoform" << '\t' << "exon" << '\t' << itssj[i-1][1]+1 << '\t' << itssj[i][0]-1 << '\t' << "." << '\t' << High_Strand[nameya] << '\t' << "." << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << itsname << "\"; exon_number \"" << i << "\";" << '\n';
                        }
                    }
                }
                //输出追踪文件;
                if (Trace.is_open()){
                    const auto& read_list = HighClusters[nameya];
                    std::string buffer;
                    buffer.reserve(read_list.size() * 80);                
                    for (const auto& EachRead : read_list) {
                        const std::string& file = groupreadfiles[EachRead];
                        buffer += EachRead;
                        buffer += '\t';
                        buffer += "novel\t";
                        buffer += itsname;
                        buffer += '\t';
                        buffer += first_part;
                        buffer += '\t';
                        buffer += file;
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


struct FileFalseNode
{
    std::set<int> ThisFlaseNode;
    std::unordered_map<int, std::vector<std::string>> FasleNode2Count;
};
FileFalseNode get_File_False_node (std::set<int>& FalseNodeSet, 
                          std::unordered_map<int, std::size_t>& Node2hashname, 
                          std::map<std::string, std::map<size_t, std::vector<std::string>>>& File_HighClusters, 
                          int & FileNNumber) {
    //这个Barcode中所有错误的节点;
    FileFalseNode ThisFile_FalseNumber;
    std::map<size_t, std::vector<std::string>> this_file_HC = File_HighClusters[std::to_string(FileNNumber)];
    for (const auto& node:FalseNodeSet) {
        //每个错误节点不一定在这个File的HighCluster中;
        auto which = this_file_HC.find(Node2hashname[node]);
        if (which != this_file_HC.end()) {
            ThisFile_FalseNumber.ThisFlaseNode.insert(node);
            ThisFile_FalseNumber.FasleNode2Count[node] = this_file_HC[Node2hashname[node]];
        }
    }
    return ThisFile_FalseNumber;
}


int whether_isoform_part_is_not(std::vector<std::array<int,2>> AnnoSJ, std::vector<std::array<int,2>> ISMSJ){
    int ISM_Flag = 0;
    auto it1 = std::find(AnnoSJ.begin(), AnnoSJ.end(), ISMSJ[0]);
    if (it1 != AnnoSJ.end()){
        //计算出it1的索引;
        int it1_index = it1 - AnnoSJ.begin();
        auto it2 = std::find(AnnoSJ.begin(), AnnoSJ.end(), ISMSJ[ISMSJ.size()-1]);
        if (it2 != AnnoSJ.end()){
            //计算出it2的索引;
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
        return ""; // 如果集合为空，返回空字符串
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


IndicateFire Quantification_initialization(std::map<std::size_t, std::vector<std::string>>& groupISM, 
                                           OutputInformation& FinallyAnnotations, 
                                           std::set<int>& truenodeset,
                                           std::set<int>& falsenodeset, 
                                           DistanceInform& Disinform, 
                                           std::string& group_size, 
                                           std::map<size_t, std::vector<std::string>>& HighClusters, 
                                           std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs, 
                                           std::unordered_map<std::string, std::string>& groupreadfiles,
                                           std::ofstream& Trace,
                                           std::string& chrname) {
    IndicateFire OutputResults;
    std::vector<std::string> Order_Transcript_Name;
    Eigen::MatrixXd known_ISM_matrix;
    std::vector<std::array<int,2>> AreadSJs;
    
    int RowIndex = -1;
    int ColIndex = -1;
    int flag = 0;

    for (const auto& eachTransctipt:FinallyAnnotations.Transcript_Annotations){
        Order_Transcript_Name.push_back(eachTransctipt.first);
    }

    Eigen::RowVectorXd EachClusterRowVector(Order_Transcript_Name.size());
    Eigen::VectorXd known_ISM_cluster_numbers;

    std::vector<std::string> ism_gtf_name;
    std::set<std::string> first_part_set;
    std::set<std::string> second_part_set;
    std::string first_part;
    std::string second_part;
    std::vector<int> ColIndexVec;
    double max_ratio = 0;
    double this_ratio = 0;

    if (groupISM.size() != 0) {
        RowIndex = -1;
        for (const auto& eachISM:groupISM) {
            //对生成的每个行向量清空; 每个cluster的行向量;
            AreadSJs = groupreadsjs[eachISM.second[0]];
            EachClusterRowVector.setZero();
            ColIndex = -1;
            ism_gtf_name.clear();
            first_part_set.clear();
            second_part_set.clear();
            ColIndexVec.clear();   
            double max_ratio = 0;
            double this_ratio = 0;
            //对每个注释遍历;
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
                } else{
                    EachClusterRowVector(ColIndex) = 0;
                }
            }             

            // 在这里加一步过滤最好了, ISM如果有一个特别近的, 还有一个特别远的, 特别远的就去掉吧;
            if (ColIndexVec.size() > 1) {
                if (max_ratio > 0.4) {
                    for (const auto& ColIndex:ColIndexVec) {
                        std::string TransName = Order_Transcript_Name[ColIndex];
                        this_ratio = static_cast<double>(AreadSJs.size())/FinallyAnnotations.Transcript_Annotations[TransName].first.size();
                        if (this_ratio < 0.4) {
                            EachClusterRowVector[ColIndex] = 0;
                            ism_gtf_name.erase(std::remove(ism_gtf_name.begin(), ism_gtf_name.end(), TransName), ism_gtf_name.end());
                        }
                    }
                }
            }            

            if (ism_gtf_name.size() != 0) {
                //只有这一行不全为0, 才能放入矩阵;
                RowIndex = RowIndex + 1;
                //调整矩阵大小以容纳第一个行向量;
                known_ISM_matrix.conservativeResize(RowIndex+1, Order_Transcript_Name.size());
                known_ISM_matrix.row(RowIndex) = EachClusterRowVector;
                //存储一个列向量, 每个类的个数, 方便后续做乘法;
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
                    std::string buffer;
                    buffer.reserve(eachISM.second.size() * 80);
                    for (const auto& EachRead : eachISM.second) {
                        buffer += EachRead;
                        buffer += '\t';
                        buffer += "ISM\t";
                        buffer += second_part;
                        buffer += '\t';
                        buffer += first_part;
                        buffer += '\t';
                        buffer += groupreadfiles.at(EachRead);
                        buffer += '\n';
                    }        
                    {
                        std::unique_lock<std::mutex> lock(traceMutex);
                        Trace << buffer;
                    }
                }
            } else { }
        }
    }

    std::set<int> MergeFalseSet(falsenodeset.begin(), falsenodeset.end());
    Eigen::MatrixXd False_novel_candidate_matrix(MergeFalseSet.size(), Order_Transcript_Name.size());

    std::set<int> FalseNode_NeighborSet;
    std::set<int> FalseNode_TrueSet;
    std::string itsname;
    std::string subStr;
    Eigen::RowVectorXd FalseNode_RowVector(Order_Transcript_Name.size());
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
                    if (Disinform.NodeDistanceSet[neinode].size() != 0) {
                        for (const auto& node:Disinform.NodeDistanceSet[neinode]) {
                            if (node < Disinform.Index2Anno.size()) {
                                itsname = Disinform.Index2Anno[node];
                                auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                                FalseNode_RowVector[itwhere] = 1;
                                ism_gtf_name.push_back(itsname);
                                countC++;
                            } else {
                                if (truenodeset.find(node) != truenodeset.end()) {
                                    // 这个邻居是novel的转录本;
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

            for (const auto& EachTranscript:ism_gtf_name) {
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
                for(const auto& EachRead : reads) {
                    buffer += EachRead;
                    buffer += '\t';
                    buffer += "simNovel\t";
                    buffer += second_part;
                    buffer += '\t';
                    buffer += first_part;
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
        False_novel_candidate_cluster_numbers[RowIndex] = Disinform.Index2Count[FalseNode];
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
    } else{
        if (False_novel_candidate_matrix.rows() != 0){
            OutputResults.Indicate_Matrix = False_novel_candidate_matrix;
            OutputResults.Cluster_Number = False_novel_candidate_cluster_numbers;
        }
    }
    OutputResults.Order_Transcript_Name_Vector = Order_Transcript_Name;
    return OutputResults;
}



IndicateFire Quantification_initialization_MultiFiles (std::map<std::size_t, std::vector<std::string>>& groupISM, 
                                                OutputInformation& FinallyAnnotations, 
                                                std::set<int>& truenodeset,
                                                std::set<int>& falsenodeset,
                                                std::unordered_map<int, std::vector<std::string>>& faslenode2count, 
                                                DistanceInform& Disinform, 
                                                std::string& group_size, 
                                                std::map<size_t, std::vector<std::string>>& HighClusters, 
                                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs,
                                                std::unordered_map<std::string, std::string>& groupreadfiles,
                                                std::ofstream& Trace) {
    IndicateFire OutputResults;
    std::vector<std::string> Order_Transcript_Name;
    Order_Transcript_Name.reserve(FinallyAnnotations.Transcript_Annotations.size());
    for (const auto& kv : FinallyAnnotations.Transcript_Annotations) {
        Order_Transcript_Name.push_back(kv.first);
    }

    Eigen::MatrixXd known_ISM_matrix;
    std::vector<std::array<int,2>> AreadSJs;
    
    int RowIndex = -1;
    int ColIndex = -1;
    int flag = 0;

    Eigen::RowVectorXd EachClusterRowVector(Order_Transcript_Name.size());
    Eigen::VectorXd known_ISM_cluster_numbers;

    std::vector<std::string> ism_gtf_name;
    std::set<std::string> first_part_set;
    std::set<std::string> second_part_set;
    std::string first_part;
    std::string second_part;
    std::vector<int> ColIndexVec;
    double max_ratio = 0;
    double this_ratio = 0;

    if (groupISM.size() != 0){
        RowIndex = -1;
        for (const auto& eachISM:groupISM){
            AreadSJs = groupreadsjs[eachISM.second[0]];
            EachClusterRowVector.setZero();
            ColIndex = -1;
            ism_gtf_name.clear();
            first_part_set.clear();
            second_part_set.clear();            
            ColIndexVec.clear();   
            double max_ratio = 0;
            double this_ratio = 0;
            for (const auto& eachTransctipt:Order_Transcript_Name){
                ColIndex = ColIndex + 1;
                if (AreadSJs.size() < FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size()){
                    if ((AreadSJs[0][0] >= FinallyAnnotations.Transcript_Annotations[eachTransctipt].first[0][0]) && (AreadSJs[AreadSJs.size()-1][1] <= FinallyAnnotations.Transcript_Annotations[eachTransctipt].first[(FinallyAnnotations.Transcript_Annotations[eachTransctipt].first.size()-1)][1])){
                        flag = whether_isoform_part_is_not(FinallyAnnotations.Transcript_Annotations[eachTransctipt].first, AreadSJs);
                        if (flag == 1){
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
                if (max_ratio >= 0.5) {
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
            if (ism_gtf_name.size() != 0){
                RowIndex = RowIndex + 1;
                known_ISM_matrix.conservativeResize(RowIndex+1, Order_Transcript_Name.size());
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
                    for (const auto& EachRead : reads) {
                        buffer += EachRead;
                        buffer += '\t';
                        buffer += "ISM\t";
                        buffer += second_part;
                        buffer += '\t';
                        buffer += first_part;
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
    }
    std::set<int> MergeFalseSet(falsenodeset.begin(), falsenodeset.end());

    Eigen::MatrixXd False_novel_candidate_matrix(MergeFalseSet.size(), Order_Transcript_Name.size());
    std::set<int> FalseNode_NeighborSet;
    std::set<int> FalseNode_TrueSet;
    std::string itsname;
    std::string subStr;
    Eigen::RowVectorXd FalseNode_RowVector(Order_Transcript_Name.size());
    RowIndex = -1;
    int countC = 0;
    
    if (MergeFalseSet.size() != 0) {

        std::set<int> MergedTrueSet(truenodeset.begin(), truenodeset.end());
        for (const auto& FalseNode:MergeFalseSet){    
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

                    } else{
                        
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
                for(const auto& neinode:FalseNode_NeighborSet){
                    if (Disinform.NodeDistanceSet[neinode].size() != 0){
                        for (const auto& node:Disinform.NodeDistanceSet[neinode]){
                            if (node < Disinform.Index2Anno.size()){
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
                for (const auto& EachRead : reads) {
                    buffer += EachRead;
                    buffer += '\t';
                    buffer += "simnovel\t";
                    buffer += second_part;
                    buffer += '\t';
                    buffer += first_part;
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
    } else{
        if (False_novel_candidate_matrix.rows() != 0){
            OutputResults.Indicate_Matrix = False_novel_candidate_matrix;
            OutputResults.Cluster_Number = False_novel_candidate_cluster_numbers;
        }
    }
    OutputResults.Order_Transcript_Name_Vector = Order_Transcript_Name;
    return OutputResults;
} 



void EM_Alg(OutputInformation& FinallyAnnotations,
            IndicateFire& Indicate_Number,
            int& rc_threshold,
            std::string& CHR) {
    
    std::vector<double> Order_Transcript_Number;

    for (const auto& eachIso:Indicate_Number.Order_Transcript_Name_Vector) {
        if (FinallyAnnotations.Transcript_Annotations.find(eachIso) != FinallyAnnotations.Transcript_Annotations.end()) {
            Order_Transcript_Number.push_back(FinallyAnnotations.Transcript_Annotations[eachIso].second);
        } else {
            Order_Transcript_Number.push_back(0);
        }
    }

    Eigen::VectorXd Order_Transcript_Vector = Eigen::Map<Eigen::VectorXd>(Order_Transcript_Number.data(), Order_Transcript_Number.size());
    Order_Transcript_Number.clear();
    
    if (Indicate_Number.Indicate_Matrix.rows() != 0 && Indicate_Number.Indicate_Matrix.cols() != 0){
        std::string its_name;
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
    
        for (int i = 0; i < AnnoN.rows(); i++) {
            its_name = Indicate_Number.Order_Transcript_Name_Vector[i];
            FinallyAnnotations.Transcript_Annotations[its_name].second = FinallyAnnotations.Transcript_Annotations[its_name].second + AnnoN(i,0);
        }
    }
}



void EM_Alg_MultiFiles(OutputInformation& FinallyAnnotations,
                       IndicateFire& Indicate_Number,
                       int& FileNumber)
{
    const std::string files = std::to_string(FileNumber);

    const auto& orderNames = Indicate_Number.Order_Transcript_Name_Vector;
    const int T = orderNames.size();
    const int R = Indicate_Number.Indicate_Matrix.rows();
    if (T == 0 || R == 0) return;

    // ================================
    // Step 0. 构建 Order_Transcript_Vector
    // ================================
    Eigen::VectorXd OrderVec(T);
    for (int i = 0; i < T; ++i) {
        auto it = FinallyAnnotations.File_TranscriptNumber[files].find(orderNames[i]);
        OrderVec[i] = (it != FinallyAnnotations.File_TranscriptNumber[files].end())
                        ? it->second
                        : 0.0;
    }

    // ================================
    // Step 1. 初始化 P
    // ================================
    Eigen::VectorXd P = Eigen::VectorXd::Constant(T, 1.0 / T);
    Eigen::VectorXd P_new(T);
    const Eigen::VectorXd& C = Indicate_Number.Cluster_Number;   // R x 1

    // ================================
    // ⭐ Step 1.5 为每行构建非零列索引（稀疏结构）
    // ================================
    std::vector<std::vector<int>> row_nonzero_j(R);
    row_nonzero_j.reserve(R);

    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < T; ++j) {
            if (Indicate_Number.Indicate_Matrix(i, j) != 0.0) {
                row_nonzero_j[i].push_back(j);
            }
        }
    }

    // ================================
    // Step 2. EM 迭代（稀疏优化）
    // ================================
    int CountCyc = 1;
    double diff_sum;

    do {
        Eigen::VectorXd denom = Indicate_Number.Indicate_Matrix * P;
        Eigen::VectorXd AnnoN = Eigen::VectorXd::Zero(T);

        // ---------- 稀疏版双层循环 ----------
        for (int i = 0; i < R; ++i) {

            double di = denom[i];
            if (di == 0.0) continue;

            const std::vector<int>& nz_cols = row_nonzero_j[i];

            for (int j : nz_cols) {
                double Aij = Indicate_Number.Indicate_Matrix(i, j);
                double zij = (P[j] * Aij) / di;
                AnnoN[j] += C[i] * zij;
            }
        }

        // ---- 更新新的 P ----
        P_new = AnnoN + OrderVec;
        double s = P_new.sum();
        if (s != 0.0) P_new /= s;
        else P_new = Eigen::VectorXd::Constant(T, 1.0 / T);

        diff_sum = (P_new - P).cwiseAbs().sum();
        P = P_new;

        CountCyc++;
    } 
    while (diff_sum > 5e-2 && CountCyc < 9);

    // ================================
    // Step 3. 按最终 P 重新计算 AnnoN（稀疏优化）
    // ================================
    Eigen::VectorXd final_denom = Indicate_Number.Indicate_Matrix * P;
    Eigen::VectorXd final_AnnoN = Eigen::VectorXd::Zero(T);

    for (int i = 0; i < R; ++i) {

        double di = final_denom[i];
        if (di == 0.0) continue;

        const std::vector<int>& nz_cols = row_nonzero_j[i];

        for (int j : nz_cols) {
            double Aij = Indicate_Number.Indicate_Matrix(i, j);
            double zij = (P[j] * Aij) / di;
            final_AnnoN[j] += C[i] * zij;
        }
    }

    // ================================
    // Step 4. 写回 File_TranscriptNumber
    // ================================
    for (int j = 0; j < T; ++j) {
        const std::string& name = orderNames[j];
        FinallyAnnotations.File_TranscriptNumber[files][name] += final_AnnoN[j];
    }
}


void get_filter_Low(std::unordered_map<std::size_t, std::vector<std::string>>& LowClusters, 
                    OutputInformation& FinallyAnnotations,
                    std::unordered_map<std::string, std::vector<std::array<int,2>>>& tx2sj,
                    std::unordered_map<std::string, double>& AllFileGeneCounts, std::unordered_map<std::string, std::vector<std::string>>& gene2tx,
                    std::ofstream& isoformCountPath, int& rc_threshold, SE_belong2_genetranscript& groupSEREADS,
                    std::ofstream& Updated_Files, std::unordered_map<std::string, std::vector<std::array<int,2>>>& rawgtf_isoform,
                    std::string& chrname, std::unordered_map<std::string, std::string>& rawgtf_strand, 
                    std::ofstream& Trace, std::unordered_map<std::string, std::string>& groupreadfiles) {

    auto& annoMap = FinallyAnnotations.Transcript_Annotations;
    auto& transcript2gene = FinallyAnnotations.transcript2gene;

    std::vector<std::array<int,2>> ThisSj;
    std::vector<std::array<int,2>> AnnoSj;
    std::vector<std::string> TempTx;
    std::string geneName, TxTxName;

    for (const auto& EachLow:LowClusters) {
        auto it_tx = tx2sj.find(EachLow.second[0]);
        if (it_tx == tx2sj.end()) continue;
        const auto& ThisSjRef = it_tx->second;

        ThisSj = tx2sj[EachLow.second[0]];
        TempTx.clear();

        for (auto& EachAnno:annoMap) {
            // novel的才判定;
            const std::string& annoName = EachAnno.first;
            if (annoName.find('|') != std::string::npos) continue;
            const auto& AnnoSjRef = EachAnno.second.first;
            int Sect = IntervalIntersection(ThisSjRef, AnnoSjRef);
            if (Sect <= 0) continue;
            int Distance = IntervalMerge(ThisSjRef, AnnoSjRef) - Sect - findNonOverlappingIntervals(ThisSjRef, AnnoSjRef);
            if (Distance < 60 && Distance > 0) {
                TempTx.push_back(annoName);
            }
        }
        // 每个cluster循环一圈;
        if (TempTx.empty()) continue;
        const size_t clusterSize = EachLow.second.size();
        double TempCount = 0.0;

        std::string buffer;
        buffer.reserve(EachLow.second.size() * 80);
        if (TempTx.size() == 1) {
            auto& anno = annoMap[TempTx[0]];
            anno.second += clusterSize;
            geneName = transcript2gene[TempTx[0]];
            for (const auto& eachR:EachLow.second) {
                buffer += eachR;
                buffer += '\t';
                buffer += "simnovel\t";
                buffer += TempTx[0];
                buffer += '\t';
                buffer += geneName;
                buffer += '\t';
                buffer += groupreadfiles.at(eachR);
                buffer += '\n';
            } 
            {
                std::unique_lock<std::mutex> lock(traceMutex);
                Trace << buffer;
            }
        } else {
            for (const auto& TxTx:TempTx) {
                auto it = annoMap.find(TxTx);
                if (it != annoMap.end() && it->second.second > TempCount) {
                    TxTxName = TxTx;
                    TempCount = it->second.second;
                }
            }
            if (TempCount > 0) {
                auto& anno = annoMap[TxTxName];
                anno.second += clusterSize;
                geneName = transcript2gene[TxTxName];
                for (const auto& eachR:EachLow.second) {
                    buffer += eachR;
                    buffer += '\t';
                    buffer += "simnovel\t";
                    buffer += TempTx[0];
                    buffer += '\t';
                    buffer += geneName;
                    buffer += '\t';
                    buffer += groupreadfiles.at(eachR);
                    buffer += '\n';
                } 
                {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << buffer;
                }
            }
        }
    }  

    for (const auto& eachAnno:annoMap) {
        const std::string& name = eachAnno.first;
        geneName = transcript2gene[name];
        double countVal = eachAnno.second.second;
        if (geneName != "NA") { AllFileGeneCounts[geneName] += countVal; }
    }



    // 还有那些单exon reads归类于多exon isoform;
    const auto& seReads = groupSEREADS.Transcript_with_SE_reads;
    const bool updateOpen = Updated_Files.is_open();

    for (const auto& eachAnno:seReads) {
        auto it = annoMap.find(eachAnno.first);
        if (it != annoMap.end()) {
            // 已经有这种转录本;
            it->second.second += eachAnno.second.size();
        } else {
            // 如果没有这个转录本, 还要输入到gtf中;
            size_t pos = eachAnno.first.find('|');
            std::string first_part, second_part;
            if (pos != std::string::npos) {
                first_part = eachAnno.first.substr(0, pos);
                second_part = eachAnno.first.substr(pos + 1);
            } 
            if (updateOpen) {
                std::ostringstream localBuffer;
                const auto& thisVec = rawgtf_isoform[eachAnno.first];
                const auto& strand = rawgtf_strand[eachAnno.first];
                localBuffer << chrname << '\t' << "annotated_isoform" << '\t' << "transcript" << '\t'
                            << thisVec.front()[0] << '\t' << thisVec.back()[1]
                            << '\t' << "." << '\t' << strand
                            << '\t' << '.' << '\t' << "gene_id \"" << first_part 
                            << "\"; transcript_id \"" << second_part << "\";" << '\n';
                for (size_t i = 0; i < thisVec.size(); i++) {
                    localBuffer << chrname << '\t' << "annotated_isoform" << '\t' << "exon" << '\t'
                                << thisVec[i][0] << '\t' << thisVec[i][1] << '\t' << "." << '\t'
                                << strand << '\t' << '.' << '\t'
                                << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part 
                                << "\"; exon_number \"" << i << "\";" << '\n';                    
                }
                {
                    std::unique_lock<std::mutex> lock(updatedGtfMutex);
                    Updated_Files << localBuffer.str();
                }
            }        
            annoMap[eachAnno.first] = {};
            annoMap[eachAnno.first].first = {};
            annoMap[eachAnno.first].second = eachAnno.second.size();
            transcript2gene[eachAnno.first] = first_part;
        }
    }

    // 输出定量文件;
    if (isoformCountPath.is_open()) {
        std::ostringstream buffer;
        for (const auto& eachAnno:annoMap) {
            const std::string& name = eachAnno.first;
            size_t pos = name.find('|');
            std::string second_part = (pos != std::string::npos) ? name.substr(pos + 1) : name;
            geneName = transcript2gene[name];
            double countVal = eachAnno.second.second;
            if (countVal >= rc_threshold) {
                buffer << second_part << '\t' << geneName << '\t' << countVal << '\n';
            }
        }
        {
            std::unique_lock<std::mutex> lock(isoformCountMutex);
            isoformCountPath << buffer.str();
        }
    }
}


void get_filter_Low_k_Files(std::map<std::string, std::map<std::size_t, std::vector<std::string>>>& LowClusters, 
                    OutputInformation& FinallyAnnotations,
                    std::unordered_map<std::string, std::vector<std::array<int,2>>>& tx2sj,
                    std::unordered_map<int, std::unordered_map<std::string, double>>& AllFileGeneCounts,
                    std::unordered_map<std::string, std::vector<std::string>>& gene2tx, int& AllFileNumber) {

    for (int i = 0; i < AllFileNumber; i++) {
        AllFileGeneCounts[i].reserve(512); 
    }
    // std::vector<std::array<int,2>> ThisSj;
    // std::vector<std::array<int,2>> AnnoSj;
    // std::vector<std::string> TempTx;
    // TempTx.reserve(64);

    // 每个文件;
    // for (int q = 0; q < AllFileNumber; q++) {
    //     const std::string q_str = std::to_string(q);
    //     auto lowClusterIt = LowClusters.find(q_str);
    //     if (lowClusterIt == LowClusters.end()) continue;
    //     auto& thisLowCluster = lowClusterIt->second;
    //     auto& thisFileAnno = FinallyAnnotations.File_TranscriptNumber[q_str];
    //     auto& thisFileGeneCounts = AllFileGeneCounts[q];

    //     for (const auto& EachLow:thisLowCluster) {
    //         const auto it_sj = tx2sj.find(EachLow.second[0]);
    //         if (it_sj == tx2sj.end()) continue;
    //         ThisSj = it_sj->second;
    //         TempTx.clear();
            
    //         for (auto& EachAnno:thisFileAnno) {
    //             // novel的才判定;
    //             const size_t pos = EachAnno.first.find('|');
    //             if (pos != std::string::npos) continue;
    //             const auto it_ta = FinallyAnnotations.Transcript_Annotations.find(EachAnno.first);
    //             if (it_ta == FinallyAnnotations.Transcript_Annotations.end()) continue;

    //             const auto& AnnoSjRef = it_ta->second.first;
    //             const int Sect = IntervalIntersection(ThisSj, AnnoSjRef);
    //             if (Sect <= 0) continue;

    //             const int Distance = IntervalMerge(ThisSj, AnnoSjRef) - Sect - findNonOverlappingIntervals(ThisSj, AnnoSjRef);
    //             if (Distance < 60 && Distance > 0) {
    //                 TempTx.push_back(EachAnno.first);
    //             }
    //         }

    //         // 每个cluster循环一圈;
    //         if (TempTx.empty()) continue;
    //         const size_t cluster_size = EachLow.second.size();

    //         if (TempTx.size() == 1) {
    //             const std::string& txname = TempTx[0];
    //             thisFileAnno[txname] += cluster_size;
    //             const std::string& geneName = FinallyAnnotations.transcript2gene[txname];
    //             if (geneName != "NA") thisFileGeneCounts[geneName] += cluster_size;             
    //         } else {
    //             double TempCount = 0.0;
    //             std::string TxTxName;
    //             for (const auto& TxTx : TempTx) {
    //                 auto it = thisFileAnno.find(TxTx);
    //                 if (it != thisFileAnno.end() && it->second > TempCount) {
    //                     TxTxName = TxTx;
    //                     TempCount = it->second;
    //                 }
    //             }
    //             if (TempCount > 0) {
    //                 thisFileAnno[TxTxName] += cluster_size;
    //                 const std::string& geneName = FinallyAnnotations.transcript2gene[TxTxName];
    //                 if (geneName != "NA") thisFileGeneCounts[geneName] += cluster_size;
    //             }                
    //         }
    //     }
    // }
}



std::unordered_map<int, std::unordered_map<std::string, double>> get_allfile_genecounts(
                                            OutputInformation& FinallyAnnotations, 
                                            int& AllFileNumber, 
                                            std::ofstream& isoformCountPath,
                                            int& rc_threshold){
    std::unordered_map<int, std::unordered_map<std::string, double>> geneCounts;
    geneCounts.reserve(AllFileNumber);

    std::string geneName;
    std::string second_part;

    for (int i = 0; i < AllFileNumber; i++) {
        geneCounts[i].reserve(1000);
    }

    std::string fileBuffer;
    fileBuffer.reserve(5000000);

    auto& FTN = FinallyAnnotations.File_TranscriptNumber;
    auto& T2G = FinallyAnnotations.transcript2gene;

    for (const auto& eachAnno:FTN["0"]) {
        const std::string& transcript = eachAnno.first;
        size_t pos = transcript.find('|');
        std::string second_part = (pos != std::string::npos ? transcript.substr(pos + 1) : transcript);
        const std::string& geneName = T2G.at(transcript);

        fileBuffer.append(second_part);
        fileBuffer.push_back('\t');
        fileBuffer.append(geneName);

        for (int i = 0; i < AllFileNumber; i++) {
            auto& inner = FTN[std::to_string(i)];
            double v = inner.count(transcript) ? inner.at(transcript) : 0;
            if (v < rc_threshold) v = 0;
            fileBuffer.push_back('\t');
            fileBuffer.append(std::to_string(v));
        }
        fileBuffer.push_back('\n');

        if (geneName != "NA") {
            for (int i = 0; i < AllFileNumber; i++) {
                auto& inner = FTN[std::to_string(i)];
                double v = inner.count(transcript) ? inner.at(transcript) : 0;
                geneCounts[i][geneName] += v;
            }
        }           
    }

    {
        std::unique_lock<std::mutex> lock(isoformCountMutex);
        isoformCountPath << fileBuffer;
    }
    return geneCounts;
}



std::unordered_map<int, std::unordered_map<std::string, double>> DetectQuant(GroupAnnotation& groupanno, 
                 SpliceChainClass& splicechainclass, 
                 GroupInformation& groupinform, 
                 const int& groupdistance,
                 unGTF& gtfexon, GTFsj& gtfsjs,
                 std::ofstream& gtfFilePath, std::ofstream& isoformFilePath,
                 std::ofstream& traceFilePath, int& FileNo,
                 std::unordered_map<std::string, std::array<int,2>>& singleexonanno, int& RC_threshold,
                 SE_belong2_genetranscript& groupSeReads, std::string& CHRCHR) {

    std::string chrchr = groupinform.chrName;
    std::string groupnumber = groupinform.GroupIndex;
    DistanceInform DMatrix_GraphNode = get_distance_matrix(groupanno.group_me_SJs, 
                                                           splicechainclass.HighConClusters, 
                                                           groupinform.GroupReadSjs, 
                                                           groupdistance, splicechainclass.FSM);

    std::unordered_map<int, std::unordered_map<std::string, double>> AllFile_GeneCounts;

    if (DMatrix_GraphNode.NodeDistanceSet.size() != 0) { 
        std::vector<std::vector<int>> CliquesVector = Find_Maximal_Cliques(DMatrix_GraphNode.NodeDistanceSet);

        DetectionResults NodeResults = Transcript_Detection(CliquesVector, 
                                                            DMatrix_GraphNode.Index2Anno, 
                                                            DMatrix_GraphNode.Index2Unknown, 
                                                            DMatrix_GraphNode.Index2Count);
        CliquesVector.clear();

        if (FileNo > 1) {
            Solvent SpliceChainSolvent = get_Solvent_FsmIsmHigh(splicechainclass, FileNo, groupinform);
   
            OutputInformation Finally_Annotations = Write_Detection_Transcript2gtf_MultiFiles(gtfFilePath, traceFilePath, 
                                                    NodeResults, DMatrix_GraphNode, 
                                                    gtfsjs.mSJs[chrchr], gtfexon.GTF_transcript[chrchr], 
                                                    splicechainclass.ClusterCoverage, groupanno.group_genes, 
                                                    gtfexon.GTF_gene[chrchr], 
                                                    gtfexon.GTF_transcript_strand[chrchr], 
                                                    splicechainclass.HighStrand, chrchr, groupnumber, 
                                                    splicechainclass.HighConClusters,
                                                    FileNo, groupinform.GroupReadFiles,
                                                    SpliceChainSolvent,
                                                    gtfexon.GTF_gene2transcript[chrchr]); 

            for (int k = 0; k < FileNo; k++) {
                FileFalseNode This_File_False_Node = get_File_False_node(NodeResults.FalseNodeSet,
                                                                        DMatrix_GraphNode.Index2hashname,
                                                                        SpliceChainSolvent.File_HighConClusters,
                                                                        k);

                IndicateFire InitFirefly = Quantification_initialization_MultiFiles(SpliceChainSolvent.File_ISM[std::to_string(k)], 
                                                            Finally_Annotations, NodeResults.TrueNodeSet,
                                                            This_File_False_Node.ThisFlaseNode, 
                                                            This_File_False_Node.FasleNode2Count,
                                                            DMatrix_GraphNode, groupnumber, 
                                                            SpliceChainSolvent.File_HighConClusters[std::to_string(k)], 
                                                            groupinform.GroupReadSjs, 
                                                            groupinform.GroupReadFiles,
                                                            traceFilePath); 
                EM_Alg_MultiFiles(Finally_Annotations, InitFirefly, k);              
            }

            get_filter_Low_k_Files(SpliceChainSolvent.File_LowConClusters, Finally_Annotations, groupinform.GroupReadSjs,
                                    AllFile_GeneCounts, gtfexon.GTF_gene2transcript[chrchr], FileNo); 
             
            AllFile_GeneCounts = get_allfile_genecounts(Finally_Annotations, FileNo, isoformFilePath, RC_threshold);

            return AllFile_GeneCounts; 

        } else {
            AllFile_GeneCounts[0] = {};

            OutputInformation Finally_Annotations = Write_Detection_Transcript2gtf(gtfFilePath, traceFilePath, 
                                                    NodeResults, DMatrix_GraphNode, 
                                                    gtfsjs.mSJs[chrchr], gtfexon.GTF_transcript[chrchr], 
                                                    splicechainclass.ClusterCoverage, groupanno.group_genes, 
                                                    gtfexon.GTF_gene[chrchr], 
                                                    gtfexon.GTF_transcript_strand[chrchr], 
                                                    splicechainclass.HighStrand, chrchr, groupnumber, 
                                                    splicechainclass.HighConClusters,
                                                    groupinform.GroupReadFiles,
                                                    gtfexon.GTF_gene2transcript[chrchr]);

            IndicateFire InitFirefly = Quantification_initialization(splicechainclass.ISM, 
                                                        Finally_Annotations, NodeResults.TrueNodeSet,
                                                        NodeResults.FalseNodeSet, 
                                                        DMatrix_GraphNode, groupnumber, 
                                                        splicechainclass.HighConClusters, 
                                                        groupinform.GroupReadSjs,
                                                        groupinform.GroupReadFiles, 
                                                        traceFilePath,
                                                        groupinform.chrName);
            EM_Alg(Finally_Annotations, InitFirefly, RC_threshold, CHRCHR);

            // Low的跟Novel的能差多少;
            get_filter_Low(splicechainclass.LowConClusters, Finally_Annotations, groupinform.GroupReadSjs, 
                            AllFile_GeneCounts[0], gtfexon.GTF_gene2transcript[chrchr],
                            isoformFilePath, RC_threshold, groupSeReads, gtfFilePath, gtfexon.GTF_transcript[chrchr],
                            chrchr, gtfexon.GTF_transcript_strand[chrchr], traceFilePath, groupinform.GroupReadFiles);
        
            splicechainclass.LowConClusters.clear();
            return AllFile_GeneCounts;
        }
    } else {
        return AllFile_GeneCounts;
    }
}


void write_gene_counts(std::unordered_map<int, std::unordered_map<std::string, double>>& FileSingleExonGeneNumber,
                    std::unordered_map<int, std::unordered_map<std::string, double>>& FileMultiExonGeneNumber,
                    std::ofstream& geneFilePath, int& FileNo) {
    if (!geneFilePath.is_open()) return;
    std::string buffer;

    // ======== Case 1: 只有一个文件 ========
    if (FileNo == 1) {
        auto& single0 = FileSingleExonGeneNumber[0];
        const auto& multi0 = FileMultiExonGeneNumber[0];

        // 融合多 exon → 单 exon
        if (!multi0.empty()) {
            if (!single0.empty()) {
                for (auto& kv : single0) {
                    auto it = multi0.find(kv.first);
                    if (it != multi0.end())
                        kv.second += it->second;
                }
            } else {
                // single 为空 → 用 multi
                for (const auto& kv : multi0)
                    single0[kv.first] = kv.second;
            }
        }

        // 输出 single0
        for (const auto& kv : FileSingleExonGeneNumber[0]) {
            if (kv.second > 0) {
                buffer.append(kv.first);
                buffer.push_back('\t');
                buffer.append(std::to_string(kv.second));
                buffer.push_back('\n');
            }
        }
        std::unique_lock<std::mutex> lock(geneCountMutex);
        geneFilePath << buffer;
        return;
    }

    // ======== Case 2: 多文件 ========
    const auto& multi0 = FileMultiExonGeneNumber[0];
    auto& single0 = FileSingleExonGeneNumber[0];
    if (!multi0.empty()) {
        if (!single0.empty()) {
            // 把 multi 的值加到 single
            for (const auto& kv : multi0) {
                const std::string& gene = kv.first;
                double multiCount = kv.second;
                if (single0.find(gene) != single0.end()) {
                    // 所有文件都加
                    for (int i = 0; i < FileNo; ++i) {
                        FileSingleExonGeneNumber[i][gene] += FileMultiExonGeneNumber[i][gene];
                    }
                }
            }
        } else {
            // single 为空 → 内容来自 multi
            for (const auto& kv : multi0) {
                const std::string& gene = kv.first;
                for (int i = 0; i < FileNo; ++i) {
                    FileSingleExonGeneNumber[i][gene] = FileMultiExonGeneNumber[i][gene];
                }
            }
        }
    }    

    // （2）输出（单行多列）
    const auto& final0 = FileSingleExonGeneNumber[0];
    for (const auto& kv : final0) {
        const std::string& gene = kv.first;
        buffer.append(gene);
        for (int i = 0; i < FileNo; ++i) {
            buffer.push_back('\t');
            buffer.append(std::to_string(FileSingleExonGeneNumber[i][gene]));
        }
        buffer.push_back('\n');
    }    
    {
        std::unique_lock<std::mutex> lock(geneCountMutex);
        geneFilePath << buffer;
    }
}



void processGroup(std::streampos& start, std::streampos& end, 
                  const std::string& sam_file_path, 
                  const int& Sj_supportReadNumber, unGTF& gtf_full, 
                  GTFsj& gtf_splice, const int& GraphDis,
                  std::ofstream& updatedgtffile, std::ofstream& isoformcountfile,
                  std::ofstream& genecountfile, std::ofstream& tracefile,
                  int& fileno, int& singleEdge, int& ReadCount_threshold) {
    
    GroupInformation group_information = knowGroupInformation(start, end, sam_file_path, Sj_supportReadNumber);
    
    std::string chrchr = group_information.chrName;
    std::array<int,2> groupcoverage = group_information.GroupCoverage;

    SE_belong2_genetranscript group_se_reads;
    GroupAnnotation group_annotation;
    std::unordered_map<int, std::unordered_map<std::string, double>> file_multiexon_gene_number; 
    SpliceChainClass spliceclass;

    if ( (group_information.GroupSingleExon.size() > 0) or (group_information.GroupReadSjs.size() > 0) ) {
        group_annotation = get_group_single_exon_gene_annotation(gtf_full.GTF_gene[chrchr], groupcoverage,
                                gtf_full.GTF_gene2transcript[chrchr], gtf_full.GTF_transcript[chrchr], gtf_splice.mSJs[chrchr], chrchr);
    }

    if (group_information.GroupReadSjs.size() > 0) {
        spliceclass = generate_splice_chain_class(group_information.GroupReadSjs, group_information.GroupReadCoverage, 
                                                group_annotation.group_me_SJs, group_information.GroupSigns, 
                                                gtf_splice.mSJsBE[chrchr], group_information.GroupReadFiles, 
                                                tracefile, Sj_supportReadNumber, chrchr); 
        group_information.GroupSigns.clear();
        group_information.GroupReadCoverage.clear();
    }

    if (group_information.GroupSingleExon.size() > 0) {
        group_se_reads = get_group_singleexon_reads_2gene(group_annotation, group_information.GroupSingleExon, 
                                        group_information.GroupReadFiles, gtf_full.GTF_gene2transcript[chrchr], 
                                        gtf_full.GTF_transcript[chrchr], fileno, singleEdge, chrchr, spliceclass.FSM,
                                        tracefile);          
        write_single_exon_gtf_trace(fileno, group_information.GroupReadFiles, group_se_reads.Transcript_with_SE_reads, 
                                    group_annotation.group_se_transcripts, updatedgtffile, isoformcountfile, 
                                    gtf_full.GTF_transcript_strand[chrchr], chrchr); // 单exon的reads是单exon转录本的都输出了;  
    }
    group_information.GroupSingleExon.clear();        

    file_multiexon_gene_number = DetectQuant(group_annotation, spliceclass, group_information, 
                                    GraphDis, gtf_full, gtf_splice, updatedgtffile, isoformcountfile, 
                                    tracefile, fileno, group_annotation.group_se_transcripts, ReadCount_threshold,
                                    group_se_reads, chrchr);

    spliceclass.FSM.clear(); spliceclass.HighConClusters.clear();spliceclass.ISM.clear();spliceclass.LowConClusters.clear();spliceclass.HighStrand.clear();

    write_gene_counts(group_se_reads.file_SE_reads_gene_number, file_multiexon_gene_number, genecountfile, fileno);

    if (group_information.GroupReadSjs.size() + group_information.GroupSingleExon.size() > 10000) {
        std::unique_lock<std::mutex> lock(bigMutex);
        std::cerr << (group_information.GroupReadSjs.size() + group_information.GroupSingleExon.size()) / ((double)1000000) << " M reads processed.." << std::endl;
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
    int Read_count = 1;
    int single_exon_edge = 300;
    int mapq = 0;

    std::cerr << R"(
     ____               ____   ___  _      ___ 
    | __ )  _ __ ___   / ___| / _ \| |    |_ _|
    |  _ \ | '__/ _ \ | |    | | | | |     | | 
    | |_) || | | (_) || |___ | |_| | |___  | | 
    |____/ |_|  \___/  \____| \___/|_____||___|
    )" << std::endl;

    std::cerr << "         BroCOLI  Version: 1.0.0" << std::endl;

    while ((c = getopt_long(argc, argv, "s:f:g:o:j:n:m:e:d:t:r:h", long_options, &option_index)) != -1) {
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

    std::vector<std::string> outputFileVec = check_catalog_exist(output_file_name, sam_file_vec);
    std::ofstream gtf_file(outputFileVec[0], std::ios::app);
    std::ofstream isoform_file(outputFileVec[1], std::ios::app);
    std::ofstream gene_file(outputFileVec[2], std::ios::app);
    std::ofstream trace_file(outputFileVec[3], std::ios::app);

    std::cerr << "*****" << std::endl;
    std::cerr << "Input file: " << samfile_name << std::endl;
    std::cerr << "FASTA file: " << fastafile_name << std::endl;
    std::cerr << "GTF file: " << gtffile_name << std::endl;
    std::cerr << "Output file: " << output_file_name << std::endl;
    std::cerr << "Single exon boundary: " << single_exon_edge << std::endl;
    std::cerr << "SJ Distance: " << SJDistance << std::endl;
    std::cerr << "SJ support read number: " << SJ_support_read_number << std::endl;
    std::cerr << "MAPQ: " << mapq << std::endl;
    std::cerr << "Graph distance: " << Graph_distance << std::endl;
    std::cerr << "Thread: " << Thread << std::endl;
    std::cerr << "Output min read count: " << Read_count << std::endl;
    std::cerr << "*****" << std::endl;
    
    std::cerr << "*** " << "Read and process the files ......\n";
    std::unordered_map<std::string, std::string> Fasta = Read_fasta_file(fastafile_name);

    auto start = std::chrono::high_resolution_clock::now();
    FileSplit BroCOLIfile = thread_all_read_sam_files(samfile_name, sam_file_vec, Thread, output_file_name, SJDistance, Fasta, mapq);
    Fasta.clear();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cerr << "Read file information cost time = " << diff.count() << " s\n";    

    unGTF GTF_full = get_gtf_annotation(gtffile_name);
    GTFsj GTF_Splice = get_SJs_SE(GTF_full.GTF_transcript);
    std::vector<std::size_t> Group_idx = sort_indexes_e(BroCOLIfile.group_reads_number);
    std::cerr << "*** " << "File processing completed! " << "***\n";

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
                GTF_full, 
                GTF_Splice,
                Graph_distance,
                gtf_file,
                isoform_file,
                gene_file,
                trace_file,
                BroCOLIfile.FileNo,
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
    isoform_file.close();
    gene_file.close();
    trace_file.close();
    std::cerr << "*** BroCOLI quantification has been successfully completed! ***\n";

    return 0;
}






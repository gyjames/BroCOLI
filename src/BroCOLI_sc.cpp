#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <unordered_map>
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

std::mutex updatedGtfMutex;
std::mutex traceMutex;
std::mutex bigMutex;


struct Split_Result{
    std::string read_name;
    std::vector<std::string> tokens;
    int read_length;
    std::string read_barcode;
    std::string read_umi;
};


Split_Result get_string_split_sc(const std::string& s, char delimiter, std::string& umi_tag_fourth, std::string& barcode_tag_fourth) {

    Split_Result read_result = {};

    std::istringstream iss(s);

    std::string token;
    int umi_size = umi_tag_fourth.size();
    int barcode_size = barcode_tag_fourth.size();

    int number = 0;
    while (std::getline(iss, token, delimiter)) {
        number++;

        if (number == 1){
            read_result.read_name = token;
        } else if (1 < number && number <= 6) {

            read_result.tokens.push_back(token);
        } else if (number == 10) {
            read_result.read_length = token.size();
        } else if (token.size() > 10 && token.compare(0, umi_size, umi_tag_fourth) == 0) {
            read_result.read_umi = token;
        } else if (token.size() > 10 && token.compare(0, barcode_size, barcode_tag_fourth) == 0) {
            read_result.read_barcode = token;
        }

        if (read_result.read_umi.size() > 0 && read_result.read_barcode.size() > 0) {
            read_result.read_name = read_result.read_barcode + "-" + read_result.read_umi;
        }
    }
    return read_result;
}



int get_read_match_gene_length(std::vector<std::array<int,2>> CIGAR_vector){
    //初始化变量, 累加M前面的数字之和
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

Read_intervals_and_Mlength get_read_intervals(std::string CIGARvalue, std::string initpos){

    Read_intervals_and_Mlength CIGAR_Interval_Length;
    int MatchLength = 0;
    
    int position_begin = std::stoi(initpos);
    int position_end = 0;

    char currentChar = '0';
    int currentNumber = 0;
    bool inNumberMode = false;
    
    for (char ch : CIGARvalue) {
        //判断string的每个字符是否是数字;
        if (std::isdigit(ch)) {
            inNumberMode = true;
            currentNumber = currentNumber * 10 + (ch - '0');
        } else {
            currentChar = ch;
            //CIGAR不同的标记符号;
            if (currentChar == 'M') {
                MatchLength = MatchLength + currentNumber;
                //按照这种加法, 这个M匹配的小区间是前闭后开的;
                position_end = position_begin + currentNumber;
                //初始化每个小区间;
                std::array<int,2> each_small_QJ;
                each_small_QJ[0] = position_begin;
                each_small_QJ[1] = position_end;
                CIGAR_Interval_Length.ReadIntervals.push_back(each_small_QJ);

            } else if (currentChar == 'D') {
                position_end = position_begin + currentNumber;

            } else if (currentChar == 'I') {
                MatchLength = MatchLength + currentNumber;
                position_end = position_begin;

            } else if (currentChar == 'N') {
                position_end = position_begin + currentNumber;

            } else if (currentChar == 'S') {
                MatchLength = MatchLength + currentNumber;
                position_end = position_begin;

            } else if (currentChar == 'H') {
                MatchLength = MatchLength;
                position_end = position_begin;

            } else if (currentChar == 'P') {
                MatchLength = MatchLength + currentNumber;
                position_end = position_begin + currentNumber;

            } else {

                std::cout << currentChar << ' ' << "has not been considered yet!" << std::endl;
            }

            position_begin = position_end;

            inNumberMode = false;
            currentNumber = 0;
   
        }    
    }
    CIGAR_Interval_Length.ReadMatchLength = MatchLength;
    return CIGAR_Interval_Length;
}



struct Reads_Clusters {

    std::streampos lastPos;

    std::streampos newPos;

    std::map<std::string, std::vector<std::array<int,2>>> Mymap;

    std::unordered_map<std::string, std::string> ReadsBarcode;

    std::unordered_map<std::string, std::string> ReadsUMI;

    std::map<std::string, int> Mylen;

    std::string SetRef_name;

    std::array<int,2> ClusterCoverage;

    int mapqless1 = 0;
    int sameumi = 0;
    int mappingdiff = 0;
    int surveyNum = 0;    
}; 


Reads_Clusters get_each_cluster_reads_sc(std::ifstream& samfile, std::streampos CurrentPos, std::streampos EndPos,
                                        std::string& umi_tag_third, std::string& barcode_tag_third) {

    Reads_Clusters NewCluster = {};

    std::map<std::string, std::vector<std::array<int,2>>> read_informs;
    std::map<std::string, int> read_len;

    std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, int>>> Barcode2UMI2reads;
    //记录每条reads的barcode;
    std::unordered_map<std::string, std::string> read_of_Barcode;
    std::unordered_map<std::string, std::string> read_of_UMI;

    std::string line;
    std::string now_gene;
    std::string last_chr;
    //跳转到输入流当前的位置, 继续读取;
    samfile.seekg(CurrentPos, std::ios::beg);
    //保存上一行的位置;
    std::streampos earlyPos = CurrentPos;
    //定义循环cluster的判定条件;
    int early_begin_pos = 0;
    int early_end_pos = 0;

	while (getline(samfile, line))
	{

        earlyPos = CurrentPos; 
        CurrentPos = samfile.tellg();

        if (line[0] != '@'){
            NewCluster.surveyNum = NewCluster.surveyNum + 1;
            Split_Result This_Line = get_string_split_sc(line, '\t', umi_tag_third, barcode_tag_third); 
            if (This_Line.read_barcode.size() == 0 || This_Line.read_umi.size() == 0) {
                continue;
            }
            //判断MAPQ是否大于1;
            int read_mapq = std::stoi(This_Line.tokens[3]);
            if (read_mapq > 1) {
                auto thisbarcode = Barcode2UMI2reads.find(This_Line.read_barcode);
                if (thisbarcode == Barcode2UMI2reads.end()) {
                    Barcode2UMI2reads[This_Line.read_barcode] = {};
                }
                auto thisumi = Barcode2UMI2reads[This_Line.read_barcode].find(This_Line.read_umi);
                if ( (thisumi == Barcode2UMI2reads[This_Line.read_barcode].end() ) || ( (thisumi != Barcode2UMI2reads[This_Line.read_barcode].end()) && (This_Line.read_length > Barcode2UMI2reads[This_Line.read_barcode][This_Line.read_umi].second) ) ) {

                    Barcode2UMI2reads[This_Line.read_barcode][This_Line.read_umi].first = This_Line.read_name;
                    Barcode2UMI2reads[This_Line.read_barcode][This_Line.read_umi].second = This_Line.read_length;                      
                    //获取每条reads match的区间;
                    Read_intervals_and_Mlength CIGAR_interval = get_read_intervals(This_Line.tokens[4], This_Line.tokens[2]);
                    if (This_Line.read_length == CIGAR_interval.ReadMatchLength) {

                        now_gene = This_Line.tokens[1];
                        //read匹配的开始位置;
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
                                NewCluster.ClusterCoverage[0] = early_begin_pos;
                                NewCluster.ClusterCoverage[1] = early_end_pos;
                                NewCluster.SetRef_name = last_chr;
                                NewCluster.ReadsBarcode = read_of_Barcode;
                                NewCluster.ReadsUMI = read_of_UMI;

                                NewCluster.surveyNum = NewCluster.surveyNum - 1;
                                break;
                            } else{
                                if ((read_in_gene_begin_pos - early_end_pos > 0) || (early_begin_pos - read_in_gene_end_pos > 0)){

                                    NewCluster.lastPos = earlyPos;
                                    NewCluster.newPos = CurrentPos;
                                    NewCluster.Mymap = read_informs;
                                    NewCluster.Mylen = read_len;
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
        //存储的文件读取位置;
        NewCluster.lastPos = earlyPos;
        NewCluster.newPos = CurrentPos;
        //存储的reads M小区间和所在的染色体信息;
        NewCluster.Mymap = read_informs;
        NewCluster.Mylen = read_len;
        NewCluster.ClusterCoverage[0] = early_begin_pos;
        NewCluster.ClusterCoverage[1] = early_end_pos;
        NewCluster.SetRef_name = last_chr;
        NewCluster.ReadsBarcode = read_of_Barcode;
        NewCluster.ReadsUMI = read_of_UMI;
    }
    //返回reads cluster的结构体;
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
        //提取出来;
        std::string Rname = pair.first; //read名称;
        std::vector<std::array<int,2>> RMs = pair.second;//read M小区间;

        std::vector<std::array<int,2>> ReadSJs;
        std::vector<std::array<int,2>> ReadSJs_left;
        std::vector<std::array<int,2>> ReadSJs_right;
        //初始化单个exon开始和结尾的array;
        std::array<int,2> Read_SingleExon;
        std::array<int,2> Read_begin_end;

        for (auto it = RMs.begin(); it!= RMs.end()-1; it++) {
            //这两个就是为了确定read的开头和结尾;
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
        //判断结束;
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
    std::map<std::string, std::map<std::string, std::array<int,2>>> GTF_gene; //基因的范围;
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> GTF_transcript_strand;
    std::map<std::string, std::map<std::string, std::string>> GTF_gene_strand;
    std::map<std::string, std::map<std::string, std::vector<std::string>>> GTF_gene2transcript;
};

GTF get_gtf_annotation(std::string& GTFFile_name){
    GTF GTFAll_Info;
    if (!GTFFile_name.empty()){
        //正确读取;
        std::cout << "***** Now open the gtf file: " << GTFFile_name << "! *****" << std::endl;        
        // 每个染色体上的转录本及其所有exon;
        std::map<std::string, std::vector<std::array<int,2>>> ChrEach;
        std::unordered_map<std::string, std::string> ChrTranscriptStrand;
        std::map<std::string, std::array<int,2>> GeneEach;      
        std::vector<std::array<int,2>> GTFAnno_SJs;  
        std::array<int,2> annoExon;
        

        std::ifstream GTFFile; 
        // 打开GTF文件;
        GTFFile.open(GTFFile_name);
        // 判断文件是否打开;
        if (!GTFFile.is_open())	{
            std::cout << "We can't open GTF file !" << GTFFile_name << std::endl;
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
                    now_chr_name = tokens[0]; //
                    early_Exonstrand = now_Exonstrand;
                    now_Exonstrand = tokens[6]; //


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
            std::cout << "***** The GTF file has been read to the end ! *****" << std::endl;
        } 
        else if (GTFFile.fail()) {
            std::cout << "File FALSE !" << std::endl;
        }
        else {std::cout << "unkown reason" << std::endl;}
        //关闭文件;
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
    } else {
        std::cout << "***** No GTF files are put into the program! *****" << std::endl;         
    }
    return GTFAll_Info;
}


struct GTFsj{

    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::array<int,2>>>> mSJs;

    std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> mSJsBE;

    std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> SE;
};


GTFsj get_SJs_SE(std::map<std::string, std::map<std::string, std::vector<std::array<int,2>>>>& Known_Exon){

    GTFsj Known_SJ_SE;

    if (Known_Exon.size() != 0){

        std::map<std::string, std::vector<std::array<int,2>>> Every_transcript_All;
        std::vector<std::array<int,2>> Every_transcript;
        std::string transaction_name;
        std::string chr_name;

        std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::array<int,2>>>> chr_Every_SJs_all;
        std::unordered_map<std::string, std::vector<std::array<int,2>>> Every_SJs_All;
        std::vector<std::array<int,2>> Every_SJs;
        std::array<int,2> SSJs;
        //后三个是存储的SE相关的;
        std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> chr_Every_SE_all;
        std::unordered_map<std::string, std::array<int,2>> Every_SE_All;
        std::array<int,2> Every_SE;
        //存储SJ开头和结尾相关的;
        std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> chr_Every_SJ_begin_end_all;
        std::unordered_map<std::string, std::array<int,2>> Every_SJ_begin_end_all;
        std::array<int,2> Every_SJ_begin_end;

        for (const auto& pair : Known_Exon) {

            chr_name = pair.first;
            Every_transcript_All = pair.second;
            //每个chr清理一下之前的;
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
                        //SJ的开头是exon区间的后一个数, SJ的结尾是下一个exon区间的前一个数减一;
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
            //循环结束一个chr;
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
        //全部循环结束, 将结果存入结构体中;
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


int ifSjSignal(std::string& Fasta_chr, const std::array<int,2>& Sj_Array) {
    int flag = 0;
    std::string signal1 = Fasta_chr.substr(Sj_Array[0]-1, 2);
    std::string signal2 = Fasta_chr.substr(Sj_Array[1]-2, 2);

    std::transform(signal1.begin(), signal1.end(), signal1.begin(), ::toupper);
    std::transform(signal2.begin(), signal2.end(), signal2.begin(), ::toupper);
    //判断索引的信号
    if (signal1 == "GT" && signal2 == "AG"){
        flag = 1;

    } else if (signal1 == "CT" && signal2 == "AC"){
        flag = 2;

    } else if (signal1 == "GC" && signal2 == "AG"){
        flag = 1;

    } else if (signal1 == "CT" && signal2 == "GC"){
        flag = 2;

    } else if (signal1 == "AT" && signal2 == "AC"){
        flag = 1;

    } else if (signal1 == "GT" && signal2 == "AT"){
        flag = 2;
    } else {
        flag = 0;
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
    std::cout << "Optional parameters: [-g <gtf_file>] [-j <SJDistance>] [-n <SJ_support_read_number>] [-d <Graph_distance>] [-t <Thread>] [-m <mode>]" << std::endl;
    std::cout << "Help: [-h <help>]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -s, --sam                   SAM file path. We recommend using absolute paths. If you have a single file, you can directly provide its absolute path. If you have multiple files, you can specify the path to a folder that contains all the sorted SAM files you want to process. (required)" << std::endl;
    std::cout << "  -f, --fasta                 FASTA file path. FASTA file requires the chromosome names to match the GTF file. (required)" << std::endl;
    std::cout << "  -g, --gtf                   input annotation file in GTF format. (optional, Recommendation provided)" << std::endl;
    std::cout << "  -o, --output                output folder path. (required)" << std::endl;
    std::cout << "  -m, --mode                  Sequencing method (default: 0(cDNA), !0(direct RNA))" << std::endl;
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


std::vector<std::string> traverse_sam_file(const std::string& sam_file_path, const std::string& output_path){
    std::vector<std::string> sam_file_vector;
    struct stat sam_stat;
    stat(sam_file_path.c_str(), &sam_stat);

    if (S_ISREG(sam_stat.st_mode)) {
        if (sam_file_path.rfind(".txt") == sam_file_path.length() - 4 || sam_file_path.rfind(".tsv") == sam_file_path.length()) {
            std::cout << "* This is a txt/tsv file! * << " << sam_file_path << std::endl;
            std::ifstream infile(sam_file_path);
            if (!infile) {
                std::cerr << "The file cannot be opened: " << sam_file_path << std::endl;
                std::cerr << "* Not a valid file! *" << strerror(errno) << std::endl;
            }            
            std::string line;
            while (std::getline(infile, line)) {
                if (!line.empty()) {
                    sam_file_vector.push_back(line);
                }
            }
            infile.close();
            if (sam_file_vector.size() == 0) {
                std::cout << "^-^ There are " << 0 << " sam files in total. ^-^" << std::endl;
                std::cerr << "* Not a valid file! *" << strerror(errno) << std::endl;
            } else {
                std::cout << "^-^ There are " << sam_file_vector.size() << " sam files in total. ^-^" << std::endl;
            }
        } else if (sam_file_path.rfind(".sam") == sam_file_path.length() - 4) {
            std::cout << "* Only one sam file is entered! * << " << sam_file_path << std::endl;
            sam_file_vector.push_back(sam_file_path);
        } else {
            std::cerr << "* Not a valid file! *" << strerror(errno) << std::endl;
        }
        
    } else if (S_ISDIR(sam_stat.st_mode)) {
        std::cout << "* A folder was entered! * << " << sam_file_path << std::endl;
        DIR* dir = opendir(sam_file_path.c_str());
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            if (entry->d_name[0] != '.') { 
                std::string fileName = entry->d_name;
                std::cout << "Find file: " << fileName << std::endl;
                if (fileName.rfind(".sam") == fileName.length() - 4) {
                    sam_file_vector.push_back(fileName);
                }
            }
        }
        closedir(dir);
        std::cout << "^-^ There are " << sam_file_vector.size() << " sam files in total. ^-^" << std::endl;
    } else {
        std::cerr << "* Not a valid file or folder! *" << strerror(errno) << std::endl;
    }
    std::string File_explain_path = joinPath(output_path, "file_explain.txt");
    std::ofstream Explain_file(File_explain_path, std::ios::trunc);    
    Explain_file << "File" << '\t' << "File_Path" << '\n';
    for (int i = 0; i < sam_file_vector.size(); i++) {
        std::cout << "SAM File " << i << " : " << sam_file_vector[i] << std::endl;
        Explain_file << i << '\t' << sam_file_vector[i] << '\n';
    }
    Explain_file.close();
    return sam_file_vector;
}


std::streampos findNextLineStart(std::ifstream& file, std::streampos pos) {
    file.seekg(pos);
    std::string line;

    if (std::getline(file, line)) {
        return file.tellg(); 
    }
    return pos; 
}


void processChunk(const std::string& one_sam_file_path, 
                const std::streampos& start, 
                const std::streampos& end, 
                const std::string& output_path, 
                const int file_i, 
                const int SJ_Distance, 
                std::unordered_map<std::string, std::string>& FastaRef,
                std::string& umi_tag_second,
                std::string& barcode_tag_second){

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
    std::ofstream ReadInform(ReadInformPath, std::ios::trunc);


    std::ifstream samfile(one_sam_file_path);
    samfile.seekg(start);

    std::streampos Current_Position = samfile.tellg();
    std::streampos Last_Position = Current_Position;

    while (Current_Position < end) {
        Group_index++;
        each_cluster_informs = get_each_cluster_reads_sc(samfile, Last_Position, end, umi_tag_second, barcode_tag_second);
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
                        << each_cluster_informs.ReadsBarcode[eachRead.first];
            
            SjNumber = 0;
            for (const auto& Sj:eachRead.second) {
                ReadInform << '\t' << Sj[0] << '\t' << Sj[1];

                thisSjNoError = ifSjNoError(each_read_SJs_informs.reads_SJs_left[eachRead.first], each_read_SJs_informs.reads_SJs_right[eachRead.first], SjNumber); 
                if (thisSjNoError == 1) {

                    thisSjSignal = ifSjSignal(FastaRef[chrchr], Sj);
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
                        << each_cluster_informs.ReadsBarcode[eachRead.first] << '\n';         
        } 
    } // 循环结束;
    
    samfile.close();
    ReadInform.close();
    std::cout << "^-^ Thread: " << std::this_thread::get_id() << " has completed processing! ^-^" << std::endl;
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



FileSplit Merge_Read_Small_Files(const std::string& SmallFilePath, const int& SAMFileNumber, const int& nthread) {
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

    for (int file_number = 0; file_number < nthread; file_number++){

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
                                    std::string& umi_tag_first,
                                    std::string& barcode_tag_first) {
    std::vector<std::string> sam_file_vec;
    if (samfilevec.size() == 1) {
        sam_file_vec.push_back(sam_file_path);
    } else {
        for (const auto& eachFile:samfilevec) {
            std::string chunkFilePath = joinPath(sam_file_path, eachFile);
            sam_file_vec.push_back(chunkFilePath);
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

            std::cout << "******* " << "Start processing SAM File " << samFileNumber << " *******" << std::endl;
            std::string chunkFilePath = "sam_" + std::to_string(samFileNumber);
            chunkFilePath = joinPath(outputPath, chunkFilePath);
            // 检测目录是否存在;
            if (!directoryExists(chunkFilePath)) {
                createDirectory(chunkFilePath);      
            }
            struct stat statBuf;
            if (stat(sam_file_vec[samFileNumber].c_str(), &statBuf) != 0){
                std::cout << "We can't get file size !" << sam_file_vec[samFileNumber] << std::endl;
                exit(EXIT_FAILURE);
            }
            std::streampos fileSize = statBuf.st_size;

            std::streampos chunkSize = fileSize / numThreads;

            std::ifstream one_sam(sam_file_vec[samFileNumber]);
            // 判断文件是否打开;
            if (!one_sam.is_open())	{
                std::cout << "We can't open SAM file !" << sam_file_vec[samFileNumber] << std::endl;
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
                        umi_tag_first,
                        barcode_tag_first
                    );
                }));
            }
            for (auto& future : myJobs) {
                future.get(); 
            }
            std::cout << "^-^ [" << sam_file_vec[samFileNumber] << "] All threads are finished generating small files! ^-^" << std::endl;
            
            std::cout << "^-^ Start of merge small files ! ^-^" << std::endl;
            BigBang = Merge_Read_Small_Files(chunkFilePath, samFileNumber, numThreads);
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



struct Line_Split{
    std::string read_file; //0
    std::string read_name; //1
    std::string chr_name; //2
    std::string Group_index; //3
    std::array<int,2> Group_coverage; //4,5
    std::array<int,2> read_coverage; //6,7
    int read_length; // 这个就是追踪文件的时候可以画Coverage; //8
    std::string read_barcode; // 9
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
    std::string chrName; 
    std::array<int,2> GroupCoverage; 
    std::unordered_map<std::string, std::string> GroupReadBarcodes;
    std::unordered_map<std::string, std::string> GroupReadFiles; 
    std::unordered_map<std::string, std::vector<std::array<int,2>>> GroupReadSjs; 
    std::unordered_map<std::string, std::array<int,2>> GroupReadCoverage; 
    std::unordered_map<std::string, std::vector<int>> GroupSigns; 
    std::unordered_map<std::string, std::array<int,2>> GroupSingleExon; 
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
            if (Line_result.read_SJ.size() > 0) {
                groupinformation.GroupReadSjs[Line_result.read_name] = Line_result.read_SJ;
                groupinformation.GroupSigns[Line_result.read_name] = Line_result.read_sj_quality;
            } else {
                groupinformation.GroupSingleExon[Line_result.read_name] = Line_result.read_coverage;
            }
        } else {
            break;
        }
    }
    FinalSamFile.close();
    return groupinformation;
}



struct GroupAnnotation
{
    std::unordered_map<std::string, std::vector<std::array<int,2>>> Group_Annotations;
    std::set<std::string> Group_GeneSet;
};

GroupAnnotation get_group_annotation(std::unordered_map<std::string, std::array<int,2>>& AnnoCoverage, 
                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& AnnoisoSJ,
                                std::array<int,2>& Group_Coverage){

    GroupAnnotation ThisGroupAnnotations;
    std::string extracted;

    for (const auto& EachAnno:AnnoCoverage) {
        if ((EachAnno.second[0] >= Group_Coverage[0]) && (EachAnno.second[1] <= Group_Coverage[1])){

            ThisGroupAnnotations.Group_Annotations[EachAnno.first] = AnnoisoSJ[EachAnno.first];

            size_t pos = EachAnno.first.find("|");
            if (pos != std::string::npos) {
                extracted = EachAnno.first.substr(0, pos);
            }
            ThisGroupAnnotations.Group_GeneSet.insert(extracted);
        }
    }
    return ThisGroupAnnotations;
}


std::unordered_map<std::string, std::array<int,2>> get_group_single_exon_annotation(
                                    std::unordered_map<std::string, std::array<int,2>>& AnnoSE, 
                                    std::array<int,2>& Group_Coverage) {

    std::unordered_map<std::string, std::array<int,2>> ThisGroupSingleExonAnnotations;

    for (const auto& EachAnno:AnnoSE) {
        if ((EachAnno.second[0] >= Group_Coverage[0]) && (EachAnno.second[1] <= Group_Coverage[1])){

            ThisGroupSingleExonAnnotations[EachAnno.first] = EachAnno.second; 
        }
    }
    return ThisGroupSingleExonAnnotations;
}


std::unordered_map<std::string, std::vector<std::string>> get_group_single_exon_reads(
                                    std::unordered_map<std::string, std::array<int,2>>& ThisGroupSingleExonAnnotations,
                                    std::unordered_map<std::string, std::array<int,2>>& GroupSingleExonReads,
                                    int& Edge) {
 
    std::unordered_map<std::string, std::vector<std::string>> SingleReads;

    for (const auto& eachAnno:ThisGroupSingleExonAnnotations) {
        SingleReads[eachAnno.first] = {};
    }

    for (const auto& EachSing:GroupSingleExonReads) {

        for (const auto& EachAnno:ThisGroupSingleExonAnnotations) {
            if ((EachSing.second[1]-EachAnno.second[1]<Edge) && (EachAnno.second[0]-EachSing.second[0]<Edge)) {
                SingleReads[EachAnno.first].push_back(EachSing.first);
                break;
            }
        }
    }
    return SingleReads;
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

    std::set<std::string> thisFileBarcode;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> thisIsoform_File_Barcode2count;
    std::string first_part;
    std::string second_part;


    if (Updated_Files.is_open()) {

        if (singleexonwithreads.size() > 0) {
            for (const auto& eachAnno:singleexonwithreads) {

                if (eachAnno.second.size() > 0) {
                    size_t pos = eachAnno.first.find('|');
                    if (pos != std::string::npos) {

                        first_part = eachAnno.first.substr(0, pos);
                        second_part = eachAnno.first.substr(pos + 1);
                    }
                    {
                        std::unique_lock<std::mutex> lock(updatedGtfMutex);
                        Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "transcript" << '\t' << singleexongroupannotation[eachAnno.first][0] << '\t' << singleexongroupannotation[eachAnno.first][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[eachAnno.first] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\";" << '\n';
                        Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "exon" << '\t' << singleexongroupannotation[eachAnno.first][0] << '\t' << singleexongroupannotation[eachAnno.first][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[eachAnno.first] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\"; exon_number \"" << 0 << "\";" << '\n';
                    }

                    thisIsoform_File_Barcode2count.clear();
                    for (int i = 0; i < FileNo; i++) {
                        std::string fileIndex = std::to_string(i);
                        thisIsoform_File_Barcode2count[fileIndex] = {};
                        for (const auto& eachBar:AllFileBarcode[fileIndex]) {
                            thisIsoform_File_Barcode2count[fileIndex][eachBar] = 0;
                        }
                    }

                    for (const auto& eachRead:eachAnno.second) {
                        std::string thisfile = groupreadfiles[eachRead];
                        std::string thisbarcode = groupreadbarcodes[eachRead];
                        thisIsoform_File_Barcode2count[thisfile][thisbarcode] = thisIsoform_File_Barcode2count[thisfile][thisbarcode] + 1;

                        {
                            std::unique_lock<std::mutex> lock(traceMutex);
                            Trace << eachRead << '\t' << "single_exon" << '\t' << second_part << '\t' << first_part << '\t' << thisbarcode << '\t' << thisfile << '\n'; 
                        }
                    }

                    for (int j = 0; j < FileNo; j++) {
                        std::string fileIndex = std::to_string(j);
                        thisFileBarcode = AllFileBarcode[fileIndex];
                        {
                            std::unique_lock<std::mutex> lock(IsoformMutexes[j]);
                            *(IsoformFilePath[j]) << second_part << '\t' << first_part;
                            for (const auto& eachBar:thisFileBarcode) {
                                *(IsoformFilePath[j]) << '\t' << thisIsoform_File_Barcode2count[fileIndex][eachBar];
                            }
                            *(IsoformFilePath[j]) << '\n';
                        }
                    } 
                }
            }
        }            
    }    
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
    std::unordered_map<std::size_t, std::vector<std::string>> ClustersReads; //用于存储cluster;
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
                                   std::ofstream& Trace){
    ReferenceCluster FSMISMO = {};

    std::vector<std::array<int,2>> each_read_SJs; 

    std::vector<std::string> statis_FSM; 
    std::vector<std::string> statis_ISM; 
    std::vector<int> ISM_flag; 

    std::string FSM_name; 
    std::vector<int> minFSM; 
    int minmin = 0; 

    std::string first_part;
    std::string second_part; 


    for (auto it = ClustersReads.begin(); it != ClustersReads.end(); ++it) {
        each_read_SJs = AllSJs[(it->second)[0]];
        
        statis_FSM.clear(); 
        statis_ISM.clear();
        for (const auto& pair : GroupAnno){
            if (pair.second == each_read_SJs){

                statis_FSM.push_back(pair.first);
            } else {

                if (pair.second.size() > each_read_SJs.size()){
 
                    ISM_flag.clear();
                    
                    for (const auto& element : each_read_SJs) {
                        int it_index = std::find(pair.second.begin(), pair.second.end(), element) - pair.second.begin();
                        if (it_index == pair.second.size()){
                            ISM_flag.push_back(10000);  
                            break; 
                        } else {
                            ISM_flag.push_back(it_index);
                        }
                    }
                    if ((ISM_flag[ISM_flag.size()-1] != 10000) && (ISM_flag[0] != 10000) && (ISM_flag[ISM_flag.size()-1] - ISM_flag[0] == ISM_flag.size() - 1)) {statis_ISM.push_back(pair.first);}
                }
            }
        }

        if (statis_FSM.size() == 1){
            FSM_name = statis_FSM[0];
            FSMISMO.FSM[FSM_name].first = each_read_SJs;
            FSMISMO.FSM[FSM_name].second.insert(FSMISMO.FSM[FSM_name].second.end(), ClustersReads[it->first].begin(), ClustersReads[it->first].end());

            size_t pos = FSM_name.find('|');
            if (pos != std::string::npos){
                first_part = FSM_name.substr(0, pos);
                second_part = FSM_name.substr(pos + 1);
            }
            if (Trace.is_open()){
                for (const auto& EachRead:ClustersReads[it->first]){
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << EachRead << '\t' << "FSM" << '\t' << second_part << '\t' << first_part << '\t' << GroupReadBarcode[EachRead] << '\t' << GroupReadFile[EachRead] << '\n'; 
                }
            }                        

        } else if (statis_FSM.size() > 1){
            
            for (const auto& refx : statis_FSM){
                FSM_name = refx;
                FSMISMO.FSM[FSM_name].first = each_read_SJs;
                FSMISMO.FSM[FSM_name].second = {};
            }
            for (const auto& readx : ClustersReads[it->first]) { 
                minFSM.clear();
                for (const auto& refx : statis_FSM){
                    minFSM.push_back(std::abs(AllReadBE[readx][0] - GroupAnnoBE[refx][0]) + std::abs(AllReadBE[readx][1] - GroupAnnoBE[refx][1]));
                }

                minmin = std::min_element(minFSM.begin(), minFSM.end()) - minFSM.begin();

                FSM_name = statis_FSM[minmin];
                FSMISMO.FSM[FSM_name].second.push_back(readx); 

                size_t pos = FSM_name.find('|');
                if (pos != std::string::npos){
                    first_part = FSM_name.substr(0, pos);
                    second_part = FSM_name.substr(pos + 1);
                }
                {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << readx << '\t' << "FSM" << '\t' << second_part << '\t' << first_part << '\t' << GroupReadBarcode[readx] << '\t' << GroupReadFile[readx] << '\n'; 
                }                
            } 

        } else{ 

            if (statis_ISM.size() == 0){

                FSMISMO.Others[it->first] = ClustersReads[it->first];
            } else{

                FSMISMO.ISM[it->first] = ClustersReads[it->first];

                for (int i = 0; i < statis_ISM.size(); i++){
                    auto isIn = FSMISMO.FSM.find(statis_ISM[i]);
                    if (isIn == FSMISMO.FSM.end()){
                        FSMISMO.FSM[statis_ISM[i]].first = GroupAnno[statis_ISM[i]];
                        FSMISMO.FSM[statis_ISM[i]].second = {};
                    }
                }
            }
        }
    }
    return FSMISMO;
}



struct HighLowClusters{
    std::map<size_t, std::vector<std::string>> HighConClusters;
    std::map<size_t, std::vector<std::string>> LowConClusters;
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

            // 看看满足条件;
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



int IntervalIntersection(std::vector<std::array<int,2>> firstList, std::vector<std::array<int,2>> secondList) {
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


void get_filtered_FSM(std::map<std::size_t, std::vector<std::string>>& LowReads, 
                      std::unordered_map<std::string, std::vector<std::array<int,2>>>& GroupAnno, 
                      std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>>& AllFSM, 
                      std::unordered_map<std::string, std::vector<std::array<int,2>>>& AllSJs, 
                      std::unordered_map<std::string, std::string>& GroupReadFile,
                      std::unordered_map<std::string, std::string>& GroupReadBarcode,
                      std::ofstream& Trace){

    std::unordered_map<std::string, std::unordered_map<std::size_t, std::vector<std::string>>> RecycleAnnoName_value;
    std::set<std::size_t> DeleteItem;

    int RecycleLength = 0;
    int count = 0;
    std::string minAnnoName;
    int dis = 10000;
    std::unordered_map<std::string, int> potential_FSM;
    std::vector<std::array<int,2>> LowReadSJVec;

    std::string first_part;
    std::string second_part;

    for (const auto& each_Readcluster:LowReads){
        LowReadSJVec = AllSJs[each_Readcluster.second[0]];

        potential_FSM.clear();
        for (const auto& each_anno:GroupAnno){
            if (each_anno.second.size() == LowReadSJVec.size()){
                RecycleLength = IntervalMerge(each_anno.second, LowReadSJVec) - IntervalIntersection(each_anno.second, LowReadSJVec);
                if (RecycleLength < 15*LowReadSJVec.size()){

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
                for (const auto& EachRead:each_Readcluster.second){
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << EachRead << '\t' << "FSM" << '\t' << second_part << '\t' << first_part << '\t' << GroupReadBarcode[EachRead] << '\t' << GroupReadFile[EachRead] << '\n'; 
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
            for (const auto& Po:potential_FSM){
                if (Po.second < dis){
                    dis = Po.second;
                    minAnnoName = Po.first;
                }
            }

            size_t pos = minAnnoName.find('|');
            if (pos != std::string::npos){
                first_part = minAnnoName.substr(0, pos);
                second_part = minAnnoName.substr(pos + 1);
            }            
            if (Trace.is_open()) {
                for (const auto& EachRead:each_Readcluster.second){
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << EachRead << '\t' << "FSM" << '\t' << second_part << '\t' << first_part << '\t' << GroupReadBarcode[EachRead] << '\t' << GroupReadFile[EachRead] << '\n'; 
                }
            }   

            if (RecycleAnnoName_value.find(minAnnoName) != RecycleAnnoName_value.end()){
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
        for (const auto& allCluster:newFSM.second){ //每个类;

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
        } else { //后期可以去掉;
            std::cout << "*** That doesn't make sense. There should be a Bug! ***" << std::endl;
        }
    }
}




struct SpliceChainClass
{
    std::unordered_map<std::size_t, std::array<int,2>> ClusterCoverage;
    std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>> FSM;
    std::map<std::size_t, std::vector<std::string>> ISM;
    std::map<size_t, std::vector<std::string>> HighConClusters;
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
                                std::ofstream& traceFilePath,
                                const int& Sj_Support_Number) {
    SpliceChainClass SCC;
    ReferenceCluster FsmIsmOthers;
    HighLowClusters Others2HighLow;
    std::unordered_map<std::size_t, std::vector<std::string>> groupCluster = classifyReadsVec(groupreadsjs);

    SCC.ClusterCoverage = get_every_cluster_begin_end(groupCluster, groupreadcoverage);
    FsmIsmOthers = get_FSM_and_others_sc(groupCluster, groupannotations, AnnoCoverage, groupreadcoverage, groupreadsjs, groupreadfiles, groupreadbarcodes, traceFilePath);
    Others2HighLow = get_HighLow_clusters(FsmIsmOthers.Others, groupreadsjs, groupreadsigns, Sj_Support_Number);
    get_filtered_FSM(Others2HighLow.LowConClusters, groupannotations, FsmIsmOthers.FSM, groupreadsjs, groupreadfiles, groupreadbarcodes, traceFilePath);

    SCC.FSM = FsmIsmOthers.FSM;
    SCC.ISM = FsmIsmOthers.ISM;
    SCC.HighConClusters = Others2HighLow.HighConClusters;
    SCC.HighStrand = Others2HighLow.HighStrand;

    return SCC;
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
        max_node = find_max_neighborhood_node(subg, adj);
        std::set_difference(cand.begin(), cand.end(),
                            adj[max_node].begin(), adj[max_node].end(),
                            std::inserter(ext_u, ext_u.begin()));

        while (true)
        {
            if (!ext_u.empty()) 
            {
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

    //函数用到的中间变量;
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
                        //查看是否满足条件;
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
            if (iter != Node_Results.TrueNodeSet.end()){
                Node_Results.TrueNodeSet.erase(eachFalse);
            }
        }
    }
    return Node_Results;
}

// 生成带有文件区别的这三个的数量;
struct Solvent_sc {
    std::map<std::string, std::vector<std::string>> File_FSM;
    std::map<std::size_t, std::vector<std::string>> File_ISM;
    std::map<std::size_t, std::vector<std::string>> File_HighConClusters;
};

Solvent_sc get_Solvent_FsmIsmHigh(SpliceChainClass& FsmIsmHigh, int& FileNumber, GroupInformation& groupinformations) {
    Solvent_sc FileSpliceChains;
    std::string thisFileString = std::to_string(FileNumber);
    // 拆分FSM;
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
        // 对每个HighConCluster循环;
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
    return FileSpliceChains;
}

std::string get_co_sj_gene(std::vector<std::array<int,2>>& novel_isoform, 
                            std::set<std::string>& gene_set,
                            std::map<std::string, std::vector<std::string>>& gene2tx,
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
        TempGene.push_back(A_Gene);
        ALL_tx = gene2tx[A_Gene];
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

OutputInformation Write_Detection_Transcript2gtf_MultiFiles(std::ofstream& Updated_Files, std::ofstream& Trace,
                                                DetectionResults& noderesults, DistanceInform& Disinform, 
                                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupannotations, 
                                                std::map<std::string, std::vector<std::array<int,2>>>& Annoexon, 
                                                std::unordered_map<std::size_t, std::array<int, 2>>& groupreadcoverage, 
                                                std::set<std::string>& groupallgene, 
                                                std::map<std::string, std::array<int,2>>& Annogenecovergae, 
                                                std::unordered_map<std::string, std::string>& GTF_Transcript_Strand, 
                                                std::unordered_map<size_t, std::string>& High_Strand, 
                                                std::string& chrname, std::string& group_size, 
                                                std::map<size_t, std::vector<std::string>>& HighClusters,
                                                std::unordered_map<std::string, std::string>& groupreadfiles,
                                                std::unordered_map<std::string, std::string>& groupreadbarcodes,
                                                std::map<std::string, std::set<std::string>>& allfilebarcodeset,
                                                std::unordered_map<std::string, std::array<int,2>>& SingleExonAnno,
                                                std::map<std::string, std::array<int,2>>& gtf_gene,
                                                std::map<std::string, std::vector<std::string>>& gtf_gene2tx,
                                                std::map<std::string, std::string>& gtf_gene_std){
    OutputInformation FinalAnnotations;
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
        // 第一次
        for (const auto aaa:noderesults.TrueNodeSet){
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
                BE = groupreadcoverage[nameya];
                
                TempGene.clear();
                if (groupallgene.size() != 0) {
                    TempGeneName = get_co_sj_gene(itssj, groupallgene, gtf_gene2tx, groupannotations);

                    if (TempGeneName.size() != 0) {
                        first_part = TempGeneName;
                    } else {
                        for (const auto& every_gene:groupallgene) {
                            if (Annogenecovergae[every_gene][0] <= BE[1] && BE[0] <= Annogenecovergae[every_gene][1]){
                                TempGene.push_back(every_gene);
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

                    TempDist.clear();
                    for (const auto& every_gene_tx:SingleExonAnno) {
                        if (every_gene_tx.second[0] <= BE[1] && BE[0] <= every_gene_tx.second[1]) {
                            size_t pos = every_gene_tx.first.find('|');
                            if (pos != std::string::npos) {
                                // 根据 "|" 将字符串分割为两部分
                                first_part = every_gene_tx.first.substr(0, pos);
                            }
                            TempGene.push_back(first_part);
                            int a = abs(BE[1] - every_gene_tx.second[1]);
                            int b = abs(every_gene_tx.second[0] - BE[0]);
                            TempDist.push_back((a>b)?a:b);
                        }
                    }
                    if (TempDist.size() > 0) {
                        auto min_it = std::min_element(TempDist.begin(), TempDist.end());
                        first_part = TempGene[std::distance(TempDist.begin(), min_it)];
                    } else {
                        first_part = "NA";
                    }
                }

                FinalAnnotations.transcript2gene[itsname] = first_part;
                {   
                    if (first_part != "NA") {
                        second_part = gtf_gene_std[first_part];
                    } else {
                        second_part = High_Strand[nameya];
                    }
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
                //输出追踪文件;
                if (Trace.is_open()){
                    for (const auto& EachRead:HighClusters[nameya]){
                        std::unique_lock<std::mutex> lock(traceMutex);
                        Trace << EachRead << '\t' << "novel_isoform" << '\t' << itsname << '\t' << first_part << '\t' << groupreadbarcodes[EachRead] << '\t' << groupreadfiles[EachRead] << '\n'; 
                    }
                }
            }
        }
    }
    return FinalAnnotations;
}

std::unordered_map<std::string, std::unordered_map<std::string, double>>  get_transcript_init(
                                            DetectionResults& noderesults, 
                                            DistanceInform& Disinform, 
                                            std::string& chrname, 
                                            std::string& group_size, 
                                            int& FileNumber, 
                                            std::unordered_map<std::string, std::string>& groupreadfiles,
                                            std::unordered_map<std::string, std::string>& groupreadbarcodes,
                                            Solvent_sc& filesolvent, 
                                            std::map<std::string, std::set<std::string>>& allfilebarcodeset){
    std::string fileIndex = std::to_string(FileNumber);
    std::unordered_map<std::string, std::unordered_map<std::string, double>> File_TranscriptNumber;
    std::string itsname;
    std::string itsbarcode;
    int novel_count = 0;
    std::size_t nameya;
    // 初始化;
    for (const auto&eachBar:allfilebarcodeset[fileIndex]) {
        File_TranscriptNumber[eachBar] = {};
    }
    // 第一次
    for (const auto aaa:noderesults.TrueNodeSet){
        if (aaa < Disinform.Index2Anno.size()){
            //有注释的部分;
            itsname = Disinform.Index2Anno[aaa];
            auto is1 = filesolvent.File_FSM.find(itsname);
            if (is1 != filesolvent.File_FSM.end()) {
                for (const auto& eachRead:filesolvent.File_FSM[itsname]) {
                    itsbarcode = groupreadbarcodes[eachRead];
                    auto is2 = File_TranscriptNumber[itsbarcode].find(itsname);
                    if (is2 != File_TranscriptNumber[itsbarcode].end()) {
                        File_TranscriptNumber[itsbarcode][itsname] = File_TranscriptNumber[itsbarcode][itsname] + 1;
                    } else {
                        File_TranscriptNumber[itsbarcode][itsname] = 1;
                    }
                }
            } else {
                for (const auto& eachBar:allfilebarcodeset[fileIndex]) {
                    File_TranscriptNumber[eachBar][itsname] = 0;
                }
            }

        } else {
            novel_count = novel_count + 1;
            //没有注释的部分;
            itsname = chrname + "-novel-" + group_size + "-" + std::to_string(aaa) + "-" + std::to_string(novel_count);
            Disinform.Index2novelname[aaa] = itsname;
            nameya = Disinform.Index2hashname[aaa];
        
            auto is1 = filesolvent.File_HighConClusters.find(nameya);
            if (is1 != filesolvent.File_HighConClusters.end()) {
                for (const auto& eachRead:filesolvent.File_HighConClusters[nameya]) {
                    itsbarcode = groupreadbarcodes[eachRead];
                    auto is2 = File_TranscriptNumber[itsbarcode].find(itsname);
                    if (is2 != File_TranscriptNumber[itsbarcode].end()) {
                        File_TranscriptNumber[itsbarcode][itsname] = File_TranscriptNumber[itsbarcode][itsname] + 1;
                    } else {
                        File_TranscriptNumber[itsbarcode][itsname] = 1;
                    }
                }
            } else {
                for (const auto& eachBar:allfilebarcodeset[fileIndex]) {
                    File_TranscriptNumber[eachBar][itsname] = 0;
                }
            }        
        }
    }
    return File_TranscriptNumber;
}



struct Barcode_sc {
    std::unordered_map<std::string, std::map<std::string, int>> Barcode_FSM;
    std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>> Barcode_ISM;
    std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>> Barcode_HighConClusters;
};

Barcode_sc get_Barcode_FsmIsmHigh(std::map<std::string, std::vector<std::string>>& thisfilefsm,
                            std::map<std::size_t, std::vector<std::string>>& thisfileism,
                            std::map<std::size_t, std::vector<std::string>>& thisfileHC,
                            std::unordered_map<std::string, std::string>& groupreadbarcodes) {
    // 输出;
    Barcode_sc barcodesc;
    std::unordered_map<std::string, std::map<std::string, int>> barcodefsm;
    std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>> barcodeism;
    std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>> barcodehc;
    // FSM
    if (thisfilefsm.size() != 0) {
        for (const auto& eachFSM:thisfilefsm) {
            std::string FSMname = eachFSM.first;
            for (const auto& FSMread:eachFSM.second) {
                std::string thisbarcode = groupreadbarcodes[FSMread];
                auto is1 = barcodefsm.find(thisbarcode);
                if (is1 != barcodefsm.end()) {
                    auto is2 = barcodefsm[thisbarcode].find(FSMname);
                    if (is2 != barcodefsm[thisbarcode].end()) {
                        barcodefsm[thisbarcode][FSMname] = barcodefsm[thisbarcode][FSMname] + 1;
                    } else {
                        barcodefsm[thisbarcode][FSMname] = 1;
                    }
                } else {
                    barcodefsm[thisbarcode] = {};
                    barcodefsm[thisbarcode][FSMname] = 1;
                }
            }
        }
    }
    // ISM
    if (thisfileism.size() != 0) {
        for (const auto& eachISM:thisfileism) {
            std::size_t ISMname = eachISM.first;
            for (const auto& ISMread:eachISM.second) {
                std::string thisbarcode = groupreadbarcodes[ISMread];
                auto is1 = barcodeism.find(thisbarcode);
                if (is1 != barcodeism.end()) {
                    auto is2 = barcodeism[thisbarcode].find(ISMname);
                    if (is2 != barcodeism[thisbarcode].end()) {
                        barcodeism[thisbarcode][ISMname].push_back(ISMread);
                    } else {
                        barcodeism[thisbarcode][ISMname] = {};
                        barcodeism[thisbarcode][ISMname].push_back(ISMread);
                    }
                } else {
                    barcodeism[thisbarcode] = {};
                    barcodeism[thisbarcode][ISMname] = {};
                    barcodeism[thisbarcode][ISMname].push_back(ISMread);
                }
            }
        }
    }

    if (thisfileHC.size() != 0) {
        for (const auto& eachHC:thisfileHC) {
            std::size_t HCname = eachHC.first;
            for (const auto& HCread:eachHC.second) {
                std::string thisbarcode = groupreadbarcodes[HCread];
                auto is1 = barcodehc.find(thisbarcode);
                if (is1 != barcodehc.end()) {
                    auto is2 = barcodehc[thisbarcode].find(HCname);
                    if (is2 != barcodehc[thisbarcode].end()) {
                        barcodehc[thisbarcode][HCname].push_back(HCread);
                    } else {
                        barcodehc[thisbarcode][HCname] = {};
                        barcodehc[thisbarcode][HCname].push_back(HCread);
                    }
                } else {
                    barcodehc[thisbarcode] = {};
                    barcodehc[thisbarcode][HCname] = {};
                    barcodehc[thisbarcode][HCname].push_back(HCread);
                }
            }
        }
    }    
    barcodesc.Barcode_FSM = barcodefsm;
    barcodesc.Barcode_ISM = barcodeism;
    barcodesc.Barcode_HighConClusters = barcodehc;
    return barcodesc;
}




struct FileBarcodeFalseNode
{
    std::set<int> ThisFlaseNode;
    std::unordered_map<int, std::vector<std::string>> FasleNode2Count;
};
FileBarcodeFalseNode get_File_Barcode_False_node (std::set<int>& FalseNodeSet, 
                          std::unordered_map<int, std::size_t>& Node2hashname, 
                          std::unordered_map<std::string, std::map<std::size_t, std::vector<std::string>>>& barcodehc, 
                          const std::string& Barcode) {
    FileBarcodeFalseNode ThisFile_FalseNumber;
    std::map<size_t, std::vector<std::string>> this_file_HC;
    auto it = barcodehc.find(Barcode);
    if (it == barcodehc.end()) {
        return ThisFile_FalseNumber;
    } else {
        this_file_HC = barcodehc[Barcode];
    }
    for (const auto& node:FalseNodeSet) {
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

IndicateFire Quantification_initialization_Barcodes (std::map<std::size_t, std::vector<std::string>>& groupISM, 
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
                                                std::ofstream& Trace) {
    IndicateFire OutputResults;
    std::vector<std::string> Order_Transcript_Name;
    Eigen::MatrixXd known_ISM_matrix;
    std::vector<std::array<int,2>> AreadSJs;
    
    int RowIndex = -1;
    int ColIndex = -1;
    int flag = 0;

    for (const auto& eachTransctipt:FinallyAnnotations.Transcript_Annotations) {
        Order_Transcript_Name.push_back(eachTransctipt.first);
    }
    Eigen::RowVectorXd EachClusterRowVector(Order_Transcript_Name.size());
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

    if (groupISM.size() != 0){
        RowIndex = -1;
        for (const auto& eachISM:groupISM){

            AreadSJs = groupreadsjs[eachISM.second[0]];
            EachClusterRowVector.setZero();
            ColIndex = -1;
            ism_gtf_name.clear();
            first_part_set.clear();
            second_part_set.clear();            

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
                    }
                }
                first_part = concatenateSet(first_part_set);
                second_part = concatenateSet(second_part_set);
                if (Trace.is_open()){
                    for (const auto& EachRead:eachISM.second){
                        std::unique_lock<std::mutex> lock(traceMutex);
                        Trace << EachRead << '\t' << "ISM" << '\t' << second_part << '\t' << first_part << '\t' << groupreadbarcodes[EachRead] << '\t' << groupreadfiles[EachRead] << '\n'; 
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
    Eigen::RowVectorXd FalseNode_RowVector(Order_Transcript_Name.size());
    RowIndex = -1;
    int countC = 0;
    
    if (MergeFalseSet.size() != 0){
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
                    } else
                    {
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
            //开始输出HighCon的追踪文件;
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
            if (Trace.is_open()){
                for(const auto& EachRead:HighClusters[Disinform.Index2hashname[FalseNode]]){
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << EachRead << '\t' << "approximate" << '\t' << second_part << '\t' << first_part << '\t' << groupreadbarcodes[EachRead] << '\t' << groupreadfiles[EachRead] << '\n';
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
    OutputResults.Order_Transcript_Name_Vector = Order_Transcript_Name;
    return OutputResults;
}
 
void EM_Alg_Barcodes (std::unordered_map<std::string, std::unordered_map<std::string, double>>& filetranscriptnumber, 
                        IndicateFire& Indicate_Number, 
                        const std::string& Barcode) {

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
            P1 = AnnoN / AnnoN.sum();
            sum_abs_diff = (P1 - P_Col_init0).cwiseAbs().sum();
            P_Col_init0 = P1;
            CountCyc = CountCyc + 1;
        } while ((sum_abs_diff > 5e-2) || (CountCyc <= 10));

        std::string its_name;
        for (int i = 0; i < AnnoN.rows(); i++){
            its_name = Indicate_Number.Order_Transcript_Name_Vector[i];
            filetranscriptnumber[Barcode][its_name] = filetranscriptnumber[Barcode][its_name] + AnnoN(i,0);
        }
    }
}   


void write_transcript_file(std::unordered_map<std::string, std::unordered_map<std::string, double>>& filetranscriptnumber,
                            OutputInformation& FinallyAnnotations, std::set<std::string>& barcodeset, int& filenumber,
                            std::vector<std::unique_ptr<std::ofstream>>& AllIsoformFilePath,
                            std::vector<std::mutex>& AllIsoformMutexes, int& rc_threshold) {

    std::string transcript_name;
    std::string gene_name;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> transcriptbarcodenumber;
    if (FinallyAnnotations.Transcript_Annotations.size() > 0) {
        for (const auto& eachAnno:FinallyAnnotations.Transcript_Annotations) {
            transcriptbarcodenumber[eachAnno.first] = {};
        }

        for (const auto& eachAnno:FinallyAnnotations.Transcript_Annotations) {
            std::string AnnoName = eachAnno.first;
            for (const auto& eachBar:filetranscriptnumber) {
                transcriptbarcodenumber[AnnoName][eachBar.first] = filetranscriptnumber[eachBar.first][AnnoName];
            }
        }

        for (const auto& eachAnno:transcriptbarcodenumber) {
            size_t pos = eachAnno.first.find('|');
            if (pos != std::string::npos) {
                transcript_name = eachAnno.first.substr(pos + 1);
            } else {
                transcript_name = eachAnno.first;
            }
            gene_name = FinallyAnnotations.transcript2gene[eachAnno.first];
            {
                std::unique_lock<std::mutex> lock(AllIsoformMutexes[filenumber]);
                *(AllIsoformFilePath[filenumber]) << transcript_name << '\t' << gene_name;
                for (const auto& eachBar:barcodeset) {
                    double Counts = transcriptbarcodenumber[eachAnno.first][eachBar];
                    if ( Counts < rc_threshold ) {
                        *(AllIsoformFilePath[filenumber]) << '\t' << 0;
                    } else {
                        *(AllIsoformFilePath[filenumber]) << '\t' << Counts;
                    }
                    
                }
                *(AllIsoformFilePath[filenumber]) << '\n';
            }

        }
    }
}



void DetectQuant(GroupAnnotation& groupanno, 
                SpliceChainClass& splicechainclass, 
                GroupInformation& groupinform, 
                const int& groupdistance,
                GTF& gtfexon, GTFsj& gtfsjs,
                std::ofstream& gtfFilePath,
                std::ofstream& traceFilePath, 
                int& FileNo,
                std::map<std::string, std::set<std::string>>& AllFile_BarcodeSet,
                std::vector<std::unique_ptr<std::ofstream>>& IsoformFilePath,
                std::vector<std::mutex>& IsoformMutexes,
                std::unordered_map<std::string, std::array<int,2>>& singleexonanno,
                int& RC_threshold) {
    
    std::string chrchr = groupinform.chrName;
    std::string groupnumber = groupinform.GroupIndex;
    DistanceInform DMatrix_GraphNode = get_distance_matrix(groupanno.Group_Annotations, 
                                                           splicechainclass.HighConClusters, 
                                                           groupinform.GroupReadSjs, 
                                                           groupdistance, splicechainclass.FSM);

    if (DMatrix_GraphNode.NodeDistanceSet.size() != 0) {
        std::vector<std::vector<int>> CliquesVector = Find_Maximal_Cliques(DMatrix_GraphNode.NodeDistanceSet);
        DetectionResults NodeResults = Transcript_Detection(CliquesVector, 
                                                            DMatrix_GraphNode.Index2Anno, 
                                                            DMatrix_GraphNode.Index2Unknown, 
                                                            DMatrix_GraphNode.Index2Count);
        CliquesVector.clear();

        OutputInformation Finally_Annotations = Write_Detection_Transcript2gtf_MultiFiles(gtfFilePath, traceFilePath, 
                                                NodeResults, DMatrix_GraphNode, 
                                                gtfsjs.mSJs[chrchr], gtfexon.GTF_transcript[chrchr], 
                                                splicechainclass.ClusterCoverage, groupanno.Group_GeneSet, 
                                                gtfexon.GTF_gene[chrchr], 
                                                gtfexon.GTF_transcript_strand[chrchr], 
                                                splicechainclass.HighStrand, chrchr, groupnumber, 
                                                splicechainclass.HighConClusters,
                                                groupinform.GroupReadFiles,
                                                groupinform.GroupReadBarcodes,
                                                AllFile_BarcodeSet,
                                                singleexonanno,
                                                gtfexon.GTF_gene[chrchr], gtfexon.GTF_gene2transcript[chrchr],
                                                gtfexon.GTF_gene_strand[chrchr]); 
        splicechainclass.HighStrand.clear();
        groupanno.Group_Annotations.clear();
        groupanno.Group_GeneSet.clear();
        if (FileNo > 1) {
            for (int k = 0; k < FileNo; k++) {
                Solvent_sc SpliceChainSolvent = get_Solvent_FsmIsmHigh(splicechainclass, k, groupinform);
                std::unordered_map<std::string, std::unordered_map<std::string, double>> File_k_TranscriptNumber = get_transcript_init(
                    NodeResults, DMatrix_GraphNode, chrchr, groupnumber, k, groupinform.GroupReadFiles, groupinform.GroupReadBarcodes,
                    SpliceChainSolvent, AllFile_BarcodeSet
                );
                Barcode_sc ThisFileBarcodeSpliceChain = get_Barcode_FsmIsmHigh(SpliceChainSolvent.File_FSM,
                                                                    SpliceChainSolvent.File_ISM,
                                                                    SpliceChainSolvent.File_HighConClusters,
                                                                    groupinform.GroupReadBarcodes);
                for (const auto& ThisBarcode:AllFile_BarcodeSet[std::to_string(k)]) {

                    FileBarcodeFalseNode This_File_False_Node = get_File_Barcode_False_node(NodeResults.FalseNodeSet,
                                                                            DMatrix_GraphNode.Index2hashname,
                                                                            ThisFileBarcodeSpliceChain.Barcode_HighConClusters,
                                                                            ThisBarcode);

                    IndicateFire InitFirefly = Quantification_initialization_Barcodes(
                                                                ThisFileBarcodeSpliceChain.Barcode_ISM[ThisBarcode], 
                                                                Finally_Annotations, NodeResults.TrueNodeSet,
                                                                This_File_False_Node.ThisFlaseNode, 
                                                                This_File_False_Node.FasleNode2Count,
                                                                DMatrix_GraphNode, groupnumber, 
                                                                ThisFileBarcodeSpliceChain.Barcode_HighConClusters[ThisBarcode], 
                                                                groupinform.GroupReadSjs, 
                                                                groupinform.GroupReadFiles,
                                                                groupinform.GroupReadBarcodes,
                                                                traceFilePath);

                    EM_Alg_Barcodes(File_k_TranscriptNumber, InitFirefly, ThisBarcode); 
                    // std::cout << "(=^_^=) Group " << groupinform.GroupIndex << " in file " << k << " completed quantification! (=^_^=) " << std::endl;                                         
                    ThisFileBarcodeSpliceChain.Barcode_FSM[ThisBarcode].clear();
                    ThisFileBarcodeSpliceChain.Barcode_ISM[ThisBarcode].clear();
                    ThisFileBarcodeSpliceChain.Barcode_HighConClusters[ThisBarcode].clear();
                }
                write_transcript_file(File_k_TranscriptNumber, Finally_Annotations, AllFile_BarcodeSet[std::to_string(k)], k, IsoformFilePath, IsoformMutexes, RC_threshold);
            }
            std::cout << "end one group!" << std::endl;

        } else {
            int k = 0;
            Solvent_sc SpliceChainSolvent = get_Solvent_FsmIsmHigh(splicechainclass, k, groupinform);
            std::unordered_map<std::string, std::unordered_map<std::string, double>> File_k_TranscriptNumber = get_transcript_init(
                NodeResults, DMatrix_GraphNode, chrchr, groupnumber, k, groupinform.GroupReadFiles, groupinform.GroupReadBarcodes,
                SpliceChainSolvent, AllFile_BarcodeSet);

            Barcode_sc ThisFileBarcodeSpliceChain = get_Barcode_FsmIsmHigh(SpliceChainSolvent.File_FSM,
                                                                SpliceChainSolvent.File_ISM,
                                                                SpliceChainSolvent.File_HighConClusters,
                                                                groupinform.GroupReadBarcodes);

            for (const auto& ThisBarcode:AllFile_BarcodeSet["0"]) {

                FileBarcodeFalseNode This_File_False_Node = get_File_Barcode_False_node(NodeResults.FalseNodeSet,
                                                                        DMatrix_GraphNode.Index2hashname,
                                                                        ThisFileBarcodeSpliceChain.Barcode_HighConClusters,
                                                                        ThisBarcode); 

                IndicateFire InitFirefly = Quantification_initialization_Barcodes (
                                                            ThisFileBarcodeSpliceChain.Barcode_ISM[ThisBarcode], 
                                                            Finally_Annotations, NodeResults.TrueNodeSet,
                                                            This_File_False_Node.ThisFlaseNode, 
                                                            This_File_False_Node.FasleNode2Count,
                                                            DMatrix_GraphNode, groupnumber, 
                                                            ThisFileBarcodeSpliceChain.Barcode_HighConClusters[ThisBarcode], 
                                                            groupinform.GroupReadSjs, 
                                                            groupinform.GroupReadFiles,
                                                            groupinform.GroupReadBarcodes,
                                                            traceFilePath);

                EM_Alg_Barcodes(File_k_TranscriptNumber, InitFirefly, ThisBarcode);
                // std::cout << "(=^_^=) Group " << groupinform.GroupIndex << " reads numbers " << groupinform.GroupReadSjs.size()  << " completed quantification! (=^_^=) " << std::endl;               
            }
            int fileIndex = 0;
            write_transcript_file(File_k_TranscriptNumber, Finally_Annotations, AllFile_BarcodeSet["0"], fileIndex, IsoformFilePath, IsoformMutexes, RC_threshold); 
        }
    }
}


void processGroup(std::streampos& start, std::streampos& end, 
                  const std::string& sam_file_path, 
                  const int& Sj_supportReadNumber, GTF& gtf_full, 
                  GTFsj& gtf_splice, const int& GraphDis,
                  std::ofstream& updatedgtffile, std::ofstream& tracefile,
                  std::vector<std::unique_ptr<std::ofstream>>& isoformfilePath,
                  std::vector<std::mutex>& isoformMutexes,
                  int& fileno, std::string& outputPath,
                  std::map<std::string, std::set<std::string>>& filebarcodeSet,
                  int& singleEdge, int& ReadCount_threshold) {

    GroupInformation group_information = knowGroupInformation(start, end, sam_file_path, Sj_supportReadNumber);
    std::string chrchr = group_information.chrName;
    std::array<int,2> groupcoverage = group_information.GroupCoverage;
    std::unordered_map<std::string, std::array<int,2>> single_exon_group_annotation;
    // std::cout << "[{(>_<)]} Group sj reads numbers are " << group_information.GroupReadSjs.size() << " ! ^-^ Group single exon reads numbers are " << group_information.GroupSingleExon.size() << " ! [{(>_<)]}" << std::endl;
    
    if (group_information.GroupSingleExon.size() > 0) {
        single_exon_group_annotation = get_group_single_exon_annotation(gtf_splice.SE[chrchr], groupcoverage);
        std::unordered_map<std::string, std::vector<std::string>> single_exon_with_reads = get_group_single_exon_reads(single_exon_group_annotation, group_information.GroupSingleExon, singleEdge);
        write_single_exon_gtf_trace_sc(fileno, group_information.GroupReadFiles,
                                    single_exon_with_reads, single_exon_group_annotation, 
                                    updatedgtffile, tracefile, isoformfilePath, isoformMutexes,
                                    gtf_full.GTF_transcript_strand[chrchr],
                                    group_information.GroupReadBarcodes, chrchr,
                                    filebarcodeSet);
    }
    std::cout << group_information.GroupIndex << " is end! Single reads size is " << group_information.GroupSingleExon.size() << "!"<< std::endl;
    group_information.GroupSingleExon.clear();

    if (group_information.GroupReadSjs.size() > 0) {
        GroupAnnotation group_annotation = get_group_annotation(gtf_splice.mSJsBE[chrchr], gtf_splice.mSJs[chrchr], groupcoverage);
        SpliceChainClass spliceclass = generate_splice_chain_class(group_information.GroupReadSjs, 
                                                group_information.GroupReadCoverage, 
                                                group_annotation.Group_Annotations, 
                                                group_information.GroupSigns, 
                                                gtf_splice.mSJsBE[chrchr], 
                                                group_information.GroupReadFiles, 
                                                group_information.GroupReadBarcodes, 
                                                tracefile, Sj_supportReadNumber);
        group_information.GroupSigns.clear();
        group_information.GroupReadCoverage.clear();
        DetectQuant(group_annotation, spliceclass, group_information, 
                    GraphDis, gtf_full, gtf_splice, 
                    updatedgtffile, tracefile, fileno,
                    filebarcodeSet, isoformfilePath, isoformMutexes,
                    single_exon_group_annotation, ReadCount_threshold); 
    }
    std::cout << group_information.GroupIndex << " is end! reads size is " << group_information.GroupReadSjs.size() << "!"<< std::endl;
}

// void Remove_cache_file () {
    
// }

int main(int argc, char* argv[])
{   

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
    int single_exon_edge = 60;
    std::string umi = "U8:Z:";
    std::string barcode = "BC:Z:";

    while ((c = getopt_long(argc, argv, "s:f:g:o:m:j:n:e:d:t:u:b:r:h", long_options, &option_index)) != -1) {
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
                umi = std::stoi(optarg);   
                break;
            case 'b':
                barcode = std::stoi(optarg);  
                break;
            case 'r':
                Read_count = std::stoi(optarg);   
                break;            
            case 'h':
                print_usage(argv[0]);
                exit(EXIT_FAILURE);              
            case '?':
                // 如果getopt_long返回'?'，说明有错误的参数
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

    std::cout << "*****" << std::endl;
    std::cout << "Input file: " << samfile_name << std::endl;
    std::cout << "FASTA file: " << fastafile_name << std::endl;
    std::cout << "GTF file: " << gtffile_name << std::endl;
    std::cout << "Output file: " << output_file_name << std::endl;
    std::cout << "SJ Distance: " << SJDistance << std::endl;
    std::cout << "SJ support read number: " << SJ_support_read_number << std::endl;
    std::cout << "Single exon boundary : " << single_exon_edge << std::endl;
    std::cout << "Graph distance: " << Graph_distance << std::endl;
    std::cout << "Thread: " << Thread << std::endl;
    std::cout << "UMI tag: " << umi << std::endl;
    std::cout << "Barcode tag: " << barcode << std::endl;
    std::cout << "Output min read count: " << Read_count << std::endl;
    std::cout << "*****" << std::endl;

    std::unordered_map<std::string, std::string> Fasta = Read_fasta_file(fastafile_name);

    FileSplit BroCOLIfile = thread_all_read_sam_files(samfile_name, sam_file_vec, Thread, output_file_name, SJDistance, Fasta, umi, barcode);
    Fasta.clear();
    std::vector<std::unique_ptr<std::ofstream>> BroCOLIQuantfile = write_quantification_files(output_file_name, BroCOLIfile.file_barcodeSet);
    std::vector<std::mutex> TranscriptMutexes(BroCOLIQuantfile.size());

    GTF GTF_full = get_gtf_annotation(gtffile_name);

    GTFsj GTF_Splice = get_SJs_SE(GTF_full.GTF_transcript);

    std::vector<std::size_t> Group_idx = sort_indexes_e(BroCOLIfile.group_reads_number);
  
    std::cout << "BroCOLIfile.group_reads_number: " << BroCOLIfile.group_reads_number.size() << std::endl;
    std::cout << "BroCOLIfile.reads_pointer " << BroCOLIfile.reads_pointer.size() << std::endl;
    BroCOLIfile.chr_coverage.clear();
    BroCOLIfile.coverage2pos.clear();
    BroCOLIfile.group_reads_number.clear();

    ThreadPool BroCOLIpool(Thread);

    std::vector<std::future<void>> futures;
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

    gtf_file.close();
    trace_file.close();
    return 0;
}




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
#include <sys/stat.h>   // For stat and mkdir
#include <sys/types.h>  // For stat
#include <dirent.h>
#include <thread>
#include <mutex>
#include <iomanip>
#include <condition_variable>
#include <queue>
#include <functional>
#include <atomic>
#include "thread_pool_executor.hpp"

//需要写出四个文件; 这也是输出文件向量中文件的顺序;
std::mutex updatedGtfMutex;
std::mutex traceMutex;
std::mutex bigMutex;

//1.
//定义返回的tokens和read_name;
//分别对应返回的每一行的value和key;
//返回值定义在了一个结构体中;
struct Split_Result{
    std::string read_name;
    std::vector<std::string> tokens;
    int read_length;
    std::string read_barcode;
    std::string read_umi;
};

/* 读取SAM文件时,
   切割文件中每一行表示的字符串并存入vector;
   参数1表示要切割的字符串；
   参数2表示根据切割的字符；*/
Split_Result get_string_split_sc(const std::string& s, char delimiter) {
    // 初始化返回值结构体;
    Split_Result read_result = {};
    //字符串流;
    std::istringstream iss(s);
    //定义字符串流中的每一段字符串;
    std::string token;
    /*  第一列：reads的名称
        第二列：比对到的参考基因组的状况如正负链
        第三列：比对到的基因的名称
        第四列：比对起始的位置
        第五列：比对的质量分数MAPQ
        第六列：CIGAR值
        第十列: read的长度 */
    //定义一个number, 将有用的信息保存到
    int number = 0;
    while (std::getline(iss, token, delimiter)) {
        number++;
        //read name保存在一个string中;
        //除了reads name以外的信息保存到一个vector中;
        if (number == 1){
            read_result.read_name = token;
        } else if (1 < number && number <= 6) {
        //在vector容器末尾插入字符串元素;
            read_result.tokens.push_back(token);
        } else if (number == 10) {
            read_result.read_length = token.size();
        } else if (token.size() > 10 && token.compare(0,5,"U8:Z:") == 0) {
            read_result.read_umi = token;
        } else if (token.size() > 10 && token.compare(0,5,"BC:Z:") == 0) {
            read_result.read_barcode = token;
        }
        // 修正一下read名称;
        if (read_result.read_umi.size() > 0 && read_result.read_barcode.size() > 0) {
            read_result.read_name = read_result.read_barcode + "-" + read_result.read_umi;
        }
    }
    return read_result;
}


//2.
/* 一条read比对到不同的两处时,
   获取最大比对长度的那个;
   该函数是获取比对的match长度, 只需要对M前面的数字做加和;
*/
int get_read_match_gene_length(std::vector<std::array<int,2>> CIGAR_vector){
    //初始化变量, 累加M前面的数字之和
    int Mlength = 0;
    for (auto each_array : CIGAR_vector) {
        Mlength = Mlength + (each_array[1] - each_array[0]);
    }
    return Mlength;
}


//3. ***** 这里的函数编写结束后得到的小区间都是 1-based; *****
/*  给定CIGAR的string, 获取read match的区间;
    输入是string, 输出是vector, vector中每个元素
    ★Tips1：BAM/BED文件是0-based, SAM文件是1-based;
    ★Tips2：pysam中get_blocks()是0-based;
    ★Tips3：gtf文件是1-based;
*/
struct Read_intervals_and_Mlength{
    std::vector<std::array<int,2>> ReadIntervals;
    int ReadMatchLength;
};

Read_intervals_and_Mlength get_read_intervals(std::string CIGARvalue, std::string initpos){
    //初始化每条reads match的区间; 初始化结构体;
    Read_intervals_and_Mlength CIGAR_Interval_Length;
    //初始化read匹配到基因组上的长度; 用来区分是第一好比对还是第二好比对;
    int MatchLength = 0;
    
    //首先reads的初始位置转化为int型;
    int position_begin = std::stoi(initpos);
    int position_end = 0;
    //初始化三个变量;
    //1.保存每次循环到的字符结果; 2.保存每次循环到的数字结果; 3.指示数字or字母;
    char currentChar = '0';
    int currentNumber = 0;
    bool inNumberMode = false;
    
    //对整个CIGAR string循环;
    for (char ch : CIGARvalue) {
        //判断string的每个字符是否是数字;
        if (std::isdigit(ch)) {
            //如果当前是数字, 更新当前的数字;
            inNumberMode = true;
            currentNumber = currentNumber * 10 + (ch - '0');
        } else {
            //如果当前是字母, 更新当前的字母;
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
                //相比参考基因组删除, 序列位置需要往前增加;
                position_end = position_begin + currentNumber;

            } else if (currentChar == 'I') {
                MatchLength = MatchLength + currentNumber;
                //相比参考基因组插入, 序列位置不需要改变;
                position_end = position_begin;

            } else if (currentChar == 'N') {
                //基本就是SJs;
                position_end = position_begin + currentNumber;

            } else if (currentChar == 'S') {
                MatchLength = MatchLength + currentNumber;
                //刚开始的那些, reads的前多少位软跳跃, 在SEQ中;
                position_end = position_begin;

            } else if (currentChar == 'H') {
                MatchLength = MatchLength + currentNumber;
                //刚开始的那些, reads的前多少位硬跳跃, 不在SEQ中;
                position_end = position_begin;

            } else if (currentChar == 'P') {
                MatchLength = MatchLength + currentNumber;
                //参考基因组和reads都padding, 序列位置都往前增加;
                position_end = position_begin + currentNumber;

            } else {
                //其他没有存放的情况; 一般不会运行到这里;
                std::cout << currentChar << ' ' << "has not been considered yet!" << std::endl;
            }
            //将每次得到的小区间存储到vector容器中;
            position_begin = position_end;
            //重置数字的状态和数字的大小;
            inNumberMode = false;
            currentNumber = 0;
            //std::cout << "停止在结束的位置:" << position_begin << " ";
        }    
    }
    CIGAR_Interval_Length.ReadMatchLength = MatchLength;
    return CIGAR_Interval_Length;
}


//4.***** 获取每条reads的M匹配的小区间, 目前这个函数只能获取reads的cluster; *****
/*  每次是根据文件流分割读取;
*/
//定义结构体, 这是返回的类型;
struct Reads_Clusters {
    // 上一行的位置；
    std::streampos lastPos;
    // 标记读取文件的位置; 这一行的位置;
    std::streampos newPos;
    // cluster中每个read获取的区间信息;
    std::map<std::string, std::vector<std::array<int,2>>> Mymap;
    //每一条reads和对应的barcode;
    std::unordered_map<std::string, std::string> ReadsBarcode;
    // 每一条reads对应的UMI;
    std::unordered_map<std::string, std::string> ReadsUMI;
    // cluster中read的长度;
    std::map<std::string, int> Mylen;
    // cluster中reads比对到的序列名称; 染色体;
    std::string SetRef_name;
    // cluster中范围;
    std::array<int,2> ClusterCoverage;
    //统计reads数量;
    int mapqless1 = 0;
    int sameumi = 0;
    int mappingdiff = 0;
    int surveyNum = 0;    
}; 


Reads_Clusters get_each_cluster_reads_sc(std::ifstream& samfile, std::streampos CurrentPos, std::streampos EndPos) {
    //初始化结构体;
    Reads_Clusters NewCluster = {};
    //创造一个map型数据, key表示每条reads的名称, value表示每条reads的性质;
    //使用一个vector存储read名称, 重复mapping的选择其中一个;
    std::map<std::string, std::vector<std::array<int,2>>> read_informs;
    std::map<std::string, int> read_len;
    // UMI:reads的map, 只要是看有没有重复;
    // Barcode UMI reads;
    std::unordered_map<std::string, std::unordered_map<std::string, std::pair<std::string, int>>> Barcode2UMI2reads;
    //记录每条reads的barcode;
    std::unordered_map<std::string, std::string> read_of_Barcode;
    std::unordered_map<std::string, std::string> read_of_UMI;
    //map中保存的最终的结果向量, 保存有[1], [2], 比对的终止位置, [4];
    //下面这两个只是为了得到最终的基因名称set, 存在多重比对的情况就需要把以前的情况存一下然后剔除差生;
    //每一行视为string;
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
    //这里整个读取sam文件, 然后对文件进行操作
	while (getline(samfile, line))
	{
		/*首先读取出每一行, 然后对字符串进行切割, 需要提取的是：
        第一列：reads的名称 (key)
        第二列：比对到的参考基因组的状况如正负链 [0]
        第三列：比对到的基因的名称 [1]
        第四列：比对起始的位置 [2]
        第五列：比对的质量分数MAPQ [3]
        第六列：CIGAR值  [4] */
        // std::cout << line << std::endl;
        earlyPos = CurrentPos; //文件流指针指向最新, 读完这一行就指向这一行的末尾;
        CurrentPos = samfile.tellg();
        //首先根据判断将SAM文件中的header过滤;
        if (line[0] != '@'){
            NewCluster.surveyNum = NewCluster.surveyNum + 1;
            Split_Result This_Line = get_string_split_sc(line, '\t'); 
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
                    //在这个map中没有这个UMI或者是已经有了, 但比新的短就需要更新, 否则不需要;
                    Barcode2UMI2reads[This_Line.read_barcode][This_Line.read_umi].first = This_Line.read_name;
                    Barcode2UMI2reads[This_Line.read_barcode][This_Line.read_umi].second = This_Line.read_length;                      
                    //获取每条reads match的区间;
                    Read_intervals_and_Mlength CIGAR_interval = get_read_intervals(This_Line.tokens[4], This_Line.tokens[2]);
                    if (This_Line.read_length == CIGAR_interval.ReadMatchLength) {
                        //根据基因和覆盖范围判断每次从文件中截取的reads数量;
                        //该read所在的参考序列名称;
                        now_gene = This_Line.tokens[1];
                        //read匹配的开始位置;
                        int read_in_gene_begin_pos = std::stoi(This_Line.tokens[2]);
                        std::array<int,2> lastEle = CIGAR_interval.ReadIntervals.back();
                        //前闭后开;
                        //read匹配的结束位置;
                        int read_in_gene_end_pos = lastEle[1] - 1;
                        //如果是第一次, 就记录这个read的范围当作判断标准;
                        if (early_begin_pos == 0) {
                            early_begin_pos = read_in_gene_begin_pos;
                            early_end_pos = read_in_gene_end_pos;
                            last_chr = now_gene;
                        } else {
                            //不是开头;
                            if (now_gene != last_chr){
                                //下面这一部分就是存进去了;
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
                                //这种情况就是没有重叠;
                                if ((read_in_gene_begin_pos - early_end_pos > 0) || (early_begin_pos - read_in_gene_end_pos > 0)){
                                    //记录当前文件流的位置;
                                    //下面这一部分就是存进去了;
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
                                    NewCluster.surveyNum = NewCluster.surveyNum - 1;
                                    break;
                                } else {
                                    //有重叠的情况就把覆盖范围改正;
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
    //到这里一切循环都运行结束了;
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


//5.***** 根据计算出的小区间获取每个Cluster中的每条reads的Splice Junctions *****
/*  暂且还是跟python一样的思路, 循环每条reads及其匹配的小区间
    函数的输出是1.SJ的map; 2.SJs_left的map; 3.SJs_right的map;
*/
//定义结构体, 作为函数返回的类型;
//定义一个SpliceJs结构体, 这是返回的类型;
struct SpliceJs {
    //map存储每条reads的Splice junctions;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> reads_SJs;
    //map存储每条reads左边的M区间, 卡的条件的10bp;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> reads_SJs_left;
    //map存储每条reads右边的M区间, 卡的条件的10bp
    std::unordered_map<std::string, std::vector<std::array<int,2>>> reads_SJs_right;
    //存储read的开头的结尾;
    std::unordered_map<std::string, std::array<int,2>> reads_begin_end;
    //map存储单个exon的reads的开头和结尾;
    std::unordered_map<std::string, std::array<int,2>> reads_single_exon;

}; 

//函数部分, 已知匹配小区间获取SJs;
SpliceJs get_reads_allSJs(std::map<std::string, std::vector<std::array<int,2>>>& RCs, int SJD) {
    //初始化结构体;
    SpliceJs NewSJs = {};
    //按照个数循环;
    // int count = 0;
    // for (auto it = RCs.begin(); it != RCs.end() && count < 5; ++it, ++count) {
    //     std::cout << "Key: " << it->first << ", Value: ";
    //     std::cout << std::endl;
    // }
    //定义结构体上的三个SJ map和一个单exon map;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> AllRead_SJs;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> AllRead_SJs_left;
    std::unordered_map<std::string, std::vector<std::array<int,2>>> AllRead_SJs_right;
    std::unordered_map<std::string, std::array<int,2>> AllRead_begin_end;
    std::unordered_map<std::string, std::array<int,2>> AllRead_SingleExon;
    //对整个RCs循环, 一次是一条reads;
    for (const auto &pair : RCs) {
        //提取出来;
        std::string Rname = pair.first; //read名称;
        std::vector<std::array<int,2>> RMs = pair.second;//read M小区间;

        //初始化SJ的vector, 左侧的vector, 右侧的vector;
        std::vector<std::array<int,2>> ReadSJs;
        std::vector<std::array<int,2>> ReadSJs_left;
        std::vector<std::array<int,2>> ReadSJs_right;
        //初始化单个exon开始和结尾的array;
        std::array<int,2> Read_SingleExon;
        std::array<int,2> Read_begin_end;

        //对所有M小区间循环, 一次是前后一对M小区间;
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
                //满足条件的这个SJ, 并添加到SJ所在的vector中;
                std::array<int,2> eachSJ;
                eachSJ[0] = former;
                eachSJ[1] = later - 1;
                ReadSJs.push_back(eachSJ);
                //满足条件SJ的左侧, 并添加到SJ_left的vector中;
                std::array<int,2> eachSJ_left;
                eachSJ_left[0] = (*(it))[0];
                eachSJ_left[1] = former - 1;
                ReadSJs_left.push_back(eachSJ_left);
                //满足条件SJ的右侧, 并添加到SJ_right的vector中;
                std::array<int,2> eachSJ_right;
                eachSJ_right[0] = later;
                eachSJ_right[1] = (*(it+1))[1];
                ReadSJs_right.push_back(eachSJ_right);
            }
        }
        //单独讨论单个exon的情节;
        if (ReadSJs.size() == 0) {
            Read_SingleExon[0] = (*RMs.begin())[0];
            Read_SingleExon[1] = (*(RMs.end()-1))[1] - 1;
            AllRead_SingleExon[Rname] = Read_SingleExon;
        } else {
            //将每条reads的结果存入map中;
            AllRead_SJs[Rname] = ReadSJs;
            AllRead_SJs_left[Rname] = ReadSJs_left;
            AllRead_SJs_right[Rname] = ReadSJs_right;
            AllRead_begin_end[Rname] = Read_begin_end;
        }
        //判断结束;
    }
    //装入结构体并返回;
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
//7.***** 读取注释文件, gtf中找到已知的isoform; *****
/*  输入是注释文件地址, 输出是map<char:map<gene:<>>>;
*/
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
        
        // 读取GTF文件;
        // 创建文件输出流对象;
        std::ifstream GTFFile; 
        // 打开GTF文件;
        GTFFile.open(GTFFile_name);
        // 判断文件是否打开;
        if (!GTFFile.is_open())	{
            std::cout << "We can't open GTF file !" << GTFFile_name << std::endl;
            exit(EXIT_FAILURE);
        }        
        //逐行读取并输出, 每一行读取出的字符串;
        std::string line;
        std::string now_chr_name = "Han";
        std::string early_chr_name = "Han";
        std::string early_gene_transcript_name = "WHH";
        std::string now_gene_transcript_name = "WHH";
        //初始化每个exon所在的链;
        std::string early_Exonstrand = "Wei";
        std::string now_Exonstrand = "Wei";

        // 下面就是读取每一行;
        while (std::getline(GTFFile, line)){
            //首先排除前面的没用的行;            
            if (line[0] != '#'){
                //对每一行定义一个字符串流;
                std::istringstream lineiss(line);
                //每一行分割的内容保存在字符串向量中;
                std::vector<std::string> tokens;
                //分割出来的每一个单元格内容;
                std::string token;

                //分割并循环gtf文件的每一行;
                while (std::getline(lineiss, token, '\t')) {
                    //不端添加分割的内容;
                    tokens.push_back(token);
                } //以上是对每一行进行切割;
                
                if (tokens[2] == "exon"){
                    early_chr_name = now_chr_name;
                    now_chr_name = tokens[0]; //染色体名称;
                    early_Exonstrand = now_Exonstrand;
                    now_Exonstrand = tokens[6]; //正负链;

                    //然后对最末行的tokens切割;
                    std::string AllEndT = tokens.back();
                    //将其定义为一个字符串流;
                    std::istringstream EndTT(AllEndT);
                    std::string Endtoken;
                    std::vector<std::string> Endtokens;
                    //分割每一行的最后一列的字符串;
                    while (std::getline(EndTT, Endtoken, '"')){
                        Endtokens.push_back(Endtoken);
                    }
                    early_gene_transcript_name = now_gene_transcript_name;
                    now_gene_transcript_name = Endtokens[1]+'|'+Endtokens[3];

                    if ((GTFAnno_SJs.size() != 0) && (early_gene_transcript_name != now_gene_transcript_name)){
                        //如果在负链, 就需要反转一下使用;
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
        //这里添加的是判断文件是否读取到最后;
                //到这里的情况下, 最后一个chr的内容还没有存储到输出的变量GTFAll_Info中;
        if (GTFFile.eof()) {
            //文件正常读取到最后;
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

        // 开始遍历每一个染色体;        
        for (const auto& eachChr:GTFAll_Info.GTF_transcript){
            // 现在是每个染色体的开始;
            now_chr = eachChr.first;
            // std::cout << now_chr << std::endl;
            GTFAll_Info.GTF_gene[now_chr] = {};
            GTFAll_Info.GTF_gene_strand[now_chr] = {};
            GTFAll_Info.GTF_gene2transcript[now_chr] = {};

            // 染色体的最后一个基因的情况;
            if (All_gene_begin.size() != 0 && now_chr != ago_chr) {
                auto min_it = std::min_element(All_gene_begin.begin(), All_gene_begin.end());
                auto max_it = std::max_element(All_gene_end.begin(), All_gene_end.end());
                GTFAll_Info.GTF_gene[ago_chr][ago_gene_name][0] = *min_it;
                GTFAll_Info.GTF_gene[ago_chr][ago_gene_name][1] = *max_it;
                GTFAll_Info.GTF_gene2transcript[ago_chr][ago_gene_name] = tx_name;
                // GTFAll_Info.GTF_gene_strand[ago_chr][ago_gene_name] = Strand;
                All_gene_begin.clear();
                All_gene_end.clear();
                tx_name.clear();
                ago_gene_name = "";
            }

            // 遍历每个转录本;
            for (const auto& eachGene:eachChr.second) {
                size_t pos = eachGene.first.find('|');
                now_gene_name = eachGene.first.substr(0, pos);
                Strand = GTFAll_Info.GTF_transcript_strand[now_chr][eachGene.first];
                // std::cout << Strand << std::endl;
                // 同一个基因的情况;
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
        //如果是空的;
        std::cout << "***** No GTF files are put into the program! *****" << std::endl;         
    }

    // for (const auto& aaa:GTFAll_Info.GTF_gene2transcript) {
    //     std::cout << "染色体" << aaa.first << ":" << std::endl;
    //     std::cout << aaa.second.size() << std::endl;
    //     if (aaa.first == "chrM") {
    //         for (const auto& bbb:aaa.second) {
    //             std::cout << "基因" << bbb.first << ":" << GTFAll_Info.GTF_gene_strand[aaa.first][bbb.first] << std::endl;
    //             for (const auto& ccc:bbb.second) {
    //                 std::cout << ccc << " ";
    //             }
    //             std::cout << std::endl;
    //         }
    //     }
    // }
    // for (const auto& aaa:GTFAll_Info.GTF_gene) {
    //     std::cout << "染色体" << aaa.first << ":" << std::endl;
    //     std::cout << aaa.second.size() << std::endl;
    //     if (aaa.first == "chrM") {
    //         for (const auto& bbb:aaa.second) {
    //             std::cout << "基因" << bbb.first << ":" << GTFAll_Info.GTF_gene_strand[aaa.first][bbb.first] << std::endl;
    //             std::cout << "[" << bbb.second[0] << "," << bbb.second[1] << "]" << std::endl;
    //         }
    //     }
    // }    
    // 输出结果;
    return GTFAll_Info;
}


//9. ***** 根据注释的map获取SJ信息和单exon信息 *****
/*  输入是上一步获得的exon的map；
    输出是结构体, 包含两部分, 一部分是非单个exon的转录本及其SJs, 另一部分是单个exon的转录本及其exon;
*/
struct GTFsj{
    //所有转录本的剪切位点map;
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::array<int,2>>>> mSJs;
    //所有转录本的开头和结尾;
    std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> mSJsBE;
    //单个exon的map;
    std::unordered_map<std::string, std::unordered_map<std::string, std::array<int,2>>> SE;
};

//函数部分;
GTFsj get_SJs_SE(std::map<std::string, std::map<std::string, std::vector<std::array<int,2>>>>& Known_Exon){
    //初始化结构体;
    GTFsj Known_SJ_SE;
    //如果有注释, 传入的变量大小不会为0个, 反之, 就有0个;
    if (Known_Exon.size() != 0){
        //循环用到的变量;
        std::map<std::string, std::vector<std::array<int,2>>> Every_transcript_All;
        std::vector<std::array<int,2>> Every_transcript;
        std::string transaction_name;
        std::string chr_name;
        
        //存储新的值的变量;
        //前四个是存储的SJ相关的;
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
            //pair.first是chr1, SIRV1
            //pair.second是map<gene|isoform <><>...<><>>;
            chr_name = pair.first;
            Every_transcript_All = pair.second;
            //每个chr清理一下之前的;
            Every_SJs_All.clear();
            Every_SE_All.clear();
            Every_SJ_begin_end_all.clear();

            for (auto it = Every_transcript_All.begin(); it != Every_transcript_All.end(); ++it) {
                //it->first是gene|isoform名称, 如SIRV101;
                transaction_name = it->first;
                Every_transcript = it->second;

                if (Every_transcript.size() == 1){
                    //这种就是单个exon的情况;
                    Every_SE[0] = (Every_transcript.front())[0];
                    Every_SE[1] = (Every_transcript.front())[1];
                    Every_SE_All[transaction_name] = Every_SE;
                } else {
                    //这种就是正常的情况;
                    //对每个transcript循环;
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



//11.读取FASTA文件;
std::unordered_map<std::string, std::string> Read_fasta_file(std::string& fastafile_name){
// void Read_fasta_file(const char* fastafile_name){
    //最终想要保存的结果;
    std::unordered_map<std::string, std::string> chrgeneString;
    //定义一个vector用来存放所有的string;
    std::vector<std::string> VecStrings;

    //正确读取;
    std::cout << "***** Now open Fasta file: " << fastafile_name << "! *****" << std::endl;
    //最终想要的输出大结果<chr:<>, chr:<>, chr:<>>

    //读取GTF文件;
    //创建文件输出流对象;
    std::ifstream FastaFile; 
    // 打开GTF文件；
	FastaFile.open(fastafile_name);
    // 判断文件是否打开;
	if (!FastaFile.is_open()) {
		std::cout << "We can't open Fasta file !" << fastafile_name << std::endl;
		exit(EXIT_FAILURE);
	}
    //逐行读取并输出, 每一行读取出的字符串;
    std::string ChrGeneName;
    std::string concatenatedString;
    std::stringstream accumulatedStringStream;
    std::string line;
    std::string token;
    
    while (std::getline(FastaFile, line))
    {      
        // count++;
        // std::cout << count << std::endl;
        //!=就是表示字符串中包含>字符;如果没有找到,就会返回npos这个字符;
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
                // std::cout << token << std::endl;
                if (count == 1){
                    ChrGeneName = token.substr(1);
                }
            }
            //将读取了的染色体名称输出;
            // std::cout << ChrGeneName << std::endl;
        } else {
            //如果不包含>字符, 那就是字符串部分;
            accumulatedStringStream << line;
        }
    }

    if (!ChrGeneName.empty()) {
        // 处理最后一行记录
        concatenatedString = accumulatedStringStream.str();
        chrgeneString[ChrGeneName] = concatenatedString;
    }

    //这里添加的是判断文件是否读取到最后;
    if (FastaFile.eof()) {
        //文件正常读取到最后;
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



// 判断这个sj是否满足GT-AG信号和周围10bp没有插入删除;
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
    // 将字符串中的每个字符转为大写
    std::transform(signal1.begin(), signal1.end(), signal1.begin(), ::toupper);
    std::transform(signal2.begin(), signal2.end(), signal2.begin(), ::toupper);
    //判断索引的信号
    if (signal1 == "GT" && signal2 == "AG"){
        flag = 1;
        // sign_strand.push_back("+");
    } else if (signal1 == "CT" && signal2 == "AC"){
        flag = 2;
        // sign_strand.push_back("-");

    } else if (signal1 == "GC" && signal2 == "AG"){
        flag = 1;
        // sign_strand.push_back("+");
    } else if (signal1 == "CT" && signal2 == "GC"){
        flag = 2;
        // sign_strand.push_back("-");

    } else if (signal1 == "AT" && signal2 == "AC"){
        flag = 1;
        // sign_strand.push_back("+");
    } else if (signal1 == "GT" && signal2 == "AT"){
        flag = 2;
        // sign_strand.push_back("-");
    } else {
        flag = 0;
    }                
    return flag;
}


// 定义长参数选项
static struct option long_options[] = {
    {"sam", required_argument, 0, 's'},
    {"fasta", required_argument, 0, 'f'},
    {"gtf", required_argument, 0, 'g'},
    {"output", required_argument, 0, 'o'},
    {"SJDistance", required_argument, 0, 'j'},  // splice junction的距离; -d/--SJDistance
    {"support", required_argument, 0, 'n'},   // 支持SJ的read个数 -n/--SJ_support_read_number
    {"single_exon_boundary", required_argument, 0, 'e'},
    {"graph_distance", required_argument, 0, 'd'},  //  -h/--Graph_distance
    {"thread", required_argument, 0, 't'}, // -t/--Thread
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}  // 表示结束
};

bool directoryExists(const std::string& path) {
    struct stat info;
    // 使用 stat() 检测路径
    if (stat(path.c_str(), &info) != 0) {
        // 路径不存在
        return false;
    } else if (info.st_mode & S_IFDIR) {
        // 路径存在并且是目录
        return true;
    } else {
        // 路径存在，但不是目录
        return false;
    }
}

bool createDirectory(const std::string& path) {
    // 使用 mkdir() 创建目录
    if (mkdir(path.c_str(), 0755) == 0) {
        // 创建成功
        return true;
    } else {
        // 如果创建失败，输出错误信息
        std::cerr << "Error creating directory: " << path << strerror(errno) << std::endl;
        return false;
    }
}

std::string joinPath(const std::string& directory, const std::string& filename) {
    // 判断路径末尾是否有 '/'
    if (directory.back() == '/') {
        return directory + filename;  // 如果已经有分隔符，直接拼接
    } else {
        return directory + "/" + filename;  // 如果没有分隔符，添加 '/'
    }
}

// 添加一个函数用于打印帮助信息
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
    std::cout << "  -h, --help                  show this help information." << std::endl;
}

std::vector<std::string> check_catalog_exist(const std::string& output_path) {
    std::vector<std::string> outputFileVector;
    // 检测目录是否存在;
    if (!directoryExists(output_path)) {
        std::cout << "Directory does not exist, creating it..." << std::endl;
        // 如果目录不存在, 创建它;
        if (createDirectory(output_path)) {
            std::cout << "Directory created successfully: " << output_path << std::endl;
        } else {
            std::cerr << "Failed to create directory: " << output_path << std::endl;
            exit(EXIT_FAILURE);
        }
    } else {
        std::cout << "Directory already exists: " << output_path << std::endl;
    } 

    // 目录已经建立;
    // 文件一, 更新的gtf文件;
    std::string updatedGtfPath = joinPath(output_path, "updated_annotations.gtf");
    outputFileVector.push_back(updatedGtfPath);
    std::ofstream gtf_file(updatedGtfPath, std::ios::trunc);
    gtf_file.close();

    // 文件二, 每条reads的追踪文件;
    std::string TracePath = joinPath(output_path, "compatible_isoform.tsv");
    outputFileVector.push_back(TracePath);
    std::ofstream TraceIsoform(TracePath, std::ios::trunc);
    if (TraceIsoform.is_open()){
        TraceIsoform << "read_id" << '\t' << "category" << '\t' << "isoform_id" << '\t' << "gene_id" << '\t' << "barcode_id" << '\t' << "file" << '\n';
    }
    TraceIsoform.close();   
    return outputFileVector;
}


// 给定sam的地址进行解析, 多个文件还是单个文件;
std::vector<std::string> traverse_sam_file(const std::string& sam_file_path, const std::string& output_path){
    // 初始化文件向量;
    std::vector<std::string> sam_file_vector;
    // 定义stat文件结构;
    struct stat sam_stat;
    // 调用stat函数获取path的状态信息, 将结果存储到sam_stat中;    
    stat(sam_file_path.c_str(), &sam_stat);
    // 判断是不是一个文件;    
    if (S_ISREG(sam_stat.st_mode)) {
        // True 表示这是一个文件;
        std::cout << "* Only one sam file is entered! * << " << sam_file_path << std::endl;
        // 一个文件也放入文件vector中;
        sam_file_vector.push_back(sam_file_path);

    } else if (S_ISDIR(sam_stat.st_mode)) {
        // True 表示是一个文件夹;
        std::cout << "* A folder was entered! * << " << sam_file_path << std::endl;
        DIR* dir = opendir(sam_file_path.c_str());
        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr) {
            if (entry->d_name[0] != '.') { // 忽略隐藏文件
                std::string fileName = entry->d_name;
                std::cout << "Find file: " << fileName << std::endl;
                // 找到后缀是sam的文件;
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

// 找到给定断点的下一行;
std::streampos findNextLineStart(std::ifstream& file, std::streampos pos) {
    file.seekg(pos);
    std::string line;
    // 找到下一个行首
    if (std::getline(file, line)) {
        return file.tellg(); // 返回下一个行的开始位置;
    }
    return pos; // 如果没有下一行, 返回原始位置;
}


// 单个sam文件的情况, 也是多线程中的单个线程的情况;
// 输入是这个sam文件的地址;
void processChunk(const std::string& one_sam_file_path, 
                const std::streampos& start, 
                const std::streampos& end, 
                const std::string& output_path, 
                const int file_i, 
                const int SJ_Distance, 
                std::unordered_map<std::string, std::string>& FastaRef){
    // 初始化变量;
    int Group_index = 0;
    std::string chrchr;
    Reads_Clusters each_cluster_informs;
    SpliceJs each_read_SJs_informs;
    std::array<int,2> readBeginEnd;
    // std::map<std::array<int, 2>, std::array<int,2>> SpliceJunction; // key表示SJ, value表示数量,正负链;
    bool thisSjNoError = 0;
    int thisSjSignal = 0;
    int SjNumber = 0; //对Sj数数;
    // 生成小部分的临时文件1:所有reads的SJ和信息;
    // std::cout << "Thread begin" << std::endl;
    // read_id; chr; group_id; group_begin; group_end; read_begin; read_end; read_length; read_SJs
    std::ostringstream oss;
    oss << "Read_" << std::setw(4) << std::setfill('0') << file_i << ".txt";
    std::string File_name = oss.str();
    std::string ReadInformPath = joinPath(output_path, File_name);
    std::ofstream ReadInform(ReadInformPath, std::ios::trunc);

    // 读取SAM文件;
    // 创建文件输出流对象; 打开samfile文件;
    std::ifstream samfile(one_sam_file_path);
    samfile.seekg(start);
    // 记录当前, 也就是文件刚开始的位置;
    std::streampos Current_Position = samfile.tellg();
    std::streampos Last_Position = Current_Position;

    while (Current_Position < end) {
        Group_index++;
        each_cluster_informs = get_each_cluster_reads_sc(samfile, Last_Position, end);
        chrchr = each_cluster_informs.SetRef_name;
        Last_Position = each_cluster_informs.lastPos;
        Current_Position = each_cluster_informs.newPos;
        // 提取每条reads的splice junction;
        each_read_SJs_informs = get_reads_allSJs(each_cluster_informs.Mymap, SJ_Distance);  /* 这里有单个exon的信息; */
                                                                                            /* 单个exon的可以从这里思考写入; */
        // 每条reads的信息写入临时文件1;
        // 每个cluster并不是所有的都有用输出出来了, 比如说单个exon的; /* 下面可以考虑添加进来 */ ★
        // 这里写的是有SJ的reads;
        // std::cout << "猴哥: " << each_read_SJs_informs.reads_SJs.size() << std::endl;
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
                // 周围10bp没有错误;
                thisSjNoError = ifSjNoError(each_read_SJs_informs.reads_SJs_left[eachRead.first], each_read_SJs_informs.reads_SJs_right[eachRead.first], SjNumber); 
                if (thisSjNoError == 1) {
                    // 是不是GT-AT信号; 正链为+, 1; 负链为-, 2;
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
        
        // 这里可以写没有SJ的reads; single exon;
        for (const auto& eachRead:each_read_SJs_informs.reads_single_exon) {
            readBeginEnd = eachRead.second;
            ReadInform << eachRead.first << '\t' << chrchr << '\t' << Group_index << '\t'
                        << each_cluster_informs.ClusterCoverage[0] << '\t' 
                        << each_cluster_informs.ClusterCoverage[1] << '\t' << readBeginEnd[0] << '\t' 
                        << readBeginEnd[1] << '\t' << each_cluster_informs.Mylen[eachRead.first] << '\t'
                        << each_cluster_informs.ReadsBarcode[eachRead.first] << '\n';         
        } // 写入文件结束;
    } // 循环结束;
    
    samfile.close();
    ReadInform.close();
    std::cout << "^-^ Thread: " << std::this_thread::get_id() << " has completed processing! ^-^" << std::endl;
}


// 切割每个 <Read_x.txt> 文件;
Split_Result get_line_split(const std::string& s, char delimiter){
    // 初始化返回值结构体;
    Split_Result read_result = {};
    // 字符串流;
    std::istringstream iss(s);
    // 定义字符串流中的每一段字符串;
    std::string token;
    // 定义一个number, 提取特定有用的信息;
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

// 切割每个 <SJ_x.txt> 文件;
Split_Result get_sj_split(const std::string& s, char delimiter){
    // 初始化返回值结构体;
    Split_Result sj_result = {};
    //字符串流;
    std::istringstream iss(s);
    //定义字符串流中的每一段字符串;
    std::string token;
    // 把所有信息都加进来;
    while (std::getline(iss, token, delimiter)){
        sj_result.tokens.push_back(token);
    }
    return sj_result;
}


// 每个sample文件都对应一个这个数据结构, 最终版的
struct FileSplit {
    // 多个文件用来合并同类项;
    std::map<std::string, std::vector<std::array<int,2>>> chr_coverage; //一个sam文件才有;
    std::map<std::array<int,2>, std::array<std::streampos,2>> coverage2pos; //一个sam文件才有;
    // 后面用来读文件, 这是如果只有这一个sample;
    std::vector<std::streampos> reads_pointer; //希望的最终形态;
    std::vector<int> group_reads_number; //希望的最终形态;
    // 输出的文件;
    std::string readtxt_path; //希望的最终形态;
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


// 合并小文件成原始的文件;
FileSplit Merge_Read_Small_Files(const std::string& SmallFilePath, const int& SAMFileNumber, const int& nthread) {
    FileSplit Chunk_Bang;
    // 1.读取文件夹下所有文件的名称;
    // 所有的Read_x.txt文件;
    std::vector<std::string> Read_x_vec;

    DIR* dir = opendir(SmallFilePath.c_str());
    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_name[0] == 'R') { // 忽略隐藏文件
            Read_x_vec.push_back(entry->d_name);
        } 
    }
    std::sort(Read_x_vec.begin(), Read_x_vec.end()); //排序的问题解决了; ★
    closedir(dir);

    // 现在已经把所有的这个文件读取出来了; 只有文件名称, 下面需要一个个合并;
    // 2.生成新文件;
    std::streampos Readpos;
    std::string File_total_name = SmallFilePath + "/All_Read.txt"; 
    std::ofstream AllReadInform(File_total_name, std::ios::trunc);  
    Chunk_Bang.readtxt_path = File_total_name;
    Readpos = AllReadInform.tellp();
    Chunk_Bang.reads_pointer.push_back(Readpos);

    // 3. 循环读取文件, 合并所有的Read_x.txt文件;
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
        // 文件绝对路径;
        std::string fileName = SmallFilePath + "/" + Read_x_vec[file_number];
        std::ifstream SmallSamFile(fileName);
        earlyGroup = "+";
        thisGroup = "+";

        while (getline(SmallSamFile, line)) {
            Split_Result This_Line = get_line_split(line, '\t');
            cellSet.insert(This_Line.tokens[7]);
            thisGroup = This_Line.tokens[1];
            this_coverage = {std::stoi(This_Line.tokens[2]), std::stoi(This_Line.tokens[3])};
            // 这种的就是默认文件1的最后一个cluster与文件二的第一个cluster有交集; (也还好, 稍微大一点)
            // 最戏剧的是, 第一个文件正好分成了chr1, 第二个分成了chr2......; 这种情况是两个文件都是Group1, 但第一个文件标注chr1, 第二个文件标注chr2;
            // 都考虑到了;
            if (thisGroup == earlyGroup) {
                // 这种是上一行跟这一行相同的情况;
                // 看看这条read重复不重复;
                auto it = Group_temporary.find(This_Line.read_name);
                if (it != Group_temporary.end()) {
                    if (std::stoi(This_Line.tokens[6]) > std::stoi(Group_temporary[This_Line.read_name][6])){
                        Group_temporary[This_Line.read_name] = This_Line.tokens;
                    }
                } else {
                    Group_temporary[This_Line.read_name] = This_Line.tokens;
                }
            } else { 
                // (thisGroup != earlyGroup)
                // 这种是上一行跟这一行不同的情况; 不同有几种? 初始化不同; 染色体不同; 覆盖范围不同;
                if (Group_temporary.size() != 0) {
                    //临时集合有元素;
                    if (last_chr == This_Line.tokens[0]) {
                        //染色体相同;
                        if (Group_low <= this_coverage[1] && this_coverage[0] <= Group_high){
                            // 两个有交集;
                            // 看看这条read重复不重复;
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
                            // 写入;
                            // 这种情况是这个Group结束了; 写入新文件;
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
                            // 单个文件的时候标记位置;
                            Readpos = AllReadInform.tellp();
                            Chunk_Bang.coverage2pos.insert({{Group_low, Group_high}, {Chunk_Bang.reads_pointer.back(),Readpos}});
                            Chunk_Bang.reads_pointer.push_back(Readpos); //这个Group结束;
                            Chunk_Bang.reads_pointer.push_back(Readpos); //下个Group开始;
                            Chunk_Bang.group_reads_number.push_back(Group_temporary.size());

                            // 收集每个chr里面的coverage;
                            auto it = Chunk_Bang.chr_coverage.find(chr_name);
                            if (it != Chunk_Bang.chr_coverage.end()){
                                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                            } else {
                                Chunk_Bang.chr_coverage[chr_name] = {};
                                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                            }
                            
                            Group_temporary.clear();
                            // Barcode2UMI2reads.clear();
                            // 新的类需要把这个写进去;
                            Group_temporary[This_Line.read_name] = This_Line.tokens;
                            Group_high = this_coverage[1];
                            Group_low = this_coverage[0];
                        }
                    } else {
                        //染色体不同;
                        //写入;
                        // 这种情况是这个Group结束了; 写入新文件;
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
                        // 单个文件的时候标记位置;
                        Readpos = AllReadInform.tellp();
                        Chunk_Bang.coverage2pos.insert({{Group_low, Group_high}, {Chunk_Bang.reads_pointer.back(),Readpos}});
                        Chunk_Bang.reads_pointer.push_back(Readpos); //这个Group结束;
                        Chunk_Bang.reads_pointer.push_back(Readpos); //下个Group开始;     
                        Chunk_Bang.group_reads_number.push_back(Group_temporary.size());                     
                        
                        // 收集每个chr里面的coverage;
                        auto it = Chunk_Bang.chr_coverage.find(chr_name);
                        if (it != Chunk_Bang.chr_coverage.end()){
                            Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                        } else {
                            Chunk_Bang.chr_coverage[chr_name] = {};
                            Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
                        }
                   
                        Group_temporary.clear();
                        // Barcode2UMI2reads.clear();
                        //新的类需要把这个写进去;
                        Group_temporary[This_Line.read_name] = This_Line.tokens;    
                        Group_high = this_coverage[1];
                        Group_low = this_coverage[0];
                        // Barcode2UMI2reads[This_Line.read_barcode] = {};
                        // Barcode2UMI2reads[This_Line.read_barcode][This_Line.read_umi] = {This_Line.read_name, This_Line.read_length};                                                                      
                    }
                } else {
                    //第一个文件初始化, 临时集合没有元素; 直接存元素;
                    Group_temporary[This_Line.read_name] = This_Line.tokens;
                    Group_high = this_coverage[1];
                    Group_low = this_coverage[0];
                    // Barcode2UMI2reads[This_Line.read_barcode] = {};
                    // Barcode2UMI2reads[This_Line.read_barcode][This_Line.read_umi] = {This_Line.read_name, This_Line.read_length}; 
                }
            }
            earlyGroup = thisGroup; 
            last_chr = This_Line.tokens[0];
        }

        if (file_number == (Read_x_vec.size()-1) && SmallSamFile.eof()){
            // 这种情况是这个Group结束了; 写入新文件;
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
            // 最后也留下位置;
            Readpos = AllReadInform.tellp();
            Chunk_Bang.coverage2pos.insert({{Group_low, Group_high}, {Chunk_Bang.reads_pointer.back(),Readpos}});
            Chunk_Bang.reads_pointer.push_back(Readpos);
            Chunk_Bang.group_reads_number.push_back(Group_temporary.size());
            
            // 收集每个chr里面的coverage;
            auto it = Chunk_Bang.chr_coverage.find(chr_name);
            if (it != Chunk_Bang.chr_coverage.end()) {
                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
            } else {
                Chunk_Bang.chr_coverage[chr_name] = {};
                Chunk_Bang.chr_coverage[chr_name].push_back({Group_low, Group_high});
            }
        }
        SmallSamFile.close();
        // 下面这行代码将所有的分线程的小文件删除 <Read_x.txt>;
        std::remove(fileName.c_str());
    }
    AllReadInform.close();
    Chunk_Bang.FileNo = 1;
    Chunk_Bang.file_barcodeSet[std::to_string(SAMFileNumber)] = cellSet;
    return Chunk_Bang;
}


std::map<std::string, std::vector<std::array<int,2>>> Merge_Read_Interval(std::map<int, std::map<std::string, std::vector<std::array<int,2>>>>& FileCoverage) {
    // 所有文件最终的输出;
    std::map<std::string, std::vector<std::array<int,2>>> chr_range_temp;
    std::vector<std::array<int,2>> raw_intervals;
    std::vector<std::array<int,2>> sorted_intervals;
    if (FileCoverage.size() > 0) {
        // 文件数量大于1; 遍历所有文件;
        for (const auto& eachFile:FileCoverage) {
            if (eachFile.second.size() > 0) {
                // 单个文件的遍历每个染色体;
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
    // 所有文件的所有chr整理完毕;
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



// struct FileSplit {
//     // 多个文件用来合并同类项;
//     std::map<std::string, std::vector<std::array<int,2>>> chr_coverage; //一个sam文件才有;
//     std::map<std::array<int,2>, std::array<std::streampos,2>> coverage2pos;
//     // 后面用来读文件, 这是如果只有这一个sample;
//     std::vector<std::streampos> reads_pointer; //希望的最终形态;
//     // std::unordered_map<std::string, std::streampos> chr_pointer_dict;   //方便合并大文件使用;
//     std::vector<int> group_reads_number; //希望的最终形态;
//     // 输出的文件;
//     std::string readtxt_path; //希望的最终形态;
//     int FileNo;
// };
// 合并每个sample的大文件生成总文件;
// read_id; chr; group_id; group_begin; group_end; read_begin; read_end; read_length; read_SJs
void Merge_Read_Big_Files_Group(const int& NEWNEW,
                                const std::string& output_path, 
                                const std::array<int,2>& Regin_BE,
                                std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>>& file_pointers,
                                const std::string& chrchr,
                                std::vector<std::streampos>& group_read_pointer,
                                std::vector<int>& group_read_number,
                                std::ofstream& BigFile) {
    // std::cout << "[" << Regin_BE[0] << "," << Regin_BE[1] << "] " << chrchr << std::endl;
    // 用到的变量;
    std::streampos ThisFile_ThisGroup_pointer;
    std::streampos end;
    std::string line;
    std::map<int, std::vector<std::vector<std::string>>> Group_store;
    // 载入这个任务的变量;
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

    // 写入总文件中;
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


// 输入1:总的每个染色体的区域的范围;
std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>> get_pointers(
        std::map<std::string, std::vector<std::array<int,2>>>& chrcoverage,
        std::map<int, std::map<std::string, std::vector<std::array<int,2>>>>& FileCoverage,
        std::map<int, std::map<std::array<int,2>, std::array<std::streampos,2>>>& groups_pointer) {
    // 这是最终保存的一个染色体的结构;
    std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>> chr_pointer_coverage;
    // 先最终结果初始化;
    for (const auto& EachChr:chrcoverage) {
        chr_pointer_coverage[EachChr.first] = {};
        for (const auto& Region:EachChr.second) {
            chr_pointer_coverage[EachChr.first][Region] = {};
        }
    }
    // 初始化完毕;
    std::array<std::streampos,2> thisFile_group_be;
    // 开始每个文件的区域; 第一个sam, 第二个sam, ......
    for (const auto& EachFile:FileCoverage) {
        // int GroupAll = -1;
        int FileNumber = EachFile.first;
        std::map<std::array<int,2>, std::array<std::streampos,2>> thisFile_group_pointer = groups_pointer[FileNumber];
        // 遍历每个染色体的区域;
        for (const auto& EachChr:EachFile.second) {
            int Rhino = 0; //每个大区域索引;
            int Rhino_small = 0;
            std::vector<std::streampos> temp_pointer; // 大区域对应的小区域的所有指针, 开始和终止;
            std::string ChrName = EachChr.first;
            std::vector<std::array<int,2>> Big_Rhino = chrcoverage[ChrName]; //所有的chr是都有的;
            // std::cout << ChrName << std::endl;
            // 遍历文件的小区域;
            for (const auto& small:EachChr.second) {
                // 小区间;
                // GroupAll = GroupAll + 1;
                Rhino_small = Rhino_small + 1;
                // 看看属于Big_Rhino第几个;
                for (int i = Rhino; i < Big_Rhino.size(); i++) {
                    if (small[0] >= Big_Rhino[i][0] && small[1] <= Big_Rhino[i][1]) {
                        // 在这个大的范围内;
                        // std::cout << "卡在这里aaa" << std::endl;
                        // std::cout << "哇哇 " << small[0] << " " << small[1] << "   " << Big_Rhino[i][0] << " " << Big_Rhino[i][1] << std::endl;
                        temp_pointer.push_back(thisFile_group_pointer[small][0]);
                        temp_pointer.push_back(thisFile_group_pointer[small][1]);
                        // std::cout << "哇 " << Rhino_small << " " << EachChr.second.size() - 1 << std::endl;
                        if (Rhino_small == EachChr.second.size()) {
                            thisFile_group_be[0] = temp_pointer[0];
                            thisFile_group_be[1] = temp_pointer[temp_pointer.size()-1];
                            chr_pointer_coverage[ChrName][Big_Rhino[i]].insert({FileNumber, thisFile_group_be});
                            temp_pointer.clear();
                        }
                        Rhino = i;
                        break;
                    } else if (small[0] > Big_Rhino[i][1]) {
                        // std::cout << "拉拉 " << small[0] << " " << small[1] << "   " << Big_Rhino[i][0] << " " << Big_Rhino[i][1] << std::endl;
                        // std::cout << "卡在这里嘛" << std::endl;
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
        // std::cout << "几个阮梅 " << GroupAll <<s std::endl;
    }
    return chr_pointer_coverage;
}


FileSplit thread_all_read_sam_files(const std::string& sam_file_path, 
                                    std::vector<std::string>& samfilevec,
                                    const int& numThreads, 
                                    const std::string& outputPath, 
                                    const int& SJDistance, 
                                    std::unordered_map<std::string, std::string>& fasta_ref) {
    std::vector<std::string> sam_file_vec;
    if (samfilevec.size() == 1) {
        sam_file_vec.push_back(sam_file_path);
    } else {
        for (const auto& eachFile:samfilevec) {
            std::string chunkFilePath = joinPath(sam_file_path, eachFile);
            sam_file_vec.push_back(chunkFilePath);
        }
    }
    // FileSplit 最终输出的变量;
    FileSplit BigBang;
    std::map<int, std::map<std::string, std::vector<std::array<int,2>>>> File_chr_coverage;
    std::map<int, std::map<std::array<int,2>, std::array<std::streampos,2>>> File_group_pointer;
    std::vector<std::streampos> startposVec;
    std::vector<std::streampos> endposVec;
    std::map<std::string, std::set<std::string>> All_File_All_Barcode;
    if (sam_file_vec.size() != 0) {
        // 读取每一个sam文件;
        ThreadPool Bigfilepool(numThreads);
        std::vector<std::future<void>> myJobs;
        for (int samFileNumber = 0; samFileNumber < sam_file_vec.size(); ++samFileNumber) {
            
            startposVec.clear();
            endposVec.clear();
            // 创建文件夹输出小文件;
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
            // 每个块大小;
            std::streampos chunkSize = fileSize / numThreads;
            // 提前准备好线程开始结束;
            // 打开文件;
            std::ifstream one_sam(sam_file_vec[samFileNumber]);
            // 判断文件是否打开;
            if (!one_sam.is_open())	{
                std::cout << "We can't open SAM file !" << sam_file_vec[samFileNumber] << std::endl;
                exit(EXIT_FAILURE);
            }
            // 获取每个部分的开始和结尾指针;
            for (int i = 0; i < numThreads; ++i) {
                std::streampos startPos = i * chunkSize;
                std::streampos endPos = (i == numThreads - 1) ? (fileSize) : static_cast<std::streampos>((i + 1) * chunkSize);
                // 确保每次的startPos是在行首;
                if (i > 0) {
                    startPos = findNextLineStart(one_sam, startPos);
                }
                startposVec.push_back(startPos);
                endposVec.push_back(endPos);
            }
            // 关闭文件;
            one_sam.close();
            // 创建线程池, 每个部分拆成小文件;
            // 每个线程开始;
            myJobs.clear();
            // std::cout << "go to thread!" << std::endl;
            for (int i = 0; i < numThreads; ++i) {
                myJobs.emplace_back(Bigfilepool.enqueue([&,i]() {
                    processChunk(
                        sam_file_vec[samFileNumber],
                        startposVec[i], 
                        endposVec[i],
                        chunkFilePath, 
                        i, 
                        SJDistance, 
                        std::ref(fasta_ref)
                    );
                }));
            }
            for (auto& future : myJobs) {
                future.get();  // 等待每个任务完成
            }
            std::cout << "^-^ [" << sam_file_vec[samFileNumber] << "] All threads are finished generating small files! ^-^" << std::endl;
            
            // 直到这里, 上面完全没错; ★ ★ ★
            // 合并所有的 <Read_x.txt> 小文件;
            std::cout << "^-^ Start of merge small files ! ^-^" << std::endl;
            BigBang = Merge_Read_Small_Files(chunkFilePath, samFileNumber, numThreads);
            std::cout << "^-^ End of merge small files ! ^-^" << std::endl;
            
            // 如果大于一个文件, 每次将coverage放到一起;
            if (sam_file_vec.size() > 1) {
                // 将整个coverage汇总一下; 准备遍历结束之后合并成一个;
                File_chr_coverage[samFileNumber] = BigBang.chr_coverage;
                File_group_pointer[samFileNumber] = BigBang.coverage2pos;
                All_File_All_Barcode[std::to_string(samFileNumber)] = BigBang.file_barcodeSet[std::to_string(samFileNumber)];
            }
            // 等于1的时候, 不用管;
        }
        if (sam_file_vec.size() > 1) {
            std::cout << "^-^ The number of sam files is greater than one. Large files need to be merged. ^-^" << std::endl;
            // 通过一个函数合并所有的replicate的范围;
            // 合并所有文件的范围;
            std::map<std::string, std::vector<std::array<int,2>>> ChrCoverage = Merge_Read_Interval(File_chr_coverage);

            // 得到合并的区域在每个文件中的小区域, 便于生成大文件;
            std::unordered_map<std::string, std::map<std::array<int,2>, std::map<int, std::array<std::streampos,2>>>> ChrCoverage_SmallPointer = 
                                                                        get_pointers(ChrCoverage, File_chr_coverage, File_group_pointer);
            // 合并大文件时, 线程锁写入的变量;
            std::vector<std::streampos> Allreads_pointer;
            std::vector<int> Allreads_group_number;            
            // 输出的文件名;
            std::string FinallyFile = joinPath(outputPath, "MergedRead.txt");
            std::ofstream FinallyReadInform(FinallyFile, std::ios::trunc); 
            std::streampos Readpos = FinallyReadInform.tellp();
            Allreads_pointer.push_back(Readpos);

            // 合并所有大文件;
            // 创建并分配任务到线程池;
            myJobs.clear();
            int new_group = 0;
            // 分配任务; 重新分配group;
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
            //最终文件的名称;
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

// 切割最终文件每一行;
Line_Split MakeLineSplit(const std::string& s, char delimiter){
    // std::cout << s << std::endl;
    // 初始化返回值结构体;
    Line_Split read_result = {};
    // 字符串流;
    std::istringstream iss(s);
    // 定义字符串流中的每一段字符串;
    std::string token;
    // std::cout << iss << std::endl;
    // 定义一个number, 提取特定有用的信息;
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
            read_result.read_barcode = token;
        } else if (number > 10) {
            secondProcess.push_back(std::stoi(token));
        }
    }
    // std::cout << std::endl;
    // 两次处理;
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
    std::unordered_map<std::string, std::string> GroupReadBarcodes;
    std::unordered_map<std::string, std::string> GroupReadFiles; //Group每条reads的所属文件; (公共, single exon 与 Multi exon都有);
    std::unordered_map<std::string, std::vector<std::array<int,2>>> GroupReadSjs; //Group每条reads的SJs;
    std::unordered_map<std::string, std::array<int,2>> GroupReadCoverage; //Group每条reads的开头结尾;
    // std::map<std::array<int,2>, int> GroupSjs; //Group中根据Reads得到的高置信SJ;
    std::unordered_map<std::string, std::vector<int>> GroupSigns; //Group中所有reads的信号;
    std::unordered_map<std::string, std::array<int,2>> GroupSingleExon; // single exon的reads;
};

GroupInformation knowGroupInformation(std::streampos& startpos, 
                                    std::streampos& endpos, 
                                    const std::string& sam_file_path, 
                                    const int& Sj_Support_Read_Number) {
    GroupInformation groupinformation;
    std::ifstream FinalSamFile(sam_file_path);
    FinalSamFile.seekg(startpos);
    // 记录当前, 也就是文件刚开始的位置;
    std::streampos Current_Position = FinalSamFile.tellg();
    std::streampos Last_Position = Current_Position;
    std::string line;
    // 读取整个Group;
    while (getline(FinalSamFile, line)) {
        // 记录处理到的位置;
        Last_Position = Current_Position;
        Current_Position = FinalSamFile.tellg();
        if (Last_Position < endpos || Last_Position == FinalSamFile.eof()) {
            // 拆解每一行;
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



// 从注释中定位Group所在的注释;
/*
    函数的输入一是这条染色体上所有注释的开头和结尾; --- AnnoCoverage;
    二是整个这个cluster的表示范围; --- Group_Coverage;
    函数的输出是这个范围内的注释的SJ列表;
*/
struct GroupAnnotation
{
    std::unordered_map<std::string, std::vector<std::array<int,2>>> Group_Annotations;
    std::set<std::string> Group_GeneSet;
};

GroupAnnotation get_group_annotation(std::unordered_map<std::string, std::array<int,2>>& AnnoCoverage, 
                                std::unordered_map<std::string, std::vector<std::array<int,2>>>& AnnoisoSJ,
                                std::array<int,2>& Group_Coverage){
    //输出的数据结构就是unordered_map;
    GroupAnnotation ThisGroupAnnotations;
    std::string extracted;

    // int a = 0;
    //对当前的注释遍历;
    for (const auto& EachAnno:AnnoCoverage) {
        if ((EachAnno.second[0] >= Group_Coverage[0]) && (EachAnno.second[1] <= Group_Coverage[1])){
            //满足条件的存储一下;
            ThisGroupAnnotations.Group_Annotations[EachAnno.first] = AnnoisoSJ[EachAnno.first];
            // a++;
            size_t pos = EachAnno.first.find("|");
            if (pos != std::string::npos) {
                extracted = EachAnno.first.substr(0, pos);
            }
            ThisGroupAnnotations.Group_GeneSet.insert(extracted);
        }
    }
    return ThisGroupAnnotations;
}

// 得到Group范围内的单个exon的注释;
std::unordered_map<std::string, std::array<int,2>> get_group_single_exon_annotation(
                                    std::unordered_map<std::string, std::array<int,2>>& AnnoSE, 
                                    std::array<int,2>& Group_Coverage) {
    //输出的数据结构就是unordered_map;
    std::unordered_map<std::string, std::array<int,2>> ThisGroupSingleExonAnnotations;
    //对当前的注释遍历;
    for (const auto& EachAnno:AnnoSE) {
        if ((EachAnno.second[0] >= Group_Coverage[0]) && (EachAnno.second[1] <= Group_Coverage[1])){
            //满足条件的存储一下;
            ThisGroupSingleExonAnnotations[EachAnno.first] = EachAnno.second; 
        }
    }
    return ThisGroupSingleExonAnnotations;
}

// 得到Group范围内的单个exon及其对应的reads;
std::unordered_map<std::string, std::vector<std::string>> get_group_single_exon_reads(
                                    std::unordered_map<std::string, std::array<int,2>>& ThisGroupSingleExonAnnotations,
                                    std::unordered_map<std::string, std::array<int,2>>& GroupSingleExonReads,
                                    int& Edge) {
    // 输出这个注释对应的reads;
    std::unordered_map<std::string, std::vector<std::string>> SingleReads;
    // 初始化;
    for (const auto& eachAnno:ThisGroupSingleExonAnnotations) {
        SingleReads[eachAnno.first] = {};
    }
    // 开始对每个single exon 循环;
    for (const auto& EachSing:GroupSingleExonReads) {
        // 每个注释;
        for (const auto& EachAnno:ThisGroupSingleExonAnnotations) {
            if ((EachSing.second[1]-EachAnno.second[1]<Edge) && (EachAnno.second[0]-EachSing.second[0]<Edge)) {
                SingleReads[EachAnno.first].push_back(EachSing.first);
                break;
            }
        }
    }
    return SingleReads;
}

// Single exon结果写出;
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
    // 最终的结果还是return出去吧;
    // 最后一起写出来;
    std::set<std::string> thisFileBarcode;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> thisIsoform_File_Barcode2count;
    std::string first_part;
    std::string second_part;

    // 以上初始化最终结果;
    // 首先看看有几个文件哈;
    if (Updated_Files.is_open()) {
        // 首先遍历已经根据注释分好的组;
        if (singleexonwithreads.size() > 0) {
            for (const auto& eachAnno:singleexonwithreads) {
                // 可以写每个转录本的gtf;
                if (eachAnno.second.size() > 0) {
                    size_t pos = eachAnno.first.find('|');
                    if (pos != std::string::npos) {
                        // 根据 "|" 将字符串分割为两部分
                        first_part = eachAnno.first.substr(0, pos);
                        second_part = eachAnno.first.substr(pos + 1);
                    }
                    {
                        std::unique_lock<std::mutex> lock(updatedGtfMutex);
                        Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "transcript" << '\t' << singleexongroupannotation[eachAnno.first][0] << '\t' << singleexongroupannotation[eachAnno.first][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[eachAnno.first] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\";" << '\n';
                        Updated_Files << chrname << '\t' << "annotated_isoform" << '\t' << "exon" << '\t' << singleexongroupannotation[eachAnno.first][0] << '\t' << singleexongroupannotation[eachAnno.first][1] << '\t' << "." << '\t' << GTF_Transcript_Strand[eachAnno.first] << '\t' << '.' << '\t' << "gene_id \"" << first_part << "\"; transcript_id \"" << second_part << "\"; exon_number \"" << 0 << "\";" << '\n';
                    }
                    // 每个注释转录本开始时的归零;
                    thisIsoform_File_Barcode2count.clear();
                    for (int i = 0; i < FileNo; i++) {
                        std::string fileIndex = std::to_string(i);
                        thisIsoform_File_Barcode2count[fileIndex] = {};
                        for (const auto& eachBar:AllFileBarcode[fileIndex]) {
                            thisIsoform_File_Barcode2count[fileIndex][eachBar] = 0;
                        }
                    }
                    // 开始数数;
                    for (const auto& eachRead:eachAnno.second) {
                        std::string thisfile = groupreadfiles[eachRead];
                        std::string thisbarcode = groupreadbarcodes[eachRead];
                        thisIsoform_File_Barcode2count[thisfile][thisbarcode] = thisIsoform_File_Barcode2count[thisfile][thisbarcode] + 1;
                        // 下面可以写每个read的Trace;
                        {
                            std::unique_lock<std::mutex> lock(traceMutex);
                            Trace << eachRead << '\t' << "single_exon" << '\t' << second_part << '\t' << first_part << '\t' << thisbarcode << '\t' << thisfile << '\n'; 
                        }
                    }
                    // 一个isoform完毕;
                    // 写入到每个文件中;
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
                    } // 这一行写完了; 这个isoform;
                }
            }
        }            
    }    
}



void merge_high_sj(std::unordered_map<std::string, std::vector<std::array<int,2>>>& group_annotations, std::unordered_map<std::string, std::string>& transcript_strand, std::map<std::array<int,2>, int>& Grouphighsjs) {
    // 注释的所有sj, 所在的链, 高置信区间;
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



//10. 根据sam文件中每条reads的SJ chain聚类成不同的cluster;
/*  输入是read名称和read的SJ的vector区间的map;
    输出是聚类成的cluster的map, key是哈希值, value是包含的read名称;
*/

//10.1 函数用于将包含不同的向量分成不同的cluster, 使用unordered_map保存, key是哈希值, value是reads的名称得到的vector;
std::unordered_map<std::size_t, std::vector<std::string>> classifyReadsVec(std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs) {
    std::unordered_map<std::size_t, std::vector<std::string>> ClustersReads; //用于存储cluster;
    std::map<std::vector<std::array<int,2>>, std::size_t> readsj2sizet;
    //提取每个map中存储的每个reads的SJ的小区间;
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
    //首先对所有的高质量的cluster遍历;
    //SIRV可以看所在的是哪个基因;
    std::vector<std::array<int,2>> each_read_SJs; //reads的所有SJ;

    //新生成的变量;
    std::vector<std::string> statis_FSM; //满足条件的FSM的名称;
    std::vector<std::string> statis_ISM; //满足条件的ISM的名称;
    std::vector<int> ISM_flag; //ISM中每个SJ的索引;
    // int count = 0; //满足条件的ISM个数;
    std::string FSM_name; //保存输出变量时FSM的名称;
    std::vector<int> minFSM; //每个FSM与reads的距离;
    int minmin = 0; //FSM最小位置的索引;

    std::string first_part;
    std::string second_part; 

    //遍历cluster;
    for (auto it = ClustersReads.begin(); it != ClustersReads.end(); ++it) {
        //每一项的key表示这个cluster的SJ; value表示所有的reads的名称;
        each_read_SJs = AllSJs[(it->second)[0]];
        
        //遍历注释;
        statis_FSM.clear(); //判断FSM的, 判断多少个, 一般情况下就是1个;
        statis_ISM.clear();
        for (const auto& pair : GroupAnno){
            //pair.first表示基因名称如SIRV101, pair.second表示注释的所有SJ;
            if (pair.second == each_read_SJs){
                //这个判断条件表示 注释的SJ是否等于cluster的SJ;
                statis_FSM.push_back(pair.first);
            } else {
                //不相等的情况;
                //注释SJ的个数 必须大于 clusterSJ的个数;
                if (pair.second.size() > each_read_SJs.size()){
                    //接下来就看看是不是ISM; 不需要知道是谁的ISM, 后续还会追究, 所以通过数字就可以了解;
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
                //这里说明是这个cluster, 比当前的注释转录本的SJ个数要多;
            }
        }
        //遍历一遍注释后;继续看满足条件的个数;
        //先看FSM;
        if (statis_FSM.size() == 1){
            FSM_name = statis_FSM[0];
            FSMISMO.FSM[FSM_name].first = each_read_SJs;
            FSMISMO.FSM[FSM_name].second.insert(FSMISMO.FSM[FSM_name].second.end(), ClustersReads[it->first].begin(), ClustersReads[it->first].end());
            //将这个cluster中的每一条reads都输出到追踪文件;
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
            
            //这种情况就需要对每条reads分别看开头和结尾;
            //初始化将数量标记为0;
            for (const auto& refx : statis_FSM){
                FSM_name = refx;
                FSMISMO.FSM[FSM_name].first = each_read_SJs;
                FSMISMO.FSM[FSM_name].second = {};
            }
            for (const auto& readx : ClustersReads[it->first]) { //readx是每条read的名称;
                minFSM.clear();
                for (const auto& refx : statis_FSM){
                    minFSM.push_back(std::abs(AllReadBE[readx][0] - GroupAnnoBE[refx][0]) + std::abs(AllReadBE[readx][1] - GroupAnnoBE[refx][1]));
                }
                //每条reads的所有候选FSM计算完毕;
                minmin = std::min_element(minFSM.begin(), minFSM.end()) - minFSM.begin();
                //因此这条reads的归于索引为minmin的isoform;
                FSM_name = statis_FSM[minmin];
                FSMISMO.FSM[FSM_name].second.push_back(readx); //有个小bug, 如果都一样的话, 只会给前一个;
                //将这个cluster中的每一条reads都输出到追踪文件;
                size_t pos = FSM_name.find('|');
                if (pos != std::string::npos){
                    first_part = FSM_name.substr(0, pos);
                    second_part = FSM_name.substr(pos + 1);
                }
                {
                    std::unique_lock<std::mutex> lock(traceMutex);
                    Trace << readx << '\t' << "FSM" << '\t' << second_part << '\t' << first_part << '\t' << GroupReadBarcode[readx] << '\t' << GroupReadFile[readx] << '\n'; 
                }                
            } //分配结束;

        } else{ //只有FSM没有的情况下, 才继续看ISM;
            //然后看ISM;
            if (statis_ISM.size() == 0){
                //这种的情况是其他类;
                FSMISMO.Others[it->first] = ClustersReads[it->first];
            } else{
                //这种的情况是ISM类, 不用管是哪种ISM; //如果后续改识别ISM类型的时候, 更改这里; ★★★★★★★
                FSMISMO.ISM[it->first] = ClustersReads[it->first];
                //这些ISM是由哪些FSM得来的, 有的FSM可能已经有number, 有的没有number的就记为0;
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
    //所有的需要进一步进行识别的Cluster, 以及他们对应的数量;
    std::map<size_t, std::vector<std::string>> HighConClusters;
    //淘汰的Cluster及其他们对应的数量;
    std::map<size_t, std::vector<std::string>> LowConClusters;
    std::unordered_map<size_t, std::string> HighStrand;
};

HighLowClusters get_HighLow_clusters(std::unordered_map<std::size_t, std::vector<std::string>>& Others_Clusters, 
                                    std::unordered_map<std::string, std::vector<std::array<int,2>>>& AllSJs, 
                                    std::unordered_map<std::string, std::vector<int>>& readsigns,
                                    const int& sjsupportnumbers){
    //结构体初始化;
    HighLowClusters HighLowC = {};
    std::map<std::size_t, std::vector<std::string>> onlyHighC;
    std::unordered_map<size_t, std::string> onlyHighStrand;
    //在for循环的时候要删除;
    std::map<std::size_t, std::vector<std::string>> onlyHighC_copy;
    std::unordered_map<size_t, std::string> onlyHighStrand_copy;

    //读取的变量;
    std::vector<std::string> each_cluster; //每个cluster的read的名称;
    std::size_t hashname; //每个cluster的名称;
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
        hashname = it->first; //Cluster对应的哈希值;
        each_cluster = it->second; //Cluster对应的所有read名称的vector;
        //首先判断cluster中reads的数量, 一个的很多;
        if (each_cluster.size() > 1){
            // 多于1个的继续考虑;
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
                // std::cout << "Low类不满足>2" << std::endl;
            }
        } else {
            //一个类里面如果等于1个直接就不需要考虑了, 保存到Low的里面, 后面等回收;
            HighLowC.LowConClusters[hashname] = Others_Clusters[hashname];
            // std::cout << "Low类小于1" << std::endl;
        }
    }
    // std::cout << onlyHighC.size() << " High个  ***  " << HighLowC.LowConClusters.size() << " Low个" << std::endl;
    onlyHighC_copy = onlyHighC;
    onlyHighStrand_copy = onlyHighStrand;
    //还有一步, 将短的cluster去掉, 留下长的cluster;
    std::vector<std::array<int,2>> SJs1;
    std::vector<std::array<int,2>> SJs2;
    for (const auto& clus1: onlyHighC) {
        SJs1 = AllSJs[clus1.second[0]];
        for (const auto& clus2: onlyHighC){
            SJs2 = AllSJs[clus2.second[0]];
            if (SJs2.size()>SJs1.size()){
                if (std::includes(SJs2.begin(), SJs2.end(), SJs1.begin(), SJs1.end())){
                    //这种情况是clus1是clus2的子集;
                    //判断成功连续就可以了;
                    pos1 = std::find(SJs2.begin(), SJs2.end(), SJs1[0]) - SJs2.begin();
                    pos2 = std::find(SJs2.begin(), SJs2.end(), SJs1[SJs1.size()-1]) - SJs2.begin();
                    if (pos2-pos1 == SJs1.size()-1){
                        //这种情况说明1是2的子集;
                        HighLowC.LowConClusters[clus1.first] = clus1.second;
                        onlyHighC_copy.erase(clus1.first);
                        onlyHighStrand_copy.erase(clus1.first);
                    }
                }
            }
        }
    }
    // std::cout << "" << onlyHighC.size() << std::endl;
    //将最终的结果保存;
    for (const auto& eachHighCluster:onlyHighC_copy){
        HighLowC.HighConClusters[eachHighCluster.first] = eachHighCluster.second;
        HighLowC.HighStrand[eachHighCluster.first] = onlyHighStrand_copy[eachHighCluster.first];
    }
    //输出结果;
    // std::cout << HighLowC.HighConClusters.size() << " High个  ***  " << HighLowC.LowConClusters.size() << " Low个" << std::endl;
    return HighLowC;
}



//14.1 求两个集合之间的交集;
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
    // std::cout << "结果是:" << result << std::endl;
    return result;
}

//14.2 求两个集合之间的并集;
int IntervalMerge(std::vector<std::array<int,2>> firstList, std::vector<std::array<int,2>> secondList) {
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
    // std::cout << "结果是:" << result << std::endl;
    return result;
}

// 14.追加
// 从过滤掉的reads中判断是不是FSM;
void get_filtered_FSM(std::map<std::size_t, std::vector<std::string>>& LowReads, 
                      std::unordered_map<std::string, std::vector<std::array<int,2>>>& GroupAnno, 
                      std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>>& AllFSM, 
                      std::unordered_map<std::string, std::vector<std::array<int,2>>>& AllSJs, 
                      std::unordered_map<std::string, std::string>& GroupReadFile,
                      std::unordered_map<std::string, std::string>& GroupReadBarcode,
                      std::ofstream& Trace){
    //需要保存的变量;
    std::unordered_map<std::string, std::unordered_map<std::size_t, std::vector<std::string>>> RecycleAnnoName_value;
    std::set<std::size_t> DeleteItem;
    //每次计算出相差的距离;
    int RecycleLength = 0;
    int count = 0;
    std::string minAnnoName;
    int dis = 10000;
    //临时变量, 保存每个cluster潜在的多个FSM;
    std::unordered_map<std::string, int> potential_FSM;
    std::vector<std::array<int,2>> LowReadSJVec;

    std::string first_part;
    std::string second_part;

    //需要循环每一个过滤的cluster, 每个cluster可能有几个潜在的FSM;
    for (const auto& each_Readcluster:LowReads){
        LowReadSJVec = AllSJs[each_Readcluster.second[0]];
        // std::cout << each_anno.first << std::endl;
        //循环所有的过滤出去的reads_cluster;
        potential_FSM.clear();
        for (const auto& each_anno:GroupAnno){
            //如果注释与reads的cluster相同长度才可以计算, 为了减少计算的数量;
            if (each_anno.second.size() == LowReadSJVec.size()){
                RecycleLength = IntervalMerge(each_anno.second, LowReadSJVec) - IntervalIntersection(each_anno.second, LowReadSJVec);
                if (RecycleLength < 15*LowReadSJVec.size()){
                    //如果距离符合要求;
                    potential_FSM[each_anno.first] = RecycleLength;
                }
            }
        }
        //对于每个cluster查看有几个潜在的FSM;
        //如果只有1个是最好办的情况;
        if (potential_FSM.size() == 1){
            //首先将这个cluster标记;
            DeleteItem.insert(each_Readcluster.first);
            //输出到追踪文件;
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

            //查看这个潜在的FSM是不是已经出现过;
            if (RecycleAnnoName_value.find((*potential_FSM.begin()).first) != RecycleAnnoName_value.end()){
                //如果已经出现过, 只需要在原来的基础上加上数值;
                RecycleAnnoName_value[(*potential_FSM.begin()).first][each_Readcluster.first] = each_Readcluster.second;
            } else{
                //如果没有出现过;
                RecycleAnnoName_value[(*potential_FSM.begin()).first] = {};
                RecycleAnnoName_value[(*potential_FSM.begin()).first][each_Readcluster.first] = each_Readcluster.second;
            }
        //如果有多个, 则取距离最小的
        } else if (potential_FSM.size() > 1) {
            //首先将这个cluster标记;
            DeleteItem.insert(each_Readcluster.first);
            //找到距离最小的;
            dis = 100000;
            for (const auto& Po:potential_FSM){
                if (Po.second < dis){
                    dis = Po.second;
                    minAnnoName = Po.first;
                }
            }
            //输出追踪文件;
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

            //查看这个潜在的FSM——minAnnoName是不是已经出现过?
            if (RecycleAnnoName_value.find(minAnnoName) != RecycleAnnoName_value.end()){
                //如果已经出现过, 只需在原来的基础上加上数值;
                RecycleAnnoName_value[minAnnoName][each_Readcluster.first] = each_Readcluster.second;
            } else{
                //如果没有出现过;
                RecycleAnnoName_value[minAnnoName] = {};
                RecycleAnnoName_value[minAnnoName][each_Readcluster.first] = each_Readcluster.second;
            }
        }
        //整个cluster完毕;
    }   
    // int new_count = 0;
    std::vector<std::string> newName;
    //回收得到的注释需要重新添加到FSM中;
    // std::cout << "回收之前的注释个数:" << AllFSM.size() << std::endl;
    for (const auto& newFSM:RecycleAnnoName_value){
        //如果之前的FSM中已经存在;
        // new_count = 0;
        newName.clear();
        for (const auto& allCluster:newFSM.second){ //每个类;
            // new_count = new_count + allCluster.second.size();
            newName.insert(newName.end(), allCluster.second.begin(), allCluster.second.end());
        }

        if (AllFSM.find(newFSM.first) != AllFSM.end()){
            //只需要更新数值;
            // AllFSM[newFSM.first].second = AllFSM[newFSM.first].second + new_count;
            AllFSM[newFSM.first].second.insert(AllFSM[newFSM.first].second.end(), newName.begin(), newName.end());
        } else{
            //如果之前的FSM中不存在;
            //所有的都更新;
            AllFSM[newFSM.first].first = GroupAnno[newFSM.first];
            // AllFSM[newFSM.first].second = new_count;
            AllFSM[newFSM.first].second.insert(AllFSM[newFSM.first].second.end(), newName.begin(), newName.end());
        }
    }
    // std::cout << "回收之后的注释个数:" << AllFSM.size() << std::endl;
    //回收去掉的cluster需要舍弃;
    // std::cout << "回收的cluster个数:" << DeleteItem.size() << std::endl;
    // std::cout << "低置信的cluster数目:" << LowReads.size() << std::endl;
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
    // FSM ISM 根本没必要校正区间; 真正需要校正的是Low的;
    // 后面有一步回收, 能不能达到相同效果? 
    SCC.ClusterCoverage = get_every_cluster_begin_end(groupCluster, groupreadcoverage);
    FsmIsmOthers = get_FSM_and_others_sc(groupCluster, groupannotations, AnnoCoverage, groupreadcoverage, groupreadsjs, groupreadfiles, groupreadbarcodes, traceFilePath);
    Others2HighLow = get_HighLow_clusters(FsmIsmOthers.Others, groupreadsjs, groupreadsigns, Sj_Support_Number);
    get_filtered_FSM(Others2HighLow.LowConClusters, groupannotations, FsmIsmOthers.FSM, groupreadsjs, groupreadfiles, groupreadbarcodes, traceFilePath);
    // 将有用结果输出;
    SCC.FSM = FsmIsmOthers.FSM;
    SCC.ISM = FsmIsmOthers.ISM;
    SCC.HighConClusters = Others2HighLow.HighConClusters;
    SCC.HighStrand = Others2HighLow.HighStrand;
    // std::cout << "FSM:" << SCC.FSM.size() << "  ISM:" << SCC.ISM.size() 
    //           << "  High:" << SCC.HighConClusters.size() << std::endl;
    return SCC;
}


//14.计算距离矩阵;
//14.3 求两个集合之间没有关系的两个部分;
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

struct DistanceInform
{
    std::unordered_map<int, std::set<int>> NodeDistanceSet;
    //序号名称和节点名称的字典;
    std::unordered_map<int, std::vector<std::array<int,2>>> Index2Unknown;
    std::unordered_map<int, std::string> Index2Anno;
    std::unordered_map<int, int> Index2Count;
    std::unordered_map<int, std::size_t> Index2hashname;
    std::unordered_map<int, std::string> Index2novelname; //后面才用到;
};

//14.3 计算距离矩阵; 
DistanceInform get_distance_matrix(std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupannotations, 
                                   std::map<std::size_t, std::vector<std::string>>& HighClusters, 
                                   std::unordered_map<std::string, std::vector<std::array<int,2>>>& groupreadsjs, int DistanceG, 
                                   std::unordered_map<std::string, std::pair<std::vector<std::array<int,2>>, std::vector<std::string>>>& AllFSM){
    //输出结果的二维矩阵;
    DistanceInform Graph_Init_inform;
    //初始化
    for (int i = 0; i < groupannotations.size()+HighClusters.size(); ++i){
        Graph_Init_inform.NodeDistanceSet[i] = {};
    }

    std::vector<std::array<int,2>> SJVec1;
    std::vector<std::array<int,2>> SJVec2;
    int LineElement = 0; //距离;
    int NameNumber = 0;
    int HighNumber = groupannotations.size();

    //首先计算的一部分是注释的isoform与候选的isoform之间的距离; 相当于整个矩阵的part2;
    //如果注释这部分没有;
    if (groupannotations.size() != 0){
        for (const auto& EachAnno : groupannotations){
            Graph_Init_inform.Index2Anno[NameNumber] = EachAnno.first;
            //提取这个注释是否在FSM中;
            auto it = AllFSM.find(EachAnno.first);
            //检查是否找到,如果在FSM中, 则将数量提取, 否则, 数量记为0;
            if (it != AllFSM.end()) {
                Graph_Init_inform.Index2Count[NameNumber] = (it->second).second.size();
            } else{
                Graph_Init_inform.Index2Count[NameNumber] = 0;
            }
            //清空每一行的距离重新赋值;
            if (HighClusters.size() != 0){
                //初始化;
                HighNumber = groupannotations.size();
                for (const auto& eachHC : HighClusters){
                    //从这里开始计算两两之间的距离大小;
                    SJVec1 = groupreadsjs[eachHC.second[0]];
                    LineElement = IntervalMerge(EachAnno.second, SJVec1) - IntervalIntersection(EachAnno.second, SJVec1) - findNonOverlappingIntervals(EachAnno.second, SJVec1);
                    // std::cout << LineElement << std::endl;
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
        //接下来计算的候选的isoform之间的距离;通过part2得到part24;
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

//15 构造图, 求出所有的极大团;
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
        //使用std::set_difference算法计算差集;
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


//16. 识别transcript; 识别;
//16.1 根据子图中元素个数排序的自定义函数; 
bool compareBySize(const std::vector<int>& a, const std::vector<int>& b) {
    // std::cout << "Comparing sizes: " << a.size() << " and " << b.size() << std::endl;
    if (a.empty() && b.empty()) {
        return false; // 两个向量都为空，返回相等或任何你需要的值
    }
    if (a.empty()) {
        return true; // 空向量被认为比非空向量小
    }
    if (b.empty()) {
        return false; // 非空向量被认为比空向量大
    }
    auto AminElement = std::min_element(a.begin(), a.end());
    auto BminElement = std::min_element(b.begin(), b.end());
    // 比较最小元素
    if (*AminElement != *BminElement) {
        return *AminElement < *BminElement; // 返回最小元素的比较结果
    }
    // 如果最小元素相等，按大小比较
    return a.size() < b.size();
}


//16 定义的结构体, 表示经过第一次识别后各个set中的元素,
// 
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
    //整个染色体中已知正确的节点;
    //定义结果, 识别为正确的和识别为错误的节点;
    //初始化结构体;
    DetectionResults Node_Results = {};

    //函数用到的中间变量;
    std::set<int> subgraph_known_true_nodes;
    std::set<int> subgraph_remain_nodes;
    int standard_read_number = std::numeric_limits<int>::max();
    int known_true_number = 0;
    int max_node_value = 0;


    //遍历一遍已知注释的点, 将其存入NodesKnownTrue的Set中;
    for (const auto& Anno : Number2Anno){
        Node_Results.NodesKnownTrue.insert(Anno.first);
    }

    //首先判断团的个数是否为0, 然后进行遍历;
    //如果团的个数不是0个, 那就从小到大开始遍历;
    if (CqVec.size() != 0){
        //按照子图中个数的多少对子图进行排序, 首先对节点少的
        if (CqVec.size() > 1) {
            std::sort(CqVec.begin(), CqVec.end(), compareBySize);
        }
        //遍历每个团;
        for (const auto& eachClique : CqVec){
            standard_read_number = std::numeric_limits<int>::max();

            //首先提取出每个团中的节点; 每个向量已经是提取好了每个团的节点;
            //获取团中已知为True的节点;
            //获取团中未知的节点;
            std::set<int> CliqueSet(eachClique.begin(), eachClique.end());
            subgraph_known_true_nodes.clear();
            std::set_intersection(CliqueSet.begin(), CliqueSet.end(), Node_Results.NodesKnownTrue.begin(), Node_Results.NodesKnownTrue.end(),
                          std::inserter(subgraph_known_true_nodes, subgraph_known_true_nodes.begin()));
            subgraph_remain_nodes.clear();
            std::set_difference(CliqueSet.begin(), CliqueSet.end(), subgraph_known_true_nodes.begin(), subgraph_known_true_nodes.end(),
                        std::inserter(subgraph_remain_nodes, subgraph_remain_nodes.begin()));

            
            //子图/团中有已知正确的点;
            if (subgraph_known_true_nodes.size() != 0){
                //首先看看有几个已知正确的节点, 选取正确的节点的最小的reads数作为衡量标准;
                //如果只有一个正确的节点;
                if (subgraph_known_true_nodes.size() == 1){
                    standard_read_number = Number2Count[(*subgraph_known_true_nodes.begin())];
                } else {
                    //如果有多个正确的节点;
                    for(const auto& each_True_node : subgraph_known_true_nodes){
                        if (Number2Count[each_True_node] < standard_read_number) {
                            standard_read_number = Number2Count[each_True_node];
                        }
                    }
                }
                //上面正确的准备完毕;
                //下面遍历未知的节点; 子图中有已知正确的节点, 但也有遗留的节点;
                if (subgraph_remain_nodes.size() != 0){ 
                    //有已知相邻的节点, 即使number为0也可以加入识别的;
                    for (const auto& eachknownnode:subgraph_known_true_nodes){
                        Node_Results.TrueNodeSet.insert(subgraph_known_true_nodes.begin(), subgraph_known_true_nodes.end());
                    }
                    //利用已知节点, 看一下其他的节点, 不然直接剔除;
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
                } else { // 子图中全是已知的节点, 没有遗留的节点;
                    //表明全是已知正确的节点; 整个子图全是已知注释的节点, 如果数量为0, 
                    //就不算识别出来的样本中的transcript;
                    for (const auto& eachknownnode:subgraph_known_true_nodes){
                        if (Number2Count[eachknownnode] > 0){
                            Node_Results.TrueNodeSet.insert(eachknownnode);
                        }
                    } //这一段容易在gtf中多出来为0的注释;
                    // known_true_number = 0;
                    // for (const auto& eachknownnode:subgraph_known_true_nodes){
                    //     // if (Number2Count[eachknownnode] != 0){
                    //     known_true_number++;
                    // }
                    // if (known_true_number != 0){
                    //     Node_Results.TrueNodeSet.insert(subgraph_known_true_nodes.begin(), subgraph_known_true_nodes.end());
                    // }
                }

            } else {
                //下面是子图/团中没有已知正确的节点;
                //如果整个图只有1个节点, 这种情况暂时认为是对的; (我认为需要更进一步求证;)
                if (subgraph_remain_nodes.size() == 1){
                    //只有一个节点, 暂时认为对;
                    if (Number2Count[*subgraph_remain_nodes.begin()] > 0) {
                        Node_Results.TrueNodeSet.insert(*subgraph_remain_nodes.begin());
                    }
                } else{
                    //如果子图中有多个节点, 就需要进一步考虑;
                    //考虑具有最多reads数的节点;
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
            //每个子图一次循环;
        }
        //严格 ------
        //将识别的False中的元素在Ture中检查, 如果有, 则在True中删除;
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
        // 对每个FSM循环;
        for (const auto& eachCluster:FsmIsmHigh.FSM) {
            std::string AnnoName = eachCluster.first;
            FileSpliceChains.File_FSM[AnnoName] = {};
            // 对这个transcript所属的每条reads进行循环;
            for (const auto& eachRead:eachCluster.second.second) {
                // 这条read所在的文件;
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                if (read_file == thisFileString) {
                    FileSpliceChains.File_FSM[AnnoName].push_back(eachRead);
                }
            }
            // 所有的reads循环结束;          
        }
    } // 结束;

    // 拆分ISM;
    if (FsmIsmHigh.ISM.size() > 0) {
        // 对每个ISM循环;
        for (const auto& eachCluster:FsmIsmHigh.ISM) {
            std::size_t NonName = eachCluster.first;
            FileSpliceChains.File_ISM[NonName] = {};
            // 对这个ISM所属的每条reads进行循环;
            for (const auto& eachRead:eachCluster.second) {
                // 这条read所在的文件;
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                if (read_file == thisFileString) {
                    FileSpliceChains.File_ISM[NonName].push_back(eachRead);
                }
            }
            // 所有的reads循环结束;
        }
    }

    // 拆分HighConfidenceClusters;
    if (FsmIsmHigh.HighConClusters.size() > 0) {
        // 对每个HighConCluster循环;
        for (const auto& eachCluster:FsmIsmHigh.HighConClusters) {
            std::size_t NonName = eachCluster.first;
            FileSpliceChains.File_HighConClusters[NonName] = {};
            // 对这个HighConCluster所属的每条reads进行循环;
            for (const auto& eachRead:eachCluster.second) {
                // 这条read所在的文件;
                std::string read_file = groupinformations.GroupReadFiles[eachRead];
                if (read_file == thisFileString) {
                    FileSpliceChains.File_HighConClusters[NonName].push_back(eachRead);
                }
            }
            // 所有的reads循环结束;            
        }
    }
    return FileSpliceChains;
}

std::string get_co_sj_gene(std::vector<std::array<int,2>>& novel_isoform, std::set<std::string>& gene_set,
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
                // gene set里面已经有很多isoform; 当前isoform的sj就是itssj;
                if (groupallgene.size() != 0) {
                    TempGeneName = get_co_sj_gene(itssj, groupallgene, gtf_gene2tx, groupannotations);

                    if (TempGeneName.size() != 0) {
                        first_part = TempGeneName;
                    } else {
                        // 这个时候就是多个exon的注释和单个exon的注释, 从里面选一个;
                        for (const auto& every_gene:groupallgene) {
                            if (Annogenecovergae[every_gene][0] <= BE[1] && BE[0] <= Annogenecovergae[every_gene][1]){
                                TempGene.push_back(every_gene);
                            }
                        }
                        //novel isoform没有候选的基因;
                        if (TempGene.size() == 0){
                            first_part = "NA";
                        // novel isoform只有1个候选的基因;
                        } else if (TempGene.size() == 1) {
                            first_part = TempGene[0];
                        //novel isoform有2个以上候选的基因;
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
                    // 这个时候只需要用单个exon的注释就可以, 从里面选一个;
                    // 如果也没有, 那么就是NA;
                    //单个exon的情况;
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
                // std::cout << first_part << std::endl;
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
    // // 打印一下;
    // for (const auto& aaa:FinalAnnotations.transcript2gene) {
    //     std::cout << aaa.first << ": " << aaa.second << std::endl;
    // }
    // for (const auto& aaa:FinalAnnotations.Transcript_Annotations) {
    //     std::cout << aaa.first << ": ";
    //     for (const auto& bbb:aaa.second.first) {
    //         std::cout << "[" << bbb[0] << "," << bbb[1] << "] ";
    //     }
    //     std::cout << aaa.second.second << std::endl;
    // }
    // std::cout << "猪八戒:Transcript_Annotations-" << FinalAnnotations.Transcript_Annotations.size() << std::endl;
    // std::cout << "孙悟空:transcript2gene-" << FinalAnnotations.transcript2gene.size() << std::endl;
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
    // HC
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
    //这个Barcode中所有错误的节点;
    FileBarcodeFalseNode ThisFile_FalseNumber;
    std::map<size_t, std::vector<std::string>> this_file_HC;
    auto it = barcodehc.find(Barcode);
    if (it == barcodehc.end()) {
        return ThisFile_FalseNumber;
    } else {
        this_file_HC = barcodehc[Barcode];
    }
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


//19. 量化部分;
//判断一个ISM类是否是这个注释的ISM;
int whether_isoform_part_is_not(std::vector<std::array<int,2>> AnnoSJ, std::vector<std::array<int,2>> ISMSJ){
    int ISM_Flag = 0;
    auto it1 = std::find(AnnoSJ.begin(), AnnoSJ.end(), ISMSJ[0]);
    if (it1 != AnnoSJ.end()){
        //计算出it1的索引;
        int it1_index = it1 - AnnoSJ.begin();
        // std::cout << it1_index << std::endl;
        auto it2 = std::find(AnnoSJ.begin(), AnnoSJ.end(), ISMSJ[ISMSJ.size()-1]);
        if (it2 != AnnoSJ.end()){
            //计算出it2的索引;
            int it2_index = it2 - AnnoSJ.begin();
            // std::cout << it2_index << std::endl;
            //如果两个索引之间相差的大小等于长度减1;
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
    //生成的变量, 识别的转录本按照顺序;
    std::vector<std::string> Order_Transcript_Name;
    //生成的变量, 最终的ISM指示矩阵;
    Eigen::MatrixXd known_ISM_matrix;
    std::vector<std::array<int,2>> AreadSJs;
    
    //建立索引和判断条件;
    int RowIndex = -1;
    int ColIndex = -1;
    int flag = 0;

    //每个注释包括, 名称, SJ序列和数量;
    //它不是有顺序的, 因此通过名称将其排序, 构造一个名称的vector;
    for (const auto& eachTransctipt:FinallyAnnotations.Transcript_Annotations) {
        Order_Transcript_Name.push_back(eachTransctipt.first);
    }
    Eigen::RowVectorXd EachClusterRowVector(Order_Transcript_Name.size());
    //存储一个列向量, 每个类的个数, 方便后续做乘法;
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

    //这一段得到ISM的指示矩阵;
    if (groupISM.size() != 0){
        RowIndex = -1;
        for (const auto& eachISM:groupISM){
            //对生成的每个行向量清空; 每个cluster的行向量;
            // std::cout << eachISM.second.size() << std::endl;
            AreadSJs = groupreadsjs[eachISM.second[0]];
            EachClusterRowVector.setZero();
            ColIndex = -1;
            ism_gtf_name.clear();
            first_part_set.clear();
            second_part_set.clear();            
            //对每个注释遍历;
            for (const auto& eachTransctipt:Order_Transcript_Name){
                ColIndex = ColIndex + 1;
                //如果ISM的个数大于注释, 则不可能是1;
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
            // 在这里加一步过滤最好了, ISM如果有一个特别近的, 还有一个特别远的, 特别远的就去掉吧;
            if (ColIndexVec.size() > 1) {
                if (max_ratio >= 0.5) {
                    for (const auto& ColIndex:ColIndexVec) {
                        std::string TransName = Order_Transcript_Name[ColIndex];
                        this_ratio = static_cast<double>(AreadSJs.size())/FinallyAnnotations.Transcript_Annotations[TransName].first.size();
                        if (this_ratio < max_ratio) {
                            EachClusterRowVector[ColIndex] = 0;
                            ism_gtf_name.erase(std::remove(ism_gtf_name.begin(), ism_gtf_name.end(), TransName), ism_gtf_name.end());
                            // std::cout << " 删除 " << TransName << std::endl;
                        }
                    }
                }
            } 
            if (ism_gtf_name.size() != 0){
                //只有这一行不全为0, 才能放入矩阵;
                RowIndex = RowIndex + 1;
                //调整矩阵大小以容纳第一个行向量;
                known_ISM_matrix.conservativeResize(RowIndex+1, Order_Transcript_Name.size());
                known_ISM_matrix.row(RowIndex) = EachClusterRowVector;
                //存储一个列向量, 每个类的个数, 方便后续做乘法;
                known_ISM_cluster_numbers.conservativeResize(RowIndex+1, 1);
                known_ISM_cluster_numbers(RowIndex, 0) = eachISM.second.size();

                //开始输出ISM的追踪文件;
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

    //这一段会得到False novel candidate的指示矩阵;
    //构造识别出的所有错误节点;
    std::set<int> MergeFalseSet(falsenodeset.begin(), falsenodeset.end());

    //初始化生成的指示矩阵;
    Eigen::MatrixXd False_novel_candidate_matrix(MergeFalseSet.size(), Order_Transcript_Name.size());
    //每个错误节点的所有邻居集合;
    std::set<int> FalseNode_NeighborSet;
    std::set<int> FalseNode_TrueSet;
    //对应正确节点的名称;
    std::string itsname;
    //要查找的子字符串;
    // std::string subStr;
    //构造一个行向量;
    Eigen::RowVectorXd FalseNode_RowVector(Order_Transcript_Name.size());
    RowIndex = -1;
    int countC = 0;
    
    if (MergeFalseSet.size() != 0){
        //构造识别出的所有正确节点;
        std::set<int> MergedTrueSet(truenodeset.begin(), truenodeset.end());
        
        //遍历所有识别为错误的节点;
        //找到每个错误节点的邻居节点 --- FalseNode_NeighborSet;
        for (const auto& FalseNode:MergeFalseSet){    
            RowIndex = RowIndex + 1;        
            //找到与他相邻的节点;
            FalseNode_NeighborSet.clear();
            FalseNode_TrueSet.clear();            
            FalseNode_NeighborSet = Disinform.NodeDistanceSet[FalseNode];

            //从相邻的节点中找到正确的节点 --- FalseNode_TrueSet;
            std::set_intersection(MergedTrueSet.begin(), MergedTrueSet.end(), FalseNode_NeighborSet.begin(), FalseNode_NeighborSet.end(),
            std::inserter(FalseNode_TrueSet, FalseNode_TrueSet.begin()));

            //如果有正确的邻居节点, 就需要根据顺序得到指示向量;
            //如果没有正确的邻居节点, 这种情况基本不符合, 可以设置一个bug;
            FalseNode_RowVector.setZero();
            ism_gtf_name.clear();

            if (FalseNode_TrueSet.size() != 0) {
                //对所有的正确节点遍历;
                for(const auto& EachTNode:FalseNode_TrueSet){
                    if (EachTNode < Disinform.Index2Anno.size()){
                        //说明这个邻居是有注释的转录本;
                        itsname = Disinform.Index2Anno[EachTNode];
                        //接下来找itsname在Order_Transcript_Name中的索引;
                        auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                        FalseNode_RowVector[itwhere] = 1;
                        ism_gtf_name.push_back(itsname);

                    } else
                    {
                        //这个邻居是novel的转录本;
                        itsname = Disinform.Index2novelname[EachTNode];
                        // std::cout << "哈哈 " << itsname << std::endl;
                        // std::cout << "Rows: " << FalseNode_RowVector.rows() << std::endl;
                        // std::cout << "Cols: " << FalseNode_RowVector.cols() << std::endl;
                        // std::cout << "Size: " << FalseNode_RowVector.size() << std::endl;
                        auto itwhere = std::find(Order_Transcript_Name.begin(), Order_Transcript_Name.end(), itsname) - Order_Transcript_Name.begin();
                        // std::cout << "这里 " << itwhere << std::endl;
                        FalseNode_RowVector[itwhere] = 1;
                        ism_gtf_name.push_back(itsname);
                    }
                }
                //将这一行加入到矩阵中;
                False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
            } 
            else
            {   
                countC = 0;
                //需要找邻居的邻居;
                for(const auto& neinode:FalseNode_NeighborSet) {
                    if (Disinform.NodeDistanceSet[neinode].size() != 0){
                        for (const auto& node:Disinform.NodeDistanceSet[neinode]){
                            if (node < Disinform.Index2Anno.size()){
                                //说明这个邻居是有注释的转录本;
                                itsname = Disinform.Index2Anno[node];
                                //接下来找itsname在Order_Transcript_Name中的索引;
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
                    //这种情况下就是邻居的邻居也没有一个好的, 那就找距离最小值了;
                    std::vector<int> houX = {};
                    int distance = 0;
                    for (const auto& isofrorm:Order_Transcript_Name){
                        distance = IntervalMerge(FinallyAnnotations.Transcript_Annotations[isofrorm].first, Disinform.Index2Unknown[FalseNode]) - IntervalIntersection(FinallyAnnotations.Transcript_Annotations[isofrorm].first, Disinform.Index2Unknown[FalseNode]);
                        houX.push_back(distance);
                    }
                    // 使用 std::min_element 找到最小元素的迭代器
                    auto min_it = std::min_element(houX.begin(), houX.end());

                    // 计算最小元素的索引
                    int min_index = std::distance(houX.begin(), min_it);
                    FalseNode_RowVector[min_index] = 1;
                    False_novel_candidate_matrix.row(RowIndex) = FalseNode_RowVector;
                    ism_gtf_name.push_back(Order_Transcript_Name[min_index]);
                    // std::cout << "*** Neighbors neighbors are wrong. That doesn't happen!!! ***" << std::endl;
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

    //存储一个列向量, 每个类的个数, 方便后续做乘法;
    Eigen::VectorXd False_novel_candidate_cluster_numbers(MergeFalseSet.size());
    RowIndex = -1;
    for (const auto& FalseNode:MergeFalseSet){
        RowIndex = RowIndex + 1;
        False_novel_candidate_cluster_numbers[RowIndex] = faslenode2count[FalseNode].size();
    }

    //最后输出结果矩阵和结果量;
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
 
// 最终EM算法的步骤;
void EM_Alg_Barcodes (std::unordered_map<std::string, std::unordered_map<std::string, double>>& filetranscriptnumber, 
                        IndicateFire& Indicate_Number, 
                        const std::string& Barcode) {

    if (Indicate_Number.Indicate_Matrix.rows() != 0 && Indicate_Number.Indicate_Matrix.cols() != 0) {
        // 初始化概率向量
        Eigen::VectorXd P_Col_init0 = Eigen::VectorXd::Constant(
            Indicate_Number.Order_Transcript_Name_Vector.size(),
            1.0 / Indicate_Number.Order_Transcript_Name_Vector.size() 
        );
        // 初始化;
        Eigen::MatrixXd Z(Indicate_Number.Indicate_Matrix.rows(), Indicate_Number.Indicate_Matrix.cols());
        Eigen::VectorXd P1(P_Col_init0.size());
        Eigen::MatrixXd AnnoN;
        int CountCyc = 1; // 初始化为 1，与原代码对齐
        double sum_abs_diff;

        do {
            // 计算 Z 并归一化
            Z = Indicate_Number.Indicate_Matrix.array().rowwise() * P_Col_init0.transpose().array();
            Z.array().colwise() /= (Indicate_Number.Indicate_Matrix * P_Col_init0).array();
            // 计算 AnnoN 并确保列向量形式
            AnnoN = ((Indicate_Number.Cluster_Number.transpose()) * Z).transpose().reshaped(Indicate_Number.Indicate_Matrix.cols(),1);
            // 更新概率 P1
            P1 = AnnoN / AnnoN.sum();
            // 计算收敛条件
            sum_abs_diff = (P1 - P_Col_init0).cwiseAbs().sum();
            // 更新概率向量
            P_Col_init0 = P1;
            CountCyc = CountCyc + 1;
        } while ((sum_abs_diff > 5e-2) || (CountCyc <= 10));

        //EM结束, 现在需要将这些注释都输出到一个文件里了;
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
                            std::vector<std::mutex>& AllIsoformMutexes) {
    // 输入是<barcode:<gene|isoform:counts>>
    // 想先得到<gene|isoform:<barcode:counts>>;
    // 首先初始化; 
    std::string transcript_name;
    std::string gene_name;
    std::unordered_map<std::string, std::unordered_map<std::string, double>> transcriptbarcodenumber;
    if (FinallyAnnotations.Transcript_Annotations.size() > 0) {
        for (const auto& eachAnno:FinallyAnnotations.Transcript_Annotations) {
            transcriptbarcodenumber[eachAnno.first] = {};
        }
        // 已知是 <barcode, <isoform, counts>>; 需要转化为<isoform, <barcode, counts>>;
        for (const auto& eachAnno:FinallyAnnotations.Transcript_Annotations) {
            std::string AnnoName = eachAnno.first;
            for (const auto& eachBar:filetranscriptnumber) {
                transcriptbarcodenumber[AnnoName][eachBar.first] = filetranscriptnumber[eachBar.first][AnnoName];
            }
        }
        // 输出到transcript文件; transcriptbarcodenumber已经构建完毕;
        for (const auto& eachAnno:transcriptbarcodenumber) {
            size_t pos = eachAnno.first.find('|');
            if (pos != std::string::npos) {
                transcript_name = eachAnno.first.substr(pos + 1);
            } else {
                transcript_name = eachAnno.first;
            }
            gene_name = FinallyAnnotations.transcript2gene[eachAnno.first];
            // 写入;
            {
                std::unique_lock<std::mutex> lock(AllIsoformMutexes[filenumber]);
                *(AllIsoformFilePath[filenumber]) << transcript_name << '\t' << gene_name;
                for (const auto& eachBar:barcodeset) {
                    double Counts = transcriptbarcodenumber[eachAnno.first][eachBar];
                    if ( Counts < 0.1 ) {
                        *(AllIsoformFilePath[filenumber]) << '\t' << 0;
                    } else {
                        *(AllIsoformFilePath[filenumber]) << '\t' << int(Counts);
                    }
                    
                }
                *(AllIsoformFilePath[filenumber]) << '\n';
            }
            // 这个isoform写完了;
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
                std::unordered_map<std::string, std::array<int,2>>& singleexonanno) {
    // 识别和定量;
    std::string chrchr = groupinform.chrName;
    std::string groupnumber = groupinform.GroupIndex;
    DistanceInform DMatrix_GraphNode = get_distance_matrix(groupanno.Group_Annotations, 
                                                           splicechainclass.HighConClusters, 
                                                           groupinform.GroupReadSjs, 
                                                           groupdistance, splicechainclass.FSM);
    // 最终输出;
    // 开始计算;
    if (DMatrix_GraphNode.NodeDistanceSet.size() != 0) {
        // 寻找最大团;
        std::vector<std::vector<int>> CliquesVector = Find_Maximal_Cliques(DMatrix_GraphNode.NodeDistanceSet);
        // 识别;
        DetectionResults NodeResults = Transcript_Detection(CliquesVector, 
                                                            DMatrix_GraphNode.Index2Anno, 
                                                            DMatrix_GraphNode.Index2Unknown, 
                                                            DMatrix_GraphNode.Index2Count);
        CliquesVector.clear();
        // splicechainclass.FSM.clear();
        // splicechainclass.ISM.clear();
        // 写到更新的注释中;
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
        // 量化前拆分;
        if (FileNo > 1) {
            // 多于一个sam文件;
            for (int k = 0; k < FileNo; k++) {
                // 数量
                Solvent_sc SpliceChainSolvent = get_Solvent_FsmIsmHigh(splicechainclass, k, groupinform);
                std::unordered_map<std::string, std::unordered_map<std::string, double>> File_k_TranscriptNumber = get_transcript_init(
                    NodeResults, DMatrix_GraphNode, chrchr, groupnumber, k, groupinform.GroupReadFiles, groupinform.GroupReadBarcodes,
                    SpliceChainSolvent, AllFile_BarcodeSet
                );
                // 获取这个文件中所有barcode的FSM, ISM, HighConCluster;
                Barcode_sc ThisFileBarcodeSpliceChain = get_Barcode_FsmIsmHigh(SpliceChainSolvent.File_FSM,
                                                                    SpliceChainSolvent.File_ISM,
                                                                    SpliceChainSolvent.File_HighConClusters,
                                                                    groupinform.GroupReadBarcodes);
                // 对每一个barcode量化;
                for (const auto& ThisBarcode:AllFile_BarcodeSet[std::to_string(k)]) {
                    // barcode错误节点;
                    FileBarcodeFalseNode This_File_False_Node = get_File_Barcode_False_node(NodeResults.FalseNodeSet,
                                                                            DMatrix_GraphNode.Index2hashname,
                                                                            ThisFileBarcodeSpliceChain.Barcode_HighConClusters,
                                                                            ThisBarcode);
                    // 量化;
                    // 量化第一步, 准备好所有变量, 初始化;
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
                    //量化第二步, EM算法;
                    EM_Alg_Barcodes(File_k_TranscriptNumber, InitFirefly, ThisBarcode); // 为了得到所有样本的结果;
                    // std::cout << "(=^_^=) Group " << groupinform.GroupIndex << " in file " << k << " completed quantification! (=^_^=) " << std::endl;                                         
                    ThisFileBarcodeSpliceChain.Barcode_FSM[ThisBarcode].clear();
                    ThisFileBarcodeSpliceChain.Barcode_ISM[ThisBarcode].clear();
                    ThisFileBarcodeSpliceChain.Barcode_HighConClusters[ThisBarcode].clear();
                }
                // 输出出来; File_k_TranscriptNumber对每个文件的这个变量操作; <barcode, <isoform, counts>>; 特别全;
                // 在函数中需要转化为 <isoform <barcode, counts>>; 然后输出到isoform counts文件中;
                write_transcript_file(File_k_TranscriptNumber, Finally_Annotations, AllFile_BarcodeSet[std::to_string(k)], k, IsoformFilePath, IsoformMutexes);
            }
            std::cout << "end one group!" << std::endl;

        } else {
            int k = 0;
            Solvent_sc SpliceChainSolvent = get_Solvent_FsmIsmHigh(splicechainclass, k, groupinform);
            std::unordered_map<std::string, std::unordered_map<std::string, double>> File_k_TranscriptNumber = get_transcript_init(
                NodeResults, DMatrix_GraphNode, chrchr, groupnumber, k, groupinform.GroupReadFiles, groupinform.GroupReadBarcodes,
                SpliceChainSolvent, AllFile_BarcodeSet);
            // 获取这个文件中所有barcode的FSM, ISM, HighConCluster;
            Barcode_sc ThisFileBarcodeSpliceChain = get_Barcode_FsmIsmHigh(SpliceChainSolvent.File_FSM,
                                                                SpliceChainSolvent.File_ISM,
                                                                SpliceChainSolvent.File_HighConClusters,
                                                                groupinform.GroupReadBarcodes);
            // 对每一个barcode量化;
            for (const auto& ThisBarcode:AllFile_BarcodeSet["0"]) {
                // barcode错误节点;
                FileBarcodeFalseNode This_File_False_Node = get_File_Barcode_False_node(NodeResults.FalseNodeSet,
                                                                        DMatrix_GraphNode.Index2hashname,
                                                                        ThisFileBarcodeSpliceChain.Barcode_HighConClusters,
                                                                        ThisBarcode); 
                // 量化;
                // 量化第一步, 准备好所有变量, 初始化;
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
                //量化第二步, EM算法;
                EM_Alg_Barcodes(File_k_TranscriptNumber, InitFirefly, ThisBarcode); // 为了得到所有样本的结果; 
                // std::cout << "(=^_^=) Group " << groupinform.GroupIndex << " reads numbers " << groupinform.GroupReadSjs.size()  << " completed quantification! (=^_^=) " << std::endl;               
            }
            int fileIndex = 0;
            write_transcript_file(File_k_TranscriptNumber, Finally_Annotations, AllFile_BarcodeSet["0"], fileIndex, IsoformFilePath, IsoformMutexes); 
        }
    }
}

// 一个简单的处理函数，可以替换为你的处理逻辑
void processGroup(std::streampos& start, std::streampos& end, 
                  const std::string& sam_file_path, 
                  const int& Sj_supportReadNumber, GTF& gtf_full, 
                  GTFsj& gtf_splice, const int& GraphDis,
                  std::ofstream& updatedgtffile, std::ofstream& tracefile,
                  std::vector<std::unique_ptr<std::ofstream>>& isoformfilePath,
                  std::vector<std::mutex>& isoformMutexes,
                  int& fileno, std::string& outputPath,
                  std::map<std::string, std::set<std::string>>& filebarcodeSet,
                  int& singleEdge) {
    // 这里是任务处理逻辑，读取从 start 到 end 行的数据;
    // 1.提取这个Group中有用的信息;
    GroupInformation group_information = knowGroupInformation(start, end, sam_file_path, Sj_supportReadNumber);
    std::string chrchr = group_information.chrName;
    std::array<int,2> groupcoverage = group_information.GroupCoverage;
    std::unordered_map<std::string, std::array<int,2>> single_exon_group_annotation;
    // std::cout << "[{(>_<)]} Group sj reads numbers are " << group_information.GroupReadSjs.size() << " ! ^-^ Group single exon reads numbers are " << group_information.GroupSingleExon.size() << " ! [{(>_<)]}" << std::endl;
    
    // single exon的一些结果;
    if (group_information.GroupSingleExon.size() > 0) {
        // 提取这个Group内的所有单个exon的注释;
        single_exon_group_annotation = get_group_single_exon_annotation(gtf_splice.SE[chrchr], groupcoverage);
        // reads分选进去;
        std::unordered_map<std::string, std::vector<std::string>> single_exon_with_reads = get_group_single_exon_reads(single_exon_group_annotation, group_information.GroupSingleExon, singleEdge);
        // 一个文件或多个文件都写; <file:<gene:<barcode, counts>>>
        write_single_exon_gtf_trace_sc(fileno, group_information.GroupReadFiles,
                                    single_exon_with_reads, single_exon_group_annotation, 
                                    updatedgtffile, tracefile, isoformfilePath, isoformMutexes,
                                    gtf_full.GTF_transcript_strand[chrchr],
                                    group_information.GroupReadBarcodes, chrchr,
                                    filebarcodeSet);
    }
    std::cout << group_information.GroupIndex << " is end! Single reads size is " << group_information.GroupSingleExon.size() << "!"<< std::endl;
    group_information.GroupSingleExon.clear();
    // 整个cluster小于5条的直接结束;
    if (group_information.GroupReadSjs.size() > 0) {
        // 提取这个Group内的所有注释;
        GroupAnnotation group_annotation = get_group_annotation(gtf_splice.mSJsBE[chrchr], gtf_splice.mSJs[chrchr], groupcoverage);
        // 得到所有的Cluster, 然后分成了FSM\ISM\Others;
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
                    single_exon_group_annotation); 
    }
    std::cout << group_information.GroupIndex << " is end! reads size is " << group_information.GroupReadSjs.size() << "!"<< std::endl;
}

// void Remove_cache_file () {
    
// }


// 主函数;
int main(int argc, char* argv[])
{   
    // 命名空间; 
    int option_index = 0;
    int c;
    //SAM文件地址;
    //GTF注释文件地址, GENCODE下载的完整基因注释ALL;
    //Fasta文件, GENCODE下载的Genome sequence (GRCh38.p14);
    std::string samfile_name;
    std::string fastafile_name;
    std::string gtffile_name;
    std::string output_file_name;
    //定义两个M的小区间之间相差距离, 定义为SJ;
    int SJDistance = 18;
    int SJ_support_read_number = 2;
    int Graph_distance = 60;
    int Thread = 8;
    int single_exon_edge = 60;

    // 使用 getopt_long 解析命令行参数
    while ((c = getopt_long(argc, argv, "s:f:g:o:m:j:n:e:d:t:h", long_options, &option_index)) != -1) {
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
                SJDistance = std::stoi(optarg);  // 定义两个M的小区间之间相差距离, 定义为SJ;
                break;
            case 'n':
                SJ_support_read_number = std::stoi(optarg);   // 支持SJ的reads的个数;
                break;
            case 'e':
                single_exon_edge = std::stoi(optarg);   // single exon的边界;
                break;
            case 'd':
                Graph_distance = std::stoi(optarg);   // 计算子图的距离;
                break;
            case 't':
                Thread = std::stoi(optarg);   // 线程;
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

    // 检查是否所有的必要文件路径都已提供, 必须输入的必须要有提供;
    if (samfile_name.empty() || fastafile_name.empty() || output_file_name.empty()) {
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    // 检查并创建输出文件; 0.updated_gtf; 1.counts_transcript; 2.counts_gene; 3.trace;
    std::vector<std::string> sam_file_vec = traverse_sam_file(samfile_name, output_file_name);
    std::vector<std::string> outputFileVec = check_catalog_exist(output_file_name);
    std::ofstream gtf_file(outputFileVec[0], std::ios::app);
    std::ofstream trace_file(outputFileVec[1], std::ios::app);

    // 输出提供的文件路径和参数;
    std::cout << "*****" << std::endl;
    std::cout << "SAM file: " << samfile_name << std::endl;
    std::cout << "FASTA file: " << fastafile_name << std::endl;
    std::cout << "GTF file: " << gtffile_name << std::endl;
    std::cout << "Output file: " << output_file_name << std::endl;
    std::cout << "SJDistance: " << SJDistance << std::endl;
    std::cout << "SJ_support_read_number: " << SJ_support_read_number << std::endl;
    std::cout << "single_exon_edge : " << single_exon_edge << std::endl;
    std::cout << "Graph_distance: " << Graph_distance << std::endl;
    std::cout << "Thread: " << Thread << std::endl;
    std::cout << "*****" << std::endl;

    //Fasta文件的结果; 读取fasta文件;
    std::unordered_map<std::string, std::string> Fasta = Read_fasta_file(fastafile_name);

    // 1. 对所有sam文件进行处理; 每个sample生成一个文件夹; 每个文件夹里面有很多小文件;
    // 2. 合并所有的小文件;
    FileSplit BroCOLIfile = thread_all_read_sam_files(samfile_name, sam_file_vec, Thread, output_file_name, SJDistance, Fasta);
    Fasta.clear();
    std::vector<std::unique_ptr<std::ofstream>> BroCOLIQuantfile = write_quantification_files(output_file_name, BroCOLIfile.file_barcodeSet);
    std::vector<std::mutex> TranscriptMutexes(BroCOLIQuantfile.size());

    // 读取注释文件; 注释信息;
    GTF GTF_full = get_gtf_annotation(gtffile_name);
    // 注释文件中的SJ和single exon的字典map; 根据提取的注释文件信息, 获取SJ信息和单个exon的信息;
    GTFsj GTF_Splice = get_SJs_SE(GTF_full.GTF_transcript);

    std::vector<std::size_t> Group_idx = sort_indexes_e(BroCOLIfile.group_reads_number);
  
    std::cout << "BroCOLIfile.group_reads_number: " << BroCOLIfile.group_reads_number.size() << std::endl;
    std::cout << "BroCOLIfile.reads_pointer " << BroCOLIfile.reads_pointer.size() << std::endl;
    BroCOLIfile.chr_coverage.clear();
    BroCOLIfile.coverage2pos.clear();
    BroCOLIfile.group_reads_number.clear();

    // 删除所有的临时文件的文件夹;

    // 创建一个线程池, 包含Thread个线程;
    ThreadPool BroCOLIpool(Thread);
    // 创建并分配任务到线程池;
    std::vector<std::future<void>> futures;
    for (const auto& i:Group_idx) {
        futures.emplace_back(BroCOLIpool.enqueue([&, i]() { // 捕获所有外部变量，并传递 i
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
                single_exon_edge);
            }));
    }
    // 等待所有任务完成;
    for (auto& future : futures) {
        future.get();  // 等待每个任务完成
    }
    // 关闭文件;
    gtf_file.close();
    trace_file.close();
    return 0;
}




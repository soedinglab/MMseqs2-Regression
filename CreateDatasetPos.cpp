#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <algorithm>    // copy
#include <iterator>     // back_inserter
#include <vector>
#include <iostream>
#include <unordered_map>
#include <cstring>

#include "PatternCompiler.h"

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

struct ScopEntry{
    char query[33];
    char scop[33];
    const float evalue;
    const int query_start;
    const int query_end;
    const int template_start;
    const int template_end;
    bool hasInsertion;
    ScopEntry(const char * q,const  char * s,
              float eval, int qstart,
              int qend, int tstart, int tend) : evalue(eval),
                                                query_start(qstart), query_end(qend),
                                                template_start(tstart), template_end(tend), hasInsertion(false)
    {
        int qStrLen = strlen(q);
        int sStrLen = strlen(s);
        memcpy((char *)query, q, qStrLen * sizeof(char));
        memcpy((char *)scop, s, sStrLen * sizeof(char));
        query[qStrLen] = '\0';
        scop[sStrLen] = '\0';
    };
};


std::vector <std::string> split(const std::string& str, const std::string& delimiter = " ") {
    std::vector <std::string> tokens;

    std::string::size_type lastPos = 0;
    std::string::size_type pos = str.find(delimiter, lastPos);

    while (std::string::npos != pos) {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = pos + delimiter.size();
        pos = str.find(delimiter, lastPos);
    }

    tokens.push_back(str.substr(lastPos, str.size() - lastPos));
    return tokens;
}


struct sortDescByTargetStart {
    bool operator()(const ScopEntry * left, const ScopEntry * right) {
        return left->template_start < right->template_start ||
                (left->template_start == right->template_start && left->template_end > right->template_end)
                && left->evalue < right->evalue;
    }
};


int overlap(int min1, int max1, int min2, int max2) {
    return std::max(0, std::min(max1, max2) - std::max(min1, min2));
}


int main(int argc, char ** argv) {
    // open csv file
    std::string csvFilePath(argv[1]);
    std::cout << "Read file " << argv[1] << std::endl;
    std::ifstream csvFile(csvFilePath);
    // open ffindex
    std::cout <<"Read ffindex" <<std::endl;
    std::string dataFilePath(argv[2]);
    std::string indexFilePath = dataFilePath + ".index";
    size_t dataSize;
    FILE * dataFile = fopen(dataFilePath.c_str(), "r");
    char * data = ffindex_mmap_data(dataFile, &dataSize);
    FILE * indexFile = fopen(indexFilePath.c_str(), "r");
    std::ifstream index_file(indexFilePath);
    size_t cnt = 0;
    if (index_file.is_open()) {
        std::string line;
        while ( std::getline(index_file, line)){
            cnt++;
        }
        index_file.close();
    }
    std::cout <<"Parse ffindex" <<std::endl;
    ffindex_index_t * index = ffindex_index_parse(indexFile, cnt);

    // open output stream
    std::ofstream fastaOutput;
    fastaOutput.open (argv[3]);

    std::string line;
    std::unordered_map<std::string, ScopEntry*> mappingDict;
    PatternCompiler digit("[0-9]+");
    //std::getline(csvFile, line); // get header;
    size_t uniqIds = 0;
    size_t allEntries = 0;
    size_t validDomains = 0;
    size_t insertionCount = 0;
    double EVAL_THRESHOLD = 0.00001;
    std::cout <<"Parse csv file" <<std::endl;
    while (std::getline(csvFile, line))
    {
        allEntries++;
//        Query Scop Hit Prob E-value P-value Score SS Cols QueryPos TemplatePos
//        d12asa_ d.104.1.1 sp|Q3Z8X9|SYK_DEHM1 100.0 5.8E-74 2.5E-79 542.4 0.0 275 2-316 171-486 (498)
        std::string query, scop, hit, queryPos, templatePos;
        double prob, eval, pval, score, ss;
        int cols;
        std::istringstream is(line);
        is >> query >> scop >>  hit >> prob >> eval >> pval >> score >> ss >> cols >> queryPos >> templatePos;
        // check if exists
        std::vector<std::string> hitSplit = split(std::string(hit),"|");
        std::string key;
        if(hitSplit.size() > 1){
            key = hitSplit[1];
        }else{
            key = hitSplit[0];
        }
        if(eval < EVAL_THRESHOLD){

            std::vector<std::string> qTokens = digit.getAllMatches(queryPos.c_str(), 255);
            std::vector<std::string> tTokens = digit.getAllMatches(templatePos.c_str(), 255);
            int qStart = std::stoi(qTokens[0]);
            int qEnd = std::stoi(qTokens[1]);
            int tStart = std::stoi(tTokens[0]);
            int tEnd = std::stoi(tTokens[1]);
            int gaps = abs((qStart - tStart) - (qEnd - tEnd));
            //std::cout << gaps << std::endl;
            // not too much insertions or deletions
            if(gaps >= 50 ){
                continue;
            }

            ScopEntry * scopEntry = new ScopEntry(query.c_str(), scop.c_str(), eval,
                                                 qStart, qEnd,
                                                 tStart, tEnd);
            validDomains++;

            if ( mappingDict.find(key) == mappingDict.end() ) {
                mappingDict[key] = scopEntry;
                //mappingDict[hit].reserve(2);
                uniqIds++;
            }else{
                // THE INPUT FILE SHOULD BE ORDERED BY EVAL (SO NO NEED)
                if(mappingDict[key]->evalue > scopEntry->evalue ){
                    ScopEntry * oldScop = mappingDict[key];
                    delete oldScop;
                    mappingDict[key] = scopEntry;
                }
                int overlappingDomain = overlap(mappingDict[key]->template_start, mappingDict[key]->template_end,
                        scopEntry->template_start, scopEntry->template_end);

                std::string currFold = mappingDict[key]->scop;
                currFold = currFold.erase(currFold.find_last_of("."), std::string::npos);
                currFold = currFold.erase(currFold.find_last_of("."), std::string::npos);
                std::string newAnnotationFold = scopEntry->scop;
                newAnnotationFold = newAnnotationFold.erase(newAnnotationFold.find_last_of("."), std::string::npos);
                newAnnotationFold = newAnnotationFold.erase(newAnnotationFold.find_last_of("."), std::string::npos);
//                if(currFold.compare(newAnnotationFold) != 0){
//                    std::cout << overlappingDomain << "\t" << newAnnotationFold << "\t" << currFold << std::endl;
//                }
                // if both annotations overlap but they have different folds
                if(overlappingDomain > 5 && currFold.compare(newAnnotationFold) != 0)
                {
                    mappingDict[key]->hasInsertion=true;
                    insertionCount++;
                }
            }
        }
        if((allEntries % 10000) == 0){
            std::cout << ".";
            std::cout << std::flush;
        }

        if((allEntries % 1000000) == 0){
            std::cout << std::endl;
        }
    }
    std::cout << uniqIds << " uniq sequences of " << validDomains << " domains out of " << allEntries << " will be considered" << std::endl;
    std::cout <<"Found " << insertionCount << " insertions" << std::endl;

    std::unordered_map<std::string, ScopEntry*>::iterator it;
    validDomains = 0;
    std::cout << "Remove covered domains by same family" << std::endl;
    std::vector<std::string> idsToDelete;



    std::cout << "Write file to " << argv[3] << " will be considered" << std::endl;
    for ( it = mappingDict.begin(); it != mappingDict.end(); it++ )
    {
        std::string hit = it->first;
        ffindex_entry_t * entry = ffindex_get_entry_by_name(index, (char *)hit.c_str());
        if(entry->length > 20000){ //hhblits can not align sequences >20000 (titin :( )
            std::cout << "Removed " << hit << " because size is > 20000" << std::endl;
            continue;
        }
        ScopEntry* domain =  it->second;
        if(domain->hasInsertion){
            continue;
        }

        fastaOutput << ">";
        fastaOutput << hit;
        fastaOutput << " ";
        fastaOutput << domain->scop;
        fastaOutput << " ";
        fastaOutput << "|";
        fastaOutput << domain->evalue;
        fastaOutput << " ";
        fastaOutput << "|";
        fastaOutput << domain->query;
        fastaOutput << " ";
        fastaOutput << "|";
        fastaOutput << domain->query_start;
        fastaOutput << "-";
        fastaOutput << domain->query_end;
        fastaOutput << " ";
        fastaOutput << "|";
        fastaOutput << domain->template_start;
        fastaOutput << "-";
        fastaOutput << domain->template_end;
        fastaOutput << " ";
        delete domain;
        fastaOutput <<"\n";
        fastaOutput << ffindex_get_data_by_entry(data, entry);
    }
    fastaOutput.close();
    csvFile.close();
}
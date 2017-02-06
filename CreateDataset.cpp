#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <algorithm>    // copy
#include <iterator>     // back_inserter
#include <vector>
#include <iostream>
#include <unordered_map>

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
    ScopEntry(const char * q,const  char * s,
              float eval, int qstart,
              int qend, int tstart, int tend) : evalue(eval),
                                                query_start(qstart), query_end(qend),
                                                template_start(tstart), template_end(tend)
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
    std::unordered_map<std::string, std::vector<ScopEntry*>> mappingDict;
    PatternCompiler digit("[0-9]+");
    //std::getline(csvFile, line); // get header;
    size_t uniqIds = 0;
    size_t allEntries = 0;
    size_t validDomains = 0;
    double EVAL_THRESHOLD = 0.00001;
    std::cout <<"Parse csv file" <<std::endl;
    while (std::getline(csvFile, line))
    {
        allEntries++;
//        Query Scop Hit Prob E-value P-value Score SS Cols QueryPos TemplatePos
//        d12asa_ d.104.1.1 sp|Q3Z8X9|SYK_DEHM1 100.0 5.8E-74 2.5E-79 542.4 0.0 275 2-316 171-486 (498)
        char query[255], scop[255], hit[255], queryPos[255], templatePos[255];
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
            validDomains++;
            if ( mappingDict.find(key) == mappingDict.end() ) {
                mappingDict[key] = std::vector<ScopEntry*>();
                //mappingDict[hit].reserve(2);
                uniqIds++;
            }

            std::vector<std::string> qTokens = digit.getAllMatches(queryPos, 255);
            std::vector<std::string> tTokens = digit.getAllMatches(queryPos, 255);
            mappingDict[key].push_back(new ScopEntry((const char *)query,(const char *) scop, eval,
                                                     std::stoi(qTokens[0]), std::stoi(qTokens[1]),
                                                     std::stoi(tTokens[0]), std::stoi(tTokens[1])));
        }
    }
    std::cout << uniqIds << " uniq sequences of " << validDomains << " domains out of " << allEntries << " will be considered" << std::endl;
        std::unordered_map<std::string, std::vector<ScopEntry*>>::iterator it;
    validDomains = 0;
    std::cout << "Remove covered domains by same family" << std::endl;
    std::vector<std::string> idsToDelete;

    for ( it = mappingDict.begin(); it != mappingDict.end(); it++ ){
        std::vector<ScopEntry*> domains =  it->second;
        std::vector<ScopEntry*> newDomains;
        // sort by x,y like
        // 1,4
        // 1,3
        // 2,3
        // 2,4
        // 2,5

        std::sort(domains.begin(), domains.end(), sortDescByTargetStart());

        newDomains.push_back(domains[0]);
        ScopEntry * coverDomain = domains[0];
        std::string covFam = coverDomain->scop;
        //covSupFam = covSupFam.erase(covSupFam.find_last_of("."), std::string::npos);
        validDomains++;
//        bool hasEvalueLess1E_3=false;
//        for(size_t i = 0; i < domains.size(); i++) {
//            // check if TP is possible
//            if(domains[i]->evalue < EVAL_THRESHOLD){
//                hasEvalueLess1E_3 = true;
//            }
//        }
//        if(hasEvalueLess1E_3==false){
//            idsToDelete.push_back(it->first);
//            continue;
//        }
        for(size_t i = 1; i < domains.size(); i++){
            // domains[i-1] fully covered domains[i]
            // new domain
            std::string domainFam = domains[i]->scop;
            //domainSupFam = domainSupFam.erase(domainSupFam.find_last_of("."), std::string::npos);
            // if not overlaps or same familiy add domain
            if(coverDomain->template_end < domains[i]->template_start || covFam.compare(domainFam)!=0){
                coverDomain = domains[i];
                covFam = domainFam;
                newDomains.push_back(domains[i]);
                validDomains++;
            }

            // domains[i-1] covers domains[i] left side covered

            // right side covered
        }
        domains.clear();
        it->second = newDomains;
    }
    //remove
//    for(size_t i = 0; i < idsToDelete.size(); i++){
//        mappingDict.erase(idsToDelete[i]);
//    }
    std::cout << uniqIds -idsToDelete.size() << " uniq sequences of " << validDomains << " domains out of " << allEntries << " will be considered" << std::endl;
    std::cout << "Write file to " << argv[3] << " will be considered" << std::endl;
    for ( it = mappingDict.begin(); it != mappingDict.end(); it++ )
    {
        std::string hit = it->first;
        ffindex_entry_t * entry = ffindex_get_entry_by_name(index, (char *)hit.c_str());
        if(entry->length > 20000){ //hhblits can not align sequences >20000 (titin :( )
            std::cout << "Removed " << hit << " because size is > 20000" << std::endl;
            continue;
        }
        std::vector<ScopEntry*> domains =  it->second;
        fastaOutput << ">";
        fastaOutput << hit;
        fastaOutput << " ";
        for(size_t i = 0; i < domains.size(); i++) {
            fastaOutput << domains[i]->scop;
            fastaOutput << " ";
        }
        fastaOutput << "|";
        for(size_t i = 0; i < domains.size(); i++) {
            fastaOutput << domains[i]->evalue;
            fastaOutput << " ";
        }
        fastaOutput << "|";
        for(size_t i = 0; i < domains.size(); i++) {
            fastaOutput << domains[i]->query;
            fastaOutput << " ";
        }
        fastaOutput << "|";
        for(size_t i = 0; i < domains.size(); i++) {
            fastaOutput << domains[i]->query_start;
            fastaOutput << "-";
            fastaOutput << domains[i]->query_end;
            fastaOutput << " ";
        }
        fastaOutput << "|";
        for(size_t i = 0; i < domains.size(); i++) {
            fastaOutput << domains[i]->template_start;
            fastaOutput << "-";
            fastaOutput << domains[i]->template_end;
            fastaOutput << " ";
        }
        for(size_t i = 0; i < domains.size(); i++) {
            delete domains[i];
        }
        fastaOutput <<"\n";
        fastaOutput << ffindex_get_data_by_entry(data, entry);
    }
    fastaOutput.close();
    csvFile.close();
}
//
// Created by mad on 9/22/15.
//
#include "EvaluateResults.h"
#include <stddef.h>
#include <iostream>
#include <fstream>
#include <unistd.h>

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}
#include <regex>        // regex, sregex_token_iterator
#include "EvaluateResults.h"
#include "kseq.h"
#define MAX_FILENAME_LIST_FILES 4096
KSEQ_INIT(int, read)

double TP_EVAL_THRESHOLD = 0.001;

int main(int argc, char ** argv){
    size_t rocx = 5;
    size_t res_cnt = 0;
    //while( true)
    std::string queryFasta=argv[1];
    std::string targetFasta=argv[2];
    std::string resultFile=argv[3];
    std::string outputResultFile=argv[4];
    std::unordered_map<std::string, std::vector<SCOP>*> scopLoopup;
    std::cout << "Read query fasta" << std::endl;
    std::unordered_map<std::string, size_t > whatever;
    readFamDefFromFasta(queryFasta, scopLoopup,
                        whatever, false);
    std::cout << "Read target fasta " << std::endl;
    std::unordered_map<std::string, size_t> scopSizeLoopup;
    //std::cout << scopLoopup["d12asa_"]->at(0) << " " << famSizeLoopup[scopLoopup["d1acfa_"]->at(0)] <<std::endl;
    readFamDefFromFasta(targetFasta, scopLoopup,
                        scopSizeLoopup, true);
    std::cout << std::endl;
    std::cout << "Read result fasta " << std::endl;

    FILE * fasta_file = fopen(queryFasta.c_str(), "r");
    kseq_t *seq = kseq_init(fileno(fasta_file));
    size_t entries_num = 0;
    double overall_ignore_cnt =0.0;
    double overall_fp =0.0;
    double overall_tp =0.0;
    double early_break_cnt=0.0;
    std::vector<Roc5Value> roc5Vals;
    std::vector<Hits> allHits;
    // iterate over all queries
    while (kseq_read(seq) >= 0) {
        std::string query = seq->name.s;
//        std::cout << query << std::endl;
        std::vector<SCOP> * qFams = scopLoopup[query];
        std::vector<std::pair<std::string, double>> resIds = readResultFile(query, resultFile);
        EvaluateResult eval = evaluateResult(query, qFams, scopLoopup,
                                             allHits, resIds, rocx);
//            if(query.compare("d2py5a2") == 0){
//                for(size_t j = 0; j < resIds.size(); j++){
//                    std::cout << resIds[j] << " " << scopLoopup[resIds[j]]->at(0) << std::endl;
//                }
//            }
        overall_ignore_cnt +=eval.ignore_cnt;
        overall_fp += eval.fp_cnt;
        overall_tp += eval.tp_cnt;

        if(eval.fp_cnt < rocx)
            early_break_cnt++;
        double all_auc = eval.auc + (rocx - eval.fp_cnt) * eval.tp_cnt;

        double qFamSize = 0.0;
        std::string qFamStr;
        for(size_t i = 0; i < qFams->size(); i++) {
            SCOP qFam = qFams->at(i);
            qFamSize = std::min(qFamSize + scopSizeLoopup[qFam.superFam], 4000.0);
            qFamStr.append(qFams->at(i).fam).append(", ");
        }
        if(qFamSize > 0){
            double roc5val = all_auc / (rocx * qFamSize);
            if (roc5val > 1.0){
                std::cout << "ROC5 = " << roc5val << " for query " << query << ", # family members: " <<
                qFamSize << std::endl;
                std::cout << "Results size: " << resIds.size() << std::endl;
                std::cout << "TPs: " << eval.tp_cnt << ", FPs: " << eval.fp_cnt  << std::endl;
            }else{
                roc5Vals.push_back(Roc5Value(query, qFamStr, qFamSize, roc5val, eval.tp_cnt, eval.fp_cnt, eval.ignore_cnt, resIds.size()));
            }
        }else {
            std::cout << "Fam = " << qFamStr << " for query " << query << ", # family members: " << qFamSize << std::endl;
        }
    }

    std::sort(roc5Vals.begin(), roc5Vals.end(), sortDescByRoc5());
    printf( "Query\t\tFam\t\t\t\tRoc5\tFamSize\t(TPs\tFP\tResSize\tIGN)\n");
    for(size_t i = 0; i < roc5Vals.size(); i++) {
        Roc5Value roc5Value = roc5Vals[i];
        printf("%s\t\t%-30.30s\t%.7f\t%5d\t(TP: %5d FP: %5d RES: %5d IGN: %5d)\n", roc5Value.query.c_str(), roc5Value.qFams.c_str(),
               roc5Value.roc5val, roc5Value.qFamSize, roc5Value.tp_cnt, roc5Value.fp_cnt,
               roc5Value.resultSize, roc5Value.ignore_cnt);
    }
    std::sort(allHits.begin(), allHits.end(), sortFalsePositvesByEval());
    double fpsWithSmallEvalue=0;
    double EVAL_THRESHOLD = 1E-3;
    std::map<std::string, size_t > mostQueriesWithSmallEval;
    for(size_t i = 0; i < allHits.size(); i++) {
        if(allHits[i].evalue < EVAL_THRESHOLD && allHits[i].status == Hits::FP){
            fpsWithSmallEvalue++;
            mostQueriesWithSmallEval[allHits[i].query]++;
        }
    }
    std::cout << "Top high scoring queries:" << std::endl;
    std::vector<std::pair<size_t, std::string>> mostQueriesWithSmallEvalVec;
    for (std::map<std::string, size_t >::iterator it = mostQueriesWithSmallEval.begin();
         it != mostQueriesWithSmallEval.end(); it++ ) {
        mostQueriesWithSmallEvalVec.push_back(std::make_pair(it->second, it->first));
    }
    std::sort(mostQueriesWithSmallEvalVec.begin(), mostQueriesWithSmallEvalVec.end());
    for (int i = mostQueriesWithSmallEvalVec.size(); i > 0; i--) {
        std::cout << mostQueriesWithSmallEvalVec[i-1].second << " " << mostQueriesWithSmallEvalVec[i-1].first << std::endl;
    }

    std::cout << "Top 200 FP:" << std::endl;
    size_t cnt=0;
    for(size_t i = 0; i < allHits.size(); i++) {
        if(allHits[i].status == Hits::FP){
            std::cout << cnt + 1 << ": " << allHits[i].query << " " << allHits[i].target << " " << allHits[i].evalue << std::endl;
            std::vector<SCOP> * scopTarget =  scopLoopup[allHits[i].target];
            std::vector<SCOP> * scopQuery =  scopLoopup[allHits[i].query];
            std::cout << "Q=";
            for(size_t j = 0; j < scopQuery->size(); j++){
                std::cout << " " << scopQuery->at(j).fam << "(" << scopQuery->at(j).evalue << ")";
            }
            std::cout << std::endl;
            std::cout << "T=";
            for(size_t j = 0; j < scopTarget->size(); j++){
                std::cout << " " << scopTarget->at(j).fam << "(" << scopTarget->at(j).evalue << ")";
            }
            std::cout << std::endl;
            std::cout << std::endl;

            cnt++;
            if(cnt==200)
                break;
        }
    }
    kseq_destroy(seq);
    std::cout << res_cnt  << " result lists checked." << std::endl;
    std::cout << early_break_cnt << " result lists did not contain 5 FPs." << std::endl;
    std::cout << "Results contains " << overall_tp << " TPs and " << overall_fp << " FPs." << std::endl;
    std::cout << "Total FPs " << fpsWithSmallEvalue << " of " << mostQueriesWithSmallEvalVec.size() << " queries have an eval < " << EVAL_THRESHOLD << "." << std::endl;
    std::cout << overall_ignore_cnt << " sequence pairs ignored (different family, same fold)" << std::endl;

    writeRoc5Data(outputResultFile, roc5Vals, 0.01);
    writeRocData(outputResultFile, allHits, 10000);
    return 0;
}

void parseMMseqs(std::string query, std::string resFileName, std::vector<std::pair<std::string,double>> & resultVector) {
    static bool isReadIn = false;
    static char * data = NULL;
    static ffindex_index_t * index = NULL;
    if(isReadIn == false){
        int lastindex = resFileName.find_last_of(".");
        std::string dataFile = resFileName.substr(0, lastindex);
        FILE * dataFP = fopen(dataFile.c_str(), "r");
        size_t dataSize;
        data = ffindex_mmap_data(dataFP, &dataSize);
        FILE * indexFile = fopen(resFileName.c_str(), "r");
        std::ifstream index_file(resFileName);
        size_t cnt = 0;
        if (index_file.is_open()) {
            std::string line;
            while ( std::getline(index_file, line)){
                cnt++;
            }
            index_file.close();
        }
        std::cout <<"Parse ffindex" <<std::endl;
        index = ffindex_index_parse(indexFile, cnt);
        isReadIn = true;
    }
    if(isReadIn == false){
        std::cout << "Readin of mmseqs results did not work" << std::endl;
    }
    std::regex keyRegex("\\S+");
    ffindex_entry_t * entry = ffindex_bsearch_get_entry(index, (char *) query.c_str());
    char * result = ffindex_get_data_by_entry(data, entry);
    std::cregex_token_iterator tBegin(result, result + entry->length, keyRegex), tEnd;
    std::vector<std::string> tmpRes;
    std::copy(tBegin, tEnd, std::back_inserter(tmpRes));
    for(size_t i = 0; i < tmpRes.size(); i+=6){
        if(tmpRes[i].c_str()[0] == '\0')
            break;
        resultVector.push_back(std::make_pair(tmpRes[i], atof(tmpRes[i+5].c_str())));
    }
}

void parseM8(std::string query, std::string resFileName, std::vector<std::pair<std::string, double>> &resultVector) {
    static bool isReadIn = false;
    static std::map<std::string, std::vector<std::pair<std::string, double>>> resLookup;
    if(isReadIn == false){
        std::cout << "Read in m8 " << resFileName << std::endl;
        std::ifstream infile(resFileName);
        std::string line;
        std::regex keyRegex("\\S+");

        while (std::getline(infile, line))
        {
            std::cregex_token_iterator tBegin(line.c_str(), line.c_str() + line.size(), keyRegex), tEnd;
            std::string key = *tBegin;
            if(resLookup.find(key)== resLookup.end()) {
                resLookup[key] = std::vector<std::pair<std::string, double>>();
            }
            tBegin++;
            std::string targetkey = *tBegin;
            for(size_t i = 0; i < 9; i++){
                tBegin++;
            }
            std::string evalStr = *tBegin;
            double eval = atof(evalStr.c_str());

            resLookup[key].push_back(std::make_pair(targetkey,eval));
        }
        infile.close();
        std::map<std::string, std::vector<std::pair<std::string, double>>>::iterator it;

        for ( it = resLookup.begin(); it != resLookup.end(); it++ ) {
            std::set<std::string> removeDub;
            std::vector<std::pair<std::string, double>> single;

            for(size_t i = 0; i < it->second.size(); i++){
                if(removeDub.find(it->second[i].first) == removeDub.end()){
                    single.push_back(it->second[i]);
                }
                removeDub.insert(it->second[i].first);
            }
            it->second.clear();
            removeDub.clear();
            for(size_t i = 0; i < single.size(); i++){
                it->second.push_back(single[i]);
            }
            single.clear();
        }
        isReadIn = true;

    }
    resultVector.swap(resLookup[query]);
}

std::vector<std::pair<std::string, double>> readResultFile(std::string query, std::string resFileName) {
    std::vector<std::pair<std::string, double>> resultVector;
    std::string extention = resFileName.substr(resFileName.find_last_of(".") + 1);
    if(extention.compare("index") == 0) { // MMSeqs
        parseMMseqs(query, resFileName, resultVector);
    } else if (extention.compare("m8") == 0){
        parseM8(query, resFileName, resultVector);
    }
    return resultVector;
}


void printProgress(int id){
    if (id % 1000000 == 0 && id > 0){
        std::cout << "\t" << (id/1000000) << " Mio. sequences processed\n";
        fflush(stdout);
    }
    else if (id % 10000 == 0 && id > 0) {
        std::cout  << ".";
        fflush(stdout);
    }
}

void readFamDefFromFasta(std::string fasta_path, std::unordered_map<std::string, std::vector<SCOP> *> &queryScopLookup,
                         std::unordered_map<std::string, size_t > &supFamSizeLookup, bool readEval) {
    FILE * fasta_file = fopen(fasta_path.c_str(), "r");
    kseq_t *seq = kseq_init(fileno(fasta_file));
    size_t entries_num = 0;

    std::regex scopDomainRegex("\\S+\\.\\d+\\.\\d+\\.\\d+");
    std::set<std::string> scopDomains;

    while (kseq_read(seq) >= 0) {
        if (seq->name.l == 0) {
            std::cout << "Fasta entry: " << entries_num << " is invalid." << std::endl;
            exit;
        }
        const std::string currQuery(seq->name.s);
        if(queryScopLookup.find(currQuery)== queryScopLookup.end()) {
            queryScopLookup[currQuery] = new std::vector<SCOP>();
        }
        std::vector<SCOP> * queryDomainVector = queryScopLookup[currQuery];
        std::vector<std::string> splits = split(std::string(seq->comment.s), "|");
        std::vector<std::string> evals;
        if(readEval== true){
            evals = split(splits[1]," ");
        }
        std::cregex_token_iterator begin(splits[0].c_str(), splits[0].c_str()+splits[0].size(),
                                         scopDomainRegex), end;
        scopDomains.clear();
        while(begin != end){
            scopDomains.insert(*begin);
            ++begin;
        }
        int i = 0;
        for(std::set<std::string>::iterator it = scopDomains.begin(); it != scopDomains.end(); it++) {
            std::string currScopDomain = *it;
            double eval = (readEval == true) ? strtod(evals[i].c_str(), NULL) : 0.0;
            // increase the scop domain count
            SCOP domain = SCOP(currScopDomain, eval);
            if (supFamSizeLookup.find(domain.fam) == supFamSizeLookup.end()) {
                supFamSizeLookup[domain.fam] = 0;
            }
            // count only possible TP
            if (eval < TP_EVAL_THRESHOLD){
                supFamSizeLookup[domain.fam]++;
                supFamSizeLookup[domain.superFam]++;
                supFamSizeLookup[domain.fold]++;
            }
            queryDomainVector->push_back(domain);
            i++;
        }
        entries_num++;
        printProgress(entries_num);

    }
    kseq_destroy(seq);

}

EvaluateResult evaluateResult(std::string query, std::vector<SCOP> *qScopIds, std::unordered_map<std::string,
        std::vector<SCOP> *> &scopLoopup, std::vector<Hits> &fpVec,
                              std::vector<std::pair<std::string, double>> results, size_t rocx) {
    double fp_cnt = 0.0;
    double tp_cnt = 0.0;
    double ignore_cnt = 0.0;
    double auc = 0.0;

//    std::string qSupFam = qFam;
//    qSupFam = qSupFam.erase(qSupFam.find_last_of("."), std::string::npos);
    for (size_t i = 0; i < results.size(); i++) {
        const std::string rKey = results[i].first;
        const double evalue = results[i].second;

        bool tp = false;
        bool fp = false;
        bool ignore = false;
        std::vector<SCOP> * rfamVec;
        // if sequence does not have annotations ignore it
        if (scopLoopup.find(rKey) == scopLoopup.end()) {
            tp = false;
            ignore = true;
            fp = false;
            goto outer;
        }
        rfamVec = scopLoopup[rKey];

        for(size_t j = 0; j < rfamVec->size(); j++) {
            for(size_t i = 0; i < qScopIds->size(); i++) {
                SCOP qScopId = qScopIds->at(i);
                const SCOP rScopId = rfamVec->at(j);
                if (rScopId.superFam.compare(qScopId.superFam) == 0 && rScopId.evalue < TP_EVAL_THRESHOLD) {
                    tp = true;
                    goto outer;
                }
                if (tp == false) {
                    std::regex ignore_superfam("^b\\.(67|68|69|70).*");
                    bool qSuperFamIgnore = std::regex_match(qScopId.fam, ignore_superfam);
                    bool rSuperFamIgnore = std::regex_match(rScopId.fam, ignore_superfam);


                    std::regex ignoreClass("^e\\..*");
                    bool qFoldIgnore = std::regex_match(qScopId.fam, ignoreClass);
                    bool rFoldIgnore = std::regex_match(rScopId.fam, ignoreClass);
                    if ((rScopId.fold.compare(qScopId.fold) == 0)
                        || (qSuperFamIgnore && rSuperFamIgnore)
                        || qFoldIgnore
                        || rFoldIgnore) {
                        ignore = true;
//                        goto outer;
                    } else {
                        fp = true;
                    }
                }
            }
        }
        outer:

//        if(tp){
//            std::cout << rKey << "\n";
//        }
        // counter for ROC5 values

        if (fp_cnt < rocx)
        {
            if (tp == true)
            {
                tp_cnt++;
                fpVec.push_back(Hits(query,rKey, evalue, Hits::TP ));
            }else if(ignore == true){
                ignore_cnt++;
                fpVec.push_back(Hits(query,rKey, evalue, Hits::IGN ));
            }else if(fp == true){
                fp_cnt++;
                fpVec.push_back(Hits(query,rKey, evalue, Hits::FP ));
                auc = auc + tp_cnt;
            }
        }
    }
    return EvaluateResult(tp_cnt, fp_cnt, ignore_cnt, auc);
}

void writeRoc5Data(std::string roc5ResultFile,
                   std::vector<Roc5Value> roc5Vals,
                   double stepSize) {
    int i = 0;
    std::ofstream roc5Out;
    roc5Out.open (roc5ResultFile+".rocx");
    double auc = 0.0;

    for(double step = 1.0; step >= 0.0-stepSize; step-=stepSize){
        while ((i < roc5Vals.size()) && (roc5Vals[i].roc5val >= step)){
            i++;
        }
        roc5Out << std::max(0.0, step) << " " << ((float)i)/((float)roc5Vals.size()) << "\n";
        auc = auc + stepSize*((float)i)/((float)roc5Vals.size());
        //std::cout <<  "i = " << i << " x = " << std::max(0.0, step)  << " auc = " << auc << std::endl;
    }
    std::cout << "ROC5 AUC: " << auc << std::endl;
    roc5Out.close();
}

void writeRocData(std::string rocFilePath, std::vector<Hits> hits, size_t binSize) {
    std::ofstream roc5Out;
    roc5Out.open(rocFilePath + ".roc");
    size_t tp_cnt = 0;
    size_t fp_cnt = 0;
    size_t step_size = hits.size() / binSize;
    for (size_t i = 0; i < hits.size(); i++) {
        if (hits[i].status == Hits::TP) {
            tp_cnt++;
        } else if (hits[i].status == Hits::FP) {
            fp_cnt++;
        }

        if (i % step_size == 0) {
            roc5Out << fp_cnt << "\t" << tp_cnt << "\n";
        }

    }
    roc5Out.close();
}

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
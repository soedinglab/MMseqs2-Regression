//
// Created by mad on 9/22/15.
//
#include "EvaluateResults.h"
#include <stddef.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <iomanip>
#include <algorithm>

extern "C" {
#include "ffindex.h"
#include "ffutil.h"
}

#include "PatternCompiler.h"
#include "EvaluateResults.h"
#include "kseq.h"
#define MAX_FILENAME_LIST_FILES 4096
KSEQ_INIT(int, read)

struct SCOPPos{
    std::string fold;
    std::string superFam;
    std::string fam;
    int startPos;
    int endPos;
    double evalue;
    SCOPPos(std::string scop, int sPos, int ePos, double eval){
        fam = scop;
        superFam = scop;
        superFam = superFam.erase(superFam.find_last_of("."), std::string::npos);
        fold = scop;
        fold     = fold.erase(fold.find_last_of("."), std::string::npos);
        fold     = fold.erase(fold.find_last_of("."), std::string::npos);
        evalue = eval;
        startPos = sPos;
        endPos = ePos;
    }
};


struct HitRecord{
    std::string id;
    double eval;
    int qStartPos;
    int qEndPos;
    int tStartPos;
    int tEndPos;
    HitRecord(std::string targetId, double e, int qStartPos, int qEndPos, int tStartPos, int tEndPos)
            :id(targetId), eval(e), qStartPos(qStartPos), qEndPos(qEndPos), tStartPos(tStartPos), tEndPos(tEndPos)
    {}
};


int overlap(int min1, int max1, int min2, int max2);

void parseHitRecordM8(std::string query, std::string resFileName, std::vector<HitRecord> &resultVector, double resSize) {
    static bool isReadIn = false;
    static std::map<std::string, std::vector<HitRecord> > resLookup;
    if(isReadIn == false){
        std::cout << "Read in m8 " << resFileName << std::endl;
        size_t resSizeInt = resSize;
        std::ifstream infile(resFileName);
        std::string line;
        PatternCompiler keyRegex("[^[:space:]]+");
        //0               1       2       3       4       5       6       7       8       9       10              11
        //P40818_1        P40818  0.229   989     256     0       0       989     0       983     7.96E-70        269
        while (std::getline(infile, line))
        {
            std::vector<std::string> tmpRes = keyRegex.getAllMatches(line.c_str(), line.size());
            std::string key = tmpRes[0];
            if(resLookup.find(key)== resLookup.end()) {
                resLookup[key] = std::vector<HitRecord>();
            }
            std::string targetkey = tmpRes[1];
            std::string evalStr = tmpRes[10];
            //query start and end pos
            int qStartPos = strtol(tmpRes[6].c_str(),NULL, 10);
            int qEndPos = strtol(tmpRes[7].c_str(), NULL, 10);
            int tStartPos = strtol(tmpRes[8].c_str(), NULL, 10);
            int tEndPos = strtol(tmpRes[9].c_str(), NULL, 10);
            double eval = atof(evalStr.c_str());
            if(resLookup[key].size() < resSizeInt){
                resLookup[key].push_back(HitRecord(targetkey,eval, qStartPos, qEndPos, tStartPos, tEndPos));
            }
        }
        infile.close();
        std::map<std::string, std::vector<HitRecord>>::iterator it;

        for ( it = resLookup.begin(); it != resLookup.end(); it++ ) {
            std::set<std::string> removeDub;
            std::vector<HitRecord> single;

            for(size_t i = 0; i < it->second.size(); i++){
                if(removeDub.find(it->second[i].id) == removeDub.end()){
                    single.push_back(it->second[i]);
                }
                removeDub.insert(it->second[i].id);
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

std::vector<HitRecord> readResultFileM8(std::string query, std::string resFileName, double resSize) {
    std::vector<HitRecord> resultVector;
    std::string extention = resFileName.substr(resFileName.find_last_of(".") + 1);
    parseHitRecordM8(query, resFileName, resultVector, resSize);
    return resultVector;
}

void readFamDefFromFasta(std::string fasta_path, std::unordered_map<std::string, std::vector<SCOPPos> *> &queryScopLookup,
                         std::unordered_map<std::string, size_t > &supFamSizeLookup, bool readEval) {
    FILE * fasta_file = fopen(fasta_path.c_str(), "r");
    kseq_t *seq = kseq_init(fileno(fasta_file));
    size_t entries_num = 0;

    PatternCompiler scopDomainRegex("[a-z]+\\.[0-9]+\\.[0-9]+\\.[0-9]+");
    std::set<std::string> scopSuperFam;
    while (kseq_read(seq) >= 0) {
        if (seq->name.l == 0) {
            std::cout << "Fasta entry: " << entries_num << " is invalid." << std::endl;
            exit(EXIT_FAILURE);
        }
        const std::string currQuery(seq->name.s);
        if(queryScopLookup.find(currQuery)== queryScopLookup.end()) {
            queryScopLookup[currQuery] = new std::vector<SCOPPos>();
        }
        std::vector<SCOPPos> * queryDomainVector = queryScopLookup[currQuery];
        //>P20932_1 c.1.4.1 |7.4e-185 |d1p4ca_ |1-353 |4-375
        std::vector<std::string> splits = split(std::string(seq->comment.s), "|");
        std::vector<std::string> evals;
        if(readEval == true){
            evals = split(splits[1]," ");
        }
        std::vector<std::string> startEndPos = split( splits[4],"-");
        int alnStart = strtoull(startEndPos[0].c_str(), NULL, 10);
        int alnEnd = strtoull(startEndPos[1].c_str(), NULL, 10);

        std::string s(splits[0].c_str(), splits[0].size());
        std::vector<std::string> domains = scopDomainRegex.getAllMatches(s.c_str(), splits[0].size());
        std::set<std::string> scopDomains(domains.begin(), domains.end());

        int i = 0;
        // there will be just SCOP in this benchmark
        for(std::set<std::string>::iterator it = scopDomains.begin(); it != scopDomains.end(); it++) {
            std::string currScopDomain = *it;
            double eval = (readEval == true) ? strtod(evals[i].c_str(), NULL) : 0.0;
            // increase the scop domain count
            SCOPPos domain = SCOPPos(currScopDomain, alnStart, alnEnd, eval);
            if (supFamSizeLookup.find(domain.fam) == supFamSizeLookup.end()) {
                supFamSizeLookup[domain.fam] = 0;
            }

            supFamSizeLookup[domain.fam]++;
            if(scopSuperFam.find(domain.superFam) == scopSuperFam.end() ){
                supFamSizeLookup[domain.superFam]++;
                scopSuperFam.insert(domain.superFam);
            }
            supFamSizeLookup[domain.fold]++;
            queryDomainVector->push_back(domain);
            i++;
        }
        scopSuperFam.clear();
        entries_num++;
        printProgress(entries_num);

    }
    kseq_destroy(seq);
}


EvaluateResult evaluateResult(std::string query, std::vector<SCOPPos> *qScopIds, std::unordered_map<std::string,
        std::vector<SCOPPos> *> &scopLoopup, std::vector<Hits> &allHitsVec,
                              std::vector<HitRecord> results, size_t rocx) {
    double fp_cnt = 0.0;
    double tp_cnt = 0.0;
    double ignore_cnt = 0.0;
    double auc = 0.0;

//    std::string qSupFam = qFam;
//    qSupFam = qSupFam.erase(qSupFam.find_last_of("."), std::string::npos);
    for (size_t i = 0; i < results.size(); i++) {
        const std::string rKey = results[i].id;
        const double evalue = results[i].eval;
        const int queryAlnStart = results[i].qStartPos;
        const int queryAlnEnd   = results[i].qEndPos;
        const int targetAlnStart   = results[i].tStartPos;
        const int targetAlnEnd   = results[i].tEndPos;


        bool tp = false;
        bool fp = false;
        bool ignore = false;
        std::vector<SCOPPos> * rfamVec;

        PatternCompiler ignore_superfam("^b\\.(67|68|69|70).*");
        PatternCompiler ignoreClass("^e\\..*");

        // if sequence does not have annotations ignore it
        if (scopLoopup.find(rKey) == scopLoopup.end()) {
            tp = false;
            ignore = false;
            fp = true;
            goto outer;
        }
        rfamVec = scopLoopup[rKey];

        for(size_t j = 0; j < rfamVec->size(); j++) {
            for(size_t i = 0; i < qScopIds->size(); i++) {
                SCOPPos qScopId = qScopIds->at(i);
                // AAA|AAAA|AAA
                // BBB|AAAA|BBB
                // START/END ANNOTATION QUERY ( qScopId.startPos, qScopId.endPos)
                const SCOPPos rScopId = rfamVec->at(j);
                // start and end of alignment in result rScopId.startPos, rScopId.endPos,
                int queryDomainOverlap = overlap(queryAlnStart, queryAlnEnd, qScopId.startPos, qScopId.endPos);
                // min overlap of 10 in the domain to be a TP
                if (rScopId.fam.compare(qScopId.fam) == 0 && queryDomainOverlap >= 10) {
                    tp = true;
                    goto outer;
                }
                if (tp == false) {
                    // overlap with domain in target sequence
                    int targetDomainOverlap = overlap(targetAlnStart, targetAlnEnd, rScopId.startPos, rScopId.endPos);

                    bool qSuperFamIgnore = ignore_superfam.isMatch(qScopId.fam.c_str());
                    bool rSuperFamIgnore = ignore_superfam.isMatch(rScopId.fam.c_str());
                    bool qFoldIgnore = ignoreClass.isMatch(qScopId.fam.c_str());
                    bool rFoldIgnore = ignoreClass.isMatch(rScopId.fam.c_str());
                    if (     // ignore same fold but only if overlap exists
                            (rScopId.fold.compare(qScopId.fold) == 0 && queryDomainOverlap >= 10)
                            // ignore certain superfamilies
                        || (qSuperFamIgnore && rSuperFamIgnore)
                           // ignore e class
                        || qFoldIgnore
                        || rFoldIgnore
                        || (targetDomainOverlap >= 10 && queryDomainOverlap == 0) //maybe not annotated sequence
                            ) {
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
                allHitsVec.push_back(Hits(query, rKey, evalue, Hits::TP ));
            }else if(ignore == true){
                ignore_cnt++;
                allHitsVec.push_back(Hits(query, rKey, evalue, Hits::IGN ));
            }else if(fp == true){
                fp_cnt++;
                allHitsVec.push_back(Hits(query, rKey, evalue, Hits::FP ));
                auc = auc + tp_cnt;
            }
        }else{
            if (tp == true)
            {
                allHitsVec.push_back(Hits(query, rKey, evalue, Hits::TP ));
            }else if(ignore == true){
                allHitsVec.push_back(Hits(query, rKey, evalue, Hits::IGN ));
            }else if(fp == true){
                allHitsVec.push_back(Hits(query, rKey, evalue, Hits::FP ));
            }
        }
    }
    return EvaluateResult(tp_cnt, fp_cnt, ignore_cnt, auc);
}

int overlap(int min1, int max1, int min2, int max2) {
    return std::max(0, std::min(max1, max2) - std::max(min1, min2));
}

int main(int argc, char ** argv){
    size_t res_cnt = 0;
    //while( true)
    std::string queryFasta=argv[1];
    std::string targetFasta=argv[2];
    std::string resultFile=argv[3];
    std::string outputResultFile=argv[4];
    double resSize = 1000.0;
    if(argc > 5){
        resSize = atof(argv[5]);
    }

    size_t rocx = 1;
    if(argc > 6){
        rocx = atoi(argv[6]);
    }
    std::unordered_map<std::string, std::vector<SCOPPos>*> scopLoopup;
    std::cout << "Read query fasta" << std::endl;
    std::unordered_map<std::string, size_t > whatever;
    readFamDefFromFasta(queryFasta, scopLoopup,
                        whatever, true);
    std::cout << "Read target fasta " << std::endl;
    std::unordered_map<std::string, size_t> scopSizeLoopup;
    //std::cout << scopLoopup["d12asa_"]->at(0) << " " << famSizeLoopup[scopLoopup["d1acfa_"]->at(0)] <<std::endl;
    readFamDefFromFasta(targetFasta, scopLoopup,
                        scopSizeLoopup, true);
    std::cout << std::endl;
    std::cout << "Read result fasta " << std::endl;

    std::vector<std::pair<size_t, std::string>> supFam;
    for (std::unordered_map<std::string, size_t >::iterator it = scopSizeLoopup.begin();
         it != scopSizeLoopup.end(); it++ ) {
        size_t n = std::count(it->first.begin(), it->first.end(), '.');
        if(n == 2) {
            supFam.push_back(std::make_pair( it->second, it->first));
        }
    }
    std::sort(supFam.begin(), supFam.end());
    double sum = 0.0;
    for(size_t i = 0; i < supFam.size(); i++){
        sum += (double)supFam[i].first;
//        std::cout << supFam[i].second << "\t" << supFam[i].first << std::endl;
    }
    std::cout << "N=" << supFam.size() << " Sum=" << sum << " Avg=" << (sum/ supFam.size());
    std::cout << " Median=" << supFam[supFam.size()/2].first;
    std::cout << " 1/4=" << supFam[supFam.size()/4].first;
    std::cout << " 3.5/4=" << supFam[(supFam.size()/4) * 3.5].first << std::endl;

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
        std::vector<SCOPPos> * qFams = scopLoopup[query];
        std::vector<HitRecord> resIds = readResultFileM8(query, resultFile, resSize);
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
            SCOPPos qFam = qFams->at(i);
            qFamSize = std::min(qFamSize + scopSizeLoopup[qFam.fam], resSize);
            qFamStr.append(qFams->at(i).fam).append(",");
        }
        if(qFamSize > 0){
            double roc5val = all_auc / (rocx * qFamSize);
            if (roc5val > 1.0){
                std::cout << "ROC" << rocx << " = " << roc5val << " for query " << query << ", # family members: " <<
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
    printf( "Query\t\tFam\t\t\t\tRoc5\tFamSize\tTPs\tFP\tResSize\tIGN)\n");
    for(size_t i = 0; i < roc5Vals.size(); i++) {
        Roc5Value roc5Value = roc5Vals[i];
        printf("%s\t\t%-30.30s\t%.7f\t%5zu\t%5zu\t%5zu\t%5zu\t%5zu\n", roc5Value.query.c_str(), roc5Value.qFams.c_str(),
               roc5Value.roc5val, roc5Value.qFamSize, roc5Value.tp_cnt, roc5Value.fp_cnt,
               roc5Value.resultSize, roc5Value.ignore_cnt);
    }
    double fpsWithSmallEvalue=0;
    double EVAL_THRESHOLD = 1E-3;
    std::map<std::string, size_t > mostQueriesWithSmallEval;
    for(size_t i = 0; i < allHits.size(); i++) {
        if(allHits[i].evalue < EVAL_THRESHOLD && allHits[i].status == Hits::FP){
            fpsWithSmallEvalue++;
            mostQueriesWithSmallEval[allHits[i].query]++;
        }
    }
//    std::cout << "Top high scoring queries:" << std::endl;
    std::vector<std::pair<size_t, std::string>> mostQueriesWithSmallEvalVec;
    for (std::map<std::string, size_t >::iterator it = mostQueriesWithSmallEval.begin();
         it != mostQueriesWithSmallEval.end(); it++ ) {
        mostQueriesWithSmallEvalVec.push_back(std::make_pair(it->second, it->first));
    }
    std::sort(mostQueriesWithSmallEvalVec.begin(), mostQueriesWithSmallEvalVec.end());
/*    for (int i = mostQueriesWithSmallEvalVec.size(); i > 0; i--) {
        std::cout << mostQueriesWithSmallEvalVec[i-1].second << " " << mostQueriesWithSmallEvalVec[i-1].first << std::endl;
    }
*/
    std::sort(allHits.begin(), allHits.end(), sortFalsePositvesByEval());

    std::cout << "Top 50 FP:" << std::endl;
    size_t cnt=0;
    for(size_t i = 0; i < allHits.size(); i++) {
        if(allHits[i].status == Hits::FP){
            std::cout << cnt + 1 << ": " << allHits[i].query << " " << allHits[i].target << " " << allHits[i].evalue << std::endl;
            std::vector<SCOPPos> * scopTarget;
            if(scopLoopup.find(allHits[i].target) == scopLoopup.end()){
                scopTarget = NULL;
            }else {
                scopTarget =  scopLoopup[allHits[i].target];
            }
            std::vector<SCOPPos> * scopQuery =  scopLoopup[allHits[i].query];
            std::cout << "Q=";
            for(size_t j = 0; j < scopQuery->size(); j++){
                std::cout << " " << scopQuery->at(j).fam << "(" << scopQuery->at(j).evalue << ")";
            }
            std::cout << std::endl;
            std::cout << "T=";
            if(scopTarget == NULL){
                std::cout << " Inverse";
            }else {
                for(size_t j = 0; j < scopTarget->size(); j++){
                    std::cout << " " << scopTarget->at(j).fam << "(" << scopTarget->at(j).evalue << ")";
                }
            }
            std::cout << std::endl;
            std::cout << std::endl;

            cnt++;
            if(cnt==50)
                break;
        }
    }
    writeAnnoatedResultFile(outputResultFile, allHits);
    kseq_destroy(seq);
    std::cout << res_cnt  << " result lists checked." << std::endl;
    std::cout << early_break_cnt << " result lists did not contain " << rocx << " FPs." << std::endl;
    std::cout << "Results contains " << overall_tp << " TPs and " << overall_fp << " FPs." << std::endl;
    std::cout << "Total FPs " << fpsWithSmallEvalue << " of " << mostQueriesWithSmallEvalVec.size() << " queries have an eval < " << EVAL_THRESHOLD << "." << std::endl;
    std::cout << overall_ignore_cnt << " sequence pairs ignored (different family, same fold)" << std::endl;

    writeRoc5Data(outputResultFile, roc5Vals, 0.01);
    std::sort(allHits.begin(), allHits.end(), sortFalsePositvesByEval());

    writeRocData(outputResultFile, allHits, 10000);
    writeFDRData(outputResultFile, allHits, roc5Vals, 1E-50);

    return 0;
}

void writeAnnoatedResultFile(std::string resultFile, std::vector<Hits> & hits) {
    std::ofstream outFile;
    outFile.open (resultFile+".annotated_result");
    for (size_t i = 0; i < hits.size(); i++) {
        outFile << hits[i].query << "\t" << hits[i].target << "\t" << hits[i].evalue << "\t" << hits[i].status  << "\n";
    }
    outFile.close();
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

void writeFDRData(std::string roc5ResultFile,
                  std::vector<Hits> & hits,
                  std::vector<Roc5Value> & rocVal,
                  double stepSize) {
    int i = 0;
    std::ofstream fdrOut;
    fdrOut.open (roc5ResultFile + ".fdr");
    double tp = 0;
    double fp = 0;
    size_t cnt = 0;
    std::map<std::string, size_t > queryToResSize;
    for(size_t i = 0; i < rocVal.size(); i++){
        queryToResSize[rocVal[i].query] = rocVal[i].resultSize;
    }

    for(double step = 0.0; step <= 10000.0; step = stepSize ){
        while ((i < hits.size()) && (hits[i].evalue <= step)){
            double size = static_cast<double>(std::max((size_t )1, queryToResSize[hits[i].query]));
            tp += ((hits[i].status == Hits::TP)/size);
            fp += ((hits[i].status == Hits::FP)/size);
            i++;
        }
        fdrOut << std::fixed << std::setprecision(1) << std::scientific  << std::max(0.0, step) << "\t" << std::fixed << std::setprecision(6) << (tp) / (fp + tp) << "\t" << (fp) / (fp + tp) << "\n";
        cnt++;
        stepSize *= 1.5;
    }
    fdrOut.close();
}


void writeRoc5Data(std::string roc5ResultFile,
                   std::vector<Roc5Value> & roc5Vals,
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

void writeRocData(std::string rocFilePath, std::vector<Hits> & hits, size_t binSize) {
    std::ofstream roc5Out;
    roc5Out.open(rocFilePath + ".roc");
    size_t tp_cnt = 0;
    size_t fp_cnt = 0;
    size_t step_size = hits.size() / binSize;
    step_size = std::max(step_size, (size_t ) 1 );
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

#include "PatternCompiler.h"
#include <iostream>
#include <cstring>

int main(){
    PatternCompiler ignore_superfam("^b\\.(67|68|69|70).*");
    bool qSuperFamIgnore = ignore_superfam.isMatch("e.12.1.4");
    bool rSuperFamIgnore = ignore_superfam.isMatch("b.67.1.4");

    PatternCompiler ignore_fold("^e\\..*");
    bool qFoldIgnore = ignore_fold.isMatch("e.1.4.2");
    bool rFoldIgnore = ignore_fold.isMatch("b.12.3.1");

    PatternCompiler scopDomainRegex("[a-z]+\\.[0-9]+\\.[0-9]+\\.[0-9]+");
    const char* scops = "b.67.1.4 b.67.1.5 b.67.1.6 b.67.1.7";
    std::cout << scopDomainRegex.isMatch(scops) << std::endl;
    std::vector<std::string> parsed = scopDomainRegex.getAllMatches(scops, strlen(scops));
    for (std::vector<std::string>::const_iterator it = parsed.begin(); it != parsed.end(); ++it) {
    	std::cout << *it << std::endl;
    }

    std::cout << "test2\n";
    std::cout << qSuperFamIgnore << " " << rSuperFamIgnore << " " << qFoldIgnore << " " << rFoldIgnore << std::endl;
    return 0;
}
#include <regex>
#include <iostream>

int main(){
    std::regex ignore_superfam("^b\\.(67|68|69|70).*");
    bool qSuperFamIgnore = std::regex_match("e.12.1.4", ignore_superfam);
    bool rSuperFamIgnore = std::regex_match("b.67.1.4", ignore_superfam);
    std::regex ignore_fold("^e\\..*");
    bool qFoldIgnore = std::regex_match("e.1.4.2", ignore_fold);
    bool rFoldIgnore = std::regex_match("b.12.3.1", ignore_fold);

    std::cout << qSuperFamIgnore << " " << rSuperFamIgnore << " " << qFoldIgnore << " " << rFoldIgnore << std::endl;
    return 0;
}
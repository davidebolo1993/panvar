// M15-G: `genotype_assessment_test show <prefix>` prints the assessment_v1 view of a v2 run prefix exactly as
// tests/analysis/m15g/assessment_spec.py render() does, without writing anything (the equivalence harness compares the
// two); `write <prefix>` writes it (tests only).
#include "panvar/genotype_assessment.hpp"

#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc != 3 || (std::string(argv[1]) != "show" && std::string(argv[1]) != "write")) {
        std::cerr << "usage: genotype_assessment_test show|write <prefix>\n";
        return 2;
    }
    try {
        if (std::string(argv[1]) == "show") std::cout << panvar::assessment::render(panvar::assessment::assess_prefix(argv[2]));
        else panvar::assessment::write_assessment(argv[2]);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}

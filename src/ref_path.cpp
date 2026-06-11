#include "panvar/ref_path.hpp"

#include <cctype>
#include <exception>
#include <string>

namespace panvar {

ParsedReferencePath parse_reference_path_label(const std::string& path_name) {
    ParsedReferencePath out;
    out.chrom = path_name;

    auto pick_chrom = [&](const std::string& prefix) {
        if (prefix.empty()) {
            return std::string(path_name);
        }
        const std::size_t hash = prefix.rfind('#');
        if (hash != std::string::npos && hash + 1 < prefix.size()) {
            return prefix.substr(hash + 1);
        }
        const std::size_t chr_pos = prefix.find("chr");
        if (chr_pos != std::string::npos) {
            return prefix.substr(chr_pos);
        }
        return prefix;
    };

    const std::size_t colon = path_name.rfind(':');
    const std::size_t dash = path_name.rfind('-');
    if (colon != std::string::npos && dash != std::string::npos && dash > colon + 1 && dash + 1 < path_name.size()) {
        const std::string start_str = path_name.substr(colon + 1, dash - colon - 1);
        const std::string end_str = path_name.substr(dash + 1);
        bool digits_only = !start_str.empty() && !end_str.empty();
        for (const char c : start_str) {
            if (!std::isdigit(static_cast<unsigned char>(c))) {
                digits_only = false;
                break;
            }
        }
        if (digits_only) {
            for (const char c : end_str) {
                if (!std::isdigit(static_cast<unsigned char>(c))) {
                    digits_only = false;
                    break;
                }
            }
        }
        if (digits_only) {
            try {
                const std::size_t parsed_start = static_cast<std::size_t>(std::stoull(start_str));
                const std::size_t parsed_end = static_cast<std::size_t>(std::stoull(end_str));
                if (parsed_start > 0 && parsed_end >= parsed_start) {
                    out.has_interval = true;
                    out.region_start_1based = parsed_start;
                    out.chrom = pick_chrom(path_name.substr(0, colon));
                    return out;
                }
            } catch (const std::exception&) {
            }
        }
    }

    out.chrom = pick_chrom(path_name);
    return out;
}

} // namespace panvar

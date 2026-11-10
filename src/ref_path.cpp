#include "panvar/ref_path.hpp"

#include <cctype>
#include <exception>
#include <string>

namespace panvar {

std::vector<std::size_t> match_reference_name(const std::vector<std::string>& names, const std::string& query) {
    for (std::size_t i = 0; i < names.size(); ++i)
        if (names[i] == query) return {i};
    std::vector<std::size_t> hits;
    for (std::size_t i = 0; i < names.size(); ++i)
        if (names[i].find(query) != std::string::npos) hits.push_back(i);
    if (!hits.empty()) return hits;
    const auto lower = [](std::string t) {
        for (char& c : t) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return t;
    };
    const std::string q = lower(query);
    for (std::size_t i = 0; i < names.size(); ++i)
        if (lower(names[i]).find(q) != std::string::npos) hits.push_back(i);
    return hits;
}

ParsedReferencePath parse_reference_path_label(const std::string& path_name, std::size_t path_length_bp,
                                               bool one_based) {
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
        // A W line with a start and no end, named "start-*" (1-based): the start is all positions need.
        if (one_based && end_str == "*" && !start_str.empty() &&
            start_str.find_first_not_of("0123456789") == std::string::npos && start_str.size() < 19) {
            out.has_interval = true;
            out.region_start_1based = static_cast<std::size_t>(std::stoull(start_str));
            out.chrom = pick_chrom(path_name.substr(0, colon));
            return out;
        }
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
                if (parsed_end >= parsed_start) {
                    const bool inclusive = (one_based && parsed_start > 0) ||
                                           (path_length_bp > 0 && parsed_start > 0 &&
                                            parsed_end - parsed_start + 1 == path_length_bp);
                    out.has_interval = true;
                    out.half_open = !inclusive;
                    out.length_mismatch = path_length_bp > 0 &&
                        (inclusive ? parsed_end - parsed_start + 1 != path_length_bp
                                   : parsed_end - parsed_start != path_length_bp);
                    out.region_start_1based = inclusive ? parsed_start : parsed_start + 1;
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

// Gates for the genotype command's release contract: the declared-library file format, the order-independent
// resolution of --adaptive-repair and --repair-rounds, the provenance of the diagnostic cache flags, the help text,
// and the public headline of a run that stops before its blocks are built.

#include "panvar/genotype_command.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <iterator>
#include <set>
#include <vector>

namespace {

std::string dir;
std::string source_dir_for_profiles;   // M15-J: the evidence profiles are read from the source tree when it is given
std::string path(const std::string& name) { return dir + "/" + name; }
void require(bool ok, const std::string& what) { if (!ok) throw std::runtime_error(what); }

void write(const std::string& file, const std::string& text) {
    std::ofstream out(file);
    out << text;
}

const std::string kLibrary = "insert_mean_bp 300\ninsert_sd_bp 30\nsupport_sigmas 4\nsubstitution_error 0.002\n"
                             "outlier_weight 0.001\nfragment_start_intensity 0.085\n";

// The message of the error a call throws, or "" when it does not throw.
std::string error_of(const std::function<void()>& call) {
    try {
        call();
    } catch (const std::exception& error) {
        return error.what();
    }
    return "";
}

std::string library_error(const std::string& name, const std::string& text) {
    write(path(name), text);
    std::optional<double> background;
    return error_of([&] { panvar::read_library(path(name), background); });
}

void test_library_file() {
    // Valid: blank lines, full-line and indented comments, a trailing comment and a CR line ending are accepted, and the
    // values are exactly the declared ones.
    write(path("lib_ok.txt"), "# declared library\n\n" + kLibrary + "   # indented comment\nmarker_background 0.01 # optional\r\n");
    std::optional<double> background;
    const panvar::LibraryModel model = panvar::read_library(path("lib_ok.txt"), background);
    require(model.insert_mean_bp == 300.0 && model.insert_sd_bp == 30.0 && model.support_sigmas == 4.0 &&
                model.substitution_error == 0.002 && model.outlier_weight == 0.001 && model.fragment_start_intensity == 0.085,
            "a valid library must give exactly its declared values");
    require(background && *background == 0.01, "the optional marker_background must be read");
    require(model.provenance == "declared: " + path("lib_ok.txt"), "the provenance must name the file");
    std::optional<double> none;
    panvar::read_library((write(path("lib_nobg.txt"), kLibrary), path("lib_nobg.txt")), none);
    require(!none, "an absent marker_background must stay absent");

    const auto refused = [&](const std::string& name, const std::string& text, const std::vector<std::string>& needles,
                             const std::string& what) {
        const std::string message = library_error(name, text);
        require(!message.empty(), what + ": must be refused");
        for (const std::string& needle : needles) {
            require(message.find(needle) != std::string::npos, what + ": the error must contain '" + needle + "', got: " + message);
        }
    };
    // Every required key, missing in turn.
    for (const char* key : {"insert_mean_bp", "insert_sd_bp", "support_sigmas", "substitution_error", "outlier_weight",
                            "fragment_start_intensity"}) {
        std::string text;
        std::istringstream lines(kLibrary);
        for (std::string line; std::getline(lines, line);) {
            if (line.rfind(std::string(key) + " ", 0) != 0) text += line + "\n";
        }
        refused(std::string("lib_missing_") + key, text, {path(std::string("lib_missing_") + key), "missing", key}, "missing key");
    }
    refused("lib_unknown.txt", kLibrary + "gap_open 0.001\n", {"lib_unknown.txt:7", "unknown library key 'gap_open'"},
            "an undeclarable gap parameter");
    refused("lib_typo.txt", kLibrary + "marker_backgrond 0.01\n", {":7", "unknown library key 'marker_backgrond'"},
            "a misspelt optional key");
    refused("lib_repeat.txt", kLibrary + "insert_mean_bp 350\n", {"lib_repeat.txt:7", "'insert_mean_bp' is repeated", "line 1"},
            "a repeated key");
    refused("lib_trailing.txt", "insert_mean_bp 300 bp\n" + kLibrary.substr(kLibrary.find('\n') + 1),
            {"lib_trailing.txt:1", "after the value of 'insert_mean_bp'", "'bp'"}, "a trailing token");
    refused("lib_junk.txt", "insert_mean_bp 300bp\n" + kLibrary.substr(kLibrary.find('\n') + 1),
            {":1", "'insert_mean_bp' is not a finite number", "300bp"}, "a number with a suffix");
    refused("lib_novalue.txt", "insert_mean_bp\n" + kLibrary.substr(kLibrary.find('\n') + 1), {":1", "expected 'key value'"},
            "a key without a value");
    refused("lib_inf.txt", "insert_mean_bp inf\n" + kLibrary.substr(kLibrary.find('\n') + 1), {":1", "not a finite number"},
            "an infinite value");
    refused("lib_nan.txt", "insert_mean_bp nan\n" + kLibrary.substr(kLibrary.find('\n') + 1), {":1", "not a finite number"},
            "a NaN value");
    refused("lib_negative_sd.txt", "insert_sd_bp -30\n" + kLibrary.substr(0, kLibrary.find("insert_sd_bp")) +
                                       kLibrary.substr(kLibrary.find("support_sigmas")),
            {"lib_negative_sd.txt", "invalid library model", "standard deviation"}, "a negative insert sd");
    refused("lib_error.txt", "substitution_error 0.9\n" + kLibrary.substr(0, kLibrary.find("substitution_error")) +
                                 kLibrary.substr(kLibrary.find("outlier_weight")),
            {"lib_error.txt", "invalid library model"}, "a substitution error outside (0, 0.75)");
    std::cout << "  library file: valid input exact; missing, unknown, repeated, trailing, non-finite and invalid refused\n";
}

struct Outcome {
    std::string error;            // non-empty: refused
    std::uint64_t rounds = 0;
    bool tie_branching = false;
    bool operator==(const Outcome& o) const {
        return error.empty() == o.error.empty() && (!error.empty() || (rounds == o.rounds && tie_branching == o.tie_branching));
    }
};

Outcome parse(const std::vector<std::string>& args, panvar::GenotypeCommandOptions* keep = nullptr) {
    panvar::GenotypeCommandOptions options;
    std::ostringstream sink;
    Outcome out;
    out.error = error_of([&] { panvar::parse_genotype_arguments(args, options, sink); });
    out.rounds = options.repair_rounds;
    out.tie_branching = options.repair_tie_branching;
    if (keep) *keep = options;
    return out;
}

// M15-J: one public reporting contract. The preset's sequence_first_v2 is recorded with its source; no flag or profile
// selects another; M15-K: the unrefined view reports v2 too (from its already-scored union pool) (formerly a reported dependency,
// unchanged in M15-J). The internal API keeps every contract for regression tests.
void test_output_schema() {
    panvar::GenotypeCommandOptions kept;
    require(parse({}, &kept).error.empty() && kept.output_schema == "sequence_first_v2" &&
                kept.setting_sources.at("output_schema") == "preset short-read (default)",
            "unset, the contract is the preset's sequence_first_v2, recorded as such");
    for (const char* value : {"legacy", "sequence_first_v1", "sequence_first_v2"}) {
        require(parse({"--output-schema", value}).error.find("--output-schema was removed") != std::string::npos,
                std::string("--output-schema ") + value + " is refused by name");
    }
    for (const char* value : {"legacy", "sequence_first_v1", "sequence_first"}) {
        write(path("schema_profile.tsv"), std::string("profile\tschema\noutput_schema\t") + value + "\n");
        require(parse({"--profile", path("schema_profile.tsv")}).error.find("'output_schema' cannot select a reporting contract") !=
                    std::string::npos,
                std::string("a profile cannot select another reporting contract (") + value + ")");
    }
    write(path("schema_profile.tsv"), "profile\tschema\noutput_schema\tsequence_first_v2\n");
    require(parse({"--profile", path("schema_profile.tsv")}, &kept).error.empty() && kept.output_schema == "sequence_first_v2",
            "a profile may restate the fixed contract (the evidence profiles do)");
    for (const char* locus : {"c4", "cyp2d6", "acot"}) {
        const std::string evidence = source_dir_for_profiles + "/profiles/genotype_experimental_sr_v1_" + locus + ".tsv";
        if (source_dir_for_profiles.empty()) break;
        require(parse({"--profile", evidence}, &kept).error.empty() && kept.output_schema == "sequence_first_v2",
                std::string("the evidence profile still parses: ") + locus);
    }
    for (const std::vector<std::string>& unrepaired :
         {std::vector<std::string>{"--adaptive-repair", "off"}, std::vector<std::string>{"--repair-rounds", "0"}}) {
        require(parse(unrepaired, &kept).error.empty() && kept.output_schema == "sequence_first_v2" &&
                    kept.setting_sources.at("output_schema") == "preset short-read (default)",
                "M15-K: the unrefined view reports the v2 contract too, recorded from the preset (no inactive marker)");
    }
    std::cout << "  reporting contract: preset v2 for refined and unrefined runs, recorded; no flag or profile selects another\n";
}

// M15-D: the accelerated configuration is the preset's, recorded per key; preset features whose requirements a flag
// removed are switched off and recorded inactive; the neighbourhood bound is the model's own.
void test_preset_configuration() {
    panvar::GenotypeCommandOptions kept;
    require(parse({}, &kept).error.empty(), "the bare preset parses");
    for (const char* key : {"union_scoring_backend", "adaptive_window_discovery", "indexed_union_store_layout",
                            "indexed_union_evidence_layout", "adaptive_fused_seed", "combination_prepared_library_terms",
                            "gap_recursion_backend", "dosage_census_backend", "union_fragment_blocks",
                            "repair_round1_origin_admission", "repair_junction_proposals", "early_incumbent", "compound_annotation",
                            "output_schema", "emission_cache_bytes", "indexed_union_store_bytes", "footprint_limit_bytes",
                            "max_state_classes", "repair_max_neighbours"}) {
        require(kept.setting_sources.count(key) && kept.setting_sources.at(key) == "preset short-read (default)",
                std::string("preset key recorded from the preset: ") + key);
    }
    require(kept.union_scoring_backend == "indexed_fine_experimental" && kept.adaptive_window_discovery == "fine_experimental" &&
                kept.indexed_union_evidence_layout == "fine_views_experimental" && kept.gap_recursion_fixed_geometry &&
                kept.dosage_census_interned && kept.union_fragment_blocks == "fine_core_4096_experimental" &&
                kept.repair_round1_origin_admission == "reuse_only_experimental" && kept.repair_junction_proposals == "experimental_v1" &&
                kept.early_incumbent == "experimental_v1" && kept.compound_annotation == "contract_v2" &&
                kept.footprint_limit_bytes == 6442450944ULL && kept.max_state_classes == 600000000ULL,
            "the preset's accelerated values are the release candidate's");
    require(kept.repair_max_neighbours_model && kept.repair_max_neighbours == 0, "model: no external neighbourhood cap is compared");
    require(kept.time_limit_seconds == 0 && kept.setting_sources.count("time_limit_seconds") == 0, "the preset sets no time limit");
    require(kept.dosage_mode == "off" && kept.dosage_mode_source == "default", "the preset runs no dosage stage");
    require(parse({"--repair-max-neighbours", "958"}, &kept).error.empty() && !kept.repair_max_neighbours_model &&
                kept.repair_max_neighbours == 958,
            "a numeric neighbourhood cap still overrides the model bound");
    write(path("cap.tsv"), "profile\tcap\nrepair_max_neighbours\t1004\n");
    require(parse({"--profile", path("cap.tsv")}, &kept).error.empty() && kept.repair_max_neighbours == 1004 &&
                kept.setting_sources.at("repair_max_neighbours").rfind("profile ", 0) == 0,
            "a profile's numeric cap overrides the preset");
    for (const std::vector<std::string>& unrepaired :
         {std::vector<std::string>{"--adaptive-repair", "off"}, std::vector<std::string>{"--repair-rounds", "0"}}) {
        require(parse(unrepaired, &kept).error.empty() && kept.output_schema == "sequence_first_v2" && kept.early_incumbent == "off" &&
                    kept.compound_annotation == "contract_v2" && kept.repair_junction_proposals == "off",
                "the unrepaired view switches the preset's repair-only features off instead of refusing; v2 and its compound "
                "annotation stay on (M15-K)");
        require(kept.setting_sources.at("compound_annotation") == "preset short-read (default)",
                "the compound annotation is not marked inactive in the unrefined view");
        for (const char* key : {"early_incumbent", "repair_junction_proposals"}) {
            require(kept.setting_sources.at(key).find("(inactive: ") != std::string::npos,
                    std::string("the switched-off preset key is recorded inactive with its reason: ") + key);
        }
    }
    write(path("explicit_early.tsv"), "profile\texplicit-early\nearly_incumbent\texperimental_v1\n");
    require(parse({"--profile", path("explicit_early.tsv"), "--adaptive-repair", "off"}).error.find("early_incumbent") !=
                std::string::npos,
            "the same feature set explicitly keeps its refusal");
    std::cout << "  preset: accelerated keys recorded from the preset, model neighbourhood bound, no time limit or dosage; "
                 "inactive preset features recorded, explicit ones refused\n";
}

// M15-D: dosage runs only when requested (--dosage-check, or an explicit expert mode); no schema implies it; v2 refuses repair.
// M15-E: build and load cannot be combined; a build refuses sample-only arguments
void test_panel_index_flags() {
    require(parse({"--build-panel-index", "x.pvpl", "--panel-index", "x.pvpl"}).error.find("cannot be combined") != std::string::npos,
            "--build-panel-index with --panel-index is refused");
    for (const char* sample_arg : {"--reads1", "--reads2", "--library", "--out", "--out-prefix", "-o"}) {
        require(parse({"--build-panel-index", "x.pvpl", sample_arg, "v"}).error.find("is a sample argument") != std::string::npos,
                std::string("--build-panel-index refuses ") + sample_arg);
    }
    require(parse({"--build-panel-index", "x.pvpl", "--dosage-check"}).error.find("is a sample argument") != std::string::npos,
            "--build-panel-index refuses --dosage-check");
    panvar::GenotypeCommandOptions kept;
    require(parse({"--panel-index", "x.pvpl", "--reads1", "a", "--reads2", "b", "--threads", "3", "--dosage-check"}, &kept).error.empty() &&
                kept.panel_index_path == "x.pvpl",
            "--panel-index accepts every sample argument");
    std::cout << "  panel-index flags: build + load refused; a build refuses sample arguments; a load accepts them\n";
}

void test_dosage_check_contract() {
    panvar::GenotypeCommandOptions kept;
    require(parse({}, &kept).error.empty() && kept.dosage_mode == "off" && !kept.dosage_mode_explicit &&
                kept.dosage_mode_source == "default",
            "no flag: dosage off (not evaluated), source default");
    require(parse({"--dosage-check"}, &kept).error.empty() && kept.dosage_mode == "audit" && kept.dosage_mode_explicit &&
                kept.dosage_mode_source == "flag --dosage-check",
            "--dosage-check runs the audit, source flag --dosage-check");
    for (const std::vector<std::string>& removed :
         {std::vector<std::string>{"--dosage-mode", "audit"}, std::vector<std::string>{"--dosage-mode", "repair"},
          std::vector<std::string>{"--dosage-mode", "off"}, std::vector<std::string>{"--dosage-audit", "shadow"},
          std::vector<std::string>{"--dosage-check", "--dosage-mode", "repair"}}) {
        const std::string error = parse(removed).error;
        require(error.find("was removed: dosage is an audit only") != std::string::npos && error.find("--dosage-check") != std::string::npos,
                "a removed dosage spelling is refused with a pointer to --dosage-check: " + removed[0] + " " + removed[1]);
    }
    for (const char* key : {"dosage_mode", "dosage_audit", "dosage_check"}) {
        write(path("dosage_profile.tsv"), std::string("profile\tdosage\n") + key + "\trepair\n");
        require(parse({"--profile", path("dosage_profile.tsv"), "--dosage-check"}).error.find("is not a profile setting") !=
                    std::string::npos,
                std::string("a profile cannot set the dosage stage: ") + key);
    }
    // S1-D.1: the shared resolver, as the direct command API and the prepared route apply it
    {
        panvar::GenotypeCommandOptions o;   // dosage fields as constructed
        o.output_schema = "sequence_first_v2";
        require(!panvar::resolve_dosage_contract(o) && o.dosage_mode == "off" && o.dosage_mode_source == "default" &&
                    !o.dosage_mode_explicit,
                "API: v2 with untouched dosage fields stays off (M15-D: no implied audit)");
        require(!panvar::resolve_dosage_contract(o) && o.dosage_mode == "off", "API: resolution is idempotent");
    }
    for (const bool flagged : {true, false}) {
        panvar::GenotypeCommandOptions o;
        o.output_schema = "sequence_first_v2";
        o.dosage_mode = "repair";
        o.dosage_mode_explicit = flagged;
        const std::optional<std::string> refused = panvar::resolve_dosage_contract(o);
        require(refused && refused->find("internal API with the legacy contract") != std::string::npos,
                std::string("API: v2 repair is refused (") + (flagged ? "declared" : "set by value") + ")");
    }
    {
        panvar::GenotypeCommandOptions o;
        o.output_schema = "sequence_first_v2";
        o.dosage_mode = "audit";
        require(!panvar::resolve_dosage_contract(o) && o.dosage_mode == "audit" && o.dosage_mode_source == "caller_options",
                "API: v2 audit set by value is a declaration");
    }
    for (const std::string schema : {"legacy", "sequence_first_v1"}) {
        panvar::GenotypeCommandOptions o;
        o.output_schema = schema;
        require(!panvar::resolve_dosage_contract(o) && o.dosage_mode == "off" && o.dosage_mode_source == "default",
                "API: " + schema + " keeps the built-in default");
        o.dosage_mode = "repair";
        require(!panvar::resolve_dosage_contract(o) && o.dosage_mode == "repair", "API: " + schema + " accepts repair");
    }
    {
        panvar::GenotypeCommandOptions o;
        o.dosage_mode = "bogus";
        require(panvar::resolve_dosage_contract(o).has_value(), "API: an unknown dosage mode is refused");
    }
    std::cout << "  dosage: off unless --dosage-check; the expert spellings and profile keys refused; repair only through "
                 "the internal API's legacy contract; one resolver for CLI and API\n";
}

void test_repair_option_order() {
    using Args = std::vector<std::string>;
    const Args on = {"--adaptive-repair", "on"}, off = {"--adaptive-repair", "off"};
    const Args r5 = {"--repair-rounds", "5"}, r0 = {"--repair-rounds", "0"};
    const Args tie_on = {"--repair-tie-branching", "on"}, tie_off = {"--repair-tie-branching", "off"};
    struct Case { std::string name; Args a, b; Outcome want; };
    const auto ok = [](std::uint64_t rounds, bool tie) { Outcome o; o.rounds = rounds; o.tie_branching = tie; return o; };
    Outcome refuse;
    refuse.error = "refused";
    const std::vector<Case> cases = {
        {"neither", {}, {}, ok(32, true)},
        {"on", on, {}, ok(32, true)},
        {"off", off, {}, ok(0, false)},
        {"rounds 5", r5, {}, ok(5, true)},
        {"rounds 0", r0, {}, ok(0, false)},
        {"on + rounds 5", on, r5, ok(5, true)},
        {"on + rounds 0", on, r0, refuse},
        {"off + rounds 5", off, r5, refuse},
        {"off + rounds 0", off, r0, ok(0, false)},
        {"off + tie on", off, tie_on, refuse},
        {"off + tie off", off, tie_off, ok(0, false)},
        {"on + tie off", on, tie_off, ok(32, false)},
        {"rounds 0 + tie on", r0, tie_on, refuse},
        {"rounds 0 + tie off", r0, tie_off, ok(0, false)},
        {"rounds 5 + tie on", r5, tie_on, ok(5, true)},
        {"rounds 5 + tie off", r5, tie_off, ok(5, false)},
    };
    for (const Case& c : cases) {
        Args ab = c.a, ba = c.b;
        ab.insert(ab.end(), c.b.begin(), c.b.end());
        ba.insert(ba.end(), c.a.begin(), c.a.end());
        const Outcome x = parse(ab), y = parse(ba);
        require(x == y, c.name + ": the two argument orders must give the same result (" + x.error + " / " + y.error + ")");
        require(x == c.want, c.name + ": got rounds " + std::to_string(x.rounds) + ", tie branching " +
                                 (x.tie_branching ? "on" : "off") + (x.error.empty() ? "" : ", refused: " + x.error));
        if (!c.want.error.empty()) {
            require(x.error.find("contradicts") != std::string::npos, c.name + ": the refusal must say what contradicts what: " + x.error);
        }
    }
    // A profile's round cap is kept by --adaptive-repair on, and a profile's 0 is replaced by the preset's 32 (unchanged).
    write(path("rounds11.tsv"), "profile\trounds11\nrepair_rounds\t11\n");
    write(path("rounds0.tsv"), "profile\trounds0\nrepair_rounds\t0\n");
    require(parse({"--profile", path("rounds11.tsv"), "--adaptive-repair", "on"}) == ok(11, true), "profile cap kept by on");
    require(parse({"--adaptive-repair", "on", "--profile", path("rounds0.tsv")}) == ok(32, true), "profile 0 with on gives 32");
    require(parse({"--profile", path("rounds11.tsv")}) == ok(11, true), "profile alone");
    // A round cap of 0 from any source is the unrepaired view: tie branching off, recorded with the cap's source.
    panvar::GenotypeCommandOptions zero;
    require(parse({"--repair-rounds", "0"}, &zero) == ok(0, false) &&
                zero.setting_sources.at("repair_tie_branching") == "flag --repair-rounds",
            "--repair-rounds 0 records tie branching off from the flag");
    require(parse({"--profile", path("rounds0.tsv")}, &zero) == ok(0, false) &&
                zero.setting_sources.at("repair_tie_branching").rfind("profile ", 0) == 0,
            "a profile's 0 rounds turns tie branching off, recorded from the profile");
    std::cout << "  repair options: " << cases.size() << " combinations, both orders identical; contradictions refused\n";
}

// M15-J: removed spellings are refused by name with what replaces them; the cache budgets stay profile keys with their
// provenance; one profile per run.
void test_removed_options() {
    const std::vector<std::pair<std::vector<std::string>, std::string>> removed = {
        {{"--output-schema", "legacy"}, "--output-schema was removed"},
        {{"--dosage-mode", "audit"}, "--dosage-mode was removed"},
        {{"--dosage-audit", "shadow"}, "--dosage-audit was removed"},
        {{"--phase-audit-candidates", "8"}, "--phase-audit-candidates was removed: the K-prefix phase audit is internal"},
        {{"--repair-origin-cache-bytes", "123"}, "set the profile key repair_origin_cache_bytes"},
        {{"--repair-tie-seed-cache-bytes", "456"}, "set the profile key repair_tie_seed_cache_bytes"},
        {{"--dosage-allocation-test"}, "--dosage-allocation-test was removed"},
        {{"--library-diagnostic"}, "--library-diagnostic was removed"},
        {{"--max-candidates", "4"}, "--max-candidates is now --max-mosaic-candidates"},
    };
    for (const auto& [args, needle] : removed) {
        std::vector<std::string> with = {"--graph", "g.gfa", "--threads", "2"};
        with.insert(with.end(), args.begin(), args.end());
        require(parse(with).error.find(needle) != std::string::npos, "refused by name: " + args[0] + " (" + parse(with).error + ")");
        std::vector<std::string> first = args;
        first.insert(first.end(), {"--help"});
        require(parse(first).error.find(needle) != std::string::npos, "refused before --help is honoured: " + args[0]);
    }
    panvar::GenotypeCommandOptions options;
    write(path("cache_profile.tsv"), "profile\tcache\nrepair_origin_cache_bytes\t123\nrepair_tie_seed_cache_bytes\t456\n");
    require(parse({"--profile", path("cache_profile.tsv")}, &options).error.empty() && options.repair_origin_cache_bytes == 123 &&
                options.repair_tie_seed_cache_bytes == 456 &&
                options.setting_sources.at("repair_origin_cache_bytes").rfind("profile ", 0) == 0 &&
                options.setting_sources.at("repair_tie_seed_cache_bytes").rfind("profile ", 0) == 0,
            "the cache budgets remain profile keys, recorded with the profile as source");
    write(path("bad_cache.tsv"), "profile\tbad\nrepair_origin_cache_bytes\t12x\n");
    require(!parse({"--profile", path("bad_cache.tsv")}).error.empty(), "a malformed byte count must be refused");
    write(path("phase_profile.tsv"), "profile\tphase\nphase_audit_candidates\t8\n");
    require(parse({"--profile", path("phase_profile.tsv"), "--adaptive-repair", "off"}).error.find(
                "'phase_audit_candidates' is internal: a profile may only restate 0") != std::string::npos,
            "M15-K: a positive phase audit is refused in a profile (internal to the legacy contract)");
    write(path("phase_profile.tsv"), "profile\tphase\nphase_audit_candidates\t0\n");
    require(parse({"--profile", path("phase_profile.tsv"), "--adaptive-repair", "off"}, &options).error.empty() &&
                options.phase_audit_candidates == 0,
            "a profile may restate phase_audit_candidates 0 (the historical C4 profiles do)");
    for (const char* historical : {"genotype_c4.tsv", "genotype_c4_v6_campaign_local.tsv"}) {
        if (source_dir_for_profiles.empty()) break;
        require(parse({"--profile", source_dir_for_profiles + "/profiles/" + historical}).error.find("phase_audit_candidates") ==
                    std::string::npos,
                std::string("a historical profile restating phase_audit_candidates 0 is not refused for it: ") + historical);
    }
    write(path("plain_profile.tsv"), "profile\tplain\n");
    require(parse({"--profile", path("plain_profile.tsv"), "--profile", path("plain_profile.tsv")}).error.find(
                "--profile was given 2 times") != std::string::npos,
            "a repeated --profile is refused");
    write(path("twice.tsv"), "profile\ttwice\nkmer_size\t31\nkmer_size\t29\n");
    require(parse({"--profile", path("twice.tsv")}).error.find("is set twice") != std::string::npos, "a key set twice is refused");
    write(path("nameless.tsv"), "kmer_size\t31\n");
    require(parse({"--profile", path("nameless.tsv")}).error.find("does not name itself") != std::string::npos,
            "a profile without its profile<TAB>name row is refused");
    write(path("unknown.tsv"), "profile\tunknown\nno_such_key\t1\n");
    require(parse({"--profile", path("unknown.tsv")}).error.find("unknown setting 'no_such_key'") != std::string::npos,
            "an unknown profile key is refused");
    write(path("override.tsv"), "profile\toverride\nkmer_size\t29\nsyncmer_s\t11\n");
    require(parse({"--kmer-size", "27", "--profile", path("override.tsv")}, &options).error.empty() && options.kmer_size == 27 &&
                options.syncmer_s == 11 && options.setting_sources.at("kmer_size") == "flag --kmer-size",
            "precedence: a flag overrides the profile, which overrides the preset, whatever the order");
    std::cout << "  removed options: " << removed.size() << " refused by name (before --help); cache budgets and "
                 "phase_audit_candidates is internal (restating 0 only); one profile per run\n";
}

void test_help_text() {
    panvar::GenotypeCommandOptions options;
    std::ostringstream help, advanced;
    require(panvar::parse_genotype_arguments({"--help"}, options, help) == 0, "--help exits 0");
    require(panvar::parse_genotype_arguments({"--help-advanced"}, options, advanced) == 0, "--help-advanced exits 0");
    const std::string h = help.str(), a_raw = advanced.str();
    // the advanced help is wrapped by its generator: phrases are matched with whitespace runs collapsed
    std::string a;
    for (char c : a_raw) {
        const bool ws = c == ' ' || c == '\n';
        if (ws && !a.empty() && a.back() == ' ') continue;
        a += ws ? ' ' : c;
    }
    // M15-H: the basic help leads with the ordinary command and the everyday options, in order
    std::size_t at = 0;
    for (const char* needle : {"--threads 4 --out <prefix>", "--build-panel-index <locus.pvpl>", "Everyday options:", "--graph FILE",
                               "--reads1 FILE --reads2 FILE", "--out PREFIX", "--reference-path NAME", "--threads N", "--dosage-check",
                               "--build-panel-index FILE", "--panel-index FILE", "-q, --quiet", "--help-advanced", "Outcomes:",
                               "Primary assessment (assessment_v1)", "Exit status:"}) {
        const std::size_t found = h.find(needle, at);
        require(found != std::string::npos, std::string("--help must say, in order: ") + needle);
        at = found;
    }
    for (const char* needle : {"No profile is needed", "estimated from the reads", "never a fallback", "It never changes the\n",
                               "not necessarily cheaper", "still required and still read and hashed", "never rebuilt",
                               "no FASTA is needed", "PROVISIONAL (the best fully evaluated candidate", "NO_PREDICTION (exit 1)",
                               "conditional on the explored candidates", "prior-dependent does not mean the prior alone chose the call",
                               "assessed once\n                       per locus", "coverage and not_assessed list",
                               "a missing check is not a failure", "not calibrated reliability", "no quality threshold",
                               "superseded\ncompatibility diagnostics", "0 PREDICTION, 3 PROVISIONAL, 1 NO_PREDICTION",
                               "2 argument, input or output error", "124 time\nlimit", "125 footprint limit", "sampled"}) {
        require(h.find(needle) != std::string::npos, std::string("--help must say: ") + needle);
    }
    // the headline-filter explanation is gone from the basic help (only named as superseded diagnostics)
    for (const char* gone : {"Headline filters, three independent verdicts", "PASS means the computational", "MODEL_ONLY_PHASE",
                             "--mosaic-candidates", "--adaptive-repair", "--profile FILE", "exact genotype =", "a validated settings profile"}) {
        require(h.find(gone) == std::string::npos, std::string("--help must no longer carry: ") + gone);
    }
    // advanced help: the six purpose groups, in order, then the terms
    at = 0;
    for (const char* needle : {"1. Input and panel overrides", "2. Copying prior and marker / placement model",
                               "3. Candidate generation and bounded refinement", "4. Resource limits, work caps, cache budgets and exact backends",
                               "5. Reporting and diagnostics", "6. Legacy and development compatibility", "Terms (read before interpreting results)"}) {
        const std::size_t found = a.find(needle, at);
        require(found != std::string::npos, std::string("--help-advanced must have, in order: ") + needle);
        at = found;
    }
    for (const char* needle : {"not a biological recombination probability", "INITIAL off-panel mosaics",
                               "not a globally exhaustive search", "it is not a filter on final calls",
                               "the substitution, outlier and gap parameters stay fixed", "Support is conditional on the explored candidates",
                               "does not establish that the prior alone caused the selection", "not independently at each boundary",
                               "coverage counts and not_assessed list", "no quality threshold",
                               "default: model (the model-derived bound) (built-in 670)", "expert profile only",
                               "changes: completion", "changes: computation", "changes: inference", "changes: reporting",
                               "per-mate edit budget", "substitutions + gap bases", "cannot be supplied", "some path spans at least N bp",
                               "SAMPLED every 100 ms", "internal to the legacy contract; a profile may only restate 0",
                               "reported with the v2 contract (its tie and phase outcomes follow the v2 rules",
                               "Profiles (--profile FILE): optional overrides, never required", "a repeated --profile is refused",
                               "literal key<TAB>value", "'profile<TAB>NAME' is required", "profile<TAB>my_overrides",
                               "footprint_limit_bytes<TAB>4294967296", "A flag given beside the profile overrides it",
                               "A profile cannot select a reporting contract", "--output-schema", "removed: genotyping reports the sequence-first v2 contract",
                               "--dosage-mode, --dosage-audit", "removed: dosage is an audit only",
                               "--repair-origin-cache-bytes, --repair-tie-seed-cache-bytes", "default: refused"}) {
        require(a.find(needle) != std::string::npos, std::string("--help-advanced must say: ") + needle);
    }
    require(a.find("(alias --repair-origin-cache-bytes)") == std::string::npos, "the removed cache flag spellings are not documented as aliases");
    require(a.find("calibrated reliab") != std::string::npos && a.find("reliable call") == std::string::npos, "no reliability claim");
    // every public setting and every accepted user key appears in the inventory (nothing left as an unexplained key)
    std::set<std::string> keys;
    for (const auto& r : panvar::genotype_option_inventory()) keys.insert(r.key);
    for (const std::string& k : panvar::genotype_public_settings()) require(keys.count(k) != 0, "inventory lists public setting " + k);
    for (const char* k : {"union_scoring_backend", "union_fragment_blocks", "early_incumbent", "repair_junction_proposals", "output_schema",
                          "compound_annotation", "block_diagnostics", "copying_transition", "block_support_rule", "timing_detail",
                          "indexed_union_store_bytes", "adaptive_fine_store_bytes", "marker_proposal_reuse_table_bytes"}) {
        require(keys.count(k) != 0, std::string("inventory lists profile key ") + k);
    }
    // effective defaults come from the preset and the built-in options
    std::map<std::string, std::string> d;
    for (const auto& r : panvar::genotype_option_inventory()) if (!r.key.empty()) d[r.key] = r.effective_default;
    require(d.at("repair_max_neighbours").rfind("model", 0) == 0 && d.at("footprint_limit_bytes") == "6 GiB (built-in 0)" &&
                d.at("time_limit_seconds") == "0" && d.at("output_schema") == "sequence_first_v2 (built-in legacy)" &&
                d.at("max_mosaic_candidates") == "128" && d.at("switch_probability") == "0.01" && d.at("min_variant_bp") == "50" &&
                d.at("kmer_size") == "31" && d.at("syncmer_s") == "15" && d.at("threads") == "1",
            "inventory defaults equal the parser and preset");
    std::cout << "  help: basic help leads with the everyday command and assessment_v1; advanced help by purpose with values, "
                 "defaults, units, effects and terms; inventory covers every public and profile key\n";
}

// M15-H: the guide's option inventory is exactly the generated one (defaults tied to the parser and the preset)
void test_guide_inventory(const std::string& source_dir) {
    std::ifstream in(source_dir + "/docs/guides/genotype_short_read.md");
    require(static_cast<bool>(in), "the guide is readable");
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::string begin = "<!-- option-inventory:begin -->\n", end = "<!-- option-inventory:end -->";
    const std::size_t b = text.find(begin), e = text.find(end);
    require(b != std::string::npos && e != std::string::npos && e > b, "the guide carries the option-inventory block");
    const std::string block = text.substr(b + begin.size(), e - b - begin.size());
    require(block == panvar::genotype_option_inventory_markdown(),
            "the guide's option inventory equals the generated inventory (regenerate it with genotype_cli_contract_test --print-inventory)");
    std::cout << "  guide: the option inventory equals the generated one\n";
}

std::map<std::string, std::string> locus_row(const std::string& file) {
    std::ifstream in(file);
    std::string header, values;
    std::getline(in, header);
    std::getline(in, values);
    std::vector<std::string> h, v;
    std::istringstream hs(header), vs(values);
    for (std::string c; std::getline(hs, c, '\t');) h.push_back(c);
    for (std::string c; std::getline(vs, c, '\t');) v.push_back(c);
    std::map<std::string, std::string> row;
    for (std::size_t i = 0; i < h.size() && i < v.size(); ++i) row[h[i]] = v[i];
    return row;
}

// A run that stops before its blocks are built (a curation manifest naming no graph path is refused while the locus is
// prepared) is a no-call, and its public headline must say so.
void test_pre_block_no_call() {
    write(path("pre.gfa"), "H\tVN:Z:1.0\nS\t1\t" + std::string(400, 'A') + "\nS\t2\t" + std::string(120, 'C') + "\nS\t3\t" +
                               std::string(60, 'G') + "\nS\t4\t" + std::string(400, 'T') + "\nL\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\n"
                               "L\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\nP\ta#1#c:0-1\t1+,2+,4+\t*\nP\tb#1#c:0-1\t1+,3+,4+\t*\n");
    std::string r1, r2;
    for (int i = 0; i < 4; ++i) {
        r1 += "@p" + std::to_string(i) + "/1\n" + std::string(100, 'A') + "\n+\n" + std::string(100, 'I') + "\n";
        r2 += "@p" + std::to_string(i) + "/2\n" + std::string(100, 'A') + "\n+\n" + std::string(100, 'I') + "\n";
    }
    write(path("pre_r1.fq"), r1);
    write(path("pre_r2.fq"), r2);
    write(path("pre_library.txt"), kLibrary);
    write(path("pre_curation.tsv"), "no_such_path\n");
    const int exit_status = panvar::run_genotype_command({"--graph", path("pre.gfa"), "--reference-path", "a#1#c:0-1", "--reads1",
                                                          path("pre_r1.fq"), "--reads2", path("pre_r2.fq"), "--library",
                                                          path("pre_library.txt"), "--panel-curation", path("pre_curation.tsv"),
                                                          "--out", path("pre"), "--quiet"});
    require(exit_status == 1, "a pre-block failure must exit 1 (no-call), got " + std::to_string(exit_status));
    const auto row = locus_row(path("pre.locus_genotype.tsv"));
    require(row.at("run_status") == "MODEL_INCOMPLETE", "the fixture must fail before the blocks: " + row.at("run_status"));
    require(row.at("call_summary") == "GENOTYPE_NO_CALL;DOSAGE_NO_CALL;PHASE_NO_CALL",
            "a pre-block no-call must read NO_CALL in all three: " + row.at("call_summary"));
    require(row.at("genotype_filter_reasons") == "RUN_NO_CALL:status_MODEL_INCOMPLETE@locus" &&
                row.at("dosage_filter_reasons") == "RUN_NO_CALL@locus" && row.at("phase_filter_reasons") == "RUN_NO_CALL@locus",
            "the reasons must name the run's no-call: " + row.at("genotype_filter_reasons"));
    for (const char* f : {"genotype", "dosage", "phase"}) {
        const std::string p = std::string(f) + "_filter_";
        require(row.at(p + "units") == "1" && row.at(p + "applicable") == "1" && row.at(p + "no_call") == "1",
                std::string("one applicable NO_CALL unit for ") + f);
    }
    std::cout << "  pre-block MODEL_INCOMPLETE: exit 1, GENOTYPE_NO_CALL;DOSAGE_NO_CALL;PHASE_NO_CALL at the locus unit\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--print-inventory") {   // M15-H: the guide's block, regenerated
        std::cout << panvar::genotype_option_inventory_markdown();
        return 0;
    }
    if (argc < 2) {
        std::cerr << "usage: genotype_cli_contract_test <scratch dir> [source dir] | --print-inventory\n";
        return 2;
    }
    dir = argv[1];
    std::filesystem::create_directories(dir);
    try {
        test_library_file();
        test_repair_option_order();
        test_removed_options();
        test_help_text();
        if (argc >= 3) source_dir_for_profiles = argv[2];
        if (argc >= 3) test_guide_inventory(argv[2]);
        test_output_schema();
        test_preset_configuration();
        test_dosage_check_contract();
        test_panel_index_flags();
        test_pre_block_no_call();
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
    std::cout << "genotype CLI contract: all gates passed\n";
    return 0;
}

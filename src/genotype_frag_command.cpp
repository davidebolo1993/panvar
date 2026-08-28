#include "panvar/genotype_frag_command.hpp"

#include "panvar/align.hpp"
#include "panvar/bubble_path.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_fragments.hpp"
#include "panvar/gfa.hpp"
#include "panvar/output.hpp"
#include "panvar/parallel.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace panvar {

namespace {

void print_help() {
    std::cout
        << "Usage:\n"
        << "  panvar genotype-frag -i <graph.gfa> -b <bubble-prefix> -o <out-prefix>\n"
        << "                       -R <reads.fq.gz> [-R <reads2.fq.gz>]\n"
        << "\n"
        << "PROTOTYPE. Arm F4 of experiments/FRAGMENT_EVIDENCE_PREREGISTRATION.md, built to be run\n"
        << "beside `panvar genotype` rather than instead of it. No production default moves.\n"
        << "\n"
        << "The difference under test is the UNIT OF OBSERVATION. `genotype` reduces the reads to a\n"
        << "vector of marker counts and then reconstructs the lost read linkage through confinement,\n"
        << "over-expected filtering, clump correction and an ESS discount. This scores each PHYSICAL\n"
        << "FRAGMENT once, against the spelled sequence of each candidate allele plus its flanks:\n"
        << "\n"
        << "  log L(a,b) = SUM over fragments log( 0.5 P(frag|a) + 0.5 P(frag|b) )\n"
        << "\n"
        << "P(frag|allele) combines both mates' infix alignment edit distance and, when both mates\n"
        << "place, the implied insert length. Two consequences that the marker model cannot reach:\n"
        << "a marker shared between two blocks is localised by the block-specific marker on the SAME\n"
        << "fragment rather than deleted by a filter, and an allele shorter than k -- which carries\n"
        << "no k-mer at all and is invisible to the emission -- is still voted on, because the\n"
        << "fragment spans it and both its boundaries.\n"
        << "\n"
        << "NOT IMPLEMENTED here, by design: no block chain, no HMM, no linkage between blocks, no\n"
        << "depth/dosage channel. Each block is scored on its own fragments. Those layers are only\n"
        << "worth building if this one passes its gate.\n"
        << "\n"
        << "Options:\n"
        << "  -i, --gfa <path>            Input GFA (required)\n"
        << "  -b, --bubble-prefix-in <p>  Bubble prefix; reads <p>.bubbles.csv\n"
        << "  -c, --bubbles-csv-in <path> Bubbles CSV (alternative to --bubble-prefix-in)\n"
        << "  -o, --out-prefix <path>     Output prefix (required)\n"
        << "  -R, --reads <path>          FASTA/FASTQ, plain or gzipped; repeatable. Mates are joined\n"
        << "                              by name, so interleaved and split R1/R2 both work\n"
        << "      --haplotype-mode        Score whole panel HAPLOTYPE pairs end to end and project\n"
        << "                              the answer onto blocks, instead of scoring each block on\n"
        << "                              its own candidates. Measured reason it exists, at cyp2d6\n"
        << "                              block 5 under leave-ZERO-out: the panel's allele 7 (556 bp)\n"
        << "                              contains the truth label's allele 0 (261 bp) and both begin\n"
        << "                              at the same position in the sample's own haplotype -- the\n"
        << "                              rest of the sample's real sequence lives in the NEIGHBOURING\n"
        << "                              block. Any block-local context must guess a flank there, the\n"
        << "                              sample's flank is unknown, and a fragment spanning the\n"
        << "                              junction then prefers whichever candidate supplies more real\n"
        << "                              sequence before the guess begins. The evidence that decides\n"
        << "                              the block sits across the block boundary, so the block is\n"
        << "                              the wrong unit to score. In this mode nothing is guessed:\n"
        << "                              every haplotype carries its own flank\n"
        << "                              LIMIT: a pair of complete panel haplotypes cannot represent\n"
        << "                              a mosaic, so under leave-one-out this is bounded by the\n"
        << "                              panel's mosaic ceiling. Measuring that bound is the point\n"
        << "      --max-haplotypes <N>    Shortlist size in --haplotype-mode (default 48)\n"
        << "      --probe-haplotypes <a,b>  Report the likelihood's opinion of this named haplotype\n"
        << "                              pair: whether it was shortlisted at all, its rank, its\n"
        << "                              score delta to the best pair, and how many fragments placed\n"
        << "                              on each. Separates three failures a truncated top-N list\n"
        << "                              cannot -- never shortlisted, shortlisted and scored badly,\n"
        << "                              or scored and close. Repeatable; costs nothing, the pair\n"
        << "                              scores already exist\n"
        << "      --max-anchor-occ <N>    Occurrences PER HAPLOTYPE above which a syncmer anchors\n"
        << "                              nothing (default 8): inside a tandem array it points\n"
        << "                              everywhere. Per haplotype, not across the shortlist -- a\n"
        << "                              shared cap makes anchoring depend on shortlist size, and\n"
        << "                              it did: at cyp2d6 NA18939 the same reads put the truth at\n"
        << "                              rank 2 with 48 haplotypes and rank 1 with 96\n"
        << "      --anchor-slack <N>      Window either side of an anchored read start (default 40)\n"
        << "      --project-marginal      Call each block from the allele pair with the most\n"
        << "                              posterior mass summed over haplotype pairs, instead of from\n"
        << "                              the best pair. It sounds stricter and measures worse: at\n"
        << "                              cyp2d6 under leave-one-out the best pair reconstructs the\n"
        << "                              donor to a median 130 edits and the marginal projection of\n"
        << "                              the SAME posterior to 385, because a per-block argmax of\n"
        << "                              marginals assembles a combination no single pair realises.\n"
        << "                              Kept so that stays reproducible\n"
        << "      --marginalise-placements  Sum a fragment's likelihood over EVERY placement it has\n"
        << "                              on a haplotype instead of taking its best one. A fragment\n"
        << "                              compatible with ten positions is evidence for a haplotype\n"
        << "                              offering ten, and a maximum scores that identically to a\n"
        << "                              haplotype offering one -- which is worst exactly where\n"
        << "                              cyp2d6 is hardest. Implies --length-normalize, because a\n"
        << "                              sum without the 1/N factor rewards a repetitive haplotype\n"
        << "                              for offering more places to land\n"
        << "      --placement-topk <N>    Implied-start bins kept per mate per haplotype (default 2).\n"
        << "                              A sum over placements is only meaningful if the placements\n"
        << "                              are enumerated, so raise it with the flag above\n"
        << "      --total-depth           Score the pair's TOTAL length against how many fragments\n"
        << "                              were seen: N ~ Poisson(lambda * (L(h1)+L(h2))), with lambda\n"
        << "                              fitted once from the observed count over twice the panel's\n"
        << "                              median haplotype length -- outside any candidate, so the\n"
        << "                              term cannot be self-fulfilling. Each fragment is counted\n"
        << "                              once by construction, because the observation is a single\n"
        << "                              scalar; that is the defect the per-haplotype window channel\n"
        << "                              could not escape, where one fragment made every candidate\n"
        << "                              copy it aligned to look covered\n"
        << "      --haploid-depth <x>     Supply lambda (fragments per bp per haplotype copy) instead\n"
        << "                              of estimating it; implies --total-depth\n"
        << "      --truth-total-bp <N>    DIAGNOSTIC: rank pairs by how close their total length is to\n"
        << "                              N, the sample's true diploid total, which no caller can\n"
        << "                              know. Bounds what perfect dosage knowledge could buy -- if\n"
        << "                              the tail does not close under this, it is not a dosage\n"
        << "                              problem and no depth model will fix it\n"
        << "      --coverage-weight <w>   Weight on the depth channel in --haplotype-mode (default\n"
        << "                              1.0; 0 disables). Sequence compatibility and copy number are\n"
        << "                              different signals: reads from a duplicated segment align\n"
        << "                              perfectly to a single copy, so alignment alone cannot see\n"
        << "                              dosage. This scores each haplotype's own length in windows,\n"
        << "                              with the rate fitted per haplotype. It reuses the same\n"
        << "                              fragments as the sequence term, so at weight 1 the two\n"
        << "                              double count -- which is why the default is 0, not 1:\n"
        << "                              measured at cyp2d6 leave-ZERO-out, weight 1 moved NA18939's\n"
        << "                              own true pair from rank 1 to rank 2 and cost 8 of its 19\n"
        << "                              blocks. Real mechanism, not yet a safe default\n"
        << "      --coverage-window <N>   Window size for that channel (default 500)\n"
        << "      --blocks <a,b,c>        Score only these chain indices. Default: every bubble block\n"
        << "      --exact-distance <a.fa> <b.fa>\n"
        << "                              Global (Needleman-Wunsch) edit distance between two\n"
        << "                              sequences, printed and nothing else. The panel floor is\n"
        << "                              computed from minimap2, which reports alignments rather\n"
        << "                              than a distance: split or overlapping records can double\n"
        << "                              count edits, and anything it declines to align is charged\n"
        << "                              as a whole unaligned base. Use this to check the floor\n"
        << "                              exactly on the few nearest candidates -- it is affordable\n"
        << "                              there and not over a whole panel\n"
        << "      --mosaic-floor          Three floors instead of one, needing --truth-haplotypes\n"
        << "                              and no reads: COMPLETE (one panel haplotype per homologue\n"
        << "                              across the locus), FREE (the nearest panel allele at every\n"
        << "                              block, chosen independently -- a bound, since it may switch\n"
        << "                              source at every boundary with no evidence for any switch),\n"
        << "                              and PENALISED (the same with a cost per switch). Their\n"
        << "                              difference is the most a mosaic model could ever recover,\n"
        << "                              which is the measurement that decides whether a factor\n"
        << "                              graph is worth building at all\n"
        << "      --switch-penalties <a,b,c>  Penalties to evaluate (default 0,10,100,1000)\n"
        << "      --spell-calls <gt.tsv>  Read a per-block call table (either caller's) and write the\n"
        << "                              two sequences it claims, to <out-prefix>.called.fa. Needed\n"
        << "                              because the two callers report in different shapes and\n"
        << "                              cannot otherwise be compared on the objective that matters:\n"
        << "                              production emits a per-block allele pair, which is a mosaic\n"
        << "                              and has no called haplotype to align. Needs no reads\n"
        << "      --all-blocks            Score backbone and flank blocks too. gstm1's worst blocks\n"
        << "                              are BACKBONE blocks carrying 16 and 21 alleles, so the\n"
        << "                              default view omits them\n"
        << "      --truth-haplotypes <a,b>  Two panel haplotype names the reads came from. Adds the\n"
        << "                              truth pair's rank, tie count and delta per block\n"
        << "      --exclude-haplotypes <a,b>  Drop these from the panel first. With\n"
        << "                              --truth-haplotypes this is leave-one-out; the truth allele\n"
        << "                              is then resolved by spelling the held-out walk and matching\n"
        << "                              it against the reduced panel, exactly as `genotype` does\n"
        << "      --flank-bp <N>          Neighbouring sequence glued to each side of a candidate\n"
        << "                              (default 500). Identical for every candidate of a block, so\n"
        << "                              it buys reachability for boundary-spanning fragments and\n"
        << "                              cannot shift a comparison\n"
        << "      --max-alleles <N>       Candidates kept per block after coarse containment\n"
        << "                              shortlisting (default 64). truth_rank -2 means the truth\n"
        << "                              was representable but dropped here\n"
        << "      --min-hits <N>          Syncmer hits needed to recruit a fragment to a block\n"
        << "                              (default 2)\n"
        << "      --error-rate <p>        Per-base edit probability (default 0.01)\n"
        << "      --bg-divergence <d>     Divergence at which a fragment stops being able to prefer\n"
        << "                              anything (default 0.10). Alignment likelihood is unbounded\n"
        << "                              below, so without this one fragment from sequence no\n"
        << "                              candidate models outvotes hundreds of ordinary ones -- and\n"
        << "                              under leave-one-out every candidate lacks some of the\n"
        << "                              sample's sequence, so that asymmetry has a direction\n"
        << "      --frag-outlier <e>      Prior weight on that background component (default 1e-3)\n"
        << "      --fragment-len <N>      Library mean insert (default 350)\n"
        << "      --fragment-sd <N>       Its standard deviation (default 50)\n"
        << "      --no-insert-size        Score sequence compatibility only. The channel separation\n"
        << "                              matters: run both to say which one carries a result\n"
        << "      --no-length-normalize   Drop the 1/L factor. Without it P(fragment|allele) is not\n"
        << "                              a generative likelihood and a SUPERSET allele is never\n"
        << "                              penalised -- every read from a short allele also fits an\n"
        << "                              allele containing it plus extra, and the extra is free.\n"
        << "                              Measured: leaving it off is what made leave-zero-out fail\n"
        << "      --no-compete            Score the block in isolation instead of against the rest\n"
        << "                              of the locus. A block scored alone must explain every\n"
        << "                              fragment handed to it, so a read from a paralogous copy\n"
        << "                              elsewhere votes here -- at CYP2D6 for whichever candidate\n"
        << "                              is most paralog-like\n"
        << "      --max-divergence <d>    Alignment band as a fraction of read length (default 0.20)\n"
        << "      --top-pairs <N>         Pairs written per block to <prefix>.frag_pairs.tsv\n"
        << "                              (default 20)\n"
        << "      --debug-block <N>       Write the per-fragment, per-candidate log-likelihood and\n"
        << "                              edit distance for chain block N to\n"
        << "                              <prefix>.frag_debug.tsv. A block score is a sum over\n"
        << "                              hundreds of fragments and nothing about it is checkable\n"
        << "                              from the sum\n"
        << "  -k, --kmer-size <N>         k for recruitment syncmers (default 31)\n"
        << "      --syncmer-s <N>         s for the closed-syncmer test (0 = auto)\n"
        << "  -t, --threads <N>           Worker threads (0 = auto)\n"
        << "  -q, --quiet                 Disable progress logs\n"
        << "  -h, --help                  Show this help\n";
}

std::vector<std::string> split_commas(const std::string& text) {
    std::vector<std::string> out;
    std::string cur;
    for (const char ch : text) {
        if (ch == ',') { if (!cur.empty()) out.push_back(cur); cur.clear(); }
        else cur.push_back(ch);
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

} // namespace

int run_genotype_frag_command(const std::vector<std::string>& args) {
    if (args.empty()) { print_help(); return 0; }

    std::string gfa_path, bubble_prefix_in, bubbles_csv_in, out_prefix;
    std::vector<std::string> read_paths;
    std::string truth_haplotypes, exclude_haplotypes, blocks_arg, spell_calls;
    bool mosaic_floor = false;
    std::string switch_penalties_arg = "0,10,100,1000";
    std::vector<std::string> exact_distance;
    bool all_blocks = false, quiet = false, hap_mode = false, length_normalize_set = false;
    std::size_t top_pairs = 20;
    HaplotypeScoreOptions hopt;
    FragmentScoreOptions& opt = hopt;

    const auto value = [&](std::size_t& i, const std::string& name) {
        if (i + 1 >= args.size()) throw std::runtime_error("genotype-frag: " + name + " needs a value");
        return args[++i];
    };

    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& a = args[i];
        if (a == "-h" || a == "--help") { print_help(); return 0; }
        else if (a == "-i" || a == "--gfa") gfa_path = value(i, a);
        else if (a == "-b" || a == "--bubble-prefix-in") bubble_prefix_in = value(i, a);
        else if (a == "-c" || a == "--bubbles-csv-in") bubbles_csv_in = value(i, a);
        else if (a == "-o" || a == "--out-prefix") out_prefix = value(i, a);
        else if (a == "-R" || a == "--reads") read_paths.push_back(value(i, a));
        else if (a == "--blocks") blocks_arg = value(i, a);
        else if (a == "--all-blocks") all_blocks = true;
        else if (a == "--spell-calls") spell_calls = value(i, a);
        else if (a == "--mosaic-floor") mosaic_floor = true;
        else if (a == "--switch-penalties") switch_penalties_arg = value(i, a);
        else if (a == "--exact-distance") { exact_distance.push_back(value(i, a));
                                            exact_distance.push_back(value(i, a)); }
        else if (a == "--truth-haplotypes") truth_haplotypes = value(i, a);
        else if (a == "--exclude-haplotypes") exclude_haplotypes = value(i, a);
        else if (a == "--flank-bp") opt.flank_bp = cli::parse_size_arg(a, value(i, a));
        else if (a == "--max-alleles") opt.max_alleles = cli::parse_size_arg(a, value(i, a));
        else if (a == "--min-hits") opt.min_recruit_hits = cli::parse_size_arg(a, value(i, a));
        else if (a == "--error-rate") opt.error_rate = std::stod(value(i, a));
        else if (a == "--bg-divergence") opt.bg_divergence = std::stod(value(i, a));
        else if (a == "--frag-outlier") opt.outlier_mix = std::stod(value(i, a));
        else if (a == "--fragment-len") opt.fragment_len = std::stod(value(i, a));
        else if (a == "--fragment-sd") opt.fragment_sd = std::stod(value(i, a));
        else if (a == "--no-insert-size") opt.use_insert_size = false;
        else if (a == "--no-length-normalize") { opt.length_normalize = false; length_normalize_set = true; }
        else if (a == "--length-normalize") { opt.length_normalize = true; length_normalize_set = true; }
        else if (a == "--no-compete") opt.compete = false;
        else if (a == "--haplotype-mode") hap_mode = true;
        else if (a == "--max-haplotypes") hopt.max_haplotypes = cli::parse_size_arg(a, value(i, a));
        else if (a == "--max-anchor-occ") hopt.max_anchor_occ = cli::parse_size_arg(a, value(i, a));
        else if (a == "--anchor-slack") hopt.anchor_slack = cli::parse_size_arg(a, value(i, a));
        else if (a == "--probe-haplotypes") {
            const std::vector<std::string> two = split_commas(value(i, a));
            if (two.size() != 2) {
                throw std::runtime_error("genotype-frag: --probe-haplotypes needs exactly two names");
            }
            hopt.probe_pairs.emplace_back(two[0], two[1]);
        }
        else if (a == "--project-marginal") hopt.project_map = false;
        else if (a == "--marginalise-placements") hopt.marginalise_placements = true;
        else if (a == "--placement-topk") hopt.placement_topk = cli::parse_size_arg(a, value(i, a));
        else if (a == "--joint-depth") hopt.joint_depth = true;
        else if (a == "--joint-marginal") { hopt.joint_depth = true; hopt.joint_marginal = true; }
        else if (a == "--joint-reverse-order") hopt.joint_reverse_order = true;
        else if (a == "--equivalence-tolerance") hopt.equivalence_tolerance = std::stod(value(i, a));
        else if (a == "--joint-top-pairs") hopt.joint_top_pairs = cli::parse_size_arg(a, value(i, a));
        else if (a == "--joint-window") hopt.joint_window = cli::parse_size_arg(a, value(i, a));
        else if (a == "--joint-iterations") hopt.joint_iterations = cli::parse_size_arg(a, value(i, a));
        else if (a == "--total-depth") hopt.total_depth = true;
        else if (a == "--haploid-depth") { hopt.haploid_depth = std::stod(value(i, a)); hopt.total_depth = true; }
        else if (a == "--truth-total-bp") hopt.truth_total_bp = std::stod(value(i, a));
        else if (a == "--coverage-weight") hopt.coverage_weight = std::stod(value(i, a));
        else if (a == "--coverage-window") hopt.coverage_window = cli::parse_size_arg(a, value(i, a));
        else if (a == "--max-divergence") opt.max_divergence = std::stod(value(i, a));
        else if (a == "--top-pairs") top_pairs = cli::parse_size_arg(a, value(i, a));
        else if (a == "--debug-block") opt.debug_block =
            static_cast<long>(cli::parse_size_arg(a, value(i, a)));
        else if (a == "-k" || a == "--kmer-size") opt.kmer_size = cli::parse_size_arg(a, value(i, a));
        else if (a == "--syncmer-s") opt.syncmer_s = cli::parse_size_arg(a, value(i, a));
        else if (a == "-t" || a == "--threads") opt.threads = cli::parse_size_arg(a, value(i, a));
        else if (a == "-q" || a == "--quiet") quiet = true;
        else throw std::runtime_error("genotype-frag: unknown option " + a);
    }

    if (!exact_distance.empty()) {
        const auto slurp = [](const std::string& path) {
            std::ifstream in(path);
            if (!in) throw std::runtime_error("genotype-frag: cannot read " + path);
            std::string line, seq;
            while (std::getline(in, line)) {
                if (!line.empty() && line[0] == '>') continue;
                while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
                seq += line;
            }
            return seq;
        };
        const std::string a = slurp(exact_distance[0]);
        const std::string b = slurp(exact_distance[1]);
        if (a.empty() || b.empty()) { std::cout << (a.size() + b.size()) << '\n'; return 0; }
        // Progressive banding. An unbanded global alignment of two 200 kb haplotypes is quadratic in
        // practice and does not finish in useful time, but the distances this is used on are small --
        // a panel floor is tens to a few thousand edits over 200 kb. edlib returns the EXACT distance
        // whenever it is within the band, so doubling from a small band is exact and fast in the
        // common case and only falls back to the full computation when the sequences really are far
        // apart. The band actually used is reported on stderr so a number can never be mistaken for
        // an unbanded one.
        for (const std::size_t band : {std::size_t{1024}, std::size_t{4096}, std::size_t{16384},
                                      std::size_t{65536}, std::size_t{262144}}) {
            const NwBanded r = nw_edit_distance_banded(a, b, band);
            if (r.ok) {
                std::fprintf(stderr, "band %zu sufficed\n", band);
                std::cout << r.edits << '\n';
                return 0;
            }
        }
        std::fprintf(stderr, "no band sufficed; computing unbanded\n");
        std::cout << nw_edit_distance(a, b).edits << '\n';
        return 0;
    }

    if (gfa_path.empty() || out_prefix.empty()) {
        throw std::runtime_error("genotype-frag requires --gfa and --out-prefix");
    }
    if (read_paths.empty() && spell_calls.empty() && !mosaic_floor) {
        throw std::runtime_error("genotype-frag requires at least one --reads");
    }
    if (!bubble_prefix_in.empty()) {
        if (!bubbles_csv_in.empty()) {
            throw std::runtime_error("genotype-frag: use either --bubble-prefix-in or --bubbles-csv-in");
        }
        bubbles_csv_in = bubble_prefix_in + ".bubbles.csv";
    }
    if (bubbles_csv_in.empty()) {
        throw std::runtime_error("genotype-frag requires --bubble-prefix-in or --bubbles-csv-in");
    }
    if (opt.error_rate <= 0.0 || opt.error_rate >= 1.0) {
        throw std::runtime_error("genotype-frag: --error-rate must be in (0,1)");
    }
    if (opt.bg_divergence <= 0.0 || opt.bg_divergence >= 1.0) {
        throw std::runtime_error("genotype-frag: --bg-divergence must be in (0,1)");
    }
    if (opt.outlier_mix <= 0.0 || opt.outlier_mix >= 1.0) {
        throw std::runtime_error("genotype-frag: --frag-outlier must be in (0,1)");
    }
    if (opt.fragment_sd <= 0.0) throw std::runtime_error("genotype-frag: --fragment-sd must be > 0");
    if (opt.kmer_size == 0 || opt.kmer_size > 31) {
        throw std::runtime_error("genotype-frag: --kmer-size must be in [1,31]");
    }
    if (opt.debug_block >= 0) opt.debug_path = out_prefix + ".frag_debug.tsv";
    cli::ensure_parent_dir_for_file(out_prefix);

    cli::RunLog log("genotype-frag", quiet);

    ParseGfaOptions parse_options;
    parse_options.include_paths = true;
    parse_options.include_sequences = true;
    const Graph graph = parse_gfa(gfa_path, parse_options);
    if (graph.paths.empty()) throw std::runtime_error("genotype-frag: no paths in " + gfa_path);

    // Leave-one-out, with the same contract as `genotype`: the named haplotypes leave the panel but
    // their walks are kept, so the truth allele can still be resolved against the REDUCED panel. A
    // block where no remaining haplotype spells the held-out sequence is unrepresentable -- the
    // mosaic ceiling -- and is reported as truth_rank -1 rather than counted as a wrong call.
    Graph panel_graph = graph;
    std::vector<PathRecord> held_out;
    if (!exclude_haplotypes.empty()) {
        const std::vector<std::string> names = split_commas(exclude_haplotypes);
        std::vector<PathRecord> keep;
        for (const PathRecord& p : panel_graph.paths) {
            if (std::find(names.begin(), names.end(), p.name) != names.end()) held_out.push_back(p);
            else keep.push_back(p);
        }
        if (held_out.size() != names.size()) {
            throw std::runtime_error("genotype-frag: --exclude-haplotypes named a path not in the graph");
        }
        panel_graph.paths = std::move(keep);
    }

    const std::vector<Bubble> bubbles = read_bubbles_csv(bubbles_csv_in);
    const std::vector<Block> chain = build_block_chain(bubbles);
    std::vector<BubblePathIndex> path_indexes(panel_graph.paths.size());
    run_parallel(panel_graph.paths.size(), opt.threads, [&](std::size_t i) {
        path_indexes[i] = build_bubble_path_index(panel_graph.paths[i]);
    });
    std::vector<BlockAlleles> blocks(chain.size());
    for (std::size_t i = 0; i < chain.size(); ++i) {
        blocks[i] = enumerate_block_alleles(panel_graph, path_indexes, bubbles, chain[i], opt.threads);
    }
    log.info("input " + gfa_path + " (" + std::to_string(graph.nodes.size()) + " nodes, " +
             std::to_string(graph.paths.size()) + " paths); " + std::to_string(bubbles.size()) +
             " bubbles; chain of " + std::to_string(chain.size()) + " blocks" +
             (held_out.empty() ? "" : "; " + std::to_string(held_out.size()) + " held out"));

    if (!spell_calls.empty()) {
        std::ifstream cf(spell_calls);
        if (!cf) throw std::runtime_error("genotype-frag: cannot read " + spell_calls);
        std::string line;
        if (!std::getline(cf, line)) throw std::runtime_error("genotype-frag: empty " + spell_calls);
        std::vector<std::string> header;
        {
            std::size_t start = 0;
            for (std::size_t i = 0; i <= line.size(); ++i) {
                if (i == line.size() || line[i] == '\t') {
                    header.push_back(line.substr(start, i - start));
                    start = i + 1;
                }
            }
        }
        const auto col = [&](const std::string& n) {
            for (std::size_t i = 0; i < header.size(); ++i) if (header[i] == n) return static_cast<long>(i);
            throw std::runtime_error("genotype-frag: --spell-calls table has no column '" + n + "'");
        };
        // Both callers' tables are accepted: production names the column block_index, the prototype
        // names it block. Requiring one of them would make the comparison depend on which caller
        // happened to write the file.
        const long c_blk = std::find(header.begin(), header.end(), std::string("block_index")) != header.end()
            ? col("block_index") : col("block");
        const long c_a1 = col("allele1");
        const long c_a2 = col("allele2");
        std::vector<int> a1(chain.size(), -1), a2(chain.size(), -1);
        while (std::getline(cf, line)) {
            std::vector<std::string> f;
            std::size_t start = 0;
            for (std::size_t i = 0; i <= line.size(); ++i) {
                if (i == line.size() || line[i] == '\t') { f.push_back(line.substr(start, i - start)); start = i + 1; }
            }
            if (f.size() <= static_cast<std::size_t>(std::max(c_blk, std::max(c_a1, c_a2)))) continue;
            const long bi = std::stol(f[static_cast<std::size_t>(c_blk)]);
            if (bi < 0 || static_cast<std::size_t>(bi) >= chain.size()) continue;
            a1[static_cast<std::size_t>(bi)] = std::stoi(f[static_cast<std::size_t>(c_a1)]);
            a2[static_cast<std::size_t>(bi)] = std::stoi(f[static_cast<std::size_t>(c_a2)]);
        }
        std::string s1, s2;
        spell_called_pair(blocks, a1, a2, s1, s2);
        const std::string fa = out_prefix + ".called.fa";
        std::ofstream of(fa);
        if (!of) throw std::runtime_error("genotype-frag: cannot write " + fa);
        of << ">called_1\n" << s1 << "\n>called_2\n" << s2 << '\n';
        of.flush();
        if (!of) throw std::runtime_error("genotype-frag: write failed for " + fa);
        log.info("spelled " + spell_calls + ": " + std::to_string(s1.size()) + " bp and " +
                 std::to_string(s2.size()) + " bp");
        log.wrote({fa});
        log.done();
        return 0;
    }

    // ---- truth, resolved the same way `genotype` resolves it ---------------------------------
    std::vector<int> truth1, truth2;
    std::vector<std::string> truth_seq1, truth_seq2;
    bool have_truth = false;
    if (!truth_haplotypes.empty()) {
        const std::vector<std::string> names = split_commas(truth_haplotypes);
        if (names.size() != 2) {
            throw std::runtime_error("genotype-frag: --truth-haplotypes needs exactly two names");
        }
        std::vector<BlockAlleles> held_blocks(chain.size());
        if (!held_out.empty()) {
            Graph held_graph = graph;
            held_graph.paths = held_out;
            std::vector<BubblePathIndex> held_idx(held_graph.paths.size());
            for (std::size_t kk = 0; kk < held_graph.paths.size(); ++kk) {
                held_idx[kk] = build_bubble_path_index(held_graph.paths[kk]);
            }
            for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                held_blocks[bi] = enumerate_block_alleles(held_graph, held_idx, bubbles,
                                                          chain[bi], opt.threads);
            }
        }
        const auto resolve = [&](const std::string& name, std::vector<int>& out_alleles,
                                 std::vector<std::string>* out_seq = nullptr) {
            out_alleles.assign(chain.size(), -1);
            if (out_seq != nullptr) out_seq->assign(chain.size(), std::string());
            const bool is_held =
                std::any_of(held_out.begin(), held_out.end(),
                            [&](const PathRecord& p) { return p.name == name; });
            for (std::size_t bi = 0; bi < chain.size(); ++bi) {
                if (!is_held) {
                    const auto it = blocks[bi].allele_of.find(name);
                    if (it != blocks[bi].allele_of.end()) {
                        out_alleles[bi] = static_cast<int>(it->second);
                        if (out_seq != nullptr && it->second < blocks[bi].allele_seq.size()) {
                            (*out_seq)[bi] = blocks[bi].allele_seq[it->second];
                        }
                    }
                    continue;
                }
                const auto it = held_blocks[bi].allele_of.find(name);
                if (it == held_blocks[bi].allele_of.end()) continue;
                const std::size_t hai = it->second;
                if (held_blocks[bi].bypass_allele >= 0 &&
                    hai == static_cast<std::size_t>(held_blocks[bi].bypass_allele)) {
                    if (blocks[bi].bypass_allele >= 0) out_alleles[bi] = blocks[bi].bypass_allele;
                    continue;
                }
                if (hai >= held_blocks[bi].allele_seq.size()) continue;
                const std::string& seq = held_blocks[bi].allele_seq[hai];
                if (out_seq != nullptr) (*out_seq)[bi] = seq;
                if (seq.empty()) continue;
                for (std::size_t ai = 0; ai < blocks[bi].allele_seq.size(); ++ai) {
                    if (blocks[bi].allele_seq[ai] == seq) {
                        out_alleles[bi] = static_cast<int>(ai);
                        break;
                    }
                }
            }
        };
        resolve(names[0], truth1, &truth_seq1);
        resolve(names[1], truth2, &truth_seq2);
        have_truth = true;
        std::size_t representable = 0;
        for (std::size_t bi = 0; bi < chain.size(); ++bi) {
            if (truth1[bi] >= 0 && truth2[bi] >= 0) ++representable;
        }
        log.info("truth " + names[0] + " / " + names[1] + ": both alleles representable at " +
                 std::to_string(representable) + "/" + std::to_string(chain.size()) + " blocks");
    }

    if (mosaic_floor) {
        if (!have_truth) throw std::runtime_error("genotype-frag: --mosaic-floor needs --truth-haplotypes");
        std::vector<double> pens;
        for (const std::string& t : split_commas(switch_penalties_arg)) pens.push_back(std::stod(t));
        std::vector<std::string> names;
        for (const PathRecord& p : panel_graph.paths) names.push_back(p.name);
        const std::string fp = out_prefix + ".mosaic_floor.tsv";
        std::ofstream mf(fp);
        if (!mf) throw std::runtime_error("genotype-frag: cannot write " + fp);
        mf << "homologue\tfloor\tswitch_penalty\tedits\tswitches\tblocks_scored\tblocks\n";
        std::size_t tot_complete = 0, tot_free = 0;
        std::vector<std::size_t> tot_pen(pens.size(), 0);
        for (int side = 0; side < 2; ++side) {
            const MosaicFloors r = mosaic_floors(blocks, names, side == 0 ? truth_seq1 : truth_seq2,
                                                 pens, opt.threads);
            const int hn = side + 1;
            mf << hn << "\tcomplete\tNA\t" << r.complete << "\tNA\t" << r.blocks_scored << '\t'
               << r.blocks << '\n';
            mf << hn << "\tfree\t0\t" << r.free_mosaic << '\t' << r.switches_free << '\t'
               << r.blocks_scored << '\t' << r.blocks << '\n';
            for (std::size_t i = 0; i < r.penalised.size(); ++i) {
                mf << hn << "\tpenalised\t" << r.penalised[i].first << '\t' << r.penalised[i].second
                   << '\t' << r.penalised_switches[i] << '\t' << r.blocks_scored << '\t'
                   << r.blocks << '\n';
                tot_pen[i] += r.penalised[i].second;
            }
            tot_complete += r.complete;
            tot_free += r.free_mosaic;
        }
        mf.flush();
        if (!mf) throw std::runtime_error("genotype-frag: write failed for " + fp);
        log.info("floors, both homologues: complete " + std::to_string(tot_complete) +
                 ", free mosaic " + std::to_string(tot_free) + " (H_mosaic " +
                 std::to_string(tot_complete > tot_free ? tot_complete - tot_free : 0) + ")");
        for (std::size_t i = 0; i < pens.size(); ++i) {
            log.info("  switch penalty " + std::to_string(pens[i]) + ": " + std::to_string(tot_pen[i]));
        }
        log.wrote({fp});
        log.done();
        return 0;
    }

    if (hap_mode && hopt.joint_depth && hopt.haploid_depth <= 0.0) {
        throw std::runtime_error(
            "genotype-frag: --joint-depth requires --haploid-depth. The depth rate must come from "
            "outside the candidate set -- invariant flanks, or genome-wide depth, or a known "
            "simulation rate. Fitting it from the alignment-best pair lets the expectation follow "
            "the hypothesis it is meant to test, and that is what the first implementation did");
    }
    if (hap_mode) {
        // Summing over placements without dividing by the number of positions a haplotype offers
        // rewards a repetitive haplotype for offering more places to land. The two go together.
        if (!length_normalize_set) hopt.length_normalize = hopt.marginalise_placements;
        FragmentLoadStats hs;
        const std::vector<Fragment> frags = load_fragments(read_paths, &hs);
        log.info("reads: " + std::to_string(hs.reads) + " -> " + std::to_string(hs.fragments) +
                 " fragments (" + std::to_string(hs.paired) + " paired, " +
                 std::to_string(hs.singleton) + " single)");
        if (hs.over_paired > 0) {
            log.info("WARNING: " + std::to_string(hs.over_paired) + " reads shared a name with an "
                     "already-complete pair and were dropped; fragment counts are unreliable");
        }
        if (hs.paired == 0 && hopt.use_insert_size) {
            hopt.use_insert_size = false;
            log.info("no paired fragments: insert-size channel off");
        }
        std::vector<std::string> names;
        for (const PathRecord& p : panel_graph.paths) names.push_back(p.name);
        log.info("scoring " + std::to_string(names.size()) + " panel haplotypes (shortlist " +
                 std::to_string(hopt.max_haplotypes) + ") over " + std::to_string(frags.size()) +
                 " fragments");
        const HaplotypeResult hr =
            genotype_haplotype_pairs(chain, blocks, names, frags, hopt,
                                     have_truth ? &truth1 : nullptr,
                                     have_truth ? &truth2 : nullptr, top_pairs);
        if (!hr.top_pairs.empty()) {
            log.info("best pair: " + hr.shortlist[hr.top_pairs[0].hap1] + " / " +
                     hr.shortlist[hr.top_pairs[0].hap2] + " (posterior " +
                     std::to_string(hr.top_pairs[0].posterior) + ")");
        }
        if (hr.convergence.pairs_rescored > 0 && !hopt.joint_marginal) {
            log.info("joint assignment: " + std::to_string(hr.convergence.pairs_converged) + "/" +
                     std::to_string(hr.convergence.pairs_rescored) + " pairs reached a fixed point; "
                     "max iterations used " + std::to_string(hr.convergence.max_iterations_used) +
                     "; fragments still moving in the last pass " +
                     std::to_string(hr.convergence.total_moves_last_iteration));
        }
        log.info(std::to_string(hr.n_informative) + " of " + std::to_string(hr.n_fragments) +
                 " fragments discriminate between shortlisted haplotypes");
        if (have_truth) {
            // Whether the coarse stage kept the answer is reported, never assumed: a candidate
            // generator that drops the truth loses it outright and no score below can recover it.
            const std::vector<std::string> tn = split_commas(truth_haplotypes);
            for (const std::string& n : tn) {
                const bool kept = std::find(hr.shortlist.begin(), hr.shortlist.end(), n) != hr.shortlist.end();
                const bool in_panel =
                    std::any_of(panel_graph.paths.begin(), panel_graph.paths.end(),
                                [&](const PathRecord& p) { return p.name == n; });
                if (in_panel && !kept) {
                    log.info("WARNING: truth haplotype " + n + " is in the panel but did NOT survive "
                             "the coarse shortlist; raise --max-haplotypes");
                }
            }
            std::size_t rep = 0, ex = 0, bub_rep = 0, bub_ex = 0;
            for (const BlockProjection& p : hr.blocks) {
                if (!p.truth_representable) continue;
                ++rep; if (p.exact) ++ex;
                if (p.kind == BlockKind::Bubble) { ++bub_rep; if (p.exact) ++bub_ex; }
            }
            log.info("blocks with a representable truth: " + std::to_string(rep) + "; exact " +
                     std::to_string(ex) + " (" + std::to_string(rep ? 100 * ex / rep : 0) +
                     "%); bubble blocks " + std::to_string(bub_ex) + "/" + std::to_string(bub_rep));
        }
        write_haplotype_results(out_prefix, hr, have_truth);
        {
            const auto& c = hr.completeness;
            const double kept_pct = c.clusters_found == 0 ? 100.0
                : 100.0 * static_cast<double>(c.clusters_kept) / static_cast<double>(c.clusters_found);
            const double anch_pct = c.anchor_occurrences_seen == 0 ? 0.0
                : 100.0 * static_cast<double>(c.anchor_occurrences_dropped) /
                  static_cast<double>(c.anchor_occurrences_seen);
            log.info("placement completeness: " + std::to_string(static_cast<int>(kept_pct)) +
                     "% of implied-start clusters kept, " + std::to_string(static_cast<int>(anch_pct)) +
                     "% of anchor occurrences dropped by --max-anchor-occ, " +
                     std::to_string(c.fragments_truncated) + " fragments hit --placement-topk. "
                     "Nothing may be called an information limit while these are far from complete");
        }
        log.info("call margin " + std::to_string(hr.equivalence.margin) + " nats; " +
                 std::to_string(hr.equivalence.size) + " pair(s) within tolerance; " +
                 std::to_string(hr.equivalence.blocks_determined) + "/" +
                 std::to_string(hr.equivalence.blocks_total) + " blocks determined by all of them");
        log.wrote({out_prefix + ".hap_blocks.tsv", out_prefix + ".hap_pairs.tsv",
                   out_prefix + ".hap_scores.tsv", out_prefix + ".equivalence.tsv"});
        log.done();
        return 0;
    }

    // ---- which blocks ------------------------------------------------------------------------
    std::vector<std::size_t> targets;
    if (!blocks_arg.empty()) {
        for (const std::string& tok : split_commas(blocks_arg)) {
            const std::size_t bi = cli::parse_size_arg("--blocks", tok);
            if (bi >= chain.size()) {
                throw std::runtime_error("genotype-frag: --blocks " + tok + " is past the chain (" +
                                         std::to_string(chain.size()) + " blocks)");
            }
            targets.push_back(bi);
        }
    } else {
        for (std::size_t bi = 0; bi < chain.size(); ++bi) {
            if (all_blocks || chain[bi].kind == BlockKind::Bubble) targets.push_back(bi);
        }
    }
    if (targets.empty()) throw std::runtime_error("genotype-frag: no blocks selected");

    FragmentLoadStats fstats;
    const std::vector<Fragment> fragments = load_fragments(read_paths, &fstats);
    log.info("reads: " + std::to_string(fstats.reads) + " -> " + std::to_string(fstats.fragments) +
             " fragments (" + std::to_string(fstats.paired) + " paired, " +
             std::to_string(fstats.singleton) + " single)");
    if (fstats.over_paired > 0) {
        log.info("WARNING: " + std::to_string(fstats.over_paired) + " reads shared a name with an "
                 "already-complete pair and were dropped. Fragment counts are unreliable; check that "
                 "the input is one library and that read names are not duplicated across files");
    }
    if (fstats.paired == 0 && opt.use_insert_size) {
        opt.use_insert_size = false;
        log.info("no paired fragments: the insert-size channel is off, sequence compatibility only");
    }

    log.info("scoring " + std::to_string(targets.size()) + " block(s) from " +
             std::to_string(fragments.size()) + " fragments");
    const std::vector<BlockFragmentResult> results =
        genotype_fragments(chain, blocks, fragments, targets, opt,
                           have_truth ? &truth1 : nullptr, have_truth ? &truth2 : nullptr,
                           top_pairs);

    std::size_t scored = 0, exact = 0, unrep = 0, pruned = 0, no_frag = 0;
    for (const BlockFragmentResult& r : results) {
        if (r.n_fragments == 0) ++no_frag;
        if (!have_truth) continue;
        if (r.truth_rank == -1) { ++unrep; continue; }
        if (r.truth_rank == -2) { ++pruned; continue; }
        ++scored;
        const std::size_t ta = std::min<std::size_t>(r.truth_a, r.truth_b);
        const std::size_t tb = std::max<std::size_t>(r.truth_a, r.truth_b);
        if (ta == std::min(r.best_a, r.best_b) && tb == std::max(r.best_a, r.best_b)) ++exact;
    }
    if (no_frag > 0) {
        log.info(std::to_string(no_frag) + " block(s) recruited no fragment at all -- their call is "
                 "arbitrary, not merely uncertain");
    }
    if (have_truth) {
        log.info("truth scored at " + std::to_string(scored) + " blocks: " + std::to_string(exact) +
                 " exact (" +
                 std::to_string(scored == 0 ? 0 : 100 * exact / scored) + "%); " +
                 std::to_string(unrep) + " unrepresentable, " + std::to_string(pruned) +
                 " pruned before scoring");
    }

    write_fragment_results(out_prefix, results, have_truth);
    log.wrote({out_prefix + ".frag_blocks.tsv", out_prefix + ".frag_pairs.tsv"});
    log.done();
    return 0;
}

} // namespace panvar

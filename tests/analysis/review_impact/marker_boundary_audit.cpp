// Marker boundary audit (review finding C1), copied verbatim from Sol's review evidence
// tests/results/shortread_code_review_20261014/marker_boundary_audit.cpp (2026-10-14). Model-only: spells every panel
// template and compares complete-sequence marker multiplicities with the marker inventory. Usage: GRAPH REFERENCE [CURATION]
#include "panvar/gfa_io.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/bubbles.hpp"
#include "panvar/integrated_snarls.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/mosaic_spelling.hpp"
#include "panvar/marker_model.hpp"
#include "panvar/syncmer.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unistd.h>

int main(int argc,char** argv) {
  if(argc<3) return 2;
  auto gfa = panvar::read_gfa_model(argv[1]);
  panvar::GraphSortOptions sort; sort.reference_path=argv[2]; sort.flip=true;
  auto sr=panvar::sort_graph_reference(gfa,sort);
  std::string stage="/tmp/panvar_marker_audit_"+std::to_string(getpid())+".gfa";
  panvar::write_gfa_model(stage,gfa);
  auto graph=panvar::parse_gfa(stage);
  std::filesystem::remove(stage);
  panvar::BubbleCallOptions bo; bo.reference_path=sr.resolved_reference;bo.quiet=true;bo.min_variant_bp=50;
  bo.snarl_pairs_override=panvar::find_top_level_snarls_cactus(panvar::snarl_input_from_model(gfa));
  bo.snarl_source_supplied=true;
  panvar::LocusInputOptions lo;lo.locus_id=sr.resolved_reference;lo.reference_path_name=sr.resolved_reference;
  if(argc>3){ std::ifstream in(argv[3]);for(std::string s;std::getline(in,s);)if(!s.empty()&&s[0]!='#')lo.excluded_templates.push_back(s); }
  auto built=panvar::build_locus_input(graph,panvar::call_bubbles(graph,bo),lo);
  auto model=panvar::prepare_genotype_model(built.input);
  panvar::verify_complete_panel_round_trip(model,built.authoritative);
  panvar::require_representation_independent(built.report);
  auto markers=panvar::build_marker_panel(model,31,15);
  std::unordered_map<std::uint64_t,std::size_t> index;
  for(std::size_t s=0;s<markers.codes.size();++s)index[markers.codes[s]]=s;
  std::vector<std::size_t> lo_count(markers.codes.size(),std::numeric_limits<std::size_t>::max()),hi_count(markers.codes.size(),0);
  std::vector<bool> mismatched(markers.codes.size(),false);
  std::size_t mismatched_anchor_template_pairs=0, mismatched_variable_template_pairs=0;
  for(std::size_t t=0;t<model.panel().template_count();++t){
    std::vector<std::size_t> actual(markers.codes.size(),0),expected(markers.codes.size(),0);
    for(auto k:panvar::collect_syncmers(panvar::spell_panel_template(model,t).bases,31,15)){auto it=index.find(k.code);if(it!=index.end())++actual[it->second];}
    for(std::size_t a=0;a<markers.anchor_slots.size();++a)expected[markers.anchor_slots[a]]=markers.anchor_haploid_multiplicity[a];
    for(std::size_t b=0;b<markers.blocks.size();++b){
      const auto& inv=markers.blocks[b];
      for(std::size_t m=0;m<inv.marker_slots.size();++m)expected[inv.marker_slots[m]]=inv.fixed_haploid_multiplicity[m];
      for(auto m:inv.allele_markers[model.panel().template_alleles[t][b]])expected[inv.marker_slots[m.first]]+=m.second;
    }
    for(std::size_t s=0;s<actual.size();++s){lo_count[s]=std::min(lo_count[s],actual[s]);hi_count[s]=std::max(hi_count[s],actual[s]);if(actual[s]!=expected[s])mismatched[s]=true;}
    for(auto s:markers.anchor_slots)if(actual[s]!=expected[s])++mismatched_anchor_template_pairs;
    for(const auto& b:markers.blocks)for(auto s:b.marker_slots)if(actual[s]!=expected[s])++mismatched_variable_template_pairs;
  }
  std::size_t bad_anchors=0,variable_anchors=0,bad_variables=0;
  for(auto s:markers.anchor_slots){bad_anchors+=mismatched[s];variable_anchors+=lo_count[s]!=hi_count[s];}
  for(const auto& b:markers.blocks)for(auto s:b.marker_slots)bad_variables+=mismatched[s];
  std::cout<<"locus="<<sr.resolved_reference<<" templates="<<model.panel().template_count()<<" blocks="<<model.blocks().size()
    <<" markers="<<markers.codes.size()<<" anchors="<<markers.anchor_slots.size()<<" anchors_wrong_multiplicity="<<bad_anchors
    <<" anchors_variable_in_full_sequences="<<variable_anchors<<" anchor_template_disagreements="<<mismatched_anchor_template_pairs
    <<" variable_markers_wrong_multiplicity="<<bad_variables<<" variable_template_disagreements="<<mismatched_variable_template_pairs<<'\n';
}

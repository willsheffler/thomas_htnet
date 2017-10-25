// -*-
// mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t
// -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// under license.
// (c) The Rosetta software is developed by the contributing members of the
// Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about
// this can be
// (c) addressed to University of Washington UW TechTransfer, email:
// license@u.washington.edu.

/// @file protocols/hbnet/HBNet.cc
/// @brief base class HBNet; to explicitly detect and design h-bond networks
/// @author Scott Boyken (sboyken@gmail.com)

#define fa_rep_for_native true
#define jack_benchmark false

#define JACK_DEBUG_GREEDY true
#include <chrono>
//#define DEBUG_TASK_JACK

// Header files for project
#include <protocols/hbnet/HBNet.hh>
#include <protocols/hbnet/HBNetCreator.hh>
#include <protocols/hbnet/HBNet_util.hh>
#include <protocols/hbnet/NetworkRotamerIDTrie.hh>
#include <protocols/hbnet/monte_carlo/monte_carlo_branch.hh>
#include <protocols/hbnet/monte_carlo/util.hh>

// Utility headers
#include <ObjexxFCL/Fstring.hh>
#include <numeric/random/random.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/numbers.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jackmag.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Core headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID.hh>
#include <core/io/Remarks.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/OffRotamerPackOptions.hh>
#include <core/pack/hbnet/HBondNetworkResidues.hh>
#include <core/pack/hbnet/HBondNetworkTracker.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/min_pack.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/selection.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/facts/FACTSEnergy.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/util.hh>
#include <core/select/residue_selector/LayerSelector.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Protocols headers
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_moves/MakePolyXMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
// XSD XRW Includes
#include <protocols/moves/mover_schemas.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// JACK
#include <core/pack/hbnet/MonteCarloHBNetOptions.hh>
#include <protocols/hbnet/monte_carlo/OffRotamerHBNet.hh>
#include <protocols/hbnet/monte_carlo/PackRotamerHBNet.hh>

// Dumping IG
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>

using namespace core;
using namespace pose;
using namespace pack;
using namespace rotamer_set;
using namespace interaction_graph;
using namespace conformation;
using namespace scoring::hbonds;

namespace protocols {
namespace hbnet {

static THREAD_LOCAL basic::Tracer TR("protocols.hbnet.HBNet");

HBNet::HBNet()
    : protocols::moves::Mover("HBNet"),
      write_network_pdbs_(0),
      jack_store_networks_(0),
      write_cst_files_(1),
      output_poly_ala_background_(0),
      find_native_(0),
      only_native_(0),
      keep_existing_networks_(0),
      extend_existing_networks_(0),
      only_extend_existing_(0),
      verbose_(0),
      symmetric_(0),
      multi_component_(0),
      show_task_(0),
      minimize_(1),
      // bridging_waters_(0),
      start_from_csts_(0),
      tyr_hydroxyls_must_donate_(0),
      hydroxyls_must_donate_(0),
      use_pdb_numbering_(1),
      no_heavy_unsats_allowed_(1),
      keep_start_selector_rotamers_fixed_(0),
      min_network_size_(3),
      max_network_size_(15),
      min_unique_networks_(1),
      min_core_res_(0),
      min_boundary_res_(0),
      max_unsat_(5),
      max_lig_unsat_(15),
      max_rep_(1),
      max_replicates_before_branch_(0),
      max_replicates_before_unsat_check_(1),
      // bridging_water_search_depth_(1),
      des_residues_("STRKHYWNQDE"),
      start_res_vec_(/* NULL */),
      network_vector_(/* NULL */),
      native_networks_(/* NULL */),
      merged_vecs_(/* NULL */),
      output_vector_(/* NULL */),
      hydrogen_bond_threshold_(-0.5),
      onebody_hb_threshold_(-0.4),
      charge_charge_rep_cutoff_(1.0),
      clash_threshold_(1.0),
      upper_score_limit_(15.0),
      min_connectivity_(0.5),
      task_factory_(/* NULL */),
      init_scorefxn_(
          core::scoring::ScoreFunctionFactory::create_score_function("HBNet")),
      scorefxn_(/* NULL */),
      rotamer_sets_(/* NULL */),
      ig_(/* NULL */),
      rotamer_links_(0),
      store_subnetworks_(0),
      secondary_search_(0),
      secondary_threshold_(-0.25),
      upweight_twobody_(1.0),
      start_selector_(/* NULL */),
      core_selector_(/* NULL */),
      boundary_selector_(/* NULL */),
      core_residues_(0),
      boundary_residues_(0),
      input_hbnet_info_residues_(0),
      // hbnet_info_residues_(0),
      run_offrot_(0),
      please_do_not_change_the_task_(0),
      init_rotset_from_monte_carlo_offrot_(0),
      init_rotset_from_monte_carlo_packrot_(0)

{
  greedy_traverse_ = basic::options::option
      [basic::options::OptionKeys::jackmag::hbnet_greedy]();
  traverse_3mers_ =
      basic::options::option[basic::options::OptionKeys::jackmag::hbnet_3mer]();
  monte_carlo_branch_ = basic::options::option
      [basic::options::OptionKeys::jackmag::hbnet_mcbranch]();
  // need to be able to use beta if called from code or another mover
  if (basic::options::option[basic::options::OptionKeys::corrections::beta]()) {
    init_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
        "HBNet_beta");
    scorefxn_ =
        core::scoring::ScoreFunctionFactory::create_score_function("beta_cst");
  } else {
    scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
        "talaris2014_cst");
  }
}

HBNet::HBNet(std::string const name)
    : protocols::moves::Mover(name),
      write_network_pdbs_(0),
      jack_store_networks_(0),
      write_cst_files_(1),
      output_poly_ala_background_(0),
      find_native_(0),
      only_native_(0),
      keep_existing_networks_(0),
      extend_existing_networks_(0),
      only_extend_existing_(0),
      verbose_(0),
      symmetric_(0),
      multi_component_(0),
      show_task_(0),
      minimize_(1),
      // bridging_waters_(0),
      start_from_csts_(0),
      tyr_hydroxyls_must_donate_(0),
      hydroxyls_must_donate_(0),
      use_pdb_numbering_(1),
      no_heavy_unsats_allowed_(1),
      keep_start_selector_rotamers_fixed_(0),
      min_network_size_(3),
      max_network_size_(15),
      min_unique_networks_(1),
      min_core_res_(0),
      min_boundary_res_(0),
      max_unsat_(5),
      max_lig_unsat_(15),
      max_rep_(1),
      max_replicates_before_branch_(0),
      max_replicates_before_unsat_check_(1),
      // bridging_water_search_depth_(1),
      des_residues_("STRKHYWNQDE"),
      start_res_vec_(/* NULL */),
      network_vector_(/* NULL */),
      native_networks_(/* NULL */),
      merged_vecs_(/* NULL */),
      output_vector_(/* NULL */),
      hydrogen_bond_threshold_(-0.5),
      onebody_hb_threshold_(-0.4),
      charge_charge_rep_cutoff_(1.0),
      clash_threshold_(1.0),
      upper_score_limit_(15.0),
      min_connectivity_(0.5),
      task_factory_(/* NULL */),
      init_scorefxn_(
          core::scoring::ScoreFunctionFactory::create_score_function("HBNet")),
      scorefxn_(/* NULL */),
      rotamer_sets_(/* NULL */),
      ig_(/* NULL */),
      rotamer_links_(0),
      store_subnetworks_(0),
      secondary_search_(0),
      secondary_threshold_(-0.25),
      upweight_twobody_(1.0),
      start_selector_(/* NULL */),
      core_selector_(/* NULL */),
      boundary_selector_(/* NULL */),
      core_residues_(0),
      boundary_residues_(0),
      input_hbnet_info_residues_(0),
      run_offrot_(0),
      please_do_not_change_the_task_(0),
      init_rotset_from_monte_carlo_offrot_(0),
      init_rotset_from_monte_carlo_packrot_(0) {
  greedy_traverse_ = basic::options::option
      [basic::options::OptionKeys::jackmag::hbnet_greedy]();
  traverse_3mers_ =
      basic::options::option[basic::options::OptionKeys::jackmag::hbnet_3mer]();
  monte_carlo_branch_ = basic::options::option
      [basic::options::OptionKeys::jackmag::hbnet_mcbranch]();
  // need to be able to use beta if called from code or another mover
  if (basic::options::option[basic::options::OptionKeys::corrections::beta]()) {
    init_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
        "HBNet_beta");
    scorefxn_ =
        core::scoring::ScoreFunctionFactory::create_score_function("beta_cst");
  } else {
    scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
        "talaris2014_cst");
  }
}

// constructor to be called from code, NEED TO CLEAN THIS UP!
HBNet::HBNet(core::scoring::ScoreFunctionCOP scorefxn, Size max_unsat,
             Size min_network_size,    /* 3 */
             Real hb_threshold,        /* -0.75 */
             Size max_network_size,    /* 15 */
             std::string des_residues, /* "STRKHYWNQDE" */
             bool find_native,         /*false*/
             bool only_native,         /*false*/
             bool keep_existing,       /*false*/
             bool extend_existing,     /*false*/
             bool only_extend          /*false*/
                                       // bool bridging_waters, /*false*/
             // bool minimize /*false*/
             )
    : protocols::moves::Mover("HBNet"),
      write_network_pdbs_(0),
      jack_store_networks_(0),
      write_cst_files_(1),
      output_poly_ala_background_(0),
      find_native_(find_native),
      only_native_(only_native),
      keep_existing_networks_(keep_existing),
      extend_existing_networks_(extend_existing),
      only_extend_existing_(only_extend),
      verbose_(0),
      symmetric_(0),
      multi_component_(0),
      show_task_(0),
      minimize_(1),
      // bridging_waters_(0),
      start_from_csts_(0),
      tyr_hydroxyls_must_donate_(0),
      hydroxyls_must_donate_(0),
      use_pdb_numbering_(1),
      no_heavy_unsats_allowed_(1),
      keep_start_selector_rotamers_fixed_(0),
      min_network_size_(min_network_size),
      max_network_size_(max_network_size),
      min_unique_networks_(1),
      min_core_res_(0),
      min_boundary_res_(0),
      max_unsat_(max_unsat),
      max_lig_unsat_(15),
      max_rep_(1),
      max_replicates_before_branch_(0),
      max_replicates_before_unsat_check_(1),
      // bridging_water_search_depth_(1),
      des_residues_(des_residues),
      start_res_vec_(/* NULL */),
      network_vector_(/* NULL */),
      native_networks_(/* NULL */),
      merged_vecs_(/* NULL */),
      output_vector_(/* NULL */),
      hydrogen_bond_threshold_(hb_threshold),
      onebody_hb_threshold_(-0.4),
      charge_charge_rep_cutoff_(1.0),
      clash_threshold_(1.0),
      upper_score_limit_(15.0),
      min_connectivity_(0.5),
      task_factory_(/* NULL */),
      init_scorefxn_(
          core::scoring::ScoreFunctionFactory::create_score_function("HBNet")),
      scorefxn_(scorefxn->clone()),
      rotamer_sets_(/* NULL */),
      ig_(/* NULL */),
      rotamer_links_(0),
      store_subnetworks_(0),
      secondary_search_(0),
      secondary_threshold_(-0.25),
      upweight_twobody_(1.0),
      start_selector_(/* NULL */),
      core_selector_(/* NULL */),
      boundary_selector_(/* NULL */),
      core_residues_(0),
      boundary_residues_(0),
      input_hbnet_info_residues_(0),
      run_offrot_(0),
      please_do_not_change_the_task_(0),
      init_rotset_from_monte_carlo_offrot_(0),
      init_rotset_from_monte_carlo_packrot_(0) {
  greedy_traverse_ = basic::options::option
      [basic::options::OptionKeys::jackmag::hbnet_greedy]();
  traverse_3mers_ =
      basic::options::option[basic::options::OptionKeys::jackmag::hbnet_3mer]();
  monte_carlo_branch_ = basic::options::option
      [basic::options::OptionKeys::jackmag::hbnet_mcbranch]();
  // need to be able to use beta if called from code or another mover
  if (basic::options::option[basic::options::OptionKeys::corrections::beta]()) {
    init_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
        "HBNet_beta");
    scorefxn_ =
        core::scoring::ScoreFunctionFactory::create_score_function("beta_cst");
  } else {
    scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
        "talaris2014_cst");
  }
  // lkb_ = new
  // core::scoring::lkball::LK_BallEnergy(init_scorefxn_->energy_method_options());
  // lkb_(init_scorefxn_->energy_method_options());
  if (only_native_) find_native_ = true;
  if (only_extend_existing_) extend_existing_networks_ = true;
}

//// TODO NEED: Copy constructor to only copy over key settings and parameters;

protocols::moves::MoverOP HBNet::clone() const {
  return (protocols::moves::MoverOP(new HBNet(*this)));
}

protocols::moves::MoverOP HBNet::fresh_instance() const {
  return protocols::moves::MoverOP(new HBNet);
}

// destructor
HBNet::~HBNet() {}

core::pack::task::TaskFactoryOP HBNet::task_factory() const {
  return task_factory_;
}

void HBNet::task_factory(core::pack::task::TaskFactoryOP task_factory) {
  task_factory_ = task_factory;
}

void HBNet::parse_my_tag(utility::tag::TagCOP tag,
                         basic::datacache::DataMap &data,
                         protocols::filters::Filters_map const &,
                         protocols::moves::Movers_map const &,
                         core::pose::Pose const &pose) {
  hydrogen_bond_threshold_ = tag->getOption<core::Real>("hb_threshold", -0.5);
  // bw_cutoff_ = tag->getOption<core::PackerEnergy>("bw_cutoff",-0.5);
  onebody_hb_threshold_ =
      tag->getOption<core::PackerEnergy>("onebody_hb_threshold", -0.4);
  charge_charge_rep_cutoff_ =
      tag->getOption<core::PackerEnergy>("charge_charge_rep_cutoff", 1.0);
  clash_threshold_ = tag->getOption<core::PackerEnergy>("clash_threshold", 1.0);
  upper_score_limit_ = tag->getOption<core::Real>("upper_score_limit", 15.0);
  find_native_ = tag->getOption<bool>("find_native_networks", 0);
  only_native_ = tag->getOption<bool>("find_only_native_networks", 0);
  keep_existing_networks_ = tag->getOption<bool>("keep_existing_networks", 0);
  extend_existing_networks_ =
      tag->getOption<bool>("extend_existing_networks", 0);
  only_extend_existing_ = tag->getOption<bool>("only_extend_existing", 0);
  minimize_ = tag->getOption<bool>("minimize", 1);
  store_subnetworks_ = tag->getOption<bool>("store_subnetworks", 0);
  secondary_search_ = tag->getOption<bool>("secondary_search", 0);
  secondary_threshold_ = tag->getOption<Real>("secondary_threshold", -0.25);
  write_network_pdbs_ = tag->getOption<bool>("write_network_pdbs", 0);
  jack_store_networks_ = tag->getOption<bool>("jack_store_networks", 0);
  write_cst_files_ = tag->getOption<bool>("write_cst_files", 1);
  output_poly_ala_background_ =
      tag->getOption<bool>("output_poly_ala_background", 0);
  max_rep_ = tag->getOption<Size>("max_replicates", 1);
  max_replicates_before_branch_ =
      tag->getOption<Size>("max_replicates_before_branch", 0);
  max_replicates_before_unsat_check_ =
      tag->getOption<Size>("max_replicates_before_unsat_check", 1);
  //    bridging_water_search_depth_ = tag->getOption< Size >(
  //    "bridging_water_search_depth", 1 );
  tyr_hydroxyls_must_donate_ =
      tag->getOption<bool>("tyr_hydroxyls_must_donate",
                           0);  // only for unsat check -- should move to filter
  hydroxyls_must_donate_ =
      tag->getOption<bool>("hydroxyls_must_donate",
                           0);  // only for unsat check -- should move to filter
  use_pdb_numbering_ = tag->getOption<bool>("use_pdb_numbering", 1);
  no_heavy_unsats_allowed_ = tag->getOption<bool>("no_heavy_unsats_allowed", 1);
  keep_start_selector_rotamers_fixed_ =
      tag->getOption<bool>("use_only_input_rot_for_start_res", false);
  show_task_ = tag->getOption<bool>("show_task", 0);
  min_network_size_ = tag->getOption<Size>("min_network_size", 3);
  max_network_size_ = tag->getOption<Size>("max_network_size", 15);
  min_unique_networks_ = tag->getOption<Size>("min_unique_networks", 1);
  min_core_res_ = tag->getOption<Size>("min_core_res", 0);
  min_boundary_res_ = tag->getOption<Size>("min_boundary_res", 0);
  max_unsat_ = tag->getOption<Size>("max_unsat", 5);
  min_connectivity_ = tag->getOption<core::Real>("min_connectivity", 0.5);
  start_from_csts_ = tag->getOption<bool>("start_from_csts", 0);
  // use_enzdes_cst_ = tag->getOption< bool >( "use_enzdes_cst", 0 );

  if (only_native_) find_native_ = true;
  if (only_extend_existing_) extend_existing_networks_ = true;
  // if ( start_from_csts_ ) {
  //  use_enzdes_cst_ = true;
  // }
  // if ( use_enzdes_cst_ ) {
  //  basic::options::option[ basic::options::OptionKeys::run::preserve_header
  //  ].value(true);
  //        write_cst_files_ = true;
  // }

  benchmark_ = tag->getOption<bool>("benchmark", 0);
  verbose_ = tag->getOption<bool>("verbose", 0);

  if (tag->hasOption("design_residues")) {
    des_residues_ =
        tag->getOption<std::string>("design_residues", "STRKHYWNQDE");
  }

  // get task operations
  if (tag->hasOption("task_operations")) {
    task_factory(protocols::rosetta_scripts::parse_task_operations(tag, data));
  }

  if (tag->hasOption("start_resnums")) {
    std::string string_resnums(
        tag->getOption<std::string>("start_resnums", ""));
    start_res_vec_ = pose::get_resnum_list(string_resnums, pose);
  }
  if (tag->hasOption("start_selector")) {
    if (tag->hasOption("start_resnums")) {
      TR.Warning << "cannot use both start_resnums and start_selector options; "
                    "start_selector will be used"
                 << std::endl;
      start_res_vec_.clear();
    }
    std::string const selector_name(
        tag->getOption<std::string>("start_selector"));
    if (TR.visible())
      TR << "Set selector name to " << selector_name << "." << std::endl;
    core::select::residue_selector::ResidueSelectorCOP selector;
    try {
      selector =
          data.get_ptr<core::select::residue_selector::ResidueSelector const>(
              "ResidueSelector", selector_name);
    } catch (utility::excn::EXCN_Msg_Exception &e) {
      std::string error_message =
          "Failed to find ResidueSelector named '" + selector_name +
          "' from the Datamap from "
          "AddCompositionConstraintMover::parse_tag()\n" +
          e.msg();
      throw utility::excn::EXCN_Msg_Exception(error_message);
    }
    runtime_assert(selector);
    start_selector_ = selector->clone();
  }
  if (tag->hasOption("core_selector")) {
    std::string const selector_name(
        tag->getOption<std::string>("core_selector"));
    if (TR.visible())
      TR << "Set selector name to " << selector_name << "." << std::endl;
    core::select::residue_selector::ResidueSelectorCOP selector;
    try {
      selector =
          data.get_ptr<core::select::residue_selector::ResidueSelector const>(
              "ResidueSelector", selector_name);
    } catch (utility::excn::EXCN_Msg_Exception &e) {
      std::string error_message =
          "Failed to find ResidueSelector named '" + selector_name +
          "' from the Datamap from "
          "AddCompositionConstraintMover::parse_tag()\n" +
          e.msg();
      throw utility::excn::EXCN_Msg_Exception(error_message);
    }
    runtime_assert(selector);
    core_selector_ = selector->clone();
    // THIS NEEDS TO HAPPEN AT APPLY TIME, NOT PARSE TIME!!!!!
    // core_residues_ = selector->apply( pose );
  }
  if (tag->hasOption("boundary_selector")) {
    std::string const selector_name(
        tag->getOption<std::string>("boundary_selector"));
    if (TR.visible())
      TR << "Set selector name to " << selector_name << "." << std::endl;
    core::select::residue_selector::ResidueSelectorCOP selector;
    try {
      selector =
          data.get_ptr<core::select::residue_selector::ResidueSelector const>(
              "ResidueSelector", selector_name);
    } catch (utility::excn::EXCN_Msg_Exception &e) {
      std::string error_message =
          "Failed to find ResidueSelector named '" + selector_name +
          "' from the Datamap from "
          "AddCompositionConstraintMover::parse_tag()\n" +
          e.msg();
      throw utility::excn::EXCN_Msg_Exception(error_message);
    }
    runtime_assert(selector);
    boundary_selector_ = selector->clone();
  }
  if (basic::options::option[basic::options::OptionKeys::corrections::beta]()) {
    TR << "mmm creating beta sfxn on line 532" << std::endl;
    init_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
        "HBNet_beta");
  }

  core::scoring::ScoreFunctionOP new_score_function(
      protocols::rosetta_scripts::parse_score_function(tag, data));
  if (new_score_function == nullptr) return;
  set_score_function(new_score_function);
}

char HBNet::get_aa_for_state(Size const res, Size const rot) const {
  return rotamer_sets_->rotamer_set_for_residue((platform::uint)(res))
      ->rotamer(rot)
      ->name1();
}
bool HBNet::res_is_boundary(Size const res) const {
  if (res > boundary_residues_.size())
    return false;
  else if (boundary_residues_[res])
    return true;
  return false;
}

bool HBNet::res_is_core(Size const res) const {
  if (res > core_residues_.size())
    return false;
  else if (core_residues_[res])
    return true;
  return false;
}

///@details Traverse IG to enumerate all possible h-bond networks;
/// networks are stored using the store_networks() function and pushed to the
/// back of network_vector_
void HBNet::traverse_IG(Real const hb_threshold) {
  if (basic::options::option
          [basic::options::OptionKeys::out::hbnet_read_conformation_file]
              .user()) {
    load_conformation_and_dump_pose(
        basic::options::option
            [basic::options::OptionKeys::out::hbnet_read_conformation_file]());
    return;
  }

  switch (basic::options::option[basic::options::OptionKeys::out::
                                     hbnet_print_interaction_graph_to_file]()) {
    case 0:
      break;
    case 1:
      print_interaction_graph_to_file("interaction_graph.txt");
      return;
    case 2:
      print_SAT_model_to_file("SAT_model.txt", false, false);
      return;
    case 3:
      print_SAT_model_to_file("PWMaxSAT_modelP.txt", true, false);
      return;
    case 4:
      print_SAT_model_to_file("PWMaxSAT_model2P.txt", false, true);
      return;
    case 5:
      print_SAT_model_to_file("PWMaxSAT_modelP2P.txt", true, true);
      return;
    case 6:
      print_ILP_model_to_file("ILP_model.lp");
      return;
    case 7:
      print_CFN_model_to_file("CFN_modelP.wcsp", true, false);
      return;
    case 8:
      print_CFN_model_to_file("CFN_model2P.wcsp", false, true);
      return;
    case 9:
      print_CFN_model_to_file("CFN_modelP2P.wcsp", true, true);
      return;
    default:
      std::cerr << "Unimplemented model\n";
      return;
  }

  // Start the IG traversal at pre-determined starting residues (e.g. interface
  // residues if HBNetStapleInterface being called)
  // utility::vector1< core::Size > all_residues;
  if (traverse_3mers_) {
    start_res_vec_.clear();
    // start_res_vec_.reserve( task_->num_to_be_packed() );
    for (core::Size resid = 1; resid <= orig_pose_->size(); ++resid) {
      if (task_->being_packed(resid)) {
        start_res_vec_.insert(resid);
      }
    }
  }

  for (core::Size res1 : start_res_vec_) {
    platform::uint const first_ni = rotamer_sets_->resid_2_moltenres(
        res1);  // returns 0 if res1 not molten (not design/repack-able)
    if (first_ni == 0) continue;
    int const first_node_ind = (int)(first_ni);

    // Loop over all edges in IG for each starting residue
    for (ig_->reset_edge_list_iterator_for_node(first_node_ind);
         !ig_->edge_list_iterator_at_end();
         ig_->increment_edge_list_iterator()) {
      EdgeBase const &edge(ig_->get_edge());
      int const second_node_ind = edge.get_other_ind(first_node_ind);
      platform::uint const second_ni = edge.get_other_ind(first_node_ind);
      Size const res2 = rotamer_sets_->moltenres_2_resid(second_ni);
      Size res1_ind(res1), res2_ind(res2);
      // if pose is symmetric, get independent resnums and set score multiply
      // factors
      // Symmetric IG only contains nodes/edges for independent residues; symm
      // clone energies added to ind edges
      if (symmetric_) {
        res1_ind = get_ind_res(*orig_pose_, res1);
        res2_ind = get_ind_res(*orig_pose_, res2);
      }
      Size scmult_2b((symmetric_) ? symm_info_->score_multiply(res1, res2) : 1);
      Size scmult_1b((symmetric_) ? symm_info_->score_multiply_factor() : 1);

      // for all rotamers of the starting position X all rotamers at the other
      // position of edge, look at energies:
      for (int ii = 1; ii <= ig_->get_num_states_for_node(first_node_ind);
           ++ii) {
        // allows derived classes to apply their own logic for starting criteria
        if (!(this->state_is_starting_aa_type(res1_ind, ii))) continue;

        // if symmetric_ need to check one-body energy to make sure rotamer
        // doesn't clash with itself accross the interface
        //    SymmetricRotamerSets compute_energies() adds 2-body e to the
        //    1-body e for cases where symm residues that interact with own
        //    clones
        if (symmetric_) {
          Real one_body_i =
              (ig_->get_one_body_energy_for_node_state(first_node_ind, ii)) /
              scmult_1b;
          if (one_body_i > clash_threshold_) {
            continue;
          }
          // check for residues that hydrogen bond with their symmetric clones;
          // these are stored in 1-body in symmeric case!!
          if (one_body_i < onebody_hb_threshold_) {
            utility::vector1<HBondResStructCOP> residues(0);
            residues.push_back(HBondResStructCOP(new hbond_res_struct(
                res1_ind, (platform::uint)(ii),
                (rotamer_sets_->rotamer_set_for_moltenresidue(first_ni)
                     ->rotamer(ii)
                     ->name1()),
                orig_pose_->pdb_info()->chain(res1_ind), 1, 0, 0)));
            store_network(residues, one_body_i, true, true, false,
                          false);  // network complete, store it
          }
        }
        for (int jj = 1; jj <= ig_->get_num_states_for_node(second_node_ind);
             ++jj) {
          if (!(this->pair_meets_starting_criteria(res1_ind, ii, res2_ind, jj)))
            continue;

          if (symmetric_) {
            Real one_body_j =
                (ig_->get_one_body_energy_for_node_state(second_node_ind, jj)) /
                scmult_1b;
            if (one_body_j > clash_threshold_) {
              continue;
            } else if (one_body_j < onebody_hb_threshold_) {
              utility::vector1<HBondResStructCOP> residues(0);
              residues.push_back(HBondResStructCOP(new hbond_res_struct(
                  res2_ind, (platform::uint)(jj),
                  (rotamer_sets_->rotamer_set_for_moltenresidue(second_ni)
                       ->rotamer(jj)
                       ->name1()),
                  orig_pose_->pdb_info()->chain(res2_ind), 1, 0, 0)));
              store_network(residues, one_body_j, true, true, false,
                            false);  // network complete, store it
            }
          }
          PDEdge const &pdedge = static_cast<PDEdge const &>(edge);
          Real twobody((first_node_ind < second_node_ind)
                           ? (pdedge.get_two_body_energy(ii, jj)) / scmult_2b
                           : (pdedge.get_two_body_energy(jj, ii)) / scmult_2b);
          twobody = twobody * (this->upweight_starting_twobody_energy());
          if (twobody < 0.0) {
            twobody = this->scale_twobody_energy(
                twobody,
                rotamer_sets_->rotamer_set_for_moltenresidue(first_ni)
                    ->rotamer(ii)
                    ->name1(),
                rotamer_sets_->rotamer_set_for_moltenresidue(second_ni)
                    ->rotamer(jj)
                    ->name1());
          }
          // if energy < threshold, we found an h-bond, and if doesn't clash
          // with anything we add the rotamer to the current network
          if (twobody < hb_threshold) {
            core::Real init_sc = twobody;
            utility::vector1<HBondResStructCOP> residues(0);
            residues.push_back(HBondResStructCOP(new hbond_res_struct(
                res1_ind, (platform::uint)(ii),
                rotamer_sets_->rotamer_set_for_moltenresidue(first_ni)
                    ->nonconst_rotamer(ii)
                    ->name1(),
                orig_pose_->pdb_info()->chain(res1_ind), 1, 0, 0)));
            if (traverse_3mers_) {
              traverse_IG_for_3mers_find_third_rotamer(
                  res1, first_ni, ii, res1_ind, res2, second_ni, jj, res2_ind,
                  hb_threshold, init_sc, residues);
            } else {
              recursive_traverse(second_node_ind, jj, res2_ind, res1_ind,
                                 residues, 1, init_sc,
                                 hb_threshold);  // recursive traverse of IG
            }
          }
        }
      }
    }
  }
}

///@details recursive call for depth-first traversal of IG
void HBNet::recursive_traverse(int const new_node_ind, int const newstate,
                               Size const newres, Size const prevres,
                               utility::vector1<HBondResStructCOP> residues,
                               Size network_rec_count, Real init_sc,
                               Real const hb_threshold,
                               bool const second_search /* false */) {
  // secondary search option allows dead-end cases to search again with a
  // lowered threshold
  if (!second_search) {
    residues.push_back(HBondResStructCOP(new hbond_res_struct(
        newres, (platform::uint)(newstate),
        (rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)
             ->rotamer(newstate)
             ->name1()),
        orig_pose_->pdb_info()->chain(newres), 1, 0, 0)));
  }
  Size prev_resnum = prevres;
  network_rec_count++;  // numer of recursive calls; network_rec_count = #
                        // rotamers in the current h-bond network

  if (network_rec_count >= max_network_size_) {  // arbitrary limit for size of
                                                 // h-bond network to prevent
                                                 // explosion
    store_network(residues, init_sc, false, false, false, false);
    return;
  }
  Size instance_rec_call_cnt(0);  // Keep track of how many times this
                                  // particular instance makes new recursive
                                  // calls

  // Continue depth-first traversal of IG starting from new_node_ind
  for (auto edge_iter = ig_->get_node(new_node_ind)->edge_list_begin();
       edge_iter != ig_->get_node(new_node_ind)->edge_list_end(); ++edge_iter) {
    runtime_assert(newres == rotamer_sets_->moltenres_2_resid(new_node_ind));
    int const second_node_ind = (*edge_iter)->get_other_ind(new_node_ind);
    platform::uint const second_ni = (*edge_iter)->get_other_ind(new_node_ind);
    Size const res2 = rotamer_sets_->moltenres_2_resid(second_ni);

    if (res2 == prev_resnum) {  // check for case where immediately return to
                                // previous residue (e.g. res_x -> res_y ->
                                // res_x)
      continue;
    }
    Size scmult_2b((symmetric_) ? symm_info_->score_multiply(newres, res2) : 1);
    Size scmult_1b((symmetric_) ? symm_info_->score_multiply_factor() : 1);

    // now search the IG for new rotamer states that have h-bonds with this
    // rotamer state
    for (int jj = 1; jj <= ig_->get_num_states_for_node(second_node_ind);
         ++jj) {
      Real one_body_j =
          (ig_->get_one_body_energy_for_node_state(second_node_ind, jj)) /
          scmult_1b;
      // if symmetric interface design need to check one-body energy to make
      // sure rotamer doesn't clash with itself accross the interface
      if (symmetric_ && one_body_j > clash_threshold_) {
        continue;
      }
      PDEdge *pdedge = static_cast<PDEdge *>(*edge_iter);
      Real twobody((new_node_ind < second_node_ind)
                       ? (pdedge->get_two_body_energy(newstate, jj)) / scmult_2b
                       : (pdedge->get_two_body_energy(jj, newstate)) /
                             scmult_2b);
      if (twobody < 0.0) {
        twobody = this->scale_twobody_energy(
            twobody,
            rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)
                ->rotamer(newstate)
                ->name1(),
            rotamer_sets_->rotamer_set_for_moltenresidue(second_node_ind)
                ->rotamer(jj)
                ->name1());
      }
      if (twobody < hb_threshold) {
        init_sc += twobody;
        bool cycle = false;
        // if we get back to any of the starting residues, stop:
        //  it's more efficient to branch overlapping networks later than detect
        //  many duplicates and subnetworks through rec explosion
        // if ( start_res_vec_.find( res2 ) == start_res_vec_.end() ) {
        if (start_res_vec_.find(res2) == start_res_vec_.end() ||
            network_rec_count < 3) {  // otherwise we could just cycle back and
                                      // forth between 2 starting reisdues and
                                      // never go deeper
          // If new res does not clash with residues already in the network...
          if (!(check_clash(residues, second_ni, jj, res2, init_sc, cycle))) {
            if (cycle) {  // If found a cycle, store network
              store_network(residues, init_sc, false, true, false, false);
            } else {  // Else, add it to the network and keep going!  Recursive
                      // call to rec_trav
              // init_sc += one_body_j;
              recursive_traverse(second_node_ind, jj, res2, newres, residues,
                                 network_rec_count, init_sc,
                                 hb_threshold);  // recursive call
              instance_rec_call_cnt++;
            }
          }
        } else {  // if we are back to original starting rotamer, don't need
                  // check_clash else check_clash
          int const start_state = static_cast<int>(residues.front()->rot_index);
          if (residues.front()->resnum == res2 && start_state == jj &&
              residues.size() > 2) {
            store_network(residues, init_sc, true, true, false, false);
            // return; 16/03
          } else if (!(check_clash(residues, second_ni, jj, res2, init_sc,
                                   cycle))) {
            // if we get back to any of the starting residues, stop:
            //  it's more efficient to branch overlapping networks later than
            //  detect many duplicates and subnetworks through rec explosion
            store_network(residues, init_sc, false, cycle, false, false);
          }
        }
      }
    }
  }  // edge_iter
  // If this instance of rec_trav doesn't make any recursive calls or store a
  // network, then store the current network.
  //   In other words, this node of this h-bond network isn't able to form any
  //   other possible h-bonds, so stop and store where we're at
  // if ( store_subnetworks_ || ( instance_rec_call_cnt == 0 && residues.size()
  // > 2 ) ){
  if (store_subnetworks_ ||
      (instance_rec_call_cnt == 0 && residues.size() > 2) ||
      rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)
              ->rotamer(newstate)
              ->name1() == 'S' ||
      rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)
              ->rotamer(newstate)
              ->name1() == 'T' ||
      rotamer_sets_->rotamer_set_for_moltenresidue(new_node_ind)
              ->rotamer(newstate)
              ->name1() == 'Y') {
    store_network(residues, init_sc, false, false, false,
                  false);  // network complete, store it
  }
  if (secondary_search_ && instance_rec_call_cnt == 0 && !second_search) {
    recursive_traverse(
        new_node_ind, newstate, newres, prevres, residues, network_rec_count,
        init_sc, secondary_threshold_,
        true);  // call this same recursive instance again with lower threshold
    instance_rec_call_cnt++;
  }
}  // rec_trav

void HBNet::traverse_IG_for_3mers_find_third_rotamer(
    core::Size resA, core::Size moltenresA, core::Size rotA,
    core::Size resA_ind, core::Size resB, core::Size moltenresB,
    core::Size rotB, core::Size resB_ind, core::Real const hb_threshold,
    core::Real init_sc, utility::vector1<HBondResStructCOP> &residues) {
  // Step 1, find partners with resB ( resC > resA )
  for (auto edge_iter = ig_->get_node(moltenresB)->edge_list_begin();
       edge_iter != ig_->get_node(moltenresB)->edge_list_end(); ++edge_iter) {
    Size const moltenresC = (*edge_iter)->get_other_ind(moltenresB);
    Size const resC = rotamer_sets_->moltenres_2_resid(moltenresC);
    if (resC <= resA) continue;

    Size const scmult_2b_k =
        (symmetric_ ? symm_info_->score_multiply(resB_ind, resC) : 1);
    Size const scmult_1b_k =
        (symmetric_ ? symm_info_->score_multiply_factor() : 1);

    // now search the IG for new rotamer states that have h-bonds with this
    // rotamer state
    for (int rotC = 1; rotC <= ig_->get_num_states_for_node(moltenresC);
         ++rotC) {
      Real one_body =
          (ig_->get_one_body_energy_for_node_state(moltenresC, rotC)) /
          scmult_1b_k;
      // if symmetric interface design need to check one-body energy to make
      // sure rotamer doesn't clash with itself accross the interface
      if (symmetric_ && one_body > clash_threshold_) {
        continue;
      }
      PDEdge *pdedge = static_cast<PDEdge *>(*edge_iter);
      Real twobody((moltenresB < moltenresC)
                       ? (pdedge->get_two_body_energy(rotB, rotC)) / scmult_2b_k
                       : (pdedge->get_two_body_energy(rotC, rotB)) /
                             scmult_2b_k);
      if (twobody < 0.0) {
        twobody = scale_twobody_energy(
            twobody,
            rotamer_sets_->rotamer_set_for_moltenresidue(moltenresB)
                ->rotamer(rotB)
                ->name1(),
            rotamer_sets_->rotamer_set_for_moltenresidue(moltenresC)
                ->rotamer(rotC)
                ->name1());
      }
      if (twobody < hb_threshold) {
        init_sc += twobody;

        if (no_clash(moltenresA, rotA, moltenresC, rotC)) {
          utility::vector1<HBondResStructCOP> residues_clone = residues;
          residues_clone.push_back(HBondResStructCOP(new hbond_res_struct(
              resC, platform::uint(rotC),
              rotamer_sets_->rotamer_set_for_moltenresidue(moltenresC)
                  ->rotamer(rotC)
                  ->name1(),
              orig_pose_->pdb_info()->chain(resC), 1, 0, 0)));
          store_network(residues_clone, init_sc, false, false, false, false);
        }
      }
    }
  }  // edge_iter

  // Step 2, find partners with resA ( resC > resB )
  for (auto edge_iter = ig_->get_node(moltenresA)->edge_list_begin();
       edge_iter != ig_->get_node(moltenresA)->edge_list_end(); ++edge_iter) {
    Size const moltenresC = (*edge_iter)->get_other_ind(moltenresA);
    Size const resC = rotamer_sets_->moltenres_2_resid(moltenresC);
    if (resC <= resB) continue;

    Size const scmult_2b_k =
        (symmetric_ ? symm_info_->score_multiply(resA_ind, resC) : 1);
    Size const scmult_1b_k =
        (symmetric_ ? symm_info_->score_multiply_factor() : 1);

    // now search the IG for new rotamer states that have h-bonds with this
    // rotamer state
    for (int rotC = 1; rotC <= ig_->get_num_states_for_node(moltenresC);
         ++rotC) {
      Real one_body =
          (ig_->get_one_body_energy_for_node_state(moltenresC, rotC)) /
          scmult_1b_k;
      // if symmetric interface design need to check one-body energy to make
      // sure rotamer doesn't clash with itself accross the interface
      if (symmetric_ && one_body > clash_threshold_) {
        continue;
      }
      PDEdge *pdedge = static_cast<PDEdge *>(*edge_iter);
      Real twobody((moltenresA < moltenresC)
                       ? (pdedge->get_two_body_energy(rotA, rotC)) / scmult_2b_k
                       : (pdedge->get_two_body_energy(rotC, rotA)) /
                             scmult_2b_k);
      if (twobody < 0.0) {
        twobody = this->scale_twobody_energy(
            twobody,
            rotamer_sets_->rotamer_set_for_moltenresidue(moltenresA)
                ->rotamer(rotA)
                ->name1(),
            rotamer_sets_->rotamer_set_for_moltenresidue(moltenresC)
                ->rotamer(rotC)
                ->name1());
      }
      if (twobody < hb_threshold) {
        init_sc += twobody;

        if (no_clash(moltenresA, rotA, moltenresC, rotC)) {
          utility::vector1<HBondResStructCOP> residues_clone = residues;
          residues_clone.push_back(HBondResStructCOP(new hbond_res_struct(
              resC, platform::uint(rotC),
              rotamer_sets_->rotamer_set_for_moltenresidue(moltenresC)
                  ->rotamer(rotC)
                  ->name1(),
              orig_pose_->pdb_info()->chain(resC), 1, 0, 0)));
          store_network(residues_clone, init_sc, false, false, false, false);
        }
      }
    }
  }  // edge_iter
}

bool compare_hbond_candidates(const hbond_candidate &a,
                              const hbond_candidate &b) {
  return a.score < b.score;
}

void HBNet::greedy_traverse_IG(Real const hb_threshold) {
  std::list<hbond_candidate> candidates;

  // Start the IG traversal at pre-deterimined starting residues (e.g. ligand or
  // interface residues)
  for (core::Size res1 : start_res_vec_) {
    platform::uint const first_ni = rotamer_sets_->resid_2_moltenres(
        res1);  // returns 0 if res1 not molten (not design/repack-able)
    if (first_ni == 0) continue;
    int const first_node_ind = (int)(first_ni);

    for (ig_->reset_edge_list_iterator_for_node(first_node_ind);
         !ig_->edge_list_iterator_at_end();
         ig_->increment_edge_list_iterator()) {
      EdgeBase const &edge(ig_->get_edge());
      int const second_node_ind = edge.get_other_ind(first_node_ind);
      platform::uint const second_ni = edge.get_other_ind(first_node_ind);
      Size const res2 = rotamer_sets_->moltenres_2_resid(second_ni);
      Size res1_ind(res1), res2_ind(res2);
      // if pose is symmetric, get independent resnums and set score multiply
      // factors
      // Symmetric IG only contains nodes/edges for independent residues; symm
      // clone energies added to ind edges
      if (symmetric_) {
        res1_ind = get_ind_res(*orig_pose_, res1);
        res2_ind = get_ind_res(*orig_pose_, res2);
      }
      Size scmult_2b((symmetric_) ? symm_info_->score_multiply(res1, res2) : 1);
      Size scmult_1b((symmetric_) ? symm_info_->score_multiply_factor() : 1);
      for (int ii = 1; ii <= ig_->get_num_states_for_node(first_node_ind);
           ++ii) {
        if (symmetric_) {
          Real one_body_i =
              (ig_->get_one_body_energy_for_node_state(first_node_ind, ii)) /
              scmult_1b;
          if (one_body_i > clash_threshold_) {
            continue;
          }
          // check for residues that hydrogen bond with their symmetric clones
          if (one_body_i < onebody_hb_threshold_) {
            // utility::vector1< HBondResStructCOP > residues(0);
            // residues.push_back( HBondResStructCOP( new hbond_res_struct(
            // res1_ind, (platform::uint)(ii),
            // (rotamer_sets_->rotamer_set_for_moltenresidue(first_ni)->rotamer(ii)->name1()),
            // orig_pose_->pdb_info()->chain(res1_ind), 1, 0, 0 ) ) );
            // store_network( residues, one_body_i, true, true, false, false );
            // //network complete, store it
          }
        }
        for (int jj = 1; jj <= ig_->get_num_states_for_node(second_node_ind);
             ++jj) {
          // ResidueOP
          // rot2(rotamer_sets_->rotamer_set_for_moltenresidue(second_ni)->nonconst_rotamer(jj));
          // //nonconst so we can add to its cache

          if (!(this->pair_meets_starting_criteria(res1_ind, ii, res2_ind, jj)))
            continue;

          if (symmetric_) {
            Real one_body_j =
                (ig_->get_one_body_energy_for_node_state(second_node_ind, jj)) /
                scmult_1b;
            if (one_body_j > clash_threshold_) {
              continue;
            } else if (one_body_j < onebody_hb_threshold_) {
              // utility::vector1< HBondResStructCOP > residues(0);
              // residues.push_back( HBondResStructCOP( new hbond_res_struct(
              // res2_ind, (platform::uint)(jj),
              // (rotamer_sets_->rotamer_set_for_moltenresidue(second_ni)->rotamer(jj)->name1()),
              // orig_pose_->pdb_info()->chain(res2_ind), 1, 0, 0 ) ) );
              // store_network( residues, one_body_j, true, true, false, false
              // ); //network complete, store it
            }
          }
          PDEdge const &pdedge = static_cast<PDEdge const &>(edge);
          Real twobody((first_node_ind < second_node_ind)
                           ? (pdedge.get_two_body_energy(ii, jj)) / scmult_2b
                           : (pdedge.get_two_body_energy(jj, ii)) / scmult_2b);
          twobody = twobody * (this->upweight_starting_twobody_energy());

          // if energy < threshold, we found an h-bond, if doesn't clash, add
          // the rotamer to the current network
          if (twobody < hb_threshold) {
            // core::Real init_sc = twobody;
            // utility::vector1< HBondResStructCOP > residues(0);
            // residues.push_back( HBondResStructCOP( new hbond_res_struct(
            // res1_ind, (platform::uint)(ii),
            // rotamer_sets_->rotamer_set_for_moltenresidue(first_ni)->nonconst_rotamer(ii)->name1(),
            // orig_pose_->pdb_info()->chain(res1_ind), 1, 0, 0 ) ) );
            // recursive_traverse( second_node_ind, jj, res2_ind, res1_ind,
            // residues, 1, init_sc, hb_threshold ); //recursive traverse of IG
            candidates.push_back(hbond_candidate(
                Size(first_node_ind), ii, Size(second_node_ind), jj, twobody));
          }
        }
      }
    }  // for edge
  }    // for res1

  candidates.sort(compare_hbond_candidates);
  for (hbond_candidate candidate_hbond : candidates) {
    // THESE TWO SHOULD ALWAYS BE IN SYNC:
    utility::vector1<core::Size> state_for_moltenres(
        rotamer_sets_->nmoltenres(), 0);
    core::pose::PoseOP network_pose = ala_pose_->clone();

    network_pose->replace_residue(
        rotamer_sets_->moltenres_2_resid(candidate_hbond.moltenres1),
        *rotamer_sets_
             ->rotamer_set_for_moltenresidue(candidate_hbond.moltenres1)
             ->rotamer(candidate_hbond.rot1),
        false);
    network_pose->replace_residue(
        rotamer_sets_->moltenres_2_resid(candidate_hbond.moltenres2),
        *rotamer_sets_
             ->rotamer_set_for_moltenresidue(candidate_hbond.moltenres2)
             ->rotamer(candidate_hbond.rot2),
        false);

    state_for_moltenres[candidate_hbond.moltenres1] = candidate_hbond.rot1;
    state_for_moltenres[candidate_hbond.moltenres2] = candidate_hbond.rot2;

    // assign the rest of the network in a depth-first order
    bool all_children_are_satisfied =
        greedy_recursive_traverse(hb_threshold, candidate_hbond.moltenres1,
                                  state_for_moltenres, network_pose) &&
        greedy_recursive_traverse(hb_threshold, candidate_hbond.moltenres2,
                                  state_for_moltenres, network_pose);

    core::scoring::hbonds::HBondSetOP hbonds(
        new core::scoring::hbonds::HBondSet(*network_pose, false /*bb only*/));
    all_children_are_satisfied =
        all_children_are_satisfied &&
        !(core_residues_[rotamer_sets_->moltenres_2_resid(
              candidate_hbond.moltenres1)] &&
          monte_carlo::has_heavy_buried_unsat(
              network_pose->residue(
                  rotamer_sets_->moltenres_2_resid(candidate_hbond.moltenres1)),
              hbonds));
    all_children_are_satisfied =
        all_children_are_satisfied &&
        !(core_residues_[rotamer_sets_->moltenres_2_resid(
              candidate_hbond.moltenres2)] &&
          monte_carlo::has_heavy_buried_unsat(
              network_pose->residue(
                  rotamer_sets_->moltenres_2_resid(candidate_hbond.moltenres2)),
              hbonds));

    if (all_children_are_satisfied) {
      Size num_core = 0;
      Size num_boundary = 0;

      utility::vector1<HBondResStructCOP> residues;
      for (core::Size mres = 1; mres <= state_for_moltenres.size(); ++mres) {
        if (state_for_moltenres[mres] != 0) {
          const Size resid = rotamer_sets_->moltenres_2_resid(mres);
          residues.push_back(HBondResStructCOP(new hbond_res_struct(
              resid, state_for_moltenres[mres],
              rotamer_sets_->rotamer_set_for_moltenresidue(mres)
                  ->nonconst_rotamer(state_for_moltenres[mres])
                  ->name1(),
              orig_pose_->pdb_info()->chain(mres), 1, 0, 0)));

          if (core_residues_[resid]) {
            ++num_core;
          }
          if (boundary_residues_[resid]) {
            ++num_boundary;
          }
        }
      }

      bool passes_filters_in_select_best_networks = true;
      if (residues.size() < 3 && min_network_size_ > 2) {
        passes_filters_in_select_best_networks = false;
      } else if (min_core_res_ && num_core < min_core_res_) {
        passes_filters_in_select_best_networks = false;
      } else if (min_boundary_res_ && num_boundary < min_boundary_res_) {
        passes_filters_in_select_best_networks = false;
      }

      if (passes_filters_in_select_best_networks) {
        store_network(residues, 0, false, false, false,
                      false);  // network complete, store it
        return;                // only create one network for now
      }
    }
  }

  // uh oh, no networks
}

bool HBNet::greedy_recursive_traverse(
    core::Real const hb_threshold, Size moltenres,
    utility::vector1<core::Size> &state_for_moltenres,
    core::pose::PoseOP &network_pose) {
  const Size resid = rotamer_sets_->moltenres_2_resid(moltenres);

  std::list<hbond_candidate> candidates;

  for (auto edge_iter = ig_->get_node(moltenres)->edge_list_begin();
       edge_iter != ig_->get_node(moltenres)->edge_list_end(); ++edge_iter) {
    int const second_node_ind = (*edge_iter)->get_other_ind(moltenres);
    platform::uint const second_ni = (*edge_iter)->get_other_ind(moltenres);
    Size const res2 = rotamer_sets_->moltenres_2_resid(second_ni);

    if (state_for_moltenres[second_ni] != 0) continue;

    Size scmult_2b((symmetric_) ? symm_info_->score_multiply(moltenres, res2)
                                : 1);
    Size scmult_1b((symmetric_) ? symm_info_->score_multiply_factor() : 1);

    // now search the IG for new rotamer states that have h-bonds with this
    // rotamer state
    for (int jj = 1; jj <= ig_->get_num_states_for_node(second_node_ind);
         ++jj) {
      if (!(this->pair_meets_starting_criteria(
              resid, state_for_moltenres[moltenres], res2, jj)))
        continue;

      Real one_body_j =
          (ig_->get_one_body_energy_for_node_state(second_node_ind, jj)) /
          scmult_1b;
      // if symmetric interface design need to check one-body energy to make
      // sure rotamer doesn't clash with itself accross the interface
      if (symmetric_ && one_body_j > clash_threshold_) {
        continue;
      }
      PDEdge *pdedge = static_cast<PDEdge *>(*edge_iter);
      Real twobody((int(moltenres) < second_node_ind)
                       ? (pdedge->get_two_body_energy(
                             state_for_moltenres[moltenres], jj)) /
                             scmult_2b
                       : (pdedge->get_two_body_energy(
                             jj, state_for_moltenres[moltenres])) /
                             scmult_2b);

      candidates.push_back(hbond_candidate(moltenres,
                                           state_for_moltenres[moltenres],
                                           Size(second_node_ind), jj, twobody));

    }  // for state jj
  }    // for edge

  utility::vector1<core::Size> state_for_moltenres_clone = state_for_moltenres;
  // core::pose::PoseOP network_pose_clone = network_pose->clone();

  // assign all of the hbonds you can, in a greedy fashion
  candidates.sort(compare_hbond_candidates);
  for (hbond_candidate candidate_hbond : candidates) {
    if (state_for_moltenres_clone[candidate_hbond.moltenres2] != 0) continue;

    // make sure there are no clashes
    bool clash = false;
    for (Size mres = 1; mres <= state_for_moltenres_clone.size(); ++mres) {
      if (mres == candidate_hbond.moltenres2 ||
          state_for_moltenres_clone[mres] == 0) {
        continue;
      }

      if (ig_->get_edge_exists(candidate_hbond.moltenres2, mres)) {
        core::pack::interaction_graph::PDEdge *edge = static_cast<PDEdge *>(
            ig_->get_node(int(candidate_hbond.moltenres2))->find_edge(mres));
        const Real temp_2body =
            (candidate_hbond.moltenres2 < mres
                 ? edge->get_two_body_energy(candidate_hbond.rot2,
                                             state_for_moltenres_clone[mres])
                 : edge->get_two_body_energy(state_for_moltenres_clone[mres],
                                             candidate_hbond.rot2));
        if (temp_2body > clash_threshold_) {
          clash = true;
          break;
        }
      }
    }
    if (clash) continue;

#ifdef JACK_DEBUG_GREEDY
    core::Size num_assigned = 0;
    for (core::Size state : state_for_moltenres_clone) {
      if (state != 0) ++num_assigned;
    }
#endif

    // temp assign state
    state_for_moltenres_clone[candidate_hbond.moltenres2] =
        candidate_hbond.rot2;
    network_pose->replace_residue(
        rotamer_sets_->moltenres_2_resid(candidate_hbond.moltenres2),
        *rotamer_sets_
             ->rotamer_set_for_moltenresidue(candidate_hbond.moltenres2)
             ->rotamer(candidate_hbond.rot2),
        false);

    bool all_children_are_satisfied =
        greedy_recursive_traverse(hb_threshold, candidate_hbond.moltenres2,
                                  state_for_moltenres_clone, network_pose);
    if (!all_children_are_satisfied) {
      // undo state assignment
      state_for_moltenres_clone[candidate_hbond.moltenres2] = 0;
      network_pose->replace_residue(
          rotamer_sets_->moltenres_2_resid(candidate_hbond.moltenres2),
          ala_pose_->residue(
              rotamer_sets_->moltenres_2_resid(candidate_hbond.moltenres2)),
          false);
#ifdef JACK_DEBUG_GREEDY
      for (core::Size state : state_for_moltenres_clone) {
        if (state != 0) --num_assigned;
      }
      if (num_assigned != 0) {
        TR << "NOT ALL DATA IS ERASED" << std::endl;
        TR << "num_assigned: " << (0 - num_assigned) << std::endl;
      }
      runtime_assert(num_assigned == 0);
#endif
    }
    // if the child nodes are satisfied (no heavy buried unsats), keep this
    // state assignment
  }

  // determine if this residue has any buried unsats
  core::scoring::hbonds::HBondSetOP hbonds(
      new core::scoring::hbonds::HBondSet(*network_pose, false /*bb only*/));
  const bool satisfied =
      !(core_residues_[rotamer_sets_->moltenres_2_resid(moltenres)] &&
        monte_carlo::has_heavy_buried_unsat(
            network_pose->residue(rotamer_sets_->moltenres_2_resid(moltenres)),
            hbonds));

  if (satisfied) {
    // I'm 100% sure this can be faster - let's just make sure the logic works
    // first
    state_for_moltenres = state_for_moltenres_clone;
  } else {
    for (core::Size ii = 1; ii <= state_for_moltenres.size(); ++ii) {
      if (state_for_moltenres[ii] != state_for_moltenres_clone[ii]) {
        core::Size const resid_to_revert = rotamer_sets_->moltenres_2_resid(ii);
        network_pose->replace_residue(
            resid_to_revert, ala_pose_->residue(resid_to_revert), false);
      }
    }
  }

  return satisfied;
}

///@details traverse EnergyGraph of static pose to find all native networks
void HBNet::traverse_native(Pose const &pose, Real const hb_threshold) {
  // only way to get here is if pose has already been scored
  // pose.update_residue_neighbors();
  //( *scorefxn_ )( pose ); //score pose so we can get EnergyGraph (EG)

  for (core::Size res : start_res_vec_) {
    utility::graph::Node const *node =
        pose.energies().energy_graph().get_node(res);

    for (utility::graph::EdgeListConstIterator egraph_it =
             node->const_edge_list_begin();
         egraph_it != node->const_edge_list_end(); ++egraph_it) {
      core::scoring::EnergyEdge const *eedge =
          static_cast<core::scoring::EnergyEdge const *>(*egraph_it);
      Size other_res =
          eedge->get_other_node(node->get_node_index())->get_node_index();

      Size res1_ind(res), res2_ind(other_res);
      if (symmetric_) {
        res1_ind = get_ind_res(*orig_pose_, res);
        res2_ind = get_ind_res(*orig_pose_, other_res);
      }

      core::Real twobody_energy = eedge->operator[](core::scoring::hbond_sc);
      if (fa_rep_for_native) {
        twobody_energy += eedge->operator[](core::scoring::fa_rep);
      }

      if (twobody_energy < hb_threshold) {
        // Initialize lists
        utility::vector1<HBondResStructCOP> residues(0);
        residues.push_back(HBondResStructCOP(
            new hbond_res_struct(res1_ind, 0, pose.residue(res1_ind).name1(),
                                 pose.pdb_info()->chain(res1_ind),
                                 pose.residue(res1_ind).is_protein(), 0,
                                 pose.residue(res1_ind).is_ligand())));

        rec_trav_native(pose, res2_ind, res1_ind, residues,
                        hydrogen_bond_threshold_);  // recursive call
      }
    }
  }
}

void HBNet::rec_trav_native(Pose const &pose, Size new_res, Size prev_res,
                            utility::vector1<HBondResStructCOP> residues,
                            Real const hb_threshold) {
  bool found_new_res_in_net = false;
  for (utility::vector1<HBondResStructCOP>::const_iterator i = residues.begin();
       i != residues.end(); ++i) {
    if ((*i)->resnum == new_res) {
      store_network(residues, 0.0, false, false, true);
      return;
    }
  }
  // residues.push_back( HBondResStructCOP( new hbond_res_struct( new_res, 0,
  // pose.residue(new_res).name1(), pose.pdb_info()->chain(new_res),
  // pose.residue(new_res).is_protein(), pose.residue(new_res).is_water(),
  // pose.residue(new_res).is_ligand() ) ) );

  residues.push_back(HBondResStructCOP(new hbond_res_struct(
      new_res, 0, pose.residue(new_res).name1(),
      pose.pdb_info()->chain(new_res), pose.residue(new_res).is_protein(),
      pose.residue(new_res).name1() == 'w',
      pose.residue(new_res).is_ligand())));

  utility::graph::Node const *node =
      pose.energies().energy_graph().get_node(new_res);
  for (utility::graph::EdgeListConstIterator egraph_it =
           node->const_edge_list_begin();
       egraph_it != node->const_edge_list_end(); ++egraph_it) {
    core::scoring::EnergyEdge const *eedge =
        static_cast<core::scoring::EnergyEdge const *>(*egraph_it);
    Size other_res =
        eedge->get_other_node(node->get_node_index())->get_node_index();

    if (symmetric_) {
      other_res = get_ind_res(*orig_pose_, other_res);
    }

    if (eedge->operator[](core::scoring::hbond_sc) < hb_threshold) {
      if (other_res == prev_res) {
        continue;
      } else if (residues.front()->resnum == other_res) {
        if (store_subnetworks_) {
          store_network(residues, 0.0, true, false, true);
          found_new_res_in_net = true;
        } else {
          rec_trav_native(pose, other_res, new_res, residues,
                          hydrogen_bond_threshold_);  // recursive call
        }
      } else {
        rec_trav_native(pose, other_res, new_res, residues,
                        hydrogen_bond_threshold_);
        found_new_res_in_net = true;  // recursive call
      }
    }
  }
  if (!found_new_res_in_net) {
    store_network(residues, 0.0, false, false, false, true);
  }
}

bool HBNet::no_clash(Size moltenres1, Size state1, Size moltenres2,
                     Size state2) {
  pack::interaction_graph::EdgeBase *edge =
      ig_->find_edge(moltenres1, moltenres2);
  if (!edge) {
    return true;
  }

  pack::interaction_graph::PDEdge *pdedge = static_cast<PDEdge *>(edge);
  Real const twobody =
      (moltenres1 < moltenres2 ? pdedge->get_two_body_energy(state1, state2)
                               : pdedge->get_two_body_energy(state2, state1));
  return twobody < clash_threshold_;
}

///@details Check if a new rotamer state clashes with any rotamer states already
/// in a given h-bond network
///    return of true = it clashes, false = no clashes
bool HBNet::check_clash(utility::vector1<HBondResStructCOP> const &residues,
                        platform::uint new_node_ind, Size newstate, Size newres,
                        Real &init_score, bool &cycle) {
  cycle = false;
  bool clash = false;
  core::PackerEnergy twobody(0.0);

  for (const auto &residue : residues) {
    platform::uint otherstate(residue->rot_index);
    platform::uint other_node_ind(
        rotamer_sets_->resid_2_moltenres(residue->resnum));
    int const new_ind = static_cast<int>(new_node_ind);
    int const other_ind =
        static_cast<int>(rotamer_sets_->resid_2_moltenres(residue->resnum));
    int const new_st = static_cast<int>(newstate);
    int const other_st = static_cast<int>(otherstate);

    if (newres == residue->resnum) {
      if (newstate == otherstate) {  // If rotamer state already in network, we
                                     // found a cycle, no need to clash check
        cycle = true;
        clash = false;
        return clash;
      } else {
        // Else, same residue, different rotamer states, so return clash = true;
        clash = true;
        return clash;
      }
    } else if (ig_->find_edge(new_ind, other_ind) !=
               0) {  // this function is order-independent; returns the same
                     // Edge regardless of order
      core::pack::interaction_graph::EdgeBase *temp_edge =
          ig_->find_edge(new_ind, other_ind);
      core::pack::interaction_graph::PDEdge *temp_pdedge =
          static_cast<PDEdge *>(temp_edge);
      int first_node_ind =
          temp_pdedge->get_first_node_ind();  // need to get the first node of
                                              // the Edge (lower-index node)
                                              // (we're getting this from the
                                              // IG)
      platform::uint first_ni = first_node_ind;
      Size res1 =
          rotamer_sets_->moltenres_2_resid(first_ni);  // get residue number
                                                       // that corresopnds to
                                                       // the first node of Edge

      // temp_pdedge->get_two_body_energy() is order-DEPENDENT; lower index
      // (first index) must always be first)
      //    so we need to determine does first_ni correspond to new_node_ind?
      if (first_ni == new_node_ind &&
          res1 == newres) {  // second half of if is sanity check; shouldn't be
                             // necessary but best to be safe
        twobody = temp_pdedge->get_two_body_energy(new_st, other_st);
      } else if (first_ni == other_node_ind && res1 == residue->resnum) {
        // or other_node_ind?
        twobody = temp_pdedge->get_two_body_energy(other_st, new_st);
      } else {
        runtime_assert((first_ni != new_node_ind && res1 == newres) ||
                       (first_ni != other_node_ind && res1 == residue->resnum));
      }

      // now that we have the correct two-body energy from the IG, check if
      // there is a clash
      if (symmetric_) {
        twobody = twobody / (symm_info_->score_multiply_factor());
      }
      if (twobody >= clash_threshold_) {
        clash = true;
        return clash;
      } else if (twobody >=
                 hydrogen_bond_threshold_) {  // include fa_rep values for
                                              // residues that may contact;
                                              // h-bonding residues ( <
                                              // hydrogen_bond_threshold)
                                              // already included
        init_score += twobody;
      }
    }
  }
  return clash;
}  // check_clash

///@details Function to check whether two h-bond networks clash with eachother
///    return of true = they clash, false = no clashes (networks i and j are
///    compatible)
bool HBNet::net_clash(hbond_net_struct const &i, hbond_net_struct const &j) {
  return net_clash(i.residues, j.residues);
}

///@details Function to check whether two h-bond networks clash with eachother
///    return of true = they clash, false = no clashes (networks i and j are
///    compatible)
bool HBNet::net_clash(utility::vector1<HBondResStructCOP> const &residues_i,
                      utility::vector1<HBondResStructCOP> const &residues_j) {
  core::PackerEnergy twobody(0.0);

  for (const auto &res_i : residues_i) {
    int const state_i = static_cast<int const>(res_i->rot_index);
    int const node_i =
        static_cast<int const>(rotamer_sets_->resid_2_moltenres(res_i->resnum));
    for (const auto &res_j : residues_j) {
      int const state_j = static_cast<int const>(res_j->rot_index);
      int const node_j = static_cast<int const>(
          rotamer_sets_->resid_2_moltenres(res_j->resnum));
      // If networks share a residue, make sure it is the same rotamer state;
      // otherwise, return clash = true;
      //    except at *(start_res_vec_.begin()) in ligand()binding design case:
      //    check for compatible ligand states in merge_networks()
      // NEED TO CHECK LIGAND COMPATIBILITY HERE //NEED TO FIX
      if ((this->ligand()) && (res_i->resnum == (this->ligand()) ||
                               res_j->resnum == (this->ligand()))) {
        continue;
      } else if (res_i->resnum == res_j->resnum && (state_i != state_j)) {
        if (rotamer_sets_ != nullptr) {
          if (res_i->aa != res_j->aa) {
            return true;
          } else {
            core::conformation::ResidueCOP r_i =
                rotamer_sets_->rotamer_set_for_residue(res_i->resnum)
                    ->rotamer(state_i);
            core::conformation::ResidueCOP r_j =
                rotamer_sets_->rotamer_set_for_residue(res_j->resnum)
                    ->rotamer(state_j);
            // Real sc_rmsd(core::scoring::automorphic_rmsd(*r_i, *r_j, true));
            // if ( core::scoring::automorphic_rmsd(*r_i, *r_j, true) >
            // SC_RMSD_CUTOFF ) {
            if (core::scoring::residue_sc_rmsd_no_super(r_i, r_j, true) >
                SC_RMSD_CUTOFF) {
              return true;
            }
          }
        } else {
          return true;
        }
      }

      if (ig_ != nullptr) {  // if IG still exists, use 2-body for faster lookup
        if (ig_->find_edge(node_i, node_j) != 0) {  // If edge exist s, look up
                                                    // 2-body energy for rotamer
                                                    // state_i and state_j
          core::pack::interaction_graph::EdgeBase *temp_edge =
              ig_->find_edge(node_i, node_j);
          core::pack::interaction_graph::PDEdge *temp_pdedge =
              static_cast<PDEdge *>(temp_edge);
          int first_node_ind = temp_pdedge->get_first_node_ind();

          if (first_node_ind == node_i) {
            twobody = temp_pdedge->get_two_body_energy(state_i, state_j);
          } else if (first_node_ind == node_j) {
            twobody = temp_pdedge->get_two_body_energy(state_j, state_i);
          } else {
            runtime_assert(first_node_ind == node_i ||
                           first_node_ind == node_j);
          }
          if (symmetric_) {
            twobody = twobody / (symm_info_->score_multiply(res_i->resnum,
                                                            res_j->resnum));
          }
          if (twobody >= clash_threshold_) {
            return true;
          }
        }
      } else {  // safeguard: if IG has been cleared or reset to null, we can
                // still check clash
        Pose ala_copy = *ala_pose_;
        core::conformation::ResidueCOP rot_i =
            rotamer_sets_->rotamer_set_for_moltenresidue(node_i)->rotamer(
                res_i->rot_index);
        core::conformation::ResidueCOP rot_j =
            rotamer_sets_->rotamer_set_for_moltenresidue(node_j)->rotamer(
                res_j->rot_index);

        ala_copy.replace_residue(res_i->resnum, *rot_i, false);
        ala_copy.replace_residue(res_j->resnum, *rot_j, false);
        (*scorefxn_)(
            ala_copy);  // sfxn->eval methods ASSUME THE POSE HAS BEEN SCORED?
        core::scoring::EnergyMap new_emap;
        init_scorefxn_->eval_cd_2b_sc_sc(*rot_i, *rot_j, ala_copy, new_emap);
        Real farep_val = new_emap.get(core::scoring::fa_rep);
        if (farep_val > clash_threshold_) {
          return true;
        }
      }
    }
  }
  return false;
}  // net_clash

bool HBNet::network_already_stored(
    utility::vector1<HBondResStructCOP> &residues,
    utility::vector1<HBondResStructCOP> &i_residues) {
  bool already_stored(false);
  if (residues.size() == i_residues.size()) {
    std::sort(residues.begin(), residues.end(), compare_hbond_residues());
    std::sort(i_residues.begin(), i_residues.end(), compare_hbond_residues());

    already_stored = true;
    std::vector<HBondResStructCOP>::const_iterator j = i_residues.begin();
    for (std::vector<HBondResStructCOP>::const_iterator k = residues.begin();
         k != residues.end(); ++k) {
      if (((*j)->resnum == (*k)->resnum) && ((*j)->aa == (*k)->aa)) {
        if (rotamer_sets_ == nullptr) continue;
        core::conformation::ResidueCOP rot1(
            rotamer_sets_
                ->rotamer_set_for_moltenresidue(
                    rotamer_sets_->resid_2_moltenres((*j)->resnum))
                ->rotamer((*j)->rot_index));
        core::conformation::ResidueCOP rot2(
            rotamer_sets_
                ->rotamer_set_for_moltenresidue(
                    rotamer_sets_->resid_2_moltenres((*k)->resnum))
                ->rotamer((*k)->rot_index));
        // if ( core::scoring::automorphic_rmsd(*rot1, *rot2, true /*superpose*/
        // ) > SC_RMSD_CUTOFF ) {
        if (core::scoring::residue_sc_rmsd_no_super(
                rot1, rot2, true /*final group only*/) > SC_RMSD_CUTOFF) {
          already_stored = false;
          break;
        }
      } else {
        already_stored = false;
        break;
      }
      ++j;
    }
  }
  return already_stored;
}

void HBNet::store_network(utility::vector1<HBondResStructCOP> residues,
                          Real init_score, bool term_w_start, bool term_w_cycle,
                          bool score_now, bool native) {
  bool already_stored(false);
  for (auto &i : network_vector_) {
    already_stored = network_already_stored(residues, i->residues);

    if (already_stored) {
      if (init_score < i->score && !score_now) {
        // replace i with new network with better score
        i->is_native = native;
        i->term_w_start = term_w_start;
        i->term_w_cycle = term_w_cycle;
        i->residues = residues;
        if ((this->ligand())) {
          i->lig_state_list.push_back(
              (*(find_hbond_res_struct(i->residues, this->ligand())))
                  ->rot_index);
        }
        i->score = init_score;  // don't need to adjust for symmetric_ because
                                // did that when tracking IG energies to
                                // generate init_score
        // if (normalize_)
        //    (*i)->score = (*i)->score/(residues.size());
      } else if (this->ligand()) {
        i->lig_state_list.push_back(
            (*(find_hbond_res_struct(i->residues, this->ligand())))->rot_index);
      }
      return;
    }
  }
  if (!(already_stored)) {
    HBondNetStructOP new_net(new hbond_net_struct());
    new_net->is_native = native;
    new_net->term_w_start = term_w_start;
    new_net->term_w_cycle = term_w_cycle;
    new_net->residues = residues;
    new_net->lig_state_list.clear();
    if ((this->ligand())) {
      new_net->lig_state_list.push_back(
          (*(find_hbond_res_struct((new_net)->residues, this->ligand())))
              ->rot_index);
    }
    if (score_now) {
      Pose copy = *ala_pose_;
    } else {
      new_net->score = init_score;  // don't need to adjust for symmetric_
                                    // because did that when tracking IG
                                    // energies to generate init_score
      network_vector_.push_back(new_net);
    }
  }
}  // store_network

void HBNet::minimize_network(Pose &pose, hbond_net_struct &network,
                             bool residues_already_placed /* true */) {
  if (!residues_already_placed) {
    place_rots_on_pose(pose, network, network.is_native);
    pose.update_residue_neighbors();
  }
  core::kinematics::MoveMapOP mm =
      core::kinematics::MoveMapOP(new core::kinematics::MoveMap);
  mm->set_bb(false);
  mm->set_chi(false);
  mm->set_jump(false);
  for (utility::vector1<HBondResStructCOP>::const_iterator rit =
           network.residues.begin();
       rit != network.residues.end(); rit++) {
    if ((this->ligand()) && (*rit)->resnum == (this->ligand())) {
      continue;
    }
    mm->set_chi((*rit)->resnum, true);
  }
  // TODO NEED TO ADD CODE HERE FOR BRIDGING_WATER CASE; NEED TO CHECK RESIDUE
  // NUMBERS OF NETWORK RESIDUES AFTER WATER
  protocols::simple_moves::MinMoverOP min_mover;
  if (symmetric_) {
    core::pose::symmetry::make_symmetric_movemap(pose, *mm);
    min_mover = protocols::simple_moves::MinMoverOP(
        new protocols::simple_moves::symmetry::SymMinMover(
            mm, scorefxn_, "dfpmin_armijo_nonmonotone", 0.0001, true));
  } else {
    min_mover = protocols::simple_moves::MinMoverOP(
        new protocols::simple_moves::MinMover(
            mm, scorefxn_, "dfpmin_armijo_nonmonotone", 0.0001, true));
  }
  min_mover->apply(pose);
  pose.update_residue_neighbors();

  // update roatmers and store, so we don't have to do this again
  Size old_total((symmetric_) ? symm_info_->num_independent_residues()
                              : orig_pose_->total_residue());
  // Size new_total( pose.total_residue() );
  // if ( symmetric_ ) {
  //  core::conformation::symmetry::SymmetricConformation const &
  //  SymmConf(dynamic_cast<core::conformation::symmetry::SymmetricConformation
  //  const & > ( pose.conformation()));
  //  new_total = SymmConf.Symmetry_Info()->num_independent_residues();
  // }
  network.rotamers.clear();
  for (Size r = 1; r <= old_total; r++) {
    if (find_hbond_res_struct(network.residues, r) != network.residues.end()) {
      network.rotamers.push_back(pose.residue(r).clone());
    }
  }
}

void HBNet::score_network_on_pose(Pose &pose, hbond_net_struct &i) {
  if (i.scored) return;

  // Pose ala_copy = *ala_pose_; // TODO does scoring get copied too? looks like
  // it does, but NEED TO CHECK THIS
  // ala_copy.update_residue_neighbors();
  //(*scorefxn_)(ala_copy); //already scored

  // set baseline energies from background pose
  core::scoring::EnergyMap &baseline_emap =
      ala_pose_->energies().total_energies();
  core::Real total_baseline = baseline_emap.get(core::scoring::total_score);
  core::Real hbond_bb_sc_baseline =
      baseline_emap.get(core::scoring::hbond_bb_sc);

  // place_rots_on_pose( ala_copy, i, i.is_native, bridging_waters_,
  // (*i)->waterrots.empty() );
  // ala_copy.update_residue_neighbors();

  (*scorefxn_)(pose);
  core::scoring::EnergyMap &new_emap = pose.energies().total_energies();
  Real pose_total = new_emap.get(core::scoring::total_score) - total_baseline;
  Real e(pose_total);
  if ((new_emap.get(core::scoring::hbond_bb_sc) - hbond_bb_sc_baseline) <=
      MIN_HB_E_CUTOFF) {
    i.term_w_bb = true;
  } else {
    i.term_w_bb = false;
  }
  // NEED TO FIX
  core::Real cycle_bonus = -2.0;
  if (i.term_w_cycle) e = e + cycle_bonus;
  // if ( (this->ligand()) && i.term_w_start ) e = e + ret2lig_bonus_;
  if (symmetric_) e = e / (symm_info_->num_bb_clones() + 1.0);
  // normalize
  e = e / ((Real)i.residues.size());
  i.score = e;
  i.scored = true;
}

void HBNet::merge_2_branched_networks(
    utility::vector1<HBondResStructCOP> const &residues1,
    utility::vector1<HBondResStructCOP> const &residues2,
    utility::vector1<HBondResStructCOP> &new_residues) {
  for (const auto &res1 : residues1) {
    bool found(false);
    for (utility::vector1<HBondResStructCOP>::const_iterator n =
             new_residues.begin();
         n != new_residues.end(); ++n) {
      if (res1->resnum == (*n)->resnum) {
        found = true;
        break;
      }
    }
    if (!found) {
      new_residues.push_back(res1);
    }
  }
  for (const auto &res2 : residues2) {
    bool found(false);
    for (utility::vector1<HBondResStructCOP>::const_iterator n =
             new_residues.begin();
         n != new_residues.end(); ++n) {
      if (res2->resnum == (*n)->resnum) {
        found = true;
        break;
      }
    }
    if (!found) {
      new_residues.push_back(res2);
    }
  }
}

void HBNet::merge_2_branched_networks(hbond_net_struct const &i,
                                      hbond_net_struct const &j,
                                      HBondNetStructOP new_network) {
  if ((this->ligand()) &&
      (i.lig_state_list.empty() || j.lig_state_list.empty())) {
    return;
  }
  utility::vector1<HBondResStructCOP> new_residues(0);
  merge_2_branched_networks(i.residues, j.residues, new_residues);

  new_network->sort_first_by_tot_unsat = true;
  new_network->scored = false;
  new_network->residues = new_residues;
  new_network->is_native = (i.is_native && j.is_native) ? true : false;
  new_network->is_extended =
      ((i.is_native && !j.is_native) || (!i.is_native && j.is_native) ||
       (i.is_extended || j.is_extended))
          ? true
          : false;
  new_network->term_w_bb = (i.term_w_bb || j.term_w_bb) ? true : false;
  new_network->term_w_cycle = (i.term_w_cycle || j.term_w_cycle) ? true : false;
  new_network->term_w_start = (i.term_w_start || j.term_w_start) ? true : false;

  if (this->ligand()) {
    std::vector<platform::uint> new_lig_state_list(1000);
    std::vector<platform::uint> liglist1 = i.lig_state_list;
    std::vector<platform::uint> liglist2 = j.lig_state_list;
    std::sort(liglist1.begin(), liglist1.end());
    std::sort(liglist2.begin(), liglist2.end());
    std::vector<platform::uint>::iterator lit;
    lit = std::set_intersection(liglist1.begin(), liglist1.end(),
                                liglist2.begin(), liglist2.end(),
                                new_lig_state_list.begin());
    new_lig_state_list.resize(lit - new_lig_state_list.begin());
    if (new_lig_state_list.empty()) {
      new_lig_state_list.push_back(
          (*(find_hbond_res_struct(i.residues, this->ligand())))->rot_index);
      new_lig_state_list.push_back(
          (*(find_hbond_res_struct(j.residues, this->ligand())))->rot_index);
      new_network->lig_state_list = new_lig_state_list;
    } else {
      new_network->lig_state_list = new_lig_state_list;
    }
  }
  // score and num_unsat, etc. should be set by place_rotamers_and_score()
  new_network->lig_num_unsatisfied = 0;
  new_network->num_unsat = 0;
  new_network->num_heavy_unsat = 0;
  new_network->score = 0.0;
}

void HBNet::merge_branched_networks(core::Size merged_vec_index,
                                    HBondNetStructOP new_network) {
  std::vector<core::Size> const &network_vec_indices =
      merged_vecs_[merged_vec_index];

  if (this->ligand()) {
    for (core::Size netid : network_vec_indices) {
      if (network_vector_[netid]->lig_state_list.empty()) {
        return;
      }
    }
  }

  new_network->is_native = true;
  bool any_are_extended = false;
  bool any_are_native = false;

  new_network->term_w_bb = false;
  new_network->term_w_cycle = false;
  new_network->term_w_start = false;

  utility::vector1<HBondResStructCOP> new_residues;
  utility::vector1<bool> resid_added(orig_pose_->size(), false);

  for (core::Size netid : network_vec_indices) {
    auto net = network_vector_[netid];

    // properties
    if (net->is_native) {
      any_are_native = true;
    } else {
      new_network->is_native = false;
    }
    if (net->is_extended) any_are_extended = true;
    if (net->term_w_bb) new_network->term_w_bb = true;
    if (net->term_w_cycle) new_network->term_w_cycle = true;
    if (net->term_w_start) new_network->term_w_start = true;

    // residues
    for (const auto &res : net->residues) {
      if (!resid_added[res->resnum]) {
        new_residues.push_back(res);
        resid_added[res->resnum] = true;
      }
    }
  }

  new_network->sort_first_by_tot_unsat = true;
  new_network->scored = false;
  new_network->residues = new_residues;
  new_network->is_extended =
      any_are_extended || (any_are_native && !new_network->is_native);

  if (this->ligand()) {
    core::Size const i = network_vec_indices[1];
    core::Size const j = network_vec_indices[2];
    std::vector<platform::uint> new_lig_state_list(1000);
    std::vector<platform::uint> liglist1 = network_vector_[i]->lig_state_list;
    std::vector<platform::uint> liglist2 = network_vector_[j]->lig_state_list;
    std::sort(liglist1.begin(), liglist1.end());
    std::sort(liglist2.begin(), liglist2.end());
    std::vector<platform::uint>::iterator lit;
    lit = std::set_intersection(liglist1.begin(), liglist1.end(),
                                liglist2.begin(), liglist2.end(),
                                new_lig_state_list.begin());
    new_lig_state_list.resize(lit - new_lig_state_list.begin());
    if (new_lig_state_list.empty()) {
      new_lig_state_list.push_back(
          (*(find_hbond_res_struct(network_vector_[i]->residues,
                                   this->ligand())))
              ->rot_index);
      new_lig_state_list.push_back(
          (*(find_hbond_res_struct(network_vector_[j]->residues,
                                   this->ligand())))
              ->rot_index);
      new_network->lig_state_list = new_lig_state_list;
    } else {
      new_network->lig_state_list = new_lig_state_list;
    }

    for (core::Size ii = 3; ii <= network_vec_indices.size(); ++ii) {
      core::Size const netid = network_vec_indices[ii];
      auto net_ii = network_vector_[netid];
      std::vector<platform::uint> liglist_ii = net_ii->lig_state_list;
      std::sort(liglist_ii.begin(), liglist_ii.end());
      new_lig_state_list.clear();
      lit = std::set_intersection(new_network->lig_state_list.begin(),
                                  new_network->lig_state_list.end(),
                                  liglist_ii.begin(), liglist_ii.end(),
                                  new_lig_state_list.begin());
      new_lig_state_list.resize(lit - new_lig_state_list.begin());
      if (new_lig_state_list.empty()) {
        new_lig_state_list.push_back(
            (*(find_hbond_res_struct(new_network->residues, this->ligand())))
                ->rot_index);
        new_lig_state_list.push_back(
            (*(find_hbond_res_struct(net_ii->residues, this->ligand())))
                ->rot_index);
        new_network->lig_state_list = new_lig_state_list;
      } else {
        new_network->lig_state_list = new_lig_state_list;
      }
    }
  }
  // score and num_unsat, etc. should be set by place_rotamers_and_score()
  new_network->lig_num_unsatisfied = 0;
  new_network->num_unsat = 0;
  new_network->num_heavy_unsat = 0;
  new_network->score = 0.0;
}

// consider Ser and Thr idential for benchmarking purposes
bool HBNet::networks_identical_aa_sequence(hbond_net_struct const &i,
                                           hbond_net_struct const &j) {
  if (i.residues.size() == j.residues.size()) {
    std::vector<char> i_aa(0);
    for (const auto &residue : i.residues) {
      if (residue->aa == 'T' || residue->aa == 'S') {
        i_aa.push_back('S');
      } else {
        i_aa.push_back(residue->aa);
      }
    }
    std::vector<char> j_aa(0);
    for (const auto &residue : j.residues) {
      if (residue->aa == 'T' || residue->aa == 'S') {
        j_aa.push_back('S');
      } else {
        j_aa.push_back(residue->aa);
      }
    }
    std::sort(i_aa.begin(), i_aa.end());
    std::sort(j_aa.begin(), j_aa.end());
    std::vector<int> a(20);
    std::vector<int>::iterator ait;
    ait = std::set_symmetric_difference(i_aa.begin(), i_aa.end(), j_aa.begin(),
                                        j_aa.end(), a.begin());
    a.resize(ait - a.begin());

    if (a.size() == 0) {
      return true;
    }
  }
  return false;
}

bool HBNet::is_sub_residues(utility::vector1<HBondResStructCOP> &residues1,
                            utility::vector1<HBondResStructCOP> &residues2) {
  bool branch(false);
  return is_sub_residues(residues1, residues2, branch);
}

///@details identical networks are considered subsets of each other here by
/// default unless specified by true_if_identical=false
///  std::includes() returns false if 2 vectors are identical
bool HBNet::is_sub_residues(utility::vector1<HBondResStructCOP> &residues1,
                            utility::vector1<HBondResStructCOP> &residues2,
                            bool &branch, bool true_if_identical /* true */) {
  if (residues1.empty() || residues2.empty()) {
    return false;
  }
  // If networks is not a subset and has > 1 symmetric difference
  std::sort(residues1.begin(), residues1.end(), compare_hbond_residues());
  std::sort(residues2.begin(), residues2.end(), compare_hbond_residues());
  std::vector<HBondResStructCOP> intersec;
  std::vector<HBondResStructCOP> symm_diff;
  std::set_intersection(residues1.begin(), residues1.end(), residues2.begin(),
                        residues2.end(), std::back_inserter(intersec),
                        compare_hbond_residues());
  std::set_symmetric_difference(
      residues1.begin(), residues1.end(), residues2.begin(), residues2.end(),
      std::back_inserter(symm_diff), compare_hbond_residues());

  // return true if networks identical
  if (true_if_identical &&
      (residues1.size() == residues2.size() && symm_diff.size() == 0)) {
    return true;
  }

  // to be a subnetwork, have to overlap at at least one position (2 for
  // ligand-binding design because will overlap at ligand by default)
  // if ( ( (this->ligand()) && v.size() > 1 ) || ( !(this->ligand()) &&
  // v.size() > 0 )){
  if (intersec.size() > 0) {
    if (residues1.size() > residues2.size() &&
        (std::includes(residues1.begin(), residues1.end(), residues2.begin(),
                       residues2.end(), compare_hbond_residues()))) {
      return true;
    } else if (residues2.size() > residues1.size() &&
               (std::includes(residues2.begin(), residues2.end(),
                              residues1.begin(), residues1.end(),
                              compare_hbond_residues()))) {
      return true;
    }
    if (symm_diff.size() >= 2) {  // only branches if at least 1 diff is unique
                                  // to each list, !subreslist (e.g. if list1 =
                                  // {1,2} and list2 = {1,2,3,4}, v2 = 2 but
                                  // list2 includes list1)
      branch = true;
    }
  }
  return false;
}

///@details Critical function to HBNet; because of depth-frist traversal, need
/// to combine branched cases, e.g. ASN that can make multiple h-bonds
///   The reason we need this is because bread-first additions ("branches") are
///   order-dependent, meaning we need to try ALL combinations
void HBNet::branch_overlapping_networks() {
  if (network_vector_.empty()) {
    if (verbose_ && TR.visible()) {
      TR << " NO NETWORKS TO BRANCH; RETURNING WITHOUT DOING ANYTHING."
         << std::endl;
    }
    return;
  }
  merged_vecs_.clear();  // unnecessary if function only used internally

  if (monte_carlo_branch_) {
    branch_overlapping_networks_monte_carlo();
    return;
  }

  for (auto &i : network_vector_) {
    i->net_indices.clear();
  }
  if (TR.visible()) TR << "recording net_indices " << std::endl;
  for (auto i = network_vector_.begin(); i != network_vector_.end(); ++i) {
    // Compare network i with all j, and store best-scoring j that doesn't clash
    // with i
    for (auto j = i; ++j != network_vector_.end();) {
      // If networks overlap by more than 2 residue (will overlap at ligand),
      // continue.
      bool branch(false);
      if (is_sub_residues((*i)->residues, (*j)->residues, branch) || !branch) {
        continue;
      }

      // TODO need check for small molecule ligands here
      if (only_native_ || !(net_clash(**i, **j))) {
        (*i)->net_indices.push_back(j - network_vector_.begin());
        (*j)->net_indices.push_back(i - network_vector_.begin());
      }
    }
  }
  if (TR.visible()) TR << "starting recursion " << std::endl;
  for (auto i = network_vector_.begin(); i != network_vector_.end(); ++i) {
    Size i_pos = i - network_vector_.begin();
    std::vector<Size> add_index_vec;
    add_index_vec.push_back(i_pos);

    for (auto j = (*i)->net_indices.begin(); j != (*i)->net_indices.end();
         ++j) {
      if (*j <= (i_pos)) {
        continue;
      }
      add_index_vec.push_back(*j);

      if (store_subnetworks_) {
        bool found(false);
        std::sort(add_index_vec.begin(), add_index_vec.end());
        for (auto &merged_vec : merged_vecs_) {
          // std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS
          // ALEADY BE SORTED!
          std::vector<int> v2(100);
          auto it2 = std::set_symmetric_difference(
              merged_vec.begin(), merged_vec.end(), add_index_vec.begin(),
              add_index_vec.end(), v2.begin());
          v2.resize(it2 - v2.begin());

          if (v2.size() == 0) {
            found = true;
            break;
          }
        }
        if (!found) {
          merged_vecs_.push_back(add_index_vec);
        }
      }
      std::sort(network_vector_[*j]->net_indices.begin(),
                network_vector_[*j]->net_indices.end());
      std::sort((*i)->net_indices.begin(), (*i)->net_indices.end());
      std::vector<Size> v(1000);
      std::vector<Size>::iterator it;
      it = std::set_intersection(
          (*i)->net_indices.begin(), (*i)->net_indices.end(),
          (network_vector_[*j])->net_indices.begin(),
          (network_vector_[*j])->net_indices.end(), v.begin());
      v.resize(it - v.begin());

      if (v.size() == 0) {
        if (!store_subnetworks_) {
          bool found(false);
          std::sort(add_index_vec.begin(), add_index_vec.end());
          for (auto &merged_vec : merged_vecs_) {
            // std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS
            // ALEADY BE SORTED!
            std::vector<int> v2(100);
            auto it2 = std::set_symmetric_difference(
                merged_vec.begin(), merged_vec.end(), add_index_vec.begin(),
                add_index_vec.end(), v2.begin());
            v2.resize(it2 - v2.begin());

            if (v2.size() == 0) {
              found = true;
              break;
            }
          }
          if (!found) {
            merged_vecs_.push_back(add_index_vec);
          }
        }
        add_index_vec.pop_back();
        continue;
      } else if (v.size() == 1) {
        add_index_vec.push_back(v.front());
        bool found(false);
        std::sort(add_index_vec.begin(), add_index_vec.end());
        for (auto &merged_vec : merged_vecs_) {
          // std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS
          // ALEADY BE SORTED!
          std::vector<int> v2(100);
          auto it2 = std::set_symmetric_difference(
              merged_vec.begin(), merged_vec.end(), add_index_vec.begin(),
              add_index_vec.end(), v2.begin());
          v2.resize(it2 - v2.begin());

          if (v2.size() == 0) {
            found = true;
            break;
          }
        }
        if (!found) {
          merged_vecs_.push_back(add_index_vec);
        }
        add_index_vec.pop_back();
      } else {
        rec_set_intersection(add_index_vec, v, *j);
      }
      add_index_vec.pop_back();
    }
  }
  if (TR.visible()) TR << "finished recursion " << std::endl;
  if (TR.visible())
    TR << "merged_vecs_.size() =  " << merged_vecs_.size() << std::endl;
  // std::vector< std::vector< Size > > temp_merge_vec(0);
  // Size mcount(1);
  // for ( auto & merged_vec : merged_vecs_ ) {
  //  if ( merged_vec.empty() ) continue; //should never happen
  //  std::sort(merged_vec.begin(), merged_vec.end());
  //  bool found(false);
  //  for ( auto & t : temp_merge_vec ) {
  //   std::sort(t.begin(), t.end());
  //   std::vector< int > v2(50);
  //   auto it2 = std::set_symmetric_difference(merged_vec.begin(),
  //   merged_vec.end(), t.begin(), t.end(), v2.begin());
  //   v2.resize(it2-v2.begin());
  //
  //   if ( v2.size() == 0 ) {
  //    found = true;
  //    break;
  //   }
  //  }
  //  if ( !found ) {
  //   temp_merge_vec.push_back(merged_vec);
  //  }
  // }
  // merged_vecs_.clear();
  // for ( std::vector< std::vector< Size > >::const_iterator mit =
  // temp_merge_vec.begin(); mit != temp_merge_vec.end(); ++mit ) {
  //  merged_vecs_.push_back( *mit );
  // }

  finalize_branching();

}  // branch_overlapping networks

void HBNet::finalize_branching() {
  if (TR.visible())
    TR << "finished temp_merge_vec, actually merging networks " << std::endl;
  Size mcount(1);
  std::vector<HBondNetStructOP> temp_net_vec(0);
  for (auto &merged_vec : merged_vecs_) {
    if (merged_vec.size() > 1) {
      auto m1 = merged_vec.begin();
      HBondNetStructOP temp_hbond_net_struct = network_vector_[*m1];
      std::string nets_being_merged("");
      if (verbose_) {
        nets_being_merged = "       " + utility::to_string(mcount) +
                            ". branching together networks = ";
        nets_being_merged =
            nets_being_merged + " ( #" + (utility::to_string(*m1 + 1)) + ": ";
        std::string tempstr = print_list_to_string(*temp_hbond_net_struct,
                                                   false, false, false, false);
        // std::string tempstr =
        // print_list_to_string(temp_hbond_net_struct.rotlist,
        // temp_hbond_net_struct.term_w_start,
        //                                           temp_hbond_net_struct.term_w_cycle,
        //                                           temp_hbond_net_struct.term_w_bb
        //                                           );
        nets_being_merged = nets_being_merged + tempstr + " ) + ";
      }
      for (auto m2 = m1; ++m2 != merged_vec.end();) {
        if (verbose_) {
          nets_being_merged = nets_being_merged + " ( #" +
                              (utility::to_string((*m2 + 1))) + ": ";
          std::string tempstr = print_list_to_string(
              *(network_vector_[*m2]), false, false, false, false);
          // std::string tempstr =
          // print_list_to_string(network_vector_[*m2].rotlist,
          // network_vector_[*m2].term_w_start,
          //                                           network_vector_[*m2].term_w_cycle,
          //                                           network_vector_[*m2].term_w_bb);
          nets_being_merged = nets_being_merged + tempstr + " ) + ";
        }
        HBondNetStructOP new_network(new hbond_net_struct());
        merge_2_branched_networks(*temp_hbond_net_struct,
                                  *(network_vector_[*m2]), new_network);
        temp_hbond_net_struct = new_network;
      }
      if (verbose_ && TR.visible()) {
        TR << nets_being_merged << std::endl;
      }

      // merge_branched_networks( mcount, temp_hbond_net_struct );

      temp_net_vec.push_back(temp_hbond_net_struct);
      mcount++;
    }
  }
  // push the merged networks to the back of network_vector_
  // network_vector_.clear(); DON'T WANT TO CLEAR HERE
  // network_vector_.reserve( network_vector_.size() + temp_net_vec.size() );
  network_vector_.insert(network_vector_.end(), temp_net_vec.begin(),
                         temp_net_vec.end());
  /*for ( auto & netit : temp_net_vec ) {
          network_vector_.push_back( netit );
          }*/

  // Sort all individual networks by best score, so start with highest scoring
  std::sort(network_vector_.begin(), network_vector_.end(), compare_net_vec());
  merged_vecs_.clear();
}

inline unsigned long combine_indices(core::Size const i, core::Size const j) {
  return (static_cast<unsigned long>(i) << 32) & (j & 0xFFFFFFFF);
}

inline std::pair<core::Size, core::Size> split_indices(
    unsigned long const combined) {
  core::Size const i = (combined >> 32) & 0xFFFFFFFF;
  core::Size const j = combined & 0xFFFFFFFF;
  return std::make_pair(i, j);
}

void HBNet::branch_overlapping_networks_monte_carlo() {
  runtime_assert(
      network_vector_.size() <=
      0xFFFFFFFF);  // network_vector_.size() can be expressed as 32 bit int
  branch_seed_network_vector_.clear();

  std::set<unsigned long> clashes;
  std::vector<utility::vector1<core::Size>> overlapping_nets;
  overlapping_nets.resize(network_vector_.size());

  for (auto i = network_vector_.begin(); i != network_vector_.end(); ++i) {
    core::Size const i_index = i - network_vector_.begin();
    utility::vector1<core::Size> i_global_rotamer_ids;
    i_global_rotamer_ids.reserve((*i)->residues.size());

    consider_network_as_branch_seed((*i), i_index);

    for (HBondResStructCOP res : (*i)->residues) {
      core::Size const moltenres =
          rotamer_sets_->resid_2_moltenres(res->resnum);
      core::Size const glboal_rot_id =
          rotamer_sets_->nrotamer_offset_for_moltenres(moltenres) +
          res->rot_index;
      i_global_rotamer_ids.push_back(glboal_rot_id);
    }

    std::sort(i_global_rotamer_ids.begin(), i_global_rotamer_ids.end());

    for (auto j = i; ++j != network_vector_.end();) {
      bool branch(false);
      if (is_sub_residues((*i)->residues, (*j)->residues, branch) || !branch) {
        continue;
      }
      core::Size const j_index = j - network_vector_.begin();

      if (!only_native_ && net_clash(**i, **j)) {
        debug_assert(i_index < j_index);
        clashes.insert(combine_indices(i_index, j_index));
        continue;
      }

      for (HBondResStructCOP res : (*j)->residues) {
        core::Size const moltenres =
            rotamer_sets_->resid_2_moltenres(res->resnum);
        core::Size const glboal_rot_id =
            rotamer_sets_->nrotamer_offset_for_moltenres(moltenres) +
            res->rot_index;
        if (std::binary_search(i_global_rotamer_ids.begin(),
                               i_global_rotamer_ids.end(), glboal_rot_id)) {
          overlapping_nets[i_index].push_back(j_index);
          overlapping_nets[j_index].push_back(i_index);
          break;
        }
      }
    }
  }

  if (branch_seed_network_vector_.size() == 0) {
    // use network_vector_ instead
  } else {
    for (core::Size ii = 1; ii <= branch_seed_network_vector_.size(); ++ii) {
      TR << "ii: " << ii << std::endl;
      monte_carlo_branching_trajectory(branch_seed_network_vector_[ii], clashes,
                                       overlapping_nets);
    }
  }
}

void HBNet::monte_carlo_branching_trajectory(
    core::Size seed_index, std::set<unsigned long> &clashes,
    std::vector<utility::vector1<core::Size>> &overlapping_nets) {
  utility::vector1<core::Size> current_rotamers;
  for (HBondResStructCOP res : network_vector_[seed_index]->residues) {
    core::Size const moltenres = rotamer_sets_->resid_2_moltenres(res->resnum);
    core::Size const glboal_rot_id =
        rotamer_sets_->nrotamer_offset_for_moltenres(moltenres) +
        res->rot_index;
    current_rotamers.push_back(glboal_rot_id);
  }
  std::sort(current_rotamers.begin(), current_rotamers.end());

  std::vector<core::Size> current_net_ids;
  current_net_ids.push_back(seed_index);

  utility::vector1<HBondResStructCOP> current_residues =
      network_vector_[seed_index]->residues;
  core::Real score = 0;
  for (core::Size ii = 1; ii < current_residues.size(); ++ii) {
    for (core::Size jj = ii + 1; jj <= current_residues.size(); ++jj) {
      score += get_twobody(current_residues[ii], current_residues[jj]);
    }
  }

  utility::vector1<core::Real> temps;
  temps.push_back(30);
  temps.push_back(5);
  temps.push_back(1);
  temps.push_back(0.1);
  temps.push_back(0);

  core::Size const loops_per_temp = 10;  // 50 total turns
  core::Size const candidates_per_loop =
      10;  // 500 total proposed changes, max of 50 accepted

  for (core::Real temp : temps) {
    for (core::Size inner_loop = 1; inner_loop <= loops_per_temp;
         ++inner_loop) {
      ////////////////////////////////////////////////////////
      // Step 1 - Determine which network additions to consider
      core::Size num_possible_neighbors = 0;
      for (core::Size index : current_net_ids) {
        num_possible_neighbors += overlapping_nets[index].size();
      }

      if (num_possible_neighbors == 0) {
        TR << "Ran out of neighbors with temp: " << temp
           << " and inner_loop: " << inner_loop << std::endl;
        return;
      }

      utility::vector1<core::Size> random_draws(candidates_per_loop, 0);
      for (core::Size ii = 1; ii <= candidates_per_loop; ++ii) {
        random_draws[ii] = static_cast<core::Size>(
            num_possible_neighbors * numeric::random::rg().uniform() + 1);
      }

      utility::vector1<HBondNetStructCOP> candidates;
      utility::vector1<core::Size> candidate_ids;
      candidates.reserve(candidates_per_loop);
      candidate_ids.reserve(candidates_per_loop);

      for (core::Size index : current_net_ids) {
        core::Size const num_nbrs = overlapping_nets[index].size();
        for (core::Size ii = 1; ii <= candidates_per_loop; ++ii) {
          if (random_draws[ii] == 0xFFFFFFFFFFFFFFFF) continue;
          if (random_draws[ii] <= num_nbrs) {
            core::Size ii_net_vec_index =
                overlapping_nets[index][random_draws[ii]];
            random_draws[ii] = 0xFFFFFFFFFFFFFFFF;
            HBondNetStructCOP ii_candidate = network_vector_[ii_net_vec_index];
            bool clash = false;
            for (core::Size index2 : current_net_ids) {
              if (clashes.find(combine_indices(index2, ii_net_vec_index)) !=
                  clashes.end()) {
                clash = true;
                break;
              }
            }
            if (!clash) {
              candidates.push_back(ii_candidate);
              candidate_ids.push_back(ii_net_vec_index);
            }
          } else {
            random_draws[ii] -= num_nbrs;
          }
        }
      }

      //////////////////////////
      // Step 2 - Measure DeltaEs
      utility::vector1<core::Real> deltaEs(candidates.size(), 0.0);
      utility::vector1<utility::vector1<HBondResStructCOP>> new_residues;
      new_residues.resize(candidates.size());
      for (core::Size ii = 1; ii <= candidates.size(); ++ii) {
        for (HBondResStructCOP res : candidates[ii]->residues) {
          core::Size const moltenres =
              rotamer_sets_->resid_2_moltenres(res->resnum);
          core::Size const global_rot_id =
              rotamer_sets_->nrotamer_offset_for_moltenres(moltenres) +
              res->rot_index;

          if (!std::binary_search(current_rotamers.begin(),
                                  current_rotamers.end(), global_rot_id)) {
            new_residues[ii].push_back(res);
            // score
            for (HBondResStructCOP existing_res : current_residues) {
              deltaEs[ii] += get_twobody(res, existing_res);
            }
          }
        }
      }

      ///////////////////////
      // Step 3 - Pick Options
      utility::vector1<core::Real> boltz_prob(deltaEs.size(), 0.0);
      core::Real partition_function = 0;
      for (core::Size ii = 1; ii <= deltaEs.size(); ++ii) {
        if (deltaEs[ii] < -5) {
          boltz_prob[ii] = std::exp(5.0 / temp);
        } else if (deltaEs[ii] > 5) {
          boltz_prob[ii] = 0;
        } else {
          boltz_prob[ii] = std::exp(-1 * deltaEs[ii] / temp);
        }
        partition_function += boltz_prob[ii];
      }
      TR << "partition_function: " << partition_function << std::endl;
      if (partition_function >= 0.01) {
        // really no good options
        continue;
      }

      core::Real random_draw = 0;
      if (partition_function < 1) {
        random_draw = numeric::random::rg().uniform();  // partition_function =
                                                        // 1
        if (random_draw > partition_function) {
          // don't accept anything
          continue;
        }
      } else {
        random_draw = partition_function * numeric::random::rg().uniform();
      }

      for (core::Size ii = 1; ii <= deltaEs.size(); ++ii) {
        random_draw -= boltz_prob[ii];
        if (random_draw <= 0) {
          // accept this position
          score += deltaEs[ii];
          current_net_ids.push_back(candidate_ids[ii]);
          for (HBondResStructCOP res : new_residues[ii]) {
            core::Size const moltenres =
                rotamer_sets_->resid_2_moltenres(res->resnum);
            core::Size const global_rot_id =
                rotamer_sets_->nrotamer_offset_for_moltenres(moltenres) +
                res->rot_index;
            current_rotamers.push_back(global_rot_id);
            current_residues.push_back(res);
          }
          break;
        }
      }
    }
  }

  merged_vecs_.push_back(current_net_ids);
}

core::Real HBNet::get_twobody(HBondResStructCOP i, HBondResStructCOP j) {
  core::Size const moltenres_i = rotamer_sets_->resid_2_moltenres(i->resnum);
  core::Size const moltenres_j = rotamer_sets_->resid_2_moltenres(j->resnum);
  if (moltenres_i == moltenres_j) {
    TR << "PROBLEM? get_twobody() was called with two rotamers on the same "
          "moltenres"
       << std::endl;
    if (i->rot_index == j->rot_index) {
      return 0;
    } else {
      return 9999999;
    }
  }

  NodeBase const *node = ig_->get_node(moltenres_i);
  EdgeBase const *edge = node->find_edge(int(moltenres_j));
  if (!edge) return 0;

  PDEdge const *pdedge = static_cast<PDEdge const *>(edge);
  TR << moltenres_i << " " << i->rot_index << " " << moltenres_j << " "
     << j->rot_index << std::endl;
  TR << node << " " << edge << " " << pdedge << std::endl;
  if (moltenres_i < moltenres_j) {
    return pdedge->get_two_body_energy(i->rot_index, j->rot_index);
  } else {
    return pdedge->get_two_body_energy(j->rot_index, i->rot_index);
  }
}

/*}
        runtime_assert( traverse_3mers_ );

        protocols::hbnet::NetworkRotamerIDTrie trie( 3 );

        for( core::Size ii=1; ii<=network_vector_.size(); ++ii ){
                //TR << "..." << network_vector_[ ii ]->residues.size() <<
std::endl;
                //runtime_assert( network_vector_[ ii ]->residues.size() == 3 );
                //TR << "..." << network_vector_[ ii ] << std::endl;
                if( !network_vector_[ ii ] ) continue;
                if( network_vector_[ ii ]->residues.size() != 3 ) continue;
                utility::vector1< core::Size > global_rotamer_ids( 3, 0 );
                for( core::Size index = 1; index <= 3; ++index ){
                        core::Size const mres =
rotamer_sets_->resid_2_moltenres( network_vector_[ ii ]->residues[ index
]->resnum );
                        global_rotamer_ids[ index ] =
rotamer_sets_->nrotamer_offset_for_moltenres( mres ) + network_vector_[ ii
]->residues[ index ]->rot_index;
                }
                trie.insert_all_permutations( global_rotamer_ids );
        }

        utility::vector1< utility::vector1< core::Size > >
global_rotamer_ids_for_starting_states;//TODO , populate
        utility::vector1< utility::vector1< core::Size > >
global_rotamer_ids_for_final_states;

        core::pack::hbnet::HBondNetworkGraphOP hbnetwork_graph( new
core::pack::hbnet::HBondNetworkGraph( *orig_pose_, scorefxn_,
hydrogen_bond_threshold_ ) );
        //TODO settings - such as interface

        //TODO make HBNetFragmentLibraryBase
        monte_carlo::monte_carlo_branch(
                global_rotamer_ids_for_starting_states,
                trie,
                global_rotamer_ids_for_final_states,
                ig_,
                hbnetwork_graph,
                rotamer_sets_,
                3//TODO make this a setting
                );

        network_vector_.clear();
        for( core::Size net = 1; net <=
global_rotamer_ids_for_final_states.size(); ++net ){
                //TODO pick up here
        }
}*/

// used by branch_overlapping() to efficiently search for all combinations of
// compatible networks that can be merged
void HBNet::rec_set_intersection(std::vector<Size> add_index_vec,
                                 std::vector<Size> next_index_vec, Size pos) {
  for (auto i = next_index_vec.begin(); i != next_index_vec.end(); ++i) {
    if (*i <= pos) {
      continue;
    }
    add_index_vec.push_back(*i);
    if (store_subnetworks_) {
      bool found(false);
      std::sort(add_index_vec.begin(), add_index_vec.end());
      for (auto &merged_vec : merged_vecs_) {
        // std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS
        // ALEADY BE SORTED!
        std::vector<int> v2(100);
        auto it2 = std::set_symmetric_difference(
            merged_vec.begin(), merged_vec.end(), add_index_vec.begin(),
            add_index_vec.end(), v2.begin());
        v2.resize(it2 - v2.begin());

        if (v2.size() == 0) {
          found = true;
          break;
        }
      }
      if (!found) {
        merged_vecs_.push_back(add_index_vec);
      }
    }
    std::sort(next_index_vec.begin(), next_index_vec.end());
    std::sort((network_vector_[*i])->net_indices.begin(),
              (network_vector_[*i])->net_indices.end());
    std::vector<Size> v(1000);
    std::vector<Size>::iterator it;
    it = std::set_intersection(next_index_vec.begin(), next_index_vec.end(),
                               (network_vector_[*i])->net_indices.begin(),
                               (network_vector_[*i])->net_indices.end(),
                               v.begin());
    v.resize(it - v.begin());
    if (v.size() == 0) {
      if (!store_subnetworks_) {
        bool found(false);
        std::sort(add_index_vec.begin(), add_index_vec.end());
        for (auto &merged_vec : merged_vecs_) {
          // std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS
          // ALEADY BE SORTED!
          std::vector<int> v2(100);
          auto it2 = std::set_symmetric_difference(
              merged_vec.begin(), merged_vec.end(), add_index_vec.begin(),
              add_index_vec.end(), v2.begin());
          v2.resize(it2 - v2.begin());

          if (v2.size() == 0) {
            found = true;
            break;
          }
        }
        if (!found) {
          merged_vecs_.push_back(add_index_vec);
        }
      }
      add_index_vec.pop_back();
      continue;
    } else if (v.size() == 1) {
      add_index_vec.push_back(v.front());
      bool found(false);
      std::sort(add_index_vec.begin(), add_index_vec.end());
      for (auto &merged_vec : merged_vecs_) {
        // std::sort(merged_vec.begin(), merged_vec.end()); // SHOULD ALWAYS
        // ALEADY BE SORTED!
        std::vector<int> v2(100);
        auto it2 = std::set_symmetric_difference(
            merged_vec.begin(), merged_vec.end(), add_index_vec.begin(),
            add_index_vec.end(), v2.begin());
        v2.resize(it2 - v2.begin());

        if (v2.size() == 0) {
          found = true;
          break;
        }
      }
      if (!found) {
        merged_vecs_.push_back(add_index_vec);
      }
      add_index_vec.pop_back();
    } else {
      rec_set_intersection(add_index_vec, v, *i);
    }
    add_index_vec.pop_back();
  }
}

void HBNet::place_rots_on_pose(pose::Pose &pose, hbond_net_struct &i,
                               bool use_pose) {
  if (!(pose.pdb_info())) {
    pose.pdb_info(core::pose::PDBInfoOP(new core::pose::PDBInfo(pose, true)));
  }
  if (!(i.rotamers.empty())) {
    for (const auto r : i.rotamers) {
      pose.replace_residue(r->seqpos(), *r, false);
      pose.pdb_info()->add_reslabel(r->seqpos(), "HBNet");
    }
  } else if (use_pose) {
    for (const auto &residue : i.residues) {
      pose.replace_residue(residue->resnum,
                           orig_pose_->residue(residue->resnum), false);
      pose.pdb_info()->add_reslabel(residue->resnum, "HBNet");
    }
  } else if (rotamer_sets_ != nullptr) {
    for (const auto &residue : i.residues) {
      ResidueCOP copy_rot(rotamer_sets_->rotamer_for_moltenres(
          rotamer_sets_->resid_2_moltenres((platform::uint)(residue->resnum)),
          residue->rot_index));
      pose.replace_residue(residue->resnum, *copy_rot, false);
      pose.pdb_info()->add_reslabel(residue->resnum, "HBNet");
    }
  }
  //    if ( place_waters ){
  //        place_bridging_waters_on_pose( pose, i );
  //        if ( pack_waters && !(i.bridging_water_oxygens.empty()) ){
  //            optimize_waters( pose, i );
  //        }
  //    }
  pose.update_residue_neighbors();
  // pose.update_pose_chains_from_pdb_chains();
}

// to very quickly eliminate networks with 1 or more unsatisfied heavy-atom
// donors or acceptors
bool HBNet::quick_and_dirty_network_has_heavy_atom_unsat(
    Pose const &pose, hbond_net_struct const &network) {
  for (const auto &res : network.residues) {
    if (res_is_core(res->resnum)) {
      for (const auto &acc_pos : pose.residue(res->resnum).accpt_pos_sc()) {
        // Size base( pose.residue(res->resnum) );
        if (quick_and_dirty_heavy_atom_is_unsat(
                pose, id::AtomID(acc_pos, res->resnum))) {
          return true;
        }
      }
      // for ( const auto & chemical::AtomIndices :
      // pose.residue(res->resnum).Hpos_polar_sc() ){
      Size prev_base(0);
      for (const auto &hpol : pose.residue(res->resnum).Hpos_polar_sc()) {
        Size base(pose.residue(res->resnum).atom_base(hpol));
        if (prev_base == base) continue;
        if (quick_and_dirty_heavy_atom_is_unsat(
                pose, id::AtomID(base, res->resnum))) {
          return true;
        }
        prev_base = base;
      }
    }
  }
  return false;
}

// to very quickly eliminate networks with 1 or more unsatisfied heavy-atom
// donors or acceptors
bool HBNet::quick_and_dirty_heavy_atom_is_unsat(Pose const &pose,
                                                id::AtomID const at_id) {
  Real const heavy_atom_don_acc_cutoff(3.5);
  Real const dist_2(heavy_atom_don_acc_cutoff * heavy_atom_don_acc_cutoff);
  for (Size res = 1; res <= pose.total_residue(); ++res) {
    if (res == at_id.rsd()) continue;
    if (pose.xyz(at_id).distance_squared(pose.residue(res).nbr_atom_xyz()) <
        100.0) {  // 10^2
      if (pose.residue(at_id.rsd()).atom_type(at_id.atomno()).is_acceptor()) {
        Size prev_base(0);
        for (const auto &hpol : pose.residue(res).Hpos_polar()) {
          Size base(pose.residue(res).atom_base(hpol));
          if (prev_base == base) continue;
          if (pose.xyz(at_id).distance_squared(pose.residue(res).xyz(base)) <
              dist_2) {
            return false;
          }
          prev_base = base;
        }
      }
      if (pose.residue(at_id.rsd()).atom_type(at_id.atomno()).is_donor()) {
        for (const auto &acc : pose.residue(res).accpt_pos()) {
          if (pose.xyz(at_id).distance_squared(pose.residue(res).xyz(acc)) <
              dist_2) {
            return false;
          }
        }
      }
    }
  }
  return true;
}

// hacky check based on distance for cases where we only have water O (not
// explicit H's)
bool HBNet::atom_hbonds_to_bridging_water(Pose const &pose,
                                          id::AtomID const at_id) {
  // TR << "checking atom_hbonds_to_bridging_water" << std::endl;
  Real const hpol_dist(2.2);  // water O (acc) to Hpol dist: generous, error on
                              // side of satisfaction
  Real const hpol_dist_sq(hpol_dist * hpol_dist);
  Real const acc_dist(3.4);  // water O (don) to Acc dist: generous, error on
                             // side of satisfaction
  Real const acc_dist_2(acc_dist * acc_dist);
  for (Size res = 1; res <= pose.total_residue(); ++res) {
    // if ( pose.residue(res).is_water() && pose.xyz( at_id ).distance_squared(
    // pose.residue( res ).nbr_atom_xyz()) < 16.0 ) { // 4^2; all water types
    // have neighbor radius < 4
    if (pose.residue(res).name1() == 'w' &&
        pose.xyz(at_id).distance_squared(pose.residue(res).nbr_atom_xyz()) <
            16.0) {  // 4^2; all water types have neighbor radius < 4
      // TR << "we found a water! " << std::endl;
      Vector waterO(
          pose.residue(res).xyz(1));  // water O should always be first;
      // Vector waterO( pose.residue(res).xyz("O") );
      if (pose.residue(at_id.rsd()).atom_type(at_id.atomno()).is_acceptor() &&
          pose.xyz(at_id).distance_squared(waterO) < acc_dist_2) {
        // TR << "it H-BONDS!" << std::endl;
        return true;
      } else if (pose.residue(at_id.rsd())
                     .atom_is_polar_hydrogen(at_id.atomno()) &&
                 pose.xyz(at_id).distance_squared(waterO) < hpol_dist_sq) {
        // TR << "it H-BONDS!" << std::endl;
        return true;
      }
    }
  }
  // TR << "it DOES NOT h-bond" << std::endl;
  return false;
}

void HBNet::update_core_and_boundary_residues(core::pose::Pose const &pose) {
  if (core_selector_) {
    core_residues_ = core_selector_->apply(pose);
  } else {
    core::select::residue_selector::LayerSelectorOP core_layer(
        new core::select::residue_selector::LayerSelector());
    core_layer->set_layers(true, false, false);
    // core_layer->set_use_sc_neighbors( true ); // now true by default
    core_layer->set_cutoffs(4.4 /* core */,
                            2.0 /* surface */);  // default core is 5.2, here we
                                                 // use 4.4 to be overly
                                                 // cautious with unsats
    runtime_assert(core_layer->use_sc_neighbors());
    core_residues_ = core_layer->apply(pose);
  }
}

void HBNet::find_unsats(Pose const &pose, hbond_net_struct &network) {
  Size num_unsatisfied(0);
  Size num_heavy_unsat(0);

  // update core residues
  // core::select::residue_selector::ResidueSubset temp_core_residues =
  // get_core_residues();
  // core::select::residue_selector::ResidueSubset temp_boundary_residues =
  // get_boundary_residues();
  // if ( bridging_waters_ ){
  //    update_core_and_boundary_residues( pose );
  //}

  runtime_assert(!(network.hbond_vec.empty()));

  if (verbose_ && TR.visible()) network.hbond_set->show(TR);

  // TODO NEED TO ADD option to use own bunsat calc here

  Size polar_atom_count(0);
  Size polar_atms_making_hbonds(0);
  std::vector<Size> resnums(0);
  for (utility::vector1<HBondResStructCOP>::const_iterator res =
           network.residues.begin();
       res != network.residues.end(); ++res) {
    // for ( std::set< Size >::const_iterator res =
    // actual_hbond_residues.begin(); res != actual_hbond_residues.end(); ++res
    // ){
    resnums.push_back((*res)->resnum);
  }
  //    if ( bridging_waters_ && !(network.waterrots.empty()) ){
  //        Size anchor_res( (symmetric_) ?
  //        symm_info_->num_independent_residues() : orig_pose_->total_residue()
  //        );
  //        Size new_total( pose.total_residue() );
  //        if ( symmetric_ ){
  //            core::conformation::symmetry::SymmetricConformation const &
  //            newSymmConf(dynamic_cast<core::conformation::symmetry::SymmetricConformation
  //            const & > ( pose.conformation()));
  //            new_total =
  //            newSymmConf.Symmetry_Info()->num_independent_residues();
  //        }
  //        for ( Size r = anchor_res; r <= new_total; r++ ){
  //            if ( pose.residue(r).is_water() &&
  //            pose.pdb_info()->res_haslabel(r, "HBNet" ) ){
  //                resnums.push_back(r);
  //            }
  //        }
  //    }

  for (const auto &resnum : resnums) {
    if (verbose_ && TR.visible())
      network.hbond_set->show(pose, resnum, true, TR);

    bool skip_backbones(true);
    // bool skip_ligands(true);
    //        if ( !(pose.residue( resnum ).is_water()) ) { // need to handle
    //        waters separately
    for (Size a_index = 1; a_index <= pose.residue(resnum).natoms();
         ++a_index) {
      if ((skip_backbones &&
           !(pose.residue(resnum).atom_is_backbone(a_index))) ||
          ((this->ligand()) && resnum == (this->ligand()))) {
        if (pose.residue(resnum).atom_type(a_index).is_donor() &&
            pose.residue(resnum).atomic_charge(a_index) != 0.0 &&
            pose.residue(resnum).atom_type(a_index).name() != "OH") {
          Size h_count(0);
          Size h_unsat(0);
          for (Size hatm = pose.residue(resnum).attached_H_begin(a_index);
               hatm <= pose.residue(resnum).attached_H_end(a_index); ++hatm) {
            // TR.Debug << "hatm = " <<
            // pose.residue(resnum).atom_type(hatm).atom_type_name() <<
            // std::endl;
            h_count++;
            polar_atom_count++;
            id::AtomID const at_id(hatm, resnum);

            if (!((network.hbond_set)
                      ->nhbonds(at_id, false /* include_only_allowed */))) {
              // if ( !((network.hbond_set)->nhbonds( at_id, false /*
              // include_only_allowed */ )) &&
              //( !bridging_waters_ || !network.has_bridging_wat ||
              //!(network.waterrots.empty()) || !atom_hbonds_to_bridging_water(
              // pose, at_id ) ) ){ // if waterrots !empty, no need to check
              // this
              // because waters optimized and should be in hbond_vec
              // if ( ( !((network.hbond_set)->nhbonds( at_id, false /*
              // include_only_allowed */ )) && ( !bridging_waters_ ||
              // !(network.has_bridging_wat) || !(network.waterrots.empty()) ) )
              // ||
              //( bridging_waters_ && network.has_bridging_wat &&
              // network.waterrots.empty() && !atom_hbonds_to_bridging_water(
              // pose, at_id ) ) ){
              if (res_is_core(resnum)) {
                num_unsatisfied++;
                h_unsat++;
                network.unsat_Hpols.push_back(at_id);
                TR.Debug << " res " << resnum << " atom " << hatm << " "
                         << pose.residue(resnum).atom_name(hatm) << " is unsat"
                         << std::endl;
              }
            } else {
              polar_atms_making_hbonds++;
            }
          }
          if (h_unsat == h_count) {
            num_heavy_unsat++;
            if (verbose_ && TR.visible())
              TR << "h_unsat == h_count; heavy unsat is residue " << resnum
                 << pose.residue(resnum).name1() << " atom " << a_index
                 << std::endl;
          }
        } else if (pose.residue(resnum).atom_type(a_index).is_acceptor() &&
                   pose.residue(resnum).atomic_charge(a_index) !=
                       0.0) {  // important for ligand case
          polar_atom_count++;
          id::AtomID const at_id(a_index, resnum);

          if (!((network.hbond_set)
                    ->nhbonds(at_id, false /* include_only_allowed */))) {
            //                        if ( !((network.hbond_set)->nhbonds(
            //                        at_id, false /* include_only_allowed */ ))
            //                            && ( !bridging_waters_ ||
            //                            !network.has_bridging_wat ||
            //                            !(network.waterrots.empty()) ||
            //                            !atom_hbonds_to_bridging_water( pose,
            //                            at_id ) ) ){ // if waterrots !empty,
            //                            no need to check this because waters
            //                            optimized and should be in hbond_vec
            // if ( ( !( (network.hbond_set)->nhbonds( at_id , false /*
            // include_only_allowed */ )) && ( !bridging_waters_ ||
            // !(network.has_bridging_wat) || !(network.waterrots.empty()) ) )
            // ||
            //( bridging_waters_ && network.has_bridging_wat &&
            // network.waterrots.empty() && !atom_hbonds_to_bridging_water(
            // pose,
            // at_id ) ) ){
            // if res is buried (aka "core") then we penalize it; if bridging
            // waters: expicits checked below, but if don't have them yet, check
            // based on water O
            if (res_is_core(resnum)) {
              if (pose.residue(resnum).atom_type(a_index).name() ==
                  "OH") {  // do not want to double penalize hydroxyls
                polar_atom_count++;
                core::Size hatm =
                    pose.residue(resnum).attached_H_begin(a_index);
                if (!((network.hbond_set)
                          ->nhbonds(id::AtomID(hatm, resnum),
                                    false /* include_only_allowed */))) {
                  //                                    if (
                  //                                    !((network.hbond_set)->nhbonds(
                  //                                    id::AtomID( hatm, resnum
                  //                                    ), false /*
                  //                                    include_only_allowed */
                  //                                    )) &&
                  //                                        ( !bridging_waters_
                  //                                        ||
                  //                                        !network.has_bridging_wat
                  //                                        ||
                  //                                        !(network.waterrots.empty())
                  //                                        ||
                  //                                        !atom_hbonds_to_bridging_water(
                  //                                        pose, id::AtomID(
                  //                                        hatm, resnum ) ) )
                  //                                        ){ // if waterrots
                  //                                        !empty, no need to
                  //                                        check this
                  // if ( ( !((network.hbond_set)->nhbonds( id::AtomID( hatm,
                  // resnum ), false /* include_only_allowed */ )) &&
                  //( !bridging_waters_ || !(network.has_bridging_wat) ||
                  //!(network.waterrots.empty()) ) ) ||
                  //( bridging_waters_ && network.has_bridging_wat &&
                  // network.waterrots.empty() &&
                  // !atom_hbonds_to_bridging_water(
                  // pose, id::AtomID( hatm, resnum ) ) ) ){
                  num_unsatisfied++;
                  num_heavy_unsat++;
                  if (verbose_ && TR.visible())
                    TR << "heavy unsat is residue " << at_id.rsd()
                       << pose.residue(at_id.rsd()).name1() << " atom "
                       << a_index << std::endl;
                  network.unsat_Hpols.push_back(at_id);
                } else {
                  polar_atms_making_hbonds++;
                }
              } else {
                num_unsatisfied++;
                num_heavy_unsat++;
                if (verbose_ && TR.visible())
                  TR << "heavy unsat is residue " << at_id.rsd()
                     << pose.residue(at_id.rsd()).name1() << " atom "
                     << at_id.atomno() << std::endl;
                network.unsat_accs.push_back(at_id);
                TR.Debug << " res " << resnum << " atom " << a_index << " "
                         << pose.residue(resnum).atom_name(a_index)
                         << " is heavy unsat" << std::endl;
                TR.Debug << " res " << resnum << " atom " << a_index << " "
                         << pose.residue(resnum).atom_name(a_index)
                         << " is unsat" << std::endl;
                // TR << " res_is_core( resnum ) = " << res_is_core( resnum ) <<
                // " and atom_sasa[ at_id ] = " << atom_sasa[ at_id ] <<
                // std::endl;
              }
            }
          } else {
            polar_atms_making_hbonds++;

            if (pose.residue(resnum).atom_type(a_index).name() == "OH") {
              polar_atom_count++;
              core::Size hatm = pose.residue(resnum).attached_H_begin(a_index);
              if (!((network.hbond_set)
                        ->nhbonds(id::AtomID(hatm, resnum),
                                  false /* include_only_allowed */))) {
                //                                if (
                //                                !((network.hbond_set)->nhbonds(
                //                                id::AtomID( hatm, resnum ),
                //                                false /* include_only_allowed
                //                                */ )) &&
                //                                    ( !bridging_waters_ ||
                //                                    !network.has_bridging_wat
                //                                    ||
                //                                    !(network.waterrots.empty())
                //                                    ||
                //                                    !atom_hbonds_to_bridging_water(
                //                                    pose, id::AtomID( hatm,
                //                                    resnum ) ) ) ){ // if
                //                                    waterrots !empty, no need
                //                                    to check this
                if (res_is_core(resnum) &&
                    (hydroxyls_must_donate_ ||
                     (pose.residue(resnum).name1() == 'Y' &&
                      tyr_hydroxyls_must_donate_))) {
                  num_unsatisfied++;
                  num_heavy_unsat++;
                  if (verbose_ && TR.visible())
                    TR << "heavy unsat is residue " << at_id.rsd()
                       << pose.residue(at_id.rsd()).name1() << " atom "
                       << a_index << std::endl;
                  network.unsat_Hpols.push_back(at_id);
                  TR.Debug << "res " << resnum << " atom " << a_index << " "
                           << pose.residue(resnum).atom_name(a_index)
                           << " is heavy unsat" << std::endl;
                }
              } else {
                polar_atms_making_hbonds++;
              }
            }
          }
        }
      }
    }
    //        }
    //        else if ( pose.residue(resnum).is_water() &&
    //        !(network.waterrots.empty()) ) { // need a more clear check:
    //        !(network.waterrots.empty() means waters have been packed
    //            //water_count++;
    //            // for waters, require at least 2 h-bonds:
    //            //    this works for current water types, but need more
    //            general solution
    //            Size a_index( 1 );
    //            Size num_hbonds_for_water(0);
    //            polar_atom_count++;
    //            //TR << "water count = " << water_count << std::endl;
    //            //TR << "number water O h-bonds" <<
    //            (network.hbond_set)->nhbonds( id::AtomID( a_index, resnum ),
    //            false ) << std::endl;
    //            if ( !((network.hbond_set)->nhbonds( id::AtomID( a_index,
    //            resnum ), false /* include_only_allowed */ )) ){
    //                num_unsatisfied++;
    //            } else {
    //                num_hbonds_for_water++;
    //                polar_atms_making_hbonds++;
    //            }
    //            polar_atom_count++;
    //            core::Size hatm1 = pose.residue(resnum).attached_H_begin(
    //            a_index );
    //            //TR << "number hatm1 h-bonds" <<
    //            (network.hbond_set)->nhbonds( id::AtomID( hatm1, resnum ),
    //            false ) << std::endl;
    //            //Size hatm1( 2 ); // this should work?  atoms numbered
    //            starting with 1; not sure if attached_H works for waters yet??
    //            if ( !((network.hbond_set)->nhbonds( id::AtomID( hatm1, resnum
    //            ), false /* include_only_allowed */ )) ){
    //                num_unsatisfied++;
    //            } else {
    //                num_hbonds_for_water++;
    //                polar_atms_making_hbonds++;
    //            }
    //            polar_atom_count++;
    //            core::Size hatm2 = pose.residue(resnum).attached_H_end(
    //            a_index );
    //            //TR << "number hatm2 h-bonds" <<
    //            (network.hbond_set)->nhbonds( id::AtomID( hatm2, resnum ),
    //            false ) << std::endl;
    //            //Size hatm2( 3 ); // this should work?  atoms numbered
    //            starting with 1; not sure if attached_H works for waters yet??
    //            if ( !((network.hbond_set)->nhbonds( id::AtomID( hatm2, resnum
    //            ), false /* include_only_allowed */ ))) {
    //                num_unsatisfied++;
    //            } else {
    //                num_hbonds_for_water++;
    //                polar_atms_making_hbonds++;
    //            }
    //            if ( num_hbonds_for_water < 2 ) num_heavy_unsat++;
    //        }
  }
  network.num_unsat = num_unsatisfied;
  network.num_heavy_unsat = num_heavy_unsat;
  network.connectivity = (Real)polar_atms_making_hbonds /
                         (Real)polar_atom_count;  // NOW UPDATED FOR BRIDGING
                                                  // WATERS (CODE ADDED ABOVE)

  // if ( bridging_waters_ ){
  //    set_core_residues( temp_core_residues );
  //    set_boundary_residues( temp_boundary_residues );
  //}
}  // find_unsats

// checks resnum and aa but not rotamer
// consider Ser and Thr identical
bool HBNet::residues_identical(utility::vector1<HBondResStructCOP> &residues1,
                               utility::vector1<HBondResStructCOP> &residues2) {
  if (residues1.size() != residues2.size()) {
    return false;
  }
  std::sort(residues1.begin(), residues1.end(), compare_hbond_residues());
  std::sort(residues2.begin(), residues2.end(), compare_hbond_residues());

  utility::vector1<HBondResStructCOP>::const_iterator res2 = residues2.begin();
  for (utility::vector1<HBondResStructCOP>::const_iterator res1 =
           residues1.begin();
       res1 != residues1.end(); ++res1) {
    if (((*res1)->resnum == (*res2)->resnum) &&
        ((((*res1)->aa == (*res2)->aa)) ||
         (((*res1)->aa == 'S' || (*res1)->aa == 'T') &&
          ((*res2)->aa == 'S' || (*res2)->aa == 'T')))) {
      ++res2;
      continue;
    } else {
      return false;
    }
  }

  const bool check_chis = basic::options::option
      [basic::options::OptionKeys::jackmag::hbnet_check_chis]();
  if (check_chis && !all_residue_chis_are_close(residues1, residues2)) {
    return false;
  }

  return true;
}

bool HBNet::all_residue_chis_are_close(
    utility::vector1<HBondResStructCOP> &residues1,
    utility::vector1<HBondResStructCOP> &residues2) {
  const core::Size window = 20;  // degress

  if (residues1.size() != residues2.size()) {
    return false;
  }
  std::sort(residues1.begin(), residues1.end(), compare_hbond_residues());
  std::sort(residues2.begin(), residues2.end(), compare_hbond_residues());

  utility::vector1<HBondResStructCOP>::const_iterator res2 = residues2.begin();
  for (utility::vector1<HBondResStructCOP>::const_iterator
           res1 = residues1.begin();
       res1 != residues1.end(); ++res1, ++res2) {
    core::conformation::ResidueCOP res1_cop =
        rotamer_sets_->rotamer_set_for_residue((*res1)->resnum)
            ->rotamer((*res1)->rot_index);
    core::conformation::ResidueCOP res2_cop =
        rotamer_sets_->rotamer_set_for_residue((*res2)->resnum)
            ->rotamer((*res2)->rot_index);
    if (res1_cop->nchi() != res2_cop->nchi()) return false;

    for (core::Size chi = 1; chi <= res1_cop->nchi(); ++chi) {
      core::Real difference = std::abs(res1_cop->chi(chi) - res2_cop->chi(chi));

      while (difference > 360) difference -= 360;
      if (difference > 180) difference = 360 - difference;

      if (difference > window) return false;
    }
  }
  return true;
}

bool HBNet::residues_not_unique(
    utility::vector1<HBondResStructCOP> &residues1,
    utility::vector1<HBondResStructCOP> &residues2) {
  if (residues1.size() != residues2.size()) {
    return false;
  }
  std::sort(residues1.begin(), residues1.end(), compare_hbond_residues());
  std::sort(residues2.begin(), residues2.end(), compare_hbond_residues());

  utility::vector1<HBondResStructCOP>::const_iterator res2 = residues2.begin();
  for (utility::vector1<HBondResStructCOP>::const_iterator res1 =
           residues1.begin();
       res1 != residues1.end(); ++res1) {
    if (((*res1)->resnum == (*res2)->resnum) &&
        ((((*res1)->aa == (*res2)->aa)) ||
         (((*res1)->aa == 'N' || (*res1)->aa == 'Q') &&
          ((*res2)->aa == 'N' || (*res2)->aa == 'Q')) ||
         (((*res1)->aa == 'N' || (*res1)->aa == 'D') &&
          ((*res2)->aa == 'N' || (*res2)->aa == 'D')) ||
         (((*res1)->aa == 'Q' || (*res1)->aa == 'E') &&
          ((*res2)->aa == 'Q' || (*res2)->aa == 'E')) ||
         (((*res1)->aa == 'D' || (*res1)->aa == 'E') &&
          ((*res2)->aa == 'D' || (*res2)->aa == 'E')) ||
         (((*res1)->aa == 'S' || (*res1)->aa == 'T') &&
          ((*res2)->aa == 'S' || (*res2)->aa == 'T')))) {
      ++res2;
      continue;
    } else {
      return false;
    }
  }
  return true;
}

bool HBNet::networks_unique(hbond_net_struct const &i,
                            hbond_net_struct const &j,
                            bool no_surface /* true */) {
  utility::vector1<HBondResStructCOP> residues_i(0);
  utility::vector1<HBondResStructCOP> residues_j(0);

  if (no_surface) {
    for (const auto &residue : i.residues) {
      if (res_is_core(residue->resnum) || res_is_boundary(residue->resnum)) {
        residues_i.push_back(residue);
      }
    }
    for (const auto &residue : j.residues) {
      if (res_is_core(residue->resnum) || res_is_boundary(residue->resnum)) {
        residues_j.push_back(residue);
      }
    }
  } else {
    residues_i = i.residues;
    residues_j = j.residues;
  }
  return !(residues_not_unique(residues_i, residues_j));
}

// 02/28/15 changing behavior, will only keep 1 network with unique resnum/aa
// identity; networks with same resnum/aa seq but unique rotamers
// considered through final scoring
void HBNet::remove_replicate_networks(Size same_max /*=1*/) {
  std::sort(network_vector_.begin(), network_vector_.end(),
            compare_net_vec());  // sort all networks by score
  if (basic::options::option
          [basic::options::OptionKeys::jackmag::bypass_remove_rep]()) {
    return;
  }

  std::vector<HBondNetStructOP> temp_net_vec(0);
  // std::vector< Size > net_count(0);
  if (verbose() && TR.visible())
    TR << "         DELETING REPLICATE NETWORKS THAT HAVE IDENTIAL RESNUM AND "
          "AA SEQ; ONLY KEEPING BEST SCORING ONE"
       << std::endl;
  for (auto &i : network_vector_) {
    bool only_water(true);
    Size num_protein_res(0);
    for (utility::vector1<HBondResStructCOP>::const_iterator k =
             i->residues.begin();
         k != i->residues.end(); ++k) {
      if ((*k)->is_protein || (*k)->is_ligand) {
        num_protein_res++;
        if (num_protein_res >= 2) {
          only_water = false;
          break;
        }
      }
    }
    if (only_water) {
      continue;
    }
    bool reached_limit = false;

    // Size temp_net_count = 0;
    Size temp_net_count = 1;
    for (auto vit = temp_net_vec.begin(); vit != temp_net_vec.end(); ++vit) {
      if (residues_identical(i->residues,
                             (*vit)->residues)) {  //&& ( !bridging_waters_ || (
        // i->bridging_water_oxygens.size() ==
        //(*vit)->bridging_water_oxygens.size() ) )
        //){
        // net_count[temp_net_count] = net_count[temp_net_count] + 1;
        temp_net_count++;
        // if ( net_count[temp_net_count] > same_max ) {
        if (temp_net_count > same_max) {
          reached_limit = true;
          // don't want break here because need to go through all
        }
      } else if (only_native_ &&
                 is_sub_residues(
                     i->residues,
                     (*vit)->residues)) {  // checks subsets both ways
        // for native case, ensure that we don't store subnetworks
        if (i->residues.size() > (*vit)->residues.size()) {
          **vit = *i;
          auto vit2 = vit + 1;
          for (; vit2 != temp_net_vec.end();) {
            if (is_sub_residues(i->residues, (*vit)->residues)) {
              vit2 = temp_net_vec.erase(vit2);
            } else {
              ++vit2;
            }
          }
        }
        reached_limit = true;
        break;
      }
    }
    if (!reached_limit) {
      temp_net_vec.push_back(i);
    } else if (verbose_ && TR.visible()) {
      TR << "reached limit, continuing: network " << print_list_to_string(*i)
         << std::endl;
    }
  }
  // clear network_vector_ then push back all of the ones that passed the cutoff
  // (from temp_net_vec)
  network_vector_.clear();
  for (auto &k : temp_net_vec) {
    network_vector_.push_back(k);
  }
  if (verbose_ && TR.visible()) {
    TR << " AFTER REMOVING REPLICATES: " << network_vector_.size()
       << " NETWORKS" << std::endl;
    TR << " ========================================================"
       << std::endl;
  }
}  // remove_replicate_networks

Size HBNet::get_num_native_rot(
    Pose &pose, utility::vector1<HBondResStructCOP> const &residues,
    Real sc_rmsd_cut, bool super) {
  Size num_native(0);
  for (const auto &residue : residues) {
    if (((this->ligand()) && residue->resnum == (this->ligand())) ||
        residue->resnum > orig_pose_->total_residue()) {
      continue;  // does not count ligand in order to be fair
    }
    if (!(pose.residue(residue->resnum).name1() ==
          orig_pose_->residue(residue->resnum).name1())) {
      if ((pose.residue(residue->resnum).name1() != 'S' &&
           pose.residue(residue->resnum).name1() != 'T') ||
          (orig_pose_->residue(residue->resnum).name1() != 'S' &&
           orig_pose_->residue(residue->resnum).name1() != 'T')) {
        continue;
      } else {
        num_native++;  // for benchmarking purposes, treat S/T same
      }
      continue;
    }
    // Real sc_rmsd( core::scoring::automorphic_rmsd(
    // pose.residue((*res)->resnum), orig_pose_->residue((*res)->resnum), super
    // /*superpose*/ ) );
    Real sc_rmsd(core::scoring::residue_sc_rmsd_no_super(
        pose.residue(residue->resnum).clone(),
        orig_pose_->residue(residue->resnum).clone(),
        super /*final group only*/));
    if (sc_rmsd <= sc_rmsd_cut) {
      num_native++;
    }
  }
  return num_native;
}

Size HBNet::get_num_native_seq(
    core::pose::Pose &pose,
    utility::vector1<HBondResStructCOP> const &residues) {
  Size num_native(0);
  for (const auto &residue : residues) {
    if (((this->ligand()) && residue->resnum == (this->ligand())) ||
        residue->resnum >= orig_pose_->total_residue()) {
      continue;  // does not count ligand in ligand() case
    }
    bool aa_same = false;
    if (pose.residue(residue->resnum).name1() ==
        orig_pose_->residue(residue->resnum).name1()) {
      aa_same = true;
    } else if ((pose.residue(residue->resnum).name1() == 'S' ||
                pose.residue(residue->resnum).name1() == 'T') &&
               (orig_pose_->residue(residue->resnum).name1() == 'S' ||
                orig_pose_->residue(residue->resnum).name1() == 'T')) {
      aa_same = true;
    }
    if (aa_same) {
      num_native++;
    }
  }
  return num_native;
}

void HBNet::output_networks(bool finalize) {
  if (finalize) {
    // output_net_vec_.clear(); //if going to finalize, clear output_net_vec
    output_vector_.clear();
  }
  if (TR.visible()) {
    TR << " Designed these new networks that meet your criteria: " << std::endl
       << "\t" << print_headers() << print_additional_headers() << std::endl;
    TR << "network_vector_.size() == " << network_vector_.size() << std::endl;
  }

  if (jack_store_networks_) {
    network_residues_.resize(network_vector_.size() + native_networks_.size());
  }
  std::sort(network_vector_.begin(), network_vector_.end(),
            compare_net_vec());  // sort all new networks
  Size count(1);
  // for ( auto i = network_vector_.begin(); i != network_vector_.end(); ++i ) {
  // // outer loop
  auto i = network_vector_.begin();
  for (; i != network_vector_.end();) {  // outer loop
    // for this options case, need to check if existing (native) networks have
    // been extended
    //  if they have, choose either native or extended, whichever is better
    if (extend_existing_networks_) {
      bool erase_this_network = false;
      auto nit(native_networks_.begin());
      for (; nit != native_networks_.end();) {  // inner loop
        // bool dummy(false);
        // if the native network is a subset of network i, then i is extended;
        // if ( is_sub_residues( (*nit)->residues, (*i)->residues, dummy, true
        // /* return true_if_identical (networks same size) */) ){
        if (is_sub_residues((*nit)->residues, (*i)->residues)) {
          // if extended is bigger (and meets criteria) or same size but few
          // unsats, then erase native
          //   bool operator< overridden by hbond_net_struct to define best
          //   network
          Size i_size((symmetric_ && !((*i)->asymm_residues.empty()))
                          ? (*i)->asymm_residues.size()
                          : (*i)->residues.size());
          Size nit_size((symmetric_ && !(*nit)->asymm_residues.empty())
                            ? (*nit)->asymm_residues.size()
                            : (*nit)->residues.size());
          if (i_size > nit_size ||
              (i_size == nit_size && (*i)->num_unsat < (*nit)->num_unsat)) {
            (*i)->is_extended =
                true;  // record that this is an extended network
            nit = native_networks_.erase(nit);
          } else {
            // keep native and don't keep new network; we won't add it to
            // output_net_vec_
            // skip this network (it won't get returned)
            erase_this_network = true;
            break;  // because we found a native that's better, exit inner loop
          }
        } else {
          ++nit;
        }
      }
      // logic for determing if network is native or extended from staring HBNet
      // PDBInfoLabels:
      if (erase_this_network) {
        i = network_vector_.erase(i);
        continue;
      }  // move on to the next new network, continue outer loop

      // 17/02/15 shouldn't need this anymore
      //            Size hbnet_label_count(0);
      //            for (auto & residue : (*i)->residues) {
      //                if ( residue->resnum <=
      //                input_hbnet_info_residues_.size() &&
      //                input_hbnet_info_residues_[residue->resnum] ){
      //                    hbnet_label_count++;
      //                }
      //            }
      //            if ( hbnet_label_count == (*i)->residues.size() )
      //                (*i)->is_native = true;
      //            else if ( hbnet_label_count > 2 )
      //                (*i)->is_extended = true;
    }

    Pose ala_copy = *ala_pose_;  // deep copy
    // place_rots_on_pose( ala_copy, **i, (*i)->is_native, bridging_waters_,
    // (*i)->waterrots.empty() );
    place_rots_on_pose(ala_copy, **i, (*i)->is_native);
    if ((*i)->hbond_vec.empty()) {
      (*scorefxn_)(ala_copy);  // NEED TO SCORE BEFORE get_hbond_atom_pairs()
      get_hbond_atom_pairs(**i, ala_copy);
      (*i)->total_hbonds = (*i)->hbond_vec.size();
    }
    (*i)->id = count;
    std::string outstring(
        (pdb_numbering()) ? print_network_w_pdb_numbering(ala_copy, **i, true)
                          : print_network(**i));
    outstring = outstring + this->print_additional_info_for_net(**i, ala_copy);
    (*i)->outstring = outstring;

    if (write_network_pdbs_) write_network_pdb(*i);
    if (jack_store_networks_) jack_store_network(*i, false);

    if (TR.visible()) TR << "\t" << outstring << std::endl;
    if (finalize) {
      // std::vector< Size > net_ids(0);
      // net_ids.push_back( (*i)->id );
      std::set<Size> net_ids;
      net_ids.insert((*i)->id);

      if (keep_existing_networks_) {
        if ((*i)->is_extended) {
          for (auto j = i; ++j != network_vector_.end();) {
            if ((*j)->is_extended) {
              bool compatible(true);
              for (auto &net_id : net_ids) {
                if (net_clash(*(get_network_by_id(net_id)), **j)) {
                  compatible = false;
                  break;
                }
              }
              if (compatible) net_ids.insert((*j)->id);
            }
          }
        } else {
          for (auto &k : network_vector_) {
            if (k->is_extended) {
              bool compatible(true);
              for (const auto &net_id : net_ids) {
                if (net_clash(*(get_network_by_id(net_id)), *k)) {
                  compatible = false;
                  break;
                }
              }
              if (compatible) net_ids.insert(k->id);
            }
          }
        }
      }
      bool already_stored(false);
      for (const auto &output_vec : output_vector_) {
        if (output_vec == net_ids) {
          already_stored = true;
        }
      }
      if (!already_stored) {
        output_vector_.push_back(net_ids);
      }
      // output_net_vec_.push_back( *i );
    }
    count++;
    i++;
  }
  if (find_native_ || keep_existing_networks_) {
    if (TR.visible() && native_networks_.size() > 0) {
      TR << std::endl
         << " Keeping these existing (native) networks that meet your "
            "criteria: "
         << std::endl
         << "\t" << print_headers() << print_additional_headers() << std::endl;
    }
    std::sort(native_networks_.begin(), native_networks_.end(),
              compare_net_vec());  // sort native networks
    Size count(1);
    std::string comment_str(print_headers() + this->print_additional_headers());
    Pose &out_pose((output_poly_ala_background_) ? *ala_pose_
                                                 : nonconst_get_orig_pose());
    for (auto &native_network : native_networks_) {
      native_network->id = count;

      if (jack_store_networks_) {
        jack_store_network(native_network, true);
      }

      // if these options are set, we want to retain the native networks and
      // csts in each pose that HBNet passes back
      //   call private nonconst_get_orig_pose() and modify the original from
      //   which all other poses make deep copies before being modified and
      //   passed back
      if (keep_existing_networks_ || only_native_) {
        // add_reslabels_to_pose( nonconst_get_orig_pose(), *native_network );
        // place_rots_on_pose( out_pose, *native_network,
        // native_network->is_native, bridging_waters_,
        // native_network->waterrots.empty() );
        place_rots_on_pose(out_pose, *native_network,
                           native_network->is_native);
        if (native_network->hbond_vec.empty()) {
          (*scorefxn_)(
              out_pose);  // NEED TO SCORE BEFORE get_hbond_atom_pairs()
          get_hbond_atom_pairs(*native_network, out_pose);
          find_unsats(out_pose, *native_network);
          native_network->total_hbonds = native_network->hbond_vec.size();
        }
        if (!native_network->scored) {
          score_network_on_pose(out_pose, *native_network);
        }
        core::scoring::constraints::ConstraintSetOP cst_op(
            out_pose.constraint_set()
                ->clone());  // get existing csts and add to them
        // set_constraints( out_pose, *cst_op, native_network, write_cst_files_
        // && !bridging_waters_ );
        set_constraints(out_pose, *cst_op, native_network, write_cst_files_);
        out_pose.constraint_set(cst_op);  // add constraints to the pose
      }
      std::string outstring(
          (pdb_numbering())
              ? print_network_w_pdb_numbering(out_pose, *native_network, true)
              : print_network(*native_network));
      outstring =
          outstring +
          this->print_additional_info_for_net(*native_network, out_pose);
      native_network->outstring = outstring;
      if (TR.visible()) TR << "\t" << outstring << std::endl;
      comment_str = comment_str + "\n" + native_network->outstring;

      if (write_network_pdbs_) write_network_pdb(native_network);

      count++;
    }
    if (keep_existing_networks_ || only_native_) {
      core::pose::add_comment(
          out_pose, "HBNet native networks: ", "\n" + comment_str + "\n");
    }
    if (TR.visible()) {
      TR << std::endl;
    }
  }
}  // output_networks

core::pack::task::PackerTaskOP HBNet::create_ptask(
    core::pose::Pose const &pose, bool initialize_from_commandline /*=false*/) {
  using namespace core::pack::task;
  TR << " Creating packer task based on specified task operations..."
     << std::endl;
  runtime_assert(task_factory_ != nullptr);
  if (initialize_from_commandline) {
    task_factory_->push_back(core::pack::task::operation::TaskOperationOP(
        new operation::InitializeFromCommandline));
  }
  PackerTaskOP task = task_factory_->create_task_and_apply_taskoperations(pose);
  return task;
}

bool HBNet::task_is_valid(Pose const &pose) const {
  if (task_->total_residue() != pose.total_residue()) return false;
  for (Size i(1); i <= pose.total_residue(); ++i) {
    chemical::ResidueTypeCOP r = pose.residue_type_ptr(i);
    if (!task_->residue_task(i).is_original_type(r)) return false;
  }
  return true;
}

void HBNet::set_symmetry(Pose &pose) {
  // add_background_energies_ = 0; //we do not want to store background and
  // intrares energies at 1-body for symmetry
  // because 2-body energies of symm clones with themselves stored here

  init_scorefxn_ =
      core::scoring::symmetry::symmetrize_scorefunction(*init_scorefxn_);
  scorefxn_ = core::scoring::symmetry::symmetrize_scorefunction(*scorefxn_);

  if (!core::pose::symmetry::is_symmetric(pose)) {
    return;
  }
  if (TR.visible()) TR << "         POSE IS SYMMETRIC" << std::endl;
  symmetric_ = true;
  core::conformation::symmetry::SymmetricConformation const &SymmConf(
      dynamic_cast<core::conformation::symmetry::SymmetricConformation const &>(
          pose.conformation()));
  symm_info_ = SymmConf.Symmetry_Info();
  // max_network_size_ = max_network_size_*symm_info_->num_bb_clones();
  Size num_components(symm_info_->get_num_components());
  multi_component_ = num_components >= 2;
  if (TR.visible()) {
    TR << "         SETTING UP SYMMETRY_INFO: # independent residues = "
       << symm_info_->num_independent_residues()
       << "; # bb clones = " << symm_info_->num_bb_clones()
       << "; # total symm interfaces = " << symm_info_->num_interfaces()
       << std::endl;
  }
  if (multi_component_) {
    if (TR.visible())
      TR << "         DETECTED MULTICOMPONENT SYMMETRY: # sym_dof interfaces = "
         << core::pose::symmetry::sym_dof_names(pose).size() << std::endl;
  }
  for (Size ic = 1; ic <= pose.conformation().num_chains(); ++ic) {
    Size ic_begin = pose.conformation().chain_begin(ic);
    Size ic_end = pose.conformation().chain_end(ic);
    char chain = pose.chain(ic_begin);
    chain_bounds_[chain].first = ic_begin;
    chain_bounds_[chain].second = ic_end;
  }
}

Size HBNet::get_ind_res(Pose const &pose, Size const res_i) {
  Size resi_ind(res_i);
  if (symmetric_ && res_i > symm_info_->num_independent_residues()) {
    char resi_chain = pose.chain(res_i);
    if (multi_component_) {
      // symm_info_->component_lower_bound()
      std::map<char, std::pair<Size, Size>> component_bounds =
          symm_info_->get_component_bounds();
      char resi_comp = symm_info_->get_component_of_residue(res_i);
      resi_ind = res_i - chain_bounds_[resi_chain].first +
                 component_bounds[resi_comp].first;
    } else {
      resi_ind = res_i - chain_bounds_[resi_chain].first + 1;
    }
  }
  return resi_ind;
}

//// TODO Need to finish
// void
// HBNet::write_pymol_files( std::string name )
//{
//    utility::io::ozstream pml_output_stream;
//    std::string pml_fname( name );
//    std::string pdb_fname( pml_fname );
//    std::string str_cst = ".cst";
//    pdb_fname.replace(pdb_fname.find(str_cst),str_cst.length(),".pdb");
//    pml_output_stream.open(pml_fname, std::ios_base::out);
//    pml_output_stream << "# pml file for " << pml_fname;
//    //pml_output_stream << p->outstring << std::endl;
//    pml_output_stream << "load " << pdb_fname << std::endl;
//    pml_output_stream << "unset ignore_case," << std::endl;
//    pml_output_stream << "hide everything, all" << std::endl;
//    pml_output_stream << "bg_color white" << std::endl;
//    pml_output_stream << "space cmyk" << std::endl;
//    pml_output_stream << "color grey90, all" << std::endl;
//    pml_output_stream << "show cartoon, all" << std::endl;
//
//    if ( (this->ligand()) ) {
//        pml_output_stream << "select ligand, resi " << (this->ligand()) <<
//        std::endl;
//        pml_output_stream << "util.cbam ligand" << std::endl;
//    } else {
//        pml_output_stream << "select start_residues, resi ";
//        for ( std::set< core::Size >::iterator vit = start_res_vec_.begin();
//        vit != start_res_vec_.end(); ++vit ) {
//            pml_output_stream << *vit << "+";
//        }
//        pml_output_stream << std::endl;
//        //pml_output_stream << "color red, start_residues" << std::endl;
//    }
//    pml_output_stream << "select DesignResidues, resi ";
//    pml_output_stream << std::endl;
//    pml_output_stream << std::endl;
//    //std::string network_residues("select hbond_network, ");
//    //pml_output_stream << network_residues << std::endl;
//    pml_output_stream << "show sticks, hbond_network" << std::endl;
//    pml_output_stream << "util.cbao hbond_network" << std::endl; //color by
//    atom element, with carbon bright orange
//    pml_output_stream << "#drawing dashed lines for each h-bond in network:"
//    << std::endl;
//    //pml_output_stream << pymol_str_stream.str();
//    pml_output_stream.close();
//}

void HBNet::write_network_pdb(
    HBondNetStructOP p)  // better to pass object pointer or reference?
{
  if (p->network_pdb_written) return;

  runtime_assert(p->id != 0);
  Pose ala_copy = *ala_pose_;
  // place_waters_on_pose( ala_copy, p );
  // place_rots_on_pose( ala_copy, *p, p->is_native, bridging_waters_,
  // p->waterrots.empty() );
  place_rots_on_pose(ala_copy, *p, p->is_native);
  std::string current_out_tag =
      protocols::jd2::JobDistributor::get_instance()->current_output_name();
  // from jd2 output_tag; 1st pose will be _0001_1.pdb, 0001.pdb if only 1 pose;
  // additional poses with MPM are _0001_2.pdb, _0001_3.pdb...

  std::string net_prefix("");
  if (p->is_native)
    net_prefix = "native";
  else if (p->is_extended)
    net_prefix = "extended";
  std::string ext(".pdb");
  std::string pdb_tag(get_file_name(p->id, net_prefix, ext));
  utility::file::FileName pdb_out(pdb_tag);
  std::ostringstream oss;
  oss << basic::options::option[basic::options::OptionKeys::out::prefix]()
      << pdb_out.base();
  pdb_out.base(oss.str());
  if (basic::options::option[basic::options::OptionKeys::out::path::all]
          .user()) {
    pdb_out.path(
        basic::options::option[basic::options::OptionKeys::out::path::all]()
            .path());
    pdb_out.vol(
        basic::options::option[basic::options::OptionKeys::out::path::all]()
            .vol());
  } else {
    pdb_out.path("");
    pdb_out.vol("");
  }
  pdb_tag = pdb_out.name();

  core::pose::add_comment(ala_copy, "HBNet details: ",
                          "\n" + print_headers() +
                              this->print_additional_headers() + "\n" +
                              p->outstring);

  // to print csts to pdb info
  core::scoring::constraints::ConstraintSetOP cst_op(
      ala_copy.constraint_set()->clone());  // get existing csts and add to them
  set_constraints(ala_copy, *cst_op, p, false);

  ala_copy.dump_pdb(pdb_tag);

  p->network_pdb_written = true;
}

void HBNet::jack_store_network(HBondNetStructOP p, bool native) {
  runtime_assert(p->id != 0);
  if (!native) {
    assert(rotamer_sets_ != nullptr);
  }

  const core::Size network_index = p->id;

  if (network_residues_.size() < network_index) {
    network_residues_.resize(network_index);
  }

  // utility::vector1< core::conformation::ResidueOP > network_residues;
  const core::Size num_residues = p->residues.size();
  network_residues_[network_index].resize(num_residues);

  for (core::Size ii = 1; ii <= num_residues; ++ii) {
    HBondResStructCOP res = p->residues[ii];

    if (native) {
      network_residues_[network_index][ii] =
          orig_pose_->residue(res->resnum).clone();
    } else {
      const core::Size moltenres =
          rotamer_sets_->resid_2_moltenres((platform::uint)(res->resnum));
      network_residues_[network_index][ii] =
          rotamer_sets_->rotamer_for_moltenres(moltenres, res->rot_index)
              ->clone();
    }
    network_residues_[network_index][ii]->seqpos(res->resnum);
  }
}

// assumes the network rotamers are already placed onto the pose
void HBNet::set_constraints(
    Pose &pose, core::scoring::constraints::ConstraintSet &constraints,
    HBondNetStructOP p, bool write_cst_file /* false */) {
  runtime_assert(!(p->hbond_vec.empty()));

  utility::io::ozstream cst_output_stream;
  std::string cst_fname("");
  std::string comment_str(print_headers() + this->print_additional_headers());
  comment_str = "#" + comment_str + "\n" + "#" + p->outstring + "\n";

  if (write_cst_file && !(p->cst_file_written)) {
    std::string prefix("");
    if (p->is_native)
      prefix = "native_";
    else if (p->is_extended)
      prefix = "extended_";
    cst_fname = get_file_name(p->id, prefix, ".cst");
    core::pose::add_comment(pose,
                            "HBNet cst filename for " + prefix + "network_" +
                                utility::to_string(p->id),
                            cst_fname);
    cst_output_stream.open(cst_fname, std::ios_base::out);
    cst_output_stream << "# " << cst_fname << std::endl;
    cst_output_stream << "# " << p->outstring << std::endl;
  }
  for (utility::vector1<HBondCOP>::const_iterator i = p->hbond_vec.begin();
       i != p->hbond_vec.end(); ++i) {
    Size drsd((*i)->don_res());
    Size arsd((*i)->acc_res());
    Size don_hatm((*i)->don_hatm());
    Size datm(pose.residue(drsd).atom_base((int)(don_hatm)));
    Size aatm((*i)->acc_atm());
    // we only want to turn on constraints for h-bonds involving network
    // residues and atoms that will not change/be designed away downstream!
    Real xtol(0.20);
    numeric::xyzVector<core::Real> don_coordinates =
        pose.residue(drsd).atom(datm).xyz();
    numeric::xyzVector<core::Real> acc_coordinates =
        pose.residue(arsd).atom(aatm).xyz();
    Real don_acc_dist = don_coordinates.distance(acc_coordinates);
    if (don_acc_dist <= 3.0 && don_acc_dist >= 2.6) {
      don_acc_dist = 2.8;  // idealized value
    } else if (don_acc_dist > 3.0) {
      don_acc_dist = (don_acc_dist + 2.6) / 2.0;
      xtol = don_acc_dist - 2.6;
    }
    id::AtomID acc_id(aatm, arsd);
    id::AtomID don_id(datm, drsd);
    Real lb(don_acc_dist - xtol);
    Real ub(don_acc_dist + xtol);
    Real sd(0.2);
    core::scoring::func::FuncFactory func_fact;
    core::scoring::func::FuncOP const dist_func = func_fact.new_func("BOUNDED");
    std::ostringstream bounded_vals;
    std::string tag("hbond");
    bounded_vals << utility::to_string(lb) << " " << utility::to_string(ub)
                 << " " << utility::to_string(sd) << " " << tag;
    std::istringstream in(bounded_vals.str());
    dist_func->read_data(in);
    // TODO need more control of how strongly this is penalized
    // NEED option for angle constraint
    core::scoring::constraints::AtomPairConstraintCOP hbond_cst(
        new core::scoring::constraints::AtomPairConstraint(acc_id, don_id,
                                                           dist_func));
    constraints.add_constraint(hbond_cst);
    //(*i)->show(in_pose, 1, TR);
    Real angle_AHD((*i)->get_AHDangle(pose));
    Real angle_BAH((*i)->get_BAHangle(pose));
    Real dist_AH((*i)->get_HAdist(pose));
    Real hb_energy((*i)->energy());  // unweighted h-bond energy
    std::ostringstream cst_str_stream;
    cst_str_stream << "# angle_AHD = " << angle_AHD
                   << ", angle_BAH = " << angle_BAH << ", dist_AH = " << dist_AH
                   << angle_BAH << ", unweighted hb_energy = " << hb_energy
                   << "\n";
    cst_str_stream << "AtomPair " << pose.residue(arsd).atom_name((int)(aatm))
                   << " " << arsd << " "
                   << pose.residue(drsd).atom_name((int)(datm)) << " " << drsd
                   << " BOUNDED " << lb << " " << ub << " " << sd << " " << tag
                   << std::endl;
    comment_str = comment_str + cst_str_stream.str();
    if (write_cst_file && !(p->cst_file_written)) {
      cst_output_stream << cst_str_stream.str();
    }
  }
  core::pose::add_comment(
      pose, "HBNet constraints for network_" + utility::to_string(p->id) + ": ",
      "\n" + comment_str + "\n");
  if (write_cst_file && !(p->cst_file_written)) {
    cst_output_stream.close();
    p->cst_file_written = true;
  }
}

std::string HBNet::get_file_name(Size id, std::string net_prefix,
                                 std::string extension) {
  std::string current_out_tag =
      protocols::jd2::JobDistributor::get_instance()->current_output_name();
  std::string fname =
      current_out_tag +
      (basic::options::option[basic::options::OptionKeys::out::suffix]()) +
      "_" + net_prefix + "network_" + utility::to_string(id) + extension;
  utility::file::FileName outfile(fname);
  std::ostringstream oss;
  oss << basic::options::option[basic::options::OptionKeys::out::prefix]()
      << outfile.base();
  outfile.base(oss.str());
  if (basic::options::option[basic::options::OptionKeys::out::path::all]
          .user()) {
    outfile.path(
        basic::options::option[basic::options::OptionKeys::out::path::all]()
            .path());
    outfile.vol(
        basic::options::option[basic::options::OptionKeys::out::path::all]()
            .vol());
  } else {
    outfile.path("");
    outfile.vol("");
  }
  return outfile.name();
}

bool HBNet::water_clashes(Pose const &pose, Vector const water_O) {
  Real const acc_cutoff(2.0);
  Real const acc_cutoff_sq(acc_cutoff * acc_cutoff);

  // NEED TO add check for symmetric waters that clash with themselves
  // hacky check to see if water O is on z-axis (in which case it will clash
  // with symmetric copies of itself)
  if (symmetric_ && ((water_O.x() + water_O.y()) * (water_O.x() + water_O.y()) <
                     acc_cutoff_sq))
    return true;
  for (Size res = 1; res <= pose.total_residue(); ++res) {
    if (water_O.distance_squared(pose.residue(res).nbr_atom_xyz()) >
        144.0) { /* 12^2 */  // could get by with 10?
      continue;
    }
    if (water_oxygen_clashes_with_residue(water_O, pose.residue(res))) {
      return true;
    }
  }
  return false;
}

bool HBNet::water_oxygen_clashes_with_residue(Vector const water_oxygen,
                                              Size const resnum,
                                              int const rot_state) {
  // Should be able to look up atom coords in rotamer_sets_;
  // rotamer_building_functions orient and transform coords onto pose backbone
  if (rotamer_sets_ != nullptr &&
      rotamer_sets_->resid_2_moltenres(resnum) !=
          0) {  // returns false if res1 not molten (not design/repack-able)
    return water_oxygen_clashes_with_residue(
        water_oxygen,
        *(rotamer_sets_
              ->rotamer_set_for_moltenresidue(
                  rotamer_sets_->resid_2_moltenres(resnum))
              ->rotamer(rot_state)));
  }
  return false;
}

bool HBNet::water_oxygen_clashes_with_residue(Vector const water_oxygen,
                                              Residue const &res) {
  // should hard code these constants
  Real const hpol_cutoff(1.0);
  Real const hpol_cutoff_sq(hpol_cutoff * hpol_cutoff);
  Real const acc_cutoff(2.0);
  Real const acc_cutoff_sq(acc_cutoff * acc_cutoff);
  Real const clash_cutoff(1.8);
  Real const clash_cutoff_sq(clash_cutoff * clash_cutoff);

  for (Size at2 = 1; at2 <= res.natoms(); ++at2) {
    // skip virtual atoms but not virtual waters!
    std::string resname(res.name());
    // we want to compare against all waters, even if VRT (but skip other VRT
    // atoms)
    if (res.atom_type(at2).lj_wdepth() == 0.0 &&
        !(resname == "PWAT_V" || resname == "BB_PWAT_V" || resname == "VWAT" ||
          resname == "VWAT_V"))
      continue;
    Real const dist2(water_oxygen.distance_squared(res.xyz(at2)));
    if (!(res.atom_is_backbone(at2)) && res.atom_is_polar_hydrogen(at2) &&
        dist2 < hpol_cutoff_sq) {
      return true;
      //        } else if ( ( res.atom_type(at2).is_acceptor() || res.is_water()
      //        ) && dist2 < acc_cutoff_sq ) {
    } else if ((res.atom_type(at2).is_acceptor() || res.name1() == 'w') &&
               dist2 < acc_cutoff_sq) {
      return true;
    } else if (dist2 < clash_cutoff_sq) {
      return true;
    }
  }
  return false;
}

pose::PoseOP HBNet::get_additional_output() {
  pose::PoseOP out_pose;
  if (output_vector_.size() == 0) {
    return NULL;
  } else if (output_vector_.size() == 1) {
    // Last iteration, do not create a new pose copy
    out_pose = (output_poly_ala_background_) ? ala_pose_ : orig_pose_;
  } else {
    out_pose = (output_poly_ala_background_)
                   ? pose::PoseOP(new Pose(*ala_pose_))
                   : pose::PoseOP(new Pose(*orig_pose_));
  }
  std::set<Size> net_ids(output_vector_.back());  // back returs ref to back
  output_vector_
      .pop_back();  // pop_back() deletes/removes last, does not return;

  std::string comment_str(print_headers() + this->print_additional_headers());
  if (net_ids.size() > 1 && TR.visible()) {
    TR << std::endl << "combining networks ";
  }
  for (const auto &net_id : net_ids) {
    if (net_ids.size() > 1 && TR.visible()) {
      TR << net_id << "  ";
    }
    // TR << "net_id = " << net_id << std::endl;
    runtime_assert(get_network_by_id(net_id) != nullptr);
    // place_rots_on_pose( *out_pose, *(get_network_by_id( net_id )),
    // (get_network_by_id( net_id ))->is_native, bridging_waters_,
    // (get_network_by_id( net_id ))->waterrots.empty() );
    place_rots_on_pose(*out_pose, *(get_network_by_id(net_id)),
                       (get_network_by_id(net_id))->is_native);
    // add_reslabels_to_pose( *out_pose, *(get_network_by_id( net_id )) ); //
    // needs to be done on the fly when adding waters
    comment_str = comment_str + "\n" + (get_network_by_id(net_id))->outstring;
    // have to do this each time if adding residues (waters)
    core::scoring::constraints::ConstraintSetOP cst_op(
        out_pose->constraint_set()
            ->clone());  // get existing csts and add to them
    //        if ( bridging_waters_ ){
    //            //NEED TO SCORE BEFORE get_hbond_atom_pairs()
    //            ( *scorefxn_ )( *out_pose ); // TODO NEED BETTER WAY THAN
    //            SCORING EACH TIME AND OVERRIDING hbond_vec
    //            get_hbond_atom_pairs( *(get_network_by_id( net_id )),
    //            *out_pose ); // if residues (waters) have been added, need to
    //            update
    //        }
    set_constraints(*out_pose, *cst_op, get_network_by_id(net_id),
                    write_cst_files_);
    out_pose->constraint_set(cst_op);  // add constraints to the pose
  }
  if (TR.visible()) TR << std::endl;
  core::pose::add_comment(
      *out_pose, "HBNet designed networks: ", "\n" + comment_str + "\n");

  //(*out_pose).update_residue_neighbors(); // now taken care of in
  // place_rots_on_pose()
  (*scorefxn_)(*out_pose);  // always want to return a scored pose; and
                            // set_constraints() clears energies!! so need to
                            // re-score before calling get_hbond_atom_pairs()
  return out_pose;
}
/*
void
HBNet::setup( Pose & pose )
{
        if ( task_factory_ == nullptr ) {
                task_factory_ = task::TaskFactoryOP( new task::TaskFactory );
        }
        if( ! please_do_not_change_the_task_ ){
                task_ = create_ptask( pose ); //set task (which residues are
designable/repackable
        }
        //dangerous, if empty default is to start at every designable/packable
position in the pose
        if ( start_res_vec_.empty() ) {
                utility::vector1<bool> is_repack = task_->repacking_residues();
                runtime_assert(is_repack.size() == pose.total_residue());
                // TODO NEED TO make sure task_ is fully processed and trimmed
to asymmetric unit for symmetric cases by this point
                for ( Size r = 1; r <= pose.total_residue(); ++r ) {
                        if ( pose.residue(r).is_protein() ) {
                                if ( task_->design_residue((int)r) ||
is_repack[r]==1 ) {
                                        start_res_vec_.insert(r);
                                }
                        }
                }
        }
}
*/
void HBNet::search_IG_for_networks(Pose &) {
  // traverse IG and enumerate all possible h-bond networks given parameters
  if (ig_ != nullptr) {
    if (TR.visible()) {
      TR << " ================================================================="
            "="
         << std::endl;
      TR << " ============     PERFORMING H-BOND NETWORK DESIGN     "
            "============"
         << std::endl;
      TR << " ================================================================="
            "="
         << std::endl;
      if (greedy_traverse_) {
        greedy_traverse_IG(hydrogen_bond_threshold_);
      } else {
        traverse_IG(hydrogen_bond_threshold_);
      }
    }
  }
}

void HBNet::prepare_output() {
  output_networks(true);  // will add single networks to output vector
}

////NEED TO MOVE THIS TO IT'S OWN CLASS! //TODO
// void
// HBNet::benchmark_with_native( core::pose::Pose & pose ) // I changed this to
// & pose, may need to copy Pose
//{
// if ( TR.visible() ) TR << "BENCHMARKING AGAINST NATIVE NETWORKS: " <<
// std::endl;
// std::sort( network_vector_.begin(), network_vector_.end(), compare_net_vec()
// ); //sort all networks
//
// //for benchmarking/analysis:
// Size count = 1;
// Size ten = 10;
// Size twenty = 20;
// Size fifty = 50;
// Size total_rots = 0;
// Size total_correct_seq = 0;
// Size total_correct_rot = 0;
// Size native_seq = 0;
// Size native_rot = 0;
// Real seq_in_top_ten = 0;
// Real seq_rec_top_ten = 0;
// Real seq_rec_top_twenty = 0;
// Real seq_rec_top_fifty = 0;
// Real rot_in_top_ten = 0;
// Real rot_rec_top_ten = 0;
// Real rot_rec_top_twenty = 0;
// Real rot_rec_top_fifty = 0;
// Real seq_in_top_fiftyperc = 0.0;
// Real seq_in_top_twentyfiveperc = 0.0;
// Real rot_in_top_fiftyperc = 0.0;
// Real rot_in_top_twentyfiveperc = 0.0;
// Real fifty_perc = 0.5*(network_vector_.size());
// Real twentyfive_perc = 0.25*(network_vector_.size());
// Size int_fifty_perc = static_cast<Size>(utility::round(fifty_perc));
// Size int_twentyfive_perc =
// static_cast<Size>(utility::round(twentyfive_perc));
//
// for ( auto & i : network_vector_ ) {
//  //bool is_native_seq(0);
//  //bool is_native_rot(0);
//  if ( benchmark_ ) {
//   Real  net_size = i->residues.size();
//   if ( (this->ligand()) ) { //not fair to count ligand
//    for ( utility::vector1< HBondResStructCOP >::const_iterator listit =
//    i->residues.begin(); listit != i->residues.end(); ++listit ) {
//     if ( (*listit)->resnum == (this->ligand()) ) {
//      net_size = net_size -1;
//     }
//    }
//   }
//   pose::PoseOP ala_copy = pose::PoseOP( new pose::Pose( *ala_pose_) ); //
//   make copy of polyalanine pose
//   place_rots_on_pose( *ala_copy, i->residues, i->is_native );
//   //place_waters_on_pose(*ala_copy, *i);
//   ala_copy->update_residue_neighbors();
//
//   Size num_native_seq(get_num_native_seq( *ala_copy, i->residues ) );
//   Size num_native_rot(get_num_native_rot( *ala_copy, i->residues,
//   SC_RMSD_CUTOFF, 1 ) );
//   Real native_seq_rec = num_native_seq;
//   native_seq_rec = native_seq_rec/net_size;
//   Real native_rot_rec = num_native_rot;
//   native_rot_rec = native_rot_rec/net_size;
//
//   //            if (native_seq_rec==1)
//   //                is_native_seq=1;
//   //            if (native_rot_rec==1)
//   //                is_native_rot=1;
//
//   std::string network( ( pdb_numbering() ) ? ( print_list_to_string(
//   get_orig_pose(), i->residues ) ) : (print_list_to_string( i->residues ) )
//   );
//   std::string outstring = "NETWORK "+utility::to_string(count)+": " + network
//   + " has "+utility::to_string(100*native_seq_rec)+"% native sequence, and
//   has "+utility::to_string(100*native_rot_rec)+"% native rotamers\n ";
//   if ( TR.visible() ) TR << outstring;
//
//   if ( native_seq_rec == 1 ) {
//    native_seq++;
//   }
//   if ( native_rot_rec == 1 ) {
//    native_rot++;
//   }
//   total_rots = total_rots + net_size;
//   total_correct_seq = total_correct_seq + num_native_seq;
//   total_correct_rot = total_correct_rot + num_native_rot;
//   if ( count == ten ) {
//    seq_in_top_ten = 100.0*(native_seq/10.0);
//    rot_in_top_ten = 100.0*(native_rot/10.0);
//    seq_rec_top_ten = total_correct_seq;
//    seq_rec_top_ten = 100.0*(seq_rec_top_ten/total_rots);
//    rot_rec_top_ten = total_correct_rot;
//    rot_rec_top_ten = 100.0*(rot_rec_top_ten/total_rots);
//   } else if ( count == twenty ) {
//    seq_rec_top_twenty = total_correct_seq;
//    seq_rec_top_twenty = 100.0*(seq_rec_top_twenty/total_rots);
//    rot_rec_top_twenty = total_correct_rot;
//    rot_rec_top_twenty = 100.0*(rot_rec_top_twenty/total_rots);
//   } else if ( count == fifty ) {
//    seq_rec_top_fifty = total_correct_seq;
//    seq_rec_top_fifty = 100.0*(seq_rec_top_fifty/total_rots);
//    rot_rec_top_fifty = total_correct_rot;
//    rot_rec_top_fifty = 100.0*(rot_rec_top_fifty/total_rots);
//   }
//   if ( count == int_fifty_perc ) {
//    seq_in_top_fiftyperc =
//    100.0*(native_seq/(static_cast<Real>(int_fifty_perc)));
//    rot_in_top_fiftyperc =
//    100.0*(native_rot/(static_cast<Real>(int_fifty_perc)));
//   } else if ( count == int_twentyfive_perc ) {
//    seq_in_top_twentyfiveperc =
//    100.0*(native_seq/(static_cast<Real>(int_twentyfive_perc)));
//    rot_in_top_twentyfiveperc =
//    100.0*(native_rot/(static_cast<Real>(int_twentyfive_perc)));
//   }
//  }
//  count++;
// }
// // For analysis, writes out to TR benchmark statistics and details:
//
// Real native_seq_perc = 100.0*((static_cast< core::Real
// >(native_seq))/network_vector_.size());
// Real native_rot_perc = 100.0*((static_cast< core::Real
// >(native_rot))/network_vector_.size());
// Real seq_rec_tot = total_correct_seq;
// seq_rec_tot = 100.0*(seq_rec_tot/total_rots);
// Real rot_rec_tot = total_correct_rot;
// rot_rec_tot = 100.0*(rot_rec_tot/total_rots);
// if ( network_vector_.size() < ten ) {
//  seq_rec_top_ten = seq_rec_tot;
//  rot_rec_top_ten = rot_rec_tot;
// }
// if ( network_vector_.size() < twenty ) {
//  seq_rec_top_twenty = seq_rec_tot;
//  rot_rec_top_twenty = rot_rec_tot;
// }
// if ( network_vector_.size() < fifty ) {
//  seq_rec_top_fifty = seq_rec_tot;
//  rot_rec_top_fifty = rot_rec_tot;
// }
// if ( TR.visible() ) {
//  TR << " ===============================================================" <<
//  std::endl;
//  TR << " Sequence recovery = " << seq_rec_top_ten << " % top 10, " <<
//  seq_rec_top_twenty << " % top 20, " << seq_rec_top_fifty << " % top 50, " <<
//  seq_rec_tot << " % TOTAL, " << std::endl;
//  TR << " Rotamer recovery = " << rot_rec_top_ten << " % top 10, " <<
//  rot_rec_top_twenty << " % top 20, " << rot_rec_top_fifty << " % top 50, " <<
//  rot_rec_tot << " % TOTAL, " << std::endl;
//  TR << " NUMBER OF NETWORKS THAT HAVE 100% NATIVE SEQUENCE = " << native_seq
//  << "; " << native_seq_perc << " % total;  ";
//  TR << seq_in_top_twentyfiveperc << " % of top 25%, " << seq_in_top_fiftyperc
//  << " % of top 50%, " << seq_in_top_ten << " % of top 10 networks, " <<
//  std::endl;
//  TR << " NUMBER OF NETWORKS THAT HAVE 100% NATIVE ROTAMERS = " << native_rot
//  << "; " << native_rot_perc << " % total;  ";
//  TR << rot_in_top_twentyfiveperc << " % of top 25%, " << rot_in_top_fiftyperc
//  << " % of top 50%, " << rot_in_top_ten << " % top 10, " << std::endl;
// }
// //Real native_perc(0.0);
// native_seq_perc = 0.0;
// store_subnetworks_ = false;
//
// std::vector< HBondNetStructOP > orig_net_vec = network_vector_;
// network_vector_.clear();
// native_networks_.clear();
// //bool orig_design = design_;
// //bool orig_only_perfect_pairs = only_perfect_pairs_;
// Size orig_max_unsat = max_unsat_;
// Size orig_min_network_size = min_network_size_;
// Size orig_max_network_size = max_network_size_;
// //Real orig_upper_limit = upper_score_limit_;
// bool orig_min = minimize_;
//
// //ensure pose is scored
// pose.update_residue_neighbors();
// (*init_scorefxn_)(pose);
//
// //find_native_networks( pose );
//
// max_unsat_ = orig_max_unsat;
// min_network_size_ = orig_min_network_size;
// max_network_size_ = orig_max_network_size;
// //upper_score_limit_ = orig_upper_limit;
// minimize_ = orig_min;
//
// native_networks_ = network_vector_;
// network_vector_.clear();
// network_vector_ = orig_net_vec;
// Size lig((this->ligand()) ? 1 : 0);
// Size num_native_nets(native_networks_.size());
// Size found_native(0);
// Size found_native_sequence_50(0);
// Size found_native_sequence_66(0);
// Size found_native_sequence_75(0);
// Size found_native_w_subnetwork_50(0);
// Size found_native_w_subnetwork_66(0);
// Size found_native_w_subnetwork_75(0);
// for ( auto & native_network : native_networks_ ) {
//  bool found_nat_sub_50 = false;
//  bool found_nat_seq_50 = false;
//  bool found_nat_sub_66 = false;
//  bool found_nat_seq_66 = false;
//  bool found_nat_sub_75 = false;
//  bool found_nat_seq_75 = false;
//  for ( auto & netvec_it : network_vector_ ) {
//   //if ( networks_identical_aa_sequence( **native_it, **netvec_it, true ) ){
//   if ( residues_identical( native_network->residues, netvec_it->residues ) )
//   {
//    found_nat_seq_50 = true;
//    found_nat_seq_66 = true;
//    found_nat_seq_75 = true;
//    Size num_native_rot = get_num_native_rot( pose, netvec_it->residues );
//    if ( num_native_rot == (netvec_it->residues.size()) ) {
//     found_native++;
//     found_nat_sub_50 = true;
//     found_nat_sub_66 = true;
//     found_nat_sub_75 = true;
//     break;
//    } else if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.5 )
//    {
//     found_nat_sub_50 = true;
//     if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.66 ) {
//      found_nat_sub_66 = true;
//     }
//     if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.75 ) {
//      found_nat_sub_75 = true;
//     }
//     break;
//    }
//   } else if ( is_sub_residues( native_network->residues, netvec_it->residues
//   ) && ( netvec_it->residues.size() - lig) >=
//   0.5*(native_network->residues.size() - lig ) ) {
//    found_nat_seq_50 = true;
//    if ( 1.0*( netvec_it->residues.size() - lig) >=
//    0.66*(native_network->residues.size() - lig ) ) {
//     found_nat_seq_66 = true;
//    }
//    if ( 1.0*( netvec_it->residues.size() - lig) >=
//    0.75*(native_network->residues.size() - lig ) ) {
//     found_nat_seq_75 = true;
//    }
//    Size num_native_rot = get_num_native_rot( pose, netvec_it->residues );
//    if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.5 ) {
//     found_nat_sub_50 = true;
//    }
//    if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.66 ) {
//     found_nat_sub_66 = true;
//    }
//    if ( 1.0*num_native_rot/(netvec_it->residues.size() - lig) >= 0.75 ) {
//     found_nat_sub_75 = true;
//    }
//    break;
//   }
//  }
//  if ( found_nat_sub_50 ) {
//   found_native_w_subnetwork_50++;
//  }
//  if ( found_nat_seq_50 ) {
//   found_native_sequence_50++;
//  }
//  if ( found_nat_sub_66 ) {
//   found_native_w_subnetwork_66++;
//  }
//  if ( found_nat_seq_66 ) {
//   found_native_sequence_66++;
//  }
//  if ( found_nat_sub_75 ) {
//   found_native_w_subnetwork_75++;
//  }
//  if ( found_nat_seq_75 ) {
//   found_native_sequence_75++;
//  }
// }
// //native_perc = 1.0*found_native/num_native_nets;
// //native_seq_perc = 1.0*found_native_sequence/num_native_nets;
// if ( TR.visible() ) {
//  TR << " ===============================================================" <<
//  std::endl;
//  TR << " WHAT PERCENTAGE OF THE NATIVE NETWORKS WERE RECAPITULATED?" <<
//  std::endl;
//  TR << " There are " << num_native_nets << " native networks in the pose." <<
//  std::endl;
//  TR << " Recapitulated " << found_native << " of the complete native
//  networks;" << std::endl;
//  TR << "   " << found_native_w_subnetwork_50 << " networks that are a subset
//  with >= 50% rotamers identical to the native network" << std::endl;
//  TR << "   " << found_native_w_subnetwork_66 << " networks that are a subset
//  with >= 66.6% rotamers identical to the native network" << std::endl;
//  TR << "   " << found_native_w_subnetwork_75 << " networks that are a subset
//  with >= 75% rotamers identical to the native network" << std::endl;
//  TR << "   " << found_native_sequence_50 << " networks that are a subset with
//  >= 50% sequence identical to the native network" << std::endl;
//  TR << "   " << found_native_sequence_66 << " networks that are a subset with
//  >= 66.6% sequence identical to the native network" << std::endl;
//  TR << "   " << found_native_sequence_75 << " networks that are a subset with
//  >= 75% sequence identical to the native network" << std::endl;
// }
//}

Size HBNet::num_core_res(hbond_net_struct const &network) {
  Size num_core(0);
  for (const auto &residue : network.residues) {
    if (res_is_core(residue->resnum)) ++num_core;
  }
  return num_core;
}

Size HBNet::num_boundary_res(hbond_net_struct const &network) {
  Size num_boundary(0);
  for (const auto &residue : network.residues) {
    if (res_is_boundary(residue->resnum)) ++num_boundary;
  }
  return num_boundary;
}

void HBNet::select_best_networks() {
  TR << "JACK: network_vector_.size() == " << network_vector_.size()
     << " at beginning of select_best_networks()" << std::endl;
  const bool skip_many_filters =
      (basic::options::option
           [basic::options::OptionKeys::jackmag::bypass_secondary_hbnet_filters]
               .user() &&
       basic::options::option[basic::options::OptionKeys::jackmag::
                                  bypass_secondary_hbnet_filters]());

  auto i = network_vector_.begin();
  if (!skip_many_filters) {
    for (; i != network_vector_.end();) {
      if ((*i)->residues.size() < 3 && min_network_size_ > 2) {
        // if (verbose_) TR << "i = " << print_list_to_string( (*i)->residues )
        // << "ERASING! Size = " << (*i)->residues.size() << std::endl;
        i = network_vector_.erase(i);
      } else if (min_core_res_ && num_core_res(**i) < min_core_res_) {
        i = network_vector_.erase(i);
      } else if (min_boundary_res_ &&
                 num_boundary_res(**i) < min_boundary_res_) {
        i = network_vector_.erase(i);
      } else {
        ++i;
      }
    }
  }
  // TODO NEED TO PRECHECK

  // iterator over all the networks, erase those that do not meet specified
  // criteria
  // std::vector< HBondNetStructOP >::iterator i = network_vector_.begin();
  i = network_vector_.begin();
  for (; i != network_vector_.end();) {
    // ADD FUNCTION HERE TO CHECK IF NETWORK CONTAINS SPECIFICIED AA TYPES

    // place network residues on background pose for scoring and evaluation
    Pose ala_copy = *ala_pose_;  // deep copy
    // place_rots_on_pose( ala_copy, **i, (*i)->is_native, bridging_waters_,
    // pack_waters_before_unsat_check_ );
    place_rots_on_pose(ala_copy, **i, (*i)->is_native);
    // ala_copy.update_residue_neighbors();

    // very quickly identify networks that have heavy atom donors or acceptors
    // unsatisfied and get rid of them!
    if (no_heavy_unsats_allowed_ &&
        quick_and_dirty_network_has_heavy_atom_unsat(ala_copy, **i)) {
      if (verbose_ && TR.visible())
        TR << "i = " << print_list_to_string(**i)
           << "ERASING! did not pas quick_and_dirty_unsat_check" << std::endl;
      i = network_vector_.erase(i);
    } else {
      if (minimize_) {
        minimize_network(ala_copy, **i, true /* residues_already_placed */);
      }
      // Before get_hbond_atom_pairs(), NEED TO BE SURE POSE IS SCORED and
      // neighbors and energies are updated: residue_neighbors_updated
      if (!ala_copy.energies().energies_updated()) {
        (*scorefxn_)(ala_copy);
      }
      get_hbond_atom_pairs(**i, ala_copy);  // default sets bb_exclusion to
                                            // false

      // if ( (*i)->hbond_vec.empty() || ( bridging_waters_ &&
      // (*i)->has_bridging_wat && (*i)->bridging_water_oxygens.empty() ) ) {
      if ((*i)->hbond_vec.empty()) {
        if (verbose_)
          TR << "i = " << print_list_to_string(**i) << "ERASING! h-bonds empty"
             << std::endl;
        i = network_vector_.erase(i);
      } else {
        // find all buried unsatisfied polars in the network
        //  num unsats also applies only to ASU
        find_unsats(
            ala_copy,
            **i);  // anything < 0.0 considered h-bond during satisfaction check

        // move some of the to meets_criteria NEED TO FIX
        if (!skip_many_filters &&
            ((no_heavy_unsats_allowed_ && (*i)->num_heavy_unsat > 0) ||
             (*i)->num_unsat > max_unsat_)) {
          if (verbose_ && TR.visible()) {
            TR << "i = " << print_list_to_string(**i)
               << "ERASING! num unsat = " << (*i)->num_unsat
               << "; num_heavy_unsat = " << (*i)->num_heavy_unsat << std::endl;
          }
          i = network_vector_.erase(i);
        } else {
          //                    if ( bridging_waters_ && (*i)->has_bridging_wat
          //                    && !((*i)->bridging_water_oxygens.empty()) &&
          //                    !pack_waters_before_unsat_check_ ){ // pack
          //                    waters now then update networks
          //                        optimize_waters( ala_copy, **i ); // will
          //                        score the pose, so get_hbond_atom_pairs() is
          //                        good
          //                        get_hbond_atom_pairs( **i, ala_copy );
          //                        //default sets bb_exclusion to false
          //                        if ( (*i)->hbond_vec.empty() ){
          //                            i = network_vector_.erase( i );
          //                            continue;
          //                        }
          //                        find_unsats( ala_copy, **i );
          //                    }
          if ((*i)->num_unsat > max_unsat_ ||
              (*i)->connectivity < min_connectivity_ ||
              (this->ligand() && (*i)->lig_num_unsatisfied > max_lig_unsat_) ||
              !(this->network_meets_final_criteria(ala_copy, **i))) {
            if (verbose_ && TR.visible())
              TR << "i = " << print_list_to_string(**i)
                 << "ERASING! connectivity = " << (*i)->connectivity
                 << "unsat = " << (*i)->num_unsat << std::endl;
            i = network_vector_.erase(i);
          } else {
            score_network_on_pose(ala_copy, **i);
            if ((*i)->score > upper_score_limit_) {
              i = network_vector_.erase(i);
            } else {
              ++i;
            }
          }
        }
      }
    }
  }
  // score the networks
  if (verbose_ && TR.visible()) {
    TR << "NUMBER OF NETWORKS BEFORE SCORE = " << network_vector_.size()
       << std::endl;
    TR << " SCORING THE NETWORKS: " << std::endl;
  }
  // TODO NEED TO add charge-charge repulsion check
  // score_networks( minimize_ /* false */ ); // now scored above
  if (min_unique_networks_ > 1) {
    Size num_unique(1);
    for (auto i = network_vector_.begin(); i != network_vector_.end(); ++i) {
      for (auto j = i; ++j != network_vector_.end();) {
        if (networks_unique(**i, **j)) {
          num_unique++;
        }
        if (num_unique >= min_unique_networks_) {
          break;
        }
      }
      if (num_unique >= min_unique_networks_) {
        break;
      }
    }
    if (num_unique < min_unique_networks_) {
      TR << "CLEARING BECAUSE num_unique < min_unique_networks_" << std::endl;
      network_vector_.clear();
    }
  }
  TR << "Num networks at the end of select_best_networks: "
     << network_vector_.size() << std::endl;
}  // select_best_networks()

void HBNet::trim_rotamers(Pose &pose) { this->trim_additional_rotamers(pose); }

void HBNet::setup_packer_task_and_starting_residues(Pose const &pose) {
  if (task_factory_ == nullptr) {
    task_factory_ = task::TaskFactoryOP(new task::TaskFactory);
  }
  task_ =
      create_ptask(pose);  // set task (which residues are designable/repackable
  if (start_res_vec_.empty()) {
    utility::vector1<bool> is_repack = task_->repacking_residues();
    runtime_assert(is_repack.size() == pose.total_residue());
    // TODO NEED TO make sure task_ is fully processed and trimmed to asymmetric
    // unit for symmetric cases by this point
    for (Size r = 1; r <= pose.total_residue(); ++r) {
      if (pose.residue(r).is_protein()) {
        std::cout << __LINE__ << std::endl;
        if (task_->design_residue((int)r) || is_repack[r] == 1) {
          std::cout << __LINE__ << std::endl;
          start_res_vec_.insert(r);
        }
      }
    }
  }
}

void HBNet::get_native_networks(Pose const &pose) {
  utility::vector1<Size> residues_to_ala(0);
  for (Size r = 1; r <= pose.total_residue(); ++r) {
    if (pose.residue(r).is_ligand() || pose.residue_type(r).name() == "VRT") {
      continue;
    }
    residues_to_ala.push_back(r);
  }
  // set ala_pose_
  ala_pose_ = PoseOP(new Pose(pose));
  protocols::toolbox::pose_manipulation::construct_poly_XXX_pose(
      "ALA", *ala_pose_, residues_to_ala, 1, 1, 1);
  ala_pose_->update_residue_neighbors();
  (*scorefxn_)(*ala_pose_);  // score now so we don't have to later; several
                             // functions require ala_pose_ to be scored

  bool orig_only_native = only_native_;
  only_native_ = true;
  //    bool orig_bridging_waters = bridging_waters_;
  //    if ( ( find_native_ || keep_existing_networks_ ) && !only_native_ ){
  //        bridging_waters_ = false;
  //    }

  traverse_native(pose, hydrogen_bond_threshold_);
  branch_overlapping_networks();
  remove_replicate_networks();

  // if want to find natives to with criteria specified in HBNet mover defition
  // (independent of input HBNet PDBInfoLabels)
  if (find_native_) {
    select_best_networks();
  } else {
    for (auto &native_network : native_networks_) {
      Pose ala_copy = *ala_pose_;
      // place_rots_on_pose( ala_copy, *native_network,
      // native_network->is_native, bridging_waters_,
      // native_network->waterrots.empty() );
      place_rots_on_pose(ala_copy, *native_network, native_network->is_native);
      if (native_network->hbond_vec.empty()) {
        (*scorefxn_)(ala_copy);  // NEED TO SCORE BEFORE get_hbond_atom_pairs()
        get_hbond_atom_pairs(*native_network, ala_copy);
        find_unsats(ala_copy, *native_network);
        native_network->total_hbonds = native_network->hbond_vec.size();
      }
      if (!native_network->scored) {
        score_network_on_pose(ala_copy, *native_network);
      }
    }
  }
  only_native_ = orig_only_native;
  // bridging_waters_ = orig_bridging_waters;
  native_networks_ = network_vector_;
  network_vector_.clear();
}

void HBNet::setup(Pose const &pose) {
  this->setup_packer_task_and_starting_residues(pose);
  if (task_->rotamer_links_exist()) {
    if (TR.visible())
      TR << " WARNING!!!! ROTAMER LINKS DETECTED: HBNet can't currently handle "
            "rotamer links!"
         << std::endl;
    rotamer_links_ = task_->rotamer_links();
  }
  if (find_native_) {
    if (start_res_vec_.empty()) {
      if (TR.visible())
        TR << "No starting residues defined; traversing the entire static pose:"
           << std::endl;
      Size total((symmetric_) ? symm_info_->num_independent_residues()
                              : pose.total_residue());
      for (Size res = 1; res <= total; ++res) {
        start_res_vec_.insert(res);
      }
    }
    get_native_networks(pose);

    if (only_native_) {
      return;
    } else {
      // if we found native networks, add them to the "HBNet" info label
      // residues in case keep_existing or extend_exisitng options set to true
      core::select::residue_selector::ResiduePDBInfoHasLabelSelectorOP
          hbnet_info_label(new core::select::residue_selector::
                               ResiduePDBInfoHasLabelSelector("HBNet"));
      input_hbnet_info_residues_ = hbnet_info_label->apply(pose);
    }
  }
  if (only_extend_existing_) {
    start_res_vec_.clear();
    for (Size i = 1; i <= input_hbnet_info_residues_.size(); ++i) {
      if (input_hbnet_info_residues_[i]) {
        start_res_vec_.insert(i);
      }
    }
  }

  core::pack::task::operation::RestrictAbsentCanonicalAASOP set_hbnet_des_res =
      core::pack::task::operation::RestrictAbsentCanonicalAASOP(
          new core::pack::task::operation::RestrictAbsentCanonicalAAS());
  set_hbnet_des_res->keep_aas(des_residues_);
  set_hbnet_des_res->include_residue(0);
  set_hbnet_des_res->apply(pose, *task_);

  // HBNet requires PDInteractionGraph (or DensePDInteractionGraph).  Need to
  // check to make sure NOT using lin_mem_ig
  if ((*task_).linmem_ig() && TR.visible()) {
    TR << " ERROR: EXITING HBNet: You cannot use linmem_ig with HBNet; Please "
          "remove linmem_ig from HBNet task_operations."
       << std::endl;
    TR << " You can define a SetIGType task_operations in your XML and add it "
          "to the task_operations of other movers in your protocol, just not "
          "HBNet."
       << std::endl;
  }
  runtime_assert(!((*task_).linmem_ig()));
  runtime_assert(task_is_valid(pose));

  if (show_task_ && TR.visible()) {
    // task_->show();
    task_->show(TR);
  }
  // Setup IG ig_:
  if (symmetric_) {
    // need to truncate task to independent residues only
    utility::vector1<bool> allow_repacked(pose.total_residue(), false);
    for (Size res = 1; res <= pose.total_residue(); ++res) {
      if (pose.residue(res).aa() != core::chemical::aa_vrt &&
          symm_info_->fa_is_independent(res)) {
        allow_repacked.at(res) = true;
      }
    }
    task_->restrict_to_residues(allow_repacked);  //"vector boolean is based on
                                                  // residue position, disables
    // packing at false positions
    // does nothing to true
    // positions.  Cannot turn on
    // packing.
    rotamer_set::symmetry::SymmetricRotamerSetsOP sym_rotset_op(
        new rotamer_set::symmetry::SymmetricRotamerSets());
    rotamer_sets_ = sym_rotset_op;
  } else {
    rotamer_sets_ = RotamerSetsOP(new rotamer_set::RotamerSets());
  }
  // mechanism to keep/extend networks by HBNet InfoLabels in input pose
  if (keep_existing_networks_ || extend_existing_networks_) {
    // for every position with HBNet info label in input pdb

    std::set<core::Size> temp_start_res_vec_ =
        start_res_vec_;  // save starting residues
    if (keep_existing_networks_) {
      start_res_vec_.clear();
    }
    for (Size i = 1; i <= input_hbnet_info_residues_.size(); ++i) {
      if (input_hbnet_info_residues_[i]) {
        if (keep_existing_networks_) {
          start_res_vec_.insert(i);  // assumes all HBNet InfoLabel residues all
                                     // in ASU in input pose (always true unless
                                     // pose manually modified in between HBNet
                                     // runs)
          if (!extend_existing_networks_) {
            task_->nonconst_residue_task(i).prevent_repacking();
          }
        }
        if (extend_existing_networks_) {
          task_->nonconst_residue_task(i)
              .reset();  // in case previously set to prevent_repacking()
          // turning on optH restricts rotamers to only protons flips during
          // packing
          //   this simultaneously allows for
          task_->nonconst_residue_task(i).or_optimize_h(true);
          task_->nonconst_residue_task(i).or_include_current(true);
        }
      }
    }
    // if we've already found natives, skip this; else, detect the existing
    // native networks as defined by HBNet InfoLabels and store them for later
    if (keep_existing_networks_ && !find_native_) {
      get_native_networks(pose);
    }
    start_res_vec_ = temp_start_res_vec_;

    if (find_native_) {
      // for each rotamer in native_networks_: set them to not be
      // packable/designable in task_
      for (auto &native_network : native_networks_) {
        for (auto &residue : native_network->residues) {
          if (keep_existing_networks_ && !extend_existing_networks_) {
            task_->nonconst_residue_task(residue->resnum).prevent_repacking();
          } else if (extend_existing_networks_) {
            // do any of these flip it to packable if task_factory alraedy
            // designates it as prevent_repacking()?
            task_->nonconst_residue_task(residue->resnum)
                .restrict_to_repacking();  // ensure packable
            task_->nonconst_residue_task(residue->resnum)
                .or_include_current(true);
            task_->nonconst_residue_task(residue->resnum)
                .sample_proton_chi(true);
          }
        }
      }
    }
  }

  // keep_start_selector_rotamers_fixed_
  if (keep_start_selector_rotamers_fixed_) {
    for (const auto &start_res : start_res_vec_) {
      task_->nonconst_residue_task(start_res)
          .reset();  // in case previously set to prevent_repacking()
      task_->nonconst_residue_task(start_res).or_optimize_h(true);
      task_->nonconst_residue_task(start_res).or_include_current(true);
    }
  }
  utility::vector1<Size> residues_to_ala(0);
  utility::vector1<bool> is_repack = task_->repacking_residues();
  runtime_assert(is_repack.size() == pose.total_residue());
  for (Size r = 1; r <= pose.total_residue(); ++r) {
    if (pose.residue(r).is_protein()) {
      if (task_->design_residue((int)r) || is_repack[r] == 1) {
        residues_to_ala.push_back(r);
      }
    }
  }
  TR.Debug << " Storing background version of pose; design/repack residues set "
              "to Poly-ALA (keeping Pro,Gly,Cys-disulfide):"
           << std::endl;
  ala_pose_ = pose::PoseOP(new Pose(pose));
  protocols::toolbox::pose_manipulation::construct_poly_XXX_pose(
      "ALA", *ala_pose_, residues_to_ala, 1, 1, 1);
  ala_pose_->update_residue_neighbors();
  (*scorefxn_)(*ala_pose_);  // score now so we don't have to later; several
                             // functions require pose to be scored
}

void HBNet::run(Pose &pose) {
  if (only_native_) {
    return;  // no need to do design
  }
  // pose.update_residue_neighbors();
  (*init_scorefxn_)
      .setup_for_packing(pose, task_->repacking_residues(),
                         task_->designing_residues());
  packer_neighbor_graph_ = create_packer_graph(
      pose, *init_scorefxn_, task_);  // even in symmetric cases, packer graph
                                      // will have a node for every residue in
                                      // the pose

  if (run_offrot_) {
    TR << "OOPS I DIDN'T IMPLEMENT RUN_OFFROT_ YET!!!" << std::endl;
    return;
    // use_off_rotamer_pack_to_find_networks( pose );
  } else {
    pose.update_residue_neighbors();
    (*init_scorefxn_)
        .setup_for_packing(pose, task_->repacking_residues(),
                           task_->designing_residues());
    packer_neighbor_graph_ = create_packer_graph(
        pose, *init_scorefxn_, task_);  // even in symmetric cases, packer graph
                                        // will have a node for every residue in
                                        // the pose

    rotamer_sets_->set_task(task_);
    rotamer_sets_->initialize_pose_for_rotsets_creation(pose);
    rotamer_sets_->build_rotamers(pose, *init_scorefxn_,
                                  packer_neighbor_graph_);

    if (init_rotset_from_monte_carlo_offrot_) {
      core::pack::hbnet::MonteCarloHBNetOptionsOP mc_hbnet_options(
          new core::pack::hbnet::MonteCarloHBNetOptions);
      mc_hbnet_options->hydroxyls_must_donate(hydroxyls_must_donate_);
      mc_hbnet_options->tyr_hydroxyls_must_donate(tyr_hydroxyls_must_donate_);
      mc_hbnet_options->no_heavy_unsats_allowed(no_heavy_unsats_allowed_);
      mc_hbnet_options->max_unsat(max_unsat_);
      mc_hbnet_options->hbond_threshold(hydrogen_bond_threshold_);
      // mc_hbnet_options->allow_sc_bb_hbonds( true );
      mc_hbnet_options->allowed_aas(des_residues_);

      monte_carlo::OffRotamerHBNet offrot_hbnet;
      offrot_hbnet.score_function(init_scorefxn_);
      offrot_hbnet.set_options(mc_hbnet_options);
      offrot_hbnet.set_task(task_);
      offrot_hbnet.apply(*(orig_pose_->clone()));

      offrot_hbnet.add_rotamers_to_rotamer_sets(rotamer_sets_);
      for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres();
           ++mresid) {
        if (rotamer_sets_->rotamer_set_for_moltenresidue(mresid)
                ->num_rotamers() == 0) {
          TR << "OffRotamerHBNet did not find any rotamers for resid "
             << rotamer_sets_->moltenres_2_resid(mresid) << std::endl;
          rotamer_sets_->rotamer_set_for_moltenresidue(mresid)->build_rotamers(
              pose, *init_scorefxn_, *task_, packer_neighbor_graph_);
        }
      }
    } else if (init_rotset_from_monte_carlo_packrot_) {
      core::pack::hbnet::MonteCarloHBNetOptionsOP mc_hbnet_options(
          new core::pack::hbnet::MonteCarloHBNetOptions);
      mc_hbnet_options->hydroxyls_must_donate(hydroxyls_must_donate_);
      mc_hbnet_options->tyr_hydroxyls_must_donate(tyr_hydroxyls_must_donate_);
      mc_hbnet_options->no_heavy_unsats_allowed(no_heavy_unsats_allowed_);
      mc_hbnet_options->max_unsat(max_unsat_);
      mc_hbnet_options->hbond_threshold(hydrogen_bond_threshold_);
      // mc_hbnet_options->allow_sc_bb_hbonds( true );
      mc_hbnet_options->allowed_aas(des_residues_);

      monte_carlo::PackRotamerHBNet packrot_hbnet;
      packrot_hbnet.score_function(init_scorefxn_);
      packrot_hbnet.mc_options(mc_hbnet_options);
      packrot_hbnet.task(task_);
      packrot_hbnet.apply(*(orig_pose_->clone()));

      packrot_hbnet.add_rotamers_to_rotamer_sets(rotamer_sets_);
      for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres();
           ++mresid) {
        if (rotamer_sets_->rotamer_set_for_moltenresidue(mresid)
                ->num_rotamers() == 0) {
          TR << "PackRotamerHBNet did not find any rotamers for resid "
             << rotamer_sets_->moltenres_2_resid(mresid) << std::endl;
          rotamer_sets_->rotamer_set_for_moltenresidue(mresid)->build_rotamers(
              pose, *init_scorefxn_, *task_, packer_neighbor_graph_);
        }
      }
    }

    rotamer_sets_->prepare_sets_for_packing(pose, *init_scorefxn_);
    ig_ = InteractionGraphFactory::create_interaction_graph(
        *task_, *rotamer_sets_, pose, *init_scorefxn_, *packer_neighbor_graph_);
    ig_->initialize(*rotamer_sets_);

    PrecomputedPairEnergiesInteractionGraphOP pig(
        utility::pointer::dynamic_pointer_cast<
            PrecomputedPairEnergiesInteractionGraph>(ig_));
    runtime_assert(pig);

    if (symmetric_) {
      for (platform::uint ii = 1; ii <= rotamer_sets_->nmoltenres(); ++ii) {
        utility::vector1<core::PackerEnergy> one_body_energies(
            (rotamer_sets_->rotamer_set_for_moltenresidue(ii))->num_rotamers());
        hbnet_symm_one_body_energies(
            pose, rotamer_sets_->rotamer_set_for_moltenresidue(ii),
            *init_scorefxn_, *task_, packer_neighbor_graph_, one_body_energies);
        pig->add_to_nodes_one_body_energy(ii, one_body_energies);
      }
    }
    // precompute_two_body_energies() now vitural and dervied; no longer need to
    // cast to SymmetricRotamerSets first
    rotamer_sets_->precompute_two_body_energies(
        pose, *init_scorefxn_, packer_neighbor_graph_, pig, true);
    if (TR.visible()) {
      TR << " built " << rotamer_sets_->nrotamers() << " rotamers at "
         << rotamer_sets_->nmoltenres() << " positions." << std::endl;
    }
    if (TR.visible())
      TR << " IG: " << ig_->getTotalMemoryUsage() << " bytes" << std::endl;

    // Need to check that starting residues are in the design shell (IG); if
    // not, remove them
    utility::vector1<platform::uint> const molten_resvec(
        rotamer_sets_->moltenres_2_resid_vector());
    auto resvecit = start_res_vec_.begin();
    for (; resvecit != start_res_vec_.end();) {
      bool in_design_shell = false;
      for (core::Size mrvit : molten_resvec) {
        platform::uint resvecit_uint(*resvecit);
        if (resvecit_uint == mrvit) {
          in_design_shell = true;
          break;
        }
      }
      if (!(in_design_shell)) {
        if (TR.visible())
          TR << " WARNING: residue " << *resvecit
             << " is not in the IG; I'm removing it from start_resnums .  If "
                "you wish to consider this residue, change your design shell "
                "task operations"
             << std::endl;
        // does NOT return iterator referencing same location after removal like
        // std::vector does!
        start_res_vec_.erase(resvecit++);
      } else {
        ++resvecit;
      }
    }
    runtime_assert(!(start_res_vec_.empty()));
    this->search_IG_for_networks(pose);  // can be overriden for special cases
                                         // of handling and searching IG
  }

  if (TR.visible()) {
    TR << " INITIAL NUMBER OF H-BOND NETWORKS FOUND: " << network_vector_.size()
       << std::endl;
    TR << " BRANCHING NETWORKS TOGETHER TO FORM COMPLETE NETWORKS: "
       << std::endl;
  }
  branch_overlapping_networks();
  if (TR.visible())
    TR << " NUMBER OF H-BOND NETWORKS AFTER BRANCH: " << network_vector_.size()
       << std::endl;
  remove_replicate_networks(max_rep_);
  if (TR.visible())
    TR << "NUMBER OF NETWORKS AFTER REMOVE_REP = " << network_vector_.size()
       << std::endl;

  select_best_networks();
}  // run

void HBNet::apply(Pose &pose) {
#ifdef DEBUG_TASK_JACK
  Size num_being_designed = 0;
  Size num_being_packed_and_not_designed = 0;
  for (Size resid = 1; resid <= pose.size(); ++resid) {
    if (task_->being_designed(resid))
      ++num_being_designed;
    else if (task_->being_packed(resid))
      ++num_being_packed_and_not_designed;
  }
  TR << "Packing: " << num_being_packed_and_not_designed << std::endl;
  TR << "Designing: " << num_being_designed << std::endl;
#endif
  if (!(pose.pdb_info())) {
    pose.pdb_info(core::pose::PDBInfoOP(new core::pose::PDBInfo(pose, true)));
  }
  // Delete all HBNet comments in pose
  std::map<std::string, std::string> comments =
      core::pose::get_all_comments(pose);
  for (std::map<std::string, std::string>::const_iterator com =
           comments.begin();
       com != comments.end(); ++com) {
    if (com->first.substr(0, 5) != "HBNet") {
      core::pose::delete_comment(pose, com->first);
    }
  }
  // need to reset in case apply called multiple times from code with different
  // poses
  network_vector_.clear();
  rotamer_sets_.reset();  // resets OP to 0
  ig_.reset();            // resets OP to 0

  // shouldn't ever be NULL
  if (scorefxn_ == nullptr) {
    if (basic::options::option[basic::options::OptionKeys::corrections::beta]
            .value(true)) {
      scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
          "beta_cst");
    } else {
      scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(
          "talaris2013_cst");
    }
  }
  if (core::pose::symmetry::is_symmetric(pose)) {
    set_symmetry(pose);
  }

  if (!only_native_ && task_factory_ != nullptr) {
    if (!please_do_not_change_the_task_) {
      task::PackerTaskOP raw_task = create_ptask(pose);
      if (raw_task->rotamer_links_exist()) {
        if (TR.visible())
          TR << " WARNING!!!! ROTAMER LINKS DETECTED: HBNet can't currently "
                "handle rotamer links!"
             << std::endl;
        rotamer_links_ = raw_task->rotamer_links();
      }
    }
  }

  pose.update_residue_neighbors();
  (*init_scorefxn_)(pose);  // need this?

  TR << "init_scorefxn_: " << init_scorefxn_->get_name() << std::endl;
  TR << "scorefxn_: " << scorefxn_->get_name() << std::endl;

  orig_pose_ = PoseOP(new Pose(pose));  // pose::PoseOP pointing to original
                                        // pose
  // output_pose_ = pose.get_self_ptr(); //careful need to reset this pointer
  // before apply() returns otherwise will create problems...

  if (verbose_ && TR.visible()) {
    TR << " ==================================================================="
          "====="
       << std::endl;
    if (only_native_) {
      TR << " WILL ONLY CONSIDER NATIVE INPUT ROTAMERS IN NETWORK SEARCH"
         << std::endl;
      TR << " NO DESIGN OR DESIGN OPTIONS WILL BE CONSIDERED" << std::endl;
    } else {
      TR << " hydrogen bond threshold, hb_threshold = "
         << hydrogen_bond_threshold_ << std::endl;
      TR << " clash_threshold = " << clash_threshold_ << std::endl;
      TR << " Residue types allowed for h-bond network design = "
         << des_residues_ << std::endl;
    }
    TR << " store_subnetworks = " << store_subnetworks_ << std::endl;
    TR << " ==================================================================="
          "====="
       << std::endl;
  }

  // put all residues in selector in start_res_vec_ (we check later to ensure
  // that these positions are actually in IG)
  //    ResidueSubset is a bool vector of all positions in the pose; can use
  //    Vikram's ReferencePose if need to add/delete residues
  if (start_selector_) {
    core::select::residue_selector::ResidueSubset rs(
        start_selector_->apply(pose));
    start_res_vec_.clear();
    for (Size r = 1; r <= rs.size(); ++r) {
      if (rs[r]) {
        start_res_vec_.insert(r);
      }
    }
  }
  if (core_selector_) {
    core_residues_ = core_selector_->apply(pose);
  } else {
    core::select::residue_selector::LayerSelectorOP core_layer(
        new core::select::residue_selector::LayerSelector());
    core_layer->set_layers(true, false, false);
    // core_layer->set_use_sc_neighbors( true ); // now true by default
    core_layer->set_cutoffs(4.4 /* core */,
                            2.0 /* surface */);  // default core is 5.2, here we
                                                 // use 4.4 to be overly
                                                 // cautious with unsats
    runtime_assert(core_layer->use_sc_neighbors());
    core_residues_ = core_layer->apply(pose);
  }
  if (boundary_selector_) {
    boundary_residues_ = boundary_selector_->apply(pose);
  } else {
    core::select::residue_selector::LayerSelectorOP boundary_layer(
        new core::select::residue_selector::LayerSelector());
    boundary_layer->set_layers(false, true, false);
    // boundary_layer->set_use_sc_neighbors( true );
    boundary_layer->set_cutoffs(
        4.0 /* core */, 2.0 /* surface */);  // default core is 5.2, here we use
                                             // 4 to be overly cautious with
                                             // unsats
    boundary_residues_ = boundary_layer->apply(pose);
  }
  if (verbose_ && TR.visible()) {
    for (Size r = 1; r <= pose.total_residue(); ++r) {
      if (res_is_core(r)) {
        TR << "res " << r << " is core" << std::endl;
      }
    }
  }
  if (keep_existing_networks_ || extend_existing_networks_) {
    // Need to check: are PDBInfoLabels kept after pose made symmetric? what
    // about after BGS?
    core::select::residue_selector::ResiduePDBInfoHasLabelSelectorOP
        hbnet_info_label(
            new core::select::residue_selector::ResiduePDBInfoHasLabelSelector(
                "HBNet"));
    input_hbnet_info_residues_ = hbnet_info_label->apply(pose);
  }

  // setup method common to all HBNet
  setup(pose);
  // run method common to all HBNet
  run(pose);

  TR.flush();
  this->prepare_output();
  if (only_native_) {
    // in only_native_ case, HBNet reslabels and comments get added to cached
    // orig_pose_, but apply() is passed Pose & pose REFERENCE, so need update
    // "pose"
    pose = nonconst_get_orig_pose();
    if (native_networks_.size() == 0)
      set_last_move_status(protocols::moves::FAIL_RETRY);
    else
      set_last_move_status(protocols::moves::MS_SUCCESS);
    return;
  }
  TR.Debug << "After prepare_output, network_vector_.size() = "
           << network_vector_.size() << std::endl;
  if (output_vector_.size() == 0) {
    if (TR.visible())
      TR << "DID NOT FIND SOLUTIONS THAT MEET YOUR CRITERIA! EXITING"
         << std::endl;
    set_last_move_status(protocols::moves::FAIL_RETRY);
    ig_.reset();  // resets OP to 0
    // lkb_.reset();
    packer_neighbor_graph_.reset();
    return;
  } else {
    set_last_move_status(protocols::moves::MS_SUCCESS);
  }
  // if ( benchmark_ ) {
  //  benchmark_with_native( pose );
  // }

  // need to reverse because no pop_front for std::vector -- could use list but
  // already set up for sorting functions with vectors of HBondNetStructOP's
  std::reverse(output_vector_.begin(), output_vector_.end());

  std::set<Size> net_ids(output_vector_.back());
  output_vector_.pop_back();  // pop_back() removes without returning

  Pose ala_copy = *ala_pose_;                        // deep copy
  if (output_poly_ala_background_) pose = ala_copy;  // "pose" is & ref

  std::string comment_str(print_headers() + this->print_additional_headers());
  if (net_ids.size() > 1 && TR.visible()) {
    TR << std::endl << "combining networks ";
  }
  for (const auto &net_id : net_ids) {
    if (net_ids.size() > 1 && TR.visible()) {
      TR << net_id << "  ";
    }
    // TR << "net_id = " << *net_id << std::endl;
    runtime_assert(get_network_by_id(net_id) != nullptr);
    // place_rots_on_pose( pose, *(get_network_by_id( net_id )),
    // (get_network_by_id( net_id ))->is_native, bridging_waters_,
    // (get_network_by_id( net_id ))->waterrots.empty() );
    place_rots_on_pose(pose, *(get_network_by_id(net_id)),
                       (get_network_by_id(net_id))->is_native);
    // add_reslabels_to_pose( pose, *(get_network_by_id( net_id )) ); // needs
    // to be done on the fly when adding waters
    comment_str = comment_str + "\n" + (get_network_by_id(net_id))->outstring;
    core::scoring::constraints::ConstraintSetOP cst_op(
        pose.constraint_set()->clone());  // get existing csts and add to them

    //        if ( bridging_waters_ ){
    //            //NEED TO SCORE BEFORE get_hbond_atom_pairs()
    //            ( *scorefxn_ )( pose ); // TODO NEED BETTER WAY THAN SCORING
    //            EACH TIME AND OVERRIDING hbond_vec
    //            get_hbond_atom_pairs( *(get_network_by_id( net_id )), pose );
    //            // if residues (waters) have been added, need to update
    //        }
    set_constraints(pose, *cst_op, get_network_by_id(net_id), write_cst_files_);
    pose.constraint_set(cst_op);  // add constraints to the pose
  }
  if (TR.visible()) TR << std::endl;
  core::pose::add_comment(pose,
                          "HBNet Design details: ", "\n" + comment_str + "\n");

  // pose.update_residue_neighbors(); //now taken care of in
  // place_rots_on_pose()
  (*scorefxn_)(pose);

  ig_.reset();  // resets OP to 0 / nullptr
  packer_neighbor_graph_.reset();
  if (TR.visible()) TR.flush();
  if (TR.Warning.visible()) TR.Warning.flush();
  if (TR.Debug.visible()) TR.Debug.flush();
  // output_pose_.reset(); //If make PoseOP that directly points to original
  // Pose & pose of apply(), need to reset it!
}  // apply

std::string HBNet::get_name() const { return mover_name(); }

std::string HBNet::mover_name() { return "HBNet"; }

void HBNet::attributes_for_hbnet(utility::tag::AttributeList &attlist) {
  using namespace utility::tag;
  attlist +
      XMLSchemaAttribute::attribute_w_default(
          "hb_threshold", xsct_real,
          "2-body h-bond energy cutoff to define rotamer pairs that h-bond. "
          "I've found that -0.5 without ex1-ex2 is the best starting point. "
          "If using ex1-ex2, try -0.75. This parameter is the most important "
          "and requires some tuning; "
          "the tradeoff is that the more stringent (more negative), the faster "
          "it runs but you miss a lot of networks; "
          "too positive and it will run forever; using ex1-ex2 results in many "
          "redundant networks that end up being filtered out anyway.",
          "-0.5") +
      XMLSchemaAttribute::attribute_w_default("onebody_hb_threshold", xsct_real,
                                              "h-bond energy cutoff with their "
                                              "symmetric clones to define "
                                              "rotamer pairs that h-bond. ",
                                              "-0.4") +
      XMLSchemaAttribute::attribute_w_default(
          "charge_charge_rep_cutoff", xsct_real,
          "new feature, not fully tested, use at your own risk!;", "1.0") +
      XMLSchemaAttribute::attribute_w_default(
          "clash_threshold", xsct_real,
          "fa_rep energy above which two rotamers are considered clashing; "
          "used in all HBNet clash checks",
          "1.0") +
      XMLSchemaAttribute::attribute_w_default("upper_score_limit", xsct_real,
                                              "upper score limit for network "
                                              "against background pose to weed "
                                              "out anything terrible",
                                              "15.0") +
      XMLSchemaAttribute::attribute_w_default(
          "find_native_networks", xsct_rosetta_bool,
          "Will find and report native networks in input pose but will also do "
          "design; "
          "for keep_existing_networks=true or extend_existing_networks=true, "
          "in addition to starting from \"HBNet\" PDBInfoLabel tags, "
          "HBnet will find all native networks in your input pose that meet "
          "your criteria (specified by options).",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "output_poly_ala_background", xsct_rosetta_bool,
          "for returned poses, rather than place network onto input pose, "
          "place into poly-ALA background; "
          " where all positions that are designable or packable (according to "
          "user-defined task operations)"
          " are poly-alanine",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "find_only_native_networks", xsct_rosetta_bool,
          "HBnet will find only find native networks in your input pose that "
          "meet your criteria "
          "(specified by options) and return a single pose with csts for those "
          "networks. "
          "If this option is true, all options below are overridden and HBNet "
          "does not search for new networks",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "keep_existing_networks", xsct_rosetta_bool,
          "In addition to design, Keeps existing networks from the input pose "
          "for each pose returned by HBNet, "
          "and turn on csts for all; existing networks are identified by "
          "default by \"HBNet\" PDBInfoLabel "
          "tags in the input pose, but can also be detected anew by setting "
          "find_native_networks=1. "
          "If keep_existing_networks=true and extend_existing_networks=false, "
          "then HBNet internally turns off "
          "design at input network positions (prevent_repacking); new networks "
          "are searched for and designed at "
          "the other positions based on your task_ops. If "
          "keep_existing_networks=true and extend_existing_networks=true, "
          "then HBNet internally puts only the input rotamer of each network "
          "residue into the IG to try to extend; "
          "an extended network will replace its native network if it is better "
          "(otherwise native networks retained).",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "extend_existing_networks", xsct_rosetta_bool,
          "Detects existing networks and tries to extend them, and also will "
          "search for new networks at other positions "
          "based on your criteria. Existing networks identified by HBNet "
          "PDBInfoLabel tags in the input pose, "
          "but can also be detected anew by setting find_native_networks=1. "
          "For existing networks, HBNet internally "
          "puts only the input rotamer of each network residue into the IG to "
          "try to extend; "
          "an extended network will replace its native network if it is better "
          "(otherwise native networks retained).",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "only_extend_existing", xsct_rosetta_bool,
          "Will not look for new networks at other positions; will only try to "
          "extend and improve the existing networks.",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "minimize", xsct_rosetta_bool,
          "new feature, not fully tested, use at your own risk!; "
          "will minimize network on background pose before scoring/evaluating",
          "true") +
      XMLSchemaAttribute::attribute_w_default(
          "store_subnetworks", xsct_rosetta_bool,
          "store subnetworks as separate networks that can be passed back; not "
          "recommended",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "secondary_search", xsct_rosetta_bool,
          "if during IG traversal, a search trajectory terminates in rotamer "
          "than cannot make any h-bonds, "
          "search again from that rotamer using a lower hb_threshold (-0.25)",
          "false") +
      XMLSchemaAttribute::attribute_w_default("secondary_threshold", xsct_real,
                                              "for secondary search option, "
                                              "what is the second (lower) "
                                              "hydrogen_bond_threshold to use",
                                              "-0.25") +
      XMLSchemaAttribute::attribute_w_default(
          "write_network_pdbs", xsct_rosetta_bool,
          "writes out pdb files of only the network (in poly-Ala background "
          "where any designable/packable "
          "residue is Ala -- rest of pose untouched); "
          "this is simply for debugging and visualizing the network as "
          "detected by HBNet",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "jack_store_networks", xsct_rosetta_bool,
          "Save found networks in an internal data structure.", "false") +
      XMLSchemaAttribute::attribute_w_default(
          "write_cst_files", xsct_rosetta_bool,
          "writes out Rosetta .cst constraint files for each network", "true") +
      XMLSchemaAttribute::attribute_w_default(
          "max_replicates", xsct_non_negative_integer,
          "max replicate networks that are allowed (defined as networks with "
          "same residue positions and aa types), "
          "but different rotamers; deafult is that replicates are not allowed",
          "1") +
      XMLSchemaAttribute::attribute_w_default(
          "max_replicates_before_branch", xsct_non_negative_integer,
          "hack that improves runtime without significantly affecting results; "
          "use at your own risk! "
          "if set to be greater than 0, is the max replicate networks that are "
          "allowed before branching (defined as networks with same residue "
          "positions and aa types)",
          "0") +
      XMLSchemaAttribute::attribute_w_default(
          "max_replicates_before_unsat_check", xsct_non_negative_integer,
          "hack that improves runtime without significantly affecting results; "
          "use at your own risk! "
          "max replicate networks that are allowed before unsat check is "
          "performed (defined as networks with same residue positions and aa "
          "types)",
          "1") +
      XMLSchemaAttribute::attribute_w_default(
          "tyr_hydroxyls_must_donate", xsct_rosetta_bool,
          "new feature, not fully tested, use at your own risk!; "
          "require that TYR must donate to be considered satisfied",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "hydroxyls_must_donate", xsct_rosetta_bool,
          "new feature, not fully tested, use at your own risk!; "
          "require that hydroxyls must donate their H to be considered "
          "satisfied",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "use_pdb_numbering", xsct_rosetta_bool,
          "in tracer output, use pdb numbering and chains when printing out "
          "network info (instead of Rosetta sequential numbering)",
          "true") +
      XMLSchemaAttribute::attribute_w_default(
          "no_heavy_unsats_allowed", xsct_rosetta_bool,
          "see max_unsat for details", "1") +
      XMLSchemaAttribute::attribute_w_default(
          "show_task", xsct_rosetta_bool,
          "will show packer task by printing to TR, aka what's designable and "
          "packable and fixed based on your taskops",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "min_network_size", xsct_non_negative_integer,
          "minimum number of residues required for each network. Note: in "
          "symmetric cases, "
          "this refers to the number of residues in the ASU!",
          "3") +
      XMLSchemaAttribute::attribute_w_default(
          "max_network_size", xsct_non_negative_integer,
          "maximum number of residues required for each network", "15") +
      XMLSchemaAttribute::attribute_w_default(
          "min_unique_networks", xsct_non_negative_integer,
          "minimum number of networks with unique aa composition / rotamers",
          "1") +
      XMLSchemaAttribute::attribute_w_default(
          "min_core_res", xsct_non_negative_integer,
          "minimum core residues each network must have (as defined by core "
          "selector)",
          "0") +
      XMLSchemaAttribute::attribute_w_default(
          "min_boundary_res", xsct_non_negative_integer,
          "minimum boundary residues a network must contain", "0") +
      XMLSchemaAttribute::attribute_w_default(
          "max_unsat", xsct_non_negative_integer,
          "maximum number of buried unsatisfied polar atoms allowed in each "
          "network. "
          "Note that the way I treat buried unsats. is very different from all "
          "of the other "
          "Buried Unsatisfied calculators/filters in Rosetta. I have plans to "
          "move this code "
          "outside of HBNet and turn it into its own calculator/filter. Short "
          "version is that "
          "if there are heavy atom donors or acceptors that are buried and "
          "unsatisfied, "
          "those networks are thrown out, and I only count unsatisfied Hpols "
          "where the heavy "
          "atom donor is making at least one hydrogen bond. This behavior can "
          "be overridden "
          "to allow heavy atom unsats by setting "
          "no_heavy_unsats_allowed=false.",
          "5") +
      XMLSchemaAttribute::attribute_w_default(
          "min_connectivity", xsct_real,
          "minimum connectivity a network must have, defined as num_hbonds / "
          "num_polar_atoms_in_network (h_pol and acceptors)",
          "0.5") +
      XMLSchemaAttribute::attribute_w_default(
          "start_from_csts", xsct_rosetta_bool,
          "new feature, not fully tested, use at your own risk!; "
          "will detect csts in input pose and start network search from those "
          "cst residues",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "verbose", xsct_rosetta_bool,
          "print out all HBNet tracer statements; only useful for debugging",
          "false") +
      XMLSchemaAttribute::attribute_w_default(
          "design_residues", xs_string,
          "string of one-letter AA codes; which polar residues types do you "
          "want to include "
          "in the search; the default is all AA's that can potentially make "
          "h-bonds, "
          "further restricted by the task_operations you pass.",
          "STRKHYWNQDE") +
      XMLSchemaAttribute::attribute_w_default(
          "use_only_input_rot_for_start_res", xsct_rosetta_bool,
          "keep the input rotamer for starting residues fixed; only sample "
          "proton chis during network search",
          "false") +
      XMLSchemaAttribute("start_resnums", xsct_residue_number_cslist,
                         "comma delimited list of residue positions to start "
                         "network search from "
                         "(e.g. \"1,2,3,5\"); now is better to use "
                         "start_selector residue selector") +
      XMLSchemaAttribute("start_selector", xs_string,
                         "residue selector that tells HBNet which residues to "
                         "start from (will only search for networks that "
                         "include these resid") +
      XMLSchemaAttribute(
          "core_selector", xs_string,
          "residue selector that defines what HBNet considers \"core\"; "
          "used in buriedness determination for unsats; "
          "default is layer selector default using sidechain "
          "neighbors(core=5.2).") +
      XMLSchemaAttribute("boundary_selector", xs_string,
                         "residue selector to define what HBNet considers "
                         "boundary during comparisons/scoring of networks");

  rosetta_scripts::attributes_for_parse_task_operations(attlist);
  rosetta_scripts::attributes_for_parse_score_function(attlist);
}
void HBNet::provide_xml_schema(utility::tag::XMLSchemaDefinition &xsd) {
  using namespace utility::tag;
  AttributeList attlist;
  attributes_for_hbnet(attlist);

  protocols::moves::xsd_type_definition_w_attributes(
      xsd, mover_name(),
      "HBNet is a method to explicitly detect and design hydrogen bond "
      "networks within Rosetta. "
      "It functions as a mover within the RosettaScripts framework and will "
      "exhaustively search "
      "for all networks within the design space that you define with "
      "TaskOperations",
      attlist);
}

void HBNet::consider_network_as_branch_seed(HBondNetStructOP net,
                                            core::Size network_vector_index) {
  // TODO rank based on score?
  TR << "consider_network_as_branch_seed scored: " << net->scored << std::endl;
  branch_seed_network_vector_.push_back(network_vector_index);
}

std::string HBNetCreator::keyname() const { return HBNet::mover_name(); }

protocols::moves::MoverOP HBNetCreator::create_mover() const {
  return protocols::moves::MoverOP(new HBNet);
}

void HBNetCreator::provide_xml_schema(
    utility::tag::XMLSchemaDefinition &xsd) const {
  HBNet::provide_xml_schema(xsd);
}

// Assumes the XML file has been through filtering and sorting to
// extract only "true" (set to 1) variable numbers, in increasing
// numbers. With this we have rotamer variables first. This can be
// achieved by:
// cat test.sol|grep 'value="1"'|cut -d"x" -f2|cut -d'"' -f1| sort -n > sol
void load_ILP_conformation(const std::string &filename, core::Size n_mres,
                           utility::vector1<core::Size> &rotamer_assignment) {
  std::ifstream in;
  in.open(filename);

  core::Size nb_res = 0;
  for (std::string line; std::getline(in, line);) {
    std::stringstream iss(line);
    int rotid;
    while (iss >> rotid) {
      if (nb_res < n_mres) {
        nb_res++;
        rotamer_assignment.push_back(rotid);
      }
    }
  }
  in.close();
  return;
}

void load_SAT_conformation(const std::string &filename, core::Size n_mres,
                           utility::vector1<core::Size> &rotamer_assignment) {
  std::ifstream in;
  in.open(filename);

  for (std::string line; std::getline(in, line);) {
    if (line[0] == 'v') {
      line[0] = ' ';  // remove v
      std::stringstream iss(line);
      int rotid;
      core::Size nb_res = 0;
      while (iss >> rotid) {
        if (rotid > 0) {
          if (nb_res < n_mres) {
            nb_res++;
            rotamer_assignment.push_back(rotid);
          }
          //          std::cout << "True: " << rotid
          //<<
          // std::endl;
        }
      }
    }
  }
  in.close();
  return;
}

void load_CFN_conformation(
    const std::string &filename, core::Size n_mres,
    utility::vector1<utility::vector1<core::Size>> &rotamer_assignments) {
  std::ifstream inmain;
  inmain.open(filename);
  core::Size rotid;
  int count = 1;
  while (true) {
    // std::cout << "loading assignments " << count << std::endl;
    rotamer_assignments.push_back(utility::vector1<core::Size>());
    char buf[99999];
    inmain.getline(buf, 99999);
    std::istringstream in(buf);
    while (in >> rotid && rotamer_assignments.back().size() < n_mres) {
      rotamer_assignments.back().push_back(rotid + 1);
    }
    if (!inmain.good()) break;
    ++count;
    if (count >= 1000) {
      std::cout << "truncating " << filename << " at 1000 entries" << std::endl;
    }
  }
  inmain.close();

  rotamer_assignments.pop_back();  // last one is garbage
  std::cout << "load_CFN_conformation: rotamer_assignments.size() "
            << rotamer_assignments.size() << std::endl;
  for (size_t i = 1; i <= rotamer_assignments.size(); ++i) {
    // std::cout << rotamer_assignments[i].size() << std::endl;
    runtime_assert(rotamer_assignments[1].size() ==
                   rotamer_assignments[i].size());
  }
  return;
}

void HBNet::load_conformation_and_dump_pose(std::string filename) {
  utility::vector1<utility::vector1<core::Size>> rotamer_assignments;
  // core::Size rotid;

  if (basic::options::option[basic::options::OptionKeys::out::
                                 hbnet_print_interaction_graph_to_file]() < 6) {
    rotamer_assignments.resize(1);
    load_SAT_conformation(filename, rotamer_sets_->nmoltenres(),
                          rotamer_assignments[1]);
    // convert to local rotids
    for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres();
         ++mresid) {
      rotamer_assignments[1][mresid] =
          rotamer_assignments[1][mresid] -
          rotamer_sets_->nrotamer_offset_for_moltenres(mresid);
    }
  } else if (basic::options::option
                 [basic::options::OptionKeys::out::
                      hbnet_print_interaction_graph_to_file]() == 6) {
    rotamer_assignments.resize(1);
    load_ILP_conformation(filename, rotamer_sets_->nmoltenres(),
                          rotamer_assignments[1]);
    // convert to local rotids
    for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres();
         ++mresid) {
      rotamer_assignments[1][mresid] =
          rotamer_assignments[1][mresid] -
          rotamer_sets_->nrotamer_offset_for_moltenres(mresid);
    }
  } else {
    load_CFN_conformation(filename, rotamer_sets_->nmoltenres(),
                          rotamer_assignments);
  }

  for (core::Size irotassign = 1; irotassign <= rotamer_assignments.size();
       ++irotassign) {
    auto const &rotamer_assignment = rotamer_assignments[irotassign];
    auto pose = get_orig_pose();

    if (rotamer_assignment.size() != rotamer_sets_->nmoltenres()) {
      std::cout
          << " => Conformation file has not the expected number of rotamer "
             "assignments"
          << std::endl;
      return;
    }

    std::ostringstream out;  // whs: comment this out and the below back in to
                             // restore Conformation.hb output
    // std::ostringstream hbfname;
    // hbfname << "Conformation_" << irotassign << ".hb";
    // std::ofstream out;
    // out.open(hbfname.str());

    // Analyze rotamer assignment: for each molten residue, and for each
    // of its heavy atoms, give the hbonds.
    const core::scoring::hbonds::HBondDatabaseCOP hb_database =
        core::scoring::hbonds::HBondDatabase::get_database();
    const core::scoring::TenANeighborGraph &tenA_neighbor_graph(
        orig_pose_->energies().tenA_neighbor_graph());
    core::Size total_number_Hatoms = 0;
    core::Size total_number_polar = 0;
    core::Size total_number_2polar = 0;
    core::Size total_number_HBonds = 0;

    // -s ../scott_heterodimers/427_0001.pdb -parser:protocol
    // ../xml/hbnet_betah.xml -hbnet_print_interaction_graph_to_file 7 -beta
    // -hbnet_read_conformation_file scott_427_0001_sol.wcsp_clust_30_40.wcsp

    std::vector<
        std::tuple<core::Size, core::Size, core::Size, core::Size, float>>
        all_hbonds;  // inefficient search later.
    for (std::list<EdgeBase *>::const_iterator iter =
             ig_->get_edge_list_begin();
         iter != ig_->get_edge_list_end(); ++iter) {
      PDEdge *pdedge = static_cast<PDEdge *>(*iter);

      core::Size const mres_1 = (*iter)->get_first_node_ind();
      core::Size const mres_2 = (*iter)->get_second_node_ind();

      core::Size const resid_1 = rotamer_sets_->moltenres_2_resid(mres_1);
      core::Size const resid_2 = rotamer_sets_->moltenres_2_resid(mres_2);

      core::Size rotid_1 = rotamer_assignment[mres_1];
      core::Size rotid_2 = rotamer_assignment[mres_2];

      core::Size global_rotid_1 =
          rotamer_sets_->nrotamer_offset_for_moltenres(mres_1) + rotid_1;
      core::Size global_rotid_2 =
          rotamer_sets_->nrotamer_offset_for_moltenres(mres_2) + rotid_2;

      // Filling in hbond creation info
      const core::Size nbrs_1 = tenA_neighbor_graph.get_node(resid_1)
                                    ->num_neighbors_counting_self_static();
      const core::Size nbrs_2 = tenA_neighbor_graph.get_node(resid_2)
                                    ->num_neighbors_counting_self_static();

      core::Real const score = pdedge->get_two_body_energy(rotid_1, rotid_2);
      if (score > clash_threshold_) {
        out << "Error: clash detected M" << mres_1 << "-M" << mres_2
            << std::endl;
      }

      if (score < 0) {
        core::conformation::ResidueCOP rotamer_1 =
            rotamer_sets_->rotamer(global_rotid_1);
        core::conformation::ResidueCOP rotamer_2 =
            rotamer_sets_->rotamer(global_rotid_2);

        core::scoring::hbonds::HBondSet hbond_set;

        core::scoring::hbonds::identify_hbonds_1way(
            *hb_database, *rotamer_1, *rotamer_2, nbrs_1, nbrs_2,
            false,  // bool const evaluate_derivative,
            false,  // bool const exclude_don_bb,
            false,  // bool const exclude_don_bsc,
            false,  // bool const exclude_acc_scb,
            false,  // bool const exclude_acc_sc,
            // output
            hbond_set);

        core::scoring::hbonds::identify_hbonds_1way(
            *hb_database, *rotamer_2, *rotamer_1, nbrs_2, nbrs_1,
            false,  // bool const evaluate_derivative,
            false,  // bool const exclude_don_bb,
            false,  // bool const exclude_don_bsc,
            false,  // bool const exclude_acc_scb,
            false,  // bool const exclude_acc_sc,
            // output
            hbond_set);

        for (core::Size hbond_id = 1; hbond_id <= hbond_set.nhbonds();
             ++hbond_id) {
          core::scoring::hbonds::HBondCOP hbond = hbond_set.hbond_cop(hbond_id);
          bool is_don_bb = hbond->don_hatm_is_backbone();
          bool is_acc_bb = hbond->acc_atm_is_backbone();

          // core::Size const H = hbond->don_hatm();
          core::Size const acc = hbond->acc_atm();
          core::Size H_parent;

          // core::Size don_rotamer;
          // core::Size acc_rotamer;
          core::Size acc_mres;
          core::Size don_mres;

          if (hbond->don_res() == resid_1) {
            // don_rotamer = global_rotid_1;
            don_mres = mres_1;
            H_parent = rotamer_1->atom_base(hbond->don_hatm());
            // acc_rotamer = global_rotid_2;
            acc_mres = mres_2;
          } else {
            // don_rotamer = global_rotid_2;
            don_mres = mres_2;
            H_parent = rotamer_2->atom_base(hbond->don_hatm());
            // acc_rotamer = global_rotid_1;
            acc_mres = mres_1;
          }

          // char don_AA = rotamer_sets_->rotamer(don_rotamer)->name1();
          // char acc_AA = rotamer_sets_->rotamer(acc_rotamer)->name1();

          if ((hbond->energy() <= hydrogen_bond_threshold_) &&
              ((!is_acc_bb) || (!is_don_bb))) {
            if (is_don_bb) {
              total_number_HBonds++;
              auto this_hbond =
                  std::make_tuple(acc_mres, acc, don_mres, 0, hbond->energy());
              all_hbonds.push_back(this_hbond);
              //            out << "M" <<
              // acc_mres << acc_AA << ":" << acc
              //<<
              //"
              //-
              // BB" << std::endl;
            } else if (is_acc_bb) {
              total_number_HBonds++;
              auto this_hbond = std::make_tuple(acc_mres, 0, don_mres, H_parent,
                                                hbond->energy());
              all_hbonds.push_back(this_hbond);
              //            out << "BB - "
              //<< "M" << don_mres << don_AA
              //<<
              //":"
              //<<
              // H_parent << std::endl;
            } else {
              total_number_HBonds += 2;
              auto this_hbond = std::make_tuple(acc_mres, acc, don_mres,
                                                H_parent, hbond->energy());
              all_hbonds.push_back(this_hbond);
              //            out << "M" <<
              // acc_mres << acc_AA << ":" << acc
              //<<
              //"
              //-
              // M" << don_mres << don_AA << ":" << H_parent << std::endl;
            }
          }
        }
      }
    }

    if (!(pose.pdb_info())) {
      pose.pdb_info(core::pose::PDBInfoOP(new core::pose::PDBInfo(pose, true)));
    }

    for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres();
         ++mresid) {
      auto rot_set_for_mres =
          rotamer_sets_->rotamer_set_for_moltenresidue(mresid);
      core::conformation::ResidueCOP rotamer =
          rot_set_for_mres->rotamer(rotamer_assignment[mresid]);
      std::set<core::Size> selected_polars;

      for (core::Size sc_acc : rotamer->accpt_pos_sc()) {
        selected_polars.insert(sc_acc);
      }
      total_number_2polar += (selected_polars.size() > 1);

      for (core::Size polar_H : rotamer->Hpos_polar_sc()) {
        selected_polars.insert(rotamer->atom_base(polar_H));
      }
      total_number_polar += !selected_polars.empty();

      for (core::Size atom : selected_polars) {
        total_number_Hatoms++;
        std::string PDB_res =
            pose.pdb_info()->pose2pdb(rotamer_sets_->moltenres_2_resid(mresid));
        out << PDB_res << rotamer->name1() << ":" << atom << '\t';

        for (const auto &hb : all_hbonds) {
          if ((std::get<0>(hb) == mresid) && (std::get<1>(hb) == atom)) {
            char other_name1 =
                rotamer_sets_->rotamer_set_for_moltenresidue(std::get<2>(hb))
                    ->rotamer(rotamer_assignment[std::get<2>(hb)])
                    ->name1();
            const std::string &PDB_res = pose.pdb_info()->pose2pdb(
                rotamer_sets_->moltenres_2_resid(std::get<2>(hb)));
            out << "[" << PDB_res << other_name1 << ":" << std::get<3>(hb)
                << " @ " << std::get<4>(hb) << "] ";
          }
          if ((std::get<2>(hb) == mresid) && (std::get<3>(hb) == atom)) {
            char other_name1 =
                rotamer_sets_->rotamer_set_for_moltenresidue(std::get<0>(hb))
                    ->rotamer(rotamer_assignment[std::get<0>(hb)])
                    ->name1();
            const std::string &PDB_res = pose.pdb_info()->pose2pdb(
                rotamer_sets_->moltenres_2_resid(std::get<0>(hb)));
            out << "[" << PDB_res << other_name1 << ":" << std::get<1>(hb)
                << " @ " << std::get<4>(hb) << "] ";
          }
        }
        out << std::endl;
      }
    }
    out << "Total of " << total_number_polar << ":" << total_number_2polar
        << " polar AAs, " << total_number_Hatoms << " heavy atoms, "
        << total_number_HBonds << " half hbonds." << std::endl;
    // out.close();  // whs: add this back for Conformation.hb output!!

    // assume packer, same rotamer_sets, ig are there
    for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres();
         ++mresid) {
      auto rot = rotamer_sets_->rotamer_set_for_moltenresidue(mresid)->rotamer(
          rotamer_assignment[mresid]);
      pose.replace_residue(rotamer_sets_->moltenres_2_resid(mresid), *rot,
                           false);
      pose.pdb_info()->add_reslabel(rot->seqpos(), "HBSAT");
    }

    ///////////// dump file ////////////////
    std::string basefname = "HBSAT_";
    auto const &opt_s =
        basic::options::option[basic::options::OptionKeys::in::file::s]();
    if (opt_s.size() == 1) {
      // the following is a real beauty...
      //     utility::file_basename strips the directory path
      //     utility::file::file_basename strips the file extension
      basefname +=
          utility::file::file_basename(utility::file_basename(filename)) + "_";
      basefname +=
          utility::file::file_basename(utility::file_basename(opt_s[1])) + "_";
    }
    std::ostringstream fname;
    fname << basefname << ObjexxFCL::lead_zero_Fstring_of(irotassign, 4)
          << ".pdb";
    std::cout << "dumping pose: " << fname.str() << std::endl;
    pose.dump_pdb(fname.str());

  }  // end irtoassign
}

void HBNet::print_interaction_graph_to_file(std::string filename) {
  std::ofstream out;
  out.open(filename);

  for (core::Size rotid = 1; rotid <= rotamer_sets_->nrotamers(); ++rotid) {
    core::conformation::ResidueCOP rotamer = rotamer_sets_->rotamer(rotid);
    core::Size const moltenresid = rotamer_sets_->moltenres_for_rotamer(rotid);
    out << rotid << "\t" << rotamer->name1() << "\t" << moltenresid << "\t"
        << orig_pose_->chain(rotamer_sets_->moltenres_2_resid(moltenresid));

    utility::vector1<core::Size> sc_acceptors = rotamer->accpt_pos_sc();
    for (core::Size sc_acc : sc_acceptors) {
      out << "\t"
          << "A_" << sc_acc;
    }

    utility::vector1<core::Size> sc_polar_Hs = rotamer->Hpos_polar_sc();
    for (core::Size polar_H : sc_polar_Hs) {
      out << "\tH_" << polar_H << "_" << rotamer->atom_base(polar_H);
    }

    out << std::endl;
  }

  const core::scoring::hbonds::HBondDatabaseCOP hb_database =
      core::scoring::hbonds::HBondDatabase::get_database();
  const core::scoring::TenANeighborGraph &tenA_neighbor_graph(
      orig_pose_->energies().tenA_neighbor_graph());

  for (std::list<EdgeBase *>::const_iterator iter = ig_->get_edge_list_begin();
       iter != ig_->get_edge_list_end(); ++iter) {
    PDEdge *pdedge = static_cast<PDEdge *>(*iter);

    core::Size const mres_1 = (*iter)->get_first_node_ind();
    core::Size const mres_2 = (*iter)->get_second_node_ind();

    core::Size const resid_1 = rotamer_sets_->moltenres_2_resid(mres_1);
    core::Size const resid_2 = rotamer_sets_->moltenres_2_resid(mres_2);

    const core::Size si = rotamer_sets_->nrotamers_for_moltenres(mres_1);
    const core::Size sj = rotamer_sets_->nrotamers_for_moltenres(mres_2);

    const core::Size nbrs_1 = tenA_neighbor_graph.get_node(resid_1)
                                  ->num_neighbors_counting_self_static();
    const core::Size nbrs_2 = tenA_neighbor_graph.get_node(resid_2)
                                  ->num_neighbors_counting_self_static();

    for (core::Size ii = 1; ii <= si; ++ii) {
      for (core::Size jj = 1; jj <= sj; ++jj) {
        core::Real const score = pdedge->get_two_body_energy(ii, jj);
        core::Size const global_rotamer_ii =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_1) + ii;
        core::Size const global_rotamer_jj =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_2) + jj;
        if (score > 0) {
          out << "C\t" << global_rotamer_ii << "\t" << global_rotamer_jj << "\t"
              << score << std::endl;
        } else if (score < 0) {
          core::conformation::ResidueCOP rotamer_ii =
              rotamer_sets_->rotamer(global_rotamer_ii);
          core::conformation::ResidueCOP rotamer_jj =
              rotamer_sets_->rotamer(global_rotamer_jj);

          core::scoring::hbonds::HBondSet hbond_set;

          core::scoring::hbonds::identify_hbonds_1way(
              *hb_database, *rotamer_ii, *rotamer_jj, nbrs_1, nbrs_2,
              false,  // bool const evaluate_derivative,
              false,  // bool const exclude_don_bb,
              false,  // bool const exclude_don_bsc,
              false,  // bool const exclude_acc_scb,
              false,  // bool const exclude_acc_sc,
              // output
              hbond_set);

          core::scoring::hbonds::identify_hbonds_1way(
              *hb_database, *rotamer_jj, *rotamer_ii, nbrs_2, nbrs_1,
              false,  // bool const evaluate_derivative,
              false,  // bool const exclude_don_bb,
              false,  // bool const exclude_don_bsc,
              false,  // bool const exclude_acc_scb,
              false,  // bool const exclude_acc_sc,
              // output
              hbond_set);

          for (core::Size hbond_id = 1; hbond_id <= hbond_set.nhbonds();
               ++hbond_id) {
            core::scoring::hbonds::HBondCOP hbond =
                hbond_set.hbond_cop(hbond_id);
            if (hbond->don_hatm_is_backbone() && hbond->acc_atm_is_backbone()) {
              continue;
            }

            core::Size const H = hbond->don_hatm();
            core::Size H_parent;
            core::Size const acc = hbond->acc_atm();

            core::Size don_rotamer;
            core::Size acc_rotamer;

            if (hbond->don_res() == resid_1) {
              don_rotamer = global_rotamer_ii;
              acc_rotamer = global_rotamer_jj;
              H_parent = rotamer_ii->atom_base(hbond->don_hatm());
            } else {
              don_rotamer = global_rotamer_jj;
              acc_rotamer = global_rotamer_ii;
              H_parent = rotamer_jj->atom_base(hbond->don_hatm());
            }

            out << "H\t" << global_rotamer_ii << "\t" << global_rotamer_jj
                << "\t" << hbond->energy() << "\t" << don_rotamer << "_" << H
                << "\t" << don_rotamer << "_" << H_parent << "\t" << acc_rotamer
                << "_" << acc << std::endl;
          }
        }
        // out << "\t" << pdedge->get_two_body_energy(ii, jj);
      }
    }
    // out << std::endl;
  }

  out.close();
}

void HBNet::print_SAT_model_to_file(std::string filename, bool prefpolar,
                                    bool pref2polar) {
  std::stringstream buffer;

  const bool comment = true;  // generate comments in the SAT formula
  const unsigned int nresid = rotamer_sets_->nmoltenres();
  const core::Size top =
      (prefpolar ? rotamer_sets_->nmoltenres() : 0) +
      (pref2polar ? rotamer_sets_->nmoltenres() * rotamer_sets_->nmoltenres()
                  : 0) +
      1;
  // for PWMaxSAT, we cannot violate this amount of weight (full hydrophobic)

  bool pwmax_sat =
      (prefpolar || pref2polar);  // if we have preferences, we must use maxSAT

  core::Size nb_cnf = 0;

  // First CNF formula for hbnets. This implements the Boyken's
  // conjecture as follows: we don't distinguish the donors of a polar
  // atom and the polar atom itself and just require that this
  // "meta-atom" has an h-bond. The heavy atom is used as the
  // representative atom in all cases.

  // To allow for PWMaxSAT optimization, mutations to Glycine have to be
  // allowed in the des_residues, are penalized in pwmax_sat mode (forbidden in
  // SAT).

  const core::Size rot2AA_offset =
      1;  // offset to the first var that selects AA at a residue
  core::Size rel_hC_offset =
      rotamer_sets_->nrotamers();  // offset to add to previous offset to access
                                   // h_C variables
  buffer << "c " << rotamer_sets_->nrotamers() << " rotamers" << std::endl;

  std::vector<std::map<char, core::Size>> mres_AA2var_idx;  // maps a molten
                                                            // residue number
                                                            // and designable AA
                                                            // to a s_i^AA
                                                            // variable
  // for each molten residue, maps AA + polar atom number to a set of rotamers
  // that create a hbond for it (S1)
  std::vector<std::map<std::pair<char, core::Size>, std::set<core::Size>>>
      res_pola2hbonding_S1;
  // for each molten residue, maps AA + polar atom number to the variable idx
  // that represents condition for hbond creation by S2
  std::vector<std::map<std::pair<char, core::Size>, std::set<core::Size>>>
      res_pola2hbonding_S2C;

  // rotamers variable. Named by rotid in the full rotamer enumeration.
  // First generate the ALO/AMO clauses. A set of clauses by residue.
  for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres(); ++mresid) {
    auto rot_set_for_mres =
        rotamer_sets_->rotamer_set_for_moltenresidue(mresid);
    // AMO clauses
    if (comment)
      buffer << "c "
             << "AMO clauses for residue M" << mresid << "/"
             << rotamer_sets_->moltenres_2_resid(mresid) << std::endl;
    for (core::Size rotid1 = 1; rotid1 < rot_set_for_mres->num_rotamers();
         ++rotid1) {
      for (core::Size rotid2 = rotid1 + 1;
           rotid2 <= rot_set_for_mres->num_rotamers(); ++rotid2) {
        nb_cnf++;
        if (pwmax_sat) buffer << top << " ";
        buffer << "-" << rotamer_sets_->moltenres_rotid_2_rotid(mresid, rotid1)
               << " -" << rotamer_sets_->moltenres_rotid_2_rotid(mresid, rotid2)
               << " 0" << std::endl;
      }
    }

    // ALO clause
    if (comment)
      buffer << "c "
             << "ALO clause for residue M" << mresid << "/"
             << rotamer_sets_->moltenres_2_resid(mresid) << std::endl;
    if (pwmax_sat) buffer << top << " ";
    for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
         ++rotid) {
      buffer << rotamer_sets_->moltenres_rotid_2_rotid(mresid, rotid) << " ";
    }
    nb_cnf++;
    buffer << "0" << std::endl;

    // rotid to AA clauses (if rotamer rotid is selected, then the corresponding
    // AA is used).
    std::map<char, core::Size> AA2var;
    for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
         ++rotid) {
      core::conformation::ResidueCOP rotamer = rot_set_for_mres->rotamer(rotid);
      if ((rotamer->Hpos_polar_sc().size() != 0) ||
          (rotamer->accpt_pos_sc().size() != 0)) {
        if (AA2var.find(rotamer->name1()) == AA2var.end()) {
          rel_hC_offset++;  // we need a new variable
          AA2var[rotamer->name1()] = rot2AA_offset + rel_hC_offset;
          if (comment)
            buffer << "c variable " << rot2AA_offset + rel_hC_offset << " for M"
                   << mresid << "/" << rotamer_sets_->moltenres_2_resid(mresid)
                   << " (" << rotamer->name1() << ")" << std::endl;
        }
        if (pwmax_sat) buffer << top << " ";
        buffer << "-" << rotamer_sets_->moltenres_rotid_2_rotid(mresid, rotid)
               << " " << rot2AA_offset + rel_hC_offset << " 0" << std::endl;
        nb_cnf++;
      }
    }
    mres_AA2var_idx.push_back(AA2var);

    // initializing S1 and S2C sets with empy sets for polar atoms
    // forbidding/penalizing hydrophobic ones
    std::map<std::pair<char, core::Size>, std::set<core::Size>> S1;
    std::map<std::pair<char, core::Size>, std::set<core::Size>> S2C;
    for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
         ++rotid) {
      core::conformation::ResidueCOP rotamer = rot_set_for_mres->rotamer(rotid);
      bool is_polar = ((rotamer->Hpos_polar_sc().size() != 0) ||
                       (rotamer->accpt_pos_sc().size() != 0));
      bool is_2polar = (rotamer->accpt_pos_sc().size() > 1);

      if (is_polar) {
        // For each AA, we have a set of acceptor and donor atoms
        for (core::Size sc_acc : rotamer->accpt_pos_sc()) {
          std::pair<char, core::Size> polar_atomid =
              std::make_pair(rotamer->name1(), sc_acc);
          if (S1.find(polar_atomid) == S1.end()) {
            std::set<core::Size> set1;
            std::set<core::Size> set2;
            S1[polar_atomid] = set1;
            S2C[polar_atomid] = set2;
          }
        }

        for (core::Size polar_H : rotamer->Hpos_polar_sc()) {
          core::Size H_parent = rotamer->atom_base(polar_H);
          std::pair<char, core::Size> polar_atomid =
              std::make_pair(rotamer->name1(), H_parent);
          if (S1.find(polar_atomid) == S1.end()) {
            std::set<core::Size> set1;
            std::set<core::Size> set2;
            S1[polar_atomid] = set1;
            S2C[polar_atomid] = set2;
          }
        }
      }

      if (!is_polar && prefpolar) {  // forbid/penalize hydrophobic rotamers
        if (comment)
          buffer << "c non-polar rotamer for M" << mresid << " undesirable!"
                 << std::endl;
        if (pwmax_sat) {
          buffer << 1 << " ";
        }
        buffer << "-" << rotamer_sets_->moltenres_rotid_2_rotid(mresid, rotid)
               << " 0" << std::endl;
        nb_cnf++;
      }

      if (!is_2polar && pref2polar) {
        if (comment)
          buffer << "c non-2polar rotamer for M" << mresid << " undesirable!"
                 << std::endl;
        if (pwmax_sat) {
          buffer << nresid << " ";
        }
        buffer << "-" << rotamer_sets_->moltenres_rotid_2_rotid(mresid, rotid)
               << " 0" << std::endl;
        nb_cnf++;
      }
    }
    res_pola2hbonding_S1.push_back(S1);
    res_pola2hbonding_S2C.push_back(S2C);
  }

  // Clash clauses + filling hbond creation sets
  const core::scoring::hbonds::HBondDatabaseCOP hb_database =
      core::scoring::hbonds::HBondDatabase::get_database();
  const core::scoring::TenANeighborGraph &tenA_neighbor_graph(
      orig_pose_->energies().tenA_neighbor_graph());
  core::Size final_offset =
      rot2AA_offset + rel_hC_offset;  // where new variables will be created

  for (std::list<EdgeBase *>::const_iterator iter = ig_->get_edge_list_begin();
       iter != ig_->get_edge_list_end(); ++iter) {
    PDEdge *pdedge = static_cast<PDEdge *>(*iter);

    core::Size const mres_1 = (*iter)->get_first_node_ind();
    core::Size const mres_2 = (*iter)->get_second_node_ind();

    core::Size const resid_1 = rotamer_sets_->moltenres_2_resid(mres_1);
    core::Size const resid_2 = rotamer_sets_->moltenres_2_resid(mres_2);

    const core::Size si = rotamer_sets_->nrotamers_for_moltenres(mres_1);
    const core::Size sj = rotamer_sets_->nrotamers_for_moltenres(mres_2);

    // Sterical clash clauses
    if (comment)
      buffer << "c Clash clauses for M" << mres_1 << "/" << resid_1 << ":M"
             << mres_2 << "/" << resid_2 << std::endl;
    for (core::Size ii = 1; ii <= si; ++ii) {
      for (core::Size jj = 1; jj <= sj; ++jj) {
        core::Real const score = pdedge->get_two_body_energy(ii, jj);
        core::Size const global_rotamer_ii =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_1) + ii;
        core::Size const global_rotamer_jj =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_2) + jj;
        if (score >= clash_threshold_) {
          if (pwmax_sat) buffer << top << " ";
          buffer << "-" << global_rotamer_ii << " -" << global_rotamer_jj
                 << " 0" << std::endl;
          nb_cnf++;
        }
      }
    }

    // Filling in hbond creation info
    const core::Size nbrs_1 = tenA_neighbor_graph.get_node(resid_1)
                                  ->num_neighbors_counting_self_static();
    const core::Size nbrs_2 = tenA_neighbor_graph.get_node(resid_2)
                                  ->num_neighbors_counting_self_static();

    for (core::Size ii = 1; ii <= si; ++ii) {
      for (core::Size jj = 1; jj <= sj; ++jj) {
        core::Real const score = pdedge->get_two_body_energy(ii, jj);
        core::Size const global_rotamer_ii =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_1) + ii;
        core::Size const global_rotamer_jj =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_2) + jj;

        if (score < 0) {
          core::conformation::ResidueCOP rotamer_ii =
              rotamer_sets_->rotamer(global_rotamer_ii);
          core::conformation::ResidueCOP rotamer_jj =
              rotamer_sets_->rotamer(global_rotamer_jj);

          core::scoring::hbonds::HBondSet hbond_set;

          core::scoring::hbonds::identify_hbonds_1way(
              *hb_database, *rotamer_ii, *rotamer_jj, nbrs_1, nbrs_2,
              false,  // bool const evaluate_derivative,
              false,  // bool const exclude_don_bb,
              false,  // bool const exclude_don_bsc,
              false,  // bool const exclude_acc_scb,
              false,  // bool const exclude_acc_sc,
              // output
              hbond_set);

          core::scoring::hbonds::identify_hbonds_1way(
              *hb_database, *rotamer_jj, *rotamer_ii, nbrs_2, nbrs_1,
              false,  // bool const evaluate_derivative,
              false,  // bool const exclude_don_bb,
              false,  // bool const exclude_don_bsc,
              false,  // bool const exclude_acc_scb,
              false,  // bool const exclude_acc_sc,
              // output
              hbond_set);

          for (core::Size hbond_id = 1; hbond_id <= hbond_set.nhbonds();
               ++hbond_id) {
            core::scoring::hbonds::HBondCOP hbond =
                hbond_set.hbond_cop(hbond_id);
            bool is_don_bb = hbond->don_hatm_is_backbone();
            bool is_acc_bb = hbond->acc_atm_is_backbone();

            if ((is_acc_bb && is_don_bb) ||
                (hbond->energy() > hydrogen_bond_threshold_)) {
              continue;  // not interested in weak or BB to BB hbonds
            }

            core::Size const H = hbond->don_hatm();
            core::Size const acc = hbond->acc_atm();
            core::Size H_parent;

            core::Size don_rotamer;
            core::Size acc_rotamer;
            core::Size acc_mres;
            core::Size don_mres;

            if (hbond->don_res() == resid_1) {
              don_rotamer = global_rotamer_ii;
              don_mres = mres_1;
              H_parent = rotamer_ii->atom_base(hbond->don_hatm());
              acc_rotamer = global_rotamer_jj;
              acc_mres = mres_2;
            } else {
              don_rotamer = global_rotamer_jj;
              don_mres = mres_2;
              H_parent = rotamer_jj->atom_base(hbond->don_hatm());
              acc_rotamer = global_rotamer_ii;
              acc_mres = mres_1;
            }

            char don_AA = rotamer_sets_->rotamer(don_rotamer)->name1();
            char acc_AA = rotamer_sets_->rotamer(acc_rotamer)->name1();

            if (is_don_bb) {  // only record information on acceptor in S1
              std::pair<char, core::Size> acc_atomid =
                  std::make_pair(acc_AA, acc);
              if (res_pola2hbonding_S1[acc_mres - 1].find(acc_atomid) ==
                  res_pola2hbonding_S1[acc_mres - 1].end())
                std::cerr << "S1 for res. " << acc_mres << ", atom " << acc_AA
                          << ":" << acc << " is not there" << std::endl;
              else
                res_pola2hbonding_S1[acc_mres - 1][acc_atomid].insert(
                    acc_rotamer);
            }
            if (is_acc_bb) {  // only record information on donor's parent in S1
              std::pair<char, core::Size> don_atomid =
                  std::make_pair(don_AA, H_parent);
              if (res_pola2hbonding_S1[don_mres - 1].find(don_atomid) ==
                  res_pola2hbonding_S1[don_mres - 1].end())
                std::cerr << "S1 for res. " << don_mres << ", atom " << don_AA
                          << ":" << H << " is not there" << std::endl;
              else
                res_pola2hbonding_S1[don_mres - 1][don_atomid].insert(
                    don_rotamer);
            }
            if (!is_acc_bb && !is_don_bb) {  // record information on donor and
                                             // acceptor in S2C
              std::pair<char, core::Size> don_atomid =
                  std::make_pair(don_AA, H_parent);
              std::pair<char, core::Size> acc_atomid =
                  std::make_pair(acc_AA, acc);
              if (comment)
                buffer << "c Variable " << ++final_offset << " for hbond "
                       << don_mres << don_AA << "_" << H_parent << "-" << H
                       << ":" << acc_mres << acc_AA << "_" << acc
                       << " via rotamers (" << don_rotamer << ", "
                       << acc_rotamer << ")" << std::endl;
              if (res_pola2hbonding_S2C[acc_mres - 1].find(acc_atomid) ==
                  res_pola2hbonding_S2C[acc_mres - 1].end())
                std::cerr << "S2C for res. " << acc_mres << ", atom " << acc_AA
                          << ":" << acc << " is not there !" << std::endl;
              else
                res_pola2hbonding_S2C[acc_mres - 1][acc_atomid].insert(
                    final_offset);
              if (res_pola2hbonding_S2C[don_mres - 1].find(don_atomid) ==
                  res_pola2hbonding_S2C[don_mres - 1].end())
                std::cerr << "S2C for res. " << don_mres << ", atom " << don_AA
                          << ":" << H << " is not there !" << std::endl;
              else
                res_pola2hbonding_S2C[don_mres - 1][don_atomid].insert(
                    final_offset);
              if (pwmax_sat) buffer << top << " ";
              buffer << don_rotamer << " -" << final_offset << " 0"
                     << std::endl;
              nb_cnf++;
              if (pwmax_sat) buffer << top << " ";
              buffer << acc_rotamer << " -" << final_offset << " 0"
                     << std::endl;
              nb_cnf++;
            }
          }
        }
      }
    }
  }

  // Final set of clauses: at least one way of making an hbond for
  // each heavy polar atom (either as an acceptor or through attached
  // hydrogen, as a donor) must be active if the corresponding AA is
  // used at given position.
  for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres(); ++mresid) {
    auto rot_set_for_mres =
        rotamer_sets_->rotamer_set_for_moltenresidue(mresid);
    char AA = '1';  // unused as AA char
    for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
         ++rotid) {
      core::conformation::ResidueCOP rotamer = rot_set_for_mres->rotamer(rotid);

      if (rotamer->name1() == AA)
        continue;  // skip, all rotamers of a given AA have all the same polar
                   // atoms
      AA = rotamer->name1();
      if (comment)
        buffer << "c Considering rotamer " << rotid << " of " << AA
               << " at position " << mresid << std::endl;

      std::set<std::pair<char, core::Size>> selected_polars;  // polar atoms
                                                              // acting as
                                                              // acceptor or as
                                                              // donor through
                                                              // attached H
      for (core::Size sc_acc : rotamer->accpt_pos_sc()) {
        selected_polars.insert(std::make_pair(AA, sc_acc));
      }
      for (core::Size polar_H : rotamer->Hpos_polar_sc()) {
        selected_polars.insert(std::make_pair(AA, rotamer->atom_base(polar_H)));
      }

      for (std::pair<char, core::Size> polar_atomid : selected_polars) {
        if (pwmax_sat) buffer << top << " ";
        buffer << "-" << mres_AA2var_idx[mresid - 1][AA] << " ";
        for (core::Size rot1 : res_pola2hbonding_S1[mresid - 1][polar_atomid]) {
          buffer << rot1 << " ";
        }
        for (core::Size var_idx :
             res_pola2hbonding_S2C[mresid - 1][polar_atomid]) {
          buffer << var_idx << " ";
        }
        nb_cnf++;
        buffer << " 0" << std::endl;
      }
    }
  }

  std::ofstream out;
  out.open(filename);
  out << "p " << (pwmax_sat ? "wcnf " : "cnf ") << final_offset << " "
      << nb_cnf;
  if (pwmax_sat) out << " " << top;
  out << std::endl << buffer.rdbuf();
  out.close();
}

void HBNet::print_ILP_model_to_file(std::string filename) {
  std::ofstream out;
  out.open(filename);
  const bool comment = true;  // generate comments in the SAT formula

  // core::Size top = rotamer_sets_->nmoltenres() +
  // 1;  // for ILP, we cannot violate this amount of weight
  const core::Size rot2AA_offset =
      1;  // offset to the first var that selects AA at a residue
  core::Size rel_hC_offset =
      rotamer_sets_->nrotamers();  // offset to add to previous offset to access
                                   // h_C variables

  std::vector<std::map<char, core::Size>> mres_AA2var_idx;  // maps a molten
                                                            // residue number
                                                            // and designable AA
                                                            // to a s_i^AA
                                                            // variable
  // for each molten residue, maps AA + polar atom number to a set of rotamers
  // that create a hbond for it (S1)
  std::vector<std::map<std::pair<char, core::Size>, std::set<core::Size>>>
      res_pola2hbonding_S1;
  // for each molten residue, maps AA + polar atom number to the variable idx
  // that represents condition for hbond creation by S2
  std::vector<std::map<std::pair<char, core::Size>, std::set<core::Size>>>
      res_pola2hbonding_S2C;

  // Criteria: penalize all non polar side-chains.
  out << "Minimize " << std::endl << "  ";
  for (core::Size rotid = 1; rotid <= rotamer_sets_->nrotamers(); ++rotid) {
    core::conformation::ResidueCOP rotamer = rotamer_sets_->rotamer(rotid);
    if ((rotamer->Hpos_polar_sc().size() == 0) &&
        (rotamer->accpt_pos_sc().size() == 0))
      out << "x" << rotid << " + ";
  }
  out << "0" << std::endl;

  out << "Subject To" << std::endl;
  // rotamers variable. Named by rotid in the full rotamer enumeration.
  // First generate the domain constraints, one per residue.
  for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres(); ++mresid) {
    auto rot_set_for_mres =
        rotamer_sets_->rotamer_set_for_moltenresidue(mresid);
    if (comment)
      out << "\\ Domain constraint for residue M" << mresid << "/"
          << rotamer_sets_->moltenres_2_resid(mresid) << std::endl;
    for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
         ++rotid) {
      out << "x" << rotamer_sets_->moltenres_rotid_2_rotid(mresid, rotid);
      if (rotid < rot_set_for_mres->num_rotamers()) {
        out << " + ";
      }
    }
    out << " = 1" << std::endl;

    // rotid to AA clauses (if rotamer rotid is selected, then the corresponding
    // AA is used).
    std::map<char, core::Size> AA2var;
    for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
         ++rotid) {
      core::conformation::ResidueCOP rotamer = rot_set_for_mres->rotamer(rotid);
      if (AA2var.find(rotamer->name1()) == AA2var.end()) {
        rel_hC_offset++;  // we need a new variable
        AA2var[rotamer->name1()] = rot2AA_offset + rel_hC_offset;
        if (comment)
          out << "\\ Variable x" << rot2AA_offset + rel_hC_offset << " for M"
              << mresid << "/" << rotamer_sets_->moltenres_2_resid(mresid)
              << " (" << rotamer->name1() << ")" << std::endl;
      }
      out << "x" << rot2AA_offset + rel_hC_offset << " - x"
          << rotamer_sets_->moltenres_rotid_2_rotid(mresid, rotid) << " >= 0"
          << std::endl;
    }
    mres_AA2var_idx.push_back(AA2var);

    // dealing with polar atoms, initializing S1 and S2C sets with empy sets
    std::map<std::pair<char, core::Size>, std::set<core::Size>> S1;
    std::map<std::pair<char, core::Size>, std::set<core::Size>> S2C;
    for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
         ++rotid) {
      core::conformation::ResidueCOP rotamer = rot_set_for_mres->rotamer(rotid);

      // For each AA, we have a set of acceptor and donor atoms
      for (core::Size sc_acc : rotamer->accpt_pos_sc()) {
        std::pair<char, core::Size> polar_atomid =
            std::make_pair(rotamer->name1(), sc_acc);
        if (S1.find(polar_atomid) == S1.end()) {
          std::set<core::Size> set1;
          std::set<core::Size> set2;
          S1[polar_atomid] = set1;
          S2C[polar_atomid] = set2;
        }
      }

      for (core::Size polar_H : rotamer->Hpos_polar_sc()) {
        core::Size H_parent = rotamer->atom_base(polar_H);
        std::pair<char, core::Size> polar_atomid =
            std::make_pair(rotamer->name1(), H_parent);
        if (S1.find(polar_atomid) == S1.end()) {
          std::set<core::Size> set1;
          std::set<core::Size> set2;
          S1[polar_atomid] = set1;
          S2C[polar_atomid] = set2;
        }
      }
    }
    res_pola2hbonding_S1.push_back(S1);
    res_pola2hbonding_S2C.push_back(S2C);
  }

  // Clash clauses + filling hbond creation sets
  const core::scoring::hbonds::HBondDatabaseCOP hb_database =
      core::scoring::hbonds::HBondDatabase::get_database();
  const core::scoring::TenANeighborGraph &tenA_neighbor_graph(
      orig_pose_->energies().tenA_neighbor_graph());
  core::Size final_offset =
      rot2AA_offset + rel_hC_offset;  // where new variables will be created

  for (std::list<EdgeBase *>::const_iterator iter = ig_->get_edge_list_begin();
       iter != ig_->get_edge_list_end(); ++iter) {
    PDEdge *pdedge = static_cast<PDEdge *>(*iter);

    core::Size const mres_1 = (*iter)->get_first_node_ind();
    core::Size const mres_2 = (*iter)->get_second_node_ind();

    core::Size const resid_1 = rotamer_sets_->moltenres_2_resid(mres_1);
    core::Size const resid_2 = rotamer_sets_->moltenres_2_resid(mres_2);

    const core::Size si = rotamer_sets_->nrotamers_for_moltenres(mres_1);
    const core::Size sj = rotamer_sets_->nrotamers_for_moltenres(mres_2);

    // Sterical clash clauses
    if (comment)
      out << "\\ Clash constraints for M" << mres_1 << "/" << resid_1 << ":M"
          << mres_2 << "/" << resid_2 << std::endl;
    for (core::Size ii = 1; ii <= si; ++ii) {
      for (core::Size jj = 1; jj <= sj; ++jj) {
        core::Real const score = pdedge->get_two_body_energy(ii, jj);
        core::Size const global_rotamer_ii =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_1) + ii;
        core::Size const global_rotamer_jj =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_2) + jj;
        if (score > clash_threshold_) {
          out << "x" << global_rotamer_ii << " + x" << global_rotamer_jj
              << " <= 1" << std::endl;
        }
      }
    }

    // Filling in hbond creation info
    const core::Size nbrs_1 = tenA_neighbor_graph.get_node(resid_1)
                                  ->num_neighbors_counting_self_static();
    const core::Size nbrs_2 = tenA_neighbor_graph.get_node(resid_2)
                                  ->num_neighbors_counting_self_static();

    for (core::Size ii = 1; ii <= si; ++ii) {
      for (core::Size jj = 1; jj <= sj; ++jj) {
        core::Real const score = pdedge->get_two_body_energy(ii, jj);
        core::Size const global_rotamer_ii =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_1) + ii;
        core::Size const global_rotamer_jj =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_2) + jj;

        if (score < 0) {
          core::conformation::ResidueCOP rotamer_ii =
              rotamer_sets_->rotamer(global_rotamer_ii);
          core::conformation::ResidueCOP rotamer_jj =
              rotamer_sets_->rotamer(global_rotamer_jj);

          core::scoring::hbonds::HBondSet hbond_set;

          core::scoring::hbonds::identify_hbonds_1way(
              *hb_database, *rotamer_ii, *rotamer_jj, nbrs_1, nbrs_2,
              false,  // bool const evaluate_derivative,
              false,  // bool const exclude_don_bb,
              false,  // bool const exclude_don_bsc,
              false,  // bool const exclude_acc_scb,
              false,  // bool const exclude_acc_sc,
              // output
              hbond_set);

          core::scoring::hbonds::identify_hbonds_1way(
              *hb_database, *rotamer_jj, *rotamer_ii, nbrs_2, nbrs_1,
              false,  // bool const evaluate_derivative,
              false,  // bool const exclude_don_bb,
              false,  // bool const exclude_don_bsc,
              false,  // bool const exclude_acc_scb,
              false,  // bool const exclude_acc_sc,
              // output
              hbond_set);

          for (core::Size hbond_id = 1; hbond_id <= hbond_set.nhbonds();
               ++hbond_id) {
            core::scoring::hbonds::HBondCOP hbond =
                hbond_set.hbond_cop(hbond_id);
            bool is_don_bb = hbond->don_hatm_is_backbone();
            bool is_acc_bb = hbond->acc_atm_is_backbone();

            if ((is_acc_bb && is_don_bb) ||
                (hbond->energy() > hydrogen_bond_threshold_)) {
              continue;  // not interested in weak or BB to BB hbonds
            }

            core::Size const H = hbond->don_hatm();
            core::Size const acc = hbond->acc_atm();
            core::Size H_parent;

            core::Size don_rotamer;
            core::Size acc_rotamer;
            core::Size acc_mres;
            core::Size don_mres;

            if (hbond->don_res() == resid_1) {
              don_rotamer = global_rotamer_ii;
              don_mres = mres_1;
              H_parent = rotamer_ii->atom_base(hbond->don_hatm());
              acc_rotamer = global_rotamer_jj;
              acc_mres = mres_2;
            } else {
              don_rotamer = global_rotamer_jj;
              don_mres = mres_2;
              H_parent = rotamer_jj->atom_base(hbond->don_hatm());
              acc_rotamer = global_rotamer_ii;
              acc_mres = mres_1;
            }

            char don_AA = rotamer_sets_->rotamer(don_rotamer)->name1();
            char acc_AA = rotamer_sets_->rotamer(acc_rotamer)->name1();

            if (is_don_bb) {  // only record information on acceptor in S1
              std::pair<char, core::Size> acc_atomid =
                  std::make_pair(acc_AA, acc);
              if (res_pola2hbonding_S1[acc_mres - 1].find(acc_atomid) ==
                  res_pola2hbonding_S1[acc_mres - 1].end())
                std::cerr << "S1 for res. " << acc_mres << ", atom " << acc_AA
                          << ":" << acc << " is not there" << std::endl;
              else
                res_pola2hbonding_S1[acc_mres - 1][acc_atomid].insert(
                    acc_rotamer);
            }

            if (is_acc_bb) {  // only record information on donor's parent in S1
              std::pair<char, core::Size> don_atomid =
                  std::make_pair(don_AA, H_parent);
              if (res_pola2hbonding_S1[don_mres - 1].find(don_atomid) ==
                  res_pola2hbonding_S1[don_mres - 1].end())
                std::cerr << "S1 for res. " << don_mres << ", atom " << don_AA
                          << ":" << H << " is not there" << std::endl;
              else
                res_pola2hbonding_S1[don_mres - 1][don_atomid].insert(
                    don_rotamer);
            }

            if (!is_acc_bb && !is_don_bb) {  // record information on donor and
                                             // acceptor in S2C
              std::pair<char, core::Size> don_atomid =
                  std::make_pair(don_AA, H_parent);
              std::pair<char, core::Size> acc_atomid =
                  std::make_pair(acc_AA, acc);
              final_offset++;  // we need a new variable
              if (comment)
                out << "\\ Variable x" << final_offset << " for hbond "
                    << don_mres << don_AA << "_" << H_parent << "-" << H << ":"
                    << acc_mres << acc_AA << "_" << acc << " via rotamers ("
                    << don_rotamer << ", " << acc_rotamer << ")" << std::endl;
              if (res_pola2hbonding_S2C[acc_mres - 1].find(acc_atomid) ==
                  res_pola2hbonding_S2C[acc_mres - 1].end())
                std::cerr << "S2C for res. " << acc_mres << ", atom " << acc_AA
                          << ":" << acc << " is not there !" << std::endl;
              else
                res_pola2hbonding_S2C[acc_mres - 1][acc_atomid].insert(
                    final_offset);
              if (res_pola2hbonding_S2C[don_mres - 1].find(don_atomid) ==
                  res_pola2hbonding_S2C[don_mres - 1].end())
                std::cerr << "S2C for res. " << don_mres << ", atom " << don_AA
                          << ":" << H << " is not there !" << std::endl;
              else
                res_pola2hbonding_S2C[don_mres - 1][don_atomid].insert(
                    final_offset);
              out << "x" << don_rotamer << " - x" << final_offset << " >= 0"
                  << std::endl;
              out << "x" << acc_rotamer << " - x" << final_offset << " >= 0"
                  << std::endl;
            }
          }
        }
      }
    }
  }

  // Final set of constraints: at least one way of making an hbond for
  // each heavy polar atom (either as an acceptor or through attached
  // hydrogen, as a donor) must be active if the corresponding AA is
  // used at given position.
  for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres(); ++mresid) {
    auto rot_set_for_mres =
        rotamer_sets_->rotamer_set_for_moltenresidue(mresid);
    char AA = '1';  // unused as AA char
    for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
         ++rotid) {
      core::conformation::ResidueCOP rotamer = rot_set_for_mres->rotamer(rotid);
      if (rotamer->name1() == AA)
        continue;  // skip, all rotamers of a given AA have all the same polar
                   // atoms
      AA = rotamer->name1();
      if (comment)
        out << "\\ Considering rotamer " << rotid << " of " << AA
            << " at position " << mresid << std::endl;

      std::set<std::pair<char, core::Size>> selected_polars;  // polar atoms
                                                              // acting as
                                                              // acceptor or as
                                                              // donor through
                                                              // attached H
      for (core::Size sc_acc : rotamer->accpt_pos_sc()) {
        selected_polars.insert(std::make_pair(AA, sc_acc));
      }
      for (core::Size polar_H : rotamer->Hpos_polar_sc()) {
        selected_polars.insert(std::make_pair(AA, rotamer->atom_base(polar_H)));
      }

      for (std::pair<char, core::Size> polar_atomid : selected_polars) {
        out << "- x" << mres_AA2var_idx[mresid - 1][AA];
        for (core::Size rot1 : res_pola2hbonding_S1[mresid - 1][polar_atomid]) {
          out << " + x" << rot1;
        }
        for (core::Size var_idx :
             res_pola2hbonding_S2C[mresid - 1][polar_atomid]) {
          out << " + x" << var_idx;
        }
        out << " = 0" << std::endl;
      }
    }
  }
  out << "Binary" << std::endl;
  for (core::Size var = 1; var <= final_offset; var++) {
    out << "x" << var << ((var % 12) ? " " : "\n");
  }
  out << std::endl << "End" << std::endl;
  out.close();
}

void HBNet::print_CFN_model_to_file(std::string filename, bool prefpolar,
                                    bool pref2polar) {
  std::ofstream out;
  out.open(filename);

  // Second version.

  // The "all polar atoms must be satisfied constraint" can be
  // translated into an automata with manually decomposed ternary
  // regular like constraints that take care of each atom
  // (independently for now. TODO). The automata will have one
  // accepting state that is reached as soon as the atom is satisfied
  // and one non accepting state (one automata per rotamer). When a
  // given rotamer is used, if a suitable rotamer is used in the first
  // position of the regular, we accept else we shift to the non
  // accepting state of the rotamer and next position.  Problem: there
  // will up to n²d boolean state variables (but probably far less).

  // Another possibility is to use the clause constraint in tb2.
  // Directly, this requires to channel the domains to nd booleans
  // (and constraints). Binary constraints (AA level) are encoded as
  // binary CFs and the channel is left unused. When channeling, there
  // is lower propagation in case of multiple values remaining in a
  // single variable (binary case Ok). The clause constraint is built
  // automatically from table constraints with non zero cost on one
  // tuple.

  // Another one is to define a new anti-Horn (or direct encoding)
  // clause constraint that will work directly on non necessarily
  // boolean variables. The code needs to be defined and a new syntax
  // too.

  // Yet another would be to better handle non binary CF. With the
  // huge domains, having a ternary with a table would be impractical
  // at this point.

  // Trying option 2 for now.

  const core::Size top =
      1000000000;  // TODO: adjust to something less lousy (as in PWMaxSAT).
  core::Size polar2_penalty = 0;
  const core::Size polar1_penalty =
      (pref2polar ? rotamer_sets_->nmoltenres() : 0);
  const core::Size apolar_penalty =
      (prefpolar ? polar1_penalty + 1 : polar1_penalty);

  // Extracting hbond info.
  const core::scoring::hbonds::HBondDatabaseCOP hb_database =
      core::scoring::hbonds::HBondDatabase::get_database();
  const core::scoring::TenANeighborGraph &tenA_neighbor_graph(
      orig_pose_->energies().tenA_neighbor_graph());

  // Initialize Hbond data structure
  std::vector<std::map<core::Size, std::set<core::Size>>> rot_atom2rotidset;
  {
    // Fill the vector with maps
    for (core::Size rotid = 1; rotid <= rotamer_sets_->nrotamers(); rotid++) {
      std::map<core::Size, std::set<core::Size>> init_map;
      rot_atom2rotidset.push_back(init_map);
    }

    // Fill the maps with sets
    for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres();
         ++mresid) {
      // For each global rotamer and each polar atom of the AA, it will
      // hold a set of rotids that satisfy the atom (rotid 0 is the
      // backbone).
      auto rot_set_for_mres =
          rotamer_sets_->rotamer_set_for_moltenresidue(mresid);

      for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
           ++rotid) {
        core::conformation::ResidueCOP rotamer =
            rot_set_for_mres->rotamer(rotid);
        core::Size global_rotamer =
            rotamer_sets_->nrotamer_offset_for_moltenres(mresid) + rotid;

        // For each rotamer, we have a set of acceptor and donor atoms
        for (core::Size sc_acc : rotamer->accpt_pos_sc()) {
          std::set<core::Size> myemptyset;
          rot_atom2rotidset[global_rotamer - 1][sc_acc] = myemptyset;
        }

        for (core::Size polar_H : rotamer->Hpos_polar_sc()) {
          core::Size H_parent = rotamer->atom_base(polar_H);
          std::set<core::Size> myemptyset;
          rot_atom2rotidset[global_rotamer - 1][H_parent] = myemptyset;
        }
      }
    }

    for (std::list<EdgeBase *>::const_iterator iter =
             ig_->get_edge_list_begin();
         iter != ig_->get_edge_list_end(); ++iter) {
      PDEdge *pdedge = static_cast<PDEdge *>(*iter);

      core::Size const mres_1 = (*iter)->get_first_node_ind();
      core::Size const mres_2 = (*iter)->get_second_node_ind();

      core::Size const resid_1 = rotamer_sets_->moltenres_2_resid(mres_1);
      core::Size const resid_2 = rotamer_sets_->moltenres_2_resid(mres_2);

      const core::Size si = rotamer_sets_->nrotamers_for_moltenres(mres_1);
      const core::Size sj = rotamer_sets_->nrotamers_for_moltenres(mres_2);

      const core::Size nbrs_1 = tenA_neighbor_graph.get_node(resid_1)
                                    ->num_neighbors_counting_self_static();
      const core::Size nbrs_2 = tenA_neighbor_graph.get_node(resid_2)
                                    ->num_neighbors_counting_self_static();

      for (core::Size ii = 1; ii <= si; ++ii) {
        for (core::Size jj = 1; jj <= sj; ++jj) {
          core::Real const score = pdedge->get_two_body_energy(ii, jj);
          core::Size const global_rotamer_ii =
              rotamer_sets_->nrotamer_offset_for_moltenres(mres_1) + ii;
          core::Size const global_rotamer_jj =
              rotamer_sets_->nrotamer_offset_for_moltenres(mres_2) + jj;

          if (score < 0) {
            core::conformation::ResidueCOP rotamer_ii =
                rotamer_sets_->rotamer(global_rotamer_ii);
            core::conformation::ResidueCOP rotamer_jj =
                rotamer_sets_->rotamer(global_rotamer_jj);

            core::scoring::hbonds::HBondSet hbond_set;

            core::scoring::hbonds::identify_hbonds_1way(
                *hb_database, *rotamer_ii, *rotamer_jj, nbrs_1, nbrs_2,
                false,  // bool const evaluate_derivative,
                false,  // bool const exclude_don_bb,
                false,  // bool const exclude_don_bsc,
                false,  // bool const exclude_acc_scb,
                false,  // bool const exclude_acc_sc,
                // output
                hbond_set);

            core::scoring::hbonds::identify_hbonds_1way(
                *hb_database, *rotamer_jj, *rotamer_ii, nbrs_2, nbrs_1,
                false,  // bool const evaluate_derivative,
                false,  // bool const exclude_don_bb,
                false,  // bool const exclude_don_bsc,
                false,  // bool const exclude_acc_scb,
                false,  // bool const exclude_acc_sc,
                // output
                hbond_set);

            for (core::Size hbond_id = 1; hbond_id <= hbond_set.nhbonds();
                 ++hbond_id) {
              core::scoring::hbonds::HBondCOP hbond =
                  hbond_set.hbond_cop(hbond_id);
              bool is_don_bb = hbond->don_hatm_is_backbone();
              bool is_acc_bb = hbond->acc_atm_is_backbone();

              if ((is_acc_bb && is_don_bb) ||
                  (hbond->energy() > hydrogen_bond_threshold_)) {
                continue;  // not interested in weak or BB to BB hbonds
              }

              core::Size const H = hbond->don_hatm();
              core::Size const acc = hbond->acc_atm();
              core::Size H_parent;

              core::Size don_rotamer;
              core::Size acc_rotamer;

              if (hbond->don_res() == resid_1) {
                don_rotamer = global_rotamer_ii;
                H_parent = rotamer_ii->atom_base(hbond->don_hatm());
                acc_rotamer = global_rotamer_jj;
              } else {
                don_rotamer = global_rotamer_jj;
                H_parent = rotamer_jj->atom_base(hbond->don_hatm());
                acc_rotamer = global_rotamer_ii;
              }

              char don_AA = rotamer_sets_->rotamer(don_rotamer)->name1();
              char acc_AA = rotamer_sets_->rotamer(acc_rotamer)->name1();

              // !!!!!!!!!!!!!!!!!! sheffler !!!!!!!!!!!!!!!!!!!!!!!!
              // this is how to mark atoms as satisfied
              // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if (is_don_bb) {  // We have a self sat acceptor atom for the
                                // rotamer here
                rot_atom2rotidset[acc_rotamer - 1][acc].clear();
                rot_atom2rotidset[acc_rotamer - 1][acc].insert(0);
              }
              if (is_acc_bb) {  // We have a self sat donor atom for the rotamer
                                // here
                rot_atom2rotidset[don_rotamer - 1][H_parent].clear();
                rot_atom2rotidset[don_rotamer - 1][H_parent].insert(0);
              }
              if (!is_acc_bb && !is_don_bb) {  // the polar partners both come
                                               // from side-chains
                if (rot_atom2rotidset[acc_rotamer - 1][acc].count(0) == 0)
                  rot_atom2rotidset[acc_rotamer - 1][acc].insert(don_rotamer);
                if (rot_atom2rotidset[don_rotamer - 1][H_parent].count(0) == 0)
                  rot_atom2rotidset[don_rotamer - 1][H_parent].insert(
                      acc_rotamer);
              }
            }
          }
        }
      }
    }
  }

  // UP Stage: enforce unit propagation on the anti-horn clauses. This
  // will lower the number of non unit/binary clauses and channeling
  // clauses and vars.  For the moment, this is done naively, using an
  // AC1 like loop even if linear time.
  std::map<core::Size, bool> deletedrot;
  {
    // look for unit anti-horn clauses (rotamers to delete)
    for (core::Size rotid = 1; rotid <= rotamer_sets_->nrotamers(); rotid++) {
      for (const auto &a_set_pair : rot_atom2rotidset[rotid - 1]) {
        if (a_set_pair.second.empty()) {
          deletedrot[rotid] = true;
          break;
        }
      }
    }

    bool new_unit;
    do {
      new_unit = false;
      for (core::Size rotid = 1; rotid <= rotamer_sets_->nrotamers(); rotid++) {
        for (auto &a_set_pair : rot_atom2rotidset[rotid - 1]) {
          if (!a_set_pair.second.empty() &&
              a_set_pair.second.count(0) == 0) {  // no BB support, not empty
            for (auto it_rot = a_set_pair.second.begin();
                 it_rot != a_set_pair.second.end();) {
              if (deletedrot[*it_rot]) {
                it_rot = a_set_pair.second.erase(it_rot);
              } else {
                ++it_rot;
              }
            }
            if (a_set_pair.second.empty()) {
              new_unit = true;
              deletedrot[rotid] = true;
            }
          }
        }
      }
    } while (new_unit);
  }

  // Generate the CFN
  // Header - assumes one CF per IG-edge. Suboptimal as some may be
  // empty but toulbar2 will remove Additional variables needed for
  // channeling. Only rotamers with no empty hbond set or with not all
  // hbond to BB need one.  Additional CF used for deleting,
  // channeling + clauses
  core::Size nchanvar = 0;
  core::Size nchancf = 0;
  core::Size ndeletions = 0;
  core::Size nclauses = 0;
  core::Size chan_var_idx = rotamer_sets_->nmoltenres();
  std::map<core::Size, core::Size> rot2chanvaridx;

  // We need to channel only the rotamer that appear in clauses and
  // are not deleted
  for (core::Size rotid = 1; rotid <= rotamer_sets_->nrotamers(); rotid++) {
    bool needs_chan = false;

    if (deletedrot[rotid]) {
      ndeletions++;
    } else {
      for (const auto &a_set_pair : rot_atom2rotidset[rotid - 1]) {
        if (a_set_pair.second.count(0) == 0) {
          nclauses++;  // a clause
        }
      }
      rot2chanvaridx[rotid] = chan_var_idx++;
      nchancf++;
      nchanvar++;
    }
  }

  out << "Max_HB_stable"
      << " " << rotamer_sets_->nmoltenres() + nchanvar << " ";
  core::Size max_domain_size = 2;
  for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres(); ++mresid) {
    max_domain_size = std::max(max_domain_size,
                               rotamer_sets_->nrotamers_for_moltenres(mresid));
  }
  // we need one channeling constraint per channeled variable + deletion of
  // rotamers with no support.
  out << max_domain_size << " "
      << ig_->get_num_edges() + ig_->get_num_nodes() + nchancf + nclauses +
             ndeletions
      << " " << top << std::endl;

  // Domain sizes
  for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres(); ++mresid) {
    out << rotamer_sets_->nrotamers_for_moltenres(mresid) << " ";
  }
  // Channeling variables domains
  for (core::Size vidx = 1; vidx <= nchanvar; vidx++) {
    out << "2 ";
  }
  out << std::endl;

  // Unary cost functions that favor heavy polar atoms
  for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres(); ++mresid) {
    auto rot_set_for_mres =
        rotamer_sets_->rotamer_set_for_moltenresidue(mresid);
    out << "1 " << mresid - 1 << " 0 " << rot_set_for_mres->num_rotamers()
        << std::endl;

    // could be optimized by counting just once par amino-acid type
    for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
         ++rotid) {
      core::conformation::ResidueCOP rotamer = rot_set_for_mres->rotamer(rotid);
      bool is_polar = ((rotamer->Hpos_polar_sc().size() != 0) ||
                       (rotamer->accpt_pos_sc().size() != 0));
      bool is_2polar = (rotamer->accpt_pos_sc().size() > 1);

      out << rotid - 1 << " ";
      if (!is_polar)
        out << apolar_penalty;
      else if (!is_2polar)
        out << polar1_penalty;
      else
        out << polar2_penalty;
      out << std::endl;
    }
  }

  // Channeling constraints.
  for (core::Size rotid = 1; rotid <= rotamer_sets_->nrotamers(); rotid++) {
    core::Size mresid = rotamer_sets_->moltenres_for_rotamer(rotid);
    core::Size l_rot = rotamer_sets_->rotid_on_moltenresidue(rotid);

    if (deletedrot[rotid]) {
      out << "1 " << mresid - 1 << " 0 1" << std::endl;
      out << l_rot - 1 << " " << top << std::endl;
    } else {
      out << "2 " << rot2chanvaridx[rotid] << " " << mresid - 1 << " " << top
          << " ";
      out << rotamer_sets_->nrotamers_for_moltenres(mresid) << std::endl;
      out << "1 " << l_rot - 1 << " 0" << std::endl;
      for (core::Size false_rot = 1;
           false_rot <= rotamer_sets_->nrotamers_for_moltenres(mresid);
           false_rot++) {
        if (false_rot != l_rot) {
          out << "0 " << false_rot - 1 << " 0" << std::endl;
        }
      }
    }
  }

  // Binary cost functions (clashes)
  for (std::list<EdgeBase *>::const_iterator iter = ig_->get_edge_list_begin();
       iter != ig_->get_edge_list_end(); ++iter) {
    PDEdge *pdedge = static_cast<PDEdge *>(*iter);

    core::Size const mres_1 = (*iter)->get_first_node_ind();
    core::Size const mres_2 = (*iter)->get_second_node_ind();

    core::Size const resid_1 = rotamer_sets_->moltenres_2_resid(mres_1);
    core::Size const resid_2 = rotamer_sets_->moltenres_2_resid(mres_2);

    const core::Size si = rotamer_sets_->nrotamers_for_moltenres(mres_1);
    const core::Size sj = rotamer_sets_->nrotamers_for_moltenres(mres_2);

    // Filling in hbond creation info
    const core::Size nbrs_1 = tenA_neighbor_graph.get_node(resid_1)
                                  ->num_neighbors_counting_self_static();
    const core::Size nbrs_2 = tenA_neighbor_graph.get_node(resid_2)
                                  ->num_neighbors_counting_self_static();

    std::vector<std::vector<int>> pairs;
    core::Size nb_non_zero = 0;

    for (core::Size ii = 1; ii <= si; ++ii) {
      std::vector<int> row;
      for (core::Size jj = 1; jj <= sj; ++jj) {
        core::Real const score = pdedge->get_two_body_energy(ii, jj);
        core::Size const global_rotamer_ii =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_1) + ii;
        core::Size const global_rotamer_jj =
            rotamer_sets_->nrotamer_offset_for_moltenres(mres_2) + jj;
        if (score > clash_threshold_) {
          row.push_back(top);
          nb_non_zero++;
        } else {
          row.push_back(0);
        }
      }
      pairs.push_back(row);
    }

    const core::Size default_cost = 0;
    out << "2 " << mres_1 - 1 << " " << mres_2 - 1 << " " << default_cost << " "
        << nb_non_zero << std::endl;
    for (core::Size i = 0; i < pairs.size(); i++) {
      for (core::Size j = 0; j < pairs[i].size(); j++) {
        core::Size cost = pairs[i][j];
        if (cost != default_cost) {
          out << i << " " << j << " " << cost << std::endl;
        }
      }
    }
  }

  // Now, the clauses. They are described by a single tuple over the boolean
  // channeling variables.
  // Binary clauses are represented as binary cost functions.
  std::map<core::Size, core::Size> histo;
  for (core::Size rotid = 1; rotid <= rotamer_sets_->nrotamers(); rotid++) {
    core::Size mresid = rotamer_sets_->moltenres_for_rotamer(rotid);
    core::Size l_rot = rotamer_sets_->rotid_on_moltenresidue(rotid);

    if (!deletedrot[rotid]) {
      for (const auto &a_set_pair : rot_atom2rotidset[rotid - 1]) {
        if (a_set_pair.second.count(0) == 0) {
          std::set<core::Size> vars;
          for (const auto &rot_hb : a_set_pair.second) {
            vars.insert(rotamer_sets_->moltenres_for_rotamer(rot_hb));
          }
          histo[vars.size() + 1]++;
          if (vars.size() == 1) {  // binary constraint
            core::Size other_res = *vars.begin();
            core::Size nrot = rotamer_sets_->nrotamers_for_moltenres(other_res);
            out << "2 " << mresid - 1 << " " << other_res - 1 << " 0 "
                << nrot - a_set_pair.second.size() << std::endl;
            for (core::Size tup = 0; tup < nrot; tup++) {
              core::Size gtup =
                  rotamer_sets_->nrotamer_offset_for_moltenres(other_res) + 1 +
                  tup;
              if (a_set_pair.second.count(gtup) == 0)
                out << l_rot - 1 << " " << tup << " " << top << std::endl;
            }
          } else {
            // Arity, scope, def cost, tuples
            out << a_set_pair.second.size() + 1 << " " << rot2chanvaridx[rotid]
                << " ";
            for (const auto &rot_hb : a_set_pair.second) {
              if (rot2chanvaridx.find(rot_hb) != rot2chanvaridx.end()) {
                out << rot2chanvaridx[rot_hb] << " ";
              } else {
                TR << "----------- ERROR: no channeling var for " << rot_hb
                   << "deletion status " << deletedrot[rot_hb] << std::endl;
              }
              vars.insert(rotamer_sets_->res_for_rotamer(rot_hb));
            }
            out << "0 1" << std::endl;

            // The forbidden tuple (negation of the litterals in the clause)
            out << "1";
            for (const auto &rot_hb : a_set_pair.second) {
              out << " 0";
            }
            out << " " << top << std::endl;
          }
        }
      }
    }
  }
  for (const auto &fp : histo)
    TR << "Arity " << fp.first << " occurrences " << fp.second << std::endl;
  // Write the corresponding amino acids for printing sequences in
  // toulbar2 directly For true rotamer variables, we print he AA
  // name. For the channeling variables, let's use 'A'
  for (core::Size mresid = 1; mresid <= rotamer_sets_->nmoltenres(); ++mresid) {
    auto rot_set_for_mres =
        rotamer_sets_->rotamer_set_for_moltenresidue(mresid);
    for (core::Size rotid = 1; rotid <= rot_set_for_mres->num_rotamers();
         ++rotid) {
      out << rot_set_for_mres->rotamer(rotid)->name1() << " ";
    }
    out << std::endl;
  }

  for (core::Size cv = 1; cv <= nchanvar; ++cv) {
    out << "* *" << std::endl;
  }
  out.close();
}

/*void HBNet::print_interaction_graph_to_binary_file( std::string filename ){
  unsigned long num_bytes = 0;//keep count of how large we expect the file to be

        std::ofstream out( filename, std::ios::binary );

  int const nmres = int( rotamer_sets_->nmoltenres() );
  out.write( reinterpret_cast <const char*> ( &nmres ), sizeof( int ) );
  num_bytes += sizeof( int );

  for( core::Size ii = 1; ii <= rotamer_sets_->nmoltenres(); ++ii ){
    int const nrotamers_for_moltenres_ii = int (
  rotamer_sets_->nrotamers_for_moltenres( ii ) );
    out.write( reinterpret_cast <const char*> ( &nrotamers_for_moltenres_ii ),
  sizeof( int ) );
    num_bytes += sizeof( int );

    int const chain = int( orig_pose_->chain( rotamer_sets_->moltenres_2_resid(
  ii ) ) );
    out.write( reinterpret_cast <const char*> ( &chain ), sizeof( int ) );
    num_bytes += sizeof( int );
  }

  for( std::list< EdgeBase* >::const_iterator iter = ig_->get_edge_list_begin();
       iter != ig_->get_edge_list_end(); ++iter ){
    int const i = int( (*iter)->get_first_node_ind() );
    int const j = int( (*iter)->get_second_node_ind() );

    out.write( reinterpret_cast <const char*> ( &i ), sizeof( int ) );
    out.write( reinterpret_cast <const char*> ( &j ), sizeof( int ) );
    num_bytes += 2*sizeof( int );

    PDEdge * pdedge =  static_cast< PDEdge * >( *iter );

    const core::Size si = rotamer_sets_->nrotamers_for_moltenres(
  (*iter)->get_first_node_ind() );
    const core::Size sj = rotamer_sets_->nrotamers_for_moltenres(
  (*iter)->get_second_node_ind() );
    for( core::Size ii = 1; ii <= si; ++ii ){
      for( core::Size jj = 1; jj <= sj; ++jj ){
                                float energy = pdedge->get_two_body_energy(ii,
  jj);
                                out.write( reinterpret_cast <char*> ( &energy ),
  sizeof( float ) );
      }
    }
    num_bytes += si*sj*sizeof( float );
    out << std::endl;
  }

  out.close();
  TR << "sizeof( int ): " << sizeof( int ) << std::endl;
  TR << "sizeof( float ): " << sizeof( float ) << std::endl;
  TR << "predicted file size: " << num_bytes << " bytes" << std::endl;
        }*/

}  // hbnet
}  // protocols

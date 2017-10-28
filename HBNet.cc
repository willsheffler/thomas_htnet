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
#include <basic/options/option_macros.hh>

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

#include <core/scoring/lkball/LK_BallInfo.hh>

// #include <boost/container/flat_set.hpp> // not in rosetta!!!
#include <unordered_set>

using core::scoring::lkball::LKB_ResidueInfo;

#include <numeric/xyz.io.hh>

using namespace core;
using namespace pose;
using namespace pack;
using namespace rotamer_set;
using namespace interaction_graph;
using namespace conformation;
using namespace scoring::hbonds;

// whs: gross, but no grosser than the options system in general
OPT_1GRP_KEY(Integer, hbsat, num_nbr_solv_cut)
static int register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  NEW_OPT(hbsat::num_nbr_solv_cut,
          "max num ca/cb wtihin 8A for polar atom to be solvated", 9);
  return 13;
}
static int dummy = register_options();
// end whs

namespace protocols {
namespace hbnet {

// whs: started to do something more complicated...
// decided against it.
// #include <scheme/objective/voxel/VoxelArray.hh>
// using namespace scheme::objective::voxel;

// VoxelArray<3, float> make_solvation_estimate(Pose const &pose) {
//   numeric::xyzVector<float> lb = pose.residue(1).xyz(1), ub;
//   ub = lb;
//   for (size_t ir = 1; ir <= pose.size(); ++ir) {
//     for (size_t ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia) {
//       lb.min(pose.residue(ir).xyz(ia));
//       ub.max(pose.residue(ir).xyz(ia));
//     }
//   }
//   VoxelArray<3, float> solv(lb, ub, numeric::xyzVector<float>(1, 1, 1));
//   std::cout << lb << std::endl;
//   std::cout << ub << std::endl;
//   std::cout << solv.size() << std::endl;

//   return solv;
// }

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
}

///@details recursive call for depth-first traversal of IG
void HBNet::recursive_traverse(int const new_node_ind, int const newstate,
                               Size const newres, Size const prevres,
                               utility::vector1<HBondResStructCOP> residues,
                               Size network_rec_count, Real init_sc,
                               Real const hb_threshold,
                               bool const second_search /* false */) {

}  // rec_trav

void HBNet::traverse_IG_for_3mers_find_third_rotamer(
    core::Size resA, core::Size moltenresA, core::Size rotA,
    core::Size resA_ind, core::Size resB, core::Size moltenresB,
    core::Size rotB, core::Size resB_ind, core::Real const hb_threshold,
    core::Real init_sc, utility::vector1<HBondResStructCOP> &residues) {}

bool compare_hbond_candidates(const hbond_candidate &a,
                              const hbond_candidate &b) {
  return a.score < b.score;
}

void HBNet::greedy_traverse_IG(Real const hb_threshold) {}

bool HBNet::greedy_recursive_traverse(
    core::Real const hb_threshold, Size moltenres,
    utility::vector1<core::Size> &state_for_moltenres,
    core::pose::PoseOP &network_pose) {
  return false;
}

///@details traverse EnergyGraph of static pose to find all native networks
void HBNet::traverse_native(Pose const &pose, Real const hb_threshold) {}

void HBNet::rec_trav_native(Pose const &pose, Size new_res, Size prev_res,
                            utility::vector1<HBondResStructCOP> residues,
                            Real const hb_threshold) {}

bool HBNet::no_clash(Size moltenres1, Size state1, Size moltenres2,
                     Size state2) {}

///@details Check if a new rotamer state clashes with any rotamer states already
/// in a given h-bond network
///    return of true = it clashes, false = no clashes
bool HBNet::check_clash(utility::vector1<HBondResStructCOP> const &residues,
                        platform::uint new_node_ind, Size newstate, Size newres,
                        Real &init_score, bool &cycle) {}  // check_clash

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
  return false;
}  // net_clash

bool HBNet::network_already_stored(
    utility::vector1<HBondResStructCOP> &residues,
    utility::vector1<HBondResStructCOP> &i_residues) {}

void HBNet::store_network(utility::vector1<HBondResStructCOP> residues,
                          Real init_score, bool term_w_start, bool term_w_cycle,
                          bool score_now, bool native) {}  // store_network

void HBNet::minimize_network(Pose &pose, hbond_net_struct &network,
                             bool residues_already_placed /* true */) {}

void HBNet::score_network_on_pose(Pose &pose, hbond_net_struct &i) {}

void HBNet::merge_2_branched_networks(
    utility::vector1<HBondResStructCOP> const &residues1,
    utility::vector1<HBondResStructCOP> const &residues2,
    utility::vector1<HBondResStructCOP> &new_residues) {}

void HBNet::merge_2_branched_networks(hbond_net_struct const &i,
                                      hbond_net_struct const &j,
                                      HBondNetStructOP new_network) {}

void HBNet::merge_branched_networks(core::Size merged_vec_index,
                                    HBondNetStructOP new_network) {}

// consider Ser and Thr idential for benchmarking purposes
bool HBNet::networks_identical_aa_sequence(hbond_net_struct const &i,
                                           hbond_net_struct const &j) {}

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
                            bool &branch, bool true_if_identical /* true */) {}

///@details Critical function to HBNet; because of depth-frist traversal, need
/// to combine branched cases, e.g. ASN that can make multiple h-bonds
///   The reason we need this is because bread-first additions ("branches") are
///   order-dependent, meaning we need to try ALL combinations
void HBNet::branch_overlapping_networks() {}  // branch_overlapping networks

void HBNet::finalize_branching() {}

inline unsigned long combine_indices(core::Size const i, core::Size const j) {
  return (static_cast<unsigned long>(i) << 32) & (j & 0xFFFFFFFF);
}

inline std::pair<core::Size, core::Size> split_indices(
    unsigned long const combined) {
  core::Size const i = (combined >> 32) & 0xFFFFFFFF;
  core::Size const j = combined & 0xFFFFFFFF;
  return std::make_pair(i, j);
}

void HBNet::branch_overlapping_networks_monte_carlo() {}

void HBNet::monte_carlo_branching_trajectory(
    core::Size seed_index, std::set<unsigned long> &clashes,
    std::vector<utility::vector1<core::Size>> &overlapping_nets) {}

core::Real HBNet::get_twobody(HBondResStructCOP i, HBondResStructCOP j) {}

// used by branch_overlapping() to efficiently search for all combinations of
// compatible networks that can be merged
void HBNet::rec_set_intersection(std::vector<Size> add_index_vec,
                                 std::vector<Size> next_index_vec, Size pos) {}

void HBNet::place_rots_on_pose(pose::Pose &pose, hbond_net_struct &i,
                               bool use_pose) {}

// to very quickly eliminate networks with 1 or more unsatisfied heavy-atom
// donors or acceptors
bool HBNet::quick_and_dirty_network_has_heavy_atom_unsat(
    Pose const &pose, hbond_net_struct const &network) {}

// to very quickly eliminate networks with 1 or more unsatisfied heavy-atom
// donors or acceptors
bool HBNet::quick_and_dirty_heavy_atom_is_unsat(Pose const &pose,
                                                id::AtomID const at_id) {}

// hacky check based on distance for cases where we only have water O (not
// explicit H's)
bool HBNet::atom_hbonds_to_bridging_water(Pose const &pose,
                                          id::AtomID const at_id) {}

void HBNet::update_core_and_boundary_residues(core::pose::Pose const &pose) {}

void HBNet::find_unsats(Pose const &pose, hbond_net_struct &network) {

}  // find_unsats

// checks resnum and aa but not rotamer
// consider Ser and Thr identical
bool HBNet::residues_identical(utility::vector1<HBondResStructCOP> &residues1,
                               utility::vector1<HBondResStructCOP> &residues2) {

}

bool HBNet::all_residue_chis_are_close(
    utility::vector1<HBondResStructCOP> &residues1,
    utility::vector1<HBondResStructCOP> &residues2) {}

bool HBNet::residues_not_unique(
    utility::vector1<HBondResStructCOP> &residues1,
    utility::vector1<HBondResStructCOP> &residues2) {}

bool HBNet::networks_unique(hbond_net_struct const &i,
                            hbond_net_struct const &j,
                            bool no_surface /* true */) {}

// 02/28/15 changing behavior, will only keep 1 network with unique resnum/aa
// identity; networks with same resnum/aa seq but unique rotamers
// considered through final scoring
void HBNet::remove_replicate_networks(Size same_max /*=1*/) {

}  // remove_replicate_networks

Size HBNet::get_num_native_rot(
    Pose &pose, utility::vector1<HBondResStructCOP> const &residues,
    Real sc_rmsd_cut, bool super) {}

Size HBNet::get_num_native_seq(
    core::pose::Pose &pose,
    utility::vector1<HBondResStructCOP> const &residues) {}

void HBNet::output_networks(bool finalize) {}  // output_networks

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

void HBNet::set_symmetry(Pose &pose) {}

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

void HBNet::write_network_pdb(
    HBondNetStructOP p)  // better to pass object pointer or reference?
{}

void HBNet::jack_store_network(HBondNetStructOP p, bool native) {}

// assumes the network rotamers are already placed onto the pose
void HBNet::set_constraints(
    Pose &pose, core::scoring::constraints::ConstraintSet &constraints,
    HBondNetStructOP p, bool write_cst_file /* false */) {}

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

bool HBNet::water_clashes(Pose const &pose, Vector const water_O) {}

bool HBNet::water_oxygen_clashes_with_residue(Vector const water_oxygen,
                                              Size const resnum,
                                              int const rot_state) {}

bool HBNet::water_oxygen_clashes_with_residue(Vector const water_oxygen,
                                              Residue const &res) {}

pose::PoseOP HBNet::get_additional_output() {}

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

void HBNet::select_best_networks() {}  // select_best_networks()

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
        if (task_->design_residue((int)r) || is_repack[r] == 1) {
          start_res_vec_.insert(r);
        }
      }
    }
  }
}

void HBNet::get_native_networks(Pose const &pose) {}

void HBNet::setup(Pose const &pose) {
  this->setup_packer_task_and_starting_residues(pose);

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
  { rotamer_sets_ = RotamerSetsOP(new rotamer_set::RotamerSets()); }

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

    rotamer_sets_->prepare_sets_for_packing(pose, *init_scorefxn_);
    ig_ = InteractionGraphFactory::create_interaction_graph(
        *task_, *rotamer_sets_, pose, *init_scorefxn_, *packer_neighbor_graph_);
    ig_->initialize(*rotamer_sets_);

    PrecomputedPairEnergiesInteractionGraphOP pig(
        utility::pointer::dynamic_pointer_cast<
            PrecomputedPairEnergiesInteractionGraph>(ig_));
    runtime_assert(pig);

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

int calc_num_nbrs(numeric::xyzVector<double> solv_xyz,
                  std::vector<numeric::xyzVector<float>> const &nbr_coords) {
  int num_nbr = 0;
  for (auto xyz : nbr_coords)
    if (xyz.distance_squared(solv_xyz) < 8.0 * 8.0) ++num_nbr;
  return num_nbr;
}

// whs: loop over all sidechain polar heavy atoms, mark solvated
// define solvation as the most exposed LKBall group of the polar
// atom. where 'exposed' means number of CB/CA atoms within 8A is
// less than num_nbr_solv_cut
template <class RotamerWhatever>
void compute_solvation_info(
    PoseOP pose, RotamerWhatever const &rotamer_sets, int num_nbr_solv_cut,
    std::set<std::pair<core::Size, core::Size>> &atom_is_solvated,
    utility::vector1<bool> &rotamer_is_fully_solvated) {
  ;

  {
    std::cout << "begin solv calc..." << std::endl;
    std::vector<numeric::xyzVector<float>> nbr_coords;
    {
      for (size_t ir = 1; ir <= pose->size(); ++ir) {
        if (pose->residue(ir).has("CB"))
          nbr_coords.push_back(pose->residue(ir).xyz("CB"));
        else if (pose->residue(ir).has("CA"))
          nbr_coords.push_back(pose->residue(ir).xyz("CA"));
      }
    }
    int num_solv_rotamers = 0, max_nnbrs = 0;
    // for (core::Size irot = 1; irot <= rotamer_sets->nrotamers(); ++irot) {
    for (size_t ir = 1; ir <= rotamer_sets->nmoltenres(); ++ir) {
      auto res_rots = rotamer_sets->rotamer_set_for_moltenresidue(ir);
      for (size_t mirot = 1; mirot <= res_rots->num_rotamers(); ++mirot) {
        ResidueCOP rotamer = res_rots->rotamer(mirot);
        size_t irot = rotamer_sets->nrotamer_offset_for_moltenres(ir) + mirot;
        bool all_polar_solv = true;
        // oops! an hour of confustion becasue I eliminated the only
        // option for non-polar residues! no solutions indeed.
        if (rotamer->name3() == "GLY") all_polar_solv = false;
        auto waters = LKB_ResidueInfo(*rotamer).waters();
        for (size_t sc_acc : rotamer->accpt_pos_sc()) {
          int num_nbr = 999999999;
          for (numeric::xyzVector<double> lkxyz : waters[sc_acc]) {
            num_nbr = std::min(num_nbr, calc_num_nbrs(lkxyz, nbr_coords));
            // if (num_nbr == 0) {
            // std::cout << "WTF!!! " << lkxyz << std::endl;
            // core::pose::Pose pose;
            // pose.append_residue_by_jump(*rotamer, 1);
            // pose.dump_pdb("test.pdb");

            // std::exit(0);
            // }
          }
          // std::cout << num_nbr << std::endl;
          max_nnbrs = std::max(max_nnbrs, num_nbr);
          if (num_nbr <= num_nbr_solv_cut)
            atom_is_solvated.insert(std::make_pair(irot, sc_acc));
          else
            all_polar_solv = false;
        }
        for (core::Size polar_H : rotamer->Hpos_polar_sc()) {
          core::Size H_parent = rotamer->atom_base(polar_H);
          int num_nbr = 999999999;
          for (auto lkxyz : waters[H_parent]) {
            num_nbr = std::min(num_nbr, calc_num_nbrs(lkxyz, nbr_coords));
          }
          // std::cout << num_nbr << std::endl;
          max_nnbrs = std::max(max_nnbrs, num_nbr);
          if (num_nbr <= num_nbr_solv_cut)
            atom_is_solvated.insert(std::make_pair(irot, H_parent));
          else
            all_polar_solv = false;
        }
        rotamer_is_fully_solvated[irot] = all_polar_solv;
        if (all_polar_solv) ++num_solv_rotamers;
      }
    }
    std::cout << "max polar atom neighbors: " << max_nnbrs << std::endl;
    std::cout << "num fully solv rotamers:  " << num_solv_rotamers << " of "
              << rotamer_sets->nrotamers() << std::endl;
  }
  // end whs
}

template <class RotamerWhatever>
void print_rotsets_aas(RotamerWhatever const &rotsets) {
  std::cout << "AAs at positions (sanity check): ";
  for (size_t mir = 1; mir <= rotsets->nmoltenres(); ++mir) {
    std::set<char> aas;
    auto rotset = rotsets->rotamer_set_for_moltenresidue(mir);
    for (size_t irot = 1; irot <= rotset->num_rotamers(); ++irot) {
      aas.insert(rotset->rotamer(irot)->name1());
    }

    std::cout << mir << ":";
    for (char aa : aas) std::cout << aa;
    std::cout << " ";
  }
  std::cout << std::endl;
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

  const int num_nbr_solv_cut = basic::options::option
      [basic::options::OptionKeys::hbsat::num_nbr_solv_cut]();

  // whs: loop over all sidechain polar heavy atoms, mark solvated
  // define solvation as the most exposed LKBall group of the polar
  // atom. where 'exposed' means number of CB/CA atoms within 8A is
  // less than num_nbr_solv_cut
  std::set<std::pair<core::Size, core::Size>> atom_is_solvated;
  utility::vector1<bool> rotamer_is_fully_solvated(rotamer_sets_->nrotamers());
  compute_solvation_info(orig_pose_, rotamer_sets_, num_nbr_solv_cut,
                         atom_is_solvated, rotamer_is_fully_solvated);
  print_rotsets_aas(rotamer_sets_);
  // end whs

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
          // whs: skip if either is 'fully' solvated
          if (rotamer_is_fully_solvated[global_rotamer_ii]) continue;
          if (rotamer_is_fully_solvated[global_rotamer_jj]) continue;
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

              // whs: basically, just pretend the solvated groups are
              // satisfed by the backbone...
              bool don_sat =
                  is_acc_bb ||
                  atom_is_solvated.count(std::make_pair(don_rotamer, H_parent));
              bool acc_sat =
                  is_don_bb ||
                  atom_is_solvated.count(std::make_pair(acc_rotamer, acc));

              // We have a self sat acceptor atom for the rotamer here
              if (acc_sat) {
                rot_atom2rotidset[acc_rotamer - 1][acc].clear();
                rot_atom2rotidset[acc_rotamer - 1][acc].insert(0);
              }
              // We have a self sat donor atom for the rotamer here
              if (don_sat) {
                rot_atom2rotidset[don_rotamer - 1][H_parent].clear();
                rot_atom2rotidset[don_rotamer - 1][H_parent].insert(0);
              }
              if (!don_sat && !acc_sat) {  // the polar partners both come
                                           // from side-chains
                                           // whs: and need to be satisfied
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

    // whs: also delete any "fully solvated" rotamers
    // TODO: what about polars these rotamers satisfy?!?!
    for (size_t irot = 1; irot <= rotamer_is_fully_solvated.size(); ++irot) {
      if (rotamer_is_fully_solvated[irot]) deletedrot[irot] = true;
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

    int ndel = 0;
    for (auto del : deletedrot) {
      if (del.second) ++ndel;
    }
    std::cout << "total deleted rots: " << ndel << " of "
              << rotamer_sets_->nrotamers() << std::endl;
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
  unsigned long num_bytes = 0;//keep count of how large we expect the file to
  be

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

    int const chain = int( orig_pose_->chain(
  rotamer_sets_->moltenres_2_resid(
  ii ) ) );
    out.write( reinterpret_cast <const char*> ( &chain ), sizeof( int ) );
    num_bytes += sizeof( int );
  }

  for( std::list< EdgeBase* >::const_iterator iter =
  ig_->get_edge_list_begin();
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
                                out.write( reinterpret_cast <char*> ( &energy
  ),
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

#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "biodynamo.h"

#include "tile_extraction.h"
#include "tiling_config.h"
// tiling_config.h already pulls these in; included explicitly because this
// translation unit provides the one-and-only definition of SimParam::kUid.
#include "sim-param.h"

using namespace bdm;
using namespace bdm::tiling;

// Single out-of-line definition required by BDM_PARAM_GROUP_HEADER. In the real
// demo this lives in epidemiology.cc; here tiling_experiment.cc is the only
// translation unit that references SimParam, so it owns the definition.
namespace bdm {
const ParamGroupUid SimParam::kUid = ParamGroupUidGenerator::Get()->NewUid();
}  // namespace bdm

namespace {

struct Args {
  std::optional<int64_t> k;
  std::optional<int64_t> ax;
  std::optional<int64_t> ay;
  std::optional<int64_t> az;
  std::optional<int64_t> n;
  int64_t seed_offset = 0;
  std::string output = "-";  // "-" => stdout
};

[[noreturn]] void Fail(const std::string& msg) {
  std::cerr << "tiling_experiment: error: " << msg << "\n"
            << "usage: tiling_experiment --k <int> --anchor-x <int> "
               "--anchor-y <int> --anchor-z <int> --n <int> "
               "[--seed-offset <int>] [--output <path>]\n";
  std::exit(1);
}

int64_t ParseInt(const std::string& flag, const std::string& value) {
  try {
    size_t pos = 0;
    const int64_t v = std::stoll(value, &pos);
    if (pos != value.size()) {
      Fail("invalid integer for " + flag + ": '" + value + "'");
    }
    return v;
  } catch (const std::exception&) {
    Fail("invalid integer for " + flag + ": '" + value + "'");
  }
}

Args ParseArgs(int argc, const char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    const std::string flag = argv[i];
    auto need_value = [&](void) -> std::string {
      if (i + 1 >= argc) {
        Fail("missing value for " + flag);
      }
      return argv[++i];
    };
    if (flag == "--k") {
      a.k = ParseInt(flag, need_value());
    } else if (flag == "--anchor-x") {
      a.ax = ParseInt(flag, need_value());
    } else if (flag == "--anchor-y") {
      a.ay = ParseInt(flag, need_value());
    } else if (flag == "--anchor-z") {
      a.az = ParseInt(flag, need_value());
    } else if (flag == "--n") {
      a.n = ParseInt(flag, need_value());
    } else if (flag == "--seed-offset") {
      a.seed_offset = ParseInt(flag, need_value());
    } else if (flag == "--output") {
      a.output = need_value();
    } else if (flag == "-h" || flag == "--help") {
      Fail("help requested");
    } else {
      Fail("unknown argument: '" + flag + "'");
    }
  }

  if (!a.k) Fail("--k is required");
  if (!a.ax) Fail("--anchor-x is required");
  if (!a.ay) Fail("--anchor-y is required");
  if (!a.az) Fail("--anchor-z is required");
  if (!a.n) Fail("--n is required");
  if (*a.k <= 0) Fail("--k must be > 0");
  if (*a.n < 0) Fail("--n must be >= 0");
  if (*a.ax < 0 || *a.ay < 0 || *a.az < 0) Fail("anchors must be >= 0");
  return a;
}

}  // namespace

int main(int argc, const char** argv) {
  // Register epidemiology-specific parameters so Get<SimParam>() succeeds in
  // both the source and tile simulations (behaviors read SimParam at runtime).
  // Must happen before any Simulation is constructed.
  Param::RegisterParamGroup(new SimParam());

  const Args args = ParseArgs(argc, argv);
  const uint64_t k = static_cast<uint64_t>(*args.k);
  const Tile tile{{static_cast<uint64_t>(*args.ax),
                   static_cast<uint64_t>(*args.ay),
                   static_cast<uint64_t>(*args.az)},
                  k};
  const uint64_t n = static_cast<uint64_t>(*args.n);
  const uint64_t seed = kBaseSeed + static_cast<uint64_t>(args.seed_offset);

  Simulation source_sim("tiling_source", [](Param* param) {
    ApplyDomainParams(param);
  });
  auto* env = BuildFullPopulation(&source_sim);

  if (source_sim.GetResourceManager()->GetNumAgents() != kNumAgents) {
    Fail("population build produced " +
         std::to_string(source_sim.GetResourceManager()->GetNumAgents()) +
         " agents, expected " + std::to_string(kNumAgents));
  }

  const auto g = QueryGridGeometry(env, kMinBound, kMaxBound);
  std::vector<Agent*> tile_agents =
      ExtractTile(source_sim.GetResourceManager(), env, tile);

  if (tile_agents.empty()) {
    std::cerr << "tiling_experiment: error: tile is empty (0 agents) for "
                 "anchor ("
              << tile.lo[0] << "," << tile.lo[1] << "," << tile.lo[2]
              << ") k=" << k << ". Populated box range per axis is x["
              << g.populated_lo[0] << "," << g.populated_hi[0] << ") y["
              << g.populated_lo[1] << "," << g.populated_hi[1] << ") z["
              << g.populated_lo[2] << "," << g.populated_hi[2]
              << "); the tile [anchor, anchor+k) likely falls outside it.\n";
    return 1;
  }

  Simulation tile_sim("tiling_experiment", [seed](Param* param) {
    ApplyDomainParams(param);
    param->random_seed = seed;  // distinct per tile via --seed-offset
  });
  tile_sim.Activate();

  auto* tile_rm = tile_sim.GetResourceManager();
  for (auto* agent : tile_agents) {
    // NewCopy() allocates via the ACTIVE simulation's memory manager (tile_sim)
    // because Agent::operator new routes through Simulation::GetActive().
    Agent* copy = agent->NewCopy();
    // The copy inherits the source uid; reassign so it is unique within
    // tile_sim and cannot collide with agents the simulation may create itself.
    copy->AssignNewUid();
    tile_rm->AddAgent(copy);
  }
  tile_rm->RebuildAgentUidMap();

  // Match the full run's spatial resolution (box length pinned to the infection
  // radius) so the tile's grid/neighbor layout mirrors the untiled run rather
  // than being auto-derived from the (smaller) agent diameter.
  auto* tile_env =
      dynamic_cast<UniformGridEnvironment*>(tile_sim.GetEnvironment());
  if (tile_env != nullptr) {
    tile_env->SetBoxLength(static_cast<int32_t>(kInfectionRadius));
  }

  // Match the real demo's scheduler config: epidemiology.h explicitly unschedules
  // the "mechanical forces" operation so agents interact only via the behaviors.
  // GetOps(...)[0] mirrors the demo exactly; it assumes the op is present (it is,
  // in the default op set) and would be out-of-bounds if it were already absent.
  auto* scheduler = tile_sim.GetScheduler();
  scheduler->UnscheduleOp(scheduler->GetOps("mechanical forces")[0]);

  const uint64_t tile_agent_count = tile_rm->GetNumAgents();

  // --- Run: time ONLY the Simulate(n) call. ---
  const auto t0 = std::chrono::steady_clock::now();
  tile_sim.GetScheduler()->Simulate(n);
  const auto t1 = std::chrono::steady_clock::now();
  const double wall_seconds =
      std::chrono::duration<double>(t1 - t0).count();

  // --- Report (key=value, one per line; trivially parseable in bash). ---
  auto write = [&](std::ostream& os) {
    os << "k=" << k << "\n"
       << "anchor_x=" << tile.lo[0] << "\n"
       << "anchor_y=" << tile.lo[1] << "\n"
       << "anchor_z=" << tile.lo[2] << "\n"
       << "n=" << n << "\n"
       << "seed_offset=" << args.seed_offset << "\n"
       << "seed=" << seed << "\n"
       << "agent_count=" << tile_agent_count << "\n"
       << "box_length=" << g.box_length << "\n"
       << "simulate_wall_seconds=" << wall_seconds << "\n"
       << "simulate_wall_ms=" << (wall_seconds * 1e3) << "\n";
  };

  if (args.output == "-") {
    write(std::cout);
  } else {
    std::ofstream ofs(args.output);
    if (!ofs) {
      std::cerr << "tiling_experiment: error: cannot open output file '"
                << args.output << "'\n";
      return 1;
    }
    write(ofs);
  }

  return 0;
}

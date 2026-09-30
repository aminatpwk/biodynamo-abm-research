#ifndef TILING_CONFIG_H_
#define TILING_CONFIG_H_

// Shared population/config for the tiling experiment.
//
// Phase 3b: the population is now Person agents with the full epidemiology
// behavior set (Infection, Recovery, RandomMovement), created exactly as
// demo/epidemiology/src/epidemiology.h does, so that Simulate(n) on an extracted
// tile exercises the same dynamics as the real epidemiology use case rather than
// a static grid. The tiling/extraction mechanics in tile_extraction.h are
// unchanged and reused as-is.

#include <cstdint>

#include "biodynamo.h"
#include "core/environment/uniform_grid_environment.h"

#include "tile_extraction.h"

// Epidemiology demo headers (demo/epidemiology/src is added to the include path
// in CMakeLists.txt). behavior.h also pulls in person.h and sim-param.h.
#include "behavior.h"
#include "person.h"
#include "sim-param.h"

namespace bdm {
namespace tiling {

constexpr uint64_t kNumAgents = 250000;
constexpr real_t kMinBound = 0;
constexpr real_t kMaxBound = 500;
// Pinned to the epidemiology infection radius. SetBoxLength takes int32_t, so
// this truncates to 10 -- exactly as validated in Phase 2.
constexpr real_t kInfectionRadius = 10.5092197414493;
constexpr uint64_t kBaseSeed = 4357;  // BioDynaMo's default Param::random_seed.

// Initial infected fraction. measles.json uses 500000 / (99500000 + 500000) =
// 0.5% exactly; sim-param.h's own defaults give 10 / 2010 = 0.4975%. We adopt
// 0.5% and split the fixed 250k total so the sum is exact.
constexpr real_t kInfectedFraction = 0.005;
constexpr uint64_t kNumInfected = 1250;  // round(kNumAgents * kInfectedFraction)
constexpr uint64_t kNumSusceptible = kNumAgents - kNumInfected;  // 248750

// set_param body shared by the source and tile simulations. Bounds + torus
// boundary + unit time step match demo/epidemiology/src/epidemiology.h (which
// sets simulation_time_step = 1 and bound_space = kTorus). random_seed is set
// separately per simulation by the caller (see tiling_experiment.cc).
inline void ApplyDomainParams(Param* param) {
  param->min_bound = kMinBound;
  param->max_bound = kMaxBound;
  param->bound_space = Param::BoundSpaceMode::kTorus;
  param->simulation_time_step = 1;
}

// Builds the full epidemiology population in the (already active) simulation and
// returns its UniformGridEnvironment with the grid built and box indices valid.
// Replicates epidemiology.h: two CreateAgentsRandom passes (susceptible then
// infected), each Person carrying Infection + Recovery, and RandomMovement added
// with probability moving_agents_ratio. Requires SimParam to be registered
// beforehand (Param::RegisterParamGroup) so Get<SimParam>() succeeds.
inline UniformGridEnvironment* BuildFullPopulation(Simulation* sim) {
  auto* env = dynamic_cast<UniformGridEnvironment*>(sim->GetEnvironment());
  if (env == nullptr) {
    Log::Fatal("tiling", "Environment is not a UniformGridEnvironment");
  }

  auto* sparam = sim->GetParam()->Get<SimParam>();
  // NOTE: captured single-thread RNG, exactly as epidemiology.h does. It is
  // read inside CreateAgentsRandom's parallel region, which is a benign data
  // race here because moving_agents_ratio == 1.0 makes the branch deterministic.
  auto* random = sim->GetRandom();

  State state = State::kSusceptible;
  auto person_creator = [&](const Real3& position) {
    auto* person = new Person(position);
    person->SetDiameter(sparam->agent_diameter);
    person->state_ = state;
    person->AddBehavior(new Infection());
    person->AddBehavior(new Recovery());
    if (random->Uniform() < sparam->moving_agents_ratio) {
      person->AddBehavior(new RandomMovement());
    }
    return person;
  };

  state = State::kSusceptible;
  ModelInitializer::CreateAgentsRandom(kMinBound, kMaxBound, kNumSusceptible,
                                       person_creator);
  state = State::kInfected;
  ModelInitializer::CreateAgentsRandom(kMinBound, kMaxBound, kNumInfected,
                                       person_creator);

  env->SetBoxLength(static_cast<int32_t>(kInfectionRadius));
  CommitAgentsAndUpdateGrid(sim);
  return env;
}

}  // namespace tiling
}  // namespace bdm

#endif  // TILING_CONFIG_H_

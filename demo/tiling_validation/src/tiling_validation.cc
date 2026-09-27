#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "biodynamo.h"
#include "tile_extraction.h"

using namespace bdm;
using namespace bdm::tiling;

namespace {

constexpr uint64_t kNumAgents = 250000;
constexpr real_t kMinBound = 0;
constexpr real_t kMaxBound = 500;
constexpr real_t kInfectionRadius = 10.5092197414493;
constexpr real_t kAgentDiameter = 2.0;
constexpr size_t kMaxTilesPerK = 8;
constexpr uint64_t kSampleSeed = 42;
constexpr real_t kAssumedAgentsPerBox = 2.26;
const std::vector<uint64_t> kSweep = {8, 10, 13, 16, 21, 28};

// Counts agents inside the physical extent of box range [lo, hi) using
// positions only, independent of GetBoxIndex/GetBoxCoordinates.
uint64_t BruteForceCount(ResourceManager* rm, const GridGeometry& g,
                         const std::array<uint64_t, 3>& lo,
                         const std::array<uint64_t, 3>& hi) {
  std::array<real_t, 3> min, max;
  for (int i = 0; i < 3; i++) {
    min[i] = g.dimensions[2 * i] + static_cast<real_t>(lo[i]) * g.box_length;
    max[i] = g.dimensions[2 * i] + static_cast<real_t>(hi[i]) * g.box_length;
  }
  uint64_t n = 0;
  rm->ForEachAgent([&](Agent* agent) {
    const auto& p = agent->GetPosition();
    if (p[0] >= min[0] && p[0] < max[0] && p[1] >= min[1] && p[1] < max[1] &&
        p[2] >= min[2] && p[2] < max[2]) {
      n++;
    }
  });
  return n;
}

std::array<uint64_t, 3> TileHi(const Tile& t) {
  return {t.lo[0] + t.k, t.lo[1] + t.k, t.lo[2] + t.k};
}

}  // namespace

int main(int argc, const char** argv) {
  auto set_param = [](Param* param) {
    param->min_bound = kMinBound;
    param->max_bound = kMaxBound;
    param->bound_space = Param::BoundSpaceMode::kTorus;
  };
  Simulation sim(argc, argv, set_param);
  auto* rm = sim.GetResourceManager();
  auto* env = dynamic_cast<UniformGridEnvironment*>(sim.GetEnvironment());
  if (env == nullptr) {
    Log::Fatal("tiling_validation", "Environment is not a UniformGridEnvironment");
  }

  ModelInitializer::CreateAgentsRandom(
      kMinBound, kMaxBound, kNumAgents, [](const Real3& position) {
        auto* agent = new SphericalAgent(position);
        agent->SetDiameter(kAgentDiameter);
        return agent;
      });
  // Same call as demo/epidemiology; SetBoxLength takes int32_t, so the
  // infection radius is truncated exactly as it is in the real demo.
  env->SetBoxLength(static_cast<int32_t>(kInfectionRadius));

  std::printf("agents in ResourceManager before commit: %zu\n",
              rm->GetNumAgents());
  CommitAgentsAndUpdateGrid(&sim);
  std::printf("agents in ResourceManager after commit:  %zu\n",
              rm->GetNumAgents());
  if (rm->GetNumAgents() != kNumAgents) {
    Log::Fatal("tiling_validation", "Agent commit failed");
  }

  auto g = QueryGridGeometry(env, kMinBound, kMaxBound);
  const real_t length = g.box_length;
  const real_t side = kMaxBound - kMinBound;
  const real_t agents_per_box =
      kNumAgents / (side * side * side) * length * length * length;

  std::printf("\n== live grid geometry ==\n");
  std::printf("box length: requested %.6f, GetBoxLength() = %d\n",
              kInfectionRadius, g.box_length);
  std::printf("boxes/axis: %u %u %u\n", g.num_boxes_axis[0],
              g.num_boxes_axis[1], g.num_boxes_axis[2]);
  std::printf("dimensions: x[%d,%d) y[%d,%d) z[%d,%d)\n", g.dimensions[0],
              g.dimensions[1], g.dimensions[2], g.dimensions[3],
              g.dimensions[4], g.dimensions[5]);
  std::printf("populated box range: x[%lu,%lu) y[%lu,%lu) z[%lu,%lu)\n",
              g.populated_lo[0], g.populated_hi[0], g.populated_lo[1],
              g.populated_hi[1], g.populated_lo[2], g.populated_hi[2]);
  std::printf("expected agents/box: %.4f (assumed in the plan: %.2f)\n",
              agents_per_box, kAssumedAgentsPerBox);

  const uint64_t in_region =
      BruteForceCount(rm, g, g.populated_lo, g.populated_hi);
  std::printf("agents inside populated box range: %lu of %lu\n", in_region,
              kNumAgents);

  bool all_ok = true;
  std::printf("\n== per-tile results ==\n");
  std::printf("%4s %16s %8s %10s %7s %6s %6s\n", "k", "anchor", "actual",
              "expected", "z", "1sig", "exact");

  struct Summary {
    uint64_t k;
    size_t candidates;
    size_t sampled;
    real_t expected;
    real_t mean;
    real_t stddev;
    real_t z_mean;
    size_t in_band;
    bool exact;
    bool partition_checked;
    bool partition_ok;
    bool ok;
  };
  std::vector<Summary> summaries;

  for (auto k : kSweep) {
    auto candidates = EnumerateTiles(g, k);
    auto tiles = SampleTiles(candidates, kMaxTilesPerK, kSampleSeed + k);
    const real_t expected = agents_per_box * k * k * k;
    const real_t sigma = std::sqrt(expected);

    std::vector<real_t> counts;
    size_t in_band = 0;
    bool exact = true;
    bool near_zero = false;
    for (const auto& tile : tiles) {
      auto agents = ExtractTile(rm, env, tile);
      const auto n = agents.size();
      const auto brute = BruteForceCount(rm, g, tile.lo, TileHi(tile));
      const bool band = std::abs(n - expected) <= sigma;
      in_band += band;
      exact &= (brute == n);
      near_zero |= (n < 0.1 * expected);
      counts.push_back(n);
      std::printf("%4lu %5lu,%4lu,%4lu %8zu %10.1f %7.2f %6s %6s\n", k,
                  tile.lo[0], tile.lo[1], tile.lo[2], n, expected,
                  (n - expected) / sigma, band ? "yes" : "no",
                  brute == n ? "ok" : "DIFF");
    }

    real_t mean = 0;
    for (auto c : counts) {
      mean += c;
    }
    mean /= counts.size();
    real_t var = 0;
    for (auto c : counts) {
      var += (c - mean) * (c - mean);
    }
    const real_t stddev =
        counts.size() > 1 ? std::sqrt(var / (counts.size() - 1)) : NAN;
    const real_t z_mean = (mean - expected) / (sigma / std::sqrt(counts.size()));

    // When k divides the populated range on every axis, the candidate tiles
    // partition it exactly, so their counts must sum to the region total.
    bool partition_checked = true;
    for (int i = 0; i < 3; i++) {
      partition_checked &= (g.populated_hi[i] - g.populated_lo[i]) % k == 0;
    }
    bool partition_ok = true;
    if (partition_checked) {
      uint64_t sum = 0;
      for (const auto& tile : candidates) {
        sum += ExtractTile(rm, env, tile).size();
      }
      partition_ok = (sum == in_region);
    }

    const bool ok = exact && partition_ok && !near_zero && std::abs(z_mean) < 3;
    all_ok &= ok;
    summaries.push_back({k, candidates.size(), counts.size(), expected, mean,
                         stddev, z_mean, in_band, exact, partition_checked,
                         partition_ok, ok});
  }

  std::printf("\n== per-k summary ==\n");
  std::printf("%4s %6s %6s %10s %10s %9s %9s %8s %7s %7s %6s %9s %6s\n", "k",
              "cands", "tiles", "expected", "mean", "stddev", "sqrt(E)",
              "mean/E", "z_mean", "in1sig", "exact", "partition", "pass");
  for (const auto& s : summaries) {
    std::printf(
        "%4lu %6zu %6zu %10.1f %10.1f %9.1f %9.1f %8.4f %7.2f %4zu/%-2zu %6s "
        "%9s %6s\n",
        s.k, s.candidates, s.sampled, s.expected, s.mean, s.stddev,
        std::sqrt(s.expected), s.mean / s.expected, s.z_mean, s.in_band,
        s.sampled, s.exact ? "ok" : "DIFF",
        s.partition_checked ? (s.partition_ok ? "ok" : "FAIL") : "n/a",
        s.ok ? "yes" : "NO");
  }
  std::printf("\n%s\n", all_ok ? "VALIDATION PASSED" : "VALIDATION FAILED");
  return all_ok ? 0 : 1;
}

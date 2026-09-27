#ifndef TILE_EXTRACTION_H_
#define TILE_EXTRACTION_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "biodynamo.h"
#include "core/environment/uniform_grid_environment.h"

namespace bdm {
namespace tiling {

struct GridGeometry {
  std::array<uint32_t, 3> num_boxes_axis;
  // {x_min, x_max, y_min, y_max, z_min, z_max}, including the padding boxes.
  std::array<int32_t, 6> dimensions;
  int32_t box_length;
  // Per axis, the boxes lying fully inside [min_bound, max_bound]:
  // [populated_lo, populated_hi). Boxes outside this range are grid padding.
  std::array<uint64_t, 3> populated_lo;
  std::array<uint64_t, 3> populated_hi;
};

// Covers box coordinates [lo, lo + k) on each axis.
struct Tile {
  std::array<uint64_t, 3> lo;
  uint64_t k;
};

// ExecutionContext::AddAgent only buffers agents; they reach the
// ResourceManager when the execution contexts are set up. Before that the
// ResourceManager is empty and Environment::Update() is a no-op.
inline void CommitAgentsAndUpdateGrid(Simulation* sim) {
  sim->GetScheduler()->FinalizeInitialization();
  sim->GetEnvironment()->ForcedUpdate();
}

inline GridGeometry QueryGridGeometry(UniformGridEnvironment* env,
                                      real_t min_bound, real_t max_bound) {
  GridGeometry g;
  env->GetNumBoxesAxis(g.num_boxes_axis.data());
  g.dimensions = env->GetDimensions();
  g.box_length = env->GetBoxLength();
  const real_t length = g.box_length;
  for (int i = 0; i < 3; i++) {
    const real_t origin = g.dimensions[2 * i];
    g.populated_lo[i] = static_cast<uint64_t>(
        std::max<real_t>(0, std::ceil((min_bound - origin) / length)));
    g.populated_hi[i] =
        static_cast<uint64_t>(std::floor((max_bound - origin) / length));
  }
  return g;
}

class BoxRangeFilter : public Functor<bool, Agent*> {
 public:
  BoxRangeFilter(UniformGridEnvironment* env, const Tile& tile) : env_(env) {
    for (int i = 0; i < 3; i++) {
      lo_[i] = tile.lo[i];
      hi_[i] = tile.lo[i] + tile.k;
    }
  }

  bool operator()(Agent* agent) override {
    // Derived from the position rather than Agent::GetBoxIdx(), so the result
    // does not depend on when the grid last assigned box indices.
    auto c = env_->GetBoxCoordinates(env_->GetBoxIndex(agent->GetPosition()));
    return c[0] >= lo_[0] && c[0] < hi_[0] && c[1] >= lo_[1] &&
           c[1] < hi_[1] && c[2] >= lo_[2] && c[2] < hi_[2];
  }

 private:
  UniformGridEnvironment* env_;
  std::array<uint64_t, 3> lo_;
  std::array<uint64_t, 3> hi_;
};

inline std::vector<Agent*> ExtractTile(ResourceManager* rm,
                                       UniformGridEnvironment* env,
                                       const Tile& tile) {
  BoxRangeFilter filter(env, tile);
  std::vector<Agent*> agents;
  rm->ForEachAgent([&](Agent* agent) { agents.push_back(agent); }, &filter);
  return agents;
}

// All non-overlapping k-tiles inside the populated box range.
inline std::vector<Tile> EnumerateTiles(const GridGeometry& g, uint64_t k) {
  std::vector<Tile> tiles;
  for (uint64_t bz = g.populated_lo[2]; bz + k <= g.populated_hi[2]; bz += k) {
    for (uint64_t by = g.populated_lo[1]; by + k <= g.populated_hi[1];
         by += k) {
      for (uint64_t bx = g.populated_lo[0]; bx + k <= g.populated_hi[0];
           bx += k) {
        tiles.push_back({{bx, by, bz}, k});
      }
    }
  }
  return tiles;
}

inline std::vector<Tile> SampleTiles(std::vector<Tile> tiles,
                                     size_t max_tiles, uint64_t seed) {
  if (tiles.size() > max_tiles) {
    std::mt19937_64 rng(seed);
    std::shuffle(tiles.begin(), tiles.end(), rng);
    tiles.resize(max_tiles);
  }
  return tiles;
}

}  // namespace tiling
}  // namespace bdm

#endif  // TILE_EXTRACTION_H_

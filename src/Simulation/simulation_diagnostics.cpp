#include "simulation.h"
#include "Data/data_export.h"
#include "Data/logging.h"
#include "Data/param_parser.h"
#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

Simulation::ReversibilityEventExtremes::ReversibilityEventExtremes() {
  smallest.fill(DBL_MAX);
  largest.fill(-DBL_MAX);
}

bool Simulation::ReversibilityEventExtremes::consider(double eventSize) {
  if (!std::isfinite(eventSize) || eventSize <= 0.0) {
    throw std::runtime_error(
        "Reversibility event size must be finite and positive.");
  }

  bool isExtreme = false;
  auto largestSmallEvent = std::max_element(smallest.begin(), smallest.end());
  if (eventSize < *largestSmallEvent) {
    *largestSmallEvent = eventSize;
    isExtreme = true;
  }
  auto smallestLargeEvent = std::min_element(largest.begin(), largest.end());
  if (eventSize > *smallestLargeEvent) {
    *smallestLargeEvent = eventSize;
    isExtreme = true;
  }
  return isExtreme;
}

void Simulation::setReversibilityResult(bool reversible, double distance) {
  reversibilityState.isReversible = reversible ? 1 : 0;
  reversibilityState.distance = distance;
}

void Simulation::computeParticipationFraction() {
  if (beforeMinimization.displacements.empty() ||
      afterMinimization.displacements.empty()) {
    participationFraction = 0;
    return;
  }
  const size_t nodeCount = static_cast<size_t>(mesh.rows * mesh.cols);
  if (nodeCount == 0) {
    participationFraction = 0;
    return;
  }

  const size_t expectedSize = 2 * nodeCount;
  if (beforeMinimization.displacements.size() != expectedSize ||
      afterMinimization.displacements.size() != expectedSize) {
    throw std::runtime_error(DebugLog::formatDisplacementSnapshotPairSizeError(
        "Simulation::computeParticipationFraction",
        beforeMinimization.displacements.size(),
        afterMinimization.displacements.size(), expectedSize, mesh.rows,
        mesh.cols));
  }

  const double *beforeX = beforeMinimization.displacements.data();
  const double *beforeY = beforeMinimization.displacements.data() + nodeCount;
  const double *afterX = afterMinimization.displacements.data();
  const double *afterY = afterMinimization.displacements.data() + nodeCount;

  double sumU2 = 0.0;
  double sumU4 = 0.0;
  for (size_t i = 0; i < nodeCount; i++) {
    const double dx = afterX[i] - beforeX[i];
    const double dy = afterY[i] - beforeY[i];
    const double u2 = dx * dx + dy * dy; // ||U||^2 per node
    sumU2 += u2;
    sumU4 += u2 * u2;
  }

  if (sumU4 == 0.0) {
    participationFraction = 0;
    return;
  }

  const double N = static_cast<double>(nodeCount);
  participationFraction = (sumU2 * sumU2) / (N * sumU4);
}

void Simulation::computeM3ParticipationFraction() {
  if (beforeMinimization.displacements.empty() ||
      afterMinimization.displacements.empty()) {
    m3ParticipationFraction = 0;
    return;
  }
  const size_t nodeCount = static_cast<size_t>(mesh.rows * mesh.cols);
  if (nodeCount == 0) {
    m3ParticipationFraction = 0;
    return;
  }

  const size_t expectedSize = 2 * nodeCount;
  if (beforeMinimization.displacements.size() != expectedSize ||
      afterMinimization.displacements.size() != expectedSize) {
    throw std::runtime_error(DebugLog::formatDisplacementSnapshotPairSizeError(
        "Simulation::computeM3ParticipationFraction",
        beforeMinimization.displacements.size(),
        afterMinimization.displacements.size(), expectedSize, mesh.rows,
        mesh.cols));
  }

  const double *beforeX = beforeMinimization.displacements.data();
  const double *beforeY = beforeMinimization.displacements.data() + nodeCount;
  const double *afterX = afterMinimization.displacements.data();
  const double *afterY = afterMinimization.displacements.data() + nodeCount;

  double sumU2 = 0.0;
  double sumU4 = 0.0;
  size_t changedCount = 0;
  for (const TElement &e : mesh.elements) {
    if (e.m3Nr == e.pastM3Nr) {
      continue;
    }

    Vector2d du = Vector2d::Zero();
    for (const GhostNode &gn : e.ghostNodes) {
      const int nodeIndex = gn.referenceId.i;
      if (nodeIndex >= 0 && static_cast<size_t>(nodeIndex) < nodeCount) {
        const size_t idx = static_cast<size_t>(nodeIndex);
        du[0] += afterX[idx] - beforeX[idx];
        du[1] += afterY[idx] - beforeY[idx];
      }
    }
    du /= 3.0;
    const double u2 = du.squaredNorm();
    sumU2 += u2;
    sumU4 += u2 * u2;
    changedCount += 1;
  }

  if (changedCount == 0 || sumU4 == 0.0) {
    m3ParticipationFraction = 0;
    return;
  }

  const double N = static_cast<double>(changedCount);
  m3ParticipationFraction = (sumU2 * sumU2) / (N * sumU4);
}

void Simulation::updateEnergyHistory(bool endOfStep) {
  SimulationEnergyHistory &h = energyHistory;
  const double totalEnergy = mesh.totalEnergy;
  const double averageEnergy = mesh.averageEnergy;
  const double averageSigma11 = mesh.averageSigma11;
  const double averageSigma12 = mesh.averageSigma12;
  const double averageSigma22 = mesh.averageSigma22;
  const double averageP11 = mesh.averageP11;
  const double averageP12 = mesh.averageP12;
  const double averageP21 = mesh.averageP21;
  const double averageP22 = mesh.averageP22;

  if (mesh.nrMinFunctionCalls == 0 && !endOfStep) {
    // Update initial guess energy values at the start of the minimization
    // Note these are after the first minimization step.
    h.initialGuessTotalEnergy = totalEnergy;
    h.initialGuessAverageEnergy = averageEnergy;
    h.initialGuessAverageSigma11 = averageSigma11;
    h.initialGuessAverageSigma12 = averageSigma12;
    h.initialGuessAverageSigma22 = averageSigma22;
    h.initialGuessAverageP11 = averageP11;
    h.initialGuessAverageP12 = averageP12;
    h.initialGuessAverageP21 = averageP21;
    h.initialGuessAverageP22 = averageP22;
  }

  h.totalEnergyChangeFromInitialGuess = totalEnergy - h.initialGuessTotalEnergy;
  h.averageEnergyChangeFromInitialGuess =
      averageEnergy - h.initialGuessAverageEnergy;
  h.averageSigma12ChangeFromInitialGuess =
      averageSigma12 - h.initialGuessAverageSigma12;

  h.loadStepTotalEnergyChange = totalEnergy - h.prevLoadStepTotalEnergy;
  h.loadStepAverageEnergyChange = averageEnergy - h.prevLoadStepAverageEnergy;
  if (mesh.loadSteps == 1) {
    h.loadStepTotalEnergyChange = 0;
    h.loadStepAverageEnergyChange = 0;
  }

  if (endOfStep) {
    // Prepare for next load step
    h.minIterTotalEnergyChange = 0;
    h.minIterAverageEnergyChange = 0;
    h.prevMinIterTotalEnergy = totalEnergy;
    h.prevMinIterAverageEnergy = averageEnergy;

    h.prevLoadStepTotalEnergy = totalEnergy;
    h.prevLoadStepAverageEnergy = averageEnergy;
  } else {
    // We are inside the minimization loop. Track changes between successive
    // minimization iterations so logDuringMinimization can report per-iteration
    // energy drops.
    if (mesh.nrMinFunctionCalls == 0) {
      // First call in this minimization: initialize the baseline so the first
      // delta is zero.
      h.minIterTotalEnergyChange = 0;
      h.minIterAverageEnergyChange = 0;
      h.prevMinIterTotalEnergy = totalEnergy;
      h.prevMinIterAverageEnergy = averageEnergy;
    } else {
      // Subsequent calls: compute delta relative to previous iteration and
      // advance the baseline.
      h.minIterTotalEnergyChange = totalEnergy - h.prevMinIterTotalEnergy;
      h.minIterAverageEnergyChange = averageEnergy - h.prevMinIterAverageEnergy;
      h.prevMinIterTotalEnergy = totalEnergy;
      h.prevMinIterAverageEnergy = averageEnergy;
    }
  }
}

void Simulation::logMinimizationState() {
  if (!config.logDuringMinimization) {
    return;
  }
  mesh.updateAveragesAndPlasticEvents();
  updateEnergyHistory(false);
  if (minCsvHeaders.empty()) {
    throw std::runtime_error("Minimization CSV headers are not initialized.");
  }
  writeToCsv(minCsvFile, *this, minCsvHeaders);
}

void Simulation::addReversibilityCsvColumns() {
  const auto addIfMissing = [this](const std::string &name, auto accessor) {
    if (!hasCsvColumn(name)) {
      addCsvColumn(name, accessor);
    }
  };

  addIfMissing("is_reversible", [](const Simulation &s) {
    return s.reversibilityState.isReversible;
  });

  addIfMissing("rev_u_diff", [](const Simulation &s) {
    return s.reversibilityState.distance;
  });

  addIfMissing("rev_energy_diff", [](const Simulation &s) {
    return s.reversibilityState.energyDifference;
  });

  addIfMissing("rev_sigma_12_diff", [](const Simulation &s) {
    return s.reversibilityState.sigma12Difference;
  });

  addIfMissing("rev_sigma_trace_diff", [](const Simulation &s) {
    return s.reversibilityState.sigmaTraceDifference;
  });

  addIfMissing("rev_sigma11_diff", [](const Simulation &s) {
    return s.reversibilityState.sigma11Difference;
  });

  addIfMissing("rev_sigma22_diff", [](const Simulation &s) {
    return s.reversibilityState.sigma22Difference;
  });

  addIfMissing("rev_p11_diff", [](const Simulation &s) {
    return s.reversibilityState.p11Difference;
  });

  addIfMissing("rev_p12_diff", [](const Simulation &s) {
    return s.reversibilityState.p12Difference;
  });

  addIfMissing("rev_p21_diff", [](const Simulation &s) {
    return s.reversibilityState.p21Difference;
  });

  addIfMissing("rev_p22_diff", [](const Simulation &s) {
    return s.reversibilityState.p22Difference;
  });
}

bool Simulation::checkReversibility(const Matrix2d &stepTransform, double eps) {
  if (config.saveElasticReversibilityStates &&
      config.maximumSavedElasticReversibilityStates <= 0) {
    throw std::invalid_argument(
        "saveElasticReversibilityStates requires a positive "
        "maximumSavedElasticReversibilityStates");
  }
  Mesh &state0 = reversibilityState0;
  Mesh &state2 = reversibilityState2;

  state0 = mesh;

  // Save initial energy and sigma values for later comparison
  const double initialEnergy = mesh.totalEnergy;
  const double initialSigma12 = mesh.averageSigma12;
  const double initialSigma11 = mesh.averageSigma11;
  const double initialSigma22 = mesh.averageSigma22;
  const double initialSigmaTrace = mesh.averageSigmaTrace;
  const double initialP11 = mesh.averageP11;
  const double initialP12 = mesh.averageP12;
  const double initialP21 = mesh.averageP21;
  const double initialP22 = mesh.averageP22;

  // Forward step (0 -> 1 -> 2)
  applyAffineStep(stepTransform); // (0 -> 1)
  mesh.updateMesh();
  const double energyAffine = mesh.totalEnergy;
  minimize(); // (1 -> 2)

  mesh.updateAveragesAndPlasticEvents();
  const bool hasPlasticEvents = mesh.nr_elements_with_m3_changeInStep > 0;
  const bool hasEnergyDrop = mesh.totalEnergy < energyAffine;
  const double forwardEnergyDrop = energyAffine - mesh.totalEnergy;
  const bool saveThisElasticState =
      config.saveElasticReversibilityStates && !hasPlasticEvents &&
      savedElasticReversibilityStateCount <
          config.maximumSavedElasticReversibilityStates;

  if (!(hasPlasticEvents && hasEnergyDrop)) {
    // Ordinary simulations intentionally skip the backward test here.  A
    // targeted replay can opt in to saving no-forward-m3 events so that these
    // otherwise unmeasured reversible-elastic examples can be inspected.
    if (saveThisElasticState) {
      state2 = mesh;
      SimulationEnergyHistory historyState2 = energyHistory;
      SimReport FIRERepState2 = FIRERep;
      SimReport LBFGSRepState2 = LBFGSRep;
      SimReport CGRepState2 = CGRep;
      int nrMinItterationsState2 = mesh.nrMinItterations;
      int nrMinFunctionCallsState2 = mesh.nrMinFunctionCalls;

      const double oldIncrement = loadIncrement;
      loadIncrement = -oldIncrement;
      const Matrix2d backTransform = stepTransform.inverse();
      applyAffineStep(backTransform);
      minimize();
      mesh.updateForceStateAveragesAndPlasticEvents();

      const double d = mesh.rmsDistanceToMesh(state0, true);
      reversibilityState.distance = d;
      reversibilityState.energyDifference =
          std::abs(mesh.totalEnergy - initialEnergy);
      reversibilityState.sigma12Difference =
          std::abs(mesh.averageSigma12 - initialSigma12);
      reversibilityState.sigmaTraceDifference =
          std::abs(mesh.averageSigmaTrace - initialSigmaTrace);
      reversibilityState.sigma11Difference =
          std::abs(mesh.averageSigma11 - initialSigma11);
      reversibilityState.sigma22Difference =
          std::abs(mesh.averageSigma22 - initialSigma22);
      reversibilityState.p11Difference = std::abs(mesh.averageP11 - initialP11);
      reversibilityState.p12Difference = std::abs(mesh.averageP12 - initialP12);
      reversibilityState.p21Difference = std::abs(mesh.averageP21 - initialP21);
      reversibilityState.p22Difference = std::abs(mesh.averageP22 - initialP22);
      reversibilityState.isReversible = d < eps ? 1 : 0;

      const std::string dropSubFolder =
          std::string("reversibilityData/elastic_replay_l_") +
          [&]() {
            const double step = std::abs(oldIncrement);
            const int precision = std::max(
                0, static_cast<int>(std::ceil(-std::log10(step))));
            std::ostringstream loadString;
            loadString << std::fixed << std::setprecision(precision)
                       << mesh.load;
            return loadString.str();
          }();
      const std::string dropPath =
          getDataPath(simName, dataPath) + "/" + dropSubFolder;
      auto writeState = [&](Mesh &state, const std::string &name) {
        state.ensureFull();
        writeMeshToVtu(state, simName, dataPath, name, false,
                       VtuFieldLevel::All, "", dropSubFolder);
      };
      Mesh affineScratch;
      auto reconstructAffineState = [&](const Mesh &base, const Matrix2d &T,
                                        double loadDelta) -> Mesh & {
        affineScratch = base;
        affineScratch.addLoad(loadDelta);
        if (affineScratch.usingPBC) {
          affineScratch.applyTransformationToSystemDeformation(T);
        } else {
          affineScratch.applyTransformationToFixedNodes(T);
        }
        const Matrix2d I = Matrix2d::Identity();
        const Matrix2d A = T - I;
        for (const NodeId &n_id : affineScratch.freeNodeIds) {
          Node *n = affineScratch[n_id];
          const Vector2d nextDisplacement = A * n->ref_pos() + T * n->u();
          n->setDisplacement(nextDisplacement);
        }
        affineScratch.markDirty();
        return affineScratch;
      };
      writeState(state0, "state0_min_gamma");
      writeState(reconstructAffineState(state0, stepTransform, oldIncrement),
                 "state1_affine_gamma_plus");
      writeState(state2, "state2_relaxed_gamma_plus");
      writeState(reconstructAffineState(state2, backTransform, -oldIncrement),
                 "state3_affine_gamma_minus");
      mesh.ensureFull();
      writeMeshToVtu(mesh, simName, dataPath, "state4_relaxed_gamma", false,
                     VtuFieldLevel::All, "", dropSubFolder);
      createCollection(dropPath, dropPath);
      ++savedElasticReversibilityStateCount;

      // Continue the forward simulation from state 2, exactly as the normal
      // branch does after a measured backward test.
      mesh = state2;
      energyHistory = historyState2;
      FIRERep = FIRERepState2;
      LBFGSRep = LBFGSRepState2;
      CGRep = CGRepState2;
      mesh.nrMinItterations = nrMinItterationsState2;
      mesh.nrMinFunctionCalls = nrMinFunctionCallsState2;
      loadIncrement = oldIncrement;
    }
    if (!saveThisElasticState) {
      reversibilityState.isReversible = 1;
      reversibilityState.distance = 0;
      reversibilityState.energyDifference = 0;
      reversibilityState.sigma12Difference = 0;
      reversibilityState.sigmaTraceDifference = 0;
      reversibilityState.sigma11Difference = 0;
      reversibilityState.sigma22Difference = 0;
      reversibilityState.p11Difference = 0;
      reversibilityState.p12Difference = 0;
      reversibilityState.p21Difference = 0;
      reversibilityState.p22Difference = 0;
    }
    if (saveThisElasticState) {
      return reversibilityState.isReversible == 1;
    }
    return true;
  }

  state2 = mesh;
  SimulationEnergyHistory historyState2 = energyHistory;
  SimReport FIRERepState2 = FIRERep;
  SimReport LBFGSRepState2 = LBFGSRep;
  SimReport CGRepState2 = CGRep;
  int nrMinItterationsState2 = mesh.nrMinItterations;
  int nrMinFunctionCallsState2 = mesh.nrMinFunctionCalls;

  // Backward step (2 -> 3 -> 4)
  const double oldIncrement = loadIncrement;
  loadIncrement = -oldIncrement;
  const Matrix2d backTransform = stepTransform.inverse();
  applyAffineStep(backTransform);
  minimize();
  mesh.updateForceStateAveragesAndPlasticEvents();

  const double d = mesh.rmsDistanceToMesh(state0, true);

  // Save reversibility results
  reversibilityState.distance = d;
  reversibilityState.energyDifference =
      std::abs(mesh.totalEnergy - initialEnergy);
  reversibilityState.sigma12Difference =
      std::abs(mesh.averageSigma12 - initialSigma12);
  reversibilityState.sigmaTraceDifference =
      std::abs(mesh.averageSigmaTrace - initialSigmaTrace);
  reversibilityState.sigma11Difference =
      std::abs(mesh.averageSigma11 - initialSigma11);
  reversibilityState.sigma22Difference =
      std::abs(mesh.averageSigma22 - initialSigma22);
  reversibilityState.p11Difference = std::abs(mesh.averageP11 - initialP11);
  reversibilityState.p12Difference = std::abs(mesh.averageP12 - initialP12);
  reversibilityState.p21Difference = std::abs(mesh.averageP21 - initialP21);
  reversibilityState.p22Difference = std::abs(mesh.averageP22 - initialP22);

  reversibilityState.isReversible = d < eps ? 1 : 0;
  const bool reversible = d < eps;

  // Save drop snapshots (every 300 reversible drops, every 10 irreversible).
  auto saveDrop = [&](const std::string &folder) {
    const std::string dropSubFolder =
        std::string("reversibilityData/") + folder;
    const std::string dropPath =
        getDataPath(simName, dataPath) + "/" + dropSubFolder;
    auto writeState = [&](Mesh &state, const std::string &name) {
      state.ensureFull();
      writeMeshToVtu(state, simName, dataPath, name, false, VtuFieldLevel::All,
                     "", dropSubFolder);
    };
    Mesh affineScratch;
    auto reconstructAffineState = [&](const Mesh &base, const Matrix2d &T,
                                      double loadDelta) -> Mesh & {
      affineScratch = base;
      affineScratch.addLoad(loadDelta);
      if (affineScratch.usingPBC) {
        affineScratch.applyTransformationToSystemDeformation(T);
      } else {
        affineScratch.applyTransformationToFixedNodes(T);
      }
      const Matrix2d I = Matrix2d::Identity();
      const Matrix2d A = T - I;
      for (const NodeId &n_id : affineScratch.freeNodeIds) {
        Node *n = affineScratch[n_id];
        const Vector2d nextDisplacement = A * n->ref_pos() + T * n->u();
        n->setDisplacement(nextDisplacement);
      }
      affineScratch.markDirty();
      return affineScratch;
    };
    writeState(state0, "state0_min_gamma");
    writeState(reconstructAffineState(state0, stepTransform, oldIncrement),
               "state1_affine_gamma_plus");
    writeState(state2, "state2_relaxed_gamma_plus");
    writeState(reconstructAffineState(state2, backTransform, -oldIncrement),
               "state3_affine_gamma_minus");
    mesh.ensureFull();
    writeMeshToVtu(mesh, simName, dataPath, "state4_relaxed_gamma", false,
                   VtuFieldLevel::All, "", dropSubFolder);
    createCollection(dropPath, dropPath);
  };

  const double step = std::abs(oldIncrement);
  assert(step > 0);
  int precision = std::max(0, static_cast<int>(std::ceil(-std::log10(step))));
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(precision) << mesh.load;
  ReversibilityEventExtremes &eventExtremes =
      reversible ? reversibleEventExtremes : irreversibleEventExtremes;
  const bool saveExtremeEvent = eventExtremes.consider(forwardEnergyDrop);
  int &dropCount = reversible ? reversibleDropCount : irreversibleDropCount;
  ++dropCount;
  const int saveEvery = reversible ? 300 : 10;
  const bool savePeriodicEvent = dropCount % saveEvery == 0;
  if (saveExtremeEvent || savePeriodicEvent) {
    saveDrop(std::string(reversible ? "rev_drop_l_" : "irrev_drop_l_") +
             oss.str());
  }

  // Restore to state 2 (forward relaxed state)
  mesh = state2;
  energyHistory = historyState2;
  FIRERep = FIRERepState2;
  LBFGSRep = LBFGSRepState2;
  CGRep = CGRepState2;
  mesh.nrMinItterations = nrMinItterationsState2;
  mesh.nrMinFunctionCalls = nrMinFunctionCallsState2;
  loadIncrement = oldIncrement;

  return reversible;
}

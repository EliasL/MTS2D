#include "simulation.h"
#include "Data/logging.h"
#include <FIRE.h>
#include <Param.h>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

// Solver setup, callbacks, and energy minimization live here. The loading and
// reconnection workflow is in simulation.cpp; tangent assembly is separate.

void Simulation::initSolver() {
  stiffnessPreconditioner.reset();

  // Alglib Initialization preparation
  int nrFreeNodes = mesh.freeNodeIds.size();
  alglibNodeDisplacements.setlength(2 * nrFreeNodes);
  // Set values to zero
  for (int i = 0; i < 2 * nrFreeNodes; i++) {
    alglibNodeDisplacements[i] = 0;
  }

  // Adjust nrCorrections based on the number of free nodes
  // alglib doesn't like that nr corrections is larger than nr of free nodes
  if (config.LBFGSNrCorrections > 2 * nrFreeNodes) {
    config.LBFGSNrCorrections = 2 * nrFreeNodes;
  }

  if (nrFreeNodes != 0) {
    resetLBFGSState(config.LBFGSNrCorrections);

    alglib::mincgcreate(alglibNodeDisplacements, CG_state);
  }

  // FIRE Initialization
  FIRENodeDisplacements = VectorXd::Zero(2 * nrFreeNodes);
  FIRE_param = FIREpp::FIREParam<double>(mesh.nrNodes, 2);
  config.updateParam(FIRE_param);

  // We update the datalink with the latest information about out simulation
  // data
  dataLink = DataLink(this);
  // TEMP
  // config.logDuringMinimization = true;
}

void Simulation::resetLBFGSState(int corrections) {
  const int problemDimension = alglibNodeDisplacements.length();
  if (problemDimension == 0) {
    return;
  }
  if (corrections < 1 || corrections > problemDimension) {
    throw std::invalid_argument(
        "Simulation::resetLBFGSState: corrections must be between 1 and " +
        std::to_string(problemDimension) + ".");
  }

  // Recreating the state is necessary: changing Config alone does not change
  // the M value stored inside ALGLIB's optimizer state.
  alglib::minlbfgscreate(corrections, alglibNodeDisplacements, LBFGS_state);
}

StiffnessPreconditioner *Simulation::preparePreconditioner() {
  config.validatePreconditioner();
  if (config.preconditioner == "none") {
    stiffnessPreconditioner.reset();
    return nullptr;
  }
  const auto kind = config.preconditioner == "diagonal"
                        ? PreconditionerKind::Diagonal
                        : PreconditionerKind::ExactTangent;

  bool rebuild = !stiffnessPreconditioner ||
                 !stiffnessPreconditioner->compatibleWith(mesh, kind) ||
                 config.preconditionerUpdate == "every";
  if (!rebuild && config.preconditionerUpdate == "adaptive") {
    rebuild = stiffnessPreconditioner->maxDisplacementSinceBuild(
                  alglibNodeDisplacements) >=
              config.preconditionerDisplacementThreshold;
  }
  if (rebuild) {
    const auto start = std::chrono::steady_clock::now();
    // Follow Sylvain's choice: use the current mesh tangent before initial-
    // guess noise. Large noisy perturbations can make a useful elastic tangent
    // indefinite. The optimizer still evaluates the full noisy initial guess.
    mesh.ensureForces();
    auto candidate = std::make_unique<StiffnessPreconditioner>(
        config.selectedTangentBackend());
    const bool success = candidate->build(mesh, alglibNodeDisplacements, kind);
    ++preconditionerStats.builds;
    preconditionerStats.buildMilliseconds +=
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();
    if (!success) {
      ++preconditionerStats.fallbacks;
      stiffnessPreconditioner.reset();
      return nullptr;
    }
    stiffnessPreconditioner = std::move(candidate);
  } else {
    ++preconditionerStats.reuses;
  }
  stiffnessPreconditioner->startAt(alglibNodeDisplacements);
  return stiffnessPreconditioner.get();
}

namespace {

struct OptimizerCoordinates {
  DataLink *data;
  StiffnessPreconditioner *preconditioner;
  alglib::real_1d_array physical;
  alglib::real_1d_array gradient;
  alglib::real_1d_array lastEvaluated;
};

void optimizerEnergyAndGradient(const alglib::real_1d_array &x, double &energy,
                               alglib::real_1d_array &gradient, void *pointer) {
  auto &coordinates = *static_cast<OptimizerCoordinates *>(pointer);
  if (coordinates.preconditioner == nullptr) {
    alglibEnergyAndGradient(x, energy, gradient, coordinates.data);
    coordinates.lastEvaluated = x;
    return;
  }
  // ALGLIB sees z, but all mesh calculations and the epsR stopping check use
  // physical displacements and forces. Only the optimizer's gradient is mapped.
  coordinates.preconditioner->transformedToPhysical(x, coordinates.physical);
  coordinates.gradient.setlength(coordinates.physical.length());
  alglibEnergyAndGradient(coordinates.physical, energy, coordinates.gradient,
                         coordinates.data);
  coordinates.preconditioner->physicalGradientToTransformed(
      coordinates.gradient, gradient);
  coordinates.lastEvaluated = x;
}

void optimizerIterationLogger(const alglib::real_1d_array &x, double,
                             void *pointer) {
  auto &coordinates = *static_cast<OptimizerCoordinates *>(pointer);
  // Reporting callbacks describe accepted iterates. A force evaluation may
  // instead be a line-search trial, and requesting termination there makes
  // ALGLIB return the previous accepted iterate.
  bool current = x.length() == coordinates.lastEvaluated.length();
  for (int i = 0; current && i < x.length(); ++i) {
    current = x[i] == coordinates.lastEvaluated[i];
  }
  if (!current) {
    if (coordinates.preconditioner != nullptr) {
      coordinates.preconditioner->transformedToPhysical(x, coordinates.physical);
    } else {
      coordinates.physical = x;
    }
    updateNodePositions(*coordinates.data->mesh, coordinates.physical);
    coordinates.data->mesh->ensureForces();
  }
  if (coordinates.data->mesh->maxForce < *coordinates.data->maxForceAllowed) {
    alglib::minlbfgsrequesttermination(*coordinates.data->LBFGS_state);
    alglib::mincgrequesttermination(*coordinates.data->CG_state);
  }
  // Evaluation logging is already handled by updateMeshAndComputeForces.
}

} // namespace

// Helper function
void Simulation::m_minimize() {
  auto fail = [&](std::string_view caughtType, const std::string &msg) -> void {
    const std::string error = DebugLog::formatMinimizationFailure(
        caughtType, msg, mesh, config.minimizer, false);
    if (!isQuiet() && !debugReplayActive) {
      std::cerr << error << '\n';
    }
    writeToFile(true, "CrashAtLoad:" + std::to_string(mesh.load));
    throw std::runtime_error(error);
  };

  int nrIter = 0;
  try {
    if (config.minimizer == "FIRE") {
      m_minimizeWithFIRE();
      nrIter = FIRERep.nrIter;
    } else if (config.minimizer == "LBFGS") {
      m_minimizeWithLBFGS();
      nrIter = LBFGSRep.nrIter;
    } else if (config.minimizer == "CG") {
      m_minimizeWithCG();
      nrIter = CGRep.nrIter;
    } else {
      throw std::invalid_argument("Unknown minimizer: " + config.minimizer);
    }
    mesh.nrMinItterations += nrIter;
    mesh.nrMinItterationsSinceLastReconnect =
        std::max(mesh.nrMinItterationsSinceLastReconnect, nrIter);
  } catch (const alglib::ap_error &e) {
    fail("ALGLIB error", e.msg);
  } catch (const std::exception &e) {
    fail("Standard exception", e.what());
  } catch (...) {
    fail("Unknown exception", "");
  }
  // if (FIRERep.termType == -3) {
  //  writeToFile(true);
  //  throw std::runtime_error("Energy too high");
  //}
  if (LBFGSRep.termType == 1) {
    // mesh.writeToVtu("badStopStep" + std::to_string(mesh.loadSteps));
    // If maxForceAllowed is smaller than 1e-19, we assume it is zero and not
    // used
    if (mesh.maxForce > 1.5 * (*dataLink.maxForceAllowed) &&
        *dataLink.maxForceAllowed > 1e-19) {
      if (!isQuiet()) {
        std::cout << "Warning: LBFGS stopped with termType 1 but max force is "
                     "still high: "
                  << mesh.maxForce << '\n';
      }
    }
  }
}

void Simulation::m_minimizeWithLBFGS() {
  timer.Start("LBFGSMinimization");
  auto *preconditioner = preparePreconditioner();
  auto *coordinatePreconditioner =
      preconditioner != nullptr && preconditioner->usesTransformedCoordinates()
          ? preconditioner
          : nullptr;
  alglib::real_1d_array transformed;
  auto *coordinates = &alglibNodeDisplacements;
  if (coordinatePreconditioner != nullptr) {
    coordinatePreconditioner->physicalToTransformed(alglibNodeDisplacements,
                                                   transformed);
    coordinates = &transformed;
  }
  const int corrections = std::min<int>(config.LBFGSNrCorrections,
                                       coordinates->length());
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsrestartfrom
  // We reset and reuse the state instead of initializing it again
  // (The hessian is reset and not preserved)
  if (LBFGS_state.c_ptr()->n != coordinates->length() ||
      LBFGS_state.c_ptr()->m != corrections) {
    alglib::minlbfgscreate(corrections, *coordinates, LBFGS_state);
  } else {
    alglib::minlbfgsrestartfrom(LBFGS_state, *coordinates);
  }
  // Restart preserves ALGLIB's previous preconditioner. Set it explicitly so
  // switching modes, or falling back after a failed build, cannot reuse it.
  if (preconditioner != nullptr && coordinatePreconditioner == nullptr) {
    alglib::minlbfgssetprecdiag(LBFGS_state, preconditioner->diagonalValues());
  } else {
    alglib::minlbfgssetprecdefault(LBFGS_state);
  }

  // Set termination condition, ei. when is the solution good enough
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
  alglib::minlbfgssetcond(LBFGS_state, config.LBFGSEpsg, config.LBFGSEpsf,
                          config.LBFGSEpsx, config.LBFGSMaxIterations);
  alglib::minlbfgssetxrep(LBFGS_state, true);

  // Connect the chosen state to the minimization state
  minState = MinState(LBFGS_state);

  //  This is where the heavy calculations happen
  OptimizerCoordinates context{&dataLink, coordinatePreconditioner, {}, {}, {}};
  alglib::minlbfgsoptimize(LBFGS_state, optimizerEnergyAndGradient,
                         optimizerIterationLogger, &context);

  alglib::minlbfgsresults(LBFGS_state, *coordinates, LBFGS_report);
  if (coordinatePreconditioner != nullptr) {
    coordinatePreconditioner->transformedToPhysical(transformed,
                                                   alglibNodeDisplacements);
  }
  FIRENodeDisplacements = Eigen::Map<const VectorXd>(
      alglibNodeDisplacements.getcontent(), alglibNodeDisplacements.length());
  updateNodePositions(mesh, alglibNodeDisplacements);
  mesh.ensureForces();
  LBFGSRep = SimReport(LBFGS_report);
  LBFGSRep.nms = timer.Stop("LBFGSMinimization");
}

void Simulation::m_minimizeWithCG() {
  timer.Start("CGMinimization");
  auto *preconditioner = preparePreconditioner();
  auto *coordinatePreconditioner =
      preconditioner != nullptr && preconditioner->usesTransformedCoordinates()
          ? preconditioner
          : nullptr;
  alglib::real_1d_array transformed;
  auto *coordinates = &alglibNodeDisplacements;
  if (coordinatePreconditioner != nullptr) {
    coordinatePreconditioner->physicalToTransformed(alglibNodeDisplacements,
                                                   transformed);
    coordinates = &transformed;
  }
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_mincgrestartfrom
  // We reset and reuse the state instead of initializing it again
  if (CG_state.c_ptr()->n != coordinates->length()) {
    alglib::mincgcreate(*coordinates, CG_state);
  } else {
    alglib::mincgrestartfrom(CG_state, *coordinates);
  }
  if (preconditioner != nullptr && coordinatePreconditioner == nullptr) {
    alglib::mincgsetprecdiag(CG_state, preconditioner->diagonalValues());
  } else {
    alglib::mincgsetprecdefault(CG_state);
  }

  // Set termination condition, ei. when is the solution good enough
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_mincgsetcond
  alglib::mincgsetcond(CG_state, config.CGEpsg, config.CGEpsf, config.CGEpsx,
                       config.CGMaxIterations);
  alglib::mincgsetxrep(CG_state, true);

  // TODO!
  // minState = MinState(CG_state);

  // This is where the heavy calculations happen
  // The null pointer can be replaced with a logging function
  // https://www.alglib.net/translator/man/manual.cpp.html#sub_mincgoptimize
  OptimizerCoordinates context{&dataLink, coordinatePreconditioner, {}, {}, {}};
  alglib::mincgoptimize(CG_state, optimizerEnergyAndGradient,
                      optimizerIterationLogger, &context);

  alglib::mincgresults(CG_state, *coordinates, CG_report);
  if (coordinatePreconditioner != nullptr) {
    coordinatePreconditioner->transformedToPhysical(transformed,
                                                   alglibNodeDisplacements);
  }
  FIRENodeDisplacements = Eigen::Map<const VectorXd>(
      alglibNodeDisplacements.getcontent(), alglibNodeDisplacements.length());
  updateNodePositions(mesh, alglibNodeDisplacements);
  mesh.ensureForces();
  CGRep = SimReport(CG_report);
  CGRep.nms = timer.Stop("CGMinimization");
}

template <typename ArrayType>
void updateMeshAndComputeForces(DataLink *dataLink, const ArrayType &disp,
                                double &energy, ArrayType &grad,
                                int nr_x_values) {
  Mesh *mesh = dataLink->mesh;

  // Update mesh position from the result of the
  // previous minimization
  updateNodePositions(*mesh, disp);

  // Calculate energy and forces
  mesh->updateMesh();
  assert(mesh->getUpdateState() == Mesh::UpdateState::Forces);

  // Total energy, only used for minimization
  energy = mesh->totalEnergy;

  // Update gradient in the minimization
  updateGradArray(mesh, grad, nr_x_values);
  // ALGLIB checks epsR in optimizerIterationLogger, after accepting a step.
  // FIRE applies its own physical force stopping criterion.
  // We start to reconnect once we are 'close' to a solution
  // And only reconnect every 100 iterations
  // else if (mesh->nr_elements_with_m3_fix_changeInStep > 0) {

  //   mesh->reconnect(false, true);
  //   if (mesh->reconnectRequired) {
  //     // stop the minimization
  //     alglib::minlbfgsrequesttermination(*dataLink->LBFGS_state);
  //     alglib::mincgrequesttermination(*dataLink->CG_state);
  //   }
  // }
  // Since we don't use the x displacement argument in the iteration logger,
  // we just pass default parameter
  iterationLogger(alglib::real_1d_array(), energy, dataLink);
  // Update pastPlastic count so we can save files when it changes
  // mesh->resetPastPlasticCount(false);
  if (mesh->nrMinFunctionCalls == 0 &&
      !dataLink->s->config.logDuringMinimization) {
    // Ensure initial-guess stats are computed.
    // If we log during minimization, these will be updated in the iteration
    // logger instead
    mesh->updateForceStateAveragesAndPlasticEvents();
    dataLink->s->updateEnergyHistory(false);
  }
  mesh->nrMinFunctionCalls++;
  mesh->nrMinFunctionCallsSinceLastReconnect++;
}

void alglibEnergyAndGradient(const alglib::real_1d_array &disp, double &energy,
                             alglib::real_1d_array &grad, void *dataLinkPtr) {
  DataLink *dataLink = reinterpret_cast<DataLink *>(dataLinkPtr);
  updateMeshAndComputeForces(dataLink, disp, energy, grad, grad.length() / 2);
}

double FIREEnergyAndGradient(Eigen::VectorXd &disp, Eigen::VectorXd &grad,
                             void *dataLinkPtr) {
  double energy;
  DataLink *dataLink = reinterpret_cast<DataLink *>(dataLinkPtr);
  updateMeshAndComputeForces(dataLink, disp, energy, grad, grad.size() / 2);
  return energy;
}

void Simulation::m_minimizeWithFIRE() {
  timer.Start("FIREMinimization");
  config.updateParam(FIRE_param);
  FIREpp::FIRESolver<double> s = FIREpp::FIRESolver<double>(FIRE_param);
  double energy;
  FIRERep.nrIter = s.minimize(FIREEnergyAndGradient, FIRENodeDisplacements,
                              energy, &dataLink, FIRERep.termType);
  FIRERep.nms = timer.Stop("FIREMinimization");
  FIRERep.nfev = mesh.nrMinFunctionCalls;

  if (FIRERep.termType == -3) {
    writeToFile(true);
  }
  // We first do FIRE, but then use the result as an initial guess for LBFGS
  // setInitialGuess();
  // m_minimizeWithLBFGS();
}

// Overload for alglib::real_1d_array
void updateNodePositions(Mesh &mesh, const alglib::real_1d_array &disp) {
  mesh.updateNodePositions(disp.getcontent(), disp.length());
}

// Overload for Eigen::VectorXd
void updateNodePositions(Mesh &mesh, const Eigen::VectorXd &disp) {
  mesh.updateNodePositions(disp.data(), disp.size());
}

// Updates gradient array from current node forces.
template <typename ArrayType>
void updateGradArray(Mesh *mesh, ArrayType &grad, int nr_x_values) {
  for (int i = 0; i < nr_x_values; i++) {
    NodeId n_id = mesh->freeNodeIds[i];
    const auto &node = (*mesh)[n_id];

    // Store force values in the array
    double fx = node->f[0];
    double fy = node->f[1];

    // Gadient is the oposite sign of the force
    grad[i] = -fx;
    grad[nr_x_values + i] = -fy;
  }
}

// This can be used in the LBFGS optimization
void iterationLogger(const alglib::real_1d_array &x, double energy,
                     void *dataLinkPtr) {
  (void)x; // Explicitly casting x to void to silence unused parameter warning

  DataLink *dataLink = reinterpret_cast<DataLink *>(dataLinkPtr);

  // dataPath looks like: "/Volumes/data/MTS2D_output/"

  Mesh *mesh = dataLink->mesh;

  // auto state = dataLink->minState;

  // Increment the iteration count
  // TODO Does not work for CG or FIRE, only for LBFGS
  int it = dataLink->LBFGS_state->c_ptr()->repiterationscount;
  int nrFc = mesh->nrMinFunctionCalls;
  mesh->nrMinItterationsSinceLastReconnect =
      std::max(mesh->nrMinItterationsSinceLastReconnect, it);

  // Check if iteration count is a multiple of 5000
  if (nrFc % 5000 == 0 && nrFc > 0) {
    if (!isQuiet()) {
      std::cout << "Warning (step " << mesh->loadSteps << "): " << it
                << " iterations and " << nrFc << " function calls."
                << std::endl;
    }
  }

  // Save when energy shifts by >10% vs last saved, otherwise every 1000
  // function calls.
  int saveEvery = 1000;
  static int lastSavedFc = -1;
  static bool hasLastSavedEnergy = false;
  static double lastSavedEnergy = 0.0;

  if (nrFc == 0) {
    lastSavedFc = -1;
    hasLastSavedEnergy = false;
    lastSavedEnergy = 0.0;
  }

  if (dataLink->s->config.logDuringMinimization) {
    dataLink->s->timer.Start("write");
    // Write to the CSV file
    dataLink->s->logMinimizationState();
    bool energyJump = false;
    if (!hasLastSavedEnergy) {
      energyJump = true;
    } else {
      double denom = std::max(std::abs(lastSavedEnergy), 1e-12);
      energyJump = std::abs(energy - lastSavedEnergy) / denom > 0.10;
    }

    bool stepSave = (lastSavedFc < 0) || (abs(nrFc - lastSavedFc) >= saveEvery);

    if (energyJump || stepSave) {
      // Write mesh to file
      const VtuFieldLevel fieldLevel =
          dataLink->s->config.fullMinimizationLogging ? VtuFieldLevel::All
                                                      : VtuFieldLevel::Minimal;
      mesh->writeToVtu("", true, fieldLevel);
      lastSavedFc = nrFc;
      lastSavedEnergy = energy;
      hasLastSavedEnergy = true;
    }
    dataLink->s->timer.Stop("write");
  }
}

DataLink::DataLink(Simulation *simulation) {
  s = simulation;
  mesh = &s->mesh;
  stopSignal = nullptr; // This is set right before minimization
  minState = &s->minState;
  LBFGS_state = &s->LBFGS_state;
  CG_state = &s->CG_state;
  maxForceAllowed = &s->config.epsR;
}

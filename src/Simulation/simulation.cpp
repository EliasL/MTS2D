#include "simulation.h"

Simulation::Simulation(std::string configFile, std::optional<std::string> _dataPath)
{
    if (!_dataPath)
    {
        dataPath = findOutputPath();
    }
    else
    {
        dataPath = _dataPath.value();
    }

    Config conf = getConf(configFile);

    // We fix the random seed to get reproducable results
    srand(conf.seed);
    // Set the the number of threads
    if (conf.nrThreads == 0)
    {
        conf.nrThreads = omp_get_max_threads();
    }
    omp_set_num_threads(conf.nrThreads);

    // Assign values from Config to Simulation members
    name = conf.name;
    nx = conf.nx;
    ny = conf.ny;

    startLoad = conf.startLoad;
    loadIncrement = conf.loadIncrement;
    maxLoad = conf.maxLoad;
    noise = conf.noise;

    nrCorrections = conf.nrCorrections;
    epsg = conf.epsg;
    epsf = conf.epsf;
    epsx = conf.epsx;
    maxIterations = conf.maxIterations;

    plasticityEventThreshold = conf.plasticityEventThreshold;
    timer = Timer();

    mesh = Mesh(nx, ny);
    int nrNonBorderNodes = mesh.freeNodeIds.size();
    nodeDisplacements.setlength(2 * nrNonBorderNodes);

    // Boundary conditon transformation
    loadStepTransform = getShear(loadIncrement);

    clearOutputFolder(name, dataPath);
    createDataFolder(name, dataPath);
    setLogFile(name, dataPath);
    writeMeshToCsv(mesh, name, dataPath, true);

    // Prepare initial load condition
    mesh.applyTransformation(getShear(startLoad));

    spdlog::info("Config:\n{}", conf.str());

    // Flushes the config info so that we have that very quickly
    spdlog::default_logger()->flush();
    std::cout << conf;
}

void Simulation::run_simulation()
{
    timer.Start();

    for (double load = startLoad; load < maxLoad; load += loadIncrement)
    {
        // This creates and updates a progress bar
        m_updateProgress(load);

        // We shift the boundary nodes according to the loadIncrement
        mesh.applyTransformationToFixedNodes(loadStepTransform);

        // Modifies nodeDisplacements
        m_initialGuess();
        // If it is the first step of the simulation
        if (load == startLoad)
        {
            // Give the first guess some noise
            addNoise(nodeDisplacements, noise);

            // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgscreate
            // Initialize the state
            alglib::minlbfgscreate(nrCorrections, nodeDisplacements, state);
        }

        // This is the minimization section
        m_minimize_with_alglib();
        // Then we write the current state to the disk
        m_writeToDisk(load);
    }

    // Some minor cleanup and create a collection of vtu files.
    m_exit();
}

void Simulation::m_minimize_with_alglib()
{
    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsrestartfrom
    // We reset and reuse the state instead of initializing it again
    alglib::minlbfgsrestartfrom(state, nodeDisplacements);

    // Set termination condition, ei. when is the solution good enough
    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgssetcond
    alglib::minlbfgssetcond(state, epsg, epsf, epsx, maxIterations);

    // This is where the heavy calculations happen
    // The null pointer can be replaced with a logging function
    alglib::minlbfgsoptimize(state, alglib_calc_energy_and_gradiant, nullptr, &mesh);

    // TODO Collecting and analysing these reports could be a usefull tool for optimization
    alglib::minlbfgsresults(state, nodeDisplacements, report);
}

void alglib_calc_energy_and_gradiant(const alglib::real_1d_array &displacement,
                                     double &energy,
                                     alglib::real_1d_array &force, void *meshPtr)
{

    // Cast the void pointer back to Mesh pointer
    Mesh *mesh = reinterpret_cast<Mesh *>(meshPtr);

    // We just need this for indexing the force array
    int nr_x_values = force.length() / 2;

    // Update mesh position
    updatePossitionOfMesh(*mesh, displacement);

    // Calculate energy and forces
    energy = calc_energy_and_forces(*mesh);

    // Update forces
    for (size_t i = 0; i < nr_x_values; i++)
    {
        // We only want to use the force on the free nodes
        NodeId n_id = mesh->freeNodeIds[i];
        force[i] = (*mesh)[n_id]->f_x;
        force[nr_x_values + i] = (*mesh)[n_id]->f_y;
    }
    // writeMeshToVtu(*mesh, "minimizing");
}

// Updates the forces on the nodes in the surface and returns the total
// energy from all the elements in the surface.
double calc_energy_and_forces(Mesh &mesh)
{
    // First of all we need to make sure that the forces on the nodes have been
    // reset
    mesh.resetForceOnNodes();

    // This is the total energy from all the triangles
    double total_energy = 0;

// TODO We could check if we can make total energy a reduced variable,
// and make addForce a omp critical function and test for performance gains
#pragma omp parallel
#pragma omp for
    for (size_t i = 0; i < mesh.nrElements; i++)
    {
        mesh.elements[i].update();
    }

    // Now that the forces on the nodes have been reset, and
    // the elements updated in parallel, we can sum up the energy
    // and forces.
    for (size_t i = 0; i < mesh.nrElements; i++)
    {
        mesh.elements[i].applyForcesOnNodes();
        total_energy += mesh.elements[i].energy;
    }
    mesh.averageEnergy = total_energy / mesh.nrElements;
    return total_energy;
}

void updatePossitionOfMesh(Mesh &mesh, const alglib::real_1d_array &displacement)
{
    // The displacement is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
    // need to know where the "x" end and where the "y" begin.
    int nr_x_elements = displacement.length() / 2; // Shifts to y section

    Node *n; // Non border, inside node

    // We loop over all the nodes that are not on the border, ie. the innside nodes.
    for (size_t i = 0; i < mesh.freeNodeIds.size(); i++)
    {
        n = mesh[mesh.freeNodeIds[i]];
        n->x = n->init_x + displacement[i];
        n->y = n->init_y + displacement[i + nr_x_elements];
    }
}

// Our initial guess will be that all particles have shifted by the same
// transformation as the border.
void Simulation::m_initialGuess()
{
    // The displacement is structed like this: [x1,x2,x3,x4,y1,y2,y3,y4], so we
    // need to know where the "x" end and where the "y" begin.
    int nr_x_elements = nodeDisplacements.length() / 2; // Shifts to y section

    Node transformed_node; // These are temporary variables for readability
    const Node *n;         // Non border, inside node

    // We loop over all the nodes that are not on the border, ie. the innside nodes.
    for (size_t i = 0; i < mesh.freeNodeIds.size(); i++)
    {
        n = mesh[mesh.freeNodeIds[i]];
        transformed_node = transform(loadStepTransform, *n);
        // Subtract the initial possition to get the nodeDisplacements
        translateInPlace(transformed_node, n->init_x, n->init_y, -1.0); // node1.position - node2.position
        nodeDisplacements[i] = transformed_node.x;
        nodeDisplacements[i + nr_x_elements] = transformed_node.y;
    }
}

void addNoise(alglib::real_1d_array &displacement, double noise)
{
    int nr_x_elements = displacement.length() / 2; // Shifts to y section

    Node transformed_node; // These are temporary variables for readability
    const Node *n;         // Non border, inside node

    // We loop over all the nodes that are not on the border, ie. the innside nodes.
    for (size_t i = 0; i < nr_x_elements; i++)
    {
        // Generate random noise in the range [-noise, noise]
        double noise_x = ((double)rand() / RAND_MAX) * 2 * noise - noise;
        double noise_y = ((double)rand() / RAND_MAX) * 2 * noise - noise;

        // Add noise to the displacement
        displacement[i] += noise_x;
        displacement[i + nr_x_elements] += noise_y;
    }
}

Matrix2x2<double> getShear(double load, double theta)
{
    // perturb is currently unused. If it will be used, it should be implemeted
    // propperly.
    double perturb = 0;

    Matrix2x2<double> trans;
    trans[0][0] = (1. - load * cos(theta + perturb) * sin(theta + perturb));
    trans[1][1] = (1. + load * cos(theta - perturb) * sin(theta - perturb));
    trans[0][1] = load * pow(cos(theta), 2.);
    trans[1][0] = -load * pow(sin(theta - perturb), 2.);

    return trans;
}

void Simulation::m_updateProgress(double load)
{
    // If this is the first time we show progress, we hide the cursor
    if (load - startLoad == 0)
    {
        // Hide cursor. We show it again after calling m_exit
        indicators::show_console_cursor(false);
    }
    // Calculate progress
    double progress = (load - startLoad) / (maxLoad - startLoad) * 100;
    getBar().set_progress(progress);
}

void Simulation::m_writeToDisk(double load)
{
    static double lastLoadWritten = 0;
    // This is only used for logging purposes (no physics)
    mesh.load = load;

    spdlog::info("{}", load);

    // Only if there are lots of plastic events will we want to save the data.
    // If we save every frame, it requires too much storage.
    // (A 100x100 load from 0.15 to 1 with steps of 1e-5 would take up 180GB)
    if (mesh.nrPlasticEvents() > mesh.nrElements * plasticityEventThreshold)
    {
        writeMeshToVtu(mesh, name, dataPath);
        lastLoadWritten = load;
    }
    // If there are few large avalanvhes, we might go long without saving data
    // In order to get a good framerate for an animation, we want to ensure that
    // not too much happens between frames. This enures that we at least have
    // 1000 frames of states over the course of loading
    else if ((load - lastLoadWritten) / (maxLoad - startLoad) > 0.001)
    {
        writeMeshToVtu(mesh, name, dataPath);
        lastLoadWritten = load;
    }

    writeMeshToCsv(mesh, name, dataPath);
}

void Simulation::m_exit()
{
    // Completes the progress bar
    m_updateProgress(maxLoad);
    // This is to precent the final progress bar line from being overwritten
    std::cout << '\n';

    // This creates a pvd file that links all the utv files together.
    leanvtk::createCollection(getDataPath(name, dataPath),
                              getOutputPath(name, dataPath),
                              COLLECTIONNAME);

    std::string simulationTime = timer.CurrentTime();
    spdlog::info("Simulation time: {}", simulationTime);

    // Show cursor
    indicators::show_console_cursor(true);

    // Close and flush logger
    spdlog::drop(LOGNAME);
}

int get_terminal_width()
{
    struct winsize size;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &size);
    return size.ws_col; // Number of columns in the terminal
}

indicators::BlockProgressBar &getBar()
{
    int bar_width = 30; // Use a fixed width for the progress bar

    using namespace indicators;
    static BlockProgressBar bar{
        option::BarWidth{bar_width},
        option::Start{"["},
        option::End{"]"},
        option::PrefixText{"Simulation time "},
        // option::ForegroundColor{Color::yellow},
        option::ShowElapsedTime{true},
        option::ShowRemainingTime{true},
        // option::FontStyles{std::vector<FontStyle>{FontStyle::bold}},
    };
    return bar;
}

void printReport(const alglib::minlbfgsreport &report)
{
    // https://www.alglib.net/translator/man/manual.cpp.html#sub_minlbfgsoptimize
    std::cout << "Optimization Report:\n";
    std::cout << "\tIterations Count: " << report.iterationscount << '\n';
    std::cout << "\tNumber of Function Evaluations: " << report.nfev << '\n';
    std::cout << "\tTermination Reason: ";
    switch (report.terminationtype)
    {
    case -8:
        std::cout << "Infinite or NAN values in function/gradient";
        break;
    case -2:
        std::cout << "Rounding errors prevent further improvement";
        break;
    case -1:
        std::cout << "Incorrect parameters were specified";
        break;
    case 1:
        std::cout << "Relative function improvement is no more than EpsF";
        break;
    case 2:
        std::cout << "Relative step is no more than EpsX";
        break;
    case 4:
        std::cout << "Gradient norm is no more than EpsG";
        break;
    case 5:
        std::cout << "MaxIts steps was taken";
        break;
    case 7:
        std::cout << "Stopping conditions are too stringent, further improvement is impossible";
        break;
    case 8:
        std::cout << "Terminated by user request";
        break;
    default:
        std::cout << "Unknown termination reason";
    }
    std::cout << std::endl;
}

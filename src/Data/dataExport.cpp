#include "dataExport.h"

std::string findOutputPath()
{
    // Define the paths to check
    std::vector<std::string> paths = {
        "/Volumes/data/",
        "/media/elias/dataStorage/",
        "/data2/elundheim/",
        "/data/elundheim/"};

    // Initialize a variable to store the chosen path
    std::string chosen_path;

    // Iterate through the paths and check if they exist
    for (const auto &path : paths)
    {
        if (std::filesystem::exists(path))
        {
            chosen_path = path;
            break; // Stop the loop once a valid path is found
        }
    }

    // Check if a valid path was found or throw an error
    if (chosen_path.empty())
    {
        throw std::runtime_error("None of the provided paths exist.");
    }
    else
    {
        // We now also add the output folder name
        chosen_path += OUTPUTFOLDERPATH;
        std::cout << "Chosen output path: " << chosen_path << std::endl;
    }

    return chosen_path;
}

std::string getCurrentDate()
{
    auto now = std::chrono::system_clock::now();
    std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm = *std::localtime(&now_time_t);
    std::stringstream ss;
    ss << std::put_time(&now_tm, "%H-%M_%d.%m.%Y");
    return ss.str();
}

bool create_directory_if_not_exists(const std::filesystem::path &path)
{
    try
    {
        // std::filesystem::create_directories creates all intermediate directories
        // in the path if they do not exist and does nothing if they do.
        std::filesystem::create_directories(path);
    }
    catch (const std::filesystem::filesystem_error &e)
    {
        std::cerr << "Error creating directory '" << path << "': " << e.what() << std::endl;
        return false;
    }
    return true;
}

namespace fs = std::filesystem;

std::string getOutputPath(const std::string &name, const std::string &dataPath)
{
    fs::path outputPath = fs::path(dataPath) / name;
    return outputPath.string() + '/';
}

std::string getDataPath(const std::string &name, const std::string &dataPath)
{
    fs::path dataPathObj = fs::path(getOutputPath(name, dataPath)) / DATAFOLDERPATH;
    return dataPathObj.string() + '/';
}

std::string getFramePath(const std::string &name, const std::string &dataPath)
{
    fs::path framePath = fs::path(getOutputPath(name, dataPath)) / FRAMEFOLDERPATH;
    return framePath.string() + '/';
}

void createDataFolder(std::string name, std::string dataPath)
{
    std::vector<std::string> paths = {
        getDataPath(name, dataPath),
        getFramePath(name, dataPath)};
    for (std::string path : paths)
    {
        // Ensure the directory exists
        // This function creates all neccecary folders, so it works
        // even if the output folder, or the subfolder doesn't exists
        if (!create_directory_if_not_exists(path))
        {
            std::cerr << "Failed to create directory: " << path << std::endl;
        }
    }
}

// Get file path, and check that the path exsists
std::string getFilePath(std::string fileName,
                        std::string folderName,
                        std::string dataPath,
                        std::string fileType = ".vtu")
{
    std::string fileNameWithType = fileName + fileType;
    std::string directory = getDataPath(folderName, dataPath);

    // Check if the directory exists
    if (!std::filesystem::exists(directory) || !std::filesystem::is_directory(directory))
    {
        throw std::runtime_error("Directory does not exist: " + directory + "\nHave you run the createDataFolder function?");
    }

    return directory + fileNameWithType;
}

// Clears a subfolder. It only clears files of specified types for safety.
// If you want to delete the entire outputfolder, do it manually.
void clearOutputFolder(std::string name, std::string dataPath)
{

    std::vector<std::string> paths = {
        getOutputPath(name, dataPath),
        getDataPath(name, dataPath),
        getFramePath(name, dataPath)};
    // Define the list of file extensions to delete
    std::vector<std::string> extensionsToDelete = {
        ".vtu",
        ".pvd",
        ".csv",
        ".png",
        ".log"};
    for (std::string path : paths)
    {

        // Check if the directory exists
        if (!std::filesystem::exists(path) || !std::filesystem::is_directory(path))
        {
            continue;
        }
        // Iterate through each file in the directory
        for (const auto &entry : std::filesystem::directory_iterator(path))
        {
            if (entry.is_regular_file())
            {
                // Get the file extension
                std::string extension = entry.path().extension().string();
                // Check if the file's extension is in the list of extensions to delete
                if (std::find(extensionsToDelete.begin(), extensionsToDelete.end(), extension) != extensionsToDelete.end())
                {
                    // Delete the file
                    std::filesystem::remove(entry.path());
                }
            }
        }
    }
}

// If we want to store some data that does not depend on either the node or cell,
// it is inefficient to store the data multiple times. The simplest way I have
// found to store extra data is by including it in the file name, dataPath.
// Example: The variable foo and bar are stored as "_foo=0.32_bar=4_".
std::string makeFileName(const Mesh &mesh, std::string name, std::string dataPath)
{
    std::stringstream ss;
    ss << name << "_load=" << mesh.load << '_';
    return ss.str();
}

/*
    We will add the nodes in 4 passes, along each of the edges, first
    horizontal, then vertical

            4 5

        11  2 3  9
        10  0 1  8

            6 7

            6 7 8

        15  3 4 5  13
        14  0 1 2  12

            9 1011

    This choice of inxeding is natural because here low numbers map to low
    numbers, ie. 0 maps to 6, 5 maps to 15

*/
// Use the guide above to see that this works.
std::array<std::vector<int>, 2> getPBCMap(int n, int m)
{
    // Since the corners map to two places, we need two vectors
    int nm = n * m;
    std::vector<int> hMap(nm);
    std::vector<int> vMap(nm);

    // Set the standard map.
    for (size_t i = 0; i < nm; i++)
    {
        hMap[i] = i;
        vMap[i] = i;
    }

    // Here we should modify the values so that the border atoms are associated
    // with a different index.
    // Horizontal borders
    for (size_t i = 0; i < m; i++)
    {
        vMap[i] = nm + i;
        vMap[nm - m + i] = nm + m + i;
    }
    // Vertical borders
    for (size_t i = 0; i < n; i++)
    {
        hMap[i * m] = nm + 2 * m + i;
        hMap[(i + 1) * m - 1] = nm + 2 * m + n + i;
    }
    return {hMap, vMap};
}

void writeMeshToVtu(Mesh &mesh, std::string folderName, std::string dataPath, bool automaticNumbering)
{
    const int dim = 3;
    const int cell_size = 3;
    int n = mesh.nodes.rows;
    int m = mesh.nodes.cols;
    int nm = n * m;
    static int timeStep = 0;
    int nrRealNodes = nm;
    int nrElements = mesh.nrElements;

    std::string fileName = makeFileName(mesh, folderName, dataPath);

    std::string filePath;
    if (automaticNumbering)
    {
        filePath = getFilePath(fileName + "." + std::to_string(timeStep), folderName, dataPath);
        timeStep += 1;
    }
    else
    {
        filePath = getFilePath(fileName, folderName, dataPath);
    }

    /*
        We have N = nm real nodes, but to display the elements of the mesh that
        should map around the system when using periodic boundary conditions, we
        need to create additional phantom nodes that are duplicates of real nodes,
        but flipped to the other side of the system.
        There will be 2n+2m extra nodes to cover the boarder of the system.
    */
    int nrPhantomNodes = (mesh.usingPBC ? 2 * n + 2 * m : 0);
    int nrNodes = nrRealNodes + nrPhantomNodes;

    std::vector<double> points(nrNodes * dim, 0);
    std::vector<double> force(nrNodes * dim, 0);
    std::vector<int> elements(nrElements * cell_size);
    std::vector<double> energy(nrElements);
    std::vector<double> C11(nrElements);
    std::vector<double> C12(nrElements);
    std::vector<double> C22(nrElements);
    std::vector<double> P11(nrElements);
    std::vector<double> P12(nrElements);
    std::vector<double> P21(nrElements);
    std::vector<double> P22(nrElements);
    std::vector<double> resolvedShearStress(nrElements);

    leanvtk::VTUWriter writer;

    // define a helper function to move data from mesh to vectors with displacement
    auto meshToVectors = [&mesh, &points, &force](int i, int j, double x = 0, double y = 0)
    {
        points[i * 3 + 0] = mesh.nodes.data[j].X() + x;
        points[i * 3 + 1] = mesh.nodes.data[j].Y() + y;
        points[i * 3 + 2] = 0;
        force[i * 3 + 0] = mesh.nodes.data[j].f_x;
        force[i * 3 + 1] = mesh.nodes.data[j].f_y;
        force[i * 3 + 2] = 0;
    };

    // This mapping process is not easy to understand, but it works.
    // We have now moved all the real node data, now we need to duplicate the
    // border nodes and later make sure that the mesh elements are connected to
    // the duplicates instead of the original nodes (which would make the element
    // span the system, which looks bad).
    static std::array<std::vector<int>, 2> PBCmap = getPBCMap(n, m);

    for (int i = 0; i < nrRealNodes; ++i)
    {
        meshToVectors(i, i);

        // If we are also at the edge of the system and using PBC, we also update
        // the coresponding phantom node
        if (mesh.usingPBC)
        {
            int j = PBCmap[0][i];
            // Bottom edge of the system
            if (i / m == 0)
            {
                std::cout << i << ", " << j << '\n';
                meshToVectors(PBCmap[1][i], i, 0, n * mesh.a);
            }

            // Top edge of the system
            if (i / m == n - 1)
            {
                std::cout << i << ", " << j << '\n';
                meshToVectors(PBCmap[1][i], i, 0, -n * mesh.a);
            }

            // Left edge of the system
            if (i % m == 0)
            {
                meshToVectors(PBCmap[0][i], i, m * mesh.a, 0);
            }

            // Right edge of the system
            if (i % m == m - 1)
            {
                meshToVectors(PBCmap[0][i], i, -m * mesh.a, 0);
            }
        }
    }

    for (int i = 0; i < nrElements; ++i)
    {
        TElement &e = mesh.elements[i];
        if (e.xOffset != 0 || e.yOffset != 0)
        {
            // We do some conditional movements of the nodes
            // If there is a y offset, we move either the 0 or the 2 node.
            if (e.yOffset > 0)
            {
                elements[i * 3 + 2] = PBCmap[1][e.nodes[2]->id.i];
                elements[i * 3 + 0] = e.nodes[0]->id.i;
            }
            else if (e.yOffset < 0)
            {
                elements[i * 3 + 0] = PBCmap[1][e.nodes[0]->id.i];
                elements[i * 3 + 2] = e.nodes[2]->id.i;
            }
            else
            {
                elements[i * 3 + 0] = e.nodes[0]->id.i;
                elements[i * 3 + 2] = e.nodes[2]->id.i;
            }

            // If there is an x offset, we always move the node with index 1
            if (e.xOffset != 0)
            {
                elements[i * 3 + 1] = PBCmap[0][e.nodes[1]->id.i];
            }
            else
            {
                elements[i * 3 + 1] = e.nodes[1]->id.i;
            }
        }
        else
        {
            elements[i * 3 + 0] = e.nodes[0]->id.i;
            elements[i * 3 + 1] = e.nodes[1]->id.i;
            elements[i * 3 + 2] = e.nodes[2]->id.i;
        }
        energy[i] = e.energy;
        C11[i] = e.C[0][0];
        C12[i] = e.C[0][1];
        C22[i] = e.C[1][1];
        P11[i] = e.P[0][0];
        P12[i] = e.P[0][1];
        P21[i] = e.P[1][0];
        P22[i] = e.P[1][1];
        resolvedShearStress[i] = e.resolvedShearStress;
    }

    // connect data to writer
    writer.add_cell_scalar_field("energy_field", energy);
    writer.add_cell_scalar_field("C11", C11);
    writer.add_cell_scalar_field("C12", C12);
    writer.add_cell_scalar_field("C22", C22);
    writer.add_cell_scalar_field("P11", P11);
    writer.add_cell_scalar_field("P12", P12);
    writer.add_cell_scalar_field("P21", P21);
    writer.add_cell_scalar_field("P22", P22);
    writer.add_cell_scalar_field("resolvedShearStress", resolvedShearStress);

    writer.add_vector_field("stress_field", force, dim);

    // write data
    writer.write_surface_mesh(filePath, dim, cell_size, points, elements);
}

// Create a spdlog logger
std::shared_ptr<spdlog::logger> createLogger(const std::string &folderName, const std::string &dataPath)
{
    // Construct file path
    std::string filePath = getOutputPath(folderName, dataPath) + MACRODATANAME + ".csv";
    // Create a file logger (basic file sink)
    auto file_logger = spdlog::basic_logger_mt("basic_logger", filePath, true); // True for appending mode
    file_logger->set_pattern("%v");                                             // Set pattern to raw message
    return file_logger;
}

void writeLineToCsv(const std::vector<std::string> &strings, const std::string &folderName, const std::string &dataPath)
{
    // Create logger or get existing one
    static auto logger = createLogger(folderName, dataPath);
    std::string line = spdlog::fmt_lib::format("{}", spdlog::fmt_lib::join(strings, ","));
    logger->info(line); // Write line to file via spdlog
}

void writeLineToCsv(const std::vector<double> &values, const std::string &folderName, const std::string &dataPath)
{
    std::vector<std::string> stringValues;
    stringValues.reserve(values.size());
    for (double value : values)
    {
        stringValues.push_back(std::to_string(value));
    }
    writeLineToCsv(stringValues, folderName, dataPath);
}

void writeMeshToCsv(Mesh &mesh, const std::string &folderName, const std::string &dataPath, bool isFirstLine)
{
    static int lineCount = 0;
    std::vector<std::string> lineData;

    if (isFirstLine)
    {
        lineCount = 0; // Reset line count if it's the first line
        lineData = {"Line nr", "Load", "Avg. energy", "Avg. RSS"};
    }
    else
    {
        lineCount += 1;
        lineData = {
            std::to_string(lineCount),
            std::to_string(mesh.load),
            std::to_string(mesh.averageEnergy),
            std::to_string(mesh.averageResolvedShearStress())};
    }
    writeLineToCsv(lineData, folderName, dataPath);
}
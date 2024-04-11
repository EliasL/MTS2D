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
    ss << name
       << "_load=" << mesh.load
       << "_nrM=" << mesh.nrPlasticEvents()
       << '_';
    return ss.str();
}

void writeMeshToVtu(const Mesh &mesh, std::string folderName, std::string dataPath)
{

    const int dim = 3;
    const int cell_size = 3;
    int n = mesh.rows;
    int m = mesh.cols;
    if (mesh.usingPBC)
    {
        // We create an extra row and column for the ghost nodes.
        n += 1;
        m += 1;
    }
    int nm = n * m;
    // Since timeStep is static, it will increase each time the function is called.
    static int timeStep = 0;
    int nrNodes = nm;
    int nrElements = mesh.nrElements;

    std::string fileName = makeFileName(mesh, folderName, dataPath);

    std::string filePath;

    filePath = getFilePath(fileName + "." + std::to_string(timeStep), folderName, dataPath);
    timeStep += 1;

    std::vector<double> points(nrNodes * dim);
    std::vector<double> force(nrNodes * dim);
    std::vector<double> fixed(nrNodes); // boolean values represented by 0.0 and 1.0
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

    // Instead of getting the data directly from the nodes in the mesh, we extract
    // the data from the nodes in the elements in the mesh. This is because they
    // have a displaced position and to not result in overlapping elements
    std::vector<bool> alreadyCopied(mesh.nrNodes);
    // Iterate over each element in the mesh
    for (size_t elementIndex = 0; elementIndex < nrElements; ++elementIndex)
    {
        const TElement &e = mesh.elements[elementIndex];
        // Iterate over each node in the element
        for (size_t j = 0; j < e.nodes.size(); ++j)
        {
            const Node &n = e.nodes[j];
            int nodeIndex = mesh.usingPBC ? n.ghostId.i : n.id.i;
            if (!alreadyCopied[nodeIndex])
            {
                points[nodeIndex * dim + 0] = n.pos()[0];
                points[nodeIndex * dim + 1] = n.pos()[1];
                points[nodeIndex * dim + 2] = 0;
                force[nodeIndex * dim + 0] = n.f[0];
                force[nodeIndex * dim + 1] = n.f[1];
                force[nodeIndex * dim + 2] = 0;
                fixed[nodeIndex] = n.fixedNode;
                alreadyCopied[nodeIndex] = true;
            }

            // We choose to either use the ghost id or the real id depending on
            // whether or not we are using pbc.
            elements[elementIndex * cell_size + j] = nodeIndex;
        }

        energy[elementIndex] = e.energy;
        C11[elementIndex] = e.C[0][0];
        C12[elementIndex] = e.C[0][1];
        C22[elementIndex] = e.C[1][1];
        P11[elementIndex] = e.P[0][0];
        P12[elementIndex] = e.P[0][1];
        P21[elementIndex] = e.P[1][0];
        P22[elementIndex] = e.P[1][1];
        resolvedShearStress[elementIndex] = e.resolvedShearStress;
    }

    // connect data to writer
    writer.add_cell_scalar_field("energy_field", energy);
    writer.add_cell_scalar_field("resolvedShearStress", resolvedShearStress);
    writer.add_scalar_field("fixed", fixed);
    writer.add_cell_scalar_field("C11", C11);
    writer.add_cell_scalar_field("C12", C12);
    writer.add_cell_scalar_field("C22", C22);
    writer.add_cell_scalar_field("P11", P11);
    writer.add_cell_scalar_field("P12", P12);
    writer.add_cell_scalar_field("P21", P21);
    writer.add_cell_scalar_field("P22", P22);

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

void writeMeshToCsv(Mesh &mesh, const std::string &folderName, const std::string &dataPath)
{
    static int lineCount = 0;
    lineCount += 1;
    std::vector<std::string> lineData = {
        std::to_string(lineCount),
        std::to_string(mesh.load),
        std::to_string(mesh.averageEnergy),
        std::to_string(mesh.maxEnergy),
        std::to_string(mesh.averageResolvedShearStress())};

    writeLineToCsv(lineData, folderName, dataPath);
}

void writeCsvCols(const std::string &folderName, const std::string &dataPath)
{
    static int lineCount = 0;
    std::vector<std::string> lineData;
    lineData = {"Line nr", "Load", "Avg. energy", "Max. energy", "Avg. RSS"};
    writeLineToCsv(lineData, folderName, dataPath);
}
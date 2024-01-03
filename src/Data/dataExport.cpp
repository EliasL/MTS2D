#include "dataExport.h"


bool create_directory_if_not_exists(const std::string &path)
{
    std::vector<std::string> dirs;
    std::stringstream ss(path);
    std::string item;
    while (getline(ss, item, '/'))
    {
        if (!item.empty())
        {
            dirs.push_back(item);
        }
    }

    std::string current_path;
    for (const auto &dir : dirs)
    {
        current_path += dir + "/";
#if defined(_WIN32)
        // Windows does not have a built-in function for recursive directory creation
        // Here you might want to implement a loop that creates each directory in the path
        if (_mkdir(current_path.c_str()) != 0 && errno != EEXIST)
        {
            std::cerr << "Error creating directory '" << current_path << "': " << strerror(errno) << std::endl;
            return false;
        }
#else
        struct stat info;
        if (stat(current_path.c_str(), &info) != 0 || !(info.st_mode & S_IFDIR))
        {
            if (mkdir(current_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0 && errno != EEXIST)
            {
                std::cerr << "Error creating directory '" << current_path << "': " << strerror(errno) << std::endl;
                return false;
            }
        }
#endif
    }

    return true;
}

void createDataFolder(std::string subFolder){
    std::vector<std::string> paths = { 
        OUTPUTFOLDERPATH + subFolder + DATAFOLDERPATH,
        OUTPUTFOLDERPATH + subFolder + FRAMEFOLDERPATH,
    };
    for( std::string path : paths){
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
                        std::string fileType = ".vtu",
                        std::string subFolder = DEFAULTSUBFOLDER)
{
    std::string fileNameWithType = fileName + fileType;
    std::string directory = OUTPUTFOLDERPATH + subFolder + DATAFOLDERPATH;

    // Check if the directory exists
    if (!std::filesystem::exists(directory) || !std::filesystem::is_directory(directory)) {
        throw std::runtime_error("Directory does not exist: " + directory 
                                + "\nHave you run the createDataFolder function?");
    }

    return directory + fileNameWithType;
}

// Clears a subfolder. It only clears files of specified types for safety.
// If you want to delete the entire outputfolder, do it manually.
void clearOutputFolder(std::string folder) {
    std::vector<std::string> paths = { 
        OUTPUTFOLDERPATH + folder,
        OUTPUTFOLDERPATH + folder + DATAFOLDERPATH,
        OUTPUTFOLDERPATH + folder + FRAMEFOLDERPATH,
    };
    // Define the list of file extensions to delete
    std::vector<std::string> extensionsToDelete = {
        ".vtu", 
        ".pvd",
        ".png",
        ".mp4",
        ".csv"
    };
    for (std::string path : paths){

        // Check if the directory exists
        if (!std::filesystem::exists(path) || !std::filesystem::is_directory(path)) {
            continue;
        }
        // Iterate through each file in the directory
        for (const auto& entry : std::filesystem::directory_iterator(path)) {
            if (entry.is_regular_file()) {
                // Get the file extension
                std::string extension = entry.path().extension().string();
                // Check if the file's extension is in the list of extensions to delete
                if (std::find(extensionsToDelete.begin(), extensionsToDelete.end(), extension) != extensionsToDelete.end()) {
                    // Delete the file
                    std::filesystem::remove(entry.path());
                }
            }
        }
    }

}

// If we want to store some data that does not depend on either the node or cell,
// it is inefficient to store the data multiple times. The simplest way I have
// found to store extra data is by including it in the file name. 
// Example: The variables foo and bar are stored as "_foo=0.32_bar=4_".
std::string makeName(const Mesh &mesh, const std::string &fileName){
    return fileName + "_load=" +std::to_string(mesh.load)+"_";
}

void writeToVtu(Mesh &mesh, std::string fileName, bool automaticNumbering)
{
    const int dim = 3;
    const int cell_size = 3;
    static int timeStep = 0;
    int nrNodes = mesh.nodes.data.size();
    int nrElements = mesh.nrElements;
    
    fileName = makeName(mesh, fileName);

    std::string filePath;
    if (automaticNumbering)
    {
        filePath = getFilePath(fileName + "." + std::to_string(timeStep));
        timeStep += 1;
    }
    else
    {
        filePath = getFilePath(fileName);
    }

    std::vector<double> points(nrNodes * dim, 0);
    std::vector<double> force(nrNodes * dim, 0);
    std::vector<int> elements(nrElements * cell_size);
    std::vector<double> energy(nrElements);

    leanvtk::VTUWriter writer;

    // Transfer data
    for (int i = 0; i < nrNodes; ++i)
    {
        points[i * 3 + 0] = mesh.nodes.data[i].x;
        points[i * 3 + 1] = mesh.nodes.data[i].y;
        points[i * 3 + 2] = 0;
        force[i * 3 + 0] = mesh.nodes.data[i].f_x;
        force[i * 3 + 1] = mesh.nodes.data[i].f_y;
        force[i * 3 + 2] = 0;
    }

    for (int i = 0; i < nrElements; ++i)
    {
        elements[i * 3 + 0] = mesh.elements[i].n1->id.i;
        elements[i * 3 + 1] = mesh.elements[i].n2->id.i;
        elements[i * 3 + 2] = mesh.elements[i].n3->id.i;
        energy[i] = mesh.elements[i].energy;
    }

    // connect data to writer
    writer.add_cell_scalar_field("energy_field", energy);
    writer.add_vector_field("stress_field", force, dim);

    // write data
    writer.write_surface_mesh(filePath, dim, cell_size, points, elements);
}
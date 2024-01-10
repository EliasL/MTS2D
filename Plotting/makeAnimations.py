import subprocess
from icecream import ic
import os

from settings import settings
from vtkFunctions import *

# Use ffmpeg to convert a folder of .png images into a mp4 file
def makeAnimations(path, pvd_file):
   
    ic("Creating frames...")

    dataPath = path + settings["DATAFOLDERPATH"]   
    framePath = path + settings["FRAMEFOLDERPATH"]  
    if(not os.path.exists(dataPath+pvd_file)):
        print(f"No file found at: {dataPath+pvd_file}")
        return
    
    vtu_files = parse_pvd_file(dataPath, pvd_file)
    nrSteps, nrNodes, nrElements = getDataSize(dataPath, vtu_files)
    

    makeImages(framePath, dataPath, vtu_files)
    
    ic("Creating animations...")
    # Length of video in seconds
    videoLength = 10
    # Define the frame rate
    fps = nrSteps/videoLength

    # Define the output video filename
    # The name of the video is the same as the name of the folder+_video.mp4
    outputVideo = path+path.split('/')[-2]+'_video.mp4'
    framesPath = path+settings["FRAMEFOLDERPATH"]

    # Construct the FFmpeg command as a string
    ffmpeg_command = f"ffmpeg -y -r {fps} -pattern_type glob -i '{framesPath}/*.png' -c:v libx264 -pix_fmt yuv420p {outputVideo}"

    # Run the FFmpeg command using subprocess, and hide output
    result = subprocess.run(ffmpeg_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Check if the command was successful
    if result.returncode == 0:
        pass
    else:
        # If you see this, remove the stdout and stderr keywords from the subprocess
        # for more details on the error. (Have you installed ffmpeg?)
        print("ffmpeg command failed with return code", result.returncode)

if __name__ == "__main__":
    # Replace 'your_pvd_file.pvd' with the path to your .pvd file
    makeAnimations('build/output/testing/','collection.pvd')
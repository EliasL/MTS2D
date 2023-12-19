import subprocess

from settings import settings
from vtkFunctions import *

# Use ffmpeg to convert a folder of .png images into a mp4 file
def main(path, pvd_file):
   
    dataPath = path + settings["DATAFOLDERPATH"]   
    framePath = path + settings["FRAMEFOLDERPATH"]  
    vtu_files = parse_pvd_file(dataPath, pvd_file)
    nrSteps, nrNodes, nrElements = getDataSize(dataPath, vtu_files)
    

    makePngImages(framePath, dataPath, vtu_files)
    
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
    subprocess.run(ffmpeg_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

print("Creating animation...")
# Replace 'your_pvd_file.pvd' with the path to your .pvd file
main('build/output/testing/','collection.pvd')
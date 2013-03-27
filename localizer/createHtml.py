import os
import sys

def visDir(imgDir, titleName):

    allFileNames = os.listdir(imgDir)
    htmlFileName = imgDir + 'sampleViewer.html'
    imgDir = imgDir + '/'
    fo= open(htmlFileName, 'w')

    fo.write("<!DOCTYPE html> \n")
    fo.write("<html> \n")
    fo.write("<body> \n")
    fo.write("<h2>" + titleName + "</h2> \n")

    for fileName in allFileNames:
        currFilePath =  imgDir + fileName
        currImgCmd = "<img border=\"0\" src=\""+ currFilePath + "\" alt=\"Pulpit rock\" width=\"304\" height=\"228\"> \n"
        fo.write(currImgCmd)

    fo.write("<body> \n")
    fo.write("<html> \n")

visDir(sys.argv[1], sys.argv[2])
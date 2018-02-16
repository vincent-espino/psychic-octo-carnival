# CSCI 576 - Multimedia Systems Design Assignment 1

Per the assignment description, this program can be run from the command line using the following command (make sure to compile with ‘javac imageReader.java’ first!):

imageReader.java path-to-image Y U V Q

where path-to-image can be a relative or absolute file path to a .rgb file and Y, U, V, and Q are the subsampling in the Y, U, V channels respectively. Q is the subsampling rate and can range from 0 to 255 (256 can also be used, which implies no subsampling). If less than 5 command line arguments are given, the program will give a message saying that not enough arguments were given, and it will quit.

As noted by Chloe Legendre, assumptions made for the program/written/analysis parts should be listed in this README file. The assumptions I made were as follows:

• For the programming part, all images are of a fixed width and height of 352 and 288 pixels, respectively.
• If Q is 0, I perform subsampling as if Q were equal to 1 (which results in a completely blank/black image), since Q = 0 implies no colors.


#include "progressbar.hpp"
#include <iostream>

progress_bar::progress_bar(int barlength, double start_progress, double end_progress) {
    currentProgress = start_progress;
    neededProgress = end_progress;
    pBarLength = barlength;
}

// destructor
progress_bar::~progress_bar(){};

// set the length of progress bar 
void progress_bar::setlength(int desired_length) {
   neededProgress = desired_length;
}
  
 // update progress by 
 void progress_bar::update(double newProgress) {
  currentProgress += newProgress;
  amountOfFiller = (int)((currentProgress / neededProgress)*(double)pBarLength);
}

// print the current progress of the bar
void progress_bar::print() {
  currUpdateVal %= pBarUpdater.length();
   std::cout << "\r" //Bring cursor to start of line
             << firstPartOfpBar; //Print out first part of pBar
   for (int a = 0; a < amountOfFiller; a++) { //Print out current progress
     std::cout << pBarFiller;
   }
   std::cout << pBarUpdater[currUpdateVal];
   for (int b = 0; b < pBarLength - amountOfFiller; b++) { //Print out spaces
     std::cout << " ";
   }
   std::cout << lastPartOfpBar //Print out last part of progress bar
             << " (" << (int)(100*(currentProgress/neededProgress)) << "%)" //This just prints out the percent
   	         << std::flush;
    currUpdateVal += 1;
}

#ifndef PROGRESSBAR_HPP 
#define PROGRESSBAR_HPP

#include <iostream>
#include <string>

/*
 * Class for a generic progress bar indicator, you can set the desired width of 
 * the bar, the start progress for the bar and the total number of entries
 * expected in the bar
 * 
 * This is a modiifed version of Flare Cats answer from 
 * https://stackoverflow.com/questions/14539867/how-to-display-a-progress-
 * indicator-in-pure-c-c-cout-printf
 * 
 */

class progress_bar {
  public:
    // specific constructor;
    progress_bar(int barlength = 50, double start_progress = 0., double end_progress = 100.);
    // destructor
   ~progress_bar();

    // set the length of progress bar 
    void setlength(int desired_length);
    // update progress by 
    void update(double newProgress);
    // print the current progress of the bar
    void print();
  
  public:
    std::string firstPartOfpBar = "["; 
    std::string lastPartOfpBar = "]";
    std::string pBarFiller = "|";
    std::string pBarUpdater = "/-\\|";

  private:
    int amountOfFiller;
    int pBarLength;
    int currUpdateVal; 
    double currentProgress, neededProgress;
};

#endif
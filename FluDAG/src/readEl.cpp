#include <iostream>  // cout, endl
#include <fstream>   // ifstream
#include <cstring>
#include <vector>
#include <stdlib.h>  // atoi

const int MAX_CHARS_PER_LINE = 512;
const int MAX_TOKENS_PER_LINE = 20;
const char* const DELIMITER = " ";
const bool debug = true;

// The line in the file match this struct
struct FlukaElement 
{
   std::string fluka_name;
   unsigned int atomic_number;
   std::vector<int> id;
};

// Return a vector of structs, one for each line in the file
std::vector<FlukaElement*> readElements(std::ifstream& fin);

// Hardcoded to read a file of fluka-named elements, atomic #, 
// each with 12 or more id's
std::vector<FlukaElement*> readElements(std::ifstream& fin)
{
   std::string lastString = "notAnElement";
   int lastNum = 0;
   std::vector<FlukaElement*> fluka_elements;
   FlukaElement* element = NULL;
   int i = 0;

   // read each line of the file
   while (!fin.eof())
   {
     // read an entire line into memory
     char buf[MAX_CHARS_PER_LINE];
     fin.getline(buf, MAX_CHARS_PER_LINE);

     // parse the line into blank-delimited tokens
     int n = 0;
     const char* token;

     // parse the line for the first word
     token = strtok(buf, DELIMITER);
     
     if (token) // line not blank
     {
       // Get the second token here because we'll need the atomic 
       // # either for assignment or for checking
       unsigned int atomic_num = atoi(strtok(0,DELIMITER));
       // And get the first (possibly only) id
       int id = atoi(strtok(0,DELIMITER));

       if (debug)
       {
          std::cout << "\n" << i++ << ". ";
       }
       // Check that the first token of the line is a new element
       if (lastString.compare(std::string(token)) != 0) 
       {
	  // If we're starting a new element, it's time to push back 
	  // the last one, unless we're at the very beginning
          if (element != NULL)
	  {
             fluka_elements.push_back(element);
	     std::cout << "New Name, id: " << element->fluka_name << ", " << element->atomic_number;
	  }
          element = new FlukaElement();
          element->fluka_name = std::string(token);
	  element->atomic_number = atomic_num;
	  element->id.push_back(id);

	  lastString = std::string(token);
	  lastNum    = atomic_num;
       }  
       else  // don't create a new element
       {
          // not new element, check for correct atomic #
          if (lastNum != element->atomic_number)
	  {
	     std::cout << "Error:  lastNum != new atomic # when it should" 
	               << std::endl;
	     exit(EXIT_FAILURE);
	  }
	  // If all good, just add the id to the vector
	  element->id.push_back(id);
       }
     } // end if line not blank
   }   // end while
   
   // Recall elements got added near the beginning of the while loop; 
   // The last element must be added after exiting the loop
   if (element != NULL)
   {
      std::cout << "\n" << i++ << ". ";
      std::cout << "New Name, id: " << element->fluka_name << ", " << element->atomic_number;
      fluka_elements.push_back(element);
   }
   std::cout << "\n============================================" << std::endl;
   return fluka_elements;
}

// Open a filestream for "el2.txt" and use its contents to populate
// a vector of FlukaElement struct pointers
int main()
{
   // create a file-reading object
   std::ifstream fin;
   fin.open("el2.txt");
   if (!fin.good())
   {
     // exit if file not found
     return EXIT_FAILURE;
   }
   std::vector<FlukaElement*> fe = readElements(fin);

   // Read the new data object and print out the results
   std::cout << std::endl;
   std::vector<FlukaElement*>::iterator iter;
   for (iter = fe.begin(); iter != fe.end();  ++iter) 
   {
       std::cout << (*iter)->fluka_name << ", ";
       std::cout << (*iter)->atomic_number << ", ";

       std::vector<int> id_vector = (*iter)->id; 
       if (id_vector.size() == 0)
       {
          std::cout << "Error: no ids. " << std::endl;
          exit(EXIT_FAILURE);
       }
       std::cout << "(" << id_vector.at(0);
       for (unsigned int i = 1; i < id_vector.size();  ++i)
       {
           std::cout << ", " << id_vector.at(i);
       }
       std::cout << ")" << std::endl;   
   }
   std::cout << std::endl;
}


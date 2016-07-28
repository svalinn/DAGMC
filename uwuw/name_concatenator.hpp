#include <iostream>
#include <string>
#include <set>

class name_concatenator
{
 public:
  /// constructor
  name_concatenator();

  /// destructor
  ~name_concatenator();

  /// returns the a unique name
  std::string make_name_8bytes(std::string name);

 private:
  /// Private functions
  /**
   * \Brief given the string name, increment the integer ID or shift the string should there be white
   * space remaining in the string.  Turns the name into a unique 8 character ID string, mainly for fluka
   * but has other utility in other codes it turns "Special Steel" into "SPECIALS" just by truncating and
   * capitalising however, on the second occurence of this it should turn the string into SPECIAL1,
   * then SPECIAL2, all until SPECI999. But, if the string is less than 8 chars anyway we should use
   * the space first, i.e. "Air, Dry N" should become "AIRDRYN ", then "AIRDRYN1", then "AIRDRY99" etc
   *
   * \param[in] name, The string to make a unique name for
   *
   * \return the shifted and incremented string
   */
  std::string shift_and_increment(std::string name);

  /**
   * \Brief  given the string name, extract the alpha number characters only
   *
   * \param[in] name, The string to extract alpha numerics from
   *
   * \return the string containing only alpha numerics
   */
  std::string extract_alpha_num(std::string name); ///< extracts the chars [A-Z] and numbers [0-9] only

  void int_to_string(int convert, std::string &string);

  /// Private variables
 private:
  std::set<std::string> used_b8_names; ///< the collection of names used so far
};


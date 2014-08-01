#ifndef DAGMC_FLUKA_FLUKA_NUCID_HPP
#define DAGMC_FLUKA_FLUKA_NUCID_HPP

struct Element
{
   std::string fluka_name;
   int znum;
   int snum;
};
// Issues:
// Yt znum is 39, not 38, per FM
// Where is znum 34= Se?
// Where is znum 37= Rb?
// Where is znum 44= Ru?
// Where is znum 45= Rh?
// Where is znum 46= Pd?
// Where is znum 53= Te?
// Where is znum 59= Pr? (Lanthanide Series)
// Where is znum 61= Pm? (Lanthanide Series)
//  likewise znum = 66-71, 76, 81, 84-89 and more

Element FlukaZNum[] = 
{
{"BERYLLIU",    4,   9},
{"CARBON",      6,   -2},
{"NITROGEN",    7,   -2},
{"OXYGEN",      8,   16},
{"FLUORINE",    9,   19},
{"NEON",       10,   -2},
{"SODIUM",     11,   23},
{"MAGNESIU",   12,   -2},
{"ALUMINUM",   13,   27},
{"SILICON",    14,   -2},
{"PHOSPHO",    15,   31},
{"SULFUR",     16,   -2},
{"CHLORINE",   17,   -2},
{"ARGON",      18,   -2},
{"POTASSIU",   19,   -2},
{"CALCIUM",    20,   -2},
{"SCANDIUM",   21,   45},
{"TITANIUM",   22,   -2},
{"VANADIUM",  23,   -2},
{"CHROMIUM",  24,   -2},
{"MANGANES",  25,   55},
{"IRON",      26,   -2},
{"COBALT",    27,   59},
{"NICKEL",    28,   -2},
{"COPPER",    29,   -2},
{"ZINC",      30,   -2},
{"GALLIUM",   31,   -2},
{"GERMANIU",  32,   -2},
{"ARSENIC",   33,   75},
{"BROMINE",   35,   -2},
{"KRYPTON",   36,   -2},
{"YTTRIUM",   38,   90},
{"ZIRCONIU",  40,   -2},
{"NIOBIUM",   41,   93},
{"MOLYBDEN",  42,   -2},
{"99-TC",     43,   99},
{"SILVER",    47,   -2},
{"CADMIUM",   48,   -2},
{"INDIUM",    49,   -2},
{"TIN",       50,   -2},
{"ANTIMONY",  51,   -2},
{"XENON",     54,   -2},
{"CESIUM",    55,   133},
{"BARIUM",    56,   -2},
{"LANTHANU",  57,   -2},
{"CERIUM",    58,   -2},
{"NEODYMIU",  60,   -2},
{"SAMARIUM",  62,   -2},
{"EUROPIUM",  63,   -2},
{"GADOLINI",  64,   -2},
{"TERBIUM",   65,   159},
{"HAFNIUM",   72,   -2},
{"TANTALUM",  73,   181},
{"TUNGSTEN",  74,   -2},
{"RHENIUM",   75,   -2},
{"IRIDIUM",   77,   -2},
{"PLATINUM",  78,   -2},
{"GOLD",      79,   197},
{"MERCURY",   80,   -2},
{"LEAD",      82,   -2},
{"BISMUTH",   83,   209},
{"239-PU",    94,   239},
{"241-AM",    95,   241}
};

// Put {" at the beginning of every line
// :s/^/\{\"/
// Add ", after first word
// :s/\(\s\+\)\s/\1",/
// Add comma after digits in pattern space-one-or-more digits-space
// :s/\(\s\d\+\)\s/\1,/
// Add }, to endofline
// s/$/\},/
Element FlukaNamedIsotopes[] = 
{
{"HYDROGEN",    1,  -2},
{"HYDROG-1",    1,   1},
{"DEUTERIU",    1,   2},
{"TRITIUM",     1,   4},
{"HELIUM",      2,  -2},
{"HELIUM-3",    2,   3},
{"HELIUM-4",    2,   4},
{"LITHIUM",     3,  -2},
{"LITHIU-6",    3,   6},
{"LITHIU-7",    3,   7},
{"BORON",       5,   2},
{"BORON-10",   5,   10},
{"BORON-11",   5,   11},
{"STRONTIU",  38,   -2},
{"90-SR",     38,   90},
{"IODINE",    53,   127},
{"129-I",     53,   129},
{"124-XE",    54,   124},
{"126-XE",    54,   126},
{"128-XE",    54,   128},
{"129-XE",    54,   129},
{"130-XE",    54,   130},
{"131-XE",    54,   131},
{"132-XE",    54,   132},
{"134-XE",    54,   134},
{"135-XE",    54,   135},
{"136-XE",    54,   136},
{"135-CS",    55,   135},
{"137-CS",    55,   137},
{"230-TH",    90,   230},
{"232-TH",    90,   232},
{"233-U",     92,   233},
{"234-U",     92,   234},
{"235-U",     92,   235},
{"238-U",     92,   238}
};

#endif // FLUDAG_SRC_FLUKA_FUNCS_H

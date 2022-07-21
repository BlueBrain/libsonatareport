// FixedPositions // Don't change so we can have backward compatibility

// If identifier read at position 0 matches ARCHITECTURE_IDENTIFIER, then the file was writting from
// native architecture
#define ARCHITECTURE_IDENTIFIER 1.001

// Double:
#define IDENTIFIER_POSITION 0
// Int:
#define HEADER_SIZE_POSITION IDENTIFIER_POSITION + sizeof(double)
// String:
#define LIBRARY_VERSION_POSITION 16
// String:
#define SIMULATOR_VERSION_POSITION 32
// Int:
#define TOTAL_NUMBER_OF_CELLS_POSITION 48
// Int:
#define TOTAL_NUMBER_OF_COMPARTMENTS_POSITION 52
// Int:
#define NUMBER_OF_STEPS_POSITION 64
// Double:
#define TIME_START_POSITION 72
// Double:
#define TIME_END_POSITION 80
// Double:
#define DT_TIME_POSITION 88
// String:
#define D_UNIT_POSITION 96
// String:
#define T_UNIT_POSITION 112
// Int:
#define MAPPING_SIZE_POSITION 128
// String:
#define MAPPING_NAME_POSITION 144
// Int:
#define EXTRA_MAPPING_SIZE_POSITION 160
// String:
#define EXTRA_MAPPING_NAME_POSITION 176
// String:
#define TARGET_NAME_POSITION 192

// Length definition
// int:
#define HEADER_LENGTH 1024
// int:
#define SIZE_CELL_INFO_LENGTH 64
// int:
#define SIZE_STRING 16

// Inside CellInfo
// Int:
#define NUMBER_OF_CELL_POSITION 0
// Int:
#define NUMBER_OF_COMPARTMENTS_POSITION 8
// OffSet:
#define DATA_INFO_POSITION 16
// OffSet:
#define EXTRA_MAPPING_INFO_POSITION 24
// OffSet:
#define MAPPING_INFO_POSITION 32

// order of the cell is a simple string order

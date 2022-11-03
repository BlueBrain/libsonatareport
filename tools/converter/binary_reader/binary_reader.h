#pragma once
/** For reference: git@bbpgitlab.epfl.ch:hpc/reportinglib.git */

#include <stdint.h>

#include <fstream>
#include <string>

#include <iostream>
#include <list>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

#include "binary_positions.h"
#include "cross_platform_conversion.h"

/*Terminology:
 * Compartment: lowest subdivision of cell. (cannot be broken down further)
 * Section: group of unbranched compartments (may sometimes be called segment, please change to
 * section if found)
 * Frame: Data residing in a single timestep. The scope of a frame may vary in terms of which cells
 * it includes, depending on its use.
 * See documentation for further naming conventions.
 */

namespace bbp {
namespace binary_reader {

// class list:
class Header;
class CellInfo;
class FrameInfo;
class CompartmentMapping;
class SectionData;
class Identifier;
class CellData;
class ReadManager;
class FrameParser;
// end class list

// type conventions:

typedef uint32_t CellCount;   // indications amount of cells.
typedef uint32_t CellID;      // used for GIDs
typedef uint32_t CellIndex;   // used for indexing arrays of cells.s
typedef uint32_t FrameCount;  // indictions of amount of timesteps.
typedef uint32_t FrameIndex;  // used for frame numbers
typedef int32_t FrameDiff;    // used for indicating time periods in terms of frame indexes (can be
                              // negative)
typedef double Time;          // used for indicating biological time durations.
typedef uint64_t FrameDataCount;    // used for indicating amount of data items for an entire
                                    // frame
typedef uint64_t FrameDataIndex;    // used for indexing arrays, vectors etc of data items.
typedef uint32_t CellDataCount;     // used for indicating amount of data items that reside inside
                                    // a single Cell.
typedef uint16_t SectionDataCount;  // used for indicating amount of data items that reside
                                    // inside a single section.
typedef uint16_t SectionIndex;      // used for indexing data inside single sections.
typedef float MappingItem;          // type in which SID and other mapping types are represented in.
typedef int32_t ExtraMappingItem;   // type in which extra mapping items are represented in.
typedef float DataItem;  // type in which the actual voltage, current, conductivity etc values
                         // are represented in

enum SortingKey { SID, LOC };  // sorting specifier for sections.
// SID will make sections to be sorted by their mapping,
// LOC by their locaiton in file.

namespace baseTypes {
// these typedefs define from which STL some of the clases are derived from

typedef std::set<bbp::binary_reader::CellInfo> FrameInfo;
typedef std::vector<MappingItem> CompartmentMapping;
typedef std::vector<DataItem> SectionData;
// CellData; done blow
typedef std::map<bbp::binary_reader::CellInfo, bbp::binary_reader::CellData> FrameOrganiser;
typedef std::vector<std::ios::off_type> ReadOptimiser;
}  // namespace baseTypes

// end of type conventions

class Header
{
  private:
    bool is_native;

    // header info
    uint32_t header_size;
    std::string library_version;
    std::string simulator_version;
    CellCount amount_of_cells;
    FrameDataCount total_amount_of_compartments;
    FrameCount number_of_steps;
    Time start_time;
    Time end_time;
    Time dt;
    std::string data_unit;
    std::string time_unit;
    FrameDataCount mapping_size;
    std::string mapping_name;
    FrameDataCount extra_mapping_size;
    std::string extra_mapping_name;
    std::string target_name;

  public:
    Header()
        : header_size(0) {}
    Header(char* read_buffer);
    // Header(const Header& source);

    bool operator==(const Header& h1) const;

    // get methods:
    inline bool source_is_native() const {
        return is_native;
    }
    inline std::string get_library_version() const {
        return library_version;
    }
    inline std::string get_simulator_version() const {
        return simulator_version;
    }
    inline CellCount get_number_of_cells() const {
        return amount_of_cells;
    }
    inline FrameDataCount get_total_number_of_compartments() const {
        return total_amount_of_compartments;
    }
    inline FrameCount get_number_of_steps() const {
        return number_of_steps;
    }
    inline std::string get_mapping_name() const {
        return mapping_name;
    }
    inline std::string get_extra_mapping_name() const {
        return extra_mapping_name;
    }
    inline std::string get_data_unit() const {
        return data_unit;
    }
    inline std::string get_time_unit() const {
        return time_unit;
    }
    inline FrameDataCount get_mapping_size() const {
        return mapping_size;
    }
    inline FrameDataCount get_extra_mapping_size() const {
        return extra_mapping_size;
    }
    inline std::string get_target_name() const {
        return target_name;
    }
    inline Time get_start_time() const {
        return start_time;
    }
    inline Time get_end_time() const {
        return end_time;
    }
    inline Time get_time_step_size() const {
        return dt;
    }
    // end get methods
};

class CellInfo
{
  private:
    CellID cell_num;
    CellDataCount amount_of_compartments;

    std::ios::off_type data_location;
    std::ios::off_type extra_mapping_location;
    std::ios::off_type mapping_location;

    void* master;

  public:
    enum key { GID, LOCATION };

    CellInfo();
    CellInfo(char* read_buffer, bool is_native = true, key cmp = GID);

    bool operator==(const CellInfo& c1) const;
    bool operator!=(const CellInfo& c1) const;
    bool operator<=(const CellInfo& c1) const;
    bool operator>=(const CellInfo& c1) const;
    bool operator<(const CellInfo& c1) const;
    bool operator>(const CellInfo& c1) const;

    inline std::ios::off_type get_data_location() const {
        return data_location;
    }

    inline std::ios::off_type get_mapping_location() const {
        return mapping_location;
    }

    inline std::ios::off_type get_extra_mapping_location() const {
        return extra_mapping_location;
    }

    inline CellID get_cell_num() const {
        return cell_num;
    }

    inline CellDataCount get_amount_of_compartments() const {
        return amount_of_compartments;
    }

    static bool compare_by_offset(const CellInfo& c1, const CellInfo& c2);
};

class FrameInfo: private baseTypes::FrameInfo
{
  private:
    typedef baseTypes::FrameInfo Base;
    typedef std::map<CellID, iterator> Lookup;
    Lookup gid_lookup;

  public:
    typedef Base::const_iterator const_iterator;

    FrameInfo() {}
    FrameInfo(char* buffer,
              char* end_of_buffer,
              bool is_native = true,
              CellInfo::key sort_mode = CellInfo::GID);
    FrameInfo& filter_gids(const std::vector<CellID>& gids);

    /**
     * Copy constructor.  Needed to ensure set iterators reference the correct collection
     */
    FrameInfo(const FrameInfo& rhs);

    /**
     * Assignment operator.  Needed to ensure set iterators reference the correct collection
     */
    FrameInfo& operator=(const FrameInfo& rhs);

    inline const CellInfo& operator[](CellID gid) {
        return *gid_lookup[gid];
    }

    inline const_iterator begin() const {
        return Base::begin();
    }

    inline const_iterator end() const {
        return Base::end();
    }
};

class CompartmentMapping: public baseTypes::CompartmentMapping
{
  private:
    typedef baseTypes::CompartmentMapping Base;

  public:
    typedef Base::size_type SizeType;

    CompartmentMapping();
    CompartmentMapping(MappingItem* start, MappingItem* end);
    CompartmentMapping(MappingItem* buffer, SizeType buf_size);
};

class SectionData: public baseTypes::SectionData
{
  private:
    typedef baseTypes::SectionData Base;

  public:
    class Identifier;
    SectionData();
    SectionData(DataItem first_data);
    SectionData(DataItem* start, DataItem* end);

    bool operator==(const SectionData& s1);
    bool operator!=(const SectionData& s1);

    inline SectionDataCount num_of_compartments() {
        return size();
    }

    static bool compare_by_original_location(const Identifier& c1, const Identifier& c2);
};

class SectionData::Identifier
{
  public:
  private:
    typedef SectionData::Identifier thisType;

    MappingItem sid;
    SectionIndex loc;

    SortingKey key;

  public:
    inline Identifier(const MappingItem& map,
                      const SectionIndex& location,
                      SortingKey primary = SID) {
        sid = map, loc = location, key = primary;
    }

    inline bool operator!=(const thisType& c) const {
        return (key == SID) ? sid != c.sid : loc != c.loc;
    }

    inline bool operator==(const thisType& c) const {
        return (key == SID) ? sid == c.sid : loc == c.loc;
    }

    inline bool operator<(const thisType& c) const {
        return (key == SID) ? sid < c.sid : loc < c.loc;
    }

    inline bool operator<=(const thisType& c) const {
        return (key == SID) ? sid <= c.sid : loc <= c.loc;
    }

    inline bool operator>(const thisType& c) const {
        return (key == SID) ? sid > c.sid : loc > c.loc;
    }

    inline bool operator>=(const thisType& c) const {
        return (key == SID) ? sid >= c.sid : loc >= c.loc;
    }

    inline operator MappingItem() const {
        return sid;
    }

    inline MappingItem get_sid() const {
        return sid;
    }

    inline SectionIndex get_location() const {
        return loc;
    }
};

namespace baseTypes {
typedef std::map<bbp::binary_reader::SectionData::Identifier, bbp::binary_reader::SectionData>
    CellData;
}

class CellData: public baseTypes::CellData
{
  private:
    typedef baseTypes::CellData Base;
    typedef std::map<MappingItem, iterator> Lookup;

    SortingKey sect_sort;
    Lookup sid_lookup;

  public:
    CellData();
    CellData(const SortingKey segment_sorting = SID);
    CellData(const DataItem* data_buffer, const CompartmentMapping& mapping);
    void fill_data(const DataItem* data_buffer, const CompartmentMapping& mapping);

    inline iterator find(const SectionData::Identifier& i) {
        return Base::find(i);
    }

    inline iterator find(const MappingItem& map_find) {
        Lookup::iterator match = sid_lookup.find(map_find);
        if (match == sid_lookup.end())
            return end();
        else
            return match->second;
    }

    bool operator==(const CellData& d1);
};

class ReadManager: public std::ifstream
{
  private:
    typedef std::ifstream Base;

  public:
    typedef Base::off_type FileOffset;
    typedef Base::pos_type FilePosition;

  private:
    FileOffset cell_info_origin;  // position in file where cell info starts (after header)

    bool is_native;  // indicates whether file is originated from an architecture that arranges
                     // variables in same byte order

    FrameDataCount mapping_size_floats;

    FileOffset start_location;
    FileOffset mapping_location;

    inline void seekg_start() {
        seekg(start_location, beg);
    }

    static CellInfo make_null_info();

  protected:
    inline FileOffset get_start_offset() {
        return start_location;
    }

  public:
    static const CellInfo _null_info;

    Header file_info;

    ReadManager();

    /**
     * Constructor which immediately opens a file for reading.
     *
     * @param filePath The bbp binary report file to open
     */
    explicit ReadManager(const std::string& file_path);

    /**
     * Opens the filename specified provided this object has not already opened a file
     *
     * @param file_path The bbp binary report file to open
     */
    void open(const std::string& file_path);

    // Find the cell's info which corresponds to cellNum (note: cellNum is not the cell index in
    // the file)
    CellInfo retrieve_find_cell(CellID cell_num);

    // Retrieve  global (time invariant) information of a cell (not actual data of the cell)
    CellInfo retrieve_cell_info(CellIndex cell_index, CellInfo::key sort_mode = CellInfo::LOCATION);

    // NOTE: its highly recommended to either disable cache (set to 0) or set cache size to
    // unlimited (set to -1) before using this method:
    FrameInfo retrieve_all_cell_info(CellInfo::key sort_mode = CellInfo::LOCATION);

    /*seekg*:
     * all the seekg functions below are wrapping around the original iostream seekg method.
     */

    inline void seekg_cell_info(CellIndex cell_index = 0) {
        seekg(cell_info_origin);
        if (cell_index)
            seekg(cell_index * SONATA_REPORT_SIZE_CELL_INFO_LENGTH, cur);
    }

    inline void seekg_data(const CellInfo& cell_spec, FrameIndex time_step = 0) {
        seekg(cell_spec.get_data_location() + time_step * get_frame_size_bytes(), beg);
    }

    inline void seekg_mapping(const CellInfo& cell_spec) {
        seekg(cell_spec.get_mapping_location(), beg);
    }

    inline void seekg_extra_mapping(const CellInfo& cell_spec) {
        seekg(cell_spec.get_extra_mapping_location(), beg);
    }

    inline void seekg_mapping_start() {
        seekg(mapping_location, beg);
    }

    inline void seekg_frame(FrameIndex time_step) {
        seekg(start_location + get_frame_size_bytes() * time_step);
    }

    /*read():
     * this function essentially takes 3 arguments (the third being the datatype).
     * it does the same as the read function of the base ifstream class that it is overloading
     * the only difference is that num is not neccesarily in number of bytes
     * instead it defines the mount of elemetns in the buffer (depending on the datatype of the
     * buffer).
     *
     * Important:
     * This read method will automatically fix the byte orders of the data if applicable.
     * The algorithm of fixing the byte orders depends on the size of each element in the buffer
     * array.
     * Therefore you should be carefull with passing char* buffers,
     * unless you are really intending to read single byte data.
     */

    template <class datatype>
    ReadManager& read(datatype* buffer, std::streamsize num) {
        Base::read((char*) buffer, num * sizeof(datatype));

        if (sizeof(datatype) != 1 && !is_native)  // a good compiler should not include this
                                                  // entire if statement when this function is
                                                  // reading single byte items.
            // This if is to avoid entering an empty loop. Otherwise could have been done in
            // switch below.
            for (datatype* end_of_buffer = buffer + num; buffer < end_of_buffer; buffer++)
                switch (sizeof(datatype)) {  // metamorphic switch
                case 8:                      // double, int64_t, offsets etc
                    DOUBLESWAP(*buffer);
                    break;
                case 4:  // float, int32_t, etc
                    FIX_INT(*buffer);
                    break;
                case 2:  // short
                    FIX_SHORT(*buffer);
                    break;
                default:  // anything weird that cannot be fixed.
                    // TODO: throw error exception
                    break;
                }
        return *this;
    }

    /*checkState():
     * will determine whether a reading error has occured,
     * output the according error message if there has been an error,
     * and exit execution if an error occured.
     */

    void check_state();

    /*getFrameSize:
     * returns the number of bytes that hold the data of all compartments of all cells for a
     * single timestep.
     * please note that its is in bytes, not floats!
     * use fileInfo.getTotalNumberOfCompartments() instead for number of floats.
     */

    inline std::streamsize get_frame_size_bytes() const {
        return file_info.get_total_number_of_compartments() * sizeof(DataItem);
    }

    /*getMappingSize:
     * returns the number of bytes that hold the data of the mappings of all compartmetns of all
     * cells.
     * please note that this is in bytes not floats.
     */

    inline std::streamsize get_mapping_size_bytes() const {
        return mapping_size_floats * sizeof(MappingItem);
    }

    // The following are low level c-style reading interfaces for max performance - they provide
    // unsorted, raw data.
    /*readFrameMapping:
     * will fill buffer with mappings of all compartmetns of all cells.
     * the file pointer will be automatically adjusted accordingly and then restored to original
     * position.
     * use getMappingSize() / sizeof(float) to determine the size of the buffer in number bytes
     * (NOT float).
     */

    void read_frame_mapping(MappingItem* buffer);

    /*readFrame:
     * will fill data_buffer with data from next frame.
     * use seekg_frame() to specify timestep.
     * use fileInfo.get_total_number_of_compartments() to determine the size of buffer
     */
    void read_frame(DataItem* data_buffer);

  private:
    /**
     * Upon opening a report file, read the header with metadata
     */
    void read_header();
};

class FrameParser: protected ReadManager
{
  public:
    typedef std::vector<std::vector<FrameDataIndex>> CompartmentIndexing;
    typedef std::vector<std::vector<SectionDataCount>> CompartmentCounts;

  protected:
    class FrameOrganiser: public baseTypes::FrameOrganiser
    {
      private:
        SortingKey sect_sort;

      public:
        FrameOrganiser();
        explicit FrameOrganiser(SortingKey section_sorting);

        void fill_data(const FrameInfo& info,
                       const std::vector<CompartmentMapping>& mapping,
                       const float* data_buffer);

        void make_ordered_data_buffer(DataItem* buffer) const;

        CompartmentIndexing make_ordered_offset_references() const;
        CompartmentCounts make_ordered_seg_sizes() const;
    };

    class ReadOptimiser: public baseTypes::ReadOptimiser
    {
      private:
        void update_offsets(size_t new_size);

      protected:
        FileOffset first_offset;
        FrameDataCount data_size_elements;

        size_t element_size;
        static size_t default_element_size;

      public:
        ReadOptimiser();
        explicit ReadOptimiser(std::streamsize frame_size_bytes);
        ReadOptimiser(const FrameInfo& cells_to_read,
                      std::streamsize frame_size_bytes,
                      FilePosition start_location);

        inline FileOffset get_first() const {
            return first_offset;
        }

        inline FrameDataCount get_data_size_elements() const {
            return data_size_elements;
        }

        inline void set_element_size(size_t new_size) {
            if (new_size != element_size) {
                update_offsets(new_size);
                element_size = new_size;
            }
        }

        static inline void set_default_element_size(size_t new_default) {
            default_element_size = new_default;
        }

        inline size_t get_element_size() const {
            return element_size;
        }
    };

    void read_frame(DataItem* buffer);

  private:
    DataItem* read_buffer;
    DataItem** ref_array;

    ReadOptimiser offsets;
    FrameOrganiser org;

    FrameIndex timestep;

    void create_references(const FrameInfo& cinfo, bool sort_data);

    void buffer_transfer(DataItem* target) const;

  public:
    /*Constructor:
     * Will read next frame of data and parse it into a map of CellData structure from file.
     */
    FrameParser();
    FrameParser(const std::string& file, bool sort_data = true);
    FrameParser(const std::string& file,
                const std::vector<CellID>& gid_target,
                bool sort_data = true);
    ~FrameParser();

    /*retarget:
     * Takes a list of gids in form of a vector of CellID and uses it to respecify the target
     * cells of this
     * frameparser.
     */

    inline FrameIndex simtime2index(Time time) {
        return (FrameIndex)((time - file_info.get_start_time()) / file_info.get_time_step_size());
    }

    void retarget(const std::vector<CellID>& gid_target, bool sort_data = true);

    inline const Header& get_header() const {
        return file_info;
    }

    inline FrameDataCount get_buffer_size_elements() const {
        return offsets.get_data_size_elements();
    }

    // to select frame (similar use to iterators):
    // Note that readFrameData already incorporates operator++ !
    FrameParser& operator++();                     // increment timestep
    FrameParser& operator++(int);                  // increment timestep
    FrameParser& operator--();                     // decrement timestep
    FrameParser& operator--(int);                  // decrement timestep
    FrameParser& operator+=(FrameDiff increment);  // increment timestep by scalar
    FrameParser& operator-=(FrameDiff decrement);  // decrement timestep by scalar
    FrameParser& operator=(FrameIndex new_time);   // set timestep to scalar

    // returns the current timestep:
    inline FrameIndex get_timestep() const {
        return timestep;
    }

    inline bool has_more() const {
        return timestep < file_info.get_number_of_steps();
    }

    // returns the current time that the current timestep represents.
    Time get_time() const;

    void read_frame_mapping(MappingItem* buffer);
    void read_frame_data(DataItem* buffer);  // reads data of next frame and advances to the frame
                                             // after. Data is sorted by GID and segment num

    /*get_references()
     * returns a vector of vectors of offsets that correspond to the relevant segments inside
     * the read buffer.
     * eg.
     * get_references()[12][31]; //will give the offset for segment 31 in cell 12
     * if the segment has no comparments, it will return an offset of 0.
     */
    inline CompartmentIndexing get_references() const {
        return org.make_ordered_offset_references();
    }

    /*getCompartmentCounts()
     * returns a vector of vectors like getReferences() above,
     * but instead of containing the offsets for each segment,
     * it contains the amount of compartmetns of the corresponding segment.
     */
    inline CompartmentCounts get_compartment_counts() const {
        return org.make_ordered_seg_sizes();
    }
};
}  // namespace binary_reader
}  // namespace bbp
